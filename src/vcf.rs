//! Minimal VCF 4.2 reader/writer and reference-consistency validation.
//!
//! Scope: single-sample-agnostic sites (CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO). Genotype
//! columns are ignored on read. Multi-allelic ALT is kept as the raw comma-separated list.

use crate::genome_store::PackedGenome;
use std::io::{BufRead, Write};

#[derive(Debug, Clone, PartialEq)]
pub struct VcfRecord {
    pub chrom: String,
    /// 1-based position, as in the VCF spec.
    pub pos: u64,
    pub id: String,
    pub reference: String,
    pub alt: String,
    pub qual: Option<f64>,
    pub filter: String,
    pub info: String,
}

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct VcfValidation {
    pub records: usize,
    pub ref_matches: usize,
    pub ref_mismatches: usize,
    pub unknown_contig: usize,
    pub out_of_range: usize,
    pub snvs: usize,
    pub transitions: usize,
    pub transversions: usize,
}

impl VcfValidation {
    pub fn all_consistent(&self) -> bool {
        self.ref_mismatches == 0 && self.unknown_contig == 0 && self.out_of_range == 0
    }

    /// Transition/transversion ratio (None when no transversions).
    pub fn ti_tv(&self) -> Option<f64> {
        (self.transversions > 0).then(|| self.transitions as f64 / self.transversions as f64)
    }
}

fn bad(msg: String) -> std::io::Error {
    std::io::Error::new(std::io::ErrorKind::InvalidData, msg)
}

/// Parses VCF text, skipping `##` / `#CHROM` header lines.
pub fn read_vcf<R: BufRead>(reader: R) -> std::io::Result<Vec<VcfRecord>> {
    let mut out = Vec::new();
    for (ln, line) in reader.lines().enumerate() {
        let line = line?;
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 5 {
            return Err(bad(format!("VCF line {}: expected >= 5 tab-separated fields", ln + 1)));
        }
        let pos: u64 = f[1].parse().map_err(|_| bad(format!("VCF line {}: bad POS", ln + 1)))?;
        if pos == 0 {
            return Err(bad(format!("VCF line {}: POS must be >= 1", ln + 1)));
        }
        let get = |i: usize| f.get(i).copied().unwrap_or(".").to_string();
        out.push(VcfRecord {
            chrom: f[0].to_string(),
            pos,
            id: f[2].to_string(),
            reference: f[3].to_ascii_uppercase(),
            alt: f[4].to_ascii_uppercase(),
            qual: f.get(5).and_then(|q| q.parse().ok()),
            filter: get(6),
            info: get(7),
        });
    }
    Ok(out)
}

/// Writes a VCF 4.2 header (with contig lengths and somatic INFO definitions) and the records.
pub fn write_vcf<W: Write>(
    mut w: W,
    genome: &PackedGenome,
    records: &[VcfRecord],
) -> std::io::Result<()> {
    writeln!(w, "##fileformat=VCFv4.2")?;
    writeln!(w, "##reference=T2T-CHM13v2.0")?;
    writeln!(
        w,
        "##INFO=<ID=VAF,Number=1,Type=Float,Description=\"Somatic Variant Allele Frequency\">"
    )?;
    writeln!(
        w,
        "##INFO=<ID=ORIGIN,Number=1,Type=String,Description=\"Mutational Origin (OsmosticDamage|ReplicationError|CrisprOfftarget)\">"
    )?;
    writeln!(
        w,
        "##INFO=<ID=ORGAN,Number=1,Type=String,Description=\"Source Organ Compartment\">"
    )?;
    for c in &genome.contigs {
        writeln!(
            w,
            "##contig=<ID={},accession={},length={},assembly=T2T-CHM13v2.0>",
            c.canonical_chr_name(),
            c.name,
            c.len
        )?;
    }
    writeln!(w, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")?;
    for r in records {
        let q = r.qual.map(|q| format!("{q}")).unwrap_or_else(|| ".".into());
        writeln!(
            w,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            r.chrom, r.pos, r.id, r.reference, r.alt, q, r.filter, r.info
        )?;
    }
    Ok(())
}

fn is_transition(a: u8, b: u8) -> bool {
    matches!((a, b), (b'A', b'G') | (b'G', b'A') | (b'C', b'T') | (b'T', b'C'))
}

/// Checks each record's REF allele against the genome and tallies SNV Ti/Tv.
pub fn validate_against_genome(genome: &PackedGenome, records: &[VcfRecord]) -> VcfValidation {
    let mut v = VcfValidation::default();
    for r in records {
        v.records += 1;
        let Some(c) = genome.find_contig(&r.chrom) else {
            v.unknown_contig += 1;
            continue;
        };
        let start = r.pos - 1;
        let end = start + r.reference.len() as u64;
        if end > c.len {
            v.out_of_range += 1;
            continue;
        }
        if c.region(start, end) == r.reference.as_bytes() {
            v.ref_matches += 1;
        } else {
            v.ref_mismatches += 1;
        }
        if r.reference.len() == 1 && r.alt.len() == 1 && r.reference != r.alt {
            v.snvs += 1;
            if is_transition(r.reference.as_bytes()[0], r.alt.as_bytes()[0]) {
                v.transitions += 1;
            } else {
                v.transversions += 1;
            }
        }
    }
    v
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::genome_store::read_fasta_streaming;
    use std::io::Cursor;

    #[test]
    fn parse_validate_and_roundtrip() {
        let g = read_fasta_streaming(Cursor::new(">c1\nACGTACGTNN\n")).unwrap();
        let vcf = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
c1\t1\trs1\tA\tG\t50\tPASS\t.\n\
c1\t2\t.\tC\tA\t.\t.\t.\n\
c1\t3\t.\tT\tC\t.\t.\t.\n\
c1\t5\t.\tAC\tA\t.\t.\t.\n\
c1\t4\t.\tA\tG\t.\t.\t.\n\
c1\t20\t.\tA\tG\t.\t.\t.\n\
cX\t1\t.\tA\tG\t.\t.\t.\n";
        let recs = read_vcf(Cursor::new(vcf)).unwrap();
        assert_eq!(recs.len(), 7);
        let v = validate_against_genome(&g, &recs);
        assert_eq!(v.records, 7);
        assert_eq!(v.ref_matches, 3);
        assert_eq!(v.ref_mismatches, 2);
        assert_eq!(v.out_of_range, 1);
        assert_eq!(v.unknown_contig, 1);
        assert!(!v.all_consistent());
        assert_eq!(v.snvs, 4);
        assert!(v.ti_tv().is_some());

        let mut out = Vec::new();
        write_vcf(&mut out, &g, &recs).unwrap();
        let again = read_vcf(Cursor::new(out)).unwrap();
        assert_eq!(again, recs);
    }

    #[test]
    fn rejects_malformed() {
        assert!(read_vcf(Cursor::new("c1\t0\t.\tA\tG\n")).is_err());
        assert!(read_vcf(Cursor::new("c1\tx\t.\tA\tG\n")).is_err());
        assert!(read_vcf(Cursor::new("c1\t1\n")).is_err());
    }
}
