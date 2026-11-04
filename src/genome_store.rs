//! Streaming FASTA ingestion, 2-bit packed genome storage, and FASTA/BED export with checksums.
//!
//! Designed for genome-scale input (≈3.1 Gbp → ≈775 MB packed). Input is streamed line by line,
//! so memory is bounded by the packed output, never by the raw text size.
//!
//! Known lossy aspects (recorded in [`ContigStats`]): IUPAC ambiguity codes are stored as `N`,
//! and soft-masking (lowercase) is not preserved.
//! The checksum is FNV-1a 64 — an integrity check, **not** a cryptographic hash.

use crate::sha256::{to_hex, Sha256};
use std::io::{BufRead, Write};

const FNV_OFFSET: u64 = 0xcbf2_9ce4_8422_2325;
const FNV_PRIME: u64 = 0x0000_0100_0000_01b3;

/// Counts of non-trivial input handling for one contig.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct ContigStats {
    pub lowercase_bases_folded: u64,
    pub iupac_ambiguous_as_n: u64,
}

/// One chromosome/contig: 2-bit packed ACGT plus a run-length list of `N` regions.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PackedContig {
    pub name: String,
    pub len: u64,
    packed: Vec<u8>,
    /// Half-open `[start, end)` runs of `N` (stored as 0 = A in `packed`).
    pub n_runs: Vec<(u64, u64)>,
    /// Half-open `[start, end)` runs of soft-masked (lowercase) bases, preserved on export.
    pub soft_runs: Vec<(u64, u64)>,
    pub stats: ContigStats,
    /// SHA-256 over the canonical uppercase sequence (N included; case/line breaks excluded).
    pub sha256: [u8; 32],
    /// FNV-1a 64 over the canonical uppercase sequence (N included).
    pub checksum: u64,
}

impl PackedContig {
    fn code_at(&self, i: u64) -> u8 {
        (self.packed[(i / 4) as usize] >> ((i % 4) * 2)) & 3
    }

    /// Base at 0-based position `i` (`A/C/G/T/N`), or `None` if out of range.
    pub fn base_at(&self, i: u64) -> Option<u8> {
        if i >= self.len {
            return None;
        }
        if self.in_n_run(i) {
            return Some(b'N');
        }
        Some(b"ACGT"[self.code_at(i) as usize])
    }

    /// Packed payload size in bytes.
    pub fn packed_bytes(&self) -> usize {
        self.packed.len()
    }

    /// GC fraction over non-N bases.
    pub fn gc_fraction(&self) -> f64 {
        let n_total: u64 = self.n_runs.iter().map(|&(s, e)| e - s).sum();
        let valid = self.len - n_total;
        if valid == 0 {
            return 0.0;
        }
        let mut gc = 0u64;
        for i in 0..self.len {
            let c = self.code_at(i);
            if (c == 1 || c == 2) && !self.in_n_run(i) {
                gc += 1;
            }
        }
        gc as f64 / valid as f64
    }

    fn in_n_run(&self, i: u64) -> bool {
        // n_runs is sorted and non-overlapping: binary search for the last run starting <= i.
        let idx = self.n_runs.partition_point(|&(s, _)| s <= i);
        idx > 0 && i < self.n_runs[idx - 1].1
    }

    /// Checks whether any genomic coordinate in `[start, end)` falls inside a soft-masked
    /// (RepeatMasker lowercase) repetitive region.
    pub fn is_repeat_masked(&self, start: u64, end: u64) -> bool {
        if self.soft_runs.is_empty() || start >= end {
            return false;
        }
        let idx = self.soft_runs.partition_point(|&(s, _)| s < end);
        if idx == 0 {
            return false;
        }
        self.soft_runs[idx - 1].1 > start
    }

    /// Recomputes the checksum from packed data (integrity verification).
    pub fn recompute_checksum(&self) -> u64 {
        let mut h = FNV_OFFSET;
        let mut run = 0usize;
        for i in 0..self.len {
            while run < self.n_runs.len() && i >= self.n_runs[run].1 {
                run += 1;
            }
            let b = if run < self.n_runs.len() && i >= self.n_runs[run].0 {
                b'N'
            } else {
                b"ACGT"[self.code_at(i) as usize]
            };
            h = (h ^ b as u64).wrapping_mul(FNV_PRIME);
        }
        h
    }

    /// Recomputes SHA-256 from packed data.
    pub fn recompute_sha256(&self) -> [u8; 32] {
        let mut sha = Sha256::new();
        let mut run = 0usize;
        let mut chunk = Vec::with_capacity(1 << 16);
        for i in 0..self.len {
            while run < self.n_runs.len() && i >= self.n_runs[run].1 {
                run += 1;
            }
            let is_n = run < self.n_runs.len() && i >= self.n_runs[run].0;
            chunk.push(if is_n { b'N' } else { b"ACGT"[self.code_at(i) as usize] });
            if chunk.len() == 1 << 16 {
                sha.update(&chunk);
                chunk.clear();
            }
        }
        sha.update(&chunk);
        sha.finalize()
    }

    /// Hex SHA-256 of the canonical sequence.
    pub fn sha256_hex(&self) -> String {
        to_hex(&self.sha256)
    }

    pub fn verify(&self) -> bool {
        self.recompute_checksum() == self.checksum && self.recompute_sha256() == self.sha256
    }

    /// Canonical uppercase bytes of `[start, end)` (clamped), `N` included.
    pub fn region(&self, start: u64, end: u64) -> Vec<u8> {
        (start.min(self.len)..end.min(self.len)).map(|i| self.base_at(i).unwrap()).collect()
    }

    /// Returns the canonical UCSC chromosome label (e.g. `"chr1"`) if `self.name` is a known
    /// `T2T-CHM13v2.0` GenBank accession (`CP068277.2`..`CP068254.1`) or already `"chr*"`.
    pub fn canonical_chr_name(&self) -> &str {
        t2t_accession_to_chr(&self.name).unwrap_or(&self.name)
    }
}

/// Maps `T2T-CHM13v2.0` GenBank accessions (and test locus names) to canonical UCSC `chr*` names.
pub fn t2t_accession_to_chr(name: &str) -> Option<&'static str> {
    match name {
        "CP068277.2" | "chr1" | "chr1_parp1_locus" => Some("chr1"),
        "CP068276.2" | "chr2" => Some("chr2"),
        "CP068275.2" | "chr3" => Some("chr3"),
        "CP068274.2" | "chr4" => Some("chr4"),
        "CP068273.2" | "chr5" | "chr5_telomere_tert_locus" => Some("chr5"),
        "CP068272.2" | "chr6" | "chr6_foxo3_locus" => Some("chr6"),
        "CP068271.2" | "chr7" => Some("chr7"),
        "CP068270.2" | "chr8" => Some("chr8"),
        "CP068269.2" | "chr9" | "chr9_cdkn2a_locus" => Some("chr9"),
        "CP068268.2" | "chr10" => Some("chr10"),
        "CP068267.2" | "chr11" | "chr11_atm_locus" => Some("chr11"),
        "CP068266.2" | "chr12" => Some("chr12"),
        "CP068265.2" | "chr13" => Some("chr13"),
        "CP068264.2" | "chr14" => Some("chr14"),
        "CP068263.2" | "chr15" => Some("chr15"),
        "CP068262.2" | "chr16" => Some("chr16"),
        "CP068261.2" | "chr17" | "chr17_tp53_locus" => Some("chr17"),
        "CP068260.2" | "chr18" => Some("chr18"),
        "CP068259.2" | "chr19" | "chr19_sirt6_locus" => Some("chr19"),
        "CP068258.2" | "chr20" => Some("chr20"),
        "CP068257.2" | "chr21" => Some("chr21"),
        "CP068256.2" | "chr22" => Some("chr22"),
        "CP068255.2" | "chrX" => Some("chrX"),
        "CP086569.2" | "chrY" => Some("chrY"),
        "CP068254.1" | "chrM" => Some("chrM"),
        _ => None,
    }
}

/// A packed multi-contig genome.
#[derive(Debug, Clone, Default)]
pub struct PackedGenome {
    pub contigs: Vec<PackedContig>,
}

struct Builder {
    name: String,
    len: u64,
    packed: Vec<u8>,
    n_runs: Vec<(u64, u64)>,
    soft_runs: Vec<(u64, u64)>,
    open_soft: Option<u64>,
    sha: Sha256,
    open_n: Option<u64>,
    stats: ContigStats,
    hash: u64,
}

impl Builder {
    fn new(name: String) -> Self {
        Self {
            name,
            len: 0,
            packed: Vec::new(),
            n_runs: Vec::new(),
            soft_runs: Vec::new(),
            open_soft: None,
            sha: Sha256::new(),
            open_n: None,
            stats: ContigStats::default(),
            hash: FNV_OFFSET,
        }
    }

    fn push(&mut self, raw: u8) {
        if raw.is_ascii_lowercase() {
            self.stats.lowercase_bases_folded += 1;
            if self.open_soft.is_none() {
                self.open_soft = Some(self.len);
            }
        } else if let Some(s) = self.open_soft.take() {
            self.soft_runs.push((s, self.len));
        }
        let up = raw.to_ascii_uppercase();
        let (code, is_n, canon) = match up {
            b'A' => (0u8, false, b'A'),
            b'C' => (1, false, b'C'),
            b'G' => (2, false, b'G'),
            b'T' => (3, false, b'T'),
            b'N' => (0, true, b'N'),
            _ => {
                self.stats.iupac_ambiguous_as_n += 1;
                (0, true, b'N')
            }
        };
        if is_n {
            if self.open_n.is_none() {
                self.open_n = Some(self.len);
            }
        } else if let Some(s) = self.open_n.take() {
            self.n_runs.push((s, self.len));
        }
        let slot = (self.len % 4) as u32;
        if slot == 0 {
            self.packed.push(0);
        }
        *self.packed.last_mut().unwrap() |= code << (slot * 2);
        self.hash = (self.hash ^ canon as u64).wrapping_mul(FNV_PRIME);
        self.sha.update(&[canon]);
        self.len += 1;
    }

    fn finish(mut self) -> PackedContig {
        if let Some(s) = self.open_n.take() {
            self.n_runs.push((s, self.len));
        }
        if let Some(s) = self.open_soft.take() {
            self.soft_runs.push((s, self.len));
        }
        PackedContig {
            name: self.name,
            len: self.len,
            packed: self.packed,
            n_runs: self.n_runs,
            soft_runs: self.soft_runs,
            stats: self.stats,
            checksum: self.hash,
            sha256: self.sha.finalize(),
        }
    }
}

/// Streams a (multi-)FASTA from any `BufRead` into a packed genome.
/// Errors on sequence data before the first header or on I/O failure.
pub fn read_fasta_streaming<R: BufRead>(mut reader: R) -> std::io::Result<PackedGenome> {
    let mut genome = PackedGenome::default();
    let mut current: Option<Builder> = None;
    let mut line = String::new();
    loop {
        line.clear();
        if reader.read_line(&mut line)? == 0 {
            break;
        }
        let t = line.trim_end_matches(['\r', '\n']);
        if let Some(h) = t.strip_prefix('>') {
            if let Some(b) = current.take() {
                genome.contigs.push(b.finish());
            }
            let name = h.split_whitespace().next().unwrap_or("").to_string();
            current = Some(Builder::new(name));
        } else if !t.is_empty() {
            let b = current.as_mut().ok_or_else(|| {
                std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    "FASTA sequence data before first '>' header",
                )
            })?;
            for &c in t.as_bytes() {
                if !c.is_ascii_whitespace() {
                    b.push(c);
                }
            }
        }
    }
    if let Some(b) = current.take() {
        genome.contigs.push(b.finish());
    }
    Ok(genome)
}

impl PackedGenome {
    pub fn total_bp(&self) -> u64 {
        self.contigs.iter().map(|c| c.len).sum()
    }

    /// Finds a contig by exact name or canonical `T2T-CHM13v2.0` chromosome alias (`"chr1"` <-> `"CP068277.2"`).
    pub fn find_contig(&self, query: &str) -> Option<&PackedContig> {
        let canonical_q = t2t_accession_to_chr(query).unwrap_or(query);
        self.contigs
            .iter()
            .find(|c| c.name == query || c.canonical_chr_name() == canonical_q)
    }

    /// Writes a line-wrapped FASTA (uppercase, `width` columns).
    pub fn write_fasta<W: Write>(&self, mut w: W, width: usize) -> std::io::Result<()> {
        let width = width.max(1);
        for c in &self.contigs {
            writeln!(w, ">{}", c.name)?;
            let mut buf = Vec::with_capacity(width + 1);
            let (mut run, mut srun) = (0usize, 0usize);
            for i in 0..c.len {
                while run < c.n_runs.len() && i >= c.n_runs[run].1 {
                    run += 1;
                }
                while srun < c.soft_runs.len() && i >= c.soft_runs[srun].1 {
                    srun += 1;
                }
                let is_n = run < c.n_runs.len() && i >= c.n_runs[run].0;
                let soft = srun < c.soft_runs.len() && i >= c.soft_runs[srun].0;
                let b = if is_n { b'N' } else { b"ACGT"[c.code_at(i) as usize] };
                buf.push(if soft { b.to_ascii_lowercase() } else { b });
                if buf.len() == width {
                    buf.push(b'\n');
                    w.write_all(&buf)?;
                    buf.clear();
                }
            }
            if !buf.is_empty() {
                buf.push(b'\n');
                w.write_all(&buf)?;
            }
        }
        Ok(())
    }

    /// Writes BED3 + name of every `N` gap (assembly gaps).
    pub fn write_n_gaps_bed<W: Write>(&self, mut w: W) -> std::io::Result<()> {
        for c in &self.contigs {
            for &(s, e) in &c.n_runs {
                writeln!(w, "{}\t{}\t{}\tN_gap", c.name, s, e)?;
            }
        }
        Ok(())
    }

    /// Writes a tab-separated manifest: name, length, N bases, GC fraction, checksum.
    pub fn write_manifest<W: Write>(&self, mut w: W) -> std::io::Result<()> {
        writeln!(
            w,
            "contig\tlength\tn_bases\tgc_fraction\tfnv1a64\tsoft_masked_bases\tiupac_as_n\tsha256"
        )?;
        for c in &self.contigs {
            let n: u64 = c.n_runs.iter().map(|&(s, e)| e - s).sum();
            writeln!(
                w,
                "{}\t{}\t{}\t{:.6}\t{:016x}\t{}\t{}\t{}",
                c.name,
                c.len,
                n,
                c.gc_fraction(),
                c.checksum,
                c.stats.lowercase_bases_folded,
                c.stats.iupac_ambiguous_as_n,
                c.sha256_hex()
            )?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    const FA: &str = ">chrT test\nACGTNN\nNNacgR\n>chrU\nGGCC\n";

    #[test]
    fn parses_packs_and_roundtrips() {
        let g = read_fasta_streaming(Cursor::new(FA)).unwrap();
        assert_eq!(g.contigs.len(), 2);
        let c = &g.contigs[0];
        assert_eq!(c.len, 12);
        assert_eq!(c.n_runs, vec![(4, 8), (11, 12)]);
        assert_eq!(c.stats.lowercase_bases_folded, 3);
        assert_eq!(c.stats.iupac_ambiguous_as_n, 1);
        assert_eq!(c.base_at(0), Some(b'A'));
        assert_eq!(c.base_at(5), Some(b'N'));
        assert_eq!(c.base_at(9), Some(b'C'));
        assert!(g.contigs.iter().all(|c| c.verify()));
        assert_eq!(c.packed_bytes(), 3);

        let mut out = Vec::new();
        g.write_fasta(&mut out, 5).unwrap();
        let g2 = read_fasta_streaming(Cursor::new(out)).unwrap();
        assert_eq!(g2.contigs[0].checksum, g.contigs[0].checksum);
        assert_eq!(g2.total_bp(), 16);
    }

    #[test]
    fn soft_mask_preserved_and_sha256_case_insensitive() {
        let a = read_fasta_streaming(Cursor::new(">s\nACgtnnAC\nacGT\n")).unwrap();
        let b = read_fasta_streaming(Cursor::new(">s\nACGTNNACACGT\n")).unwrap();
        let c = &a.contigs[0];
        assert_eq!(c.soft_runs, vec![(2, 6), (8, 10)]);
        assert_eq!(c.sha256, b.contigs[0].sha256);
        assert!(c.verify());
        let mut out = Vec::new();
        a.write_fasta(&mut out, 60).unwrap();
        assert_eq!(String::from_utf8(out).unwrap(), ">s\nACgtnnACacGT\n");
    }

    #[test]
    fn bed_and_manifest_and_errors() {
        let g = read_fasta_streaming(Cursor::new(FA)).unwrap();
        let mut bed = Vec::new();
        g.write_n_gaps_bed(&mut bed).unwrap();
        assert_eq!(String::from_utf8(bed).unwrap(), "chrT\t4\t8\tN_gap\nchrT\t11\t12\tN_gap\n");
        let mut m = Vec::new();
        g.write_manifest(&mut m).unwrap();
        assert!(String::from_utf8(m).unwrap().contains("chrU\t4\t0\t1.000000"));
        assert!(read_fasta_streaming(Cursor::new("ACGT\n")).is_err());
    }
}
