//! FM-index (BWT + C table + sampled Occ + sampled suffix array) over a 4-letter ACGT text.
//!
//! Construction uses prefix-doubling suffix-array sorting, `O(n log² n)` time and ~`16n` bytes of
//! temporary memory — intended for loci and chromosome *regions* (up to a few tens of Mbp), not
//! whole chromosomes. Whole-genome indexing needs a linear-time SA construction (SA-IS) and is not
//! implemented here.

const OCC_STEP: usize = 64;
const SA_STEP: usize = 16;

/// FM-index over `ACGT` text. Non-ACGT input is rejected.
#[derive(Debug, Clone)]
pub struct FmIndex {
    n: usize,
    /// BWT over codes 0..=3 plus sentinel code 4.
    bwt: Vec<u8>,
    c: [usize; 5],
    /// Occ checkpoints: for each block start, counts of 0..=3 before it.
    occ: Vec<[u32; 4]>,
    /// Suffix-array samples at rows where `row % SA_STEP == 0`.
    sa_samples: Vec<u32>,
}

fn code(b: u8) -> Option<u8> {
    match b {
        b'A' | b'a' => Some(0),
        b'C' | b'c' => Some(1),
        b'G' | b'g' => Some(2),
        b'T' | b't' => Some(3),
        _ => None,
    }
}

fn build_suffix_array(text: &[u8]) -> Vec<u32> {
    // Sentinel (4) maps to the unique smallest symbol 0; A,C,G,T map to 1..=4.
    let mapped: Vec<u32> = text.iter().map(|&c| if c == 4 { 0 } else { c as u32 + 1 }).collect();
    crate::sais::sa_is(&mapped, 4)
}

impl FmIndex {
    /// Builds an index. Returns `Err(position)` of the first non-ACGT byte.
    pub fn build(seq: &[u8]) -> Result<Self, usize> {
        let mut text = Vec::with_capacity(seq.len() + 1);
        for (i, &b) in seq.iter().enumerate() {
            text.push(code(b).ok_or(i)?);
        }
        text.push(4);
        let n = text.len();
        let sa = build_suffix_array(&text);

        let mut bwt = Vec::with_capacity(n);
        let mut counts = [0usize; 5];
        for &s in &sa {
            let ch = if s == 0 { 4 } else { text[s as usize - 1] };
            bwt.push(ch);
            counts[ch as usize] += 1;
        }
        let mut c = [0usize; 5];
        // Row 0 is the sentinel suffix; letters follow in code order.
        c[0] = 1;
        for i in 1..5 {
            c[i] = c[i - 1] + counts[i - 1];
        }

        let mut occ = Vec::with_capacity(n / OCC_STEP + 2);
        let mut run = [0u32; 4];
        for (i, &ch) in bwt.iter().enumerate() {
            if i % OCC_STEP == 0 {
                occ.push(run);
            }
            if ch < 4 {
                run[ch as usize] += 1;
            }
        }
        occ.push(run);

        let sa_samples = sa.iter().step_by(SA_STEP).copied().collect();
        Ok(Self { n, bwt, c, occ, sa_samples })
    }

    /// Text length (without sentinel).
    pub fn len(&self) -> usize {
        self.n - 1
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    fn occ_before(&self, ch: u8, i: usize) -> usize {
        let block = i / OCC_STEP;
        let mut cnt = self.occ[block][ch as usize] as usize;
        for &b in &self.bwt[block * OCC_STEP..i] {
            if b == ch {
                cnt += 1;
            }
        }
        cnt
    }

    fn range(&self, pattern: &[u8]) -> Option<(usize, usize)> {
        let (mut lo, mut hi) = (0usize, self.n);
        for &b in pattern.iter().rev() {
            let ch = code(b)?;
            lo = self.c[ch as usize] + self.occ_before(ch, lo);
            hi = self.c[ch as usize] + self.occ_before(ch, hi);
            if lo >= hi {
                return Some((0, 0));
            }
        }
        Some((lo, hi))
    }

    /// Number of occurrences of `pattern` (0 if it contains non-ACGT).
    pub fn count(&self, pattern: &[u8]) -> usize {
        self.range(pattern).map(|(l, h)| h - l).unwrap_or(0)
    }

    fn locate_row(&self, row: usize) -> usize {
        let (mut r, mut steps) = (row, 0usize);
        loop {
            if r % SA_STEP == 0 {
                return self.sa_samples[r / SA_STEP] as usize + steps;
            }
            let ch = self.bwt[r];
            if ch == 4 {
                return steps;
            }
            r = self.c[ch as usize] + self.occ_before(ch, r);
            steps += 1;
        }
    }

    /// Sorted 0-based start positions of every occurrence of `pattern`.
    pub fn locate(&self, pattern: &[u8]) -> Vec<usize> {
        let Some((lo, hi)) = self.range(pattern) else {
            return Vec::new();
        };
        let mut out: Vec<usize> = (lo..hi).map(|r| self.locate_row(r)).collect();
        out.sort_unstable();
        out
    }

    /// BWT backward-backtracking search returning `(0-based position, mismatch_count)` for all
    /// occurrences of `pattern` with `min_mismatches <= mismatches <= max_mismatches`.
    pub fn locate_within_mismatches(
        &self,
        pattern: &[u8],
        min_mismatches: u8,
        max_mismatches: u8,
    ) -> Vec<(usize, u8)> {
        if pattern.is_empty() {
            return Vec::new();
        }
        let mut pat_codes = Vec::with_capacity(pattern.len());
        for &b in pattern {
            let Some(c) = code(b) else {
                return Vec::new();
            };
            pat_codes.push(c);
        }

        let mut hits = Vec::new();
        // Stack holds `(pattern_idx_remaining, lo, hi, mismatches_so_far)`
        let mut stack = vec![(pat_codes.len(), 0usize, self.n, 0u8)];
        while let Some((idx, lo, hi, mm)) = stack.pop() {
            if lo >= hi {
                continue;
            }
            if idx == 0 {
                if mm >= min_mismatches && mm <= max_mismatches {
                    for row in lo..hi {
                        hits.push((self.locate_row(row), mm));
                    }
                }
                continue;
            }
            let target = pat_codes[idx - 1];
            for ch in 0..4u8 {
                let next_mm = mm + u8::from(ch != target);
                if next_mm <= max_mismatches {
                    let nlo = self.c[ch as usize] + self.occ_before(ch, lo);
                    let nhi = self.c[ch as usize] + self.occ_before(ch, hi);
                    if nlo < nhi {
                        stack.push((idx - 1, nlo, nhi, next_mm));
                    }
                }
            }
        }
        hits.sort_unstable_by_key(|&(pos, mm)| (pos, mm));
        hits.dedup_by_key(|&mut (pos, _)| pos);
        hits
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn naive(text: &[u8], pat: &[u8]) -> Vec<usize> {
        (0..=text.len().saturating_sub(pat.len()))
            .filter(|&i| text.len() >= pat.len() && &text[i..i + pat.len()] == pat)
            .collect()
    }

    #[test]
    fn matches_naive_search_on_pseudorandom_text() {
        let mut x = 12345u64;
        let text: Vec<u8> = (0..3000)
            .map(|_| {
                x = x.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
                b"ACGT"[(x >> 33) as usize % 4]
            })
            .collect();
        let fm = FmIndex::build(&text).unwrap();
        assert_eq!(fm.len(), 3000);
        for (s, l) in [(0, 1), (10, 4), (100, 8), (500, 12), (2990, 10), (0, 30)] {
            let pat = &text[s..s + l];
            assert_eq!(fm.locate(pat), naive(&text, pat), "pattern at {s} len {l}");
            assert_eq!(fm.count(pat), naive(&text, pat).len());
        }
        assert_eq!(fm.count(b"ACGTACGTACGTACGTACGTACGTAAAAAAAAAAAAAAAAAAAAAAAA"), 0);
        assert_eq!(fm.count(b"ACN"), 0);
    }

    #[test]
    fn repeats_and_rejects_non_acgt() {
        let fm = FmIndex::build(b"ACACACAC").unwrap();
        assert_eq!(fm.locate(b"AC"), vec![0, 2, 4, 6]);
        assert_eq!(fm.locate(b"CAC"), vec![1, 3, 5]);
        assert_eq!(FmIndex::build(b"ACNT").unwrap_err(), 2);
    }

    #[test]
    fn approximate_mismatch_search_matches_naive() {
        let text = b"ACGTACGTGGCCTTAAAGGCCCTGATCGATCGATCGATCGGGCAT";
        let fm = FmIndex::build(text).unwrap();
        let pat = b"ATCGATCGATCGGGCAT";
        let hits = fm.locate_within_mismatches(pat, 0, 3);
        let mut expected = Vec::new();
        for i in 0..=text.len() - pat.len() {
            let mm = text[i..i + pat.len()]
                .iter()
                .zip(pat.iter())
                .filter(|(a, b)| a != b)
                .count() as u8;
            if mm <= 3 {
                expected.push((i, mm));
            }
        }
        assert_eq!(hits, expected);
    }
}
