//! Linear-time suffix-array construction by induced sorting (SA-IS, Nong–Zhang–Chan 2009),
//! using `u32` indices (texts up to 2³²−2 symbols).
//!
//! Memory is roughly `4n` (SA) `+ 4n` (LMS map) `+ n/8..n` (type flags) bytes plus the input text,
//! so a 250 Mbp chromosome needs on the order of 2–3 GB. A whole 3.1 Gbp genome would need
//! ≈30+ GB and should be indexed per chromosome instead.

const NONE: u32 = u32::MAX;

/// Suffix array of `s`, whose symbols must lie in `0..=upper`. No sentinel is required.
pub fn sa_is(s: &[u32], upper: usize) -> Vec<u32> {
    let n = s.len();
    assert!(n < NONE as usize, "text too long for u32 SA-IS");
    match n {
        0 => return vec![],
        1 => return vec![0],
        2 => return if s[0] < s[1] { vec![0, 1] } else { vec![1, 0] },
        _ => {}
    }
    let mut sa = vec![0u32; n];
    // ls[i] == true  <=>  suffix i is S-type.
    let mut ls = vec![false; n];
    for i in (0..n - 1).rev() {
        ls[i] = if s[i] == s[i + 1] { ls[i + 1] } else { s[i] < s[i + 1] };
    }
    let mut sum_l = vec![0u32; upper + 1];
    let mut sum_s = vec![0u32; upper + 1];
    for i in 0..n {
        if !ls[i] {
            sum_s[s[i] as usize] += 1;
        } else {
            sum_l[s[i] as usize + 1] += 1;
        }
    }
    for i in 0..=upper {
        sum_s[i] += sum_l[i];
        if i < upper {
            sum_l[i + 1] += sum_s[i];
        }
    }

    let induce = |lms: &[u32], sa: &mut [u32]| {
        sa.iter_mut().for_each(|x| *x = NONE);
        let mut buf = sum_s.clone();
        for &d in lms {
            if d as usize == n {
                continue;
            }
            let c = s[d as usize] as usize;
            sa[buf[c] as usize] = d;
            buf[c] += 1;
        }
        buf.copy_from_slice(&sum_l);
        let c = s[n - 1] as usize;
        sa[buf[c] as usize] = (n - 1) as u32;
        buf[c] += 1;
        for i in 0..n {
            let v = sa[i];
            if v != NONE && v >= 1 && !ls[v as usize - 1] {
                let c = s[v as usize - 1] as usize;
                sa[buf[c] as usize] = v - 1;
                buf[c] += 1;
            }
        }
        buf.copy_from_slice(&sum_l);
        for i in (0..n).rev() {
            let v = sa[i];
            if v != NONE && v >= 1 && ls[v as usize - 1] {
                let c = s[v as usize - 1] as usize + 1;
                buf[c] -= 1;
                sa[buf[c] as usize] = v - 1;
            }
        }
    };

    let mut lms_map = vec![NONE; n + 1];
    let mut m = 0u32;
    for i in 1..n {
        if !ls[i - 1] && ls[i] {
            lms_map[i] = m;
            m += 1;
        }
    }
    let m = m as usize;
    let mut lms: Vec<u32> = Vec::with_capacity(m);
    for i in 1..n {
        if !ls[i - 1] && ls[i] {
            lms.push(i as u32);
        }
    }
    induce(&lms, &mut sa);

    if m > 0 {
        let mut sorted_lms: Vec<u32> = Vec::with_capacity(m);
        for &v in &sa {
            if lms_map[v as usize] != NONE {
                sorted_lms.push(v);
            }
        }
        let mut rec_s = vec![0u32; m];
        let mut rec_upper = 0usize;
        rec_s[lms_map[sorted_lms[0] as usize] as usize] = 0;
        for i in 1..m {
            let mut l = sorted_lms[i - 1] as usize;
            let mut r = sorted_lms[i] as usize;
            let end_l = if (lms_map[l] as usize) + 1 < m { lms[lms_map[l] as usize + 1] as usize } else { n };
            let end_r = if (lms_map[r] as usize) + 1 < m { lms[lms_map[r] as usize + 1] as usize } else { n };
            let mut same = true;
            if end_l - l != end_r - r {
                same = false;
            } else {
                while l < end_l {
                    if s[l] != s[r] {
                        break;
                    }
                    l += 1;
                    r += 1;
                }
                if l == n || s[l] != s[r] {
                    same = false;
                }
            }
            if !same {
                rec_upper += 1;
            }
            rec_s[lms_map[sorted_lms[i] as usize] as usize] = rec_upper as u32;
        }
        let rec_sa = sa_is(&rec_s, rec_upper);
        for i in 0..m {
            sorted_lms[i] = lms[rec_sa[i] as usize];
        }
        induce(&sorted_lms, &mut sa);
    }
    sa
}

#[cfg(test)]
mod tests {
    use super::*;

    fn brute(s: &[u32]) -> Vec<u32> {
        let mut sa: Vec<u32> = (0..s.len() as u32).collect();
        sa.sort_by(|&a, &b| s[a as usize..].cmp(&s[b as usize..]));
        sa
    }

    #[test]
    fn small_edge_cases() {
        assert_eq!(sa_is(&[], 3), Vec::<u32>::new());
        assert_eq!(sa_is(&[2], 3), vec![0]);
        assert_eq!(sa_is(&[1, 0], 3), vec![1, 0]);
        let banana = [1u32, 0, 2, 0, 2, 0]; // b a n a n a
        assert_eq!(sa_is(&banana, 2), brute(&banana));
        let same = vec![2u32; 50];
        assert_eq!(sa_is(&same, 3), brute(&same));
    }

    #[test]
    fn matches_brute_force_on_random_and_repetitive_texts() {
        let mut x = 987654321u64;
        let mut next = || {
            x = x.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            (x >> 33) as u32
        };
        for len in [3usize, 4, 5, 8, 31, 100, 500, 2000] {
            for alpha in [2u32, 4] {
                let s: Vec<u32> = (0..len).map(|_| next() % alpha).collect();
                assert_eq!(sa_is(&s, alpha as usize - 1), brute(&s), "random len {len} alpha {alpha}");
                let rep: Vec<u32> = (0..len).map(|i| (i % 3) as u32 % alpha).collect();
                assert_eq!(sa_is(&rep, alpha as usize - 1), brute(&rep), "periodic len {len}");
            }
        }
    }
}
