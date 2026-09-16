//! Strand-symmetric syncmers or window-dependent minimizers, followed by
//! low-complexity masking and FracMinHash sampling of mixed canonical hashes.
//!
//! Syncmer membership depends only on the k-mer itself, so intact k-mers are
//! reselected across read boundaries and flanking errors. Minimizers depend on
//! neighbouring k-mers and do not have this property. Plain FracMinHash
//! (`SelectorKind::All`) shares the syncmer's content-defined robustness but
//! not its spacing guarantee; it is provided so the two can be compared at
//! matched density. Scale D retains about 1/D of the selected distinct hashes;
//! it does not impose a fixed memory limit.

use std::collections::VecDeque;

use bincode::{Decode, Encode};
use serde::Serialize;

#[inline]
fn base_code(b: u8) -> u16
{
    match b
    {
        b'C' | b'c' => 1,
        b'G' | b'g' => 2,
        b'T' | b't' | b'U' | b'u' => 3,
        _ => 0,
    }
}

/// Number of distinct dinucleotides in a k-mer — a cheap linguistic-complexity
/// measure. A homopolymer has 1, a 2-base repeat 2, a codon repeat ~3, while
/// random sequence has ~12 of the 16 possible.
#[inline]
fn distinct_dinucleotides(kmer: &[u8]) -> u32
{
    let mut mask: u16 = 0;
    let mut prev = base_code(kmer[0]);
    for &b in &kmer[1..]
    {
        let c = base_code(b);
        mask |= 1 << ((prev << 2) | c);
        prev = c;
    }
    mask.count_ones()
}

/// Whether a k-mer is low-complexity (homopolymer or short-period repeat).
///
/// Poly-A tails, microsatellites, and nanopore homopolymer-length errors are
/// low-complexity and match many references spuriously, so this crate excludes
/// them in both references and reads when they fall below the threshold. The threshold
/// (fewer than k/4 distinct dinucleotides) masks repeats while leaving normal
/// sequence (~12 distinct) untouched.
#[inline]
pub fn is_low_complexity(kmer: &[u8]) -> bool
{
    if kmer.len() < 4
    {
        return false;
    }
    distinct_dinucleotides(kmer) * 4 < kmer.len() as u32
}

/// Which selection scheme and its parameters. Stored in signatures and the
/// index so query uses exactly the same selection as the database build.
#[derive(Encode, Decode, Serialize, Clone, Copy, PartialEq, Eq, Debug)]
pub enum SelectorKind
{
    /// Syncmer with accepted offsets `{offset, k-s-offset}` (zero-based).
    /// Offset zero is the closed form; a single central offset is open.
    /// Equal minima prefer the position nearest the centre; mirrored ties have
    /// the same acceptance decision. Masking, ties and scaling can remove the
    /// classical closed-syncmer window guarantee.
    ///
    /// The legacy variant name is retained for source compatibility.
    OpenSyncmer
    {
        s: u8, offset: u8
    },
    /// Minimizer: within each window of `w` consecutive k-mers keep the
    /// smallest. Included for benchmarking against syncmers.
    Minimizer
    {
        w: u8
    },
    /// Every k-mer, subject only to low-complexity masking and the scale
    /// threshold: plain FracMinHash. Content-defined like a syncmer, without
    /// the spacing guarantee. Included for benchmarking against syncmers.
    All,
}

/// A fully-specified selection configuration.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct SelectConfig
{
    pub k: usize,
    pub scale: u64,
    pub kind: SelectorKind,
}

impl SelectConfig
{
    /// Validates parameter ranges, returning a human-readable error otherwise.
    pub fn validate(&self) -> Result<(), String>
    {
        if self.k == 0 || self.k > 32
        {
            return Err(format!("k must be in 1..=32, got {}", self.k));
        }
        if self.scale == 0
        {
            return Err("scale must be >= 1".into());
        }
        match self.kind
        {
            SelectorKind::OpenSyncmer { s, offset } =>
            {
                let (s, offset) = (s as usize, offset as usize);
                if s == 0 || s >= self.k
                {
                    return Err(format!("syncmer s must be in 1..k ({}), got {s}", self.k));
                }
                let win = self.k - s + 1;
                if offset >= win
                {
                    return Err(format!("syncmer offset must be < k-s+1 ({win}), got {offset}"));
                }
            }
            SelectorKind::Minimizer { w } =>
            {
                if w == 0
                {
                    return Err("minimizer w must be >= 1".into());
                }
            }
            SelectorKind::All => {}
        }
        Ok(())
    }

    /// Expected fraction of a random sequence's k-mers this configuration
    /// selects: the selector's density times the FracMinHash retention `1/D`
    /// (the sampling hash is mixed, so it is uniform and the threshold keeps
    /// exactly that share). A k-mer shared by two sequences is selected in both
    /// or in neither, so chance coincidences between two selected sets live in
    /// a universe of `4^k / 2` canonical k-mers times this density — the
    /// quantity a chance null must be sized by.
    pub fn expected_density(&self) -> f64
    {
        let selector = match self.kind
        {
            SelectorKind::OpenSyncmer { s, offset } =>
            {
                let w = self.k - s as usize + 1;
                let mirror = w - 1 - offset as usize;
                if mirror == offset as usize { 1.0 / w as f64 } else { 2.0 / w as f64 }
            }
            SelectorKind::Minimizer { w } => 2.0 / (w as f64 + 1.0),
            SelectorKind::All => 1.0,
        };
        selector / self.scale.max(1) as f64
    }

    /// FracMinHash threshold: a hash is kept iff `hash <= max_hash()`.
    #[inline]
    fn max_hash(&self) -> u64
    {
        if self.scale <= 1
        {
            u64::MAX
        }
        else
        {
            u64::MAX / self.scale
        }
    }

    /// Emits each selected, scale-filtered canonical hash of `seq` (with
    /// possible duplicates; callers dedupe). Returns an error before emitting
    /// anything if the configuration is invalid.
    pub fn select(&self, seq: &[u8], emit: &mut impl FnMut(u64)) -> Result<(), String>
    {
        self.validate()?;
        let maxh = self.max_hash();
        match self.kind
        {
            SelectorKind::OpenSyncmer { s, offset } =>
            {
                self.select_syncmer(seq, s as usize, offset as usize, maxh, emit)
            }
            SelectorKind::Minimizer { w } => self.select_minimizer(seq, w as usize, maxh, emit),
            SelectorKind::All => self.select_all(seq, maxh, emit),
        }
        Ok(())
    }

    fn select_all(&self, seq: &[u8], maxh: u64, emit: &mut impl FnMut(u64))
    {
        crate::hash::for_each_canonical(seq, self.k, |p, kh| {
            let hash = crate::hash::sampling_hash(kh);
            if hash <= maxh && !is_low_complexity(&seq[p..p + self.k])
            {
                emit(hash);
            }
        });
    }

    fn select_syncmer(
        &self,
        seq: &[u8],
        s: usize,
        offset: usize,
        maxh: u64,
        emit: &mut impl FnMut(u64),
    )
    {
        let win = self.k - s + 1;
        let mut sh = crate::hash::canonical_iter(seq, s);
        let mut minima = VecDeque::with_capacity(win);
        let mut next_s = 0usize;
        for (p, kh) in crate::hash::canonical_iter(seq, self.k).enumerate()
        {
            expire(&mut minima, p);
            while next_s < p + win
            {
                if let Some(Some(h)) = sh.next()
                {
                    push_minimum(&mut minima, next_s, h);
                }
                next_s += 1;
            }
            let Some(kh) = kh
            else
            {
                continue;
            };
            let hash = crate::hash::sampling_hash(kh);
            if hash > maxh
            {
                continue;
            }
            // Keep equal minima in order, then choose the one closest to the
            // centre. Only ties require a scan (bounded by k <= 32).
            let min_hash = minima.front().expect("valid k-mer has valid s-mers").1;
            let best = minima
                .iter()
                .take_while(|&&(_, h)| h == min_hash)
                .min_by_key(|&&(pos, _)| (2 * (pos - p)).abs_diff(win - 1))
                .unwrap()
                .0
                - p;
            if (best == offset || best == win - 1 - offset)
                && !is_low_complexity(&seq[p..p + self.k])
            {
                emit(hash);
            }
        }
    }

    fn select_minimizer(&self, seq: &[u8], w: usize, maxh: u64, emit: &mut impl FnMut(u64))
    {
        let mut minima = VecDeque::with_capacity(w);
        let mut valid_run = 0usize;
        let mut last = None;
        for (p, kh) in crate::hash::canonical_iter(seq, self.k).enumerate()
        {
            let Some(kh) = kh
            else
            {
                minima.clear();
                valid_run = 0;
                last = None;
                continue;
            };
            valid_run += 1;
            expire(&mut minima, (p + 1).saturating_sub(w));
            push_minimum(&mut minima, p, kh);
            if valid_run < w
            {
                continue;
            }
            let &(pos, kh) = minima.front().unwrap();
            if last != Some(pos)
            {
                let hash = crate::hash::sampling_hash(kh);
                if hash <= maxh && !is_low_complexity(&seq[pos..pos + self.k])
                {
                    emit(hash);
                }
                last = Some(pos);
            }
        }
    }
}

// Monotone queues keep raw ntHash ordering separate from the mixed sampling
// hash. Retaining equal values preserves leftmost minimizer tie-breaking and
// lets the syncmer choose a centre-nearest minimum.
fn push_minimum(queue: &mut VecDeque<(usize, u64)>, pos: usize, hash: u64)
{
    while queue.back().is_some_and(|&(_, h)| h > hash)
    {
        queue.pop_back();
    }
    queue.push_back((pos, hash));
}

fn expire(queue: &mut VecDeque<(usize, u64)>, start: usize)
{
    while queue.front().is_some_and(|&(pos, _)| pos < start)
    {
        queue.pop_front();
    }
}

#[cfg(test)]
mod tests
{
    use super::*;

    fn select_set(cfg: &SelectConfig, seq: &[u8]) -> std::collections::HashSet<u64>
    {
        let mut set = std::collections::HashSet::new();
        cfg.select(seq, &mut |h| {
            set.insert(h);
        })
        .unwrap();
        set
    }

    #[test]
    fn all_selector_is_strand_consistent_and_matches_syncmer_conservation()
    {
        // Plain FracMinHash must be strandless, and — being content-defined —
        // must reselect exactly the intact k-mers a syncmer would: every hash
        // the syncmer selects at scale 1 is a k-mer `All` selects at scale 1.
        let all = SelectConfig { k: 15, scale: 1, kind: SelectorKind::All };
        let sync = SelectConfig {
            k: 15,
            scale: 1,
            kind: SelectorKind::OpenSyncmer { s: 7, offset: 0 },
        };
        let seq = b"ACGTTGCAACGTACGTAACCGGTTACGTTTGGCCAATGCATGCAACGTACGATCGATCGGCTA";
        let rc: Vec<u8> = seq
            .iter()
            .rev()
            .map(|&b| match b
            {
                b'A' => b'T',
                b'C' => b'G',
                b'G' => b'C',
                b'T' => b'A',
                x => x,
            })
            .collect();
        assert_eq!(select_set(&all, seq), select_set(&all, &rc));
        assert!(select_set(&sync, seq).is_subset(&select_set(&all, seq)));
        assert!(select_set(&sync, seq).len() < select_set(&all, seq).len());
        let scaled = SelectConfig { scale: 3, ..all };
        assert!(select_set(&scaled, seq).is_subset(&select_set(&all, seq)));
    }

    #[test]
    fn expected_density_matches_measured()
    {
        let mut x = 11u64;
        let seq: Vec<u8> = (0..400_000)
            .map(|_| {
                x = x.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
                b"ACGT"[(x >> 62) as usize]
            })
            .collect();
        let n = (seq.len() - 20) as f64;
        for (kind, scale) in [
            (SelectorKind::OpenSyncmer { s: 11, offset: 0 }, 1),
            (SelectorKind::OpenSyncmer { s: 11, offset: 5 }, 1),
            (SelectorKind::OpenSyncmer { s: 11, offset: 0 }, 7),
            (SelectorKind::Minimizer { w: 11 }, 1),
            (SelectorKind::All, 9),
        ]
        {
            let cfg = SelectConfig { k: 21, scale, kind };
            let got = select_set(&cfg, &seq).len() as f64 / n;
            let want = cfg.expected_density();
            assert!(
                (got - want).abs() < 0.08 * want,
                "{kind:?} scale {scale}: measured {got:.4}, expected {want:.4}"
            );
        }
    }

    #[test]
    fn syncmer_selection_is_strand_consistent()
    {
        // The selected hash set must be identical for a sequence and its
        // reverse complement — selection plus canonical hashing is strandless.
        let cfg = SelectConfig {
            k: 15,
            scale: 1,
            kind: SelectorKind::OpenSyncmer { s: 7, offset: 0 },
        };
        let seq = b"ACGTTGCAACGTACGTAACCGGTTACGTTTGGCCAATGCATGCAACGT";
        let rc: Vec<u8> = seq
            .iter()
            .rev()
            .map(|&b| match b
            {
                b'A' => b'T',
                b'C' => b'G',
                b'G' => b'C',
                b'T' => b'A',
                x => x,
            })
            .collect();
        assert_eq!(select_set(&cfg, seq), select_set(&cfg, &rc));
        assert!(!select_set(&cfg, seq).is_empty());
    }

    #[test]
    fn low_complexity_is_detected_and_not_selected()
    {
        // Homopolymer and short-period repeats are low-complexity.
        assert!(is_low_complexity(b"AAAAAAAAAAAAAAAAAAAAA"));
        assert!(is_low_complexity(b"ATATATATATATATATATATA"));
        assert!(is_low_complexity(b"CAGCAGCAGCAGCAGCAGCAG"));
        // Normal sequence is not.
        assert!(!is_low_complexity(b"ACGTTGCAACGTACGTAACCG"));

        // A poly-A stretch must contribute no selected k-mers.
        let cfg = SelectConfig {
            k: 21,
            scale: 1,
            kind: SelectorKind::OpenSyncmer { s: 11, offset: 0 },
        };
        assert!(select_set(&cfg, &[b'A'; 200]).is_empty());
    }

    #[test]
    fn scale_reduces_selection()
    {
        let seq = b"ACGTTGCAACGTACGTAACCGGTTACGTTTGGCCAATGCATGCAACGTACGATCGATCGGCTA";
        let base = SelectConfig {
            k: 15,
            scale: 1,
            kind: SelectorKind::OpenSyncmer { s: 7, offset: 0 },
        };
        let scaled = SelectConfig { scale: 4, ..base };
        assert!(select_set(&scaled, seq).len() <= select_set(&base, seq).len());
    }

    #[test]
    fn validate_rejects_bad_params()
    {
        assert!(
            SelectConfig {
                k: 10,
                scale: 1,
                kind: SelectorKind::OpenSyncmer { s: 10, offset: 0 }
            }
            .validate()
            .is_err()
        );
        assert!(
            SelectConfig {
                k: 10,
                scale: 1,
                kind: SelectorKind::OpenSyncmer { s: 5, offset: 99 }
            }
            .validate()
            .is_err()
        );
        assert!(
            SelectConfig {
                k: 15,
                scale: 1,
                kind: SelectorKind::OpenSyncmer { s: 7, offset: 0 }
            }
            .validate()
            .is_ok()
        );
    }

    fn generated_sequence(len: usize, alphabet: &[u8]) -> Vec<u8>
    {
        let mut state = 42u64;
        (0..len)
            .map(|_| {
                state ^= state << 13;
                state ^= state >> 7;
                state ^= state << 17;
                alphabet[state as usize % alphabet.len()]
            })
            .collect()
    }

    // Independent, deliberately simple full-window oracle for the rolling
    // queues. Compare emissions, including their order and duplicates.
    fn naive(cfg: &SelectConfig, seq: &[u8]) -> Vec<u64>
    {
        let kh = crate::hash::canonical_hashes(seq, cfg.k);
        let mut out = Vec::new();
        let mut emit = |pos: usize, raw: u64| {
            let hash = crate::hash::sampling_hash(raw);
            if hash <= u64::MAX / cfg.scale && !is_low_complexity(&seq[pos..pos + cfg.k])
            {
                out.push(hash);
            }
        };
        match cfg.kind
        {
            SelectorKind::OpenSyncmer { s, offset } =>
            {
                let win = cfg.k - s as usize + 1;
                let sh = crate::hash::canonical_hashes(seq, s as usize);
                for (p, kh) in kh.iter().enumerate()
                {
                    if let Some(h) = kh
                    {
                        let best = (0..win)
                            .min_by_key(|&t| (sh[p + t].unwrap(), (2 * t).abs_diff(win - 1)))
                            .unwrap();
                        if best == offset as usize || best == win - 1 - offset as usize
                        {
                            emit(p, *h);
                        }
                    }
                }
            }
            SelectorKind::Minimizer { w } =>
            {
                let mut last = None;
                for (p, window) in kh.windows(w as usize).enumerate()
                {
                    if window.iter().any(Option::is_none)
                    {
                        continue;
                    }
                    let (t, h) = window
                        .iter()
                        .enumerate()
                        .min_by_key(|(_, h)| h.unwrap())
                        .unwrap();
                    let pos = p + t;
                    if last != Some(pos)
                    {
                        emit(pos, h.unwrap());
                        last = Some(pos);
                    }
                }
            }
            SelectorKind::All =>
            {
                for (p, kh) in kh.iter().enumerate()
                {
                    if let Some(h) = kh
                    {
                        emit(p, *h);
                    }
                }
            }
        }
        out
    }

    #[test]
    fn rolling_selectors_match_window_oracle()
    {
        for seq in [
            generated_sequence(400, b"ACGT"),
            generated_sequence(400, b"AACCGGTUNacgt"),
            b"ACACACACACACACACACACACACACACACACACACACAC".to_vec(),
        ]
        {
            for k in [2, 5, 15, 21, 32]
            {
                let mut kinds = vec![
                    SelectorKind::Minimizer { w: 1 },
                    SelectorKind::Minimizer { w: 11 },
                    SelectorKind::Minimizer { w: 255 },
                ];
                for s in [1, k / 2, k - 1]
                {
                    for offset in 0..=k - s
                    {
                        kinds.push(SelectorKind::OpenSyncmer {
                            s: s as u8,
                            offset: offset as u8,
                        });
                    }
                }
                for kind in kinds
                {
                    for scale in [1, 7]
                    {
                        let cfg = SelectConfig { k, scale, kind };
                        let mut got = Vec::new();
                        cfg.select(&seq, &mut |h| got.push(h)).unwrap();
                        assert_eq!(got, naive(&cfg, &seq), "{cfg:?}");
                    }
                }
            }
        }
    }

    #[test]
    fn rna_and_dna_selection_match()
    {
        let dna = b"AGGTTTTTATTTGAAATTTGT";
        let rna: Vec<_> = dna
            .iter()
            .map(|&b| {
                if b == b'T'
                {
                    b'u'
                }
                else
                {
                    b.to_ascii_lowercase()
                }
            })
            .collect();
        assert_eq!(is_low_complexity(dna), is_low_complexity(&rna));
        for kind in [
            SelectorKind::Minimizer { w: 1 },
            SelectorKind::OpenSyncmer { s: 11, offset: 0 },
            SelectorKind::OpenSyncmer { s: 11, offset: 5 },
        ]
        {
            let cfg = SelectConfig {
                k: 21,
                scale: 1,
                kind,
            };
            assert_eq!(select_set(&cfg, dna), select_set(&cfg, &rna));
        }
        let cfg = SelectConfig {
            k: 21,
            scale: 1,
            kind: SelectorKind::Minimizer { w: 1 },
        };
        assert_eq!(select_set(&cfg, dna).len(), 1);
    }

    #[test]
    fn scale_retains_expected_fraction_and_nested_sets()
    {
        let seq = generated_sequence(300_000, b"ACGT");
        for kind in [
            SelectorKind::OpenSyncmer { s: 11, offset: 0 },
            SelectorKind::OpenSyncmer { s: 11, offset: 5 },
            SelectorKind::Minimizer { w: 11 },
        ]
        {
            let cfg = SelectConfig {
                k: 21,
                scale: 1,
                kind,
            };
            let full = select_set(&cfg, &seq);
            let mut previous = full.clone();
            for scale in [2, 10, 100]
            {
                let scaled = select_set(&SelectConfig { scale, ..cfg }, &seq);
                assert!(scaled.is_subset(&previous));
                let fraction = scaled.len() as f64 / full.len() as f64;
                let expected = 1.0 / scale as f64;
                assert!(
                    (fraction - expected).abs() < 0.01,
                    "{kind:?} scale={scale} fraction={fraction}"
                );
                previous = scaled;
            }
        }
    }

    #[test]
    fn invalid_config_returns_error_without_emission()
    {
        for cfg in [
            SelectConfig {
                k: 0,
                scale: 1,
                kind: SelectorKind::Minimizer { w: 1 },
            },
            SelectConfig {
                k: 33,
                scale: 1,
                kind: SelectorKind::Minimizer { w: 1 },
            },
            SelectConfig {
                k: 21,
                scale: 0,
                kind: SelectorKind::Minimizer { w: 1 },
            },
            SelectConfig {
                k: 21,
                scale: 1,
                kind: SelectorKind::Minimizer { w: 0 },
            },
            SelectConfig {
                k: 21,
                scale: 1,
                kind: SelectorKind::OpenSyncmer { s: 0, offset: 0 },
            },
            SelectConfig {
                k: 21,
                scale: 1,
                kind: SelectorKind::OpenSyncmer { s: 22, offset: 0 },
            },
            SelectConfig {
                k: 21,
                scale: 1,
                kind: SelectorKind::OpenSyncmer { s: 11, offset: 11 },
            },
        ]
        {
            assert!(
                cfg.select(&[b'A'; 100], &mut |_| panic!("invalid config emitted a hash"))
                    .is_err()
            );
        }
    }

    #[test]
    fn strand_symmetry_including_ties_and_all_offsets()
    {
        let seq = generated_sequence(1000, b"AACCGGTTN");
        let rc: Vec<_> = seq
            .iter()
            .rev()
            .map(|b| match b
            {
                b'A' => b'T',
                b'C' => b'G',
                b'G' => b'C',
                b'T' => b'A',
                x => *x,
            })
            .collect();
        for s in [1, 3, 7, 11]
        {
            for offset in 0..=15 - s
            {
                let cfg = SelectConfig {
                    k: 15,
                    scale: 1,
                    kind: SelectorKind::OpenSyncmer { s, offset },
                };
                assert_eq!(select_set(&cfg, &seq), select_set(&cfg, &rc), "{cfg:?}");
            }
        }
    }

    #[test]
    fn syncmers_reselect_intact_kmers_minimizers_need_context()
    {
        let seq = generated_sequence(200, b"ACGT");
        let cfg = SelectConfig {
            k: 21,
            scale: 1,
            kind: SelectorKind::OpenSyncmer { s: 11, offset: 0 },
        };
        let whole = select_set(&cfg, &seq);
        let mut isolated = std::collections::HashSet::new();
        for kmer in seq.windows(cfg.k)
        {
            isolated.extend(select_set(&cfg, kmer));
        }
        assert_eq!(whole, isolated);
        assert!(!whole.is_empty());
        let min = SelectConfig {
            kind: SelectorKind::Minimizer { w: 11 },
            ..cfg
        };
        assert!(!select_set(&min, &seq).is_empty());
        for kmer in seq.windows(min.k)
        {
            assert!(select_set(&min, kmer).is_empty());
        }
        assert!(select_set(&cfg, b"ACGT").is_empty());
        assert!(select_set(&min, b"ACGT").is_empty());
    }
}
