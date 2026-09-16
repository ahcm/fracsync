//! Canonical ntHash rolling hash.
//!
//! We hash every k-mer to a 64-bit value that is *canonical* (a k-mer and its
//! reverse complement hash to the same value) and computed by a rolling
//! recurrence so a whole read is hashed in O(n) regardless of k. We do not aim
//! to reproduce any external tool's exact values — only an internally
//! consistent canonical hash. The raw minimum is not uniformly distributed;
//! selectors mix it before fractional sampling. Consistency is what matters for
//! containment estimation: a given k-mer must hash identically wherever it
//! occurs, in a reference or in a read.
//!
//! Windows that contain a non-ACGTU base are skipped (no hash emitted for them).

// ntHash seed table (forward seed per base; the complement seed is the seed of
// the complementary base).
const SEED_A: u64 = 0x3c8b_fbb3_95c6_0474;
const SEED_C: u64 = 0x3193_c185_62a0_2b4c;
const SEED_G: u64 = 0x2032_3ed0_8257_2324;
const SEED_T: u64 = 0x2955_49f5_4be2_4456;

/// Returns `(forward_seed, complement_seed)` for a base, or `None` for non-ACGTU.
///
/// Uracil (`U`) is treated as thymine, so RNA-alphabet sequences (e.g. RNAcentral
/// FASTA, direct-RNA reads) hash identically to their DNA spelling instead of
/// being silently dropped.
#[inline]
fn seeds(b: u8) -> Option<(u64, u64)>
{
    match b
    {
        b'A' | b'a' => Some((SEED_A, SEED_T)),
        b'C' | b'c' => Some((SEED_C, SEED_G)),
        b'G' | b'g' => Some((SEED_G, SEED_C)),
        b'T' | b't' | b'U' | b'u' => Some((SEED_T, SEED_A)),
        _ => None,
    }
}

/// Calls `f(pos, canonical_hash)` for every valid k-mer of `seq`, in order.
///
/// `pos` is the 0-based start offset of the k-mer. K-mers overlapping a
/// non-ACGTU base are skipped. Uses the standard ntHash forward/reverse rolling
/// recurrence; the canonical value is `min(forward, reverse)`.
pub fn for_each_canonical(seq: &[u8], k: usize, mut f: impl FnMut(usize, u64))
{
    for (pos, hash) in canonical_iter(seq, k).enumerate()
    {
        if let Some(hash) = hash
        {
            f(pos, hash);
        }
    }
}

/// Positional rolling hashes without a per-sequence allocation.
pub(crate) fn canonical_iter(seq: &[u8], k: usize) -> impl Iterator<Item = Option<u64>> + '_
{
    let mut fwd = 0u64;
    let mut rev = 0u64;
    let mut filled = 0usize;
    let mut i = 0usize;
    std::iter::from_fn(move || {
        if k == 0 || seq.len() < k
        {
            return None;
        }
        while i < seq.len()
        {
            match seeds(seq[i])
            {
                None =>
                {
                    filled = 0;
                    fwd = 0;
                    rev = 0;
                }
                Some((sf, sr)) =>
                {
                    if filled < k
                    {
                        fwd = fwd.rotate_left(1) ^ sf;
                        rev ^= sr.rotate_left(filled as u32);
                        filled += 1;
                    }
                    else
                    {
                        let (sf_out, sr_out) = seeds(seq[i - k]).expect("previous base valid");
                        fwd = fwd.rotate_left(1) ^ sf_out.rotate_left(k as u32) ^ sf;
                        rev = rev.rotate_right(1)
                            ^ sr_out.rotate_right(1)
                            ^ sr.rotate_left((k - 1) as u32);
                    }
                }
            }
            i += 1;
            if i >= k
            {
                return Some((filled == k).then_some(fwd.min(rev)));
            }
        }
        None
    })
}

/// SplitMix64's bijective finalizer disperses the canonical minimum before
/// thresholding. This preserves strand symmetry and does not add collisions.
/// It provides approximately uniform sampling, not a cryptographic hash.
#[inline]
pub(crate) fn sampling_hash(mut hash: u64) -> u64
{
    hash = (hash ^ (hash >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
    hash = (hash ^ (hash >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
    hash ^ (hash >> 31)
}

/// Returns `Some(hash)` per k-mer start position, `None` where a window is
/// invalid (contains non-ACGTU). Length is `seq.len() + 1 - k` (empty if shorter
/// than `k`). Used by selectors that need positional alignment between k-mer and
/// s-mer hash arrays.
pub fn canonical_hashes(seq: &[u8], k: usize) -> Vec<Option<u64>>
{
    canonical_iter(seq, k).collect()
}

#[cfg(test)]
mod tests
{
    use super::*;

    fn revcomp(seq: &[u8]) -> Vec<u8>
    {
        seq.iter()
            .rev()
            .map(|&b| match b
            {
                b'A' => b'T',
                b'C' => b'G',
                b'G' => b'C',
                b'T' => b'A',
                other => other,
            })
            .collect()
    }

    fn collect(seq: &[u8], k: usize) -> Vec<(usize, u64)>
    {
        let mut v = Vec::new();
        for_each_canonical(seq, k, |p, h| v.push((p, h)));
        v
    }

    #[test]
    fn uracil_hashes_as_thymine()
    {
        // An RNA-alphabet sequence must hash identically to its DNA spelling.
        let dna = b"ACGTTGCAACGTACGTAACCGGTTACGT";
        let rna: Vec<u8> = dna
            .iter()
            .map(|&b| if b == b'T' { b'U' } else { b })
            .collect();
        assert_eq!(collect(dna, 11), collect(&rna, 11));
        assert!(!collect(&rna, 11).is_empty());
    }

    #[test]
    fn rolling_matches_recompute()
    {
        // Rolling the hash across a sequence must equal independently hashing
        // each k-mer window in isolation.
        let seq = b"ACGTTGCAACGTACGTAACCGGTTACGT";
        let k = 7;
        let rolled = collect(seq, k);
        for &(pos, h) in &rolled
        {
            let single = collect(&seq[pos..pos + k], k);
            assert_eq!(single.len(), 1);
            assert_eq!(single[0].1, h, "mismatch at pos {pos}");
        }
    }

    #[test]
    fn canonical_is_strand_symmetric()
    {
        // A k-mer and its reverse complement must produce the same canonical
        // hash, so reads from either strand map to the same reference k-mers.
        let seq = b"ACGTTGCAACGTACGTAACCGGTTACGT";
        let k = 11;
        let fwd = collect(seq, k);
        let rc = collect(&revcomp(seq), k);
        let fwd_set: std::collections::HashSet<u64> = fwd.iter().map(|x| x.1).collect();
        let rc_set: std::collections::HashSet<u64> = rc.iter().map(|x| x.1).collect();
        assert_eq!(fwd_set, rc_set);
    }

    #[test]
    fn skips_non_acgt()
    {
        // An N at position 4 invalidates every 4-mer window covering it.
        let seq = b"ACGTNACGT";
        let k = 4;
        let v = collect(seq, k);
        let positions: Vec<usize> = v.iter().map(|x| x.0).collect();
        // Valid 4-mers: [0..4) ACGT, then [5..9) ACGT. Positions 1,2,3,4 cover N.
        assert_eq!(positions, vec![0, 5]);
    }

    #[test]
    fn short_sequence_emits_nothing()
    {
        assert!(collect(b"ACG", 7).is_empty());
        assert!(canonical_hashes(b"ACG", 7).is_empty());
    }

    #[test]
    fn positional_iterator_handles_invalid_runs_and_large_k()
    {
        let seq = b"nACGTUUacgtNACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
        for k in [0, 1, 2, 21, 32, 64, 65, 200]
        {
            let got = canonical_hashes(seq, k);
            let expected: Vec<_> = if k == 0
            {
                vec![]
            }
            else
            {
                seq.windows(k)
                    .map(|window| {
                        if window.iter().any(|&b| seeds(b).is_none())
                        {
                            return None;
                        }
                        let mut fwd = 0u64;
                        let mut rev = 0u64;
                        for (j, &base) in window.iter().enumerate()
                        {
                            let (sf, sr) = seeds(base).unwrap();
                            fwd ^= sf.rotate_left((k - 1 - j) as u32);
                            rev ^= sr.rotate_left(j as u32);
                        }
                        Some(fwd.min(rev))
                    })
                    .collect()
            };
            assert_eq!(got, expected, "k={k}");
        }
    }
}
