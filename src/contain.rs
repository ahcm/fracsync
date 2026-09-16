//! Containment scoring.
//!
//! Given a reference's selected hash set `R` and a query (sample) hash set `Q`,
//! *containment* of the reference in the sample is
//!
//! ```text
//! C(R in Q) = |R ∩ Q| / |R|
//! ```
//!
//! the fraction of the reference's content present in the sample. This is the
//! quantity measured over the selected, complexity-filtered hash sets, rather
//! than all genomic k-mers. With ideal uniform hash sampling and a nonempty
//! reference sketch, the sampled ratio estimates this selected-set containment.
//! Empty sketches are reported as zero; no empty-sketch bias correction or
//! confidence interval is implemented.
//!
//! From containment we derive an average-nucleotide-identity estimate
//! `ANI ≈ C^(1/k)` (the Mash/sourmash containment-ANI point estimate): if a
//! fraction `p` of bases match, an intact k-mer needs all `k` of its bases
//! intact, so the surviving k-mer fraction is `p^k`, inverted here.

/// Size of the intersection of sorted, deduplicated slices. Uses binary
/// searches for very unequal sizes, otherwise a linear merge.
pub fn intersection_size(a: &[u64], b: &[u64]) -> usize
{
    let (small, large) = if a.len() <= b.len() { (a, b) } else { (b, a) };
    if small.len() <= large.len() / 64
    {
        return small
            .iter()
            .filter(|h| large.binary_search(h).is_ok())
            .count();
    }
    let (mut i, mut j, mut n) = (0usize, 0usize, 0usize);
    while i < a.len() && j < b.len()
    {
        match a[i].cmp(&b[j])
        {
            std::cmp::Ordering::Less => i += 1,
            std::cmp::Ordering::Greater => j += 1,
            std::cmp::Ordering::Equal =>
            {
                n += 1;
                i += 1;
                j += 1;
            }
        }
    }
    n
}

/// A reference scored against a query sample.
pub struct Score
{
    pub name: String,
    /// Selected hashes in the reference.
    pub ref_size: usize,
    /// Selected hashes shared with the sample.
    pub shared: usize,
    /// `shared / ref_size`: fraction of the reference present in the sample.
    pub containment: f64,
    /// `shared / query_size`: fraction of the sample explained by this reference.
    pub query_containment: f64,
    /// `containment^(1/k)`: average-nucleotide-identity estimate.
    pub ani: f64,
}

/// Scores sorted, deduplicated reference and query hashes selected with the
/// same configuration. `k` must be positive. ANI is uncorrected for coverage.
pub fn score(name: &str, reference: &[u64], query: &[u64], k: usize) -> Score
{
    let shared = intersection_size(reference, query);
    let containment = if reference.is_empty()
    {
        0.0
    }
    else
    {
        shared as f64 / reference.len() as f64
    };
    let query_containment = if query.is_empty()
    {
        0.0
    }
    else
    {
        shared as f64 / query.len() as f64
    };
    let ani = if containment > 0.0
    {
        containment.powf(1.0 / k as f64)
    }
    else
    {
        0.0
    };
    Score {
        name: name.to_string(),
        ref_size: reference.len(),
        shared,
        containment,
        query_containment,
        ani,
    }
}

#[cfg(test)]
mod tests
{
    use super::*;

    #[test]
    fn intersection_is_correct()
    {
        assert_eq!(intersection_size(&[1, 2, 3, 4], &[2, 4, 6]), 2);
        assert_eq!(intersection_size(&[], &[1, 2]), 0);
        assert_eq!(intersection_size(&[1, 2, 3], &[1, 2, 3]), 3);
    }

    #[test]
    fn full_containment_gives_ani_one()
    {
        let r = vec![10, 20, 30];
        let q = vec![5, 10, 20, 30, 99];
        let s = score("r", &r, &q, 21);
        assert_eq!(s.shared, 3);
        assert!((s.containment - 1.0).abs() < 1e-12);
        assert!((s.ani - 1.0).abs() < 1e-12);
    }

    #[test]
    fn partial_containment_ani_between_zero_and_one()
    {
        let r: Vec<u64> = (0..100).collect();
        let q: Vec<u64> = (0..90).collect(); // 90% contained
        let s = score("r", &r, &q, 21);
        assert!((s.containment - 0.9).abs() < 1e-12);
        assert!(s.ani > 0.9 && s.ani < 1.0);
    }

    #[test]
    fn unequal_intersections_match_set_membership()
    {
        let large: Vec<u64> = (0..10_000).map(|x| x * 3).collect();
        for small in [vec![], vec![0, 1, 3, 6, 30_000], (0..500).collect()]
        {
            let expected = small.iter().filter(|h| large.contains(h)).count();
            assert_eq!(intersection_size(&small, &large), expected);
            assert_eq!(intersection_size(&large, &small), expected);
        }
    }
}
