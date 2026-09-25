//! Fragment length distribution: empirical CDF construction and sampling.

use rand::Rng;

use crate::types::ReadPair;

/// Max fragment length `sample_in_range` will return. Synthetic tiling and
/// depth-copy generation both cap the sampled fragment length here, so a
/// read length above this (e.g. a long-read library) can never fit and must
/// be rejected before generation starts — see `validate_read_length` in
/// main.rs, which uses this same constant.
pub const MAX_FRAGMENT_LEN: i64 = 1500;

/// Empirical fragment length distribution for sampling.
pub struct FragmentDist {
    pub mean: f64,
    #[allow(dead_code)]
    pub stddev: f64,
    /// Sorted fragment lengths for CDF-based sampling.
    pub lengths: Vec<i64>,
}

impl FragmentDist {
    /// Build from observed insert sizes in extracted read pairs.
    ///
    /// Only uses positive template lengths (one mate per pair convention),
    /// and only those in `[read_length, MAX_FRAGMENT_LEN]` -- the range every
    /// generator call samples in (`simulate.rs`, `synth.rs`). A length outside
    /// it is never emitted, so admitting it would leave `mean` -- which
    /// `compute_tiling_count` normalises the fragment count by -- describing a
    /// wider distribution than the one generated (CR-FRAG). That also makes
    /// the range its own outlier filter: the old `< 10_000` cut is subsumed by
    /// the 1500bp cap.
    pub fn from_read_pairs(pairs: &[ReadPair], read_length: usize) -> Self {
        let min = read_length as i64;
        let sizes: Vec<i64> = pairs
            .iter()
            .map(|p| p.insert_size.abs())
            .filter(|&s| s > 0 && s >= min && s <= MAX_FRAGMENT_LEN)
            .collect();

        if sizes.is_empty() {
            log::warn!(
                "no insert size in [{}, {}] among {} read pair(s) -- the only range \
                 spike samples fragment lengths from -- so the library's own \
                 fragment lengths cannot be used; falling back to the default \
                 distribution (mean=400, sd=80) restricted to that range",
                min,
                MAX_FRAGMENT_LEN,
                pairs.len(),
            );
            return Self::default_dist_in_range(min);
        }

        let dist = Self::from_lengths(sizes);
        log::info!(
            "Fragment distribution: mean={:.1}, stddev={:.1}, n={}",
            dist.mean,
            dist.stddev,
            dist.lengths.len(),
        );
        dist
    }

    /// Mean, stddev and sorted CDF over a non-empty set of fragment lengths.
    ///
    /// Callers guarantee non-emptiness: `from_read_pairs` checks it before
    /// calling, and `default_dist_in_range` builds a fixed 10,000.
    fn from_lengths(mut lengths: Vec<i64>) -> Self {
        debug_assert!(
            !lengths.is_empty(),
            "from_lengths needs at least one length"
        );
        lengths.sort_unstable();

        let mean = lengths.iter().sum::<i64>() as f64 / lengths.len() as f64;
        let variance = lengths
            .iter()
            .map(|&x| (x as f64 - mean).powi(2))
            .sum::<f64>()
            / lengths.len() as f64;

        Self {
            mean,
            stddev: variance.sqrt(),
            lengths,
        }
    }

    /// Build from explicit mean and stddev: the fallback and test fixtures.
    /// Generates a synthetic sorted distribution by sampling from Normal.
    pub fn from_stats(mean: f64, stddev: f64) -> Self {
        use rand::rngs::StdRng;
        use rand::SeedableRng;

        let mut rng = StdRng::seed_from_u64(0);
        let mut lengths: Vec<i64> = (0..10_000)
            .map(|_| {
                let z: f64 = sample_normal(&mut rng);
                (mean + z * stddev).round().max(50.0) as i64
            })
            .collect();
        lengths.sort_unstable();

        Self {
            mean,
            stddev,
            lengths,
        }
    }

    /// Default distribution for typical Illumina libraries, held inside the
    /// range the generator samples.
    ///
    /// The stand-in for a pool with no usable insert size has to be a
    /// distribution the generator can actually draw from. A plain 400 +/- 80
    /// normal is not: its tails reach below a read length and (for a short
    /// enough library) the fallback would hand `sample_in_range` lengths it
    /// rejects 1000 times before clamping, while `mean` would normalise the
    /// tiling count by lengths never emitted -- the mismatch this range
    /// exists to close. So the default is clamped into `[min,
    /// MAX_FRAGMENT_LEN]` and its mean and stddev are recomputed over the
    /// clamped lengths, which is what a caller reading `mean` gets. Clamping
    /// rather than dropping keeps the set non-empty for every read length
    /// `validate_read_length` admits (`min <= MAX_FRAGMENT_LEN`). The caller
    /// warns first: this is a fallback, not a library measurement.
    fn default_dist_in_range(min: i64) -> Self {
        let lengths: Vec<i64> = Self::from_stats(400.0, 80.0)
            .lengths
            .into_iter()
            .map(|l| l.clamp(min, MAX_FRAGMENT_LEN))
            .collect();
        Self::from_lengths(lengths)
    }

    /// Sample a fragment length from the empirical CDF.
    pub fn sample<R: Rng>(&self, rng: &mut R) -> i64 {
        if self.lengths.is_empty() {
            return 400; // fallback
        }
        let idx = rng.gen_range(0..self.lengths.len());
        self.lengths[idx]
    }

    /// Sample with rejection: only accept lengths in [min, max].
    ///
    /// After 1000 failed attempts, clamps to the nearest valid value.
    ///
    /// Panics if `min > max`: callers pass a read length as `min`, so an
    /// empty range means the read length exceeds the max fragment length
    /// this function can return (e.g. a long-read library's reads can't fit
    /// in a short-read-sized fragment). Failing loudly with a clear message
    /// here is better than std's opaque `clamp` panic, and callers are
    /// expected to reject that read length before generation starts (see
    /// `validate_read_length` in main.rs).
    pub fn sample_in_range<R: Rng>(&self, rng: &mut R, min: i64, max: i64) -> i64 {
        assert!(
            min <= max,
            "sample_in_range: min ({min}) > max ({max}); the read length is longer \
             than the max supported fragment length",
        );
        for _ in 0..1000 {
            let len = self.sample(rng);
            if len >= min && len <= max {
                return len;
            }
        }
        // Fallback: clamp.
        self.sample(rng).clamp(min, max)
    }
}

/// Sample from standard normal distribution using Box-Muller transform.
fn sample_normal<R: Rng>(rng: &mut R) -> f64 {
    let u1: f64 = rng.gen_range(1e-10..1.0_f64);
    let u2: f64 = rng.gen_range(0.0..std::f64::consts::TAU);
    (-2.0 * u1.ln()).sqrt() * u2.cos()
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::rngs::StdRng;
    use rand::SeedableRng;

    #[test]
    fn test_fragment_dist_from_stats() {
        let dist = FragmentDist::from_stats(400.0, 80.0);
        assert_eq!(dist.lengths.len(), 10_000);
        assert!((dist.mean - 400.0).abs() < 1.0);

        let mut rng = StdRng::seed_from_u64(42);
        let sampled: Vec<i64> = (0..1000).map(|_| dist.sample(&mut rng)).collect();
        let mean = sampled.iter().sum::<i64>() as f64 / sampled.len() as f64;
        assert!((mean - 400.0).abs() < 30.0, "sampled mean was {}", mean);
    }

    #[test]
    fn test_sample_in_range() {
        let dist = FragmentDist::from_stats(400.0, 80.0);
        let mut rng = StdRng::seed_from_u64(42);
        for _ in 0..100 {
            let v = dist.sample_in_range(&mut rng, 300, 500);
            assert!(v >= 300 && v <= 500, "got {}", v);
        }
    }

    /// A donor pair carrying nothing but the insert size the model reads.
    fn pair_with_insert_size(name: &str, insert_size: i64) -> ReadPair {
        ReadPair {
            name: name.to_string(),
            seq1: b"ACGT".to_vec(),
            qual1: b"IIII".to_vec(),
            seq2: b"ACGT".to_vec(),
            qual2: b"IIII".to_vec(),
            ref_start: 1000,
            ref_end: 1000 + insert_size.unsigned_abs(),
            insert_size,
            chrom: "chr20".to_string(),
        }
    }

    #[test]
    fn test_fragment_model_keeps_only_lengths_the_generator_can_sample() {
        // CR-FRAG: every generator call samples in [read_length,
        // MAX_FRAGMENT_LEN], so an insert size outside that range is never
        // emitted and must not enter the model whose `mean` normalises the
        // tiling count. 100 is shorter than one read; 5000 and 9000 are past
        // the cap but inside the old (0, 10_000) window.
        let read_length = 150usize;
        let sizes = [400i64, 100, 9_000, 200, 5_000, 600];
        let pairs: Vec<ReadPair> = sizes
            .iter()
            .enumerate()
            .map(|(i, &s)| pair_with_insert_size(&format!("p{}", i), s))
            .collect();

        let dist = FragmentDist::from_read_pairs(&pairs, read_length);

        assert_eq!(
            dist.lengths,
            vec![200, 400, 600],
            "the model must hold only lengths in [{}, {}]",
            read_length,
            MAX_FRAGMENT_LEN,
        );
        // The mean of the in-range sizes alone, not of all six.
        assert!(
            (dist.mean - 400.0).abs() < 1e-9,
            "mean over the sampled range is 400, got {}",
            dist.mean,
        );
    }

    #[test]
    fn test_fragment_model_fallback_stays_inside_the_sampled_range() {
        // The empty case: a pool whose every insert size lies outside the
        // range the generator samples leaves nothing to model. Whatever
        // stands in must still be a distribution the generator can draw
        // from -- otherwise `sample_in_range` rejects 1000 draws and clamps,
        // and `mean` describes lengths that are never emitted.
        let read_length = 150usize;
        let pairs: Vec<ReadPair> = [3_000i64, 4_000, 5_000]
            .iter()
            .enumerate()
            .map(|(i, &s)| pair_with_insert_size(&format!("p{}", i), s))
            .collect();

        let dist = FragmentDist::from_read_pairs(&pairs, read_length);

        let min = read_length as i64;
        assert!(
            dist.lengths
                .iter()
                .all(|&l| l >= min && l <= MAX_FRAGMENT_LEN),
            "fallback holds lengths outside [{}, {}]: min {:?}, max {:?}",
            min,
            MAX_FRAGMENT_LEN,
            dist.lengths.iter().min(),
            dist.lengths.iter().max(),
        );
        assert!(
            dist.mean >= min as f64 && dist.mean <= MAX_FRAGMENT_LEN as f64,
            "fallback mean {} is outside [{}, {}]",
            dist.mean,
            min,
            MAX_FRAGMENT_LEN,
        );
    }

    #[test]
    #[should_panic(expected = "sample_in_range: min (2000) > max (1500)")]
    fn test_sample_in_range_min_gt_max_panics_with_clear_message() {
        // L5: a long-read library's read length (min) can exceed the max
        // supported fragment length (max). This must fail loudly with a
        // clear message, not with std's opaque `clamp` panic.
        let dist = FragmentDist::from_stats(400.0, 80.0);
        let mut rng = StdRng::seed_from_u64(42);
        dist.sample_in_range(&mut rng, 2000, 1500);
    }
}
