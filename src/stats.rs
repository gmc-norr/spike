//! Fragment length distribution: empirical CDF construction and sampling.

use rand::Rng;

use crate::types::ReadPair;

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
    /// Only uses positive template lengths (one mate per pair convention).
    pub fn from_read_pairs(pairs: &[ReadPair]) -> Self {
        let mut sizes: Vec<i64> = pairs
            .iter()
            .map(|p| p.insert_size.abs())
            .filter(|&s| s > 0 && s < 10_000) // filter obvious outliers
            .collect();

        if sizes.is_empty() {
            log::warn!("No valid insert sizes found, using default distribution (mean=400, sd=80)");
            return Self::default_dist();
        }

        sizes.sort_unstable();

        let mean = sizes.iter().sum::<i64>() as f64 / sizes.len() as f64;
        let variance = sizes
            .iter()
            .map(|&x| (x as f64 - mean).powi(2))
            .sum::<f64>()
            / sizes.len() as f64;
        let stddev = variance.sqrt();

        log::info!(
            "Fragment distribution: mean={:.1}, stddev={:.1}, n={}",
            mean,
            stddev,
            sizes.len(),
        );

        Self {
            mean,
            stddev,
            lengths: sizes,
        }
    }

    /// Build from explicit mean and stddev (e.g. from BamStats).
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

    /// Default distribution for typical Illumina libraries.
    fn default_dist() -> Self {
        Self::from_stats(400.0, 80.0)
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
    pub fn sample_in_range<R: Rng>(&self, rng: &mut R, min: i64, max: i64) -> i64 {
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
}
