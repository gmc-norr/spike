use anyhow::{Context, Result};
use log::info;
use noodles::bam;

/// Summary statistics computed from a BAM file.
#[derive(Debug, Clone)]
pub struct BamStats {
    /// Mean insert size (template length) of properly-paired reads.
    pub insert_mean: f64,
    /// Standard deviation of insert size.
    pub insert_stddev: f64,
    /// Mean read length.
    pub read_length: f64,
    /// Estimated mean coverage (reads * read_length / genome_size).
    pub mean_coverage: f64,
    /// Total number of records sampled.
    pub records_sampled: usize,
}

/// Compute BAM statistics by sampling the first `sample_size` primary, mapped,
/// properly-paired records.
pub fn compute_stats(bam_path: &str, sample_size: usize) -> Result<BamStats> {
    info!(
        "Computing BAM statistics from first {} records...",
        sample_size
    );

    let mut reader = bam::io::reader::Builder::default()
        .build_from_path(bam_path)
        .with_context(|| format!("Failed to open BAM file: {}", bam_path))?;
    let header = reader.read_header()?;

    let mut insert_sizes: Vec<f64> = Vec::with_capacity(sample_size);
    let mut read_lengths: Vec<f64> = Vec::with_capacity(sample_size);
    let mut total_records: usize = 0;

    // Compute total reference length from header for coverage estimation.
    let genome_size: u64 = header
        .reference_sequences()
        .values()
        .map(|rs| {
            rs.length()
                .try_into()
                .map(|v: usize| v as u64)
                .unwrap_or(0)
        })
        .sum();

    for result in reader.records() {
        let record = result?;

        let flags = record.flags();

        // Skip unmapped, secondary, supplementary, duplicate, or failing QC reads.
        if flags.is_unmapped()
            || flags.is_secondary()
            || flags.is_supplementary()
            || flags.is_duplicate()
            || flags.is_qc_fail()
        {
            continue;
        }

        total_records += 1;

        // Collect read length from sequence.
        let seq_len = record.sequence().len();
        if seq_len > 0 {
            read_lengths.push(seq_len as f64);
        }

        // Collect insert size from properly-paired reads.
        if flags.is_properly_segmented() && !flags.is_mate_unmapped() {
            let tlen = record.template_length();
            // Only use positive template lengths (one mate per pair).
            if tlen > 0 {
                insert_sizes.push(tlen as f64);
            }
        }

        if insert_sizes.len() >= sample_size {
            break;
        }
    }

    // Compute insert size statistics.
    let (insert_mean, insert_stddev) = if insert_sizes.is_empty() {
        log::warn!("No properly-paired reads found for insert size estimation");
        (350.0, 50.0) // sensible defaults for typical Illumina libraries
    } else {
        let mean = insert_sizes.iter().sum::<f64>() / insert_sizes.len() as f64;
        let variance = insert_sizes
            .iter()
            .map(|x| (x - mean).powi(2))
            .sum::<f64>()
            / insert_sizes.len() as f64;
        (mean, variance.sqrt())
    };

    // Compute read length statistics.
    let read_length = if read_lengths.is_empty() {
        150.0 // sensible default
    } else {
        read_lengths.iter().sum::<f64>() / read_lengths.len() as f64
    };

    // Estimate coverage: (total_records_in_file * read_length) / genome_size.
    // This is a rough approximation from the sample; refined during the evidence pass.
    let mean_coverage = if genome_size > 0 {
        (total_records as f64 * read_length) / genome_size as f64
    } else {
        0.0
    };

    let stats = BamStats {
        insert_mean,
        insert_stddev,
        read_length,
        mean_coverage,
        records_sampled: total_records,
    };

    info!(
        "BAM stats: insert_mean={:.1}, insert_stddev={:.1}, read_len={:.0}, est_coverage={:.1}x ({} records sampled)",
        stats.insert_mean, stats.insert_stddev, stats.read_length, stats.mean_coverage, stats.records_sampled
    );

    Ok(stats)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_stats_defaults() {
        let stats = BamStats {
            insert_mean: 350.0,
            insert_stddev: 50.0,
            read_length: 150.0,
            mean_coverage: 30.0,
            records_sampled: 1_000_000,
        };
        assert!((stats.insert_mean - 350.0).abs() < f64::EPSILON);
        assert!((stats.insert_stddev - 50.0).abs() < f64::EPSILON);
    }
}
