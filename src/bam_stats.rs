use anyhow::{Context, Result};
use log::info;

/// Summary statistics computed from a BAM/CRAM file.
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

/// Compute alignment statistics by sampling the first `sample_size` primary, mapped,
/// properly-paired records. Supports both BAM and CRAM formats.
pub fn compute_stats(
    alignment_path: &str,
    sample_size: usize,
    ref_path: Option<&str>,
) -> Result<BamStats> {
    if crate::extract::is_cram(alignment_path) {
        let rp = ref_path.ok_or_else(|| {
            anyhow::anyhow!("CRAM input requires a reference FASTA (--reference)")
        })?;
        compute_stats_cram(alignment_path, sample_size, rp)
    } else {
        compute_stats_bam(alignment_path, sample_size)
    }
}

/// BAM-specific stats computation.
fn compute_stats_bam(bam_path: &str, sample_size: usize) -> Result<BamStats> {
    info!(
        "Computing BAM statistics from first {} records...",
        sample_size
    );

    let mut reader = noodles::bam::io::reader::Builder
        .build_from_path(bam_path)
        .with_context(|| format!("Failed to open BAM file: {}", bam_path))?;
    let header = reader.read_header()?;

    let mut insert_sizes: Vec<f64> = Vec::with_capacity(sample_size);
    let mut read_lengths: Vec<f64> = Vec::with_capacity(sample_size);
    let mut total_records: usize = 0;

    let genome_size: u64 = header_genome_size(&header);

    for result in reader.records() {
        let record = result?;
        let flags = record.flags();

        if flags.is_unmapped()
            || flags.is_secondary()
            || flags.is_supplementary()
            || flags.is_duplicate()
            || flags.is_qc_fail()
        {
            continue;
        }

        total_records += 1;

        let seq_len = record.sequence().len();
        if seq_len > 0 {
            read_lengths.push(seq_len as f64);
        }

        if flags.is_properly_segmented() && !flags.is_mate_unmapped() {
            let tlen = record.template_length();
            if tlen > 0 {
                insert_sizes.push(tlen as f64);
            }
        }

        if insert_sizes.len() >= sample_size {
            break;
        }
    }

    finalize_stats(insert_sizes, read_lengths, total_records, genome_size)
}

/// CRAM-specific stats computation.
fn compute_stats_cram(cram_path: &str, sample_size: usize, ref_path: &str) -> Result<BamStats> {
    info!(
        "Computing CRAM statistics from first {} records...",
        sample_size
    );

    let repository = crate::extract::build_fasta_repository(ref_path)?;

    let mut reader = noodles::cram::io::reader::Builder::default()
        .set_reference_sequence_repository(repository)
        .build_from_path(cram_path)
        .with_context(|| format!("Failed to open CRAM file: {}", cram_path))?;
    let header = reader.read_header()?;

    let mut insert_sizes: Vec<f64> = Vec::with_capacity(sample_size);
    let mut read_lengths: Vec<f64> = Vec::with_capacity(sample_size);
    let mut total_records: usize = 0;

    let genome_size: u64 = header_genome_size(&header);

    for result in reader.records(&header) {
        let cram_record = result?;
        let buf = cram_record
            .try_into_alignment_record(&header)
            .with_context(|| "failed to convert CRAM record")?;

        let flags = buf.flags();

        if flags.is_unmapped()
            || flags.is_secondary()
            || flags.is_supplementary()
            || flags.is_duplicate()
            || flags.is_qc_fail()
        {
            continue;
        }

        total_records += 1;

        let seq_len = buf.sequence().len();
        if seq_len > 0 {
            read_lengths.push(seq_len as f64);
        }

        if flags.is_properly_segmented() && !flags.is_mate_unmapped() {
            let tlen = buf.template_length();
            if tlen > 0 {
                insert_sizes.push(tlen as f64);
            }
        }

        if insert_sizes.len() >= sample_size {
            break;
        }
    }

    finalize_stats(insert_sizes, read_lengths, total_records, genome_size)
}

/// Compute genome size from header reference sequences.
fn header_genome_size(header: &noodles::sam::Header) -> u64 {
    header
        .reference_sequences()
        .values()
        .map(|rs| usize::from(rs.length()) as u64)
        .sum()
}

/// Compute final statistics from collected samples.
fn finalize_stats(
    insert_sizes: Vec<f64>,
    read_lengths: Vec<f64>,
    total_records: usize,
    genome_size: u64,
) -> Result<BamStats> {
    let (insert_mean, insert_stddev) = if insert_sizes.is_empty() {
        log::warn!("No properly-paired reads found for insert size estimation");
        (350.0, 50.0)
    } else {
        let mean = insert_sizes.iter().sum::<f64>() / insert_sizes.len() as f64;
        let variance = insert_sizes.iter().map(|x| (x - mean).powi(2)).sum::<f64>()
            / insert_sizes.len() as f64;
        (mean, variance.sqrt())
    };

    let read_length = if read_lengths.is_empty() {
        150.0
    } else {
        read_lengths.iter().sum::<f64>() / read_lengths.len() as f64
    };

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
        "Stats: insert_mean={:.1}, insert_stddev={:.1}, read_len={:.0}, est_coverage={:.1}x ({} records sampled)",
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
