//! Gzipped paired FASTQ writer.

use anyhow::{Context, Result};
use flate2::write::GzEncoder;
use flate2::Compression;
use std::io::Write;
use std::path::Path;

use crate::types::ReadPair;

/// Write paired FASTQ files from a set of read pairs.
///
/// Outputs:
///   `{output_dir}/R1.fq.gz`
///   `{output_dir}/R2.fq.gz`
///
/// Returns the paths to both files.
pub fn write_paired_fastq(pairs: &[ReadPair], output_dir: &str) -> Result<(String, String)> {
    let r1_path = Path::new(output_dir).join("R1.fq.gz");
    let r2_path = Path::new(output_dir).join("R2.fq.gz");

    let r1_file = std::fs::File::create(&r1_path)
        .with_context(|| format!("failed to create {}", r1_path.display()))?;
    let r2_file = std::fs::File::create(&r2_path)
        .with_context(|| format!("failed to create {}", r2_path.display()))?;

    // Use fast compression — these are intermediate files.
    let mut r1_gz = GzEncoder::new(std::io::BufWriter::new(r1_file), Compression::fast());
    let mut r2_gz = GzEncoder::new(std::io::BufWriter::new(r2_file), Compression::fast());

    for pair in pairs {
        // Read 1.
        write!(r1_gz, "@{}/1\n", pair.name)?;
        r1_gz.write_all(&pair.seq1)?;
        write!(r1_gz, "\n+\n")?;
        r1_gz.write_all(&pair.qual1)?;
        write!(r1_gz, "\n")?;

        // Read 2.
        write!(r2_gz, "@{}/2\n", pair.name)?;
        r2_gz.write_all(&pair.seq2)?;
        write!(r2_gz, "\n+\n")?;
        r2_gz.write_all(&pair.qual2)?;
        write!(r2_gz, "\n")?;
    }

    r1_gz.finish()?;
    r2_gz.finish()?;

    let r1_str = r1_path.to_string_lossy().to_string();
    let r2_str = r2_path.to_string_lossy().to_string();

    log::info!(
        "Wrote {} read pairs to {} and {}",
        pairs.len(),
        r1_str,
        r2_str,
    );

    Ok((r1_str, r2_str))
}
