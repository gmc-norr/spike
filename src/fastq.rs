//! Gzipped paired FASTQ writer.

use anyhow::{Context, Result};
use flate2::write::GzEncoder;
use flate2::Compression;
use std::io::Write;
use std::path::Path;

use crate::types::ReadPair;

/// Printable Phred+33 range: `!` (Q0) to `~` (Q93).
const MIN_QUAL_BYTE: u8 = 33;
const MAX_QUAL_BYTE: u8 = 126;

/// Refuse to write a quality string containing a byte outside the printable
/// Phred+33 range 33..=126.
///
/// This is the last line of defense against invalid FASTQ (M14): a BAM/CRAM
/// record with no stored quality is skipped at extraction, and sampled
/// synthetic quality is clamped to this same range, but this check catches
/// any byte that slipped past both — e.g. from a future code path — instead
/// of silently emitting a corrupt FASTQ file.
fn validate_quality(qual: &[u8], read_name: &str, mate: u8) -> Result<()> {
    if let Some(&bad) = qual.iter().find(|&&b| !(MIN_QUAL_BYTE..=MAX_QUAL_BYTE).contains(&b)) {
        anyhow::bail!(
            "read {}/{} has quality byte {} outside the printable Phred+33 range \
             {}-{}; refusing to write invalid FASTQ",
            read_name,
            mate,
            bad,
            MIN_QUAL_BYTE,
            MAX_QUAL_BYTE,
        );
    }
    Ok(())
}

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
        validate_quality(&pair.qual1, &pair.name, 1)?;
        validate_quality(&pair.qual2, &pair.name, 2)?;

        // Read 1.
        writeln!(r1_gz, "@{}/1", pair.name)?;
        r1_gz.write_all(&pair.seq1)?;
        write!(r1_gz, "\n+\n")?;
        r1_gz.write_all(&pair.qual1)?;
        writeln!(r1_gz)?;

        // Read 2.
        writeln!(r2_gz, "@{}/2", pair.name)?;
        r2_gz.write_all(&pair.seq2)?;
        write!(r2_gz, "\n+\n")?;
        r2_gz.write_all(&pair.qual2)?;
        writeln!(r2_gz)?;
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

#[cfg(test)]
mod tests {
    use super::*;

    fn pair_with_qual(qual1: Vec<u8>, qual2: Vec<u8>) -> ReadPair {
        ReadPair {
            name: "read".to_string(),
            seq1: vec![b'A'; qual1.len()],
            qual1,
            seq2: vec![b'A'; qual2.len()],
            qual2,
            ref_start: 0,
            ref_end: 100,
            insert_size: 100,
            chrom: "chr1".to_string(),
        }
    }

    // --- M14: writing must refuse quality bytes outside Phred+33 33..=126 ---

    #[test]
    fn test_write_paired_fastq_refuses_quality_byte_outside_phred33_range() {
        let dir = std::env::temp_dir().join(format!(
            "spike_fastq_test_{}",
            std::process::id()
        ));
        std::fs::create_dir_all(&dir).unwrap();

        // Byte 32 (space) is exactly what `s.wrapping_add(33)` produced from
        // a missing (0xFF) donor quality before the extraction-time fix.
        let pairs = vec![pair_with_qual(vec![32; 10], vec![b'!' + 30; 10])];

        let result = write_paired_fastq(&pairs, dir.to_str().unwrap());

        assert!(
            result.is_err(),
            "writing a quality byte outside 33-126 must be refused"
        );

        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn test_write_paired_fastq_accepts_valid_quality_range() {
        let dir = std::env::temp_dir().join(format!(
            "spike_fastq_test_valid_{}",
            std::process::id()
        ));
        std::fs::create_dir_all(&dir).unwrap();

        // Boundary values 33 ('!') and 126 ('~') are both valid and must be
        // accepted.
        let pairs = vec![pair_with_qual(vec![33; 10], vec![126; 10])];

        let result = write_paired_fastq(&pairs, dir.to_str().unwrap());

        assert!(result.is_ok(), "valid boundary quality bytes must be accepted: {:?}", result.err());

        std::fs::remove_dir_all(&dir).ok();
    }
}
