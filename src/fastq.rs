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

/// Refuse to write a quality string that is not a valid FASTQ quality line
/// for its read.
///
/// Two ways it can be invalid, and the length is the one that actually
/// escapes: a quality shorter than SEQ contains no out-of-range byte at all,
/// so on the CRAM path — where noodles hands spike an *empty* quality buffer
/// for a SAM `*` — a 151 bp SEQ line paired with a zero-length QUAL line
/// passed the byte check untouched and made `samtools import` abort with
/// "truncated file" (M14).
///
/// This is the last line of defense: a record with unusable quality is
/// dropped at extraction and sampled synthetic quality is clamped to this
/// same range, but this check catches anything that slipped past both — e.g.
/// from a future code path — instead of emitting a corrupt FASTQ file.
fn validate_quality(qual: &[u8], seq: &[u8], read_name: &str, mate: u8) -> Result<()> {
    if qual.len() != seq.len() {
        anyhow::bail!(
            "read {}/{} has {} quality byte(s) for {} base(s); refusing to write \
             invalid FASTQ",
            read_name,
            mate,
            qual.len(),
            seq.len(),
        );
    }
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

type FastqGz = GzEncoder<std::io::BufWriter<std::fs::File>>;

/// A gzipped FASTQ stream into `file`.
fn fastq_gz(file: std::fs::File) -> FastqGz {
    // Use fast compression — these are intermediate files.
    GzEncoder::new(std::io::BufWriter::new(file), Compression::fast())
}

/// Write `pair`'s record for `mate` (1 or 2).
fn write_record(gz: &mut FastqGz, pair: &ReadPair, mate: u8) -> std::io::Result<()> {
    let (seq, qual) = if mate == 1 { (&pair.seq1, &pair.qual1) } else { (&pair.seq2, &pair.qual2) };
    writeln!(gz, "@{}/{}", pair.name, mate)?;
    gz.write_all(seq)?;
    write!(gz, "\n+\n")?;
    gz.write_all(qual)?;
    writeln!(gz)
}

/// Finish `gz`, the stream to `path`, after writing its records gave `records`.
fn finish(gz: FastqGz, records: std::io::Result<()>, path: &Path) -> Result<()> {
    // `finish()` only flushes flate2's own internal buffer into the
    // `BufWriter` it returns; small output can still be sitting unwritten in
    // that `BufWriter`'s buffer. Without an explicit `flush()` here, a write
    // error (e.g. a full disk) surfaces only when the `BufWriter` is dropped,
    // where `Drop::flush` errors are silently discarded -- so a failed write
    // would be reported as `Ok` (L4). It is finished even after a failed
    // record, for the same reason.
    let finished = gz.finish().and_then(|mut w| w.flush());
    records
        .with_context(|| format!("failed to write {}", path.display()))
        .and(finished.with_context(|| format!("failed to finish writing {}", path.display())))
}

/// Write one mate's records of `pairs`, in order, to `file` as gzipped FASTQ.
fn write_mate(file: std::fs::File, pairs: &[ReadPair], mate: u8, path: &Path) -> Result<()> {
    let mut gz = fastq_gz(file);
    let records = pairs.iter().try_for_each(|pair| write_record(&mut gz, pair, mate));
    finish(gz, records, path)
}

/// Both mates' records, a pair at a time: on one thread each pair is read
/// once. A file stops at its own first failed record, as in [`write_mate`].
fn write_mates_in_turn(
    (r1_file, r1_path): (std::fs::File, &Path),
    (r2_file, r2_path): (std::fs::File, &Path),
    pairs: &[ReadPair],
) -> (Result<()>, Result<()>) {
    let (mut r1, mut r2) = (fastq_gz(r1_file), fastq_gz(r2_file));
    let (mut r1_records, mut r2_records) = (Ok(()), Ok(()));
    for pair in pairs {
        if r1_records.is_ok() {
            r1_records = write_record(&mut r1, pair, 1);
        }
        if r2_records.is_ok() {
            r2_records = write_record(&mut r2, pair, 2);
        }
    }
    (finish(r1, r1_records, r1_path), finish(r2, r2_records, r2_path))
}

/// Write paired FASTQ files from a set of read pairs.
///
/// Outputs:
///   `{output_dir}/R1.fq.gz`
///   `{output_dir}/R2.fq.gz`
///
/// Returns the paths to both files.
pub fn write_paired_fastq(pairs: &[ReadPair], output_dir: &str) -> Result<(String, String)> {
    // Validate everything before creating any file, so "refuse" means "wrote
    // nothing" rather than two truncated .fq.gz files in the user's --output.
    for pair in pairs {
        validate_quality(&pair.qual1, &pair.seq1, &pair.name, 1)?;
        validate_quality(&pair.qual2, &pair.seq2, &pair.name, 2)?;
    }

    let r1_path = Path::new(output_dir).join("R1.fq.gz");
    let r2_path = Path::new(output_dir).join("R2.fq.gz");

    let r1_file = std::fs::File::create(&r1_path)
        .with_context(|| format!("failed to create {}", r1_path.display()))?;
    let r2_file = std::fs::File::create(&r2_path)
        .with_context(|| format!("failed to create {}", r2_path.display()))?;

    // Each mate's file is its own stream, so on more than one thread the two
    // are written side by side; each holds the same bytes as when they are
    // written in turn. Both are finished and flushed before either error is
    // propagated: an early `?` on R1 would leave R2 unfinished, and
    // `GzEncoder::drop` discards its own error (L4).
    let (r1_result, r2_result) = if crate::extract::one_thread() {
        write_mates_in_turn((r1_file, &r1_path), (r2_file, &r2_path), pairs)
    } else {
        rayon::join(
            || write_mate(r1_file, pairs, 1, &r1_path),
            || write_mate(r2_file, pairs, 2, &r2_path),
        )
    };

    match (r1_result, r2_result) {
        (Ok(()), Ok(())) => {}
        (Err(e), Ok(())) | (Ok(()), Err(e)) => return Err(e),
        // Both failed: report both, so R1 failing first never hides that
        // R2 also failed.
        (Err(e1), Err(e2)) => return Err(e1.context(format!("{:#}", e2))),
    }

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

    fn pair_named(name: &str, seq1_len: usize, qual1: Vec<u8>, qual2: Vec<u8>) -> ReadPair {
        ReadPair {
            name: name.to_string(),
            seq1: vec![b'A'; seq1_len],
            qual1,
            seq2: vec![b'A'; qual2.len()],
            qual2,
            ref_start: 0,
            ref_end: 100,
            insert_size: 100,
            chrom: "chr1".to_string(),
            align: None,
        }
    }

    fn pair_with_qual(qual1: Vec<u8>, qual2: Vec<u8>) -> ReadPair {
        pair_named("read", qual1.len(), qual1, qual2)
    }

    fn scratch_dir(tag: &str) -> std::path::PathBuf {
        let dir = std::env::temp_dir().join(format!(
            "spike_fastq_test_{}_{}",
            tag,
            std::process::id()
        ));
        std::fs::remove_dir_all(&dir).ok();
        std::fs::create_dir_all(&dir).unwrap();
        dir
    }

    // --- M14: writing must refuse quality bytes outside Phred+33 33..=126 ---

    #[test]
    fn test_write_paired_fastq_refuses_quality_byte_outside_phred33_range() {
        let dir = scratch_dir("range");

        // Byte 32 (space) is exactly what `s.wrapping_add(33)` produced from
        // a missing (0xFF) donor quality before the extraction-time fix.
        let pairs = vec![pair_named("badqual", 10, vec![32; 10], vec![b'!' + 30; 10])];

        let result = write_paired_fastq(&pairs, dir.to_str().unwrap());

        let err = format!(
            "{:#}",
            result.expect_err("writing a quality byte outside 33-126 must be refused")
        );
        assert!(err.contains("32"), "message must name the offending byte: {}", err);
        assert!(err.contains("badqual"), "message must name the read: {}", err);
        // "Refuse" must mean "wrote nothing": no half-finished gzip left in
        // the user's --output directory.
        assert!(
            !dir.join("R1.fq.gz").exists() && !dir.join("R2.fq.gz").exists(),
            "a refused write must leave no output file behind"
        );

        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn test_write_paired_fastq_refuses_quality_shorter_than_sequence() {
        let dir = scratch_dir("length");

        // The CRAM shape of a missing quality: a full-length SEQ line and an
        // empty QUAL line. Every byte present is in range, so only a length
        // check catches it.
        let pairs = vec![pair_named("cramstar", 151, Vec::new(), vec![b'!' + 30; 151])];

        let result = write_paired_fastq(&pairs, dir.to_str().unwrap());

        let err = format!(
            "{:#}",
            result.expect_err("a quality shorter than SEQ must be refused")
        );
        assert!(err.contains("cramstar"), "message must name the read: {}", err);
        assert!(
            !dir.join("R1.fq.gz").exists() && !dir.join("R2.fq.gz").exists(),
            "a refused write must leave no output file behind"
        );

        std::fs::remove_dir_all(&dir).ok();
    }

    // --- L4: a write error on the final flush must not be silently ignored ---

    #[test]
    #[cfg(unix)]
    fn test_write_paired_fastq_reports_error_writing_to_dev_full() {
        let dir = scratch_dir("devfull");

        // /dev/full always fails a write with ENOSPC. Symlinking R1.fq.gz to
        // it means File::create opens the device itself, so any byte that
        // actually reaches the OS write() call errors — the same failure a
        // real full disk would give partway through the final gzip flush.
        std::os::unix::fs::symlink("/dev/full", dir.join("R1.fq.gz")).unwrap();

        let pairs = vec![pair_with_qual(vec![b'!' + 30; 10], vec![b'!' + 30; 10])];

        for threads in [1, 4] {
            let result = write_on(threads, &pairs, &dir);

            assert!(
                result.is_err(),
                "writing the final FASTQ to a full disk must be reported as an \
                 error, not returned as Ok (at {} thread(s))",
                threads
            );
        }

        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    #[cfg(unix)]
    fn test_write_paired_fastq_r2_flush_error_is_attributed_to_r2() {
        let dir = scratch_dir("r2devfull");

        // Only R2 is on /dev/full; R1 is a normal file. This is the path
        // the original L4 test never exercised: nothing previously proved
        // R2's own `.flush()?` actually runs and is attributed to R2 (not
        // silently reported as a bare, unattributed I/O error).
        std::os::unix::fs::symlink("/dev/full", dir.join("R2.fq.gz")).unwrap();

        let pairs = vec![pair_with_qual(vec![b'!' + 30; 10], vec![b'!' + 30; 10])];

        for threads in [1, 4] {
            let result = write_on(threads, &pairs, &dir);

            let err = format!(
                "{:#}",
                result.expect_err("writing R2's final FASTQ to a full disk must be reported as an error")
            );
            assert!(
                err.contains("R2.fq.gz"),
                "error must name R2 as the stream that failed (at {} thread(s)): {}",
                threads,
                err
            );
            assert!(
                !err.contains("R1.fq.gz"),
                "R1 succeeded and must not be blamed for R2's failure (at {} thread(s)): {}",
                threads,
                err
            );
        }

        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    #[cfg(unix)]
    fn test_write_paired_fastq_r1_failure_does_not_swallow_r2_failure() {
        let dir = scratch_dir("bothdevfull");

        // Both streams are on /dev/full. `r1_gz.finish()?.flush()?;` would
        // short-circuit on R1's error before `r2_gz` is ever finished,
        // dropping it unfinished: `GzEncoder::drop` calls `try_finish()`
        // and discards its error, and whatever it pushes into R2's
        // `BufWriter` then hits `BufWriter::drop`'s own swallowed flush —
        // silently losing R2's failure even though the function still
        // (correctly, but only by accident) returns `Err` for R1's. Both
        // streams' errors must be computed and both must be visible in the
        // reported error, not just R1's.
        std::os::unix::fs::symlink("/dev/full", dir.join("R1.fq.gz")).unwrap();
        std::os::unix::fs::symlink("/dev/full", dir.join("R2.fq.gz")).unwrap();

        let pairs = vec![pair_with_qual(vec![b'!' + 30; 10], vec![b'!' + 30; 10])];

        for threads in [1, 4] {
            let result = write_on(threads, &pairs, &dir);

            let err = format!(
                "{:#}",
                result.expect_err("writing to a full disk must be reported as an error")
            );
            assert!(
                err.contains("R1.fq.gz"),
                "error must still name R1 as one of the streams that failed (at {} thread(s)): {}",
                threads,
                err
            );
            assert!(
                err.contains("R2.fq.gz"),
                "R2's failure must not be swallowed just because R1 failed first (at {} thread(s)): {}",
                threads,
                err
            );
        }

        std::fs::remove_dir_all(&dir).ok();
    }

    /// [`write_paired_fastq`] into `dir` on a pool of `threads` threads.
    fn write_on(threads: usize, pairs: &[ReadPair], dir: &Path) -> Result<(String, String)> {
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .unwrap()
            .install(|| write_paired_fastq(pairs, dir.to_str().unwrap()))
    }

    /// The gzip bytes one mate's file must hold: every pair's record for that
    /// mate, in order, through the same encoder at the same level.
    fn expected_gz(pairs: &[ReadPair], mate: u8) -> Vec<u8> {
        let mut gz = GzEncoder::new(Vec::new(), Compression::fast());
        for pair in pairs {
            let (seq, qual) = if mate == 1 { (&pair.seq1, &pair.qual1) } else { (&pair.seq2, &pair.qual2) };
            writeln!(gz, "@{}/{}", pair.name, mate).unwrap();
            gz.write_all(seq).unwrap();
            write!(gz, "\n+\n").unwrap();
            gz.write_all(qual).unwrap();
            writeln!(gz).unwrap();
        }
        gz.finish().unwrap()
    }

    #[test]
    fn test_each_file_holds_the_same_bytes_at_any_thread_count() {
        // R1 and R2 may be written side by side; each must still be exactly
        // what one encoder writing that mate's records in order produces.
        let pairs: Vec<ReadPair> = (0..3000)
            .map(|i| {
                let mut p = pair_named(&format!("r{}", i), 151, vec![b'!' + (i % 40) as u8; 151], vec![b'!' + 30; 151]);
                p.seq2 = vec![b"ACGT"[i % 4]; 151];
                p
            })
            .collect();
        for threads in [1, 4] {
            let dir = scratch_dir(&format!("bytes{}", threads));
            rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .unwrap()
                .install(|| write_paired_fastq(&pairs, dir.to_str().unwrap()))
                .unwrap();
            assert!(std::fs::read(dir.join("R1.fq.gz")).unwrap() == expected_gz(&pairs, 1), "R1 at {} thread(s)", threads);
            assert!(std::fs::read(dir.join("R2.fq.gz")).unwrap() == expected_gz(&pairs, 2), "R2 at {} thread(s)", threads);
            std::fs::remove_dir_all(&dir).ok();
        }
    }

    #[test]
    fn test_write_paired_fastq_accepts_valid_quality_range() {
        let dir = scratch_dir("valid");

        // Boundary values 33 ('!') and 126 ('~') are both valid and must be
        // accepted.
        let pairs = vec![pair_with_qual(vec![33; 10], vec![126; 10])];

        let result = write_paired_fastq(&pairs, dir.to_str().unwrap());

        assert!(result.is_ok(), "valid boundary quality bytes must be accepted: {:?}", result.err());

        std::fs::remove_dir_all(&dir).ok();
    }
}
