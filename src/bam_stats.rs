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
    let mut saw_segmented = false;

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
        saw_segmented |= flags.is_segmented();

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

        if scan_is_complete(insert_sizes.len(), total_records, saw_segmented, sample_size) {
            break;
        }
    }

    reject_single_end(bam_path, saw_segmented, total_records, sample_size)?;
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
    let mut saw_segmented = false;

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
        saw_segmented |= flags.is_segmented();

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

        if scan_is_complete(insert_sizes.len(), total_records, saw_segmented, sample_size) {
            break;
        }
    }

    reject_single_end(cram_path, saw_segmented, total_records, sample_size)?;
    finalize_stats(insert_sizes, read_lengths, total_records, genome_size)
}

/// Whether the sampling scan has seen enough and may stop.
///
/// A paired library stops once `sample_size` insert sizes are in hand, which
/// is what keeps the fragment-length estimate representative -- a cap on the
/// *record* count would truncate that sample on a BAM whose proper pairs are
/// sparse. A single-end library yields no insert size at all, so that
/// condition alone leaves end-of-file as the only stop and the scan reads the
/// whole BAM/CRAM (L19). Hence the second clause, which caps the record count
/// once the data is known to be single-end.
///
/// Single-end is *detected*, not inferred from an empty insert-size sample:
/// SAM flag `0x1` says the read's template had more than one segment, and a
/// paired library sets it on every read whether or not the pair aligned
/// properly. So a record carrying `0x1` keeps the scan going, and a window of
/// `sample_size` primary records without one is single-end data.
fn scan_is_complete(
    insert_sizes: usize,
    total_records: usize,
    saw_segmented: bool,
    sample_size: usize,
) -> bool {
    insert_sizes >= sample_size || (!saw_segmented && total_records >= sample_size)
}

/// Refuse a single-end library, once the scan has established it is one.
///
/// spike builds every simulated read from an extracted *pair* -- extraction
/// keeps only `is_properly_segmented` records, and the quality profile is
/// trained on R1/R2 -- so a single-end BAM yields no donor material at all.
/// Before this check it ran to completion anyway, emitting a handful of
/// synthetic pairs off an empty quality profile (all-zero qualities) and a
/// default fragment distribution. Failing here, after `sample_size` records
/// rather than after the whole file, says why.
///
/// The message distinguishes the two ways the scan can have ended, because
/// they carry different weight. At end of file `total_records` is every
/// primary record there is and the verdict is certain. At the cap it is a
/// window, and a file that does hold pairs but puts `sample_size` consecutive
/// primary records without `0x1` at its head would be reported the same way.
/// That is bounded rather than impossible: spike's own input is a
/// coordinate-sorted indexed BAM/CRAM, where a mixed library's single- and
/// paired-end reads interleave at every locus, so no supported input is known
/// to reach it -- but the message says what the cap cannot rule out instead
/// of claiming certainty it does not have.
fn reject_single_end(
    path: &str,
    saw_segmented: bool,
    total_records: usize,
    sample_size: usize,
) -> Result<()> {
    if total_records > 0 && !saw_segmented {
        let scanned = if total_records >= sample_size {
            format!(
                "all {} primary records examined (scan capped at {}); a paired library whose \
                 first {} primary records all lacked 0x1 would look the same from here, so \
                 this does not rule out pairs elsewhere in the file",
                total_records, sample_size, sample_size,
            )
        } else {
            format!("all {} primary records in the file", total_records)
        };
        anyhow::bail!(
            "{} holds no paired reads: SAM flag 0x1 is unset on {}. spike simulates from \
             extracted read pairs, so a single-end library gives it nothing to work with.",
            path,
            scanned,
        );
    }
    Ok(())
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

/// Choose the sample name for the simulated read group from the `@RG` `SM`
/// values of the original BAM, in header order.
///
/// Rule: the first read group that carries an `SM` wins. A BAM whose read
/// groups disagree is already multi-sample; picking one keeps the merged BAM
/// from gaining yet another sample, and the mismatch is logged. Characters
/// that are not safe in a `@RG` line or in the generated shell scripts are
/// replaced by `_`.
fn pick_sample_name(rg_samples: &[String]) -> Option<String> {
    let first = rg_samples.iter().find(|sm| !sm.is_empty())?;
    if rg_samples.iter().any(|sm| !sm.is_empty() && sm != first) {
        log::warn!(
            "BAM read groups carry more than one SM; tagging simulated reads with the first ({})",
            first,
        );
    }
    // align.sh quotes the sample name, so punctuation and spaces survive it
    // untouched -- and must, or the simulated read group carries an SM that
    // no longer matches the BAM's own and merged.bam is two-sample again.
    // A tab ends the SM field and a newline ends the @RG line, so no amount
    // of shell quoting lets either through. A backslash is no safer even
    // though it is not a control character: bwa-mem2 and minimap2 unescape
    // `\t`/`\n` inside the -R string themselves, so `LAB\tech01` becomes a
    // truncated SM plus a bogus extra field once the aligner, not the shell,
    // does the unescaping -- shell quoting has no say over that at all.
    let safe: String = first
        .chars()
        .map(|c| if c.is_control() || c == '\\' { '_' } else { c })
        .collect();
    if safe != *first {
        log::warn!("sample name {} rewritten to {} for the @RG line", first, safe);
    }
    Some(safe)
}

/// Read the sample name (`@RG` `SM`) of an alignment file, for reuse as the
/// sample of the simulated read group.
pub fn sample_name(alignment_path: &str, ref_path: Option<&str>) -> Result<Option<String>> {
    use noodles::sam::header::record::value::map::read_group::tag;

    let header = if crate::extract::is_cram(alignment_path) {
        let rp = ref_path.ok_or_else(|| {
            anyhow::anyhow!("CRAM input requires a reference FASTA (--reference)")
        })?;
        let repository = crate::extract::build_fasta_repository(rp)?;
        noodles::cram::io::reader::Builder::default()
            .set_reference_sequence_repository(repository)
            .build_from_path(alignment_path)
            .with_context(|| format!("Failed to open CRAM file: {}", alignment_path))?
            .read_header()?
    } else {
        noodles::bam::io::reader::Builder
            .build_from_path(alignment_path)
            .with_context(|| format!("Failed to open BAM file: {}", alignment_path))?
            .read_header()?
    };

    let rg_samples: Vec<String> = header
        .read_groups()
        .values()
        .map(|rg| {
            rg.other_fields()
                .get(&tag::SAMPLE)
                .map(|sm| String::from_utf8_lossy(sm).into_owned())
                .unwrap_or_default()
        })
        .collect();

    Ok(pick_sample_name(&rg_samples))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_pick_sample_name_takes_the_first_read_group() {
        assert_eq!(pick_sample_name(&[]), None);
        assert_eq!(
            pick_sample_name(&["HG002".to_string()]),
            Some("HG002".to_string())
        );
        // Multi-sample BAM: the first @RG wins, deterministically.
        assert_eq!(
            pick_sample_name(&["NA18488".to_string(), "HG002".to_string()]),
            Some("NA18488".to_string())
        );
        // An empty SM is not a sample name.
        assert_eq!(
            pick_sample_name(&[String::new(), "HG002".to_string()]),
            Some("HG002".to_string())
        );
    }

    #[test]
    fn test_pick_sample_name_strips_characters_unsafe_in_the_generated_scripts() {
        // A tab ends the SM field and a newline ends the @RG line, so no
        // amount of shell quoting lets either through align.sh's
        // -R '@RG\t...' argument.
        assert_eq!(
            pick_sample_name(&["HG\t002\n".to_string()]),
            Some("HG_002_".to_string())
        );
    }

    #[test]
    fn test_pick_sample_name_strips_a_backslash() {
        // bwa-mem2 and minimap2 unescape `\t`/`\n` inside the -R string they
        // are handed, so a literal backslash in the sample name forms a new
        // escape together with whatever letter follows it -- e.g. `LAB\tech01`
        // becomes an @RG line with SM truncated to "LAB" and a bogus "ech01"
        // field once the aligner unescapes its own `\t`. Quoting the shell
        // word (sh_quote) cannot prevent this: the aligner, not the shell,
        // does the unescaping. Measured with real minimap2 in the report.
        assert_eq!(
            pick_sample_name(&["LAB\\tech01".to_string()]),
            Some("LAB_tech01".to_string())
        );
    }

    #[test]
    fn test_pick_sample_name_keeps_a_name_the_scripts_can_quote() {
        // align.sh quotes the sample name, so a space or an apostrophe needs
        // no rewriting. Rewriting one would give the simulated reads an SM
        // that no longer matches the original read groups -- a two-sample
        // merged BAM, which reusing the BAM's own SM exists to prevent.
        assert_eq!(
            pick_sample_name(&["Patient 123's".to_string()]),
            Some("Patient 123's".to_string())
        );
    }

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

    /// Write a BAM of `n` primary, mapped, 100 bp records, every one carrying
    /// `flags` and `template_length`, into a single 10 kb `chrT`.
    ///
    /// `flags` is what each test is about: `0x0` is a single-end library,
    /// `0x41` a paired one whose reads never aligned as a proper pair, and
    /// `0x63` an ordinary proper pair.
    fn write_flat_bam(path: &std::path::Path, n: usize, flags: u16, template_length: i32) {
        use noodles::sam::alignment::io::Write as _;
        use std::num::NonZeroUsize;

        const READ_LEN: usize = 100;

        let header = noodles::sam::Header::builder()
            .add_reference_sequence(
                "chrT",
                noodles::sam::header::record::value::Map::<
                    noodles::sam::header::record::value::map::ReferenceSequence,
                >::new(NonZeroUsize::try_from(10_000).unwrap()),
            )
            .build();

        let mut writer = noodles::bam::io::writer::Builder
            .build_from_path(path)
            .unwrap();
        writer.write_header(&header).unwrap();

        for i in 0..n {
            let record = noodles::sam::alignment::RecordBuf::builder()
                .set_name(format!("r{}", i))
                .set_flags(noodles::sam::alignment::record::Flags::from(flags))
                .set_reference_sequence_id(0)
                .set_alignment_start(noodles::core::Position::new(1 + i % 1000).unwrap())
                .set_mapping_quality(
                    noodles::sam::alignment::record::MappingQuality::new(60).unwrap(),
                )
                .set_mate_reference_sequence_id(0)
                .set_mate_alignment_start(noodles::core::Position::new(1 + i % 1000).unwrap())
                .set_template_length(template_length)
                .set_sequence(noodles::sam::alignment::record_buf::Sequence::from(
                    vec![b'A'; READ_LEN],
                ))
                .set_quality_scores(
                    noodles::sam::alignment::record_buf::QualityScores::from(vec![40u8; READ_LEN]),
                )
                .build();
            writer.write_alignment_record(&header, &record).unwrap();
        }
        writer.try_finish().unwrap();
    }

    fn scratch_dir(name: &str) -> std::path::PathBuf {
        let dir = std::env::temp_dir().join(format!("spike_test_{}_{}", name, std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        dir
    }

    #[test]
    fn test_single_end_bam_stops_at_the_cap_instead_of_reading_the_whole_file() {
        // A single-end library never yields an insert size, so the only stop
        // condition the scan had was end-of-file: on a WGS BAM that is the
        // whole file (L19). The scan must stop after `sample_size` primary
        // records and say so, naming the number it examined.
        let dir = scratch_dir("bam_stats_single_end");
        let bam = dir.join("single_end.bam");
        write_flat_bam(&bam, 500, 0x0, 0);

        let err = match compute_stats(bam.to_str().unwrap(), 50, None) {
            Ok(stats) => panic!(
                "single-end BAM accepted; the scan read {} records of 500",
                stats.records_sampled
            ),
            Err(e) => e.to_string(),
        };
        assert!(
            err.contains("50 primary records"),
            "error must name the 50 records the cap allowed, not the file's 500: {}",
            err
        );
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn test_paired_bam_without_proper_pairs_is_not_mistaken_for_single_end() {
        // Detection keys on SAM flag 0x1, which every read of a paired library
        // carries whether or not its pair aligned properly -- not on "we
        // scanned a lot and found no insert sizes", which would cap a paired
        // BAM's fragment-length estimate at an unrepresentative sample (L15).
        let dir = scratch_dir("bam_stats_paired_no_proper");
        let bam = dir.join("paired_improper.bam");
        write_flat_bam(&bam, 500, 0x41, 0); // 0x1 | 0x40: paired, read 1, not proper

        let stats = compute_stats(bam.to_str().unwrap(), 50, None)
            .expect("a paired BAM is usable even when no pair aligned properly");
        assert_eq!(
            stats.records_sampled, 500,
            "the single-end cap must not fire on paired reads"
        );
        assert!(
            (stats.insert_mean - 350.0).abs() < f64::EPSILON,
            "no proper pair means the documented 350/50 fallback, got {}",
            stats.insert_mean
        );
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn test_paired_bam_still_stops_once_the_insert_size_sample_is_full() {
        // The paired cap is unchanged: 50 insert sizes is 50 records here.
        let dir = scratch_dir("bam_stats_paired_proper");
        let bam = dir.join("paired_proper.bam");
        write_flat_bam(&bam, 500, 0x63, 300); // 0x1 | 0x2 | 0x20 | 0x40

        let stats = compute_stats(bam.to_str().unwrap(), 50, None).unwrap();
        assert_eq!(stats.records_sampled, 50);
        assert!((stats.insert_mean - 300.0).abs() < f64::EPSILON);
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn test_single_end_bam_shorter_than_the_cap_does_not_claim_the_cap_bound() {
        // The scan reached end of file after 100 records, so "scan capped at
        // 50000" is not what happened and the 100 is the whole file, not a
        // sample of it. Saying the cap bound when it did not sends the reader
        // looking for the other 49900 records.
        let dir = scratch_dir("bam_stats_short_single_end");
        let bam = dir.join("short_single_end.bam");
        write_flat_bam(&bam, 100, 0x0, 0);

        let err = compute_stats(bam.to_str().unwrap(), 50_000, None)
            .expect_err("a single-end BAM is unusable however short it is")
            .to_string();
        assert!(
            err.contains("100 primary records"),
            "the error must name the 100 records it examined: {}",
            err
        );
        assert!(
            !err.contains("cap"),
            "the scan hit end of file, not the cap; the message must not claim otherwise: {}",
            err
        );
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn test_capped_single_end_scan_says_the_verdict_is_bounded_by_the_cap() {
        // When the cap *does* bind, the verdict rests on a window rather than
        // on the whole file, and the message has to admit that: 50000
        // consecutive primary records without 0x1 in a file that does hold
        // pairs would be reported as "holds no paired reads" too.
        let dir = scratch_dir("bam_stats_capped_single_end");
        let bam = dir.join("capped_single_end.bam");
        write_flat_bam(&bam, 500, 0x0, 0);

        let err = compute_stats(bam.to_str().unwrap(), 50, None)
            .expect_err("a single-end BAM is unusable")
            .to_string();
        assert!(
            err.contains("scan capped at 50"),
            "the cap bound here and the message must say so: {}",
            err
        );
        assert!(
            err.contains("elsewhere in the file"),
            "a capped verdict must say what it cannot rule out: {}",
            err
        );
        let _ = std::fs::remove_dir_all(&dir);
    }

    /// The shared two-contig CRAM fixture, written with the given per-record
    /// SAM flags, in a scratch directory of its own.
    fn two_contig_cram(tag: &str, first_flags: u16, last_flags: u16) -> (std::path::PathBuf, String, String) {
        let dir = scratch_dir(tag);
        let (fasta, cram) = crate::extract::test_fixtures::write_two_contig_cram_with_flags(
            &dir,
            first_flags,
            last_flags,
        );
        (dir, fasta, cram)
    }

    #[test]
    fn test_compute_stats_accepts_a_paired_cram() {
        // CRAM is a first-class input and its scan is a second loop, not a
        // shared one: M15 (index pruning) and L2 (container leakage) both
        // found it behaving differently from the BAM path, and nothing here
        // was covered either way.
        let (dir, fasta, cram) = two_contig_cram("bam_stats_cram_paired", 0x63, 0x93);

        let stats = compute_stats(&cram, 50, Some(&fasta)).expect("a paired CRAM is usable");
        assert_eq!(stats.records_sampled, 10, "5 pairs = 10 primary records");
        assert!(
            (stats.insert_mean - 300.0).abs() < f64::EPSILON,
            "every pair's positive template length is 300, got {}",
            stats.insert_mean
        );
        assert!(
            (stats.read_length - 100.0).abs() < f64::EPSILON,
            "every read is 100 bp, got {}",
            stats.read_length
        );
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn test_compute_stats_rejects_a_single_end_cram() {
        // The same rejection as the BAM path, through the CRAM loop's own
        // copy of the scan and its own `reject_single_end` call.
        let (dir, fasta, cram) = two_contig_cram("bam_stats_cram_single_end", 0x0, 0x0);

        let err = compute_stats(&cram, 50, Some(&fasta))
            .expect_err("a single-end CRAM gives spike nothing to extract from")
            .to_string();
        assert!(
            err.contains("holds no paired reads"),
            "the error must say why the CRAM is unusable: {}",
            err
        );
        assert!(
            err.contains("10 primary records") && !err.contains("cap"),
            "10 records is the whole file, well under the 50-record cap: {}",
            err
        );
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn test_scan_is_complete_caps_only_the_single_end_case() {
        // The BAM and CRAM loops share this decision.
        // Paired: keep going until the insert-size sample is full.
        assert!(!scan_is_complete(49, 10_000, true, 50));
        assert!(scan_is_complete(50, 60, true, 50));
        // Single-end: no 0x1 anywhere, so stop on the record count instead.
        assert!(!scan_is_complete(0, 49, false, 50));
        assert!(scan_is_complete(0, 50, false, 50));
    }
}
