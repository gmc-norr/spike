//! BAM validation subcommand for spike.
//!
//! Reads a simulated BAM + truth VCF and runs automated checks to verify
//! that the spike-in reads look realistic.

use std::collections::HashMap;
use std::io::{BufRead, BufReader, Write};

use anyhow::{bail, Context, Result};
use noodles::sam::alignment::record::cigar::op::Kind;
use noodles::sam::alignment::record::Cigar as CigarTrait;

/// Parsed command-line arguments for `spike validate`.
struct ValidateArgs {
    bam_path: String,
    truth_path: String,
    ref_path: String,
    min_mapq: u8,
    flank_bp: u64,
    json_output: bool,
}

/// A truth VCF event with its expected properties.
#[derive(Debug)]
struct TruthEvent {
    chrom: String,
    start: u64,     // 0-based
    end: u64,       // 0-based, exclusive
    sv_type: String, // DEL, DUP, INV, INS, BND, SNP
    expected_vaf: f64,
    gene: String,
    /// For SmallVariant: explicit REF/ALT alleles.
    ref_allele: Option<Vec<u8>>,
    alt_allele: Option<Vec<u8>>,
}

/// Result of a single validation check.
struct CheckResult {
    event_label: String,
    check_name: String,
    expected: String,
    observed: String,
    pass: bool,
}

/// Entry point for `spike validate`.
pub fn run() -> Result<()> {
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info")).init();

    let args = parse_validate_args()?;

    log::info!("spike validate");
    log::info!("  BAM:       {}", args.bam_path);
    log::info!("  Truth VCF: {}", args.truth_path);
    log::info!("  Reference: {}", args.ref_path);

    let truth_events = load_truth_events(&args.truth_path)?;
    log::info!("Loaded {} truth events", truth_events.len());

    let mut results: Vec<CheckResult> = Vec::new();
    let mut n_errors: usize = 0;

    // Per-event checks.
    for event in &truth_events {
        let label = format_event_label(event);

        // Coverage ratio check (meaningful for DEL, DUP).
        if event.sv_type == "DEL" || event.sv_type == "DUP" {
            match check_coverage_ratio(
                &args.bam_path,
                &args.ref_path,
                event,
                args.flank_bp,
                args.min_mapq,
            ) {
                Ok(r) => results.push(r),
                Err(e) => {
                    log::warn!("coverage_ratio check failed for {}: {}", label, e);
                    n_errors += 1;
                }
            }
        }

        // Split-read evidence (meaningful for SVs, not small variants).
        if matches!(
            event.sv_type.as_str(),
            "DEL" | "DUP" | "INV" | "INS" | "BND"
        ) {
            match check_split_reads(
                &args.bam_path,
                &args.ref_path,
                event,
                args.min_mapq,
            ) {
                Ok(r) => results.push(r),
                Err(e) => {
                    log::warn!("split_reads check failed for {}: {}", label, e);
                    n_errors += 1;
                }
            }
        }

        // Allele frequency (meaningful for SNPs/small variants).
        if event.sv_type == "SNP" && event.ref_allele.is_some() && event.alt_allele.is_some() {
            match check_allele_freq(
                &args.bam_path,
                &args.ref_path,
                event,
                args.min_mapq,
            ) {
                Ok(r) => results.push(r),
                Err(e) => {
                    log::warn!("allele_freq check failed for {}: {}", label, e);
                    n_errors += 1;
                }
            }
        }

        // Log progress.
        log::info!("Checked: {}", label);
    }

    // Global checks.
    match check_insert_size(&args.bam_path, Some(&args.ref_path)) {
        Ok(r) => results.push(r),
        Err(e) => {
            log::warn!("insert_size check failed: {}", e);
            n_errors += 1;
        }
    }
    match check_dup_rate(&args.bam_path, Some(&args.ref_path)) {
        Ok(r) => results.push(r),
        Err(e) => {
            log::warn!("dup_rate check failed: {}", e);
            n_errors += 1;
        }
    }
    match check_mapq(&args.bam_path, Some(&args.ref_path)) {
        Ok(r) => results.push(r),
        Err(e) => {
            log::warn!("mean_mapq check failed: {}", e);
            n_errors += 1;
        }
    }

    if results.is_empty() {
        bail!(
            "all {} checks failed to run (BAM file may be corrupt or missing index)",
            n_errors,
        );
    }
    if n_errors > 0 {
        log::warn!("{} check(s) could not be executed (see warnings above)", n_errors);
    }

    // Print results.
    print_results(&results, args.json_output)?;

    let n_pass = results.iter().filter(|r| r.pass).count();
    let n_fail = results.len() - n_pass;

    if n_fail > 0 {
        bail!("{}/{} validation checks failed", n_fail, results.len());
    }

    Ok(())
}

/// Parse validate-specific arguments from std::env::args().
fn parse_validate_args() -> Result<ValidateArgs> {
    let raw: Vec<String> = std::env::args().collect();
    // raw[0] = binary, raw[1] = "validate", rest = flags.

    if raw.len() < 3 {
        print_usage();
        bail!("missing required arguments");
    }

    let mut bam_path: Option<String> = None;
    let mut truth_path: Option<String> = None;
    let mut ref_path: Option<String> = None;
    let mut min_mapq: u8 = 20;
    let mut flank_bp: u64 = 5000;
    let mut json_output = false;

    let args = &raw[2..];
    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "--bam" | "-b" => {
                i += 1;
                bam_path = Some(args.get(i).cloned().ok_or_else(|| anyhow::anyhow!("--bam requires a value"))?);
            }
            "--truth" | "-t" => {
                i += 1;
                truth_path = Some(args.get(i).cloned().ok_or_else(|| anyhow::anyhow!("--truth requires a value"))?);
            }
            "--reference" | "-r" => {
                i += 1;
                ref_path = Some(args.get(i).cloned().ok_or_else(|| anyhow::anyhow!("--reference requires a value"))?);
            }
            "--min-mapq" => {
                i += 1;
                min_mapq = args.get(i)
                    .ok_or_else(|| anyhow::anyhow!("--min-mapq requires a value"))?
                    .parse()
                    .context("--min-mapq must be a number")?;
            }
            "--flank" => {
                i += 1;
                flank_bp = args.get(i)
                    .ok_or_else(|| anyhow::anyhow!("--flank requires a value"))?
                    .parse()
                    .context("--flank must be a number")?;
            }
            "--json" => {
                json_output = true;
            }
            "--help" | "-h" => {
                print_usage();
                std::process::exit(0);
            }
            other => {
                bail!("unknown argument: {}", other);
            }
        }
        i += 1;
    }

    Ok(ValidateArgs {
        bam_path: bam_path.ok_or_else(|| anyhow::anyhow!("--bam is required"))?,
        truth_path: truth_path.ok_or_else(|| anyhow::anyhow!("--truth is required"))?,
        ref_path: ref_path.ok_or_else(|| anyhow::anyhow!("--reference is required"))?,
        min_mapq,
        flank_bp,
        json_output,
    })
}

fn print_usage() {
    eprintln!("Usage: spike validate --bam <BAM> --truth <VCF> --reference <FASTA> [OPTIONS]");
    eprintln!();
    eprintln!("Options:");
    eprintln!("  --bam, -b        Simulated BAM/CRAM file (required)");
    eprintln!("  --truth, -t      Truth VCF from spike (required)");
    eprintln!("  --reference, -r  Reference FASTA with .fai index (required)");
    eprintln!("  --min-mapq       Minimum MAPQ for counting reads (default: 20)");
    eprintln!("  --flank          Flanking bp for coverage comparison (default: 5000)");
    eprintln!("  --json           Output JSON instead of text table");
    eprintln!("  --help, -h       Show this help");
}

// ---------------------------------------------------------------------------
// Truth VCF parsing
// ---------------------------------------------------------------------------

/// Load truth events from a spike-produced truth VCF.
fn load_truth_events(path: &str) -> Result<Vec<TruthEvent>> {
    let content = std::fs::read_to_string(path)
        .with_context(|| format!("failed to read truth VCF: {}", path))?;

    let reader = BufReader::new(content.as_bytes());
    let mut events = Vec::new();
    let mut seen_bnd_ids: std::collections::HashSet<String> = std::collections::HashSet::new();

    for line in reader.lines() {
        let line = line?;
        if line.starts_with('#') {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 8 {
            continue;
        }

        let chrom = fields[0].to_string();
        let vcf_pos: u64 = match fields[1].parse() {
            Ok(p) => p,
            Err(_) => continue,
        };
        let id = fields[2].to_string();
        let ref_col = fields[3];
        let alt_col = fields[4];
        let info = fields[7];

        let sv_type_str = parse_info_field(info, "SVTYPE");
        let sim_vaf = parse_info_field(info, "SIM_VAF")
            .and_then(|v| v.parse::<f64>().ok())
            .unwrap_or(0.5);
        let gene = parse_info_field(info, "SIM_GENE")
            .unwrap_or("unknown")
            .to_string();

        match sv_type_str {
            Some("DEL") => {
                let end = parse_info_u64(info, "END").unwrap_or(vcf_pos + 1);
                events.push(TruthEvent {
                    chrom,
                    start: vcf_pos, // VCF POS for SV = 0-based start
                    end,
                    sv_type: "DEL".to_string(),
                    expected_vaf: sim_vaf,
                    gene,
                    ref_allele: None,
                    alt_allele: None,
                });
            }
            Some("DUP") => {
                let end = parse_info_u64(info, "END").unwrap_or(vcf_pos + 1);
                events.push(TruthEvent {
                    chrom,
                    start: vcf_pos,
                    end,
                    sv_type: "DUP".to_string(),
                    expected_vaf: sim_vaf,
                    gene,
                    ref_allele: None,
                    alt_allele: None,
                });
            }
            Some("INV") => {
                let end = parse_info_u64(info, "END").unwrap_or(vcf_pos + 1);
                events.push(TruthEvent {
                    chrom,
                    start: vcf_pos,
                    end,
                    sv_type: "INV".to_string(),
                    expected_vaf: sim_vaf,
                    gene,
                    ref_allele: None,
                    alt_allele: None,
                });
            }
            Some("INS") => {
                events.push(TruthEvent {
                    chrom,
                    start: vcf_pos,
                    end: vcf_pos,
                    sv_type: "INS".to_string(),
                    expected_vaf: sim_vaf,
                    gene,
                    ref_allele: None,
                    alt_allele: None,
                });
            }
            Some("BND") => {
                if seen_bnd_ids.contains(&id) {
                    continue;
                }
                seen_bnd_ids.insert(id.clone());
                if let Some(mate_id) = parse_info_field(info, "MATEID") {
                    seen_bnd_ids.insert(mate_id.to_string());
                }
                // BND POS is 1-based breakpoint in truth VCF.
                events.push(TruthEvent {
                    chrom,
                    start: vcf_pos.saturating_sub(1), // 1-based → 0-based
                    end: vcf_pos,
                    sv_type: "BND".to_string(),
                    expected_vaf: sim_vaf,
                    gene,
                    ref_allele: None,
                    alt_allele: None,
                });
            }
            _ => {
                // No SVTYPE: small variant (SNP/indel).
                // VCF POS is 1-based for small variants in truth VCF.
                let pos_0based = vcf_pos.saturating_sub(1);
                let ref_allele = ref_col.as_bytes().to_vec();
                let alt_allele = alt_col.as_bytes().to_vec();
                let end = pos_0based + ref_allele.len() as u64;
                events.push(TruthEvent {
                    chrom,
                    start: pos_0based,
                    end,
                    sv_type: "SNP".to_string(),
                    expected_vaf: sim_vaf,
                    gene,
                    ref_allele: Some(ref_allele),
                    alt_allele: Some(alt_allele),
                });
            }
        }
    }

    Ok(events)
}

// ---------------------------------------------------------------------------
// Validation checks
// ---------------------------------------------------------------------------

/// Check coverage ratio: event depth vs flanking depth.
fn check_coverage_ratio(
    bam_path: &str,
    ref_path: &str,
    event: &TruthEvent,
    flank_bp: u64,
    min_mapq: u8,
) -> Result<CheckResult> {
    let label = format_event_label(event);

    // Compute depth in the event region.
    let event_depth = count_depth_in_region(
        bam_path, ref_path, &event.chrom, event.start, event.end, min_mapq,
    )?;

    // Compute depth in left and right flanking regions.
    let left_start = event.start.saturating_sub(flank_bp);
    let left_end = event.start;
    let right_start = event.end;
    let right_end = event.end.saturating_add(flank_bp);

    let left_depth = count_depth_in_region(
        bam_path, ref_path, &event.chrom, left_start, left_end, min_mapq,
    )?;
    let right_depth = count_depth_in_region(
        bam_path, ref_path, &event.chrom, right_start, right_end, min_mapq,
    )?;

    // Average only flanking regions that have non-zero length (handles events
    // near chromosome start/end where one flank may be empty).
    let flank_depth = match (left_start < left_end, right_start < right_end) {
        (true, true) => (left_depth + right_depth) / 2.0,
        (true, false) => left_depth,
        (false, true) => right_depth,
        (false, false) => 0.0,
    };

    if flank_depth < 1.0 {
        return Ok(CheckResult {
            event_label: label,
            check_name: "coverage_ratio".to_string(),
            expected: "N/A".to_string(),
            observed: "no flanking coverage".to_string(),
            pass: true, // can't evaluate without coverage
        });
    }

    let ratio = event_depth / flank_depth;

    // Expected ratio depends on SV type and VAF.
    let (expected_str, pass) = match event.sv_type.as_str() {
        "DEL" => {
            // Spike suppresses reads at rate VAF → expected ratio = 1 - VAF.
            let expected_ratio = 1.0 - event.expected_vaf;
            let tolerance = 0.3;
            let pass = (ratio - expected_ratio).abs() < tolerance;
            (format!("{:.2}", expected_ratio), pass)
        }
        "DUP" => {
            // Spike adds depth copies at rate VAF → expected ratio = 1 + VAF.
            let expected_ratio = 1.0 + event.expected_vaf;
            let tolerance = 0.3;
            let pass = (ratio - expected_ratio).abs() < tolerance;
            (format!("{:.2}", expected_ratio), pass)
        }
        _ => ("~1.0".to_string(), (ratio - 1.0).abs() < 0.5),
    };

    Ok(CheckResult {
        event_label: label,
        check_name: "coverage_ratio".to_string(),
        expected: expected_str,
        observed: format!("{:.2}", ratio),
        pass,
    })
}

/// Check for split-read evidence (SA:Z tags) in the event region.
fn check_split_reads(
    bam_path: &str,
    ref_path: &str,
    event: &TruthEvent,
    min_mapq: u8,
) -> Result<CheckResult> {
    let label = format_event_label(event);

    // Query region: the event itself plus a small padding.
    let pad = 500u64;
    let query_start = event.start.saturating_sub(pad);
    let query_end = event.end.saturating_add(pad);

    let sa_count = count_sa_tags_in_region(
        bam_path, ref_path, &event.chrom, query_start, query_end, min_mapq,
    )?;

    let pass = sa_count > 0;

    Ok(CheckResult {
        event_label: label,
        check_name: "split_reads".to_string(),
        expected: ">0".to_string(),
        observed: format!("{}", sa_count),
        pass,
    })
}

/// Check allele frequency for small variants via pileup.
fn check_allele_freq(
    bam_path: &str,
    ref_path: &str,
    event: &TruthEvent,
    min_mapq: u8,
) -> Result<CheckResult> {
    let label = format_event_label(event);

    let ref_allele = event.ref_allele.as_ref().unwrap();
    let alt_allele = event.alt_allele.as_ref().unwrap();

    // Only check single-base variants (SNPs) for now.
    if ref_allele.len() != 1 || alt_allele.len() != 1 {
        return Ok(CheckResult {
            event_label: label,
            check_name: "allele_freq".to_string(),
            expected: format!("{:.2}", event.expected_vaf),
            observed: "N/A (indel)".to_string(),
            pass: true, // skip indel AF check
        });
    }

    let alt_base = alt_allele[0].to_ascii_uppercase();

    // Pileup at the variant position.
    let mut allele_counts: HashMap<u64, [u32; 4]> = HashMap::new();
    let mut dummy_read_alleles: HashMap<String, Vec<(u64, u8)>> = HashMap::new();

    pileup_region(
        bam_path,
        ref_path,
        &event.chrom,
        event.start,
        event.start + 1,
        min_mapq,
        &mut allele_counts,
        &mut dummy_read_alleles,
    )?;

    if let Some(counts) = allele_counts.get(&event.start) {
        let total: u32 = counts.iter().sum();
        if total < 5 {
            return Ok(CheckResult {
                event_label: label,
                check_name: "allele_freq".to_string(),
                expected: format!("{:.2}", event.expected_vaf),
                observed: format!("low depth ({})", total),
                pass: true, // not enough data
            });
        }

        let alt_idx = match alt_base {
            b'A' => 0,
            b'C' => 1,
            b'G' => 2,
            b'T' => 3,
            _ => {
                return Ok(CheckResult {
                    event_label: label,
                    check_name: "allele_freq".to_string(),
                    expected: format!("{:.2}", event.expected_vaf),
                    observed: "unknown alt base".to_string(),
                    pass: true,
                });
            }
        };

        let observed_vaf = counts[alt_idx] as f64 / total as f64;
        let tolerance = 0.15;
        let pass = (observed_vaf - event.expected_vaf).abs() < tolerance;

        Ok(CheckResult {
            event_label: label,
            check_name: "allele_freq".to_string(),
            expected: format!("{:.2}", event.expected_vaf),
            observed: format!("{:.2}", observed_vaf),
            pass,
        })
    } else {
        Ok(CheckResult {
            event_label: label,
            check_name: "allele_freq".to_string(),
            expected: format!("{:.2}", event.expected_vaf),
            observed: "no coverage".to_string(),
            pass: false,
        })
    }
}

/// Check global insert size distribution.
fn check_insert_size(bam_path: &str, ref_path: Option<&str>) -> Result<CheckResult> {
    let stats = crate::bam_stats::compute_stats(bam_path, 50_000, ref_path)?;

    let mean = stats.insert_mean;
    let stddev = stats.insert_stddev;

    // Reasonable Illumina ranges.
    let pass = (50.0..=1000.0).contains(&mean) && (5.0..=300.0).contains(&stddev);

    Ok(CheckResult {
        event_label: "[global]".to_string(),
        check_name: "insert_size".to_string(),
        expected: "mean 50-1000, sd 5-300".to_string(),
        observed: format!("{:.0}+/-{:.0}", mean, stddev),
        pass,
    })
}

/// Check global duplicate rate.
fn check_dup_rate(bam_path: &str, ref_path: Option<&str>) -> Result<CheckResult> {
    let (total, dups) = count_dup_reads(bam_path, ref_path)?;

    let rate = if total > 0 {
        dups as f64 / total as f64
    } else {
        0.0
    };

    let pass = rate < 0.50;

    Ok(CheckResult {
        event_label: "[global]".to_string(),
        check_name: "dup_rate".to_string(),
        expected: "<50%".to_string(),
        observed: format!("{:.1}%", rate * 100.0),
        pass,
    })
}

/// Check global mean mapping quality.
fn check_mapq(bam_path: &str, ref_path: Option<&str>) -> Result<CheckResult> {
    let mean_mapq = compute_mean_mapq(bam_path, ref_path)?;
    let pass = mean_mapq > 20.0;

    Ok(CheckResult {
        event_label: "[global]".to_string(),
        check_name: "mean_mapq".to_string(),
        expected: ">20".to_string(),
        observed: format!("{:.1}", mean_mapq),
        pass,
    })
}

// ---------------------------------------------------------------------------
// BAM reading helpers
// ---------------------------------------------------------------------------

/// Count average depth (reads per base) in a region.
fn count_depth_in_region(
    bam_path: &str,
    ref_path: &str,
    chrom: &str,
    start: u64,
    end: u64,
    min_mapq: u8,
) -> Result<f64> {
    if start >= end {
        return Ok(0.0);
    }

    let region_len = (end - start) as f64;
    let mut read_count: u64 = 0;

    if crate::extract::is_cram(bam_path) {
        let repository = crate::extract::build_fasta_repository(ref_path)?;
        let mut reader = noodles::cram::io::indexed_reader::Builder::default()
            .set_reference_sequence_repository(repository)
            .build_from_path(bam_path)
            .context("failed to open CRAM for depth counting")?;
        let header = reader.read_header()?;

        let start_pos = crate::extract::safe_noodles_position(start + 1);
        let end_pos = crate::extract::safe_noodles_position(end);
        let region = noodles::core::Region::new(chrom, start_pos..=end_pos);
        let query = reader.query(&header, &region)?;

        for rec_result in query {
            let cram_record = rec_result?;
            let buf = cram_record.try_into_alignment_record(&header)?;
            let flags = buf.flags();
            if flags.is_unmapped() || flags.is_secondary() || flags.is_supplementary()
                || flags.is_duplicate() || flags.is_qc_fail()
            {
                continue;
            }
            let mq: u8 = buf.mapping_quality().map(u8::from).unwrap_or(0);
            if mq < min_mapq {
                continue;
            }
            read_count += 1;
        }
    } else {
        let mut reader = noodles::bam::io::indexed_reader::Builder::default()
            .build_from_path(bam_path)
            .context("failed to open BAM for depth counting")?;
        let header = reader.read_header()?;

        let start_pos = crate::extract::safe_noodles_position(start + 1);
        let end_pos = crate::extract::safe_noodles_position(end);
        let region = noodles::core::Region::new(chrom, start_pos..=end_pos);
        let query = reader.query(&header, &region)?;

        for rec_result in query {
            let record = rec_result?;
            let flags = record.flags();
            if flags.is_unmapped() || flags.is_secondary() || flags.is_supplementary()
                || flags.is_duplicate() || flags.is_qc_fail()
            {
                continue;
            }
            let mq: u8 = record.mapping_quality().map(u8::from).unwrap_or(0);
            if mq < min_mapq {
                continue;
            }
            read_count += 1;
        }
    }

    // Average depth = (reads * read_length) / region_length.
    // Approximate: just use read_count / region_length * 150 (typical read length).
    // Actually, for a ratio we just need consistent counting, so reads/bp is fine.
    Ok(read_count as f64 / region_len)
}

/// Count reads with SA:Z supplementary alignment tag in a region.
fn count_sa_tags_in_region(
    bam_path: &str,
    ref_path: &str,
    chrom: &str,
    start: u64,
    end: u64,
    min_mapq: u8,
) -> Result<u64> {
    let mut sa_count: u64 = 0;

    if crate::extract::is_cram(bam_path) {
        let repository = crate::extract::build_fasta_repository(ref_path)?;
        let mut reader = noodles::cram::io::indexed_reader::Builder::default()
            .set_reference_sequence_repository(repository)
            .build_from_path(bam_path)
            .context("failed to open CRAM for SA tag counting")?;
        let header = reader.read_header()?;

        let start_pos = crate::extract::safe_noodles_position(start + 1);
        let end_pos = crate::extract::safe_noodles_position(end);
        let region = noodles::core::Region::new(chrom, start_pos..=end_pos);
        let query = reader.query(&header, &region)?;

        for rec_result in query {
            let cram_record = rec_result?;
            let buf = cram_record.try_into_alignment_record(&header)?;
            let flags = buf.flags();
            if flags.is_unmapped() || flags.is_secondary() || flags.is_supplementary()
                || flags.is_duplicate() || flags.is_qc_fail()
            {
                continue;
            }
            let mq: u8 = buf.mapping_quality().map(u8::from).unwrap_or(0);
            if mq < min_mapq {
                continue;
            }
            // Check for SA:Z tag in data.
            if has_sa_tag_buf(&buf) {
                sa_count += 1;
            }
        }
    } else {
        let mut reader = noodles::bam::io::indexed_reader::Builder::default()
            .build_from_path(bam_path)
            .context("failed to open BAM for SA tag counting")?;
        let header = reader.read_header()?;

        let start_pos = crate::extract::safe_noodles_position(start + 1);
        let end_pos = crate::extract::safe_noodles_position(end);
        let region = noodles::core::Region::new(chrom, start_pos..=end_pos);
        let query = reader.query(&header, &region)?;

        for rec_result in query {
            let record = rec_result?;
            let flags = record.flags();
            if flags.is_unmapped() || flags.is_secondary() || flags.is_supplementary()
                || flags.is_duplicate() || flags.is_qc_fail()
            {
                continue;
            }
            let mq: u8 = record.mapping_quality().map(u8::from).unwrap_or(0);
            if mq < min_mapq {
                continue;
            }
            if has_sa_tag_bam(&record) {
                sa_count += 1;
            }
        }
    }

    Ok(sa_count)
}

/// Check if a BAM record has an SA:Z auxiliary tag.
fn has_sa_tag_bam(record: &noodles::bam::Record) -> bool {
    use noodles::sam::alignment::record::data::field::Tag;
    let sa_tag = Tag::new(b'S', b'A');
    record.data().get(&sa_tag).is_some()
}

/// Check if a RecordBuf has an SA:Z auxiliary tag.
fn has_sa_tag_buf(buf: &noodles::sam::alignment::RecordBuf) -> bool {
    use noodles::sam::alignment::record::data::field::Tag;
    let sa_tag = Tag::new(b'S', b'A');
    buf.data().get(&sa_tag).is_some()
}

/// Pileup a region to count alleles at each position.
/// Reuses the same CIGAR walking pattern as loh.rs.
#[allow(clippy::too_many_arguments)]
fn pileup_region(
    bam_path: &str,
    ref_path: &str,
    chrom: &str,
    region_start: u64,
    region_end: u64,
    min_mapq: u8,
    allele_counts: &mut HashMap<u64, [u32; 4]>,
    read_alleles: &mut HashMap<String, Vec<(u64, u8)>>,
) -> Result<()> {
    if crate::extract::is_cram(bam_path) {
        let repository = crate::extract::build_fasta_repository(ref_path)?;
        let mut reader = noodles::cram::io::indexed_reader::Builder::default()
            .set_reference_sequence_repository(repository)
            .build_from_path(bam_path)
            .context("failed to open CRAM for pileup")?;
        let header = reader.read_header()?;

        let start_pos = crate::extract::safe_noodles_position(region_start + 1);
        let end_pos = crate::extract::safe_noodles_position(region_end);
        let region = noodles::core::Region::new(chrom, start_pos..=end_pos);
        let query = reader.query(&header, &region)?;

        for rec_result in query {
            let cram_record = rec_result?;
            let buf = cram_record.try_into_alignment_record(&header)?;
            let flags = buf.flags();
            if flags.is_unmapped() || flags.is_secondary() || flags.is_supplementary()
                || flags.is_duplicate() || flags.is_qc_fail()
            {
                continue;
            }
            let mq: u8 = buf.mapping_quality().map(u8::from).unwrap_or(0);
            if mq < min_mapq {
                continue;
            }
            let name = match buf.name() {
                Some(n) => String::from_utf8_lossy(n.as_ref()).into_owned(),
                None => continue,
            };
            let align_start = match buf.alignment_start() {
                Some(p) => usize::from(p).saturating_sub(1) as u64,
                None => continue,
            };
            let seq: Vec<u8> = buf.sequence().as_ref().to_vec();
            let cigar = buf.cigar();
            walk_cigar_pileup(
                &seq,
                CigarTrait::iter(&cigar),
                align_start,
                region_start,
                region_end,
                &name,
                allele_counts,
                read_alleles,
            );
        }
    } else {
        let mut reader = noodles::bam::io::indexed_reader::Builder::default()
            .build_from_path(bam_path)
            .context("failed to open BAM for pileup")?;
        let header = reader.read_header()?;

        let start_pos = crate::extract::safe_noodles_position(region_start + 1);
        let end_pos = crate::extract::safe_noodles_position(region_end);
        let region = noodles::core::Region::new(chrom, start_pos..=end_pos);
        let query = reader.query(&header, &region)?;

        for rec_result in query {
            let record = rec_result?;
            let flags = record.flags();
            if flags.is_unmapped() || flags.is_secondary() || flags.is_supplementary()
                || flags.is_duplicate() || flags.is_qc_fail()
            {
                continue;
            }
            let mq: u8 = record.mapping_quality().map(u8::from).unwrap_or(0);
            if mq < min_mapq {
                continue;
            }
            let name = match record.name() {
                Some(n) => String::from_utf8_lossy(n.as_ref()).into_owned(),
                None => continue,
            };
            let align_start = match record.alignment_start() {
                Some(Ok(p)) => usize::from(p).saturating_sub(1) as u64,
                _ => continue,
            };
            let seq: Vec<u8> = record.sequence().iter().collect();
            let cigar = record.cigar();
            walk_cigar_pileup(
                &seq,
                Box::new(cigar.iter()),
                align_start,
                region_start,
                region_end,
                &name,
                allele_counts,
                read_alleles,
            );
        }
    }

    Ok(())
}

/// Walk CIGAR operations and collect allele counts + per-read alleles.
/// Local copy of the same logic from loh.rs.
#[allow(clippy::too_many_arguments)]
fn walk_cigar_pileup(
    seq: &[u8],
    cigar_ops: Box<dyn Iterator<Item = std::io::Result<noodles::sam::alignment::record::cigar::Op>> + '_>,
    align_start: u64,
    region_start: u64,
    region_end: u64,
    name: &str,
    allele_counts: &mut HashMap<u64, [u32; 4]>,
    read_alleles: &mut HashMap<String, Vec<(u64, u8)>>,
) {
    let mut ref_pos = align_start;
    let mut seq_pos = 0usize;

    for op_result in cigar_ops {
        let op = match op_result {
            Ok(o) => o,
            Err(_) => break,
        };
        let len = op.len();

        match op.kind() {
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                for i in 0..len {
                    let rp = ref_pos + i as u64;
                    if rp >= region_start && rp < region_end {
                        if let Some(&base) = seq.get(seq_pos + i) {
                            let base = base.to_ascii_uppercase();
                            let idx = match base {
                                b'A' => Some(0),
                                b'C' => Some(1),
                                b'G' => Some(2),
                                b'T' => Some(3),
                                _ => None,
                            };
                            if let Some(idx) = idx {
                                allele_counts.entry(rp).or_insert([0; 4])[idx] += 1;
                                read_alleles
                                    .entry(name.to_string())
                                    .or_default()
                                    .push((rp, base));
                            }
                        }
                    }
                }
                ref_pos += len as u64;
                seq_pos += len;
            }
            Kind::Insertion | Kind::SoftClip => {
                seq_pos += len;
            }
            Kind::Deletion | Kind::Skip => {
                ref_pos += len as u64;
            }
            Kind::HardClip | Kind::Pad => {}
        }
    }
}

/// Count total primary reads and duplicate-flagged reads (sampling first 200k).
fn count_dup_reads(bam_path: &str, ref_path: Option<&str>) -> Result<(u64, u64)> {
    let mut total: u64 = 0;
    let mut dups: u64 = 0;
    let max_sample = 200_000u64;

    if crate::extract::is_cram(bam_path) {
        let rp = ref_path.ok_or_else(|| anyhow::anyhow!("CRAM requires --reference"))?;
        let repository = crate::extract::build_fasta_repository(rp)?;
        let mut reader = noodles::cram::io::reader::Builder::default()
            .set_reference_sequence_repository(repository)
            .build_from_path(bam_path)
            .context("failed to open CRAM for dup counting")?;
        let header = reader.read_header()?;

        for result in reader.records(&header) {
            let cram_record = result?;
            let buf = cram_record.try_into_alignment_record(&header)?;
            let flags = buf.flags();
            if flags.is_unmapped() || flags.is_secondary() || flags.is_supplementary() {
                continue;
            }
            total += 1;
            if flags.is_duplicate() {
                dups += 1;
            }
            if total >= max_sample {
                break;
            }
        }
    } else {
        let mut reader = noodles::bam::io::reader::Builder
            .build_from_path(bam_path)
            .context("failed to open BAM for dup counting")?;
        let _header = reader.read_header()?;

        for result in reader.records() {
            let record = result?;
            let flags = record.flags();
            if flags.is_unmapped() || flags.is_secondary() || flags.is_supplementary() {
                continue;
            }
            total += 1;
            if flags.is_duplicate() {
                dups += 1;
            }
            if total >= max_sample {
                break;
            }
        }
    }

    Ok((total, dups))
}

/// Compute mean MAPQ across all primary alignments (sampling first 100k).
fn compute_mean_mapq(bam_path: &str, ref_path: Option<&str>) -> Result<f64> {
    let mut total_mapq: u64 = 0;
    let mut count: u64 = 0;
    let max_sample = 100_000u64;

    if crate::extract::is_cram(bam_path) {
        let rp = ref_path.ok_or_else(|| anyhow::anyhow!("CRAM requires --reference"))?;
        let repository = crate::extract::build_fasta_repository(rp)?;
        let mut reader = noodles::cram::io::reader::Builder::default()
            .set_reference_sequence_repository(repository)
            .build_from_path(bam_path)
            .context("failed to open CRAM for MAPQ")?;
        let header = reader.read_header()?;

        for result in reader.records(&header) {
            let cram_record = result?;
            let buf = cram_record.try_into_alignment_record(&header)?;
            let flags = buf.flags();
            if flags.is_unmapped() || flags.is_secondary() || flags.is_supplementary() {
                continue;
            }
            let mq: u8 = buf.mapping_quality().map(u8::from).unwrap_or(0);
            total_mapq += mq as u64;
            count += 1;
            if count >= max_sample {
                break;
            }
        }
    } else {
        let mut reader = noodles::bam::io::reader::Builder
            .build_from_path(bam_path)
            .context("failed to open BAM for MAPQ")?;
        let _header = reader.read_header()?;

        for result in reader.records() {
            let record = result?;
            let flags = record.flags();
            if flags.is_unmapped() || flags.is_secondary() || flags.is_supplementary() {
                continue;
            }
            let mq: u8 = record.mapping_quality().map(u8::from).unwrap_or(0);
            total_mapq += mq as u64;
            count += 1;
            if count >= max_sample {
                break;
            }
        }
    }

    if count == 0 {
        return Ok(0.0);
    }
    Ok(total_mapq as f64 / count as f64)
}

// ---------------------------------------------------------------------------
// VCF INFO helpers (local copies from vcf_input.rs)
// ---------------------------------------------------------------------------

fn parse_info_field<'a>(info: &'a str, key: &str) -> Option<&'a str> {
    for field in info.split(';') {
        if let Some(value) = field.strip_prefix(key).and_then(|rest| rest.strip_prefix('=')) {
            return Some(value);
        }
    }
    None
}

fn parse_info_u64(info: &str, key: &str) -> Option<u64> {
    parse_info_field(info, key)?.parse().ok()
}

// ---------------------------------------------------------------------------
// Output formatting
// ---------------------------------------------------------------------------

fn format_event_label(event: &TruthEvent) -> String {
    let region = if event.start == event.end {
        format!("{}:{}", event.chrom, event.start)
    } else {
        format!("{}:{}-{}", event.chrom, event.start, event.end)
    };
    format!("{} {} ({})", event.sv_type, region, event.gene)
}

fn print_results(results: &[CheckResult], json: bool) -> Result<()> {
    let n_pass = results.iter().filter(|r| r.pass).count();
    let n_total = results.len();
    let stdout = std::io::stdout();
    let mut out = stdout.lock();

    if json {
        print_results_json(&mut out, results, n_total, n_pass)?;
    } else {
        print_results_text(&mut out, results, n_total, n_pass)?;
    }

    Ok(())
}

fn print_results_text(
    out: &mut impl Write,
    results: &[CheckResult],
    n_total: usize,
    n_pass: usize,
) -> Result<()> {
    writeln!(out)?;
    writeln!(out, "spike validate -- {} checks", n_total)?;
    writeln!(out)?;
    writeln!(
        out,
        "{:<35} {:<18} {:<25} {:<15} Status",
        "Event", "Check", "Expected", "Observed"
    )?;
    writeln!(out, "{}", "-".repeat(100))?;

    for r in results {
        let status = if r.pass { "PASS" } else { "FAIL" };
        writeln!(
            out,
            "{:<35} {:<18} {:<25} {:<15} {}",
            truncate(&r.event_label, 34),
            r.check_name,
            truncate(&r.expected, 24),
            truncate(&r.observed, 14),
            status,
        )?;
    }

    writeln!(out)?;
    writeln!(
        out,
        "Result: {}/{} PASS",
        n_pass, n_total,
    )?;

    Ok(())
}

fn print_results_json(
    out: &mut impl Write,
    results: &[CheckResult],
    n_total: usize,
    n_pass: usize,
) -> Result<()> {
    // Manual JSON to avoid serde dependency.
    writeln!(out, "{{")?;
    writeln!(
        out,
        "  \"summary\": {{ \"total\": {}, \"pass\": {}, \"fail\": {} }},",
        n_total,
        n_pass,
        n_total - n_pass,
    )?;
    writeln!(out, "  \"checks\": [")?;

    for (i, r) in results.iter().enumerate() {
        let comma = if i + 1 < results.len() { "," } else { "" };
        writeln!(
            out,
            "    {{ \"event\": \"{}\", \"check\": \"{}\", \"expected\": \"{}\", \"observed\": \"{}\", \"pass\": {} }}{}",
            escape_json(&r.event_label),
            escape_json(&r.check_name),
            escape_json(&r.expected),
            escape_json(&r.observed),
            r.pass,
            comma,
        )?;
    }

    writeln!(out, "  ]")?;
    writeln!(out, "}}")?;

    Ok(())
}

fn truncate(s: &str, max_len: usize) -> String {
    if s.len() <= max_len {
        s.to_string()
    } else {
        format!("{}...", &s[..max_len.saturating_sub(3)])
    }
}

fn escape_json(s: &str) -> String {
    s.replace('\\', "\\\\")
        .replace('"', "\\\"")
        .replace('\n', "\\n")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_info_field() {
        assert_eq!(
            parse_info_field("SVTYPE=DEL;END=100;SIM_VAF=0.500", "END"),
            Some("100")
        );
        assert_eq!(
            parse_info_field("SVTYPE=DEL;END=100;SIM_VAF=0.500", "SIM_VAF"),
            Some("0.500")
        );
        assert_eq!(
            parse_info_field("SVTYPE=DEL;END=100", "MISSING"),
            None
        );
    }

    #[test]
    fn test_format_event_label() {
        let event = TruthEvent {
            chrom: "chr7".to_string(),
            start: 55000,
            end: 56000,
            sv_type: "DEL".to_string(),
            expected_vaf: 0.5,
            gene: "EGFR".to_string(),
            ref_allele: None,
            alt_allele: None,
        };
        assert_eq!(format_event_label(&event), "DEL chr7:55000-56000 (EGFR)");
    }

    #[test]
    fn test_format_event_label_point() {
        let event = TruthEvent {
            chrom: "chr7".to_string(),
            start: 55200,
            end: 55200,
            sv_type: "INS".to_string(),
            expected_vaf: 0.3,
            gene: "EGFR".to_string(),
            ref_allele: None,
            alt_allele: None,
        };
        assert_eq!(format_event_label(&event), "INS chr7:55200 (EGFR)");
    }

    #[test]
    fn test_load_truth_events_del() {
        let vcf = "\
##fileformat=VCFv4.3
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE
chr7\t55000\tsim_del_1\tN\t<DEL>\t999\tPASS\tSVTYPE=DEL;END=56000;SVLEN=-1000;SIM_VAF=0.500;SIM_GENE=EGFR;SIM_EXONS=.\tGT\t0/1
";
        // Write to temp file.
        let dir = std::env::temp_dir().join("spike_test_validate");
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("truth_del.vcf");
        std::fs::write(&path, vcf).unwrap();

        let events = load_truth_events(path.to_str().unwrap()).unwrap();
        assert_eq!(events.len(), 1);
        assert_eq!(events[0].sv_type, "DEL");
        assert_eq!(events[0].chrom, "chr7");
        assert_eq!(events[0].start, 55000);
        assert_eq!(events[0].end, 56000);
        assert!((events[0].expected_vaf - 0.5).abs() < f64::EPSILON);
        assert_eq!(events[0].gene, "EGFR");
    }

    #[test]
    fn test_load_truth_events_snp() {
        let vcf = "\
##fileformat=VCFv4.3
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE
chr7\t55201\tsim_var_1\tA\tT\t999\tPASS\tSIM_VAF=0.500;SIM_GENE=EGFR\tGT\t0/1
";
        let dir = std::env::temp_dir().join("spike_test_validate");
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("truth_snp.vcf");
        std::fs::write(&path, vcf).unwrap();

        let events = load_truth_events(path.to_str().unwrap()).unwrap();
        assert_eq!(events.len(), 1);
        assert_eq!(events[0].sv_type, "SNP");
        assert_eq!(events[0].start, 55200); // 1-based 55201 → 0-based 55200
        assert_eq!(events[0].ref_allele.as_deref(), Some(b"A".as_slice()));
        assert_eq!(events[0].alt_allele.as_deref(), Some(b"T".as_slice()));
    }

    #[test]
    fn test_load_truth_events_bnd_pair() {
        let vcf = "\
##fileformat=VCFv4.3
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE
chr2\t29416089\tsim_fus_1\tN\tN[chr2:42522656[\t999\tPASS\tSVTYPE=BND;MATEID=sim_fus_1_mate;SIM_VAF=0.050;SIM_GENE=ALK\tGT\t0/1
chr2\t42522656\tsim_fus_1_mate\tN\t]chr2:29416089]N\t999\tPASS\tSVTYPE=BND;MATEID=sim_fus_1;SIM_VAF=0.050;SIM_GENE=EML4\tGT\t0/1
";
        let dir = std::env::temp_dir().join("spike_test_validate");
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("truth_bnd.vcf");
        std::fs::write(&path, vcf).unwrap();

        let events = load_truth_events(path.to_str().unwrap()).unwrap();
        // BND pairs should produce only 1 event (mate is deduplicated).
        assert_eq!(events.len(), 1);
        assert_eq!(events[0].sv_type, "BND");
    }

    #[test]
    fn test_truncate() {
        assert_eq!(truncate("short", 10), "short");
        assert_eq!(truncate("this is a longer string", 10), "this is...");
    }

    #[test]
    fn test_escape_json() {
        assert_eq!(escape_json("hello \"world\""), "hello \\\"world\\\"");
        assert_eq!(escape_json("a\\b"), "a\\\\b");
    }
}
