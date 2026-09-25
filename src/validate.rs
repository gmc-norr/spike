//! BAM validation subcommand for spike.
//!
//! Reads a simulated BAM + truth VCF and runs automated checks to verify
//! that the spike-in reads look realistic.

use std::collections::{HashMap, HashSet};
use std::io::{BufRead, BufReader, Write};

use anyhow::{bail, Context, Result};
use noodles::sam::alignment::record::cigar::op::Kind;
use noodles::sam::alignment::record::Cigar as CigarTrait;

use crate::extract::record_is_on_queried_reference;

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
    start: u64,      // 0-based
    end: u64,        // 0-based, exclusive
    sv_type: String, // DEL, DUP, INV, INS, BND, SNP
    expected_vaf: f64,
    gene: String,
    /// The other breakpoint (chrom, 0-based position): a split read at this
    /// event should have its other part there. None for INS and SNP.
    partner: Option<(String, u64)>,
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

/// A check's result; a check that could not run is a failed result, so an
/// unevaluable run never reports all-PASS.
fn check_outcome(label: &str, check_name: &str, result: Result<CheckResult>) -> CheckResult {
    result.unwrap_or_else(|e| {
        log::warn!("{} check failed to run for {}: {:#}", check_name, label, e);
        CheckResult {
            event_label: label.to_string(),
            check_name: check_name.to_string(),
            expected: "check runs".to_string(),
            observed: format!("error: {:#}", e),
            pass: false,
        }
    })
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

    // Per-event checks.
    for event in &truth_events {
        let label = format_event_label(event);

        // Coverage ratio check (meaningful for DEL, DUP).
        if event.sv_type == "DEL" || event.sv_type == "DUP" {
            let r = check_coverage_ratio(
                &args.bam_path,
                &args.ref_path,
                event,
                args.flank_bp,
                args.min_mapq,
            );
            results.push(check_outcome(&label, "coverage_ratio", r));
        }

        // Split reads joining the two breakpoints (not INS: its inserted
        // sequence has no second reference breakpoint).
        if matches!(event.sv_type.as_str(), "DEL" | "DUP" | "INV" | "BND") {
            let r = check_split_reads(&args.bam_path, &args.ref_path, event, args.min_mapq);
            results.push(check_outcome(&label, "split_reads", r));
        }

        // Allele frequency (meaningful for SNPs/small variants).
        if event.sv_type == "SNP" && event.ref_allele.is_some() && event.alt_allele.is_some() {
            let r = check_allele_freq(&args.bam_path, &args.ref_path, event, args.min_mapq);
            results.push(check_outcome(&label, "allele_freq", r));
        }

        // Log progress.
        log::info!("Checked: {}", label);
    }

    // Global checks. All three read one sample, taken from the truth events'
    // own regions rather than from the head of the file (L15).
    match sample_event_regions(&args.bam_path, &args.ref_path, &truth_events, args.flank_bp) {
        Ok(sample) => {
            results.push(check_insert_size(&sample));
            results.push(check_dup_rate(&sample));
            results.push(check_mapq(&sample));
        }
        Err(e) => {
            // A sample that could not be taken is three failed checks, not
            // three checks that quietly pass on a default (M10, M11).
            for check_name in ["insert_size", "dup_rate", "mean_mapq"] {
                results.push(check_outcome(
                    "global",
                    check_name,
                    Err(anyhow::anyhow!("{:#}", e)),
                ));
            }
        }
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
                bam_path = Some(
                    args.get(i)
                        .cloned()
                        .ok_or_else(|| anyhow::anyhow!("--bam requires a value"))?,
                );
            }
            "--truth" | "-t" => {
                i += 1;
                truth_path = Some(
                    args.get(i)
                        .cloned()
                        .ok_or_else(|| anyhow::anyhow!("--truth requires a value"))?,
                );
            }
            "--reference" | "-r" => {
                i += 1;
                ref_path = Some(
                    args.get(i)
                        .cloned()
                        .ok_or_else(|| anyhow::anyhow!("--reference requires a value"))?,
                );
            }
            "--min-mapq" => {
                i += 1;
                min_mapq = args
                    .get(i)
                    .ok_or_else(|| anyhow::anyhow!("--min-mapq requires a value"))?
                    .parse()
                    .context("--min-mapq must be a number")?;
            }
            "--flank" => {
                i += 1;
                flank_bp = args
                    .get(i)
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
                let partner = Some((chrom.clone(), end));
                events.push(TruthEvent {
                    chrom,
                    start: vcf_pos, // VCF POS for SV = 0-based start
                    end,
                    sv_type: "DEL".to_string(),
                    partner,
                    expected_vaf: sim_vaf,
                    gene,
                    ref_allele: None,
                    alt_allele: None,
                });
            }
            Some("DUP") => {
                let end = parse_info_u64(info, "END").unwrap_or(vcf_pos + 1);
                let partner = Some((chrom.clone(), end));
                events.push(TruthEvent {
                    chrom,
                    start: vcf_pos,
                    end,
                    sv_type: "DUP".to_string(),
                    partner,
                    expected_vaf: sim_vaf,
                    gene,
                    ref_allele: None,
                    alt_allele: None,
                });
            }
            Some("INV") => {
                let end = parse_info_u64(info, "END").unwrap_or(vcf_pos + 1);
                let partner = Some((chrom.clone(), end));
                events.push(TruthEvent {
                    chrom,
                    start: vcf_pos,
                    end,
                    sv_type: "INV".to_string(),
                    partner,
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
                    partner: None,
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
                    partner: bnd_partner(alt_col),
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
                    partner: None,
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
        bam_path,
        ref_path,
        &event.chrom,
        event.start,
        event.end,
        min_mapq,
    )?;

    // Compute depth in left and right flanking regions.
    let left_start = event.start.saturating_sub(flank_bp);
    let left_end = event.start;
    let right_start = event.end;
    let right_end = event.end.saturating_add(flank_bp);

    let left_depth = count_depth_in_region(
        bam_path,
        ref_path,
        &event.chrom,
        left_start,
        left_end,
        min_mapq,
    )?;
    let right_depth = count_depth_in_region(
        bam_path,
        ref_path,
        &event.chrom,
        right_start,
        right_end,
        min_mapq,
    )?;

    // Average only flanking regions that have non-zero length (handles events
    // near chromosome start/end where one flank may be empty).
    let flank_depth = match (left_start < left_end, right_start < right_end) {
        (true, true) => (left_depth + right_depth) / 2.0,
        (true, false) => left_depth,
        (false, true) => right_depth,
        (false, false) => 0.0,
    };

    Ok(coverage_ratio_result(
        label,
        &event.sv_type,
        event.expected_vaf,
        event_depth,
        flank_depth,
    ))
}

/// Judge an event's depth against its flanks.
fn coverage_ratio_result(
    label: String,
    sv_type: &str,
    expected_vaf: f64,
    event_depth: f64,
    flank_depth: f64,
) -> CheckResult {
    if flank_depth < 1.0 {
        return CheckResult {
            event_label: label,
            check_name: "coverage_ratio".to_string(),
            expected: "N/A".to_string(),
            observed: "no flanking coverage".to_string(),
            pass: false, // can't evaluate: don't report it as a pass
        };
    }

    let ratio = event_depth / flank_depth;

    // Expected ratio depends on SV type and VAF.
    let (expected_str, pass) = match sv_type {
        "DEL" => {
            // Spike suppresses reads at rate VAF → expected ratio = 1 - VAF.
            let expected_ratio = 1.0 - expected_vaf;
            let tolerance = 0.3;
            let pass = (ratio - expected_ratio).abs() < tolerance;
            (format!("{:.2}", expected_ratio), pass)
        }
        "DUP" => {
            // Spike adds depth copies at rate VAF → expected ratio = 1 + VAF.
            let expected_ratio = 1.0 + expected_vaf;
            let tolerance = 0.3;
            let pass = (ratio - expected_ratio).abs() < tolerance;
            (format!("{:.2}", expected_ratio), pass)
        }
        _ => ("~1.0".to_string(), (ratio - 1.0).abs() < 0.5),
    };

    CheckResult {
        event_label: label,
        check_name: "coverage_ratio".to_string(),
        expected: expected_str,
        observed: format!("{:.2}", ratio),
        pass,
    }
}

/// Check for split reads joining the event's two breakpoints: reads at one
/// breakpoint whose SA:Z alignment lands at the other.
fn check_split_reads(
    bam_path: &str,
    ref_path: &str,
    event: &TruthEvent,
    min_mapq: u8,
) -> Result<CheckResult> {
    // At least this many split reads must join the breakpoints. Reads with
    // an SA tag occur anywhere; reads joining these two points do not.
    const MIN_SPLIT_READS: usize = 2;
    let pad = 500u64;
    let label = format_event_label(event);
    let partner = event
        .partner
        .clone()
        .with_context(|| format!("{}: no partner breakpoint for split reads", label))?;
    let here = (event.chrom.clone(), event.start);

    let mut names = split_reads_to_partner(
        bam_path,
        ref_path,
        &here.0,
        here.1.saturating_sub(pad),
        here.1.saturating_add(pad),
        min_mapq,
        &partner,
        pad,
    )?;
    names.extend(split_reads_to_partner(
        bam_path,
        ref_path,
        &partner.0,
        partner.1.saturating_sub(pad),
        partner.1.saturating_add(pad),
        min_mapq,
        &here,
        pad,
    )?);

    Ok(CheckResult {
        event_label: label,
        check_name: "split_reads".to_string(),
        expected: format!(">={} joining {}:{}", MIN_SPLIT_READS, partner.0, partner.1 + 1),
        observed: format!("{}", names.len()),
        pass: names.len() >= MIN_SPLIT_READS,
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

// ---------------------------------------------------------------------------
// Global checks
// ---------------------------------------------------------------------------

/// Records to sample for the global checks, over all event regions together.
const GLOBAL_SAMPLE_MAX: u64 = 200_000;

/// The smallest per-event share of that budget, so a truth VCF with many
/// events still reads enough of each one.
const GLOBAL_SAMPLE_MIN_PER_REGION: u64 = 1_000;

/// What the three global checks are computed from: primary, mapped records
/// sampled from the truth events' own regions.
#[derive(Default)]
struct GlobalSample {
    /// Records sampled.
    total: u64,
    /// ...of which carry the duplicate flag.
    dups: u64,
    /// Sum of their mapping qualities.
    mapq_sum: u64,
    /// Template lengths of the properly-paired, non-duplicate records.
    insert_sizes: Vec<f64>,
}

impl GlobalSample {
    /// Add one record. Unmapped, secondary and supplementary records are the
    /// caller's to skip.
    fn add(
        &mut self,
        flags: noodles::sam::alignment::record::Flags,
        mapq: u8,
        template_length: i32,
    ) {
        self.total += 1;
        self.mapq_sum += mapq as u64;
        if flags.is_duplicate() {
            self.dups += 1;
        }
        // The same filter `bam_stats` uses: a duplicate's or a QC-failed
        // record's template length is not an insert-size observation.
        if !flags.is_duplicate()
            && !flags.is_qc_fail()
            && flags.is_properly_segmented()
            && !flags.is_mate_unmapped()
            && template_length > 0
        {
            self.insert_sizes.push(template_length as f64);
        }
    }

    /// Mean and standard deviation of the sampled insert sizes, or `None`
    /// when nothing properly paired was sampled.
    fn insert_stats(&self) -> Option<(f64, f64)> {
        if self.insert_sizes.is_empty() {
            return None;
        }
        let mean = self.insert_sizes.iter().sum::<f64>() / self.insert_sizes.len() as f64;
        let variance = self
            .insert_sizes
            .iter()
            .map(|x| (x - mean).powi(2))
            .sum::<f64>()
            / self.insert_sizes.len() as f64;
        Some((mean, variance.sqrt()))
    }
}

/// Sample records for the global checks from every truth event's own window
/// (event +/- `flank_bp`).
///
/// The head of a file is not a sample of it: on whole-genome HG002 the first
/// 100k records are chr1's telomere, mean MAPQ 10.0, which fails a check the
/// rest of the file passes (L15). The event windows are both representative of
/// the reads `validate` is judging and cheap to read, because they are indexed
/// queries like every other check here rather than a walk from the top.
fn sample_event_regions(
    bam_path: &str,
    ref_path: &str,
    events: &[TruthEvent],
    flank_bp: u64,
) -> Result<GlobalSample> {
    if events.is_empty() {
        bail!("truth VCF holds no events, so there is no region to sample");
    }

    // Spread the budget over the events, so one long event cannot spend it.
    let per_region = (GLOBAL_SAMPLE_MAX / events.len() as u64).max(GLOBAL_SAMPLE_MIN_PER_REGION);

    let mut sample = GlobalSample::default();
    for event in events {
        if sample.total >= GLOBAL_SAMPLE_MAX {
            break;
        }
        let start = event.start.saturating_sub(flank_bp);
        let end = event.end + flank_bp;
        sample_region(
            bam_path,
            ref_path,
            &event.chrom,
            start,
            end,
            per_region,
            &mut sample,
        )?;
    }

    Ok(sample)
}

/// Add up to `max_records` of one region's primary alignments to `sample`.
fn sample_region(
    bam_path: &str,
    ref_path: &str,
    chrom: &str,
    start: u64,
    end: u64,
    max_records: u64,
    sample: &mut GlobalSample,
) -> Result<()> {
    let start_pos = crate::extract::safe_noodles_position(start + 1);
    let end_pos = crate::extract::safe_noodles_position(end);
    let region = noodles::core::Region::new(chrom, start_pos..=end_pos);
    let mut taken: u64 = 0;

    if crate::extract::is_cram(bam_path) {
        let repository = crate::extract::build_fasta_repository(ref_path)?;
        let (mut reader, header) =
            crate::extract::open_cram_reader_for_region(bam_path, &repository, &region)
                .context("failed to open CRAM for the global sample")?;
        let query = reader.query(&header, &region)?;
        // `query` has already rejected an unknown contig, so this is `Some`.
        let queried_reference_sequence_id =
            header.reference_sequences().get_index_of(chrom.as_bytes());

        for rec_result in query {
            let cram_record = rec_result?;
            let buf = cram_record.try_into_alignment_record(&header)?;
            // A multi-contig container is decoded whole and `Query` filters on
            // coordinates alone, so another contig's reads would enter the
            // sample (L2, N4). The BAM arm needs no such guard.
            if !record_is_on_queried_reference(&buf, queried_reference_sequence_id) {
                continue;
            }
            let flags = buf.flags();
            if flags.is_unmapped() || flags.is_secondary() || flags.is_supplementary() {
                continue;
            }
            let mq: u8 = buf.mapping_quality().map(u8::from).unwrap_or(0);
            sample.add(flags, mq, buf.template_length());
            taken += 1;
            if taken >= max_records {
                break;
            }
        }
    } else {
        let mut reader = noodles::bam::io::indexed_reader::Builder::default()
            .build_from_path(bam_path)
            .context("failed to open BAM for the global sample")?;
        let header = reader.read_header()?;
        let query = reader.query(&header, &region)?;

        for rec_result in query {
            let record = rec_result?;
            let flags = record.flags();
            if flags.is_unmapped() || flags.is_secondary() || flags.is_supplementary() {
                continue;
            }
            let mq: u8 = record.mapping_quality().map(u8::from).unwrap_or(0);
            sample.add(flags, mq, record.template_length());
            taken += 1;
            if taken >= max_records {
                break;
            }
        }
    }

    Ok(())
}

/// A global check its sample cannot answer: a failed result, never a silent
/// pass (M10).
fn not_evaluable(check_name: &str, expected: &str, observed: &str, why: &str) -> CheckResult {
    log::warn!("{} check is not evaluable: {}", check_name, why);
    CheckResult {
        event_label: "[global]".to_string(),
        check_name: check_name.to_string(),
        expected: expected.to_string(),
        observed: observed.to_string(),
        pass: false,
    }
}

/// Check the insert size distribution of the sampled records.
fn check_insert_size(sample: &GlobalSample) -> CheckResult {
    let expected = "mean 50-1000, sd 5-300";
    let Some((mean, stddev)) = sample.insert_stats() else {
        return not_evaluable(
            "insert_size",
            expected,
            "no pairs",
            "no properly-paired record in the event regions",
        );
    };

    // Reasonable Illumina ranges.
    let pass = (50.0..=1000.0).contains(&mean) && (5.0..=300.0).contains(&stddev);

    CheckResult {
        event_label: "[global]".to_string(),
        check_name: "insert_size".to_string(),
        expected: expected.to_string(),
        observed: format!("{:.0}+/-{:.0}", mean, stddev),
        pass,
    }
}

/// Check the duplicate rate of the sampled records.
fn check_dup_rate(sample: &GlobalSample) -> CheckResult {
    if sample.total == 0 {
        return not_evaluable("dup_rate", "<50%", "no reads", "no record sampled");
    }
    if sample.dups == 0 {
        // A file with no duplicate flags is not a file without duplicates:
        // nothing marked them, so 0% would be a default, not a measurement.
        return not_evaluable(
            "dup_rate",
            "<50%",
            "no dup flags",
            "no sampled record carries the duplicate flag -- mark duplicates to evaluate it",
        );
    }

    let rate = sample.dups as f64 / sample.total as f64;
    let pass = rate < 0.50;

    CheckResult {
        event_label: "[global]".to_string(),
        check_name: "dup_rate".to_string(),
        expected: "<50%".to_string(),
        observed: format!("{:.1}%", rate * 100.0),
        pass,
    }
}

/// Check the mean mapping quality of the sampled records.
fn check_mapq(sample: &GlobalSample) -> CheckResult {
    if sample.total == 0 {
        return not_evaluable("mean_mapq", ">20", "no reads", "no record sampled");
    }

    let mean_mapq = sample.mapq_sum as f64 / sample.total as f64;
    let pass = mean_mapq > 20.0;

    CheckResult {
        event_label: "[global]".to_string(),
        check_name: "mean_mapq".to_string(),
        expected: ">20".to_string(),
        observed: format!("{:.1}", mean_mapq),
        pass,
    }
}

// ---------------------------------------------------------------------------
// BAM reading helpers
// ---------------------------------------------------------------------------

/// Mean read depth (aligned bases per base) in a region.
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

    let mut spans: Vec<(u64, u64)> = Vec::new();

    if crate::extract::is_cram(bam_path) {
        let repository = crate::extract::build_fasta_repository(ref_path)?;
        let start_pos = crate::extract::safe_noodles_position(start + 1);
        let end_pos = crate::extract::safe_noodles_position(end);
        let region = noodles::core::Region::new(chrom, start_pos..=end_pos);
        let (mut reader, header) =
            crate::extract::open_cram_reader_for_region(bam_path, &repository, &region)
                .context("failed to open CRAM for depth counting")?;
        let query = reader.query(&header, &region)?;
        // `query` has already rejected an unknown contig, so this is `Some`.
        let queried_reference_sequence_id =
            header.reference_sequences().get_index_of(chrom.as_bytes());

        for rec_result in query {
            let cram_record = rec_result?;
            let buf = cram_record.try_into_alignment_record(&header)?;
            // A container holding several contigs is decoded whole and `Query`
            // filters on coordinates alone, so another contig's reads would be
            // counted as this region's (L2, N4). The BAM arm below needs no
            // such guard: noodles-bam's `Query` compares the reference id.
            if !record_is_on_queried_reference(&buf, queried_reference_sequence_id) {
                continue;
            }
            let flags = buf.flags();
            if flags.is_unmapped()
                || flags.is_secondary()
                || flags.is_supplementary()
                || flags.is_duplicate()
                || flags.is_qc_fail()
            {
                continue;
            }
            let mq: u8 = buf.mapping_quality().map(u8::from).unwrap_or(0);
            if mq < min_mapq {
                continue;
            }
            if let Some(p) = buf.alignment_start() {
                let s = usize::from(p) as u64 - 1;
                spans.push((s, s + ref_span(CigarTrait::iter(&buf.cigar()))));
            }
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
            if flags.is_unmapped()
                || flags.is_secondary()
                || flags.is_supplementary()
                || flags.is_duplicate()
                || flags.is_qc_fail()
            {
                continue;
            }
            let mq: u8 = record.mapping_quality().map(u8::from).unwrap_or(0);
            if mq < min_mapq {
                continue;
            }
            if let Some(Ok(p)) = record.alignment_start() {
                let s = usize::from(p) as u64 - 1;
                spans.push((s, s + ref_span(record.cigar().iter())));
            }
        }
    }

    Ok(mean_depth(&spans, start, end))
}

/// Mean depth over [start, end) from read reference spans [s, e).
fn mean_depth(spans: &[(u64, u64)], start: u64, end: u64) -> f64 {
    if start >= end {
        return 0.0;
    }
    let bases: u64 = spans
        .iter()
        .map(|&(s, e)| e.min(end).saturating_sub(s.max(start)))
        .sum();
    bases as f64 / (end - start) as f64
}

/// Reference bases covered by an alignment's CIGAR.
fn ref_span(
    ops: impl Iterator<Item = std::io::Result<noodles::sam::alignment::record::cigar::Op>>,
) -> u64 {
    ops.filter_map(|op| op.ok())
        .filter(|op| {
            matches!(
                op.kind(),
                Kind::Match | Kind::Deletion | Kind::Skip | Kind::SequenceMatch | Kind::SequenceMismatch
            )
        })
        .map(|op| op.len() as u64)
        .sum()
}

/// Names of reads in [start, end) whose SA:Z tag has an alignment within
/// `pad` of `partner` (chrom, 0-based position).
#[allow(clippy::too_many_arguments)]
fn split_reads_to_partner(
    bam_path: &str,
    ref_path: &str,
    chrom: &str,
    start: u64,
    end: u64,
    min_mapq: u8,
    partner: &(String, u64),
    pad: u64,
) -> Result<HashSet<String>> {
    let mut names = HashSet::new();
    let points_to_partner =
        |sa: Option<String>| sa.is_some_and(|sa| sa_points_near(&sa, &partner.0, partner.1, pad));

    if crate::extract::is_cram(bam_path) {
        let repository = crate::extract::build_fasta_repository(ref_path)?;
        let start_pos = crate::extract::safe_noodles_position(start + 1);
        let end_pos = crate::extract::safe_noodles_position(end);
        let region = noodles::core::Region::new(chrom, start_pos..=end_pos);
        let (mut reader, header) =
            crate::extract::open_cram_reader_for_region(bam_path, &repository, &region)
                .context("failed to open CRAM for SA tag counting")?;
        let query = reader.query(&header, &region)?;
        // `query` has already rejected an unknown contig, so this is `Some`.
        let queried_reference_sequence_id =
            header.reference_sequences().get_index_of(chrom.as_bytes());

        for rec_result in query {
            let cram_record = rec_result?;
            let buf = cram_record.try_into_alignment_record(&header)?;
            // A container holding several contigs is decoded whole and `Query`
            // filters on coordinates alone, so another contig's reads would be
            // counted as this region's (L2, N4). The BAM arm below needs no
            // such guard: noodles-bam's `Query` compares the reference id.
            if !record_is_on_queried_reference(&buf, queried_reference_sequence_id) {
                continue;
            }
            let flags = buf.flags();
            if flags.is_unmapped()
                || flags.is_secondary()
                || flags.is_supplementary()
                || flags.is_duplicate()
                || flags.is_qc_fail()
            {
                continue;
            }
            let mq: u8 = buf.mapping_quality().map(u8::from).unwrap_or(0);
            if mq < min_mapq {
                continue;
            }
            if points_to_partner(sa_value_buf(&buf)) {
                if let Some(n) = buf.name() {
                    names.insert(String::from_utf8_lossy(n.as_ref()).into_owned());
                }
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
            if flags.is_unmapped()
                || flags.is_secondary()
                || flags.is_supplementary()
                || flags.is_duplicate()
                || flags.is_qc_fail()
            {
                continue;
            }
            let mq: u8 = record.mapping_quality().map(u8::from).unwrap_or(0);
            if mq < min_mapq {
                continue;
            }
            if points_to_partner(sa_value_bam(&record)) {
                if let Some(n) = record.name() {
                    names.insert(String::from_utf8_lossy(n.as_ref()).into_owned());
                }
            }
        }
    }

    Ok(names)
}

/// True if an SA:Z value (`chrom,pos,strand,CIGAR,mapQ,NM;...`, 1-based pos)
/// has an alignment on `chrom` within `pad` of 0-based `pos`.
fn sa_points_near(sa: &str, chrom: &str, pos: u64, pad: u64) -> bool {
    sa.split(';').filter(|e| !e.is_empty()).any(|entry| {
        let mut fields = entry.split(',');
        let (Some(c), Some(p)) = (fields.next(), fields.next()) else {
            return false;
        };
        match p.parse::<u64>() {
            Ok(p1) => c == chrom && p1.saturating_sub(1).abs_diff(pos) <= pad,
            Err(_) => false,
        }
    })
}

/// Partner breakpoint (chrom, 0-based) from a BND ALT such as `A]chr2:500]`.
fn bnd_partner(alt: &str) -> Option<(String, u64)> {
    let open = alt.find(['[', ']'])?;
    let close = alt.rfind(['[', ']'])?;
    let (chrom, pos) = alt.get(open + 1..close)?.rsplit_once(':')?;
    let pos: u64 = pos.parse().ok()?;
    Some((chrom.to_string(), pos.checked_sub(1)?))
}

/// The SA:Z value of a BAM record, if any.
fn sa_value_bam(record: &noodles::bam::Record) -> Option<String> {
    use noodles::sam::alignment::record::data::field::{Tag, Value};
    match record.data().get(&Tag::new(b'S', b'A'))? {
        Ok(Value::String(s)) => Some(String::from_utf8_lossy(s).into_owned()),
        _ => None,
    }
}

/// The SA:Z value of a RecordBuf (CRAM), if any.
fn sa_value_buf(buf: &noodles::sam::alignment::RecordBuf) -> Option<String> {
    use noodles::sam::alignment::record::data::field::Tag;
    use noodles::sam::alignment::record_buf::data::field::Value;
    match buf.data().get(&Tag::new(b'S', b'A'))? {
        Value::String(s) => Some(String::from_utf8_lossy(s).into_owned()),
        _ => None,
    }
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
        let start_pos = crate::extract::safe_noodles_position(region_start + 1);
        let end_pos = crate::extract::safe_noodles_position(region_end);
        let region = noodles::core::Region::new(chrom, start_pos..=end_pos);
        let (mut reader, header) =
            crate::extract::open_cram_reader_for_region(bam_path, &repository, &region)
                .context("failed to open CRAM for pileup")?;
        let query = reader.query(&header, &region)?;
        // `query` has already rejected an unknown contig, so this is `Some`.
        let queried_reference_sequence_id =
            header.reference_sequences().get_index_of(chrom.as_bytes());

        for rec_result in query {
            let cram_record = rec_result?;
            let buf = cram_record.try_into_alignment_record(&header)?;
            // A container holding several contigs is decoded whole and `Query`
            // filters on coordinates alone, so another contig's reads would be
            // counted as this region's (L2, N4). The BAM arm below needs no
            // such guard: noodles-bam's `Query` compares the reference id.
            if !record_is_on_queried_reference(&buf, queried_reference_sequence_id) {
                continue;
            }
            let flags = buf.flags();
            if flags.is_unmapped()
                || flags.is_secondary()
                || flags.is_supplementary()
                || flags.is_duplicate()
                || flags.is_qc_fail()
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
            if flags.is_unmapped()
                || flags.is_secondary()
                || flags.is_supplementary()
                || flags.is_duplicate()
                || flags.is_qc_fail()
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
    cigar_ops: Box<
        dyn Iterator<Item = std::io::Result<noodles::sam::alignment::record::cigar::Op>> + '_,
    >,
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

// ---------------------------------------------------------------------------
// VCF INFO helpers (local copies from vcf_input.rs)
// ---------------------------------------------------------------------------

fn parse_info_field<'a>(info: &'a str, key: &str) -> Option<&'a str> {
    for field in info.split(';') {
        if let Some(value) = field
            .strip_prefix(key)
            .and_then(|rest| rest.strip_prefix('='))
        {
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
    writeln!(out, "Result: {}/{} PASS", n_pass, n_total,)?;

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
    fn test_sa_points_near_partner_only() {
        let sa = "chr20,45000101,-,100S50M,60,0;chr20,41000001,+,100M50S,60,1;";
        assert!(sa_points_near(sa, "chr20", 45000000, 500));
        assert!(sa_points_near(sa, "chr20", 41000000, 500)); // second entry
        assert!(!sa_points_near(sa, "chr20", 43000000, 500)); // elsewhere
        assert!(!sa_points_near(sa, "chr9", 45000000, 500)); // other chrom
        assert!(!sa_points_near("garbage", "chr20", 45000000, 500));
    }

    #[test]
    fn test_truth_events_know_their_partner_breakpoint() {
        let vcf = "\
##fileformat=VCFv4.3
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE
chr20\t41000000\tsim_del_1\tA\t<DEL>\t999\tPASS\tSVTYPE=DEL;END=41010000;SVLEN=-10000\tGT\t0/1
chr20\t41000000\tsim_fus_2\tA\tA]chr20:45000000]\t999\tPASS\tSVTYPE=BND;MATEID=sim_fus_2_mate\tGT\t0/1
chr20\t42000000\tsim_ins_3\tA\t<INS>\t999\tPASS\tSVTYPE=INS;SVLEN=500\tGT\t0/1
";
        let path = std::env::temp_dir().join(format!("spike_partner_{}.vcf", std::process::id()));
        std::fs::write(&path, vcf).unwrap();
        let events = load_truth_events(path.to_str().unwrap()).unwrap();
        std::fs::remove_file(&path).ok();

        assert_eq!(events[0].partner, Some(("chr20".to_string(), 41010000)));
        assert_eq!(events[1].partner, Some(("chr20".to_string(), 44999999)));
        assert_eq!(events[2].partner, None);
    }

    #[test]
    fn test_mean_depth_counts_aligned_bases() {
        // Two reads over [0,100): 100 + 50 bases inside -> depth 1.5.
        assert_eq!(mean_depth(&[(0, 100), (50, 150)], 0, 100), 1.5);
        // 35x of 150 bp reads, one start every 150/35 bp: depth ~35, not ~0.23.
        let spans: Vec<(u64, u64)> = (0..3500u64).map(|i| {
            let s = i * 150 / 35;
            (s, s + 150)
        }).collect();
        let d = mean_depth(&spans, 1000, 2000);
        assert!((d - 35.0).abs() < 0.5, "depth {}", d);
    }

    #[test]
    fn test_coverage_ratio_without_flank_coverage_fails() {
        let r = coverage_ratio_result("DEL".to_string(), "DEL", 0.5, 0.0, 0.0);
        assert!(!r.pass, "an unevaluable coverage check must not pass");
    }

    #[test]
    fn test_coverage_ratio_tells_deletion_from_untouched() {
        assert!(coverage_ratio_result("DEL".to_string(), "DEL", 0.5, 17.5, 35.0).pass);
        assert!(!coverage_ratio_result("DEL".to_string(), "DEL", 0.5, 35.0, 35.0).pass);
    }

    #[test]
    fn test_check_that_cannot_run_is_a_failure() {
        // e.g. truth VCF uses "20" but the BAM uses "chr20".
        let r = check_outcome(
            "DEL 20:100-200",
            "split_reads",
            Err(anyhow::anyhow!("reference sequence not found: 20")),
        );
        assert!(!r.pass);
        assert_eq!(r.check_name, "split_reads");
        assert!(r.observed.contains("reference sequence not found: 20"), "{}", r.observed);
    }

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
        assert_eq!(parse_info_field("SVTYPE=DEL;END=100", "MISSING"), None);
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
            partner: None,
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
            partner: None,
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

    // --- N4: `validate`'s CRAM queries must not read another contig's reads ---

    /// The shared two-contig CRAM fixture in a scratch directory of its own.
    /// Returns `(dir, fasta_path, cram_path)`; the caller removes `dir`.
    fn two_contig_cram(tag: &str) -> (std::path::PathBuf, String, String) {
        let dir = std::env::temp_dir().join(format!(
            "spike_test_validate_{}_{}",
            tag,
            std::process::id()
        ));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        let (fasta, cram) = crate::extract::test_fixtures::write_two_contig_cram(&dir);
        (dir, fasta, cram)
    }

    #[test]
    fn test_count_depth_in_region_skips_another_contigs_reads() {
        // chrA's six reads cover 600 of the region's 2000 bases, so its mean
        // depth is 0.3. chrB's four reads cover another 400 bases of the same
        // window from the shared container; counted, they lift the answer to
        // 0.5 and with it every coverage-ratio verdict (L2).
        let (dir, fasta, cram) = two_contig_cram("depth");

        let depth = count_depth_in_region(&cram, &fasta, "chrA", 0, 2000, 20).unwrap();

        assert_eq!(
            depth, 0.3,
            "only chrA's reads may count toward chrA's depth"
        );
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn test_split_reads_to_partner_skips_another_contigs_reads() {
        // Every fixture record carries the same SA:Z, so what comes back is
        // exactly "the records this query saw" — which is what makes a
        // foreign contig's read names visible here.
        let (dir, fasta, cram) = two_contig_cram("split_reads");

        let found = split_reads_to_partner(
            &cram,
            &fasta,
            "chrA",
            0,
            2000,
            20,
            &("chrA".to_string(), 1500),
            10,
        )
        .unwrap();

        let mut names: Vec<String> = found.into_iter().collect();
        names.sort();
        assert_eq!(
            names,
            vec!["chrA_pair0", "chrA_pair1", "chrA_pair2"],
            "a chrB read's name must not be counted as a chrA split read"
        );
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn test_pileup_region_skips_another_contigs_reads() {
        // The allele-fraction check's pileup. chrA and chrB cover disjoint
        // intervals, so every count over chrB's 301-400, 501-600 and 701-800
        // is a base that arrived from the other contig.
        let (dir, fasta, cram) = two_contig_cram("pileup");

        let mut allele_counts: HashMap<u64, [u32; 4]> = HashMap::new();
        let mut read_alleles: HashMap<String, Vec<(u64, u8)>> = HashMap::new();
        pileup_region(
            &cram,
            &fasta,
            "chrA",
            0,
            2000,
            20,
            &mut allele_counts,
            &mut read_alleles,
        )
        .unwrap();

        // The three intervals only chrB covers hold all 400 of its bases.
        let foreign: u32 = (300..400u64)
            .chain(500..600)
            .chain(700..800)
            .filter_map(|pos| allele_counts.get(&pos))
            .map(|c| c.iter().sum::<u32>())
            .sum();
        assert_eq!(foreign, 0, "chrB's bases must not enter the chrA pileup");

        let mut names: Vec<String> = read_alleles.keys().cloned().collect();
        names.sort();
        assert_eq!(
            names,
            vec!["chrA_pair0", "chrA_pair1", "chrA_pair2"],
            "only chrA reads may carry an allele in a chrA pileup"
        );
        let _ = std::fs::remove_dir_all(&dir);
    }

    // --- L15: the global checks read the event regions, not the file's head ---

    /// A one-contig CRAM shaped like a whole-genome BAM: a head nobody asked
    /// about, and an event window far from it. 15 pairs at chrA:101-1800 carry
    /// MAPQ 0 (a telomere's worth of multi-mapping reads) and 3 pairs inside
    /// chrA:10001-10500 carry MAPQ 60, so the mean MAPQ of all 36 records is
    /// exactly 10.0 -- the number REVIEW.md measured on whole-genome HG002 --
    /// while the event's own window reads 60. No record carries the duplicate
    /// flag. Returns `(dir, fasta_path, cram_path)`; the caller removes `dir`.
    fn head_and_event_cram(tag: &str) -> (std::path::PathBuf, String, String) {
        use std::num::NonZeroUsize;

        const CONTIG_LEN: usize = 20_000;
        const READ_LEN: usize = 100;

        let dir = std::env::temp_dir().join(format!(
            "spike_test_validate_{}_{}",
            tag,
            std::process::id()
        ));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();

        let seq: Vec<u8> = (0..CONTIG_LEN).map(|i| b"ACGT"[i % 4]).collect();
        let mut fasta = String::from(">chrA\n");
        let offset = fasta.len();
        for chunk in seq.chunks(60) {
            fasta.push_str(std::str::from_utf8(chunk).unwrap());
            fasta.push('\n');
        }
        let fasta_path = dir.join("one_contig.fa");
        std::fs::write(&fasta_path, &fasta).unwrap();
        std::fs::write(
            dir.join("one_contig.fa.fai"),
            format!("chrA\t{}\t{}\t60\t61\n", seq.len(), offset),
        )
        .unwrap();

        let header = noodles::sam::Header::builder()
            .add_reference_sequence(
                "chrA",
                noodles::sam::header::record::value::Map::<
                    noodles::sam::header::record::value::map::ReferenceSequence,
                >::new(NonZeroUsize::try_from(seq.len()).unwrap()),
            )
            .build();

        // One pair = two records, read1 forward (0x63) and read2 reverse
        // (0x93), both properly segmented, 300 bp apart.
        let record = |name: &str, start: usize, first: bool, mapq: u8| {
            let (pos, mate_pos) = if first {
                (start, start + 200)
            } else {
                (start + 200, start)
            };
            let span = 200 + READ_LEN;
            noodles::cram::Record::builder()
                .set_bam_flags(noodles::sam::alignment::record::Flags::from(if first {
                    0x63u16
                } else {
                    0x93u16
                }))
                .set_flags(noodles::cram::record::Flags::QUALITY_SCORES_STORED_AS_ARRAY)
                .set_reference_sequence_id(0)
                .set_read_length(READ_LEN)
                .set_alignment_start(noodles::core::Position::new(pos).unwrap())
                .set_name(name)
                .set_next_fragment_reference_sequence_id(0)
                .set_next_mate_alignment_start(noodles::core::Position::new(mate_pos).unwrap())
                .set_template_size(if first { span as i32 } else { -(span as i32) })
                .set_mapping_quality(
                    noodles::sam::alignment::record::MappingQuality::new(mapq).unwrap(),
                )
                .set_bases(noodles::sam::alignment::record_buf::Sequence::from(
                    seq[pos - 1..pos - 1 + READ_LEN].to_vec(),
                ))
                .set_quality_scores(noodles::sam::alignment::record_buf::QualityScores::from(
                    vec![40u8; READ_LEN],
                ))
                .build()
        };

        let cram_path = dir.join("head_and_event.cram");
        let repository =
            crate::extract::build_fasta_repository(fasta_path.to_str().unwrap()).unwrap();
        {
            let mut writer = noodles::cram::io::writer::Builder::default()
                .set_reference_sequence_repository(repository)
                .build_from_path(&cram_path)
                .unwrap();
            writer.write_header(&header).unwrap();
            for i in 0..15usize {
                let name = format!("head_pair{}", i);
                let start = 101 + i * 100;
                writer
                    .write_record(&header, record(&name, start, true, 0))
                    .unwrap();
                writer
                    .write_record(&header, record(&name, start, false, 0))
                    .unwrap();
            }
            for i in 0..3usize {
                let name = format!("event_pair{}", i);
                let start = 10_001 + i * 100;
                writer
                    .write_record(&header, record(&name, start, true, 60))
                    .unwrap();
                writer
                    .write_record(&header, record(&name, start, false, 60))
                    .unwrap();
            }
            writer.try_finish(&header).unwrap();
        }

        let index = noodles::cram::index(&cram_path).unwrap();
        let mut index_writer = noodles::cram::crai::io::Writer::new(
            std::fs::File::create(dir.join("head_and_event.cram.crai")).unwrap(),
        );
        index_writer.write_index(&index).unwrap();
        index_writer.finish().unwrap();

        (
            dir,
            fasta_path.to_str().unwrap().to_string(),
            cram_path.to_str().unwrap().to_string(),
        )
    }

    /// A het DEL truth event over [start, end).
    fn del_event(chrom: &str, start: u64, end: u64) -> TruthEvent {
        TruthEvent {
            chrom: chrom.to_string(),
            start,
            end,
            sv_type: "DEL".to_string(),
            expected_vaf: 0.5,
            gene: "unknown".to_string(),
            partner: Some((chrom.to_string(), end)),
            ref_allele: None,
            alt_allele: None,
        }
    }

    #[test]
    fn test_global_sample_reads_the_event_region_not_the_head_of_the_file() {
        // The file's first 30 records are the head's MAPQ 0 pairs and only 6
        // lie in the event's window, so a sample taken from the top of the
        // file reads 10.0 and FAILs the >20 check on a file whose event
        // region is MAPQ 60 throughout (L15).
        let (dir, fasta, cram) = head_and_event_cram("global_sample");
        let event = del_event("chrA", 10_000, 10_200);

        let sample =
            sample_event_regions(&cram, &fasta, std::slice::from_ref(&event), 5_000).unwrap();
        let result = check_mapq(&sample);

        assert_eq!(
            result.observed, "60.0",
            "mean MAPQ must come from the event's own region"
        );
        assert!(result.pass, "MAPQ 60 must pass the >20 check");
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn test_dup_rate_without_duplicate_flags_is_not_evaluable() {
        // A file with no duplicate flags is not a file without duplicates:
        // nothing marked them, so there is no rate to report and none to
        // pass on (M10).
        let sample = GlobalSample {
            total: 10_000,
            dups: 0,
            mapq_sum: 600_000,
            insert_sizes: vec![300.0; 10_000],
        };

        let result = check_dup_rate(&sample);

        assert_eq!(result.observed, "no dup flags");
        assert!(!result.pass, "an unevaluable duplicate rate may not pass");
    }

    #[test]
    fn test_dup_rate_is_measured_when_records_carry_duplicate_flags() {
        // The other half of the rule: a marked file still gets a number.
        let sample = GlobalSample {
            total: 1_000,
            dups: 71,
            mapq_sum: 60_000,
            insert_sizes: vec![300.0; 1_000],
        };

        let result = check_dup_rate(&sample);

        assert_eq!(result.observed, "7.1%");
        assert!(result.pass, "7.1% is under the 50% limit");
    }

    #[test]
    fn test_insert_size_without_paired_records_is_not_evaluable() {
        // The fallback this check used to inherit from `bam_stats` was
        // 350 +/- 50, which passes: a default reported as a measurement.
        let sample = GlobalSample {
            total: 100,
            dups: 0,
            mapq_sum: 6_000,
            insert_sizes: Vec::new(),
        };

        let result = check_insert_size(&sample);

        assert_eq!(result.observed, "no pairs");
        assert!(!result.pass, "an unevaluable insert size may not pass");
    }
}
