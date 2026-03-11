//! spike: Haplotype-based read spike-in simulator for genomic variants.
//!
//! Creates synthetic chimeric reads from real BAM data for any variant type
//! (SVs, fusions, SNPs, indels).

mod bam_stats;
mod exon;
mod extract;
mod fastq;
mod haplotype;
mod loh;
mod reference;
mod simulate;
mod stats;
mod synth;
mod truth;
mod types;
mod validate;
mod vcf_input;

use anyhow::{bail, Result};
use clap::Parser;
use rand::rngs::StdRng;
use rand::Rng;
use rand::SeedableRng;
use rand_distr::{Beta, Distribution};
use std::path::Path;

use exon::AfSpec;
use haplotype::VariantHaplotype;
use types::{ReadPair, ReadPool, SimConfig, SimEvent};

#[derive(Parser)]
#[command(
    name = "spike",
    about = "Haplotype-based read spike-in simulator for genomic variants",
    version
)]
struct Args {
    /// Input BAM file (coordinate-sorted, indexed).
    #[arg(short, long)]
    bam: String,

    /// Reference FASTA (with .fai index).
    #[arg(short, long)]
    reference: String,

    /// Event specification(s). Can be repeated.
    /// Formats:
    ///   --event "del:chr20:30000000-30005000"           (coordinate-based)
    ///   --event "del:GENE:exon4-exon8"                  (gene-based, requires --exon-bed)
    ///   --event "dup:GENE:exon4-exon8"                  (gene-based duplication)
    ///   --event "inv:GENE:exon4-exon8"                  (gene-based inversion)
    ///   --event "fusion:GENEA:exon14:GENEB:exon2"       (fusion, requires --exon-bed)
    ///   --event "fusion:GENEA:exon14:GENEB:exon2:inv"   (inverted fusion)
    ///   --event "dup:chr20:30000000-30005000"
    ///   --event "inv:chr20:30000000-30005000"
    ///   --event "ins:chr20:30000000:500"                (random insertion sequence)
    ///   --event "ins:chr20:30000000:ACGTACGT"           (explicit insertion sequence)
    ///   --event "snp:chr20:30000000:A:T"                (SNP/small variant, POS is 1-based)
    ///   --event "snp:chr20:30000000:A>T"                (alternate syntax, POS is 1-based)
    ///   --event "snp:chr20:30000000:ACG:A"              (small deletion, POS is 1-based)
    ///   --event "snp:chr20:30000000:A:ACGT"             (small insertion, POS is 1-based)
    /// Per-event AF (appended with ;):
    ///   --event "del:GENE:exon4-exon8;af=0.15"
    ///   --event "fusion:GENEA:exon14:GENEB:exon2;af=het"
    ///   af=<number>: exact AF, af=het: Beta(40,40)~0.5, af=hom: 1.0
    #[arg(short, long)]
    event: Vec<String>,

    /// Input VCF file with variant records. Supports DEL, INS, DUP, INV, BND,
    /// and standard SNP/indel records (no SVTYPE, explicit REF/ALT alleles).
    /// Can be combined with --event. At least one of --event or --vcf required.
    #[arg(long)]
    vcf: Option<String>,

    /// Exon BED file. Required when using gene-based --event specs (e.g. "del:GENE:exon4-exon8").
    #[arg(long)]
    exon_bed: Option<String>,

    /// Target allele fraction (0.0-1.0).
    #[arg(long, default_value_t = 0.5)]
    allele_fraction: f64,

    /// Output directory.
    #[arg(short, long, default_value = "/tmp/spike")]
    output: String,

    /// Random seed for reproducibility.
    #[arg(long, default_value_t = 42)]
    seed: u64,

    /// Number of threads for BAM reading.
    #[arg(short, long, default_value_t = 4)]
    threads: usize,

    /// Read extraction region (e.g. "chr19:11080000-11140000").
    /// When specified, reads are extracted from this region instead of event ± flank.
    /// Use this to ensure the output BAM covers the full gene/region of interest.
    #[arg(long)]
    region: Option<String>,

    /// Flanking region (bp) to include around events (used when --region is not set).
    #[arg(long, default_value_t = 10000)]
    flank: u64,

    /// Minimum mapping quality for donor reads.
    #[arg(long, default_value_t = 20)]
    min_mapq: u8,

    /// Aligner for alignment script. Presets: "bwa-mem2" (default),
    /// "minimap2", "bowtie2", or a custom command that accepts
    /// <ref> <r1.fq.gz> <r2.fq.gz> and produces SAM on stdout.
    #[arg(long, default_value = "bwa-mem2")]
    aligner: String,

    /// Path to samtools binary (used for sort/index in align script).
    #[arg(long, default_value = "samtools")]
    samtools: String,

    /// Automatically run alignment after FASTQ generation.
    #[arg(long)]
    align: bool,

    /// Indel error rate per base in synthetic reads (fraction of total error
    /// that is indel rather than substitution). Default 0.0 means substitution-only.
    /// Typical Illumina: 0.0 to 0.05.
    #[arg(long, default_value_t = 0.0)]
    indel_error_rate: f64,

    /// Optional gVCF/VCF with SNP calls for LOH simulation.
    /// When provided, het SNP positions are extracted from this file to
    /// determine which reads belong to the deleted haplotype (more accurate
    /// than the default pileup-based approach). Supports .vcf and .vcf.gz
    /// (requires bcftools in PATH for .vcf.gz).
    #[arg(long)]
    gvcf: Option<String>,

    /// Allow overlapping events on the same chromosome.
    ///
    /// By default, overlapping events are rejected to keep event effects
    /// independent and truth interpretation unambiguous.
    #[arg(long)]
    allow_overlap: bool,

    /// Duplication model: "full" (default) builds a full tandem haplotype
    /// with duplicated region appearing twice, producing both junction reads
    /// and correct depth increase from a single tiling pass. "junction" uses
    /// the legacy junction-only haplotype with separate depth copies.
    #[arg(long, default_value = "full")]
    dup_model: String,
}

/// Parsed extraction region from --region flag.
struct ExtractionRegion {
    chrom: String,
    start: u64,
    end: u64,
}

/// Compute read extraction bounds for an event.
///
/// If `--region` is set and the event is on the same chromosome, use the region bounds
/// (expanded to also include event ± flank if the event extends beyond the region).
/// Otherwise, fall back to event ± flank.
fn extraction_bounds(
    event_chrom: &str,
    event_start: u64,
    event_end: u64,
    flank: u64,
    region: &Option<ExtractionRegion>,
) -> (u64, u64) {
    if let Some(r) = region {
        if r.chrom == event_chrom {
            // Union of region and event+flank to ensure we cover both.
            let start = r.start.min(event_start.saturating_sub(flank));
            let end = r.end.max(event_end.saturating_add(flank));
            return (start, end);
        }
    }
    // Fallback: event ± flank.
    (
        event_start.saturating_sub(flank),
        event_end.saturating_add(flank),
    )
}

/// Parse a region string like "chr19:11080000-11140000" into (chrom, start, end).
/// Coordinates are 1-based inclusive (like samtools), converted to 0-based half-open internally.
fn parse_region(s: &str) -> Result<ExtractionRegion> {
    let (chrom, coords) = s
        .split_once(':')
        .ok_or_else(|| anyhow::anyhow!("invalid region '{}', expected chr:start-end", s))?;
    let (start_s, end_s) = coords
        .split_once('-')
        .ok_or_else(|| anyhow::anyhow!("invalid region '{}', expected chr:start-end", s))?;
    let start: u64 = start_s
        .replace(',', "")
        .parse()
        .map_err(|_| anyhow::anyhow!("invalid start in region '{}'", s))?;
    let end: u64 = end_s
        .replace(',', "")
        .parse()
        .map_err(|_| anyhow::anyhow!("invalid end in region '{}'", s))?;
    if start == 0 {
        anyhow::bail!("region start must be >= 1 (1-based), got 0 in '{}'", s);
    }
    if start > end {
        anyhow::bail!(
            "region start > end ({} > {}) in '{}'; check your interval",
            start, end, s
        );
    }
    // Convert from 1-based inclusive to 0-based half-open.
    Ok(ExtractionRegion {
        chrom: chrom.to_string(),
        start: start.saturating_sub(1),
        end,
    })
}

fn main() -> Result<()> {
    // Intercept "validate" subcommand before clap parses (backward-compatible).
    let raw_args: Vec<String> = std::env::args().collect();
    if raw_args.get(1).is_some_and(|s| s == "validate") {
        match validate::run() {
            Ok(()) => return Ok(()),
            Err(e) => {
                // Log the error (table output already printed), exit non-zero.
                log::error!("{}", e);
                std::process::exit(1);
            }
        }
    }

    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info")).init();

    let args = Args::parse();

    // Validate inputs.
    if args.allele_fraction <= 0.0 || args.allele_fraction > 1.0 {
        bail!("allele-fraction must be in (0.0, 1.0]");
    }

    // Create output directory.
    std::fs::create_dir_all(&args.output)?;

    // Load gene targets (only needed for gene-based --event specs).
    let genes = if let Some(bed_path) = &args.exon_bed {
        exon::parse_exon_bed(bed_path)?
    } else {
        Vec::new()
    };

    // Parse event specifications from --event flags.
    let parsed: Vec<(SimEvent, Option<AfSpec>)> = args
        .event
        .iter()
        .map(|spec| exon::parse_event_spec(spec, &genes))
        .collect::<Result<Vec<_>>>()?;

    // Load events from --vcf if provided.
    let vcf_events: Vec<SimEvent> = if let Some(vcf_path) = &args.vcf {
        vcf_input::load_events_from_vcf(vcf_path)?
    } else {
        Vec::new()
    };

    if parsed.is_empty() && vcf_events.is_empty() {
        bail!("no events specified (use --event or --vcf)");
    }

    log::info!(
        "Parsed {} --event spec(s) + {} VCF event(s)",
        parsed.len(),
        vcf_events.len(),
    );

    // Validate --dup-model.
    if args.dup_model != "full" && args.dup_model != "junction" {
        bail!(
            "invalid --dup-model '{}', expected 'full' or 'junction'",
            args.dup_model,
        );
    }

    // Compute BAM stats for read length.
    let bam_stats =
        crate::bam_stats::compute_stats(&args.bam, 50_000, Some(args.reference.as_str()))?;
    let read_length = bam_stats.read_length.round() as usize;

    let config = SimConfig {
        bam_path: args.bam.clone(),
        ref_path: args.reference.clone(),
        allele_fraction: args.allele_fraction,
        flank_bp: args.flank,
        read_length,
        min_mapq: args.min_mapq,
        gvcf_path: args.gvcf.clone(),
        indel_error_rate: args.indel_error_rate,
        dup_model: args.dup_model.clone(),
    };

    // Parse --region if provided.
    let extraction_region = if let Some(ref region_str) = args.region {
        let r = parse_region(region_str)?;
        log::info!(
            "Extraction region: {}:{}-{} (0-based half-open)",
            r.chrom,
            r.start,
            r.end,
        );
        Some(r)
    } else {
        None
    };

    let mut rng = StdRng::seed_from_u64(args.seed);

    // Resolve per-event AF specs into concrete values.
    let mut events: Vec<SimEvent> = parsed
        .into_iter()
        .map(|(mut event, af_spec)| {
            let resolved_af = match af_spec {
                Some(AfSpec::Exact(v)) => Some(v),
                Some(AfSpec::Het) => {
                    let beta = Beta::new(40.0, 40.0).unwrap();
                    let v = beta.sample(&mut rng);
                    log::info!("  Het AF sampled: {:.3}", v);
                    Some(v)
                }
                Some(AfSpec::Hom) => Some(1.0),
                None => None, // will use global default
            };
            event.set_allele_fraction(resolved_af);
            event
        })
        .collect();

    // Append VCF-sourced events (AF already embedded from VCF INFO).
    events.extend(vcf_events);

    // Load reference sequence for synthetic read generation.
    let chroms_needed: Vec<String> = events
        .iter()
        .flat_map(|e| match e {
            SimEvent::Deletion { chrom, .. }
            | SimEvent::Duplication { chrom, .. }
            | SimEvent::Inversion { chrom, .. }
            | SimEvent::Insertion { chrom, .. }
            | SimEvent::SmallVariant { chrom, .. } => vec![chrom.clone()],
            SimEvent::Fusion {
                chrom_a, chrom_b, ..
            } => vec![chrom_a.clone(), chrom_b.clone()],
        })
        .collect::<std::collections::HashSet<_>>()
        .into_iter()
        .collect();
    let chrom_refs: Vec<&str> = chroms_needed.iter().map(|s| s.as_str()).collect();
    let shared_ref = crate::reference::SharedReference::load(&config.ref_path, &chrom_refs)?;

    // Validate event coordinates against reference chromosome lengths.
    validate_event_coordinates(&events, &shared_ref)?;
    validate_event_overlaps(&events, args.allow_overlap)?;

    let mut all_output_pairs: Vec<ReadPair> = Vec::new();

    // Per-event stats collected for README/log output.
    struct EventStat {
        vaf: f64,
        kept: usize,
        chimeric: usize,
        suppressed: usize,
    }
    let mut event_stats: Vec<EventStat> = Vec::new();

    // Haplotype flank: must be >= max expected fragment length so reads near
    // haplotype edges form complete pairs. 2000bp is conservative.
    let hap_flank: u64 = 2000;

    // Process each event using the unified haplotype + tiling approach.
    for (i, event) in events.iter().enumerate() {
        log::info!("Processing event {}/{}: {:?}", i + 1, events.len(), event);
        let vaf = event.allele_fraction().unwrap_or(config.allele_fraction);
        log::info!("  Using VAF={:.3} for this event", vaf);

        // Extract reads and build pool.
        let (pool, _extraction_chrom) = extract_pool_for_event(event, &config, &extraction_region)?;

        // Build quality profile and synth generator.
        let quality_profile =
            synth::QualityProfile::from_read_pairs(&pool.pairs, config.read_length);
        let mut synth_gen = synth::SynthReadGenerator::new(
            quality_profile,
            &shared_ref,
            config.read_length,
            config.indel_error_rate,
        );

        // Build variant haplotype.
        let mut haplotype =
            build_haplotype(event, &shared_ref, hap_flank, &config.dup_model, &mut rng)?;

        // Simulate: suppress reads + tile synthetic reads across haplotype.
        let output = simulate::simulate_event(
            i + 1,
            event,
            &pool,
            &mut haplotype,
            &config,
            &mut synth_gen,
            vaf,
            &mut rng,
        )?;

        log::info!(
            "Event {}: {} kept + {} chimeric, {} suppressed",
            i + 1,
            output.kept_originals.len(),
            output.chimeric_pairs.len(),
            output.suppressed_count,
        );

        event_stats.push(EventStat {
            vaf,
            kept: output.kept_originals.len(),
            chimeric: output.chimeric_pairs.len(),
            suppressed: output.suppressed_count,
        });

        all_output_pairs.extend(output.kept_originals);
        all_output_pairs.extend(output.chimeric_pairs);
    }

    // Deduplicate by read name (in case of overlapping regions).
    let before_dedup = all_output_pairs.len();
    dedup_by_name(&mut all_output_pairs);
    if all_output_pairs.len() < before_dedup {
        log::info!(
            "Deduplicated: {} -> {} pairs",
            before_dedup,
            all_output_pairs.len(),
        );
    }

    // Write FASTQ.
    let (r1_path, r2_path) = fastq::write_paired_fastq(&all_output_pairs, &args.output)?;

    // Write truth VCF.
    let truth_path = Path::new(&args.output).join("truth.vcf");
    truth::write_truth_vcf(
        &events,
        config.allele_fraction, // default AF for events without per-event override
        &truth_path.to_string_lossy(),
        &args.reference,
    )?;

    // Write alignment convenience script.
    write_align_script(
        &args.output,
        &args.reference,
        args.threads,
        &args.aligner,
        &args.samtools,
    )?;

    // Write events BED (used by merge.sh).
    write_event_bed(&args.output, &events, args.flank)?;

    // Write merge script.
    write_merge_script(&args.output, &args.bam, &args.reference, args.threads, &args.samtools)?;

    // Collect per-event stats as (vaf, kept, chimeric, suppressed) tuples for README.
    let stats_tuples: Vec<(f64, usize, usize, usize)> = event_stats
        .iter()
        .map(|s| (s.vaf, s.kept, s.chimeric, s.suppressed))
        .collect();

    // Write README.md.
    let cmdline = std::env::args().collect::<Vec<_>>().join(" ");
    write_readme(
        &args.output,
        &cmdline,
        &args.bam,
        &args.reference,
        &events,
        &stats_tuples,
        all_output_pairs.len(),
        args.flank,
    )?;

    // Summary.
    log::info!("=== spike complete ===");
    log::info!("Output directory: {}", args.output);
    log::info!("FASTQ: {} and {}", r1_path, r2_path);
    log::info!("Truth VCF: {}", truth_path.display());
    log::info!("Total read pairs: {}", all_output_pairs.len());
    log::info!("README: {}/README.md", args.output);
    log::info!(
        "Next steps: bash {}/align.sh  →  bash {}/merge.sh",
        args.output,
        args.output
    );

    // Auto-align if requested.
    if args.align {
        run_alignment(&args.output, &args.reference, args.threads)?;
    }

    Ok(())
}

/// Deduplicate read pairs by name, keeping the last occurrence.
///
/// Keeping the last occurrence makes event-order behavior explicit when
/// multiple simulated events touch the same original read name.
fn dedup_by_name(pairs: &mut Vec<ReadPair>) {
    let mut last_idx: std::collections::HashMap<String, usize> =
        std::collections::HashMap::with_capacity(pairs.len());
    for (i, p) in pairs.iter().enumerate() {
        last_idx.insert(p.name.clone(), i);
    }
    let mut out = Vec::with_capacity(last_idx.len());
    for (i, p) in pairs.drain(..).enumerate() {
        if last_idx.get(&p.name).copied() == Some(i) {
            out.push(p);
        }
    }
    *pairs = out;
}

/// Validate overlap relationships between single-region events.
///
/// Overlaps are rejected by default to keep multi-event simulations independent.
/// When `allow_overlap` is true, overlaps are allowed but logged as warnings.
fn validate_event_overlaps(events: &[SimEvent], allow_overlap: bool) -> Result<()> {
    let mut overlaps: Vec<(usize, usize, String, u64, u64, u64, u64)> = Vec::new();

    for i in 0..events.len() {
        let regions_i = event_regions_for_overlap(&events[i]);
        for j in (i + 1)..events.len() {
            let regions_j = event_regions_for_overlap(&events[j]);
            for (chrom_i, start_i, end_i) in &regions_i {
                for (chrom_j, start_j, end_j) in &regions_j {
                    if chrom_i != chrom_j {
                        continue;
                    }
                    let has_overlap = start_i < end_j && start_j < end_i;
                    if has_overlap {
                        overlaps.push((
                            i + 1,
                            j + 1,
                            chrom_i.to_string(),
                            *start_i,
                            *end_i,
                            *start_j,
                            *end_j,
                        ));
                    }
                }
            }
        }
    }

    if overlaps.is_empty() {
        return Ok(());
    }

    if allow_overlap {
        for (i, j, chrom, start_i, end_i, start_j, end_j) in overlaps {
            log::warn!(
                "Events {} and {} overlap on {} ({}-{} vs {}-{}). \
                 Overlap composition is approximate.",
                i,
                j,
                chrom,
                start_i,
                end_i,
                start_j,
                end_j,
            );
        }
        return Ok(());
    }

    let mut msg = String::from(
        "overlapping events detected (default is to reject overlaps).\n\
         Use --allow-overlap to override.\n",
    );
    for (i, j, chrom, start_i, end_i, start_j, end_j) in overlaps.iter().take(10) {
        msg.push_str(&format!(
            "  - events {} and {} overlap on {} ({}-{} vs {}-{})\n",
            i, j, chrom, start_i, end_i, start_j, end_j
        ));
    }
    if overlaps.len() > 10 {
        msg.push_str(&format!(
            "  ... and {} more overlap(s)\n",
            overlaps.len() - 10
        ));
    }

    bail!("{}", msg.trim_end());
}

/// Return one or more non-empty regions used for overlap checks.
///
/// Coordinates are 0-based half-open. Point events use a 1bp window.
fn event_regions_for_overlap(event: &SimEvent) -> Vec<(String, u64, u64)> {
    match event {
        SimEvent::Deletion {
            chrom,
            del_start,
            del_end,
            ..
        } => vec![(chrom.clone(), *del_start, *del_end)],
        SimEvent::Duplication {
            chrom,
            dup_start,
            dup_end,
            ..
        } => vec![(chrom.clone(), *dup_start, *dup_end)],
        SimEvent::Inversion {
            chrom,
            inv_start,
            inv_end,
            ..
        } => vec![(chrom.clone(), *inv_start, *inv_end)],
        SimEvent::Insertion { chrom, pos, .. } => {
            vec![(chrom.clone(), *pos, pos.saturating_add(1))]
        }
        SimEvent::SmallVariant {
            chrom,
            pos,
            ref_allele,
            ..
        } => vec![(
            chrom.clone(),
            *pos,
            pos.saturating_add(ref_allele.len() as u64),
        )],
        SimEvent::Fusion {
            chrom_a,
            bp_a,
            chrom_b,
            bp_b,
            ..
        } => vec![
            (chrom_a.clone(), *bp_a, bp_a.saturating_add(1)),
            (chrom_b.clone(), *bp_b, bp_b.saturating_add(1)),
        ],
    }
}

/// Write a convenience shell script for alignment.
///
/// Supports presets: bwa-mem2, minimap2, bowtie2, or a custom command.
fn write_align_script(
    output_dir: &str,
    ref_path: &str,
    threads: usize,
    aligner: &str,
    samtools: &str,
) -> Result<()> {
    let script_path = Path::new(output_dir).join("align.sh");

    let align_cmd = match aligner {
        "bwa-mem2" => "bwa-mem2 mem -t \"$THREADS\" \\\n    -R '@RG\\tID:sim\\tSM:SIM\\tPL:ILLUMINA' \\\n    \"$REF\" \\\n    \"$DIR/R1.fq.gz\" \"$DIR/R2.fq.gz\" \\\n    2>\"$DIR/align.log\"".to_string(),
        "minimap2" => "minimap2 -a -x sr -t \"$THREADS\" \\\n    -R '@RG\\tID:sim\\tSM:SIM\\tPL:ILLUMINA' \\\n    \"$REF\" \\\n    \"$DIR/R1.fq.gz\" \"$DIR/R2.fq.gz\" \\\n    2>\"$DIR/align.log\"".to_string(),
        "bowtie2" => "bowtie2 -x \"$REF\" \\\n    -1 \"$DIR/R1.fq.gz\" -2 \"$DIR/R2.fq.gz\" \\\n    -p \"$THREADS\" \\\n    --rg-id sim --rg SM:SIM --rg PL:ILLUMINA \\\n    2>\"$DIR/align.log\"".to_string(),
        custom => format!(
            "{custom} \"$REF\" \"$DIR/R1.fq.gz\" \"$DIR/R2.fq.gz\" \\\n    2>\"$DIR/align.log\"",
        ),
    };

    let script = format!(
        r#"#!/bin/bash
set -euo pipefail
# Align simulated reads and sort.
# Usage: bash align.sh [REF] [THREADS]
REF="${{1:-{ref_path}}}"
THREADS="${{2:-{threads}}}"
SAMTOOLS="{samtools}"
DIR="$(cd "$(dirname "$0")" && pwd)"

echo "Aligning $DIR/R1.fq.gz + R2.fq.gz ({aligner}, $THREADS threads)..."
{align_cmd} | \
    "$SAMTOOLS" sort -@ "$THREADS" -o "$DIR/sim.bam" -

"$SAMTOOLS" index "$DIR/sim.bam"

TOTAL=$("$SAMTOOLS" view -c "$DIR/sim.bam" 2>/dev/null || echo "?")
SA_COUNT=$("$SAMTOOLS" view "$DIR/sim.bam" | grep -c "SA:Z:" 2>/dev/null || echo "0")
echo "Done: $DIR/sim.bam ($TOTAL reads, $SA_COUNT with SA tags)"
"#
    );

    std::fs::write(&script_path, script)?;

    // Make executable.
    #[cfg(unix)]
    {
        use std::os::unix::fs::PermissionsExt;
        std::fs::set_permissions(&script_path, std::fs::Permissions::from_mode(0o755))?;
    }

    Ok(())
}

/// Run alignment.
fn run_alignment(output_dir: &str, ref_path: &str, threads: usize) -> Result<()> {
    log::info!("Running alignment...");

    let script_path = Path::new(output_dir).join("align.sh");
    let status = std::process::Command::new("bash")
        .arg(&script_path)
        .arg(ref_path)
        .arg(threads.to_string())
        .status()?;

    if !status.success() {
        bail!("alignment failed with exit code {:?}", status.code());
    }

    let bam_path = Path::new(output_dir).join("sim.bam");
    log::info!("Aligned BAM: {}", bam_path.display());
    Ok(())
}

/// Validate that all event coordinates fall within reference chromosome bounds.
fn validate_event_coordinates(
    events: &[SimEvent],
    reference: &crate::reference::SharedReference,
) -> Result<()> {
    for event in events {
        match event {
            SimEvent::Deletion {
                chrom,
                del_start,
                del_end,
                ..
            } => {
                validate_range(reference, chrom, *del_start, *del_end, "DEL")?;
            }
            SimEvent::Duplication {
                chrom,
                dup_start,
                dup_end,
                ..
            } => {
                validate_range(reference, chrom, *dup_start, *dup_end, "DUP")?;
            }
            SimEvent::Inversion {
                chrom,
                inv_start,
                inv_end,
                ..
            } => {
                validate_range(reference, chrom, *inv_start, *inv_end, "INV")?;
            }
            SimEvent::Insertion { chrom, pos, .. } => {
                validate_point(reference, chrom, *pos, "INS")?;
            }
            SimEvent::SmallVariant {
                chrom,
                pos,
                ref_allele,
                ..
            } => {
                let end = *pos + ref_allele.len() as u64;
                validate_range(reference, chrom, *pos, end, "SNP/indel")?;
                validate_ref_allele(reference, chrom, *pos, ref_allele)?;
            }
            SimEvent::Fusion {
                chrom_a,
                bp_a,
                chrom_b,
                bp_b,
                ..
            } => {
                validate_point(reference, chrom_a, *bp_a, "Fusion-A")?;
                validate_point(reference, chrom_b, *bp_b, "Fusion-B")?;
            }
        }
    }
    Ok(())
}

fn validate_range(
    reference: &crate::reference::SharedReference,
    chrom: &str,
    start: u64,
    end: u64,
    label: &str,
) -> Result<()> {
    if start >= end {
        bail!(
            "{} event on {} has invalid interval: start ({}) >= end ({})",
            label, chrom, start, end,
        );
    }
    if let Some(chrom_len) = reference.chromosome_length(chrom) {
        if start >= chrom_len {
            bail!(
                "{} event start on {} is at or beyond chromosome length ({} >= {})",
                label,
                chrom,
                start,
                chrom_len,
            );
        }
        if end > chrom_len {
            bail!(
                "{} event on {} extends beyond chromosome length ({} > {})",
                label,
                chrom,
                end,
                chrom_len,
            );
        }
    }
    Ok(())
}

/// Validate a single point coordinate (for INS and Fusion breakpoints).
fn validate_point(
    reference: &crate::reference::SharedReference,
    chrom: &str,
    pos: u64,
    label: &str,
) -> Result<()> {
    if let Some(chrom_len) = reference.chromosome_length(chrom) {
        if pos >= chrom_len {
            bail!(
                "{} event position on {} is at or beyond chromosome length ({} >= {})",
                label, chrom, pos, chrom_len,
            );
        }
    }
    Ok(())
}

/// Validate that the user-specified REF allele matches the actual reference sequence.
fn validate_ref_allele(
    reference: &crate::reference::SharedReference,
    chrom: &str,
    pos: u64,
    ref_allele: &[u8],
) -> Result<()> {
    let end = pos + ref_allele.len() as u64;
    let actual = reference.fetch_sequence(chrom, pos, end)?;
    let actual_upper: Vec<u8> = actual.iter().map(|b| b.to_ascii_uppercase()).collect();
    let expected_upper: Vec<u8> = ref_allele.iter().map(|b| b.to_ascii_uppercase()).collect();
    if actual_upper != expected_upper {
        bail!(
            "REF allele mismatch at {}:{}-{}: specified '{}' but reference has '{}'. \
             Check that the position is correct (1-based in event spec) and matches the reference genome.",
            chrom,
            pos + 1, // display as 1-based
            pos + ref_allele.len() as u64,
            String::from_utf8_lossy(&expected_upper),
            String::from_utf8_lossy(&actual_upper),
        );
    }
    Ok(())
}

/// Extract reads and build a pool for a given event.
///
/// For fusions, extracts from both gene regions and merges pools.
fn extract_pool_for_event(
    event: &SimEvent,
    config: &SimConfig,
    extraction_region: &Option<ExtractionRegion>,
) -> Result<(ReadPool, String)> {
    if let SimEvent::Fusion {
        chrom_a,
        bp_a,
        chrom_b,
        bp_b,
        ..
    } = event
    {
        // Fusion: extract from both gene regions and merge pools.
        let (region_a_start, region_a_end) =
            extraction_bounds(chrom_a, *bp_a, *bp_a, config.flank_bp, extraction_region);
        let (region_b_start, region_b_end) =
            extraction_bounds(chrom_b, *bp_b, *bp_b, config.flank_bp, extraction_region);

        let pairs_a = extract::extract_read_pairs(
            &config.bam_path,
            chrom_a,
            region_a_start,
            region_a_end,
            config.min_mapq,
            Some(config.ref_path.as_str()),
        )?;
        let pairs_b = extract::extract_read_pairs(
            &config.bam_path,
            chrom_b,
            region_b_start,
            region_b_end,
            config.min_mapq,
            Some(config.ref_path.as_str()),
        )?;

        let mut all_pairs = pairs_a;
        all_pairs.extend(pairs_b);
        let frag_dist = stats::FragmentDist::from_read_pairs(&all_pairs);
        let pool = extract::build_read_pool(all_pairs, frag_dist);
        return Ok((pool, chrom_a.clone()));
    }

    // Single-region events (DEL, DUP, INV, INS).
    let (chrom, start, end) = event.primary_region().unwrap();
    let (region_start, region_end) =
        extraction_bounds(chrom, start, end, config.flank_bp, extraction_region);
    let pairs = extract::extract_read_pairs(
        &config.bam_path,
        chrom,
        region_start,
        region_end,
        config.min_mapq,
        Some(config.ref_path.as_str()),
    )?;
    let frag_dist = stats::FragmentDist::from_read_pairs(&pairs);
    let pool = extract::build_read_pool(pairs, frag_dist);
    Ok((pool, chrom.to_string()))
}

/// Build a VariantHaplotype for a given event.
fn build_haplotype(
    event: &SimEvent,
    reference: &crate::reference::SharedReference,
    flank: u64,
    dup_model: &str,
    rng: &mut StdRng,
) -> Result<VariantHaplotype> {
    match event {
        SimEvent::Deletion {
            chrom,
            del_start,
            del_end,
            ..
        } => VariantHaplotype::from_deletion(reference, chrom, *del_start, *del_end, flank),
        SimEvent::Duplication {
            chrom,
            dup_start,
            dup_end,
            ..
        } => {
            if dup_model == "junction" {
                VariantHaplotype::from_duplication(reference, chrom, *dup_start, *dup_end, flank)
            } else {
                VariantHaplotype::from_tandem_duplication(
                    reference, chrom, *dup_start, *dup_end, flank,
                )
            }
        }
        SimEvent::Inversion {
            chrom,
            inv_start,
            inv_end,
            ..
        } => VariantHaplotype::from_inversion(reference, chrom, *inv_start, *inv_end, flank),
        SimEvent::Insertion {
            chrom,
            pos,
            ins_seq,
            ins_len,
            ..
        } => {
            // Generate insertion sequence if not provided.
            let seq: Vec<u8> = if let Some(s) = ins_seq {
                s.clone()
            } else {
                let bases = [b'A', b'C', b'G', b'T'];
                (0..*ins_len).map(|_| bases[rng.gen_range(0..4)]).collect()
            };
            VariantHaplotype::from_insertion(reference, chrom, *pos, &seq, flank)
        }
        SimEvent::Fusion {
            chrom_a,
            bp_a,
            chrom_b,
            bp_b,
            inverted,
            ..
        } => VariantHaplotype::from_fusion(
            reference, chrom_a, *bp_a, chrom_b, *bp_b, flank, *inverted,
        ),
        SimEvent::SmallVariant {
            chrom,
            pos,
            ref_allele,
            alt_allele,
            ..
        } => VariantHaplotype::from_small_variant(
            reference, chrom, *pos, ref_allele, alt_allele, flank,
        ),
    }
}

/// Return BED regions (chrom, start, end) covering an event ± flank.
///
/// Used for events.bed (merge step) and README. Coordinates are 0-based half-open.
fn event_extraction_regions(event: &SimEvent, flank: u64) -> Vec<(String, u64, u64)> {
    match event {
        SimEvent::Deletion {
            chrom,
            del_start,
            del_end,
            ..
        } => vec![(
            chrom.clone(),
            del_start.saturating_sub(flank),
            del_end.saturating_add(flank),
        )],
        SimEvent::Duplication {
            chrom,
            dup_start,
            dup_end,
            ..
        } => vec![(
            chrom.clone(),
            dup_start.saturating_sub(flank),
            dup_end.saturating_add(flank),
        )],
        SimEvent::Inversion {
            chrom,
            inv_start,
            inv_end,
            ..
        } => vec![(
            chrom.clone(),
            inv_start.saturating_sub(flank),
            inv_end.saturating_add(flank),
        )],
        SimEvent::Insertion { chrom, pos, .. } => vec![(
            chrom.clone(),
            pos.saturating_sub(flank),
            pos.saturating_add(flank),
        )],
        SimEvent::SmallVariant {
            chrom,
            pos,
            ref_allele,
            ..
        } => {
            let end = pos.saturating_add(ref_allele.len() as u64);
            vec![(
                chrom.clone(),
                pos.saturating_sub(flank),
                end.saturating_add(flank),
            )]
        }
        SimEvent::Fusion {
            chrom_a,
            bp_a,
            chrom_b,
            bp_b,
            ..
        } => vec![
            (
                chrom_a.clone(),
                bp_a.saturating_sub(flank),
                bp_a.saturating_add(flank),
            ),
            (
                chrom_b.clone(),
                bp_b.saturating_sub(flank),
                bp_b.saturating_add(flank),
            ),
        ],
    }
}

/// Write events.bed: one line per extraction region (event ± flank), 0-based half-open.
///
/// Used by merge.sh to identify which reads in the original BAM to replace.
fn write_event_bed(output_dir: &str, events: &[SimEvent], flank: u64) -> Result<()> {
    use std::io::Write as IoWrite;
    let bed_path = Path::new(output_dir).join("events.bed");
    let mut f = std::fs::File::create(&bed_path)?;
    for event in events {
        for (chrom, start, end) in event_extraction_regions(event, flank) {
            writeln!(f, "{}\t{}\t{}", chrom, start, end)?;
        }
    }
    Ok(())
}

/// Write merge.sh: merges sim.bam with the original BAM, replacing event-region reads.
///
/// After running align.sh to produce sim.bam, run merge.sh to produce merged.bam,
/// which is the original BAM with the spiked reads substituted in the event regions.
fn write_merge_script(
    output_dir: &str,
    original_bam: &str,
    ref_path: &str,
    threads: usize,
    samtools: &str,
) -> Result<()> {
    let script_path = Path::new(output_dir).join("merge.sh");

    let script = format!(
        r#"#!/bin/bash
set -euo pipefail
# Merge sim.bam (spiked reads) into the original BAM.
#
# Reads in the event regions (events.bed ± flank) are replaced by the spiked
# reads from sim.bam. All other reads are kept from the original BAM.
#
# Usage: bash merge.sh [ORIGINAL_BAM] [REFERENCE_FASTA] [THREADS]
#
# REFERENCE_FASTA is required when ORIGINAL_BAM is a CRAM file.
# Requires: samtools (>= 1.13 for -U flag support)
ORIGINAL="${{1:-{original_bam}}}"
REF="${{2:-{ref_path}}}"
THREADS="${{3:-{threads}}}"
SAMTOOLS="{samtools}"
DIR="$(cd "$(dirname "$0")" && pwd)"

if [ ! -f "$DIR/sim.bam" ]; then
    echo "Error: $DIR/sim.bam not found. Run align.sh first." >&2
    exit 1
fi

echo "Extracting reads outside event regions from $ORIGINAL..."
"$SAMTOOLS" view -b -T "$REF" -L "$DIR/events.bed" -U "$DIR/outside.bam" "$ORIGINAL" -o /dev/null

echo "Merging spiked reads with outside-region originals..."
"$SAMTOOLS" merge -f -@ "$THREADS" "$DIR/merged_tmp.bam" "$DIR/sim.bam" "$DIR/outside.bam"
rm -f "$DIR/outside.bam"

echo "Sorting and indexing..."
"$SAMTOOLS" sort -@ "$THREADS" -o "$DIR/merged.bam" "$DIR/merged_tmp.bam"
"$SAMTOOLS" index "$DIR/merged.bam"
rm -f "$DIR/merged_tmp.bam"

TOTAL=$("$SAMTOOLS" view -c "$DIR/merged.bam" 2>/dev/null || echo "?")
echo "Done: $DIR/merged.bam ($TOTAL reads)"
"#
    );

    std::fs::write(&script_path, script)?;

    #[cfg(unix)]
    {
        use std::os::unix::fs::PermissionsExt;
        std::fs::set_permissions(&script_path, std::fs::Permissions::from_mode(0o755))?;
    }

    Ok(())
}

/// Return a short human-readable description of an event.
fn event_label(event: &SimEvent) -> String {
    match event {
        SimEvent::Deletion {
            chrom,
            del_start,
            del_end,
            ..
        } => format!(
            "DEL  {}:{}-{} ({}bp)",
            chrom,
            del_start + 1,
            del_end,
            del_end - del_start
        ),
        SimEvent::Duplication {
            chrom,
            dup_start,
            dup_end,
            ..
        } => format!(
            "DUP  {}:{}-{} ({}bp)",
            chrom,
            dup_start + 1,
            dup_end,
            dup_end - dup_start
        ),
        SimEvent::Inversion {
            chrom,
            inv_start,
            inv_end,
            ..
        } => format!(
            "INV  {}:{}-{} ({}bp)",
            chrom,
            inv_start + 1,
            inv_end,
            inv_end - inv_start
        ),
        SimEvent::Insertion {
            chrom,
            pos,
            ins_len,
            ..
        } => format!("INS  {}:{} ({}bp)", chrom, pos + 1, ins_len),
        SimEvent::SmallVariant {
            chrom,
            pos,
            ref_allele,
            alt_allele,
            ..
        } => format!(
            "SNV  {}:{} {}>{}",
            chrom,
            pos + 1,
            String::from_utf8_lossy(ref_allele),
            String::from_utf8_lossy(alt_allele),
        ),
        SimEvent::Fusion {
            chrom_a,
            bp_a,
            chrom_b,
            bp_b,
            inverted,
            ..
        } => format!(
            "FUSION  {}:{}>>{}:{}{}",
            chrom_a,
            bp_a + 1,
            chrom_b,
            bp_b + 1,
            if *inverted { " (inv)" } else { "" }
        ),
    }
}

/// Write README.md documenting this spike run.
fn write_readme(
    output_dir: &str,
    cmdline: &str,
    input_bam: &str,
    reference: &str,
    events: &[SimEvent],
    event_stats: &[(f64, usize, usize, usize)],
    total_pairs: usize,
    flank: u64,
) -> Result<()> {
    use std::fmt::Write as FmtWrite;
    use std::io::Write as IoWrite;

    let now = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap_or_default()
        .as_secs();
    // Format timestamp as YYYY-MM-DD HH:MM:SS UTC (simple, no external crate).
    let secs_per_day = 86400u64;
    let days_since_epoch = now / secs_per_day;
    let time_of_day = now % secs_per_day;
    let hh = time_of_day / 3600;
    let mm = (time_of_day % 3600) / 60;
    let ss = time_of_day % 60;
    // Compute calendar date from days_since_epoch (1970-01-01 = day 0).
    let (year, month, day) = days_to_ymd(days_since_epoch);
    let timestamp = format!("{:04}-{:02}-{:02} {:02}:{:02}:{:02} UTC", year, month, day, hh, mm, ss);

    let version = env!("CARGO_PKG_VERSION");

    let mut md = String::new();
    writeln!(md, "# spike run log")?;
    writeln!(md)?;
    writeln!(md, "**spike version:** {}  ", version)?;
    writeln!(md, "**Date:** {}  ", timestamp)?;
    writeln!(md)?;
    writeln!(md, "## Command")?;
    writeln!(md)?;
    writeln!(md, "```")?;
    writeln!(md, "{}", cmdline)?;
    writeln!(md, "```")?;
    writeln!(md)?;
    writeln!(md, "## Input")?;
    writeln!(md)?;
    writeln!(md, "| Field | Value |")?;
    writeln!(md, "|-------|-------|")?;
    writeln!(md, "| BAM | `{}` |", input_bam)?;
    writeln!(md, "| Reference | `{}` |", reference)?;
    writeln!(md, "| Flank | {}bp |", flank)?;
    writeln!(md)?;
    writeln!(md, "## Events")?;
    writeln!(md)?;
    writeln!(md, "| # | Event | VAF | Kept reads | Chimeric reads | Suppressed reads |")?;
    writeln!(md, "|---|-------|-----|-----------|----------------|-----------------|")?;
    for (i, event) in events.iter().enumerate() {
        let label = event_label(event);
        let (vaf, kept, chimeric, suppressed) = event_stats.get(i).copied().unwrap_or((0.0, 0, 0, 0));
        writeln!(
            md,
            "| {} | {} | {:.3} | {} | {} | {} |",
            i + 1,
            label,
            vaf,
            kept,
            chimeric,
            suppressed,
        )?;
    }
    writeln!(md)?;
    writeln!(md, "## Output files")?;
    writeln!(md)?;
    writeln!(md, "| File | Description |")?;
    writeln!(md, "|------|-------------|")?;
    writeln!(md, "| `R1.fq.gz`, `R2.fq.gz` | Simulated read pairs (total: {}) |", total_pairs)?;
    writeln!(md, "| `truth.vcf` | Ground-truth VCF of introduced variants |")?;
    writeln!(md, "| `events.bed` | Extraction regions (event ± {}bp flank) used to build the spike-in |", flank)?;
    writeln!(md, "| `align.sh` | Aligns R1/R2 → `sim.bam` (event regions ± {}bp flank) |", flank)?;
    writeln!(md, "| `merge.sh` | Merges `sim.bam` into the original BAM → `merged.bam` (full genome) |")?;
    writeln!(md)?;
    writeln!(md, "## Workflow")?;
    writeln!(md)?;
    writeln!(md, "```bash")?;
    writeln!(md, "# Step 1: align the simulated reads")?;
    writeln!(md, "bash align.sh")?;
    writeln!(md)?;
    writeln!(md, "# Step 2a: use sim.bam directly (event regions ± {}bp flank)", flank)?;
    writeln!(md, "#   → useful for targeted analysis of the introduced variants")?;
    writeln!(md)?;
    writeln!(md, "# Step 2b: produce a full modified BAM (original + spiked reads)")?;
    writeln!(md, "bash merge.sh  # produces merged.bam")?;
    writeln!(md, "```")?;

    let readme_path = Path::new(output_dir).join("README.md");
    let mut f = std::fs::File::create(&readme_path)?;
    f.write_all(md.as_bytes())?;
    Ok(())
}

/// Convert days since Unix epoch (1970-01-01) to (year, month, day).
fn days_to_ymd(mut days: u64) -> (u64, u64, u64) {
    // Gregorian calendar calculation.
    let mut year = 1970u64;
    loop {
        let leap = (year % 4 == 0 && year % 100 != 0) || year % 400 == 0;
        let days_in_year = if leap { 366 } else { 365 };
        if days < days_in_year {
            break;
        }
        days -= days_in_year;
        year += 1;
    }
    let leap = (year % 4 == 0 && year % 100 != 0) || year % 400 == 0;
    let month_days: [u64; 12] = [31, if leap { 29 } else { 28 }, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];
    let mut month = 1u64;
    for &md in &month_days {
        if days < md {
            break;
        }
        days -= md;
        month += 1;
    }
    (year, month, days + 1)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn del(chrom: &str, start: u64, end: u64) -> SimEvent {
        SimEvent::Deletion {
            chrom: chrom.to_string(),
            del_start: start,
            del_end: end,
            gene: "G".to_string(),
            exons: Vec::new(),
            allele_fraction: None,
        }
    }

    fn ins(chrom: &str, pos: u64) -> SimEvent {
        SimEvent::Insertion {
            chrom: chrom.to_string(),
            pos,
            ins_seq: Some(vec![b'A']),
            ins_len: 1,
            gene: "G".to_string(),
            allele_fraction: None,
        }
    }

    fn fusion(chrom_a: &str, bp_a: u64, chrom_b: &str, bp_b: u64) -> SimEvent {
        SimEvent::Fusion {
            chrom_a: chrom_a.to_string(),
            bp_a,
            gene_a: "A".to_string(),
            chrom_b: chrom_b.to_string(),
            bp_b,
            gene_b: "B".to_string(),
            allele_fraction: None,
            inverted: false,
        }
    }

    fn pair(name: &str, ref_start: u64) -> ReadPair {
        ReadPair {
            name: name.to_string(),
            seq1: vec![b'A'; 10],
            qual1: vec![b'!' + 30; 10],
            seq2: vec![b'T'; 10],
            qual2: vec![b'!' + 30; 10],
            ref_start,
            ref_end: ref_start + 20,
            insert_size: 20,
            chrom: "chr1".to_string(),
        }
    }

    #[test]
    fn test_overlap_policy_rejects_range_overlap() {
        let events = vec![del("chr1", 100, 200), del("chr1", 150, 250)];
        assert!(validate_event_overlaps(&events, false).is_err());
    }

    #[test]
    fn test_overlap_policy_allows_touching_boundaries() {
        let events = vec![del("chr1", 100, 200), del("chr1", 200, 300)];
        assert!(validate_event_overlaps(&events, false).is_ok());
    }

    #[test]
    fn test_overlap_policy_rejects_point_event_overlaps() {
        // Insertion point overlaps deletion interval.
        let events = vec![del("chr1", 100, 200), ins("chr1", 150)];
        assert!(validate_event_overlaps(&events, false).is_err());

        // Fusion breakpoint overlaps deletion interval on chr1.
        let events = vec![del("chr1", 140, 160), fusion("chr1", 150, "chr2", 500)];
        assert!(validate_event_overlaps(&events, false).is_err());
    }

    #[test]
    fn test_overlap_policy_allow_flag() {
        let events = vec![del("chr1", 100, 200), del("chr1", 150, 250)];
        assert!(validate_event_overlaps(&events, true).is_ok());
    }

    #[test]
    fn test_dedup_by_name_keeps_last_occurrence() {
        let mut pairs = vec![pair("dup", 10), pair("keep", 50), pair("dup", 90)];
        dedup_by_name(&mut pairs);

        assert_eq!(pairs.len(), 2);
        assert_eq!(pairs[0].name, "keep");
        assert_eq!(pairs[1].name, "dup");
        assert_eq!(pairs[1].ref_start, 90);
    }
}
