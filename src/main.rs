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
use std::collections::BTreeSet;
use std::path::Path;

use exon::AfSpec;
use haplotype::VariantHaplotype;
use types::{FusionJoin, ReadPair, ReadPool, SimConfig, SimEvent};

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

    /// Extra read extraction region (e.g. "chr19:11080000-11140000"), on top of
    /// event ± flank -- it does not replace the event window. Use this to ensure
    /// the output BAM covers the full gene/region of interest. For an event on
    /// the same chromosome, the region and the event ± flank are merged into one
    /// query when they overlap or touch, and kept as two queries when they do
    /// not, so a distant fusion partner costs one extra event-sized window
    /// rather than every read in between. A region on another chromosome than
    /// the event is ignored for that event.
    #[arg(long)]
    region: Option<String>,

    /// Flanking region (bp) to include around events. Always defines the event
    /// window; when --region is also set, the region adds another window
    /// alongside it (see --region) rather than replacing this one.
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

/// Reference flank on each side of an event in its variant haplotype. Must be
/// >= the longest expected fragment so reads near the haplotype edges form
/// complete pairs.
const HAP_FLANK: u64 = 2000;

/// Check the input BAM's mean read length against spike's max supported
/// fragment length (L5). spike's tiling and depth-copy generation both
/// sample a fragment length in `[read_length, MAX_FRAGMENT_LEN]`; if the
/// read length itself exceeds that max, no valid fragment exists and
/// `FragmentDist::sample_in_range` panics deep inside per-event processing.
/// Reject it here instead, before any work starts, with a message that
/// explains why: spike simulates fixed-length paired-end reads and does not
/// support long-read (PacBio/ONT) libraries.
fn validate_read_length(read_length: usize) -> Result<()> {
    let max = crate::stats::MAX_FRAGMENT_LEN;
    if read_length as i64 > max {
        bail!(
            "input BAM's mean read length ({read_length}bp) exceeds spike's max \
             supported fragment length ({max}bp); spike simulates fixed-length \
             paired-end reads and does not support long-read (PacBio/ONT) libraries",
        );
    }
    Ok(())
}

/// Check `--flank`: originals are only suppressed inside the extracted
/// window (event ± flank), but synthetic reads cover event ± HAP_FLANK.
fn validate_flank(flank: u64) -> Result<()> {
    if flank < HAP_FLANK {
        bail!(
            "--flank {} is too small: it must be at least {} so every original \
             read replaced by synthetic reads is extracted",
            flank,
            HAP_FLANK
        );
    }
    Ok(())
}

/// Parsed extraction region from --region flag.
struct ExtractionRegion {
    chrom: String,
    start: u64,
    end: u64,
}

/// Compute read extraction windows for an event, in ascending order.
///
/// If `--region` is set and the event is on the same chromosome, the region and
/// event ± flank are merged into one window only when they overlap or touch;
/// otherwise both are returned and the gap between them is never read. Taking
/// the min/max instead turns a region far from the event into one span over
/// everything in between (M8).
/// Otherwise, fall back to event ± flank.
fn extraction_bounds(
    event_chrom: &str,
    event_start: u64,
    event_end: u64,
    flank: u64,
    region: &Option<ExtractionRegion>,
) -> Vec<(u64, u64)> {
    let event_window = (
        event_start.saturating_sub(flank),
        event_end.saturating_add(flank),
    );
    if let Some(r) = region {
        if r.chrom == event_chrom {
            // Half-open windows that touch (r.end == event_window.0) cover a
            // gapless span, so they merge too -- one window keeps each pair in
            // a single extraction.
            if r.start <= event_window.1 && event_window.0 <= r.end {
                return vec![(r.start.min(event_window.0), r.end.max(event_window.1))];
            }
            if r.start < event_window.0 {
                return vec![(r.start, r.end), event_window];
            }
            return vec![event_window, (r.start, r.end)];
        }
    }
    // Fallback: event ± flank.
    vec![event_window]
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
    validate_flank(args.flank)?;

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
    validate_read_length(read_length)?;

    // The simulated read group reuses the original sample, so merged.bam stays
    // single-sample. Read here, not in align.sh: align.sh never sees the BAM.
    let sample = crate::bam_stats::sample_name(&args.bam, Some(args.reference.as_str()))?
        .unwrap_or_else(|| {
            log::warn!("no @RG SM in {}; tagging simulated reads SM:SIM", args.bam);
            "SIM".to_string()
        });
    log::info!("Simulated read group: @RG ID:sim SM:{}", sample);

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

    let mut event_outputs = Vec::with_capacity(events.len());

    // Per-event stats collected for README/log output.
    struct EventStat {
        vaf: f64,
        kept: usize,
        chimeric: usize,
        suppressed: usize,
    }
    let mut event_stats: Vec<EventStat> = Vec::new();
    // M14: pairs whose stored quality is unusable never reach a pool, so they
    // are in neither kept_originals nor suppressed_names.
    let mut unusable_qual_names: BTreeSet<String> = BTreeSet::new();


    // Process each event using the unified haplotype + tiling approach.
    for (i, event) in events.iter().enumerate() {
        log::info!("Processing event {}/{}: {:?}", i + 1, events.len(), event);
        let vaf = event.allele_fraction().unwrap_or(config.allele_fraction);
        log::info!("  Using VAF={:.3} for this event", vaf);

        // Extract reads and build pool.
        let (pool, _extraction_chrom) = extract_pool_for_event(
            event,
            &config,
            &extraction_region,
            &mut unusable_qual_names,
        )?;

        // Build quality profile and synth generator.
        let quality_profile =
            synth::QualityProfile::from_read_pairs(&pool.pairs, config.read_length);
        let synth_gen = synth::SynthReadGenerator::new(
            quality_profile,
            &shared_ref,
            config.read_length,
            config.indel_error_rate,
        );

        // Build variant haplotype.
        let mut haplotype =
            build_haplotype(event, &shared_ref, HAP_FLANK, &config.dup_model, &mut rng)?;

        // Simulate: suppress reads + tile synthetic reads across haplotype.
        let output = simulate::simulate_event(
            i + 1,
            event,
            &pool,
            &mut haplotype,
            &config,
            &synth_gen,
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

        event_outputs.push(output);
    }

    // Names of the originals spike took out of the BAM; merge.sh removes
    // exactly these.
    let mut replaced_names = simulate::consumed_original_names(&event_outputs);
    // M14: a pair dropped for unusable quality cannot be replaced -- spike has
    // no quality string to write for it -- but it must still be removed.
    // Left in place it would sit inside every simulated event as
    // un-suppressible reference support and dilute the realised VAF; removed,
    // it costs only its own depth, uniformly across the extraction window.
    let dropped_unreplaced = unusable_qual_names
        .iter()
        .filter(|name| !replaced_names.contains(*name))
        .count();
    replaced_names.extend(unusable_qual_names);

    let all_output_pairs = simulate::combine_event_outputs(event_outputs);

    // Write FASTQ.
    let (r1_path, r2_path) = fastq::write_paired_fastq(&all_output_pairs, &args.output)?;

    // Write truth VCF.
    let truth_path = Path::new(&args.output).join("truth.vcf");
    truth::write_truth_vcf(
        &events,
        config.allele_fraction, // default AF for events without per-event override
        &truth_path.to_string_lossy(),
        &args.reference,
        &shared_ref,
        &reference::fasta_contigs(&args.reference)?,
    )?;

    // Write alignment convenience script.
    write_align_script(
        &args.output,
        &args.reference,
        args.threads,
        &args.aligner,
        &args.samtools,
        &sample,
    )?;

    // Write events BED (extraction regions, for inspection).
    write_event_bed(&args.output, &events, args.flank)?;

    // Write the read names merge.sh removes from the original BAM.
    write_replaced_reads(&args.output, &replaced_names)?;
    log::info!(
        "Originals removed from the BAM: {} read pairs (replaced_reads.txt); {} of \
         them were dropped for unusable quality and are not replaced",
        replaced_names.len(),
        dropped_unreplaced,
    );

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
    sample: &str,
) -> Result<()> {
    let script_path = Path::new(output_dir).join("align.sh");

    let align_cmd = match aligner {
        "bwa-mem2" => format!("bwa-mem2 mem -t \"$THREADS\" \\\n    -R '@RG\\tID:sim\\tSM:{sample}\\tPL:ILLUMINA' \\\n    \"$REF\" \\\n    \"$DIR/R1.fq.gz\" \"$DIR/R2.fq.gz\" \\\n    2>\"$DIR/align.log\""),
        "minimap2" => format!("minimap2 -a -x sr -t \"$THREADS\" \\\n    -R '@RG\\tID:sim\\tSM:{sample}\\tPL:ILLUMINA' \\\n    \"$REF\" \\\n    \"$DIR/R1.fq.gz\" \"$DIR/R2.fq.gz\" \\\n    2>\"$DIR/align.log\""),
        "bowtie2" => format!("bowtie2 -x \"$REF\" \\\n    -1 \"$DIR/R1.fq.gz\" -2 \"$DIR/R2.fq.gz\" \\\n    -p \"$THREADS\" \\\n    --rg-id sim --rg SM:{sample} --rg PL:ILLUMINA \\\n    2>\"$DIR/align.log\""),
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
# `grep -c` exits 1 when it counts nothing but still prints "0", so the
# `|| echo "0"` this used to end with produced the two-line string "0\n0" on
# every sim.bam with no supplementary alignments -- which is most of them.
SA_COUNT=$("$SAMTOOLS" view "$DIR/sim.bam" | grep -c "SA:Z:" || true)
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
    unusable_qual_names: &mut BTreeSet<String>,
) -> Result<(ReadPool, String)> {
    let mut all_pairs: Vec<ReadPair> = Vec::new();
    let pool_chrom;

    if let SimEvent::Fusion {
        chrom_a,
        bp_a,
        chrom_b,
        bp_b,
        ..
    } = event
    {
        // Fusion: extract from both gene regions and merge pools.
        let windows_a =
            extraction_bounds(chrom_a, *bp_a, *bp_a, config.flank_bp, extraction_region);
        let windows_b =
            extraction_bounds(chrom_b, *bp_b, *bp_b, config.flank_bp, extraction_region);

        extract_windows(
            config,
            chrom_a,
            &windows_a,
            &mut all_pairs,
            unusable_qual_names,
        )?;
        extract_windows(
            config,
            chrom_b,
            &windows_b,
            &mut all_pairs,
            unusable_qual_names,
        )?;
        pool_chrom = chrom_a.clone();
    } else {
        // Single-region events (DEL, DUP, INV, INS).
        let (chrom, start, end) = event.primary_region().unwrap();
        let windows = extraction_bounds(chrom, start, end, config.flank_bp, extraction_region);
        extract_windows(config, chrom, &windows, &mut all_pairs, unusable_qual_names)?;
        pool_chrom = chrom.to_string();
    }

    // Windows can share reads, and the same fragment must not enter the pool
    // -- or the fragment distribution -- twice.
    extract::dedup_pairs_by_name(&mut all_pairs);
    let frag_dist = stats::FragmentDist::from_read_pairs(&all_pairs);
    let pool = extract::build_read_pool(all_pairs, frag_dist);
    Ok((pool, pool_chrom))
}

/// Extract read pairs from every window of one event side into `pairs`.
fn extract_windows(
    config: &SimConfig,
    chrom: &str,
    windows: &[(u64, u64)],
    pairs: &mut Vec<ReadPair>,
    unusable_qual_names: &mut BTreeSet<String>,
) -> Result<()> {
    for &(start, end) in windows {
        let extracted = extract::extract_read_pairs(
            &config.bam_path,
            chrom,
            start,
            end,
            config.min_mapq,
            Some(config.ref_path.as_str()),
        )?;
        unusable_qual_names.extend(extracted.unusable_qual_names);
        pairs.extend(extracted.pairs);
    }
    Ok(())
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
            join,
            ..
        } => VariantHaplotype::from_fusion(
            reference, chrom_a, *bp_a, chrom_b, *bp_b, flank, *join,
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
/// Used for events.bed (written for inspection) and README. Coordinates are 0-based half-open.
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

/// Write events.bed: one line per event ± flank window, 0-based half-open
/// (a fusion writes one line per breakpoint).
///
/// For inspection only, and incomplete as a record of what was extracted:
/// it never reflects --region, and when --region adds an extra window that
/// is actually extracted (see extraction_bounds), that window is not written
/// here either. merge.sh does not read this file -- it selects originals to
/// replace by read name (replaced_reads.txt), not by region.
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

/// Write replaced_reads.txt: the names of the originals spike took out of the BAM.
///
/// merge.sh feeds this to `samtools view -N` to drop exactly these records from
/// the original, one name per line, sorted.
fn write_replaced_reads(output_dir: &str, names: &BTreeSet<String>) -> Result<()> {
    use std::io::Write as IoWrite;
    let path = Path::new(output_dir).join("replaced_reads.txt");
    let mut f = std::io::BufWriter::new(std::fs::File::create(&path)?);
    for name in names {
        writeln!(f, "{}", name)?;
    }
    f.flush()?;
    Ok(())
}

/// Write merge.sh: merges sim.bam with the original BAM, replacing the reads spike took.
///
/// After running align.sh to produce sim.bam, run merge.sh to produce merged.bam,
/// which is the original BAM with the spiked reads substituted for the originals
/// listed in replaced_reads.txt.
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
# The originals spike extracted are listed by read name in replaced_reads.txt;
# sim.bam holds their replacements. Removing them by name (rather than by event
# region) keeps every record spike did not take -- duplicates, non-proper pairs,
# low-MAPQ or orphaned mates -- and also removes the out-of-region mates spike
# did take, so no record is lost and none appears twice.
#
# Usage: bash merge.sh [ORIGINAL_BAM] [REFERENCE_FASTA] [THREADS]
#
# ORIGINAL_BAM must be the exact BAM spike was run on: replaced_reads.txt names
# the reads spike extracted from it, and -N below only removes names it finds.
# A different BAM shares essentially no read names, so this script checks for
# that and aborts rather than silently merging sim.bam onto full, un-thinned
# original depth.
#
# REFERENCE_FASTA is required when ORIGINAL_BAM is a CRAM file.
# Requires: samtools (>= 1.13 for -N and -U flag support)
ORIGINAL="${{1:-{original_bam}}}"
REF="${{2:-{ref_path}}}"
THREADS="${{3:-{threads}}}"
SAMTOOLS="{samtools}"
DIR="$(cd "$(dirname "$0")" && pwd)"

if [ ! -f "$DIR/sim.bam" ]; then
    echo "Error: $DIR/sim.bam not found. Run align.sh first." >&2
    exit 1
fi

if [ ! -f "$DIR/replaced_reads.txt" ]; then
    echo "Error: $DIR/replaced_reads.txt not found. Re-run spike." >&2
    exit 1
fi

echo "Removing the replaced originals from $ORIGINAL..."
"$SAMTOOLS" view -b -T "$REF" -N "$DIR/replaced_reads.txt" -U "$DIR/outside.bam" "$ORIGINAL" -o "$DIR/removed.bam"

# Sanity check: every name in replaced_reads.txt came from a read spike found
# in ORIGINAL, so a correct ORIGINAL always yields at least one matching record
# per name (usually two, for a pair, plus any secondary/supplementary records),
# so ACTUAL should land well above NAMES. If ACTUAL comes in below NAMES --
# even by a single record -- that floor is broken and ORIGINAL is most likely
# not the BAM spike was run on -- samtools does not error on read names it
# cannot find, so without this check the merge below would silently combine
# sim.bam with (near) full original depth instead of the thinned original.
NAMES=$(wc -l < "$DIR/replaced_reads.txt")
ACTUAL=$("$SAMTOOLS" view -c "$DIR/removed.bam")
rm -f "$DIR/removed.bam"
if [ "$NAMES" -gt 0 ] && [ "$ACTUAL" -lt "$NAMES" ]; then
    echo "Error: $ORIGINAL yielded only $ACTUAL matching record(s) for the $NAMES replaced read names listed in replaced_reads.txt." >&2
    echo "ORIGINAL must be the exact BAM spike was run on -- point merge.sh at a" >&2
    echo "different BAM only if it is that same BAM (e.g. moved to a new path)." >&2
    rm -f "$DIR/outside.bam"
    exit 1
fi

echo "Merging spiked reads with the untouched originals..."
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
            join,
            ..
        } => format!(
            "FUSION  {}:{}>>{}:{}{}",
            chrom_a,
            bp_a + 1,
            chrom_b,
            bp_b + 1,
            match join {
                FusionJoin::Forward => "",
                FusionJoin::LeftLeft => " (left-left, B reversed)",
                FusionJoin::RightRight => " (right-right, A reversed)",
            }
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
    writeln!(md, "| `replaced_reads.txt` | Names of the originals spike took out of the BAM, plus any pair dropped for unusable quality; `merge.sh` removes exactly these |")?;
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
            join: FusionJoin::Forward,
        }
    }

    fn region(chrom: &str, start: u64, end: u64) -> Option<ExtractionRegion> {
        Some(ExtractionRegion {
            chrom: chrom.to_string(),
            start,
            end,
        })
    }

    #[test]
    fn test_validate_read_length_rejects_long_read_length() {
        // L5: a long-read (PacBio/ONT) BAM's mean read length can exceed
        // spike's max supported fragment length (1500bp), which used to
        // reach an unguarded `clamp` panic deep in tiling. Reject it here,
        // early and with a clear message, instead.
        let err = validate_read_length(2000).expect_err("2000bp should be rejected");
        let msg = err.to_string();
        assert!(msg.contains("2000"), "message should cite the read length: {msg}");
        assert!(msg.contains("1500"), "message should cite the max: {msg}");
    }

    #[test]
    fn test_validate_read_length_accepts_normal_illumina_length() {
        assert!(validate_read_length(151).is_ok());
    }

    #[test]
    fn test_extraction_bounds_keeps_distant_region_and_event_apart() {
        // M8: --region chr20:30490000-30510000 with a fusion partner at 35 Mb
        // used to take the min/max, extracting every read in the 4.5 Mb
        // between them.
        let r = region("chr20", 30_489_999, 30_510_000);

        let windows = extraction_bounds("chr20", 35_000_000, 35_000_000, 10_000, &r);

        assert_eq!(
            windows,
            vec![(30_489_999, 30_510_000), (34_990_000, 35_010_000)],
        );
    }

    #[test]
    fn test_extraction_bounds_merges_overlapping_region_and_event() {
        // Event window runs past both ends of the region: one window covering
        // the union, so no pair is extracted twice.
        let r = region("chr20", 30_489_999, 30_510_000);

        let windows = extraction_bounds("chr20", 30_495_000, 30_520_000, 10_000, &r);

        assert_eq!(windows, vec![(30_485_000, 30_530_000)]);
    }

    #[test]
    fn test_extraction_bounds_merges_windows_that_only_touch() {
        // Region ends exactly where the event window starts. Half-open windows
        // that touch cover a gapless span, so merging them extracts the same
        // reads while keeping each pair in one window.
        let r = region("chr20", 30_000_000, 30_490_000);

        let windows = extraction_bounds("chr20", 30_500_000, 30_500_000, 10_000, &r);

        assert_eq!(windows, vec![(30_000_000, 30_510_000)]);
    }

    #[test]
    fn test_extraction_bounds_splits_windows_one_bp_apart() {
        // One base of gap is still a gap: two queries, not one span.
        let r = region("chr20", 30_000_000, 30_489_999);

        let windows = extraction_bounds("chr20", 30_500_000, 30_500_000, 10_000, &r);

        assert_eq!(
            windows,
            vec![(30_000_000, 30_489_999), (30_490_000, 30_510_000)],
        );
    }

    #[test]
    fn test_extraction_bounds_ignores_region_on_another_chromosome() {
        let r = region("chr19", 11_080_000, 11_140_000);

        let windows = extraction_bounds("chr20", 30_500_000, 30_500_000, 10_000, &r);

        assert_eq!(windows, vec![(30_490_000, 30_510_000)]);
    }

    fn scratch_dir(name: &str) -> std::path::PathBuf {
        let dir = std::env::temp_dir().join(format!("spike_test_{}_{}", name, std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        dir
    }

    #[test]
    fn test_replaced_reads_file_lists_each_consumed_name_once() {
        let dir = scratch_dir("replaced_reads");
        let names: BTreeSet<String> = ["b", "a", "b"].iter().map(|n| n.to_string()).collect();

        write_replaced_reads(dir.to_str().unwrap(), &names).unwrap();

        let text = std::fs::read_to_string(dir.join("replaced_reads.txt")).unwrap();
        assert_eq!(text, "a\nb\n");
    }

    #[test]
    fn test_merge_script_removes_originals_by_read_name() {
        // Removing by BED region loses in-region records spike never extracted
        // (duplicates, non-proper pairs, low-MAPQ or orphaned mates) and keeps
        // out-of-region mates that spike did extract, so those end up twice.
        let dir = scratch_dir("merge_script");

        write_merge_script(dir.to_str().unwrap(), "orig.bam", "ref.fa", 4, "samtools").unwrap();

        let script = std::fs::read_to_string(dir.join("merge.sh")).unwrap();
        assert!(
            script.contains("-N \"$DIR/replaced_reads.txt\""),
            "merge.sh must select the originals to drop by read name:\n{}",
            script
        );
        assert!(
            !script.contains("-L \"$DIR/events.bed\""),
            "merge.sh must not drop originals by BED region:\n{}",
            script
        );
    }

    /// A stub "samtools" for exercising merge.sh's own shell logic (the
    /// shortfall guard, control flow, exit codes) without real BAM files or
    /// real samtools. It only understands the exact invocations
    /// `write_merge_script` emits:
    /// - `view -c PATH`: for a path ending in `removed.bam`, exits 1 if that
    ///   file does not exist (pins merge.sh materialising it before counting),
    ///   else echoes `$SAMTOOLS_FAKE_ACTUAL` (default 0); any other path
    ///   echoes `0`.
    /// - `view ... -o OUT -U UN ... IN`: touches OUT and UN (no real filtering).
    /// - `merge -f -@ N OUT IN...` / `sort -@ N -o OUT IN`: touches OUT.
    /// - `index FILE`: no-op.
    ///
    /// Everything exits 0 so `set -euo pipefail` never trips on the stub itself
    /// -- only the guard logic under test can fail the script.
    fn write_stub_samtools(dir: &std::path::Path) -> std::path::PathBuf {
        let path = dir.join("fake_samtools.sh");
        let script = r#"#!/bin/bash
set -u
case "$1" in
  view)
    if [ "$2" = "-c" ]; then
      case "$3" in
        */removed.bam) [ -f "$3" ] || exit 1; echo "${SAMTOOLS_FAKE_ACTUAL:-0}" ;;
        *) echo "0" ;;
      esac
      exit 0
    fi
    outfile=""
    unfile=""
    shift
    while [ $# -gt 0 ]; do
      case "$1" in
        -o) outfile="$2"; shift 2 ;;
        -U) unfile="$2"; shift 2 ;;
        -T|-N) shift 2 ;;
        -b) shift ;;
        *) shift ;;
      esac
    done
    [ -n "$outfile" ] && : > "$outfile"
    [ -n "$unfile" ] && : > "$unfile"
    exit 0
    ;;
  merge)
    shift
    out=""
    while [ $# -gt 0 ]; do
      case "$1" in
        -f) shift ;;
        -@) shift 2 ;;
        *) out="$1"; break ;;
      esac
    done
    [ -n "$out" ] && : > "$out"
    exit 0
    ;;
  sort)
    shift
    out=""
    while [ $# -gt 0 ]; do
      case "$1" in
        -@) shift 2 ;;
        -o) out="$2"; shift 2 ;;
        *) shift ;;
      esac
    done
    [ -n "$out" ] && : > "$out"
    exit 0
    ;;
  index)
    exit 0
    ;;
  *)
    exit 0
    ;;
esac
"#;
        std::fs::write(&path, script).unwrap();
        #[cfg(unix)]
        {
            use std::os::unix::fs::PermissionsExt;
            std::fs::set_permissions(&path, std::fs::Permissions::from_mode(0o755)).unwrap();
        }
        path
    }

    #[test]
    fn test_merge_script_aborts_when_original_bam_does_not_match_replaced_reads() {
        // `bash merge.sh /other.bam` lets a caller override ORIGINAL_BAM. If that
        // BAM is not the one spike ran on, none of its read names are in
        // replaced_reads.txt, so `-N` matches nothing and the unguarded script
        // would silently merge sim.bam onto full, un-thinned original depth.
        let dir = scratch_dir("merge_guard_abort");
        let samtools = write_stub_samtools(&dir);

        write_merge_script(
            dir.to_str().unwrap(),
            "orig.bam",
            "ref.fa",
            4,
            samtools.to_str().unwrap(),
        )
        .unwrap();
        std::fs::write(dir.join("sim.bam"), b"").unwrap();
        std::fs::write(dir.join("replaced_reads.txt"), "r1\nr2\nr3\nr4\nr5\n").unwrap();

        let output = std::process::Command::new("bash")
            .arg(dir.join("merge.sh"))
            .env("SAMTOOLS_FAKE_ACTUAL", "0") // wrong-BAM case: none of the 5 names matched
            .output()
            .unwrap();

        assert!(
            !output.status.success(),
            "merge.sh must abort when far fewer records matched by name than \
             replaced_reads.txt lists (a wrong ORIGINAL_BAM):\nstdout: {}\nstderr: {}",
            String::from_utf8_lossy(&output.stdout),
            String::from_utf8_lossy(&output.stderr)
        );
        let stderr = String::from_utf8_lossy(&output.stderr);
        assert!(
            stderr.contains("ORIGINAL"),
            "abort message should point at ORIGINAL_BAM as the likely cause:\n{}",
            stderr
        );
        assert!(
            !dir.join("merged.bam").exists(),
            "merge.sh must not produce merged.bam once the shortfall check fails"
        );
    }

    #[test]
    fn test_merge_script_proceeds_when_original_bam_matches_replaced_reads() {
        let dir = scratch_dir("merge_guard_ok");
        let samtools = write_stub_samtools(&dir);

        write_merge_script(
            dir.to_str().unwrap(),
            "orig.bam",
            "ref.fa",
            4,
            samtools.to_str().unwrap(),
        )
        .unwrap();
        std::fs::write(dir.join("sim.bam"), b"").unwrap();
        std::fs::write(dir.join("replaced_reads.txt"), "r1\nr2\nr3\nr4\nr5\n").unwrap();

        let output = std::process::Command::new("bash")
            .arg(dir.join("merge.sh"))
            .env("SAMTOOLS_FAKE_ACTUAL", "10") // healthy case: ~2 records per name
            .output()
            .unwrap();

        assert!(
            output.status.success(),
            "merge.sh must not false-positive when the expected records were found:\nstdout: {}\nstderr: {}",
            String::from_utf8_lossy(&output.stdout),
            String::from_utf8_lossy(&output.stderr)
        );
        assert!(
            dir.join("merged.bam").exists(),
            "merge.sh should have produced merged.bam"
        );
    }

    #[test]
    fn test_align_script_tags_reads_with_the_bam_sample_name() {
        for aligner in ["bwa-mem2", "minimap2", "bowtie2"] {
            let dir = scratch_dir(&format!("align_script_{}", aligner));

            write_align_script(dir.to_str().unwrap(), "ref.fa", 4, aligner, "samtools", "HG002")
                .unwrap();

            let script = std::fs::read_to_string(dir.join("align.sh")).unwrap();
            assert!(
                script.contains("SM:HG002"),
                "{} align.sh must tag reads with the original sample:\n{}",
                aligner,
                script
            );
            assert!(
                !script.contains("SM:SIM"),
                "{} align.sh must not invent a second sample:\n{}",
                aligner,
                script
            );
        }
    }

    fn validate_pipeline_script() -> std::path::PathBuf {
        std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("scripts/validate_pipeline.sh")
    }

    #[test]
    fn test_validate_pipeline_counts_a_header_only_vcf_as_one_zero() {
        // `grep -c` exits 1 when it counts nothing, so the harness's old
        // `$(grep -vc '^#' f || echo 0)` yielded the two-line string "0\n0",
        // and every `[[ $n -lt 5 ]]` guard built on it errored out instead of
        // stopping the run.
        let dir = scratch_dir("validate_pipeline_count");
        let vcf = dir.join("headers_only.vcf");
        std::fs::write(&vcf, "##fileformat=VCFv4.2\n#CHROM\tPOS\n").unwrap();

        let output = std::process::Command::new("bash")
            .arg("-c")
            .arg(r#"script="$1"; vcf="$2"; shift 2; source "$script"; count_records "$vcf""#)
            .arg("_")
            .arg(validate_pipeline_script())
            .arg(&vcf)
            .output()
            .unwrap();

        assert!(
            output.status.success(),
            "sourcing validate_pipeline.sh must define its helpers without running \
             the pipeline:\nstdout: {}\nstderr: {}",
            String::from_utf8_lossy(&output.stdout),
            String::from_utf8_lossy(&output.stderr)
        );
        assert_eq!(
            String::from_utf8_lossy(&output.stdout),
            "0",
            "count_records must return a single integer for a header-only VCF"
        );
    }

    #[test]
    fn test_validate_pipeline_verdict_needs_the_spike_in_to_beat_the_background() {
        // The harness's verdict gate. The truth DELs are common HG002
        // variants, so the background sample carries several of them and a
        // caller recovers those with no spike-in at all: a "TP > 0" gate
        // passes on the background alone, and the run then prints
        // "VALIDATION PASSED" having validated nothing. The gate compares
        // against the background control instead, and is on by default.
        let output = std::process::Command::new("bash")
            .arg("-c")
            .arg(
                r#"script="$1"; shift; source "$script"
echo "min_gain=$MIN_GAIN"
for case in "4 3 1" "3 3 1" "0 0 1" "2 3 1" "4 NONE 1"; do
    set -- $case
    control="$2"
    if [ "$control" = NONE ]; then control=""; fi
    rc=0
    beats_background_control "$1" "$control" "$3" || rc=$?
    echo "tp=$1 control=$2 gain=$3 rc=$rc"
done"#,
            )
            .arg("_")
            .arg(validate_pipeline_script())
            .output()
            .unwrap();

        let stdout = String::from_utf8_lossy(&output.stdout);
        let stderr = String::from_utf8_lossy(&output.stderr);
        assert!(
            output.status.success(),
            "sourcing validate_pipeline.sh must define beats_background_control:\nstdout: {}\nstderr: {}",
            stdout,
            stderr
        );
        // Spike attribution must be the default, not an opt-in flag.
        assert!(
            stdout.contains("min_gain=1"),
            "the attribution gate must be on by default (MIN_GAIN=1):\n{}",
            stdout
        );
        // Recovering more than the background does: that is the spike-in.
        assert!(
            stdout.contains("tp=4 control=3 gain=1 rc=0"),
            "beating the control by 1 must pass:\n{}",
            stdout
        );
        // Matching the background, at any level, is not a spike-in result.
        for line in [
            "tp=3 control=3 gain=1 rc=1",
            "tp=0 control=0 gain=1 rc=1",
            "tp=2 control=3 gain=1 rc=1",
        ] {
            assert!(
                stdout.contains(line),
                "a run that does not beat the background control must fail the \
                 gate ({}):\n{}",
                line,
                stdout
            );
        }
        // No control result at all is a failure, never a pass.
        assert!(
            stdout.contains("tp=4 control=NONE gain=1 rc=2"),
            "a missing control must be reported as uncomparable, not as a pass:\n{}",
            stdout
        );
    }

    #[test]
    fn test_validate_pipeline_verdict_fails_when_the_highest_vaf_has_no_truvari() {
        // The test above pins the gate's truth table; this one pins that the
        // verdict actually consults it. Step 8 used to evaluate the gate only
        // *inside* `if [[ -f summary.json ]]`, so an outdir whose highest VAF
        // has no truvari output at all -- what `--skip-to 8` over an
        // incomplete run looks like -- printed a row of N/A and then
        // "VALIDATION PASSED". Showing that needs no caller, benchmarker or
        // aligner: fabricate the outdir and drive steps 8 and 9 over it.
        let dir = scratch_dir("validate_pipeline_wiring");

        let truth_vcf = "##fileformat=VCFv4.2\n\
                         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
                         chr20\t1000\tsim_del_1\tA\t<DEL>\t.\tPASS\tSVTYPE=DEL\n\
                         chr20\t2000\tsim_del_2\tA\t<DEL>\t.\tPASS\tSVTYPE=DEL\n";
        let summary = |tp: u32| {
            format!(
                "{{\"TP-base\": {}, \"FP\": 1, \"FN\": 4, \
                 \"recall\": 0.5, \"precision\": 0.8, \"f1\": 0.6154}}",
                tp
            )
        };

        // Two fabricated outdirs, alike but for the highest VAF's truvari
        // output: "with" has it and beats the control by 1, "without" has
        // none. Everything else -- the lower VAF, the control -- is identical,
        // so only the missing summary can change the verdict.
        for case in ["with", "without"] {
            let out = dir.join(case);
            for vaf in ["0.5", "0.25"] {
                let spike_out = out.join(format!("spike_vaf_{}", vaf));
                std::fs::create_dir_all(spike_out.join("truvari")).unwrap();
                std::fs::write(spike_out.join("truth.vcf"), truth_vcf).unwrap();
                if case == "without" && vaf == "0.5" {
                    std::fs::remove_dir_all(spike_out.join("truvari")).unwrap();
                    continue;
                }
                let tp = if vaf == "0.5" { 4 } else { 3 };
                std::fs::write(spike_out.join("truvari/summary.json"), summary(tp)).unwrap();
            }
            let control = out.join("background_control/truvari");
            std::fs::create_dir_all(&control).unwrap();
            std::fs::write(control.join("summary.json"), summary(3)).unwrap();
        }

        let output = std::process::Command::new("bash")
            .arg("-c")
            .arg(
                r#"script="$1"; root="$2"; shift 2; source "$script"
VAFS=(0.5 0.25)
MIN_GAIN=1
for case in with without; do
    OUTDIR="$root/$case"
    FAILURES=()
    rc=0
    msg=$( { step8_summarize >/dev/null; step9_verdict; } 2>&1 ) || rc=$?
    echo "$case rc=$rc"
    printf '%s\n' "$msg" | sed "s/^/$case /"
done"#,
            )
            .arg("_")
            .arg(validate_pipeline_script())
            .arg(&dir)
            .output()
            .unwrap();

        let stdout = String::from_utf8_lossy(&output.stdout);
        let stderr = String::from_utf8_lossy(&output.stderr);
        assert!(
            output.status.success(),
            "driving step8_summarize + step9_verdict over a fabricated outdir \
             must not error out:\nstdout: {}\nstderr: {}",
            stdout,
            stderr
        );
        // A complete run that beats the control still passes.
        assert!(
            stdout.contains("with rc=0"),
            "a highest VAF that beats the control by --min-gain must pass:\n{}",
            stdout
        );
        // The one under review: no truvari output for the highest VAF means
        // nothing was measured, so nothing can be attributed to the spike-in.
        assert!(
            stdout.contains("without rc=1"),
            "a highest VAF with no truvari summary must fail the verdict, \
             whichever steps ran:\n{}",
            stdout
        );
        assert!(
            stdout
                .lines()
                .any(|l| l.starts_with("without ") && l.contains("truvari")),
            "the verdict must say the highest VAF has no truvari result:\n{}",
            stdout
        );
    }

    #[test]
    fn test_validate_pipeline_aborts_when_its_data_files_are_missing() {
        // A harness that cannot fail is worthless: missing inputs must stop the
        // run, not let it continue and report an empty result as a pass.
        let dir = scratch_dir("validate_pipeline_prereq");
        let output = std::process::Command::new("bash")
            .arg(validate_pipeline_script())
            .arg("--giab-dir")
            .arg(dir.join("no_such_giab_dir"))
            .arg("--outdir")
            .arg(dir.join("out"))
            .output()
            .unwrap();

        let stderr = String::from_utf8_lossy(&output.stderr);
        assert!(
            !output.status.success(),
            "validate_pipeline.sh must exit non-zero when its data files are \
             missing:\nstdout: {}\nstderr: {}",
            String::from_utf8_lossy(&output.stdout),
            stderr
        );
        assert!(
            stderr.contains("Prerequisite check failed"),
            "it must say which prerequisite failed, and must accept --giab-dir:\n{}",
            stderr
        );
    }

    #[test]
    fn test_flank_smaller_than_haplotype_flank_is_rejected() {
        // With --flank 500, originals 500-2000 bp from the event are never
        // extracted, so never suppressed, while synthetic reads still cover
        // them: depth there would be 1 + VAF.
        assert!(validate_flank(500).is_err());
        assert!(validate_flank(HAP_FLANK - 1).is_err());
        assert!(validate_flank(HAP_FLANK).is_ok());
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

    /// Every region-bounded CRAM query must go through
    /// [`extract::open_cram_reader_for_region`], which is the only place that
    /// prunes the `.crai` to the region's slices (M15, N3). A reader built
    /// straight from `indexed_reader::Builder` walks every container on the
    /// chromosome instead — and no test on the records a query returns can see
    /// that, because the unpruned path returns exactly the same records, only
    /// slowly. So the check is on where the readers are built.
    #[test]
    fn test_cram_readers_are_built_only_by_the_pruning_opener() {
        let builder = "cram::io::indexed_reader::Builder";

        for (name, src) in [
            ("loh.rs", include_str!("loh.rs")),
            ("validate.rs", include_str!("validate.rs")),
        ] {
            assert_eq!(
                src.matches(builder).count(),
                0,
                "{name} builds an indexed CRAM reader itself; it must call \
                 extract::open_cram_reader_for_region so the index is pruned"
            );
        }

        // The opener itself builds two: one to read the header (the region's
        // name cannot be mapped to a reference id without it) and one with the
        // pruned index. Nothing else in extract.rs may build one.
        assert_eq!(
            include_str!("extract.rs").matches(builder).count(),
            2,
            "extract.rs builds an indexed CRAM reader outside the pruning opener"
        );
    }
}
