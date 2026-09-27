//! spike: Haplotype-based read spike-in simulator for genomic variants.
//!
//! Creates synthetic chimeric reads from real BAM data for any variant type
//! (SVs, fusions, SNPs, indels).

mod bam_stats;
mod census;
mod exon;
mod extract;
mod fastq;
mod haplotype;
mod loh;
mod origin;
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
    ///   --event "fusion:GENEA:exon14:GENEB:exon2"       (fusion: adds one junction, not a balanced translocation; requires --exon-bed)
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
    ///   af=<number>: exact AF, af=het: 0.5, af=hom: 1.0
    #[arg(short, long)]
    event: Vec<String>,

    /// Input VCF file with variant records. Supports DEL, INS, DUP, INV, BND,
    /// and standard SNP/indel records (no SVTYPE, explicit REF/ALT alleles).
    /// Can be combined with --event. At least one of --event or --vcf required.
    #[arg(long)]
    vcf: Option<String>,

    /// Read INFO/AF from --vcf records as the allele fraction to simulate.
    ///
    /// Off by default: in a population VCF (gnomAD, 1000G) AF is the allele
    /// frequency in the population, not the fraction of this sample's reads
    /// that should carry the allele, so using it silently produces a truth
    /// set at the wrong VAF. Without this flag only SIM_VAF and VAF are
    /// read, and a record that carries AF but has no usable SIM_VAF or VAF
    /// falls back to --allele-fraction, with a count of how many did on
    /// stderr.
    #[arg(long)]
    vcf_info_af: bool,

    /// Exon BED file. Required when using gene-based --event specs (e.g. "del:GENE:exon4-exon8").
    #[arg(long)]
    exon_bed: Option<String>,

    /// Target allele fraction, in (0.0, 1.0] -- above 0 and at most 1.
    ///
    /// 0 is rejected, not treated as "plant nothing": an event asked for at
    /// AF 0 would still be written to the truth VCF. NaN is rejected for the
    /// same reason.
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
    ///
    /// Minimum 2000, and a smaller value is rejected: synthetic reads cover
    /// event +/- 2000bp, so a narrower extraction window would leave original
    /// reads the synthetic ones are meant to replace outside it.
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

    /// Allow events on the same chromosome that overlap or come within 7000bp
    /// of each other.
    ///
    /// By default they are rejected to keep event effects independent and truth
    /// interpretation unambiguous. The check is on each event's replacement
    /// footprint, not its span: an event replaces reads -- removes the
    /// originals and tiles synthetic ones -- across its span grown by 3500bp on
    /// each side (2000bp of haplotype flank + 1500bp of fragment). So two spans
    /// must be at least 7000bp apart, or each event's synthetic flank writes
    /// plain reference over the other event's edit and cancels it.
    #[arg(long)]
    allow_overlap: bool,

    /// Simulate an event even when more than half the reads over it are ones
    /// spike cannot edit.
    ///
    /// Reads below --min-mapq, not in a proper pair, or with a mate that fails
    /// a filter stay in the merged BAM as they are, so the event that reaches
    /// the reads is about VAF x (1 - that share). Above one half, spike refuses
    /// by default rather than write a truth record the reads cannot back.
    /// truth.vcf records the share as SIM_RESIST either way.
    #[arg(long)]
    allow_resistant: bool,

    /// Duplication model: "full" (default) builds a full tandem haplotype
    /// with duplicated region appearing twice, producing both junction reads
    /// and correct depth increase from a single tiling pass. "junction" uses
    /// the legacy junction-only haplotype with separate depth copies.
    #[arg(long, default_value = "full")]
    dup_model: String,

    /// Which original reads an event replaces. "clean" (default): only the
    /// donor pool's pairs (both mates at --min-mapq or above, a proper pair,
    /// no duplicate, secondary, supplementary or QC-fail flag) inside the
    /// event's footprint. "origin" (experimental): every primary read at the
    /// event and at its look-alikes, each removed by its chance of having
    /// come from the edited copy, read from its MAPQ and its XA tag. It
    /// needs the aligner's XA tags (bwa-mem and bwa-mem2 write them).
    #[arg(long, default_value = "clean")]
    edit_model: String,
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

/// Check `--allele-fraction`: must be in (0.0, 1.0].
fn validate_allele_fraction(af: f64) -> Result<()> {
    // Negated so NaN (for which both `> 0.0` and `<= 1.0` are false) is
    // rejected rather than silently let through (L8).
    if !(af > 0.0 && af <= 1.0) {
        bail!("allele-fraction must be in (0.0, 1.0]");
    }
    Ok(())
}

/// Check `--edit-model`: "clean" or "origin".
fn validate_edit_model(model: &str) -> Result<()> {
    if model != "clean" && model != "origin" {
        bail!("invalid --edit-model '{}', expected 'clean' or 'origin'", model);
    }
    Ok(())
}

/// Check `--flank`: originals are only suppressed inside the extracted
/// window (event ± flank), but synthetic reads cover event ± HAP_FLANK.
/// Turns an event's `;af=` spec into the fraction it is simulated at. `None`
/// means the event takes the global `--allele-fraction`.
///
/// `af=het` is one copy of two, so exactly 0.5 (CR7(c)). It used to draw
/// from Beta(40,40), which moved the event itself, not just what the reads
/// show: above 0.5 the other copy lost reads too.
fn resolve_af_spec(af_spec: Option<AfSpec>) -> Option<f64> {
    match af_spec {
        Some(AfSpec::Exact(v)) => Some(v),
        Some(AfSpec::Het) => Some(0.5),
        Some(AfSpec::Hom) => Some(1.0),
        None => None, // will use global default
    }
}

/// What fusion mode is, said once per run that has a fusion (CR6). It keeps
/// every original read and adds reads across one junction, so it is not a
/// balanced translocation: nothing on either partner is replaced, and there
/// is no reciprocal junction.
fn fusion_mode_warning(events: &[SimEvent]) -> Option<String> {
    let n = events
        .iter()
        .filter(|e| matches!(e, SimEvent::Fusion { .. }))
        .count();
    if n == 0 {
        return None;
    }
    Some(format!(
        "{} fusion event(s): fusion mode keeps every original read and adds junction \
         reads for one join on top. It is not a balanced translocation: no read on \
         either partner is replaced and there is no reciprocal junction, so the \
         junction evidence is stronger than a copy-neutral sample would show. See \
         CR6 in CLINICAL_SV_DESIGN_NOTES.md.",
        n
    ))
}

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
///
/// `contigs` are the reference's contig names, used to keep a contig name that
/// itself contains ':' whole (L13); an empty slice just splits on ':'.
fn parse_region(s: &str, contigs: &[String]) -> Result<ExtractionRegion> {
    let (chrom, coords) = reference::split_contig(s, contigs)
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
    validate_allele_fraction(args.allele_fraction)?;
    validate_flank(args.flank)?;

    // Create output directory.
    std::fs::create_dir_all(&args.output)?;

    // Contig names and lengths from the reference .fai. Read before the
    // --event and --region specs are parsed because a contig name may itself
    // contain ':' (L13), so the specs cannot be split without them.
    let ref_contigs = reference::fasta_contigs(&args.reference)?;
    let contig_names: Vec<String> = ref_contigs.iter().map(|(name, _)| name.clone()).collect();

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
        .map(|spec| exon::parse_event_spec(spec, &genes, &contig_names))
        .collect::<Result<Vec<_>>>()?;

    // Load events from --vcf if provided.
    let vcf_events: Vec<SimEvent> = if let Some(vcf_path) = &args.vcf {
        vcf_input::load_events_from_vcf(vcf_path, args.vcf_info_af)?
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

    validate_edit_model(&args.edit_model)?;

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
        let r = parse_region(region_str, &contig_names)?;
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
            event.set_allele_fraction(resolve_af_spec(af_spec));
            event
        })
        .collect();

    // Append VCF-sourced events (AF already embedded from VCF INFO).
    events.extend(vcf_events);

    if let Some(warning) = fusion_mode_warning(&events) {
        log::warn!("{}", warning);
    }

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
    // Per event, in the same order: the fraction its tiled fragments actually
    // plant when the additive cap or the two-fragment floor moved it off the
    // request, `None` when the request stands. The truth VCF records it.
    let mut adjusted_afs: Vec<Option<f64>> = Vec::with_capacity(events.len());
    // Per event, in the same order: the share of the reads over it spike
    // could not edit (CR4). The truth VCF records it as SIM_RESIST.
    let mut resistant: Vec<Option<f64>> = Vec::with_capacity(events.len());
    // Per event, in the same order: how far the donor's depth departs from
    // the one depth its fragments are scaled by (CR2): SIM_DEPTH_FOLD.
    let mut depth_folds: Vec<Option<f64>> = Vec::with_capacity(events.len());

    let mut event_stats: Vec<EventStat> = Vec::new();
    // M14: pairs whose stored quality is unusable never reach a pool, so they
    // are in neither kept_originals nor suppressed_names.
    let mut unusable_qual_names: BTreeSet<String> = BTreeSet::new();
    // RF8: one line per event whose reads cannot back its truth record.
    let mut refusals: Vec<String> = Vec::new();

    let edit_origin = args.edit_model == "origin";
    if edit_origin {
        // R8: once per run, on the file. A spot whose MAPQ 0 reads lack XA is
        // normal (more than 5 hits); a file with no XA at all was stripped.
        let first = origin::require_xa(&args.bam, &args.reference)?;
        log::info!(
            "--edit-model origin (experimental): reads are removed by their chance of having \
             come from each event's edited copy, at the event and at its look-alikes. The BAM \
             keeps XA tags (the first is on record {}).",
            first,
        );
    }
    // Under origin: each event's label and site, for the removal log after
    // the one draw over every event.
    let mut origin_sites: Vec<(String, origin::OriginSite)> = Vec::new();

    // Process each event using the unified haplotype + tiling approach.
    let n_events = events.len();
    for (i, event) in events.iter_mut().enumerate() {
        log::info!("Processing event {}/{}: {:?}", i + 1, n_events, event);
        let vaf = event.allele_fraction().unwrap_or(config.allele_fraction);
        log::info!("  Using VAF={:.3} for this event", vaf);

        // Extract reads and build pool.
        let (pool, _extraction_chrom, dropped_unusable_qual) = extract_pool_for_event(
            event,
            &config,
            &extraction_region,
            &mut unusable_qual_names,
        )?;

        // --edit-model origin: every primary read at the event and at its
        // look-alikes. Additive events remove nothing, so they get no site.
        let site = if edit_origin && !simulate::is_additive(event, &config.dup_model) {
            let chrom = event
                .primary_region()
                .map(|(chrom, _, _)| chrom.to_string())
                .expect("an event that removes reads has one region");
            let contig_len = shared_ref
                .chromosome_length(&chrom)
                .ok_or_else(|| anyhow::anyhow!("{} is not in the loaded reference", chrom))?;
            let footprint = origin_footprint(event, contig_len)
                .expect("an event that removes reads has one region");
            let site = origin::gather(
                &config.bam_path,
                &config.ref_path,
                &footprint,
                config.read_length,
                &pool,
            )?;
            log::info!(
                "  origin: {} fragment(s) could have come from {}; look-alike region(s): {}",
                site.fragments().len(),
                footprint,
                if site.lookalikes.is_empty() {
                    "none".to_string()
                } else {
                    site.lookalikes
                        .iter()
                        .map(|s| s.to_string())
                        .collect::<Vec<_>>()
                        .join(", ")
                },
            );
            Some(site)
        } else {
            None
        };

        // CR4: the reads over the event the pool does not hold survive it
        // untouched. Unusable-quality pairs are not among them: merge.sh
        // removes those by name.
        let editable: std::collections::HashSet<String> = pool
            .pairs
            .iter()
            .map(|p| p.name.clone())
            .chain(unusable_qual_names.iter().cloned())
            .chain(site.iter().flat_map(|s| s.removable_names()))
            .collect();
        let event_census = census::count_resistant(
            &config.bam_path,
            &config.ref_path,
            &census::census_spans(event),
            &editable,
        )?;
        log::info!(
            "  Reads over the event spike cannot edit: {} of {} ({:.3})",
            event_census.resistant,
            event_census.counted,
            event_census.fraction(),
        );
        if let Some(warning) = census::warning(&event_label(event), &event_census) {
            log::warn!("{}", warning);
        }
        if let Some(refusal) =
            census::refusal(&event_label(event), &event_census, args.allow_resistant)
        {
            refusals.push(refusal);
        }
        // RF8: once any event is refused nothing is written, so the rest are
        // only counted, not simulated -- the error names every one of them.
        if !refusals.is_empty() {
            continue;
        }

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
        let output = match site.as_ref() {
            Some(site) => simulate::simulate_event_origin(
                i + 1,
                event,
                &pool,
                &mut haplotype,
                &config,
                &synth_gen,
                vaf,
                Some(site),
                &mut rng,
            )?,
            // Keep this call to simulate_event, a pure forwarder that passes
            // None: main is its only production caller, so collapsing this
            // match orphans it and adds a clippy warning this branch forbids.
            None => simulate::simulate_event(
                i + 1,
                event,
                &pool,
                &mut haplotype,
                &config,
                &synth_gen,
                vaf,
                &mut rng,
            )?,
        };

        log::info!(
            "Event {}: {} kept + {} chimeric, {} suppressed",
            i + 1,
            output.kept_originals.len(),
            output.chimeric_pairs.len(),
            output.suppressed_count,
        );

        event_stats.push(EventStat {
            vaf,
            adjusted_vaf: output.adjusted_vaf,
            kept: output.kept_originals.len(),
            chimeric: output.chimeric_pairs.len(),
            suppressed: output.suppressed_count,
            dropped_unusable_qual,
            uncovered_breakpoint_sides: output.uncovered_breakpoint_sides.clone(),
            census: event_census,
            depth_fold: output.depth_fold.clone(),
        });

        if let Some(warning) = census::depth_fold_warning(&event_label(event), &output.depth_fold)
        {
            log::warn!("{}", warning);
        }

        adjusted_afs.push(output.adjusted_vaf);
        resistant.push(Some(event_census.fraction()));
        depth_folds.push(Some(output.depth_fold.fold));
        event_outputs.push(output);

        if let Some(site) = site {
            origin_sites.push((event_label(event), site));
        }
    }

    if let Some(message) = census::refusal_message(&refusals) {
        anyhow::bail!("{}", message);
    }

    // --edit-model origin: one draw per fragment family over every event's
    // chances (R5). Removed pool pairs move to their event's suppressed
    // names; every removed name goes to replaced_reads.txt.
    let mut origin_removed: BTreeSet<String> = BTreeSet::new();
    if edit_origin {
        let chances: Vec<origin::Chance> = event_outputs
            .iter()
            .flat_map(|o| o.origin_chances.iter().cloned())
            .collect();
        origin_removed = origin::decide(&chances, &mut rng);
        simulate::apply_removals(&mut event_outputs, &origin_removed);
        for (stat, output) in event_stats.iter_mut().zip(&event_outputs) {
            // apply_removals moved the removed pairs out of kept_originals into
            // suppressed_names, so both halves of this stat went stale at once.
            stat.kept = output.kept_originals.len();
            stat.suppressed = output.suppressed_count;
        }
        for (label, site) in &origin_sites {
            let (mut at_spot, mut elsewhere) = (0usize, 0usize);
            for f in site.fragments() {
                if origin_removed.contains(f.name) {
                    if f.at_spot(&site.footprint) {
                        at_spot += 1;
                    } else {
                        elsewhere += 1;
                    }
                }
            }
            log::info!(
                "{}: origin removed {} fragment(s) at the spot and {} at its look-alikes",
                label,
                at_spot,
                elsewhere,
            );
        }
    }

    // Names of the originals spike took out of the BAM; merge.sh removes
    // exactly these.
    let mut replaced_names = simulate::consumed_original_names(&event_outputs);
    replaced_names.extend(origin_removed.iter().cloned());
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
        &adjusted_afs,
        &resistant,
        &depth_folds,
        config.allele_fraction, // default AF for events without per-event override
        &truth_path.to_string_lossy(),
        &args.reference,
        &shared_ref,
        &ref_contigs,
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

    // Write README.md.
    let cmdline = std::env::args().collect::<Vec<_>>().join(" ");
    write_readme(
        &args.output,
        &cmdline,
        &args.bam,
        &args.reference,
        &events,
        &event_stats,
        all_output_pairs.len(),
        args.flank,
        dropped_unreplaced,
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

/// Validate overlap relationships between events.
///
/// Events are compared on their *replacement footprints*, not their spans: each
/// event is simulated independently against the original donor and replaces
/// reads across its whole footprint, so two events whose footprints intersect
/// each lay reference sequence over the other's edit and cancel it (CR1), even
/// though their spans never touch. Overlaps are rejected by default to keep
/// multi-event simulations independent. When `allow_overlap` is true, overlaps
/// are allowed but logged as warnings.
fn validate_event_overlaps(events: &[SimEvent], allow_overlap: bool) -> Result<()> {
    let mut overlaps: Vec<FootprintOverlap> = Vec::new();

    for i in 0..events.len() {
        let regions_i = event_footprints_for_overlap(&events[i]);
        for j in (i + 1)..events.len() {
            let regions_j = event_footprints_for_overlap(&events[j]);
            for region_i in &regions_i {
                for region_j in &regions_j {
                    if region_i.chrom != region_j.chrom {
                        continue;
                    }
                    let has_overlap = region_i.footprint.0 < region_j.footprint.1
                        && region_j.footprint.0 < region_i.footprint.1;
                    if has_overlap {
                        overlaps.push(FootprintOverlap {
                            event_i: i + 1,
                            event_j: j + 1,
                            chrom: region_i.chrom.clone(),
                            span_i: region_i.span,
                            footprint_i: region_i.footprint,
                            span_j: region_j.span,
                            footprint_j: region_j.footprint,
                        });
                    }
                }
            }
        }
    }

    if overlaps.is_empty() {
        return Ok(());
    }

    if allow_overlap {
        for overlap in overlaps {
            log::warn!(
                "Events {} and {} have intersecting replacement footprints on \
                 {}: spans {} and {} (footprints {} and {}). Overlap \
                 composition is approximate.",
                overlap.event_i,
                overlap.event_j,
                overlap.chrom,
                format_range(overlap.span_i),
                format_range(overlap.span_j),
                format_range(overlap.footprint_i),
                format_range(overlap.footprint_j),
            );
        }
        return Ok(());
    }

    let mut msg = format!(
        "overlapping events detected (default is to reject overlaps).\n\
         Each event replaces reads across its span grown by {}bp on each side \
         ({}bp of haplotype flank + {}bp of fragment), and two such replacement \
         footprints may not intersect: two spans on one chromosome must be at \
         least {}bp apart.\n\
         Use --allow-overlap to override.\n",
        FOOTPRINT_MARGIN,
        HAP_FLANK,
        crate::stats::MAX_FRAGMENT_LEN,
        2 * FOOTPRINT_MARGIN,
    );
    // Each line names the spans as well as the footprints: the spans are the
    // numbers the user typed (or that a --vcf record carries), and a footprint
    // appears nowhere in the input, so footprints alone leave a 50-record --vcf
    // unsearchable for the offending pair.
    for overlap in overlaps.iter().take(10) {
        msg.push_str(&format!(
            "  - events {} and {} have intersecting replacement footprints on \
             {}: spans {} and {} (footprints {} and {})\n",
            overlap.event_i,
            overlap.event_j,
            overlap.chrom,
            format_range(overlap.span_i),
            format_range(overlap.span_j),
            format_range(overlap.footprint_i),
            format_range(overlap.footprint_j),
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

/// Extra range on each side of an event's span that the event still replaces
/// reads over -- its *replacement footprint* beyond the span itself.
///
/// `build_haplotype` grows every event by `HAP_FLANK` on each side (each
/// `Haplotype::from_*` constructor fetches `start - flank .. end + flank`), and
/// `simulate_event` suppresses every original pair inside `ref_range()` and
/// tiles synthetic reads over the same range. On top of that a pair whose
/// fragment is up to `MAX_FRAGMENT_LEN` long can start that far outside the
/// range and still reach into it, so the longest fragment is added as well.
const FOOTPRINT_MARGIN: u64 = HAP_FLANK + crate::stats::MAX_FRAGMENT_LEN as u64;

/// One region of an event, with both the span the user asked for and the
/// replacement footprint around it. The span is kept alongside the footprint so
/// the rejection message can name the coordinates that appear in the user's
/// input, not only the derived ones.
struct EventFootprint {
    chrom: String,
    /// The event's own interval, 0-based half-open.
    span: (u64, u64),
    /// `span` grown by `FOOTPRINT_MARGIN` on each side.
    footprint: (u64, u64),
}

/// A pair of event regions whose replacement footprints intersect, with
/// everything the message needs: both 1-based event numbers, the chromosome,
/// and each event's span beside its footprint.
struct FootprintOverlap {
    event_i: usize,
    event_j: usize,
    chrom: String,
    span_i: (u64, u64),
    footprint_i: (u64, u64),
    span_j: (u64, u64),
    footprint_j: (u64, u64),
}

/// Format an interval the way the overlap messages print it.
fn format_range((start, end): (u64, u64)) -> String {
    format!("{}-{}", start, end)
}

/// Return the replacement footprint of each of an event's regions: the range
/// over which that event removes originals and lays synthetic reads down.
///
/// Multi-region events (`SimEvent::Fusion`) contribute one footprint per
/// breakpoint. A fusion is additive (`is_additive`, `simulate.rs`) so it
/// suppresses nothing, but it still writes synthetic reference sequence across
/// its footprint -- sequence another event may have deleted or inverted -- which
/// corrupts that event the same way, so fusions are checked like every other
/// event.
fn event_footprints_for_overlap(event: &SimEvent) -> Vec<EventFootprint> {
    event_regions_for_overlap(event)
        .into_iter()
        .map(|(chrom, start, end)| EventFootprint {
            chrom,
            span: (start, end),
            // Saturating at 0 as the haplotype constructors do at a chromosome
            // start; a footprint past a chromosome end is harmless here.
            footprint: (
                start.saturating_sub(FOOTPRINT_MARGIN),
                end.saturating_add(FOOTPRINT_MARGIN),
            ),
        })
        .collect()
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

/// Quote a value as a single shell word for the generated scripts. Single
/// quotes suppress every expansion; `'` is the one character they cannot
/// contain, so it is spliced back in as `'\''`.
fn sh_quote(value: &str) -> String {
    format!("'{}'", value.replace('\'', r"'\''"))
}

/// Resolve a path for the generated scripts. They are run from whatever
/// directory the user happens to be in, so a relative `--reference` or
/// `--bam` would resolve against the wrong one -- or nothing at all. Falls
/// back to a lexically absolute path when the file is not there to
/// canonicalize, and to the argument itself if even that fails.
///
/// `canonicalize` (over a merely-absolute path) is deliberate: both scripts
/// need the exact file spike ran on, and a directory-level symlink such as
/// `data/giab_hg38` still resolves to a reference that keeps its `.fai` (and
/// bwa/`.gzi`) indexes alongside it in the real directory, so the resolved
/// path still finds them. The hazard is the other direction: if the FASTA
/// *itself* is a symlink and the aligner's index was built beside the link
/// rather than beside the real file, resolving it orphans the index --
/// the same index/name coupling README already warns about for a bgzipped
/// `--reference`.
fn script_path(path: &str) -> String {
    std::fs::canonicalize(path)
        .or_else(|_| std::path::absolute(path))
        .map(|p| p.to_string_lossy().into_owned())
        .unwrap_or_else(|_| path.to_string())
}

/// `--samtools` is either a bare command name to look up on `$PATH` (the
/// default) or a path to a binary; only the latter needs resolving. Tell them
/// apart the way the shell itself does.
fn script_command(command: &str) -> String {
    if command.contains('/') {
        script_path(command)
    } else {
        command.to_string()
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
    let script_file = Path::new(output_dir).join("align.sh");

    // The sample comes from the BAM's own @RG, so it can hold anything the
    // SAM spec allows -- a space, an apostrophe -- and must reach the aligner
    // as one word with those characters intact.
    let rg_arg = sh_quote(&format!("@RG\\tID:sim\\tSM:{sample}\\tPL:ILLUMINA"));
    let rg_sm = sh_quote(&format!("SM:{sample}"));
    let align_cmd = match aligner {
        "bwa-mem2" => format!("bwa-mem2 mem -t \"$THREADS\" \\\n    -R {rg_arg} \\\n    \"$REF\" \\\n    \"$DIR/R1.fq.gz\" \"$DIR/R2.fq.gz\" \\\n    2>\"$DIR/align.log\""),
        "minimap2" => format!("minimap2 -a -x sr -t \"$THREADS\" \\\n    -R {rg_arg} \\\n    \"$REF\" \\\n    \"$DIR/R1.fq.gz\" \"$DIR/R2.fq.gz\" \\\n    2>\"$DIR/align.log\""),
        "bowtie2" => format!("bowtie2 -x \"$REF\" \\\n    -1 \"$DIR/R1.fq.gz\" -2 \"$DIR/R2.fq.gz\" \\\n    -p \"$THREADS\" \\\n    --rg-id sim --rg {rg_sm} --rg PL:ILLUMINA \\\n    2>\"$DIR/align.log\""),
        custom => format!(
            "{custom} \"$REF\" \"$DIR/R1.fq.gz\" \"$DIR/R2.fq.gz\" \\\n    2>\"$DIR/align.log\"",
        ),
    };

    // Defaults are baked into the script, so they must be absolute (align.sh
    // is run from anywhere) and quoted (a path may hold `}`, `"`, `$`, a
    // backtick or a space). `VAR=${1:-'...'}` is safe unquoted: an assignment
    // right-hand side is never word-split or globbed.
    let ref_default = sh_quote(&script_path(ref_path));
    let samtools_default = sh_quote(&script_command(samtools));
    // The aligner name/command line stays verbatim in align_cmd (a custom
    // --aligner is a command line, quoting it there would break the
    // documented usage) but this echo banner is a second, undocumented ride
    // on the same value -- and unlike align_cmd it sits inside a double-quoted
    // string, so it needs quoting even for a value that is already a
    // correctly-quoted shell command line.
    let aligner_echo = sh_quote(aligner);

    let script = format!(
        r#"#!/bin/bash
set -euo pipefail
# Align simulated reads and sort.
# Usage: bash align.sh [REF] [THREADS]
REF=${{1:-{ref_default}}}
THREADS="${{2:-{threads}}}"
SAMTOOLS={samtools_default}
DIR="$(cd "$(dirname "$0")" && pwd)"

echo "Aligning $DIR/R1.fq.gz + R2.fq.gz ("{aligner_echo}", $THREADS threads)..."
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

    std::fs::write(&script_file, script)?;

    // Make executable.
    #[cfg(unix)]
    {
        use std::os::unix::fs::PermissionsExt;
        std::fs::set_permissions(&script_file, std::fs::Permissions::from_mode(0o755))?;
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
) -> Result<(ReadPool, String, usize)> {
    let mut all_pairs: Vec<ReadPair> = Vec::new();
    let mut windows_searched: Vec<String> = Vec::new();
    let unusable_before = unusable_qual_names.len();
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
        windows_searched.extend(window_labels(chrom_a, &windows_a));
        windows_searched.extend(window_labels(chrom_b, &windows_b));
        pool_chrom = chrom_a.clone();
    } else {
        // Single-region events (DEL, DUP, INV, INS). The two rules that key
        // on how many loci an event is drawn from live in two files, so tie
        // them together here: a multi-locus event that reached this branch
        // would fill its pool from one window and then be graded by
        // `donor_coverage_for_tiling`'s permissive rule.
        debug_assert!(
            !event.is_multi_locus(),
            "a multi-locus event must extract one window per locus: {:?}",
            event.primary_region()
        );
        let (chrom, start, end) = event.primary_region().unwrap();
        let windows = extraction_bounds(chrom, start, end, config.flank_bp, extraction_region);
        extract_windows(config, chrom, &windows, &mut all_pairs, unusable_qual_names)?;
        windows_searched.extend(window_labels(chrom, &windows));
        pool_chrom = chrom.to_string();
    }

    // Names this event's windows added: `unusable_qual_names` is global, so
    // the delta is what *this* event lost to unreadable qualities.
    let dropped_unusable_qual = unusable_qual_names.len() - unusable_before;
    let pool = finish_donor_pool(
        all_pairs,
        event,
        &windows_searched,
        dropped_unusable_qual,
        config.read_length,
    )?;
    Ok((pool, pool_chrom, dropped_unusable_qual))
}

/// Format one side's extraction windows as `chrom:start-end`, for messages.
fn window_labels(chrom: &str, windows: &[(u64, u64)]) -> Vec<String> {
    windows
        .iter()
        .map(|(start, end)| format!("{}:{}-{}", chrom, start, end))
        .collect()
}

/// One event's contribution to the run, for the log line and the run README.
struct EventStat {
    /// The fraction that was *requested* for this event (`af=`, or
    /// `--allele-fraction`). `truth.vcf` records it as `SIM_REQ_VAF`.
    vaf: f64,
    /// The fraction this event's fragments actually plant, when the additive
    /// 0.95 cap or the two-fragment floor moved the count off `vaf`; `None`
    /// when the request stands. `truth.vcf` records it as `SIM_VAF`, so the
    /// run README prints it beside the request rather than the request
    /// alone -- printing only the request made the two files disagree on a
    /// capped or floored run.
    adjusted_vaf: Option<f64>,
    kept: usize,
    chimeric: usize,
    suppressed: usize,
    /// Pairs this event's extraction dropped for unusable base qualities
    /// (M14). They are in neither `kept` nor `suppressed`: `merge.sh` removes
    /// them from the merged BAM and nothing replaces them, so they are a
    /// depth dip in this event's window and nowhere else.
    dropped_unusable_qual: usize,
    /// Breakpoint sides of this event the donor pool had no reads over. The
    /// event is kept -- one bare side is the far edge of a sliced or panel
    /// BAM -- but the tiled fragments that land there were scaled by depth
    /// measured somewhere else, so the run README says which sides they are.
    uncovered_breakpoint_sides: Vec<String>,
    /// The reads over this event and how many of them spike could not edit
    /// (CR4). `truth.vcf` records the share as `SIM_RESIST`.
    census: census::Census,
    /// How far the donor's depth, where this event's fragments were drawn,
    /// departs from the one depth they were scaled by (CR2). `truth.vcf`
    /// records the fold as `SIM_DEPTH_FOLD`.
    depth_fold: types::DepthFold,
}

/// Smallest donor pool spike will simulate one event from.
///
/// Below this the output is invention rather than simulation: everything
/// spike puts in a read comes from this pool -- the quality profile, the
/// fragment lengths, the coverage the tiling count is scaled by -- and each
/// of those silently substitutes a constant when the pool runs out.
///
/// 30 is `synth.rs`'s own `MIN_BASE_OBS`, the observation count it requires
/// before it will sample from a quality bin, and a pool of *n* pairs puts
/// exactly *n* observations in each cycle-only bin -- the profile's final
/// fallback, and the only level an ordinary pool always reaches. So 30 pairs
/// is the smallest pool at which any bin of the profile can meet the
/// threshold the profile itself sets; below it every bin is unusable by that
/// rule and, at zero pairs, sampling returns `synth.rs`'s last-resort
/// constant Q20 byte for every base. This is a floor on "measured from this
/// library at all", not a claim that 30 pairs is enough coverage for a good
/// simulation.
const MIN_DONOR_PAIRS: usize = 30;

/// Dedup one event's extracted pairs and turn them into its donor pool.
///
/// Split out of `extract_pool_for_event` so the pool the simulation runs on
/// can be checked without a BAM behind it.
fn finish_donor_pool(
    mut all_pairs: Vec<ReadPair>,
    event: &SimEvent,
    windows_searched: &[String],
    dropped_unusable_qual: usize,
    read_length: usize,
) -> Result<ReadPool> {
    // Windows can share reads, and the same fragment must not enter the pool
    // -- or the fragment distribution -- twice. Dedup before counting: a
    // fragment handed in by two overlapping windows is not extra material.
    extract::dedup_pairs_by_name(&mut all_pairs);

    if all_pairs.len() < MIN_DONOR_PAIRS {
        bail!(
            "event {} has too few usable donor reads: {} read pair(s) extracted from \
             {}, fewer than the {} spike needs ({} read pair(s) in those windows \
             were dropped for unusable base qualities and are not in that count). \
             Every \
             simulated read is built from this pool -- its base qualities, its \
             fragment lengths and the coverage the tiling count is scaled by all come \
             from it -- so spike would invent reads rather than simulate them, and \
             still write a truth VCF beside them. Check that the event lies in a \
             covered region of the BAM, that --region (if given) covers it, that \
             --min-mapq is not filtering the window out, and -- if the drop count \
             above accounts for the shortfall -- that the file's records carry base \
             qualities spike can read (a CRAM storing them as read features does not).",
            event_label(event),
            all_pairs.len(),
            windows_searched.join(", "),
            MIN_DONOR_PAIRS,
            dropped_unusable_qual,
        );
    }

    // The model is built over the same range the generator samples in, so
    // `read_length` has to reach it -- see `FragmentDist::from_read_pairs`.
    let frag_dist = stats::FragmentDist::from_read_pairs(&all_pairs, read_length);
    Ok(extract::build_read_pool(all_pairs, frag_dist))
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

/// The reference range an event's haplotype covers: its region grown by
/// `HAP_FLANK` on each side, stopped at the contig's ends, as
/// `VariantHaplotype::ref_range` reports it once the haplotype is built.
/// `simulate` checks the two agree. `None` for a fusion, which has no single
/// region.
fn origin_footprint(event: &SimEvent, contig_len: u64) -> Option<origin::Span> {
    let (chrom, start, end) = event.primary_region()?;
    Some(origin::Span::new(
        chrom,
        start.saturating_sub(HAP_FLANK),
        end.saturating_add(HAP_FLANK).min(contig_len),
    ))
}

/// Build a VariantHaplotype for a given event, which it may mutate: an
/// insertion with no supplied sequence has the sequence generated here
/// stored back into its `ins_seq`, so the truth VCF written after the
/// haplotype loop can spell out the bases the reads were cut from. No
/// other event type and no other field is touched.
fn build_haplotype(
    event: &mut SimEvent,
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
            let hap = VariantHaplotype::from_insertion(reference, chrom, *pos, &seq, flank);
            // CR7: keep the sequence on the event. The truth VCF is written
            // after this loop, so this is what lets its ALT spell out the
            // bases the reads were cut from instead of a symbolic <INS>.
            *ins_seq = Some(seq);
            hap
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
    let script_file = Path::new(output_dir).join("merge.sh");

    // As in align.sh: baked-in defaults must be absolute (merge.sh is run
    // from anywhere) and quoted (a path may hold `}`, `"`, `$`, a backtick or
    // a space).
    let original_default = sh_quote(&script_path(original_bam));
    let ref_default = sh_quote(&script_path(ref_path));
    let samtools_default = sh_quote(&script_command(samtools));

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
ORIGINAL=${{1:-{original_default}}}
REF=${{2:-{ref_default}}}
THREADS="${{3:-{threads}}}"
SAMTOOLS={samtools_default}
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

    std::fs::write(&script_file, script)?;

    #[cfg(unix)]
    {
        use std::os::unix::fs::PermissionsExt;
        std::fs::set_permissions(&script_file, std::fs::Permissions::from_mode(0o755))?;
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
    event_stats: &[EventStat],
    total_pairs: usize,
    flank: u64,
    dropped_unreplaced: usize,
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
    // Two VAF columns, the same pair truth.vcf records as SIM_REQ_VAF and
    // SIM_VAF: printing the request alone made this table contradict the
    // truth VCF written beside it on every capped or floored run.
    writeln!(
        md,
        "| # | Event | Requested VAF | Simulated VAF | Kept reads | Chimeric reads | \
         Suppressed reads | Dropped (unusable quality) | Resistant reads | Depth fold |"
    )?;
    writeln!(
        md,
        "|---|-------|---------------|---------------|-----------|----------------|\
         -----------------|---------------------------|-----------------|------------|"
    )?;
    for (i, event) in events.iter().enumerate() {
        let label = event_label(event);
        let stat = event_stats.get(i);
        writeln!(
            md,
            "| {} | {} | {:.3} | {:.3} | {} | {} | {} | {} | {} | {:.2} |",
            i + 1,
            label,
            stat.map_or(0.0, |s| s.vaf),
            stat.map_or(0.0, |s| s.adjusted_vaf.unwrap_or(s.vaf)),
            stat.map_or(0, |s| s.kept),
            stat.map_or(0, |s| s.chimeric),
            stat.map_or(0, |s| s.suppressed),
            stat.map_or(0, |s| s.dropped_unusable_qual),
            stat.map_or("-".to_string(), |s| format!(
                "{} of {} ({:.0}%)",
                s.census.resistant,
                s.census.counted,
                s.census.fraction() * 100.0
            )),
            stat.map_or(0.0, |s| s.depth_fold.fold),
        )?;
    }
    writeln!(md)?;
    // Only for an event above the warning threshold, like the log's warning.
    let resistant: Vec<String> = events
        .iter()
        .enumerate()
        .filter_map(|(i, event)| {
            let stat = event_stats.get(i)?;
            (stat.census.fraction() > census::WARN_ABOVE).then(|| {
                format!("{} ({:.0}%)", event_label(event), stat.census.fraction() * 100.0)
            })
        })
        .collect();
    if !resistant.is_empty() {
        writeln!(
            md,
            "**Reads spike could not edit:** {}. The resistant column counts the \
             primary, non-duplicate reads over each event that are not in its donor pool: \
             below `--min-mapq`, not a proper pair, or a mate that fails a filter. They \
             stay in the merged BAM as they are, so these events are weaker than requested \
             by about that share. `truth.vcf` records it per event as `SIM_RESIST`.",
            resistant.join("; ")
        )?;
        writeln!(md)?;
    }
    // Only for an event above the warning threshold, like the log's warning.
    let off: Vec<String> = events
        .iter()
        .enumerate()
        .filter_map(|(i, event)| {
            let fold = &event_stats.get(i)?.depth_fold;
            (fold.fold > census::DEPTH_FOLD_WARN_ABOVE).then(|| {
                format!(
                    "{} ({:.1}x over {} against {:.1}x, {:.2}-fold)",
                    event_label(event),
                    fold.worst_depth,
                    fold.worst_bin,
                    fold.scaled_by,
                    fold.fold
                )
            })
        })
        .collect();
    if !off.is_empty() {
        writeln!(
            md,
            "**Donor depth off the scaling depth:** {}. spike scales every fragment an \
             event tiles by one depth, measured at one of its breakpoints. The depth fold \
             column is the largest fold between that and the donor's own depth in any ~1 kb \
             bin the fragments are drawn from; where the two differ, the event's depth there \
             is off by about that fold. `truth.vcf` records it per event as \
             `SIM_DEPTH_FOLD`.",
            off.join("; ")
        )?;
        writeln!(md)?;
    }
    // Only when the two columns differ somewhere: on an ordinary run they are
    // the same number on every row and an unconditional paragraph would train
    // the reader to skip it.
    let moved: Vec<String> = events
        .iter()
        .enumerate()
        .filter_map(|(i, event)| {
            let stat = event_stats.get(i)?;
            let simulated = stat.adjusted_vaf?;
            Some(format!(
                "{} ({:.3} requested, {:.3} simulated)",
                event_label(event),
                stat.vaf,
                simulated
            ))
        })
        .collect();
    if !moved.is_empty() {
        writeln!(
            md,
            "**Requested and simulated VAF differ:** {}. A request above the \
             additive 0.95 cap is simulated at the cap; a request too thin to \
             ask for two tiled fragments is simulated at what two of them \
             plant. The simulated column is what the reads carry, and \
             `truth.vcf` records the same two numbers per event as `SIM_VAF` \
             (simulated) and `SIM_REQ_VAF` (requested).",
            moved.join("; ")
        )?;
        writeln!(md)?;
    }
    let dropped_total: usize = event_stats.iter().map(|s| s.dropped_unusable_qual).sum();
    writeln!(
        md,
        "**Dropped (unusable quality):** {} read pair(s) across all events had no \
         readable base qualities, so they are in neither the kept nor the suppressed \
         column. {} of them are removed from the merged BAM by `merge.sh` with \
         nothing put back in their place -- a depth dip confined to the event \
         windows, which a depth-based caller can read as signal. 0 in both columns \
         means the input's qualities were all readable.",
        dropped_total, dropped_unreplaced
    )?;
    writeln!(md)?;
    // Only when it happened: on ordinary input every side is covered and an
    // unconditional paragraph would train the reader to skip it.
    let bare: Vec<String> = events
        .iter()
        .enumerate()
        .filter_map(|(i, event)| {
            let sides = &event_stats.get(i)?.uncovered_breakpoint_sides;
            (!sides.is_empty()).then(|| format!("{} ({})", event_label(event), sides.join(", ")))
        })
        .collect();
    if !bare.is_empty() {
        writeln!(
            md,
            "**Breakpoint sides with no donor coverage:** {}. The event was kept -- \
             a bare side is the far edge of a sliced or panel BAM, not a reason to \
             refuse -- but its haplotype spans every side, so the fragments tiled \
             across the bare part were scaled by donor depth measured at the *other* \
             side and land where the input BAM has no read. Expect a coverage island \
             there, and reads whose fragment lengths and qualities came from \
             somewhere else in the genome.",
            bare.join("; ")
        )?;
        writeln!(md)?;
    }
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
    fn test_validate_read_length_accepts_exact_max() {
        // L5 boundary: MAX_FRAGMENT_LEN (1500) itself must be accepted — the
        // guard rejects lengths *above* the max, not lengths *at* it. A
        // future `>` -> `>=` typo would reject this and fail silently
        // otherwise, since only far-below (151) and far-above (2000) values
        // were previously exercised.
        assert!(validate_read_length(1500).is_ok());
    }

    #[test]
    fn test_validate_read_length_rejects_one_above_max() {
        // L5 boundary: one bp past MAX_FRAGMENT_LEN (1501) must be rejected.
        // Pins the exact threshold, complementing the 1500-accepts case.
        let err = validate_read_length(1501).expect_err("1501bp should be rejected");
        let msg = err.to_string();
        assert!(msg.contains("1501"), "message should cite the read length: {msg}");
        assert!(msg.contains("1500"), "message should cite the max: {msg}");
    }

    #[test]
    fn test_validate_allele_fraction_rejects_nan_and_infinities() {
        // L8: `af <= 0.0 || af > 1.0` lets NaN through (both comparisons are
        // false for NaN), so `--allele-fraction NaN` used to be accepted,
        // suppressing every read at simulate time while the truth VCF
        // recorded SIM_VAF=NaN.
        assert!(validate_allele_fraction(f64::NAN).is_err());
        assert!(validate_allele_fraction(f64::INFINITY).is_err());
        assert!(validate_allele_fraction(f64::NEG_INFINITY).is_err());
    }

    #[test]
    fn test_help_states_the_bounds_the_code_enforces() {
        // The help said `--allele-fraction (0.0-1.0)` while the code rejects
        // 0, and said nothing at all about `--flank`'s 2000 minimum -- the
        // likeliest accidental refusal. A user who reads --help and is
        // refused anyway has been told the wrong thing.
        use clap::CommandFactory;
        let help = Args::command().render_long_help().to_string();

        assert!(
            help.contains("(0.0, 1.0]"),
            "--allele-fraction's help must state the interval the code enforces"
        );
        assert!(
            validate_allele_fraction(0.0).is_err(),
            "the help above claims 0 is refused"
        );

        assert!(
            help.contains(&HAP_FLANK.to_string()),
            "--flank's help must state its {} minimum",
            HAP_FLANK
        );
        assert!(
            validate_flank(HAP_FLANK - 1).is_err(),
            "the help above claims below {} is refused",
            HAP_FLANK
        );

        // CR1 widened the overlap check from event spans to replacement
        // footprints, so --allow-overlap now also gates pairs that do not
        // overlap at all: two spans on one chromosome closer than
        // 2 * FOOTPRINT_MARGIN are refused as well. A user reading "overlapping
        // events are rejected" and then refused for two events 1kb apart has
        // been told the wrong thing.
        let footprint_gap = 2 * FOOTPRINT_MARGIN;
        assert!(
            help.contains(&footprint_gap.to_string()),
            "--allow-overlap's help must state the {footprint_gap}bp gap the code enforces"
        );
        assert!(
            validate_event_overlaps(
                &[
                    del("chr1", 100_000, 200_000),
                    del("chr1", 200_000 + footprint_gap - 1, 300_000),
                ],
                false,
            )
            .is_err(),
            "the help above claims spans closer than {footprint_gap}bp are refused"
        );
    }

    #[test]
    fn test_validate_allele_fraction_rejects_zero() {
        assert!(validate_allele_fraction(0.0).is_err());
    }

    #[test]
    fn test_validate_allele_fraction_accepts_upper_boundary() {
        // L8 boundary: 1.0 itself must stay accepted (af is in `(0.0, 1.0]`).
        assert!(validate_allele_fraction(1.0).is_ok());
    }

    #[test]
    fn test_validate_allele_fraction_rejects_above_one() {
        assert!(validate_allele_fraction(1.5).is_err());
    }

    #[test]
    fn test_validate_allele_fraction_rejects_negative_zero() {
        // L8 boundary coverage: `-0.0 > 0.0` is false under IEEE-754, so
        // `-0.0` was already rejected by both the old and new forms of the
        // check. Added because the fix brief asked for this exact value to
        // be pinned, not because it was ever broken.
        assert!(validate_allele_fraction(-0.0).is_err());
    }

    #[test]
    fn test_validate_allele_fraction_rejects_just_above_one() {
        // L8 boundary coverage: `test_validate_allele_fraction_rejects_above_one`'s
        // 1.5 exercises the same `<= 1.0` branch as any value clearly over
        // 1; this pins a value just barely over the boundary instead.
        assert!(validate_allele_fraction(1.0000001).is_err());
    }

    #[test]
    fn test_af_het_is_exactly_one_half() {
        // CR7(c): af=het drew the event's fraction from Beta(40,40), so a
        // "heterozygote" was built as a dosage mixture: above 0.5 the other
        // copy lost reads too. A het event is one copy of two.
        // Before the fix the draw came from the run's --seed; seed 0 gave
        // 0.4605. The helper takes no RNG now, so no seed can move it.
        assert_eq!(resolve_af_spec(Some(AfSpec::Het)), Some(0.5));
        assert_eq!(resolve_af_spec(Some(AfSpec::Exact(0.15))), Some(0.15));
        assert_eq!(resolve_af_spec(Some(AfSpec::Hom)), Some(1.0));
        assert_eq!(resolve_af_spec(None), None);
    }

    #[test]
    fn test_a_fusion_run_warns_that_it_adds_one_junction_only() {
        // CR6: fusion mode keeps every original read and adds junction reads
        // on top. That is not a balanced translocation, and a user reading
        // "fusion" could take it for one.
        let warning = fusion_mode_warning(&[del("chr1", 100, 200), fusion("chr1", 5, "chr2", 9)])
            .expect("a run with a fusion must say what fusion mode is");
        assert!(warning.contains("1 fusion event"), "{}", warning);
        assert!(warning.contains("adds junction reads"), "{}", warning);
        assert!(warning.contains("not a balanced translocation"), "{}", warning);
        assert!(warning.contains("no reciprocal junction"), "{}", warning);

        assert_eq!(fusion_mode_warning(&[del("chr1", 100, 200), ins("chr1", 500)]), None);
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
if [ -n "${SAMTOOLS_ARGV_OUT:-}" ]; then
  printf '%s\n' "$@" >> "$SAMTOOLS_ARGV_OUT"
fi
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

    #[test]
    fn test_align_script_stays_valid_with_a_quoted_custom_aligner() {
        // A correctly-quoted custom --aligner command line is left verbatim
        // in the aligner-invocation position by design (README documents
        // this). But the same value is also interpolated raw into the echo
        // banner -- a second, undocumented ride on the value -- where even a
        // well-formed command line broke the script's own syntax.
        let dir = scratch_dir("align_custom_aligner_echo");
        let aligner = r#"al --preset 'a"b'"#;

        write_align_script(dir.to_str().unwrap(), "ref.fa", 4, aligner, "samtools", "HG002")
            .unwrap();

        let output = std::process::Command::new("bash")
            .arg("-n")
            .arg(dir.join("align.sh"))
            .output()
            .unwrap();

        assert!(
            output.status.success(),
            "align.sh must stay syntactically valid when a correctly-quoted \
             --aligner value reaches the echo line:\nstderr: {}",
            String::from_utf8_lossy(&output.stderr)
        );
    }

    /// A file name holding every character that has broken these generated
    /// scripts -- a space, `}` (which ends a `${...}` expansion early), `"`,
    /// `$`, a backtick and a single quote -- plus `*` and a tab. The
    /// backtick runs `touch pwned` in the script's working directory if the
    /// shell ever evaluates the path, which is what the tests below watch
    /// for. `*` and the tab are not things that have broken the scripts;
    /// they pin the fix's riskiest-looking claim end-to-end -- that
    /// `VAR=${1:-'name'}` never word-splits or globs its default, whatever
    /// the default contains -- rather than leaving it covered only by
    /// comment and by manual testing during development. A newline is not
    /// included: this file's `argv_lines` helper records a stub's argv one
    /// element per line, so an argv element that itself contains a newline
    /// would (correctly) come back looking like two elements -- a test-
    /// harness limitation, not a bug in the generated scripts.
    const HOSTILE_NAME: &str = "we ird}\"$HOME`touch pwned`'x*\ty";

    /// A stub aligner that dumps its argv, one argument per line, to
    /// `$ARGV_OUT`. Its (empty) stdout is what `samtools sort` consumes.
    fn write_stub_aligner(path: &std::path::Path) {
        let script = "#!/bin/bash\nprintf '%s\\n' \"$@\" > \"$ARGV_OUT\"\n";
        std::fs::write(path, script).unwrap();
        #[cfg(unix)]
        {
            use std::os::unix::fs::PermissionsExt;
            std::fs::set_permissions(path, std::fs::Permissions::from_mode(0o755)).unwrap();
        }
    }

    /// Run a generated script with no arguments, so it uses the defaults it
    /// was written with, from a working directory that is not spike's.
    fn run_generated_script(
        dir: &std::path::Path,
        name: &str,
        envs: &[(&str, &str)],
    ) -> std::process::Output {
        let mut cmd = std::process::Command::new("bash");
        cmd.arg(dir.join(name)).current_dir(dir);
        for (key, value) in envs {
            cmd.env(key, value);
        }
        cmd.env(
            "PATH",
            format!(
                "{}:{}",
                dir.join("bin").display(),
                std::env::var("PATH").unwrap_or_default()
            ),
        );
        cmd.output().unwrap()
    }

    fn argv_lines(path: &std::path::Path) -> Vec<String> {
        std::fs::read_to_string(path)
            .unwrap_or_default()
            .lines()
            .map(|l| l.to_string())
            .collect()
    }

    /// Set up a scratch dir with `bin/<aligner>` and a stub samtools, and
    /// return (dir, aligner argv file, samtools path).
    fn align_script_fixture(
        name: &str,
        aligner: &str,
    ) -> (std::path::PathBuf, std::path::PathBuf, std::path::PathBuf) {
        let dir = scratch_dir(name);
        let bin = dir.join("bin");
        std::fs::create_dir_all(&bin).unwrap();
        write_stub_aligner(&bin.join(aligner));
        let samtools = write_stub_samtools(&dir);
        let argv_out = dir.join("aligner_argv.txt");
        let _ = std::fs::remove_file(&argv_out);
        let _ = std::fs::remove_file(dir.join("pwned"));
        (dir, argv_out, samtools)
    }

    #[test]
    fn test_align_script_defaults_survive_shell_metacharacters_in_paths() {
        let (dir, argv_out, samtools) =
            align_script_fixture("align_hostile_ref", "bwa-mem2");
        let ref_path = dir.join(format!("{}.fa", HOSTILE_NAME));
        std::fs::write(&ref_path, b">chr1\nACGT\n").unwrap();

        write_align_script(
            dir.to_str().unwrap(),
            ref_path.to_str().unwrap(),
            4,
            "bwa-mem2",
            samtools.to_str().unwrap(),
            "Patient 123",
        )
        .unwrap();

        let output = run_generated_script(&dir, "align.sh", &[("ARGV_OUT", argv_out.to_str().unwrap())]);

        assert!(
            output.status.success(),
            "align.sh must run when --reference holds shell metacharacters:\nstdout: {}\nstderr: {}",
            String::from_utf8_lossy(&output.stdout),
            String::from_utf8_lossy(&output.stderr)
        );
        assert!(
            !dir.join("pwned").exists(),
            "align.sh ran the backtick in the reference path as a command"
        );
        let argv = argv_lines(&argv_out);
        let want = std::fs::canonicalize(&ref_path).unwrap();
        assert!(
            argv.iter().any(|a| std::path::Path::new(a) == want),
            "the aligner must be handed the reference path unmangled; wanted {}, got {:?}",
            want.display(),
            argv
        );
        assert!(
            argv.iter()
                .any(|a| a == "@RG\\tID:sim\\tSM:Patient 123\\tPL:ILLUMINA"),
            "the aligner must be handed the sample name unmangled; got {:?}",
            argv
        );
    }

    #[test]
    fn test_align_script_default_reference_is_absolute() {
        // A relative --reference is relative to spike's working directory;
        // align.sh is normally run from somewhere else entirely, so a
        // relative default resolves against the wrong directory or not at all.
        let (dir, argv_out, samtools) =
            align_script_fixture("align_relative_ref", "bwa-mem2");

        write_align_script(
            dir.to_str().unwrap(),
            "ref.fa",
            4,
            "bwa-mem2",
            samtools.to_str().unwrap(),
            "HG002",
        )
        .unwrap();

        let output = run_generated_script(&dir, "align.sh", &[("ARGV_OUT", argv_out.to_str().unwrap())]);

        assert!(
            output.status.success(),
            "align.sh failed:\nstdout: {}\nstderr: {}",
            String::from_utf8_lossy(&output.stdout),
            String::from_utf8_lossy(&output.stderr)
        );
        let argv = argv_lines(&argv_out);
        let want = std::env::current_dir().unwrap().join("ref.fa");
        assert!(
            argv.iter().any(|a| std::path::Path::new(a) == want),
            "a relative --reference must be resolved against spike's working \
             directory; wanted {}, got {:?}",
            want.display(),
            argv
        );
    }

    #[test]
    fn test_align_script_passes_a_sample_name_with_a_space_as_one_argument() {
        // bowtie2's `--rg SM:...` is the one read-group argument that is not
        // already inside single quotes.
        let (dir, argv_out, samtools) = align_script_fixture("align_bowtie2_sm", "bowtie2");

        write_align_script(
            dir.to_str().unwrap(),
            "ref.fa",
            4,
            "bowtie2",
            samtools.to_str().unwrap(),
            "Patient 123",
        )
        .unwrap();

        let output = run_generated_script(&dir, "align.sh", &[("ARGV_OUT", argv_out.to_str().unwrap())]);

        assert!(
            output.status.success(),
            "align.sh failed:\nstdout: {}\nstderr: {}",
            String::from_utf8_lossy(&output.stdout),
            String::from_utf8_lossy(&output.stderr)
        );
        let argv = argv_lines(&argv_out);
        assert!(
            argv.iter().any(|a| a == "SM:Patient 123"),
            "bowtie2 must receive SM: and the sample name as one argument; got {:?}",
            argv
        );
    }

    #[test]
    fn test_align_script_passes_a_sample_name_with_a_single_quote_through_the_rg_argument() {
        // rg_arg (bwa-mem2/minimap2's -R) is the one place that runs the
        // `'\''` splice from sh_quote through a real shell -- the path tests
        // exercise that splice on file paths, but no test had ever pushed a
        // sample name with a single quote through it.
        let (dir, argv_out, samtools) = align_script_fixture("align_bwa_mem2_quote_sm", "bwa-mem2");

        write_align_script(
            dir.to_str().unwrap(),
            "ref.fa",
            4,
            "bwa-mem2",
            samtools.to_str().unwrap(),
            "O'Brien",
        )
        .unwrap();

        let output = run_generated_script(&dir, "align.sh", &[("ARGV_OUT", argv_out.to_str().unwrap())]);

        assert!(
            output.status.success(),
            "align.sh failed:\nstdout: {}\nstderr: {}",
            String::from_utf8_lossy(&output.stdout),
            String::from_utf8_lossy(&output.stderr)
        );
        let argv = argv_lines(&argv_out);
        assert!(
            argv.iter().any(|a| a == "@RG\\tID:sim\\tSM:O'Brien\\tPL:ILLUMINA"),
            "bwa-mem2 must receive the sample name's apostrophe unmangled and \
             as one -R argument; got {:?}",
            argv
        );
    }

    #[test]
    fn test_merge_script_defaults_survive_shell_metacharacters_in_paths() {
        let dir = scratch_dir("merge_hostile_paths");
        let samtools = write_stub_samtools(&dir);
        let original = dir.join(format!("{}.bam", HOSTILE_NAME));
        let ref_path = dir.join(format!("{}.fa", HOSTILE_NAME));
        std::fs::write(&original, b"").unwrap();
        std::fs::write(&ref_path, b"").unwrap();
        let argv_out = dir.join("samtools_argv.txt");
        let _ = std::fs::remove_file(&argv_out);
        let _ = std::fs::remove_file(dir.join("pwned"));

        write_merge_script(
            dir.to_str().unwrap(),
            original.to_str().unwrap(),
            ref_path.to_str().unwrap(),
            4,
            samtools.to_str().unwrap(),
        )
        .unwrap();
        std::fs::write(dir.join("sim.bam"), b"").unwrap();
        std::fs::write(dir.join("replaced_reads.txt"), "r1\nr2\n").unwrap();

        let output = run_generated_script(
            &dir,
            "merge.sh",
            &[
                ("SAMTOOLS_ARGV_OUT", argv_out.to_str().unwrap()),
                // Healthy case: the guard is not what this test is about.
                ("SAMTOOLS_FAKE_ACTUAL", "10"),
            ],
        );

        assert!(
            output.status.success(),
            "merge.sh must run when the BAM and reference paths hold shell \
             metacharacters:\nstdout: {}\nstderr: {}",
            String::from_utf8_lossy(&output.stdout),
            String::from_utf8_lossy(&output.stderr)
        );
        assert!(
            !dir.join("pwned").exists(),
            "merge.sh ran the backtick in a path as a command"
        );
        let argv = argv_lines(&argv_out);
        for want in [
            std::fs::canonicalize(&original).unwrap(),
            std::fs::canonicalize(&ref_path).unwrap(),
        ] {
            assert!(
                argv.iter().any(|a| std::path::Path::new(a) == want),
                "samtools must be handed {} unmangled; got {:?}",
                want.display(),
                argv
            );
        }
    }

    #[test]
    fn test_align_script_accepts_a_positional_reference_holding_a_space_and_a_dollar_sign() {
        // scripts/validate_pipeline.sh (and spike's own run_alignment) never
        // rely on align.sh's baked-in defaults -- they always call it
        // positionally: `bash align.sh "$REFERENCE" "$THREADS"`. The quoting
        // and absolutising this task added is only for the defaults; this
        // pins that the already-positional path -- untouched by this task --
        // still carries a value holding a space and a `$` through as one
        // unmangled argument, the way validate_pipeline.sh depends on it to.
        let (dir, argv_out, samtools) =
            align_script_fixture("align_positional_ref", "bwa-mem2");
        let positional_ref = dir.join("pos arg $HOME.fa");
        std::fs::write(&positional_ref, b">chr1\nACGT\n").unwrap();

        write_align_script(
            dir.to_str().unwrap(),
            "ref.fa", // baked default; irrelevant, the positional argument overrides it
            4,
            "bwa-mem2",
            samtools.to_str().unwrap(),
            "HG002",
        )
        .unwrap();

        let output = std::process::Command::new("bash")
            .arg(dir.join("align.sh"))
            .arg(&positional_ref)
            .arg("4")
            .current_dir(&dir)
            .env("ARGV_OUT", argv_out.to_str().unwrap())
            .env(
                "PATH",
                format!(
                    "{}:{}",
                    dir.join("bin").display(),
                    std::env::var("PATH").unwrap_or_default()
                ),
            )
            .output()
            .unwrap();

        assert!(
            output.status.success(),
            "align.sh must run when REFERENCE is passed positionally holding a \
             space and a $, exactly how validate_pipeline.sh invokes it:\nstdout: {}\nstderr: {}",
            String::from_utf8_lossy(&output.stdout),
            String::from_utf8_lossy(&output.stderr)
        );
        let argv = argv_lines(&argv_out);
        assert!(
            argv.iter().any(|a| std::path::Path::new(a) == positional_ref),
            "the aligner must be handed the positional reference unmangled; \
             wanted {}, got {:?}",
            positional_ref.display(),
            argv
        );
    }

    #[test]
    fn test_merge_script_accepts_positional_arguments_holding_a_space_and_a_dollar_sign() {
        // validate_pipeline.sh's other call site: `bash merge.sh "$BG_BAM"
        // "$REFERENCE" "$THREADS"`, always positional. Same pin as above, for
        // merge.sh's ORIGINAL_BAM and REFERENCE.
        let dir = scratch_dir("merge_positional_paths");
        let samtools = write_stub_samtools(&dir);
        let original = dir.join("pos arg $HOME.bam");
        let ref_path = dir.join("pos ref $USER.fa");
        std::fs::write(&original, b"").unwrap();
        std::fs::write(&ref_path, b"").unwrap();
        let argv_out = dir.join("samtools_argv_positional.txt");
        let _ = std::fs::remove_file(&argv_out);

        write_merge_script(
            dir.to_str().unwrap(),
            "orig.bam", // baked default; irrelevant, the positional arguments override it
            "ref.fa",
            4,
            samtools.to_str().unwrap(),
        )
        .unwrap();
        std::fs::write(dir.join("sim.bam"), b"").unwrap();
        std::fs::write(dir.join("replaced_reads.txt"), "r1\nr2\n").unwrap();

        let output = std::process::Command::new("bash")
            .arg(dir.join("merge.sh"))
            .arg(&original)
            .arg(&ref_path)
            .arg("4")
            .current_dir(&dir)
            .env("SAMTOOLS_ARGV_OUT", argv_out.to_str().unwrap())
            // Healthy case: the shortfall guard is not what this test is about.
            .env("SAMTOOLS_FAKE_ACTUAL", "10")
            .output()
            .unwrap();

        assert!(
            output.status.success(),
            "merge.sh must run when ORIGINAL_BAM and REFERENCE are passed \
             positionally holding a space and a $, exactly how \
             validate_pipeline.sh invokes it:\nstdout: {}\nstderr: {}",
            String::from_utf8_lossy(&output.stdout),
            String::from_utf8_lossy(&output.stderr)
        );
        let argv = argv_lines(&argv_out);
        for want in [&original, &ref_path] {
            assert!(
                argv.iter().any(|a| std::path::Path::new(a) == want.as_path()),
                "samtools must be handed the positional {} unmangled; got {:?}",
                want.display(),
                argv
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
    fn test_validate_pipeline_counts_only_the_non_advisory_checks() {
        // The harness's one `spike validate` canary is step 5's
        // `passed == 0`, and `spike validate` exits 1 into a `|| true`, so that
        // count is the whole signal. Two advisory rows (`resistant`,
        // `depth_fold`) are read back from the truth VCF and never touch the
        // BAM, so they pass whatever the BAM is: on a valid, indexed,
        // header-only merged.bam -- what a `merge.sh` that silently produced
        // nothing leaves -- a one-DEL truth VCF yields exactly the nine rows
        // below, two of them passing. Counting every row made `passed` 2 and
        // the guard went quiet; counting the non-advisory rows makes it 0
        // again.
        let dir = scratch_dir("validate_pipeline_advisory");
        let report = dir.join("spike_validate.json");
        std::fs::write(
            &report,
            r#"{
  "summary": { "total": 9, "pass": 2, "fail": 7, "counted_total": 5, "counted_pass": 0, "counted_fail": 5, "strict": false },
  "checks": [
    { "event": "DEL chr20:1000-2000 (x)", "check": "coverage_ratio", "expected": "0.50", "observed": "0.00", "pass": false, "advisory": false },
    { "event": "DEL chr20:1000-2000 (x)", "check": "coverage_any_mapq", "expected": "0.50", "observed": "0.00", "pass": false, "advisory": true },
    { "event": "DEL chr20:1000-2000 (x)", "check": "split_reads", "expected": ">=2", "observed": "0", "pass": false, "advisory": false },
    { "event": "DEL chr20:1000-2000 (x)", "check": "split_reads_each_end", "expected": ">=2 at each end", "observed": "0/0", "pass": false, "advisory": true },
    { "event": "DEL chr20:1000-2000 (x)", "check": "resistant", "expected": "<=0.100", "observed": "0.010", "pass": true, "advisory": true },
    { "event": "DEL chr20:1000-2000 (x)", "check": "depth_fold", "expected": "<=1.50", "observed": "1.11", "pass": true, "advisory": true },
    { "event": "[global]", "check": "insert_size", "expected": "mean 50-1000, sd 5-300", "observed": "no reads", "pass": false, "advisory": false },
    { "event": "[global]", "check": "dup_rate", "expected": "<50%", "observed": "no reads", "pass": false, "advisory": false },
    { "event": "[global]", "check": "mean_mapq", "expected": ">20", "observed": "no reads", "pass": false, "advisory": false }
  ]
}
"#,
        )
        .unwrap();

        // An older spike's report: no `advisory` key anywhere, so every row is
        // counted exactly as it always was.
        let old = dir.join("old_spike_validate.json");
        std::fs::write(
            &old,
            r#"{
  "summary": { "total": 3, "pass": 1, "fail": 2 },
  "checks": [
    { "event": "DEL chr20:1000-2000 (x)", "check": "coverage_ratio", "expected": "0.50", "observed": "0.50", "pass": true },
    { "event": "DEL chr20:1000-2000 (x)", "check": "split_reads", "expected": ">=2", "observed": "0", "pass": false },
    { "event": "[global]", "check": "mean_mapq", "expected": ">20", "observed": "no reads", "pass": false }
  ]
}
"#,
        )
        .unwrap();

        let output = std::process::Command::new("bash")
            .arg("-c")
            .arg(
                r#"script="$1"; report="$2"; old="$3"; shift 3; source "$script"
echo "advisory_only=$(spike_validate_counts "$report")"
read -r passed total adv <<< "$(spike_validate_counts "$report")"
if [[ "$total" -eq 0 || "$passed" -eq 0 ]]; then echo "guard=FIRES"; else echo "guard=quiet"; fi
echo "old_spike=$(spike_validate_counts "$old")"
echo "missing=[$(spike_validate_counts "$report.nope" || true)]"
printf '%s' '{"checks": [' > "$report.bad"
echo "truncated=[$(spike_validate_counts "$report.bad" || true)]""#,
            )
            .arg("_")
            .arg(validate_pipeline_script())
            .arg(&report)
            .arg(&old)
            .output()
            .unwrap();

        let stdout = String::from_utf8_lossy(&output.stdout);
        let stderr = String::from_utf8_lossy(&output.stderr);
        assert!(
            output.status.success(),
            "sourcing validate_pipeline.sh and calling spike_validate_counts \
             must work:\nstdout: {}\nstderr: {}",
            stdout,
            stderr
        );
        // Five non-advisory rows, none passing; four advisory rows beside them.
        assert!(
            stdout.contains("advisory_only=0 5 4"),
            "a report whose only passing rows are advisory must count 0 of 5 \
             passing:\n{}",
            stdout
        );
        assert!(
            stdout.contains("guard=FIRES"),
            "step 5's `passed == 0` guard must fire on that report:\n{}",
            stdout
        );
        assert!(
            stdout.contains("old_spike=1 3 0"),
            "an older spike's JSON has no `advisory` key and every row of it is \
             counted:\n{}",
            stdout
        );
        // A missing or truncated report prints nothing, which step 5 reads as
        // "no parseable JSON" -- never as a pass.
        assert!(
            stdout.contains("missing=[]") && stdout.contains("truncated=[]"),
            "a missing or unparseable report must yield no counts:\n{}",
            stdout
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

    /// Truvari's base-side output for one run: `tp-base.vcf.gz` and
    /// `fn.vcf.gz` under `dir`, each holding the named `sim_del_*` records.
    /// Written as plain text; the caller gzips them in its bash snippet.
    fn write_truvari_base(dir: &std::path::Path, tp: &[u32], fn_: &[u32]) {
        std::fs::create_dir_all(dir).unwrap();
        let vcf = |ids: &[u32]| {
            let mut s = String::from(
                "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
            );
            for id in ids {
                s.push_str(&format!(
                    "chr20\t{}\tsim_del_{}\tA\t<DEL>\t.\tPASS\tSVTYPE=DEL\n",
                    id * 1000,
                    id
                ));
            }
            s
        };
        std::fs::write(dir.join("tp-base.vcf"), vcf(tp)).unwrap();
        std::fs::write(dir.join("fn.vcf"), vcf(fn_)).unwrap();
    }

    #[test]
    fn test_validate_pipeline_recall_leaves_out_what_the_background_recovers() {
        // RF12: the truth DELs are common, so the background sample carries
        // several of them and the caller recovers those with no spike-in at
        // all. Recall is scored over the rest: an event the background control
        // recovers leaves the numerator *and* the denominator, whichever side
        // the spiked run put it on. Spiked TP {1,2,3}, FN {4,5}; control TP
        // {1,4} -> 3 events outside the background, 2 of them recovered.
        let dir = scratch_dir("validate_pipeline_outside_bg");
        write_truvari_base(&dir.join("spiked"), &[1, 2, 3], &[4, 5]);
        write_truvari_base(&dir.join("control"), &[1, 4], &[2, 3, 5]);
        // A control that recovers every truth event leaves nothing to score.
        write_truvari_base(&dir.join("control_all"), &[1, 2, 3, 4, 5], &[]);

        let output = std::process::Command::new("bash")
            .arg("-c")
            .arg(
                r#"script="$1"; d="$2"; shift 2; source "$script"
gzip -f "$d"/*/*.vcf
echo "scored=$(recall_outside_background "$d/spiked" "$d/control")"
echo "all_in_bg=$(recall_outside_background "$d/spiked" "$d/control_all")"
rc=0; out=$(recall_outside_background "$d/spiked" "$d/no_control") || rc=$?
echo "missing=[$out] rc=$rc""#,
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
            "sourcing validate_pipeline.sh and calling recall_outside_background \
             must work:\nstdout: {}\nstderr: {}",
            stdout,
            stderr
        );
        assert!(
            stdout.contains("scored=3 2 0.6667"),
            "events the control recovers must leave both sides of the recall:\n{}",
            stdout
        );
        assert!(
            stdout.contains("all_in_bg=0 0 N/A"),
            "nothing outside the background is no recall, not 0 or 1:\n{}",
            stdout
        );
        // A missing control is not an empty exclusion set: that would score
        // every background event as the spike-in's own.
        assert!(
            stdout.contains("missing=[] rc=1"),
            "a missing control must yield nothing and a non-zero status:\n{}",
            stdout
        );
    }

    #[test]
    fn test_validate_pipeline_summary_scores_recall_outside_the_background() {
        // Step 8 writes the recall outside the background beside the old
        // columns, and --min-recall judges that number. The spiked run's
        // overall recall (0.8, from summary.json) clears 0.7; its recall
        // outside the background (2 of 3) does not.
        let dir = scratch_dir("validate_pipeline_outside_bg_summary");
        let truth_vcf = "##fileformat=VCFv4.2\n\
                         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
                         chr20\t1000\tsim_del_1\tA\t<DEL>\t.\tPASS\tSVTYPE=DEL\n";
        let summary = |tp: u32, recall: &str| {
            format!(
                "{{\"TP-base\": {}, \"FP\": 0, \"FN\": 1, \"recall\": {}, \
                 \"precision\": 1.0, \"f1\": 0.9}}",
                tp, recall
            )
        };
        let spike_out = dir.join("spike_vaf_0.5");
        write_truvari_base(&spike_out.join("truvari"), &[1, 2, 3], &[4, 6]);
        std::fs::write(spike_out.join("truth.vcf"), truth_vcf).unwrap();
        std::fs::write(spike_out.join("truvari/summary.json"), summary(4, "0.8")).unwrap();
        let control = dir.join("background_control/truvari");
        write_truvari_base(&control, &[1, 4], &[2, 3, 6]);
        std::fs::write(control.join("summary.json"), summary(2, "0.4")).unwrap();

        let output = std::process::Command::new("bash")
            .arg("-c")
            .arg(
                r#"script="$1"; root="$2"; shift 2; source "$script"
gzip -f "$root"/spike_vaf_0.5/truvari/*.vcf "$root"/background_control/truvari/*.vcf
VAFS=(0.5)
MIN_GAIN=1
OUTDIR="$root"
for floor in 0.7 0.6; do
    MIN_RECALL=$floor
    FAILURES=()
    rc=0
    msg=$( { step8_summarize >/dev/null; step9_verdict; } 2>&1 ) || rc=$?
    echo "floor=$floor rc=$rc"
    printf '%s\n' "$msg" | sed "s/^/floor=$floor /"
done
cat "$root/validation_summary.tsv""#,
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
            "driving step 8 and 9 over a fabricated outdir must not error \
             out:\nstdout: {}\nstderr: {}",
            stdout,
            stderr
        );
        assert!(
            stdout.contains("\tN_outside_bg\tTP_outside_bg\tRecall_outside_bg"),
            "the summary header must name the three new columns:\n{}",
            stdout
        );
        // Spiked TP {1,2,3}, FN {4,6}; control TP {1,4}: 3 outside, 2 found.
        assert!(
            stdout.lines().any(|l| l.starts_with("0.5\t") && l.ends_with("\t3\t2\t0.6667")),
            "the VAF row must end with the recall outside the background:\n{}",
            stdout
        );
        assert!(
            stdout.lines().any(|l| l.starts_with("background\t") && l.ends_with("\tn/a\tn/a\tn/a")),
            "the background row has nothing outside itself to score:\n{}",
            stdout
        );
        assert!(
            stdout.contains("floor=0.7 rc=1") && stdout.contains("floor=0.6 rc=0"),
            "--min-recall must judge the recall outside the background (0.6667), \
             not the overall 0.8:\n{}",
            stdout
        );
    }

    #[test]
    fn test_validate_pipeline_records_each_tools_version() {
        // Truvari's matching flags have changed meaning between versions, so
        // a result is only reproducible with the tool versions beside it
        // (NF7). Three stand-ins: one answering --version, one answering
        // only `version` with noise on stderr (as bwa-mem2 does), one
        // answering neither.
        let dir = scratch_dir("validate_pipeline_versions");
        let stub = |name: &str, body: &str| {
            let path = dir.join(name);
            std::fs::write(&path, format!("#!/bin/sh\n{}\n", body)).unwrap();
            std::fs::set_permissions(&path, std::os::unix::fs::PermissionsExt::from_mode(0o755)).unwrap();
            path
        };
        let flag = stub("flag", r#"[ "$1" = --version ] && echo "flagtool 1.2" || exit 1"#);
        let word = stub("word", r#"echo "Looking to launch" >&2; [ "$1" = version ] && echo "2.2.1"; exit 0"#);
        let mute = stub("mute", "exit 1");
        let out = dir.join("tool_versions.tsv");

        let output = std::process::Command::new("bash")
            .arg("-c")
            .arg(r#"script="$1"; out="$2"; shift 2; source "$script"; record_tool_versions "$out""#)
            .arg("_")
            .arg(validate_pipeline_script())
            .arg(&out)
            .env("SPIKE", &flag)
            .env("SAMTOOLS", &flag)
            .env("BWAMEM2", &word)
            .env("BCFTOOLS", &flag)
            .env("DELLY", &mute)
            .env("TRUVARI", &word)
            .env("BGZIP", &flag)
            .env("TABIX", &flag)
            .output()
            .unwrap();

        assert!(
            output.status.success(),
            "stdout: {}\nstderr: {}",
            String::from_utf8_lossy(&output.stdout),
            String::from_utf8_lossy(&output.stderr)
        );
        let text = std::fs::read_to_string(&out).unwrap_or_default();
        for want in [
            "spike\tflagtool 1.2",
            "bwa-mem2\t2.2.1",
            "delly\tunknown",
            "truvari\t2.2.1",
            "tabix\tflagtool 1.2",
        ] {
            assert!(text.lines().any(|l| l == want), "want the line {:?} in:\n{}", want, text);
        }
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
        // Half-open: one region's end equal to the next one's start is not an
        // overlap. The regions compared are replacement footprints, so the case
        // is two spans exactly 2 * FOOTPRINT_MARGIN apart, which puts their
        // footprints end to start. (Before CR1 this test used spans that touched
        // each other -- 100-200 and 200-300 -- whose footprints overlap almost
        // completely and which spike now rejects.)
        let events = vec![
            del("chr1", 100_000, 200_000),
            del("chr1", 200_000 + 2 * FOOTPRINT_MARGIN, 300_000),
        ];
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
    fn test_overlap_policy_rejects_intersecting_replacement_footprints() {
        // The review's reproduction: two 1 kb deletions 1 kb apart. Their spans
        // do not overlap, but each event replaces reads across its span grown by
        // HAP_FLANK (the haplotype's reference flank) and MAX_FRAGMENT_LEN (the
        // longest fragment a pair can span), so the two replacements land on top
        // of each other and cancel each other's depth.
        let margin = HAP_FLANK + crate::stats::MAX_FRAGMENT_LEN as u64;
        let events = vec![del("chrT", 10_000, 11_000), del("chrT", 12_000, 13_000)];
        let err = validate_event_overlaps(&events, false)
            .expect_err("events 1 kb apart replace reads over the same range")
            .to_string();
        // The message must name both events and both footprints.
        assert!(
            err.contains("events 1 and 2"),
            "the error must name both events:\n{err}"
        );
        for (start, end) in [(10_000u64, 11_000u64), (12_000, 13_000)] {
            let footprint = format!("{}-{}", start - margin, end + margin);
            assert!(
                err.contains(&footprint),
                "the error must name the footprint {footprint}:\n{err}"
            );
            // ... and the span, which is what the user actually typed. A
            // footprint appears nowhere in the input, so a 50-record --vcf
            // whose message names only footprints cannot be searched for the
            // offending pair without subtracting the margin by hand.
            let span = format!("{start}-{end}");
            assert!(
                err.contains(&span),
                "the error must name the span {span} the user typed:\n{err}"
            );
        }
    }

    #[test]
    fn test_overlap_policy_rejects_one_below_the_footprint_boundary() {
        // The reject side of the boundary that
        // `test_overlap_policy_allows_touching_boundaries` pins from the accept
        // side: one bp closer than 2 * FOOTPRINT_MARGIN and the two footprints
        // share a base, so the pair must be refused.
        let events = vec![
            del("chr1", 100_000, 200_000),
            del("chr1", 200_000 + 2 * FOOTPRINT_MARGIN - 1, 300_000),
        ];
        assert!(validate_event_overlaps(&events, false).is_err());
    }

    #[test]
    fn test_overlap_policy_allows_footprints_that_do_not_intersect() {
        // Far enough apart that neither event touches the other's replacement
        // footprint: still accepted, exactly as before CR1.
        let events = vec![del("chrT", 10_000, 11_000), del("chrT", 31_000, 32_000)];
        assert!(validate_event_overlaps(&events, false).is_ok());
    }

    #[test]
    fn test_overlap_policy_allow_flag_covers_footprints() {
        // --allow-overlap keeps its meaning: intersecting footprints are then a
        // warning, not an error.
        let events = vec![del("chrT", 10_000, 11_000), del("chrT", 12_000, 13_000)];
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

    /// `scripts/validate_pipeline.sh` clusters the truth VCF's DELs on
    /// `FOOTPRINT_GAP`, a hand-written copy of `2 * FOOTPRINT_MARGIN`: the shell
    /// cannot read the Rust constant, so nothing but this test stops the two
    /// drifting apart. If the constant grows, the script keeps records spike
    /// then rejects and the run dies in step 3 rather than in the filter that
    /// exists to prevent it.
    #[test]
    fn test_validate_pipeline_script_mirrors_the_footprint_margin() {
        let script = include_str!("../scripts/validate_pipeline.sh");
        let assignment = script
            .lines()
            .find_map(|line| line.trim().strip_prefix("FOOTPRINT_GAP="))
            .expect("validate_pipeline.sh must set FOOTPRINT_GAP");
        let literal = assignment.split_whitespace().next().unwrap_or(assignment);

        assert_eq!(
            literal.parse::<u64>().ok(),
            Some(2 * FOOTPRINT_MARGIN),
            "validate_pipeline.sh's FOOTPRINT_GAP ({}) must equal \
             2 * FOOTPRINT_MARGIN ({})",
            literal,
            2 * FOOTPRINT_MARGIN
        );
    }

    #[test]
    fn test_parse_region_on_a_contig_name_containing_colons() {
        // 525 of GRCh38's contigs are HLA alleles whose names contain ':',
        // so --region cannot be split at the first one (L13).
        let contigs = vec!["HLA-A*01:01:01:01".to_string()];
        let r = parse_region("HLA-A*01:01:01:01:1000-2000", &contigs).unwrap();
        assert_eq!(r.chrom, "HLA-A*01:01:01:01");
        assert_eq!((r.start, r.end), (999, 2000));
    }

    #[test]
    fn test_parse_region_without_contigs_splits_at_the_first_colon() {
        let r = parse_region("chr20:1000-2000", &[]).unwrap();
        assert_eq!(r.chrom, "chr20");
        assert_eq!((r.start, r.end), (999, 2000));
        assert!(parse_region("chr20", &[]).is_err());
    }

    // --- N5: a donor pool too small to simulate from must fail loudly ---

    /// Read length for the donor-pool fixtures: the fragment model keeps
    /// only insert sizes in [read_length, MAX_FRAGMENT_LEN], and every
    /// fixture pair's 400bp insert lies inside this one.
    const TEST_READ_LENGTH: usize = 150;

    /// One usable donor pair, 100 bp, at `start`.
    fn donor_pair(name: &str, start: u64) -> ReadPair {
        ReadPair {
            name: name.to_string(),
            seq1: vec![b'A'; 100],
            qual1: vec![b'!' + 35; 100],
            seq2: vec![b'T'; 100],
            qual2: vec![b'!' + 35; 100],
            ref_start: start,
            ref_end: start + 400,
            insert_size: 400,
            chrom: "chr20".to_string(),
        }
    }

    fn donor_pairs(n: usize) -> Vec<ReadPair> {
        (0..n)
            .map(|i| donor_pair(&format!("p{}", i), 1000 + i as u64))
            .collect()
    }

    #[test]
    fn test_empty_donor_pool_is_rejected_instead_of_simulated_from() {
        // An event in a zero-coverage region, an off-target panel BAM or a
        // mistyped --region all yield 0 extracted pairs. spike used to log
        // "Built read pool: 0 pairs" and carry on: the quality profile had
        // no observation in any bin, so every base got synth.rs's
        // last-resort constant Q20 byte, and the run exited 0 with a truth
        // VCF and 2 invented read pairs beside it.
        let err = match finish_donor_pool(
            Vec::new(),
            &del("chr20", 30_000_000, 30_010_000),
            &["chr20:29990000-30020000".to_string()],
            0,
            TEST_READ_LENGTH,
        ) {
            Ok(pool) => panic!(
                "empty donor pool accepted; the run would write a truth VCF off {} pairs",
                pool.pairs.len()
            ),
            Err(e) => e.to_string(),
        };
        // The user has to be able to see which event starved, and where it
        // looked: a multi-event run reports only the one that failed.
        assert!(
            err.contains("DEL  chr20:30000001-30010000 (10000bp)"),
            "the error must name the event: {}",
            err
        );
        assert!(
            err.contains("chr20:29990000-30020000"),
            "the error must name the window it searched: {}",
            err
        );
    }

    #[test]
    fn test_donor_pool_below_the_minimum_is_rejected() {
        // "Empty" is not the whole defect: one pair trains no bin of the
        // quality profile either, and its single insert size becomes the
        // whole fragment distribution.
        let err = match finish_donor_pool(
            donor_pairs(MIN_DONOR_PAIRS - 1),
            &del("chr20", 30_000_000, 30_010_000),
            &["chr20:29990000-30020000".to_string()],
            0,
            TEST_READ_LENGTH,
        ) {
            Ok(_) => panic!("{} donor pairs accepted", MIN_DONOR_PAIRS - 1),
            Err(e) => e.to_string(),
        };
        assert!(
            err.contains(&format!("{}", MIN_DONOR_PAIRS - 1)),
            "the error must say how many pairs it found: {}",
            err
        );
    }

    #[test]
    fn test_donor_pool_at_the_minimum_is_accepted() {
        // The floor is a floor, not a coverage requirement: a pool that
        // reaches it is simulated from, unchanged.
        let pool = finish_donor_pool(
            donor_pairs(MIN_DONOR_PAIRS),
            &del("chr20", 30_000_000, 30_010_000),
            &["chr20:29990000-30020000".to_string()],
            0,
            TEST_READ_LENGTH,
        )
        .expect("a pool at the minimum must be usable");
        assert_eq!(pool.pairs.len(), MIN_DONOR_PAIRS);
    }

    #[test]
    fn test_donor_pool_is_counted_after_deduplication() {
        // Overlapping windows can hand the same fragment in twice; the
        // duplicate is not extra donor material and must not count towards
        // the floor (the check runs after dedup_pairs_by_name).
        let mut pairs = donor_pairs(MIN_DONOR_PAIRS - 1);
        pairs.push(donor_pair("p0", 1000)); // same name, second window
        assert_eq!(pairs.len(), MIN_DONOR_PAIRS);
        assert!(
            finish_donor_pool(
                pairs,
                &del("chr20", 30_000_000, 30_010_000),
                &["chr20:29990000-30020000".to_string()],
                0,
                TEST_READ_LENGTH,
            )
            .is_err(),
            "a duplicated fragment must not lift a pool over the floor"
        );
    }

    #[test]
    fn test_donor_pool_error_names_the_quality_drops_and_does_not_contradict_itself() {
        // Two problems in one message. (1) A CRAM that stores qualities via
        // read features has every record dropped by `quality_is_missing`, so
        // the pool is 0 and this guard fires -- the right outcome -- but the
        // message lists three causes and none of them is quality. (2) It says
        // "has no usable donor reads: 12 read pair(s) extracted", which
        // contradicts itself whenever the pool is non-empty.
        let err = match finish_donor_pool(
            donor_pairs(12),
            &del("chr20", 30_000_000, 30_010_000),
            &["chr20:29990000-30020000".to_string()],
            40,
            TEST_READ_LENGTH,
        ) {
            Ok(_) => panic!("12 donor pairs accepted"),
            Err(e) => e.to_string(),
        };
        assert!(
            !err.contains("has no usable donor reads"),
            "12 pairs were extracted, so \"no usable donor reads\" is false: {}",
            err
        );
        assert!(
            err.contains("40") && err.contains("qualit"),
            "the error must say how many read pairs were dropped for unusable \
             quality, since that is one way the pool empties: {}",
            err
        );
        // `dropped_unusable_qual` counts `UnusableQualTally::pair_names()`,
        // which is deduplicated across mates, so the number is pairs. Calling
        // it records here and "read pair(s)" in the run README made one number
        // read as two different quantities.
        assert!(
            err.contains("40 read pair(s) in those windows"),
            "the drop count is a pair count everywhere else, so it must say so \
             here too: {}",
            err
        );
    }

    #[test]
    fn test_run_readme_reports_the_pairs_removed_without_a_replacement() {
        // M14 drops unusable-quality pairs from the pool and adds their names
        // to replaced_reads.txt so merge.sh deletes them, but nothing in any
        // output file distinguishes them from replaced pairs -- the only
        // surface is a stderr INFO line. A localized depth dip the run README
        // does not mention is a CNV caller's signal.
        let dir = std::env::temp_dir().join(format!("spike_readme_{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();

        let events = vec![del("chr20", 38_412_500, 38_422_500)];
        let stats = vec![EventStat {
            vaf: 0.5,
            adjusted_vaf: None,
            kept: 4000,
            chimeric: 300,
            suppressed: 500,
            dropped_unusable_qual: 457,
            uncovered_breakpoint_sides: vec!["chr20:38409999".to_string()],
            census: census::Census::default(),
            depth_fold: types::DepthFold::default(),
        }];
        write_readme(
            dir.to_str().unwrap(),
            "spike -b x.bam",
            "x.bam",
            "ref.fa",
            &events,
            &stats,
            4300,
            10_000,
            412,
        )
        .unwrap();
        let md = std::fs::read_to_string(dir.join("README.md")).unwrap();
        let _ = std::fs::remove_dir_all(&dir);

        assert!(
            md.contains("457"),
            "the per-event count of pairs dropped for unusable quality must be \
             in the run README:\n{}",
            md
        );
        assert!(
            md.contains("412"),
            "the total removed without a replacement must be in the run README:\n{}",
            md
        );
    }

    #[test]
    fn test_run_readme_names_a_breakpoint_side_with_no_donor_coverage() {
        // `99f1a8e` keeps an event whose footprint is only partly inside the
        // donor data, which is right, but it kept it silently: the tiled
        // fragments over the bare side are scaled by depth measured at the
        // other one and land where the input BAM has no read. Measured on the
        // chr20 37.5-41.5 Mb slice, `del:chr20:37400000-37510000` puts 258 of
        // its 516 synthetic records over chr20:37,398,000-37,400,000, taking
        // that window from 0x in the input to 19.4x in the output. Nothing in
        // any output file said so.
        let dir = std::env::temp_dir().join(format!("spike_readme_bare_{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();

        let events = vec![del("chr20", 37_400_000, 37_510_000)];
        let stats = vec![EventStat {
            vaf: 0.5,
            adjusted_vaf: None,
            kept: 3000,
            chimeric: 258,
            suppressed: 66,
            dropped_unusable_qual: 0,
            uncovered_breakpoint_sides: vec!["chr20:37397999".to_string()],
            census: census::Census::default(),
            depth_fold: types::DepthFold::default(),
        }];
        write_readme(
            dir.to_str().unwrap(), "spike -b x.bam", "x.bam", "ref.fa",
            &events, &stats, 3258, 10_000, 0,
        )
        .unwrap();
        let md = std::fs::read_to_string(dir.join("README.md")).unwrap();
        let _ = std::fs::remove_dir_all(&dir);

        assert!(
            md.contains("chr20:37397999"),
            "the run README must name the breakpoint side with no donor \
             coverage:\n{}",
            md
        );
    }

    #[test]
    fn test_run_readme_says_nothing_when_every_side_is_covered() {
        // The paragraph may not appear on ordinary input, or it trains the
        // reader to skip it.
        let dir = std::env::temp_dir().join(format!("spike_readme_ok_{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();

        let events = vec![del("chr20", 38_412_500, 38_422_500)];
        let stats = vec![EventStat {
            vaf: 0.5,
            adjusted_vaf: None,
            kept: 4000,
            chimeric: 300,
            suppressed: 500,
            dropped_unusable_qual: 0,
            uncovered_breakpoint_sides: Vec::new(),
            census: census::Census::default(),
            depth_fold: types::DepthFold::default(),
        }];
        write_readme(
            dir.to_str().unwrap(), "spike -b x.bam", "x.bam", "ref.fa",
            &events, &stats, 4300, 10_000, 0,
        )
        .unwrap();
        let md = std::fs::read_to_string(dir.join("README.md")).unwrap();
        let _ = std::fs::remove_dir_all(&dir);

        assert!(
            !md.contains("no donor coverage"),
            "a fully covered run must not carry the warning:\n{}",
            md
        );
    }

    #[test]
    fn test_run_readme_prints_the_resistant_reads_per_event() {
        // CR4: the reads spike could not edit, per event, beside what it did.
        // The paragraph only appears when an event is above the threshold.
        let dir =
            std::env::temp_dir().join(format!("spike_readme_resist_{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();

        let events = vec![
            del("chr20", 38_412_500, 38_422_500),
            del("chr20", 39_000_000, 39_010_000),
        ];
        let stat = |counted, resistant| EventStat {
            vaf: 1.0,
            adjusted_vaf: None,
            kept: 0,
            chimeric: 300,
            suppressed: 500,
            dropped_unusable_qual: 0,
            uncovered_breakpoint_sides: Vec::new(),
            census: census::Census { counted, resistant },
            depth_fold: types::DepthFold::default(),
        };
        let stats = vec![stat(100, 50), stat(200, 4)];
        write_readme(
            dir.to_str().unwrap(), "spike -b x.bam", "x.bam", "ref.fa",
            &events, &stats, 800, 10_000, 0,
        )
        .unwrap();
        let md = std::fs::read_to_string(dir.join("README.md")).unwrap();
        let _ = std::fs::remove_dir_all(&dir);

        assert!(md.contains("| Resistant reads |"), "no column:\n{}", md);
        assert!(md.contains("| 50 of 100 (50%) |"), "first row:\n{}", md);
        assert!(md.contains("| 4 of 200 (2%) |"), "second row:\n{}", md);
        let paragraph = md
            .lines()
            .find(|l| l.starts_with("**Reads spike could not edit:**"))
            .unwrap_or_else(|| panic!("no paragraph for the event above 10%:\n{}", md));
        assert!(paragraph.contains("chr20:38412501"), "{}", paragraph);
        assert!(!paragraph.contains("chr20:39000001"), "{}", paragraph);
    }

    #[test]
    fn test_run_readme_prints_the_depth_fold_per_event() {
        // CR2: the fold per event, and a paragraph naming only an event above
        // 1.5, with the bin and the two depths.
        let dir = std::env::temp_dir().join(format!("spike_readme_fold_{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();

        let events = vec![
            del("chr20", 38_412_500, 38_422_500),
            del("chr20", 39_000_000, 39_010_000),
        ];
        let stat = |fold: f64, bin: &str, depth: f64| EventStat {
            vaf: 0.5,
            adjusted_vaf: None,
            kept: 0,
            chimeric: 300,
            suppressed: 500,
            dropped_unusable_qual: 0,
            uncovered_breakpoint_sides: Vec::new(),
            census: census::Census::default(),
            depth_fold: types::DepthFold {
                fold,
                scaled_by: 40.0,
                worst_bin: bin.to_string(),
                worst_depth: depth,
            },
        };
        let stats = vec![stat(3.73, "chr20:38410000-38411000", 10.0), stat(1.12, "", 40.0)];
        write_readme(
            dir.to_str().unwrap(), "spike -b x.bam", "x.bam", "ref.fa",
            &events, &stats, 600, 10_000, 0,
        )
        .unwrap();
        let md = std::fs::read_to_string(dir.join("README.md")).unwrap();
        let _ = std::fs::remove_dir_all(&dir);

        assert!(md.contains("| Depth fold |"), "no column:\n{}", md);
        assert!(md.contains("| 3.73 |"), "first row:\n{}", md);
        assert!(md.contains("| 1.12 |"), "second row:\n{}", md);
        let paragraph = md
            .lines()
            .find(|l| l.starts_with("**Donor depth off the scaling depth:**"))
            .unwrap_or_else(|| panic!("no paragraph for the event above 1.5:\n{}", md));
        assert!(paragraph.contains("chr20:38410000-38411000"), "{}", paragraph);
        assert!(paragraph.contains("10.0x"), "{}", paragraph);
        assert!(!paragraph.contains("chr20:39000001"), "{}", paragraph);
    }

    #[test]
    fn test_run_readme_prints_the_simulated_fraction_beside_the_request() {
        // `4f06842` made truth.vcf's SIM_VAF the fraction that was actually
        // simulated and moved the request into SIM_REQ_VAF, but the run
        // README's VAF column kept printing the request alone: a capped run's
        // README said 0.990 where its own truth.vcf said SIM_VAF=0.950, and a
        // floored one said 0.030 where truth.vcf said 0.058. Two files from
        // the same run disagreed about what the reads carry.
        let dir = std::env::temp_dir().join(format!(
            "spike_readme_vaf_{}_{}",
            std::process::id(),
            "floored"
        ));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();

        let events = vec![
            del("chr20", 38_412_500, 38_422_500),
            del("chr20", 39_000_000, 39_010_000),
        ];
        let stats = vec![
            // The two-fragment floor planted 0.058 where 0.030 was asked for.
            EventStat {
                vaf: 0.03,
                adjusted_vaf: Some(0.058),
                kept: 53,
                chimeric: 2,
                suppressed: 1,
                dropped_unusable_qual: 0,
                uncovered_breakpoint_sides: Vec::new(),
                census: census::Census::default(),
                depth_fold: types::DepthFold::default(),
            },
            // Neither mechanism touched this one, so both numbers are the
            // request -- truth.vcf writes SIM_VAF=SIM_REQ_VAF=0.500 for it.
            EventStat {
                vaf: 0.5,
                adjusted_vaf: None,
                kept: 4000,
                chimeric: 300,
                suppressed: 500,
                dropped_unusable_qual: 0,
                uncovered_breakpoint_sides: Vec::new(),
                census: census::Census::default(),
                depth_fold: types::DepthFold::default(),
            },
        ];
        write_readme(
            dir.to_str().unwrap(), "spike -b x.bam", "x.bam", "ref.fa",
            &events, &stats, 4355, 10_000, 0,
        )
        .unwrap();
        let md = std::fs::read_to_string(dir.join("README.md")).unwrap();
        let _ = std::fs::remove_dir_all(&dir);

        assert!(
            md.contains("| 0.030 | 0.058 |"),
            "the floored event's row must carry the request and the fraction \
             that was simulated, in that order, so it agrees with truth.vcf's \
             SIM_REQ_VAF and SIM_VAF:\n{}",
            md
        );
        assert!(
            md.contains("| 0.500 | 0.500 |"),
            "an event neither the cap nor the floor touched must show the \
             request in both columns:\n{}",
            md
        );
        assert!(
            md.contains("| Requested VAF | Simulated VAF |"),
            "the table must say which number is the request and which is what \
             was simulated:\n{}",
            md
        );
    }

    #[test]
    fn test_donor_pool_error_names_a_fusion_by_both_breakpoints() {
        // A fusion draws from two windows; either side can be the starved
        // one, so the message carries both.
        let err = match finish_donor_pool(
            Vec::new(),
            &fusion("chr2", 42_000_000, "chr2", 29_000_000),
            &[
                "chr2:41990000-42010000".to_string(),
                "chr2:28990000-29010000".to_string(),
            ],
            0,
            TEST_READ_LENGTH,
        ) {
            Ok(_) => panic!("empty fusion donor pool accepted"),
            Err(e) => e.to_string(),
        };
        assert!(
            err.contains("FUSION  chr2:42000001>>chr2:29000001"),
            "the error must name both breakpoints: {}",
            err
        );
        assert!(
            err.contains("chr2:41990000-42010000") && err.contains("chr2:28990000-29010000"),
            "the error must name both windows: {}",
            err
        );
    }

    #[test]
    fn test_build_haplotype_stores_generated_insertion_sequence() {
        let reference = crate::reference::SharedReference::from_sequences(
            [("chr1".to_string(), b"ACGT".repeat(100))].into(),
        );
        let mut event = SimEvent::Insertion {
            chrom: "chr1".to_string(),
            pos: 200,
            ins_seq: None,
            ins_len: 30,
            gene: "G".to_string(),
            allele_fraction: None,
        };
        let mut rng = StdRng::seed_from_u64(7);
        let hap = build_haplotype(&mut event, &reference, 50, "junction", &mut rng).unwrap();
        let stored = match &event {
            SimEvent::Insertion { ins_seq, .. } => ins_seq.clone(),
            _ => unreachable!(),
        };
        let stored = stored.expect("generated sequence kept on the event for truth to write");
        // The bases kept are the ones the haplotype -- and so every read cut
        // from it -- carries, not just the same count of them.
        assert_eq!(stored.len(), 30);
        assert_eq!(hap.get_sequence(50, 30), stored.as_slice());
    }

    #[test]
    fn test_validate_edit_model_accepts_clean_and_origin() {
        assert!(validate_edit_model("clean").is_ok());
        assert!(validate_edit_model("origin").is_ok());
    }

    #[test]
    fn test_validate_edit_model_rejects_anything_else_and_names_both() {
        let err = validate_edit_model("Origin").unwrap_err().to_string();
        assert!(err.contains("'Origin'") && err.contains("clean") && err.contains("origin"), "{}", err);
    }

    #[test]
    fn test_edit_model_defaults_to_clean() {
        let args = Args::try_parse_from(["spike", "--bam", "x.bam", "--reference", "x.fa"]).unwrap();
        assert_eq!(args.edit_model, "clean");
    }

    #[test]
    fn test_origin_footprint_is_the_haplotypes_reference_range() {
        let del = |start, end| SimEvent::Deletion {
            chrom: "chr1".into(),
            del_start: start,
            del_end: end,
            gene: "G".into(),
            exons: vec![],
            allele_fraction: None,
        };
        let span = |s, e| Some(origin::Span::new("chr1", s, e));
        assert_eq!(
            origin_footprint(&del(5000, 6000), 100_000),
            span(3000, 8000)
        );
        // The flanks stop at the contig's ends, as `VariantHaplotype::from_deletion` does.
        assert_eq!(
            origin_footprint(&del(99_000, 99_500), 100_000),
            span(97_000, 100_000)
        );
        assert_eq!(origin_footprint(&del(500, 600), 100_000), span(0, 2600));
    }
}
