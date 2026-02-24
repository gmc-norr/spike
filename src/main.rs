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
    ///   --event "fusion:GENEA:exon14:GENEB:exon2"       (fusion, requires --exon-bed)
    ///   --event "dup:chr20:30000000-30005000"
    ///   --event "inv:chr20:30000000-30005000"
    ///   --event "ins:chr20:30000000:500"
    ///   --event "snp:chr20:30000000:A:T"                (SNP/small variant)
    ///   --event "snp:chr20:30000000:A>T"                (alternate syntax)
    ///   --event "snp:chr20:30000000:ACG:A"              (small deletion)
    ///   --event "snp:chr20:30000000:A:ACGT"             (small insertion)
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

    /// Automatically run bwa-mem2 alignment after FASTQ generation.
    #[arg(long)]
    align: bool,

    /// Optional gVCF/VCF with SNP calls for LOH simulation.
    /// When provided, het SNP positions are extracted from this file to
    /// determine which reads belong to the deleted haplotype (more accurate
    /// than the default pileup-based approach). Supports .vcf and .vcf.gz
    /// (requires bcftools in PATH for .vcf.gz).
    #[arg(long)]
    gvcf: Option<String>,
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
    // Convert from 1-based inclusive to 0-based half-open.
    Ok(ExtractionRegion {
        chrom: chrom.to_string(),
        start: start.saturating_sub(1),
        end,
    })
}

fn main() -> Result<()> {
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

    // Compute BAM stats for read length.
    let bam_stats = crate::bam_stats::compute_stats(&args.bam, 50_000)?;
    let read_length = bam_stats.read_length.round() as usize;

    let config = SimConfig {
        bam_path: args.bam.clone(),
        ref_path: args.reference.clone(),
        allele_fraction: args.allele_fraction,
        flank_bp: args.flank,
        read_length,
        min_mapq: args.min_mapq,
        gvcf_path: args.gvcf.clone(),
    };

    // Parse --region if provided.
    let extraction_region = if let Some(ref region_str) = args.region {
        let r = parse_region(region_str)?;
        log::info!(
            "Extraction region: {}:{}-{} (0-based half-open)",
            r.chrom, r.start, r.end,
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
    let shared_ref =
        crate::reference::SharedReference::load(&config.ref_path, &chrom_refs)?;

    let mut all_output_pairs: Vec<ReadPair> = Vec::new();

    // Haplotype flank: must be >= max expected fragment length so reads near
    // haplotype edges form complete pairs. 2000bp is conservative.
    let hap_flank: u64 = 2000;

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
        )?;

        // Build quality profile and synth generator.
        let quality_profile =
            synth::QualityProfile::from_read_pairs(&pool.pairs, config.read_length);
        let mut synth_gen =
            synth::SynthReadGenerator::new(quality_profile, &shared_ref, config.read_length);

        // Build variant haplotype.
        let haplotype = build_haplotype(event, &shared_ref, hap_flank, &mut rng)?;

        // Simulate: suppress reads + tile synthetic reads across haplotype.
        let output = simulate::simulate_event(
            event,
            &pool,
            &haplotype,
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
    write_align_script(&args.output, &args.reference, args.threads)?;

    // Summary.
    log::info!("=== spike complete ===");
    log::info!("Output directory: {}", args.output);
    log::info!("FASTQ: {} and {}", r1_path, r2_path);
    log::info!("Truth VCF: {}", truth_path.display());
    log::info!("Total read pairs: {}", all_output_pairs.len());

    // Auto-align if requested.
    if args.align {
        run_alignment(&args.output, &args.reference, args.threads)?;
    } else {
        log::info!(
            "To align: bash {}/align.sh",
            args.output,
        );
    }

    Ok(())
}

/// Deduplicate read pairs by name, keeping the first occurrence.
fn dedup_by_name(pairs: &mut Vec<ReadPair>) {
    let mut seen = std::collections::HashSet::new();
    pairs.retain(|p| seen.insert(p.name.clone()));
}

/// Write a convenience shell script for bwa-mem2 alignment.
fn write_align_script(output_dir: &str, ref_path: &str, threads: usize) -> Result<()> {
    let script_path = Path::new(output_dir).join("align.sh");
    let script = format!(
        r#"#!/bin/bash
set -euo pipefail
# Align simulated reads with bwa-mem2 and sort.
# Usage: bash align.sh [REF] [THREADS]
REF="${{1:-{ref_path}}}"
THREADS="${{2:-{threads}}}"
DIR="$(cd "$(dirname "$0")" && pwd)"

echo "Aligning $DIR/R1.fq.gz + R2.fq.gz with bwa-mem2 ($THREADS threads)..."
pixi run bwa-mem2 mem -t "$THREADS" \
    -R '@RG\tID:sim\tSM:SIM\tPL:ILLUMINA' \
    "$REF" \
    "$DIR/R1.fq.gz" "$DIR/R2.fq.gz" \
    2>"$DIR/bwamem2.log" | \
    pixi run samtools sort -@ 4 -o "$DIR/sim.bam" -

pixi run samtools index "$DIR/sim.bam"

TOTAL=$(pixi run samtools view -c "$DIR/sim.bam" 2>/dev/null || echo "?")
SA_COUNT=$(pixi run samtools view "$DIR/sim.bam" | grep -c "SA:Z:" 2>/dev/null || echo "0")
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

/// Run bwa-mem2 alignment.
fn run_alignment(output_dir: &str, ref_path: &str, threads: usize) -> Result<()> {
    log::info!("Running bwa-mem2 alignment...");

    let script_path = Path::new(output_dir).join("align.sh");
    let status = std::process::Command::new("bash")
        .arg(&script_path)
        .arg(ref_path)
        .arg(threads.to_string())
        .status()?;

    if !status.success() {
        bail!("bwa-mem2 alignment failed with exit code {:?}", status.code());
    }

    let bam_path = Path::new(output_dir).join("sim.bam");
    log::info!("Aligned BAM: {}", bam_path.display());
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
        )?;
        let pairs_b = extract::extract_read_pairs(
            &config.bam_path,
            chrom_b,
            region_b_start,
            region_b_end,
            config.min_mapq,
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
        } => VariantHaplotype::from_duplication(reference, chrom, *dup_start, *dup_end, flank),
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
                (0..*ins_len)
                    .map(|_| bases[rng.gen_range(0..4)])
                    .collect()
            };
            VariantHaplotype::from_insertion(reference, chrom, *pos, &seq, flank)
        }
        SimEvent::Fusion {
            chrom_a,
            bp_a,
            chrom_b,
            bp_b,
            ..
        } => VariantHaplotype::from_fusion(reference, chrom_a, *bp_a, chrom_b, *bp_b, flank),
        SimEvent::SmallVariant {
            chrom,
            pos,
            ref_allele,
            alt_allele,
            ..
        } => VariantHaplotype::from_small_variant(reference, chrom, *pos, ref_allele, alt_allele, flank),
    }
}
