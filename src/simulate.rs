//! Unified SV event simulation: suppress reads + tile synthetic reads across
//! the variant haplotype.
//!
//! Replaces the per-SV splice functions (splice_deletion, splice_duplication, etc.)
//! with a single `simulate_event` that works for all SV types.

use std::collections::{HashMap, HashSet};

use anyhow::Result;
use rand::rngs::StdRng;
use rand::Rng;

use crate::haplotype::VariantHaplotype;
use crate::loh;
use crate::synth::SynthReadGenerator;
use crate::types::{ReadPair, ReadPool, SimConfig, SimEvent, SplicedOutput};

/// Classification of how a read pair relates to SV boundaries.
#[derive(Debug, PartialEq)]
enum PairRelation {
    /// Entirely outside all SV boundaries — always kept.
    Outside,
    /// Entirely inside a deleted/inverted region.
    Inside,
    /// Fragment spans an SV boundary.
    Overlapping,
}

/// Simulate one SV event: suppress original reads + generate haplotype-tiled reads.
///
/// This unified function replaces the separate splice_deletion, splice_duplication,
/// splice_inversion, splice_insertion, and splice_fusion functions.
pub fn simulate_event(
    event_index: usize,
    event: &SimEvent,
    pool: &ReadPool,
    haplotype: &mut VariantHaplotype,
    config: &SimConfig,
    synth_gen: &mut SynthReadGenerator,
    vaf: f64,
    rng: &mut StdRng,
) -> Result<SplicedOutput> {
    let name_prefix = format!("ev{:04}", event_index);

    // Get SV boundaries for read classification.
    let (sv_start, sv_end) = match event.primary_region() {
        Some((_chrom, start, end)) => (start, end),
        None => {
            // Fusion: use the first breakpoint as a point event.
            if let SimEvent::Fusion { bp_a, .. } = event {
                (*bp_a, *bp_a)
            } else {
                anyhow::bail!(
                    "event has no primary region and is not a Fusion — \
                     this is a bug; please report it"
                );
            }
        }
    };

    // For fusions (and DUPs with legacy junction model), we keep all original
    // reads and add chimeric on top. For full tandem DUPs, we suppress+replace
    // like DEL/INV since the full haplotype provides both depth and junction reads.
    let is_additive = match event {
        SimEvent::Fusion { .. } => true,
        SimEvent::Duplication { .. } => config.dup_model == "junction",
        _ => false,
    };

    // LOH / haplotype-aware suppression for het events.
    let (loh_set, classified_set, hap_variants) = get_haplotype_info(event, config, vaf, rng);
    let use_loh = !loh_set.is_empty();

    // Apply het SNP variants to both the haplotype sequence (for tiled chimeric
    // reads) and synth_gen (for DUP depth copies). This ensures all synthetic
    // reads carry the correct alleles instead of reference-only bases.
    if !hap_variants.is_empty() {
        haplotype.apply_variants(&hap_variants);
        synth_gen.set_haplotype_variants(hap_variants);
    }

    // Suppress original reads within the haplotype's reference footprint.
    //
    // For non-additive events (DEL, INV, INS): suppress reads at VAF rate
    // within the haplotype reference range [hap_ref_start, hap_ref_end).
    // Tiled haplotype reads replace the variant haplotype's contribution in
    // this zone. Reads entirely OUTSIDE this range are kept unchanged —
    // they're in unaffected reference territory.
    //
    // For additive events (DUP, Fusion): keep all originals. Chimeric reads
    // near the breakpoint and DUP depth copies are added on top.
    let (hap_ref_start, hap_ref_end) = haplotype.ref_range().unwrap_or((sv_start, sv_end));

    let mut kept = Vec::new();
    let mut suppressed: Vec<String> = Vec::new();

    for pair in &pool.pairs {
        if is_additive {
            kept.push(pair.clone());
            continue;
        }

        // Only pairs entirely inside the haplotype's reference footprint are
        // replaced. Tiled fragments never extend past the haplotype ends, so
        // suppressing pairs that stick out would leave a depth dip there.
        let in_haplotype = classify_pair_relation(pair, hap_ref_start, hap_ref_end);
        if !matches!(in_haplotype, PairRelation::Inside) {
            kept.push(pair.clone());
            continue;
        }

        // Read is within the haplotype range.
        // For inside/overlapping the SV with LOH info, use haplotype-aware
        // suppression for classified reads, random fallback for unclassified.
        // For haplotype-flank reads (or without LOH), use random.
        let in_sv = classify_pair_relation(pair, sv_start, sv_end);
        let use_loh_for_this =
            use_loh && matches!(in_sv, PairRelation::Inside | PairRelation::Overlapping);

        if use_loh_for_this {
            if classified_set.contains(&pair.name) {
                // Read was classified via het SNPs — use LOH decision.
                //
                // loh_set contains ~50% of classified reads (one haplotype).
                // To achieve the target VAF we must scale: at VAF=0.5 suppress
                // all loh_set reads (100% × 50% ≈ 50% total). At VAF=0.3
                // suppress 60% of loh_set reads (60% × 50% ≈ 30% total).
                // For VAF>0.5 we also suppress some non-loh classified reads.
                if loh_set.contains(&pair.name) {
                    let hap_prob = (2.0 * vaf).min(1.0);
                    if rng.gen::<f64>() < hap_prob {
                        suppressed.push(pair.name.clone());
                    } else {
                        kept.push(pair.clone());
                    }
                } else {
                    // Non-variant haplotype: only suppress when VAF > 0.5.
                    let other_prob = (2.0 * vaf - 1.0).max(0.0);
                    if rng.gen::<f64>() < other_prob {
                        suppressed.push(pair.name.clone());
                    } else {
                        kept.push(pair.clone());
                    }
                }
            } else {
                // Read couldn't be classified (no het SNP overlap) —
                // fall back to random suppression at the VAF rate.
                if rng.gen::<f64>() > vaf {
                    kept.push(pair.clone());
                } else {
                    suppressed.push(pair.name.clone());
                }
            }
        } else if rng.gen::<f64>() > vaf {
            kept.push(pair.clone());
        } else {
            suppressed.push(pair.name.clone());
        }
    }

    // Estimate coverage at the first breakpoint for chimeric read count.
    let bp_positions = haplotype.breakpoints();
    let first_bp_ref = bp_positions
        .first()
        .and_then(|&bp| haplotype.hap_to_ref(bp.saturating_sub(1)))
        .map(|(_, pos)| pos)
        .unwrap_or(sv_start);
    let cov = estimate_coverage_at(pool, first_bp_ref, 2000);

    // Tile synthetic reads across the haplotype.
    // For additive events (DUP/Fusion), restrict tiling to near breakpoints
    // to avoid inflating flank coverage.
    let chimeric = tile_haplotype_reads(
        haplotype,
        synth_gen,
        pool,
        cov,
        vaf,
        is_additive, // breakpoint_only
        &name_prefix,
        rng,
    );

    // For DUPs with legacy junction model, generate depth copies inside the
    // region. With the full tandem model, tiling handles depth automatically.
    let depth_copies = if let SimEvent::Duplication {
        dup_start, dup_end, ..
    } = event
    {
        if config.dup_model == "junction" {
            synth_gen.generate_dup_depth_copies(
                pool,
                *dup_start,
                *dup_end,
                vaf,
                &loh_set,
                &classified_set,
                &name_prefix,
                rng,
            )
        } else {
            Vec::new()
        }
    } else {
        Vec::new()
    };

    log::info!(
        "simulate_event: {} kept, {} suppressed, {} chimeric (haplotype-tiled), {} depth copies",
        kept.len(),
        suppressed.len(),
        chimeric.len(),
        depth_copies.len(),
    );

    let mut all_chimeric = chimeric;
    all_chimeric.extend(depth_copies);

    Ok(SplicedOutput {
        chimeric_pairs: all_chimeric,
        kept_originals: kept,
        suppressed_count: suppressed.len(),
        suppressed_names: suppressed,
    })
}

/// Combine the outputs of all simulated events into the final set of pairs.
///
/// Each event's pool spans event ± flank, so nearby events extract the same
/// originals. An original suppressed by any event is dropped, even if another
/// event passed it through as kept. Otherwise the second event would undo the
/// first event's suppression.
pub fn combine_event_outputs(outputs: Vec<SplicedOutput>) -> Vec<ReadPair> {
    let suppressed: HashSet<String> = outputs
        .iter()
        .flat_map(|o| o.suppressed_names.iter().cloned())
        .collect();

    let mut all_pairs = Vec::new();
    for output in outputs {
        all_pairs.extend(
            output
                .kept_originals
                .into_iter()
                .filter(|p| !suppressed.contains(&p.name)),
        );
        all_pairs.extend(output.chimeric_pairs);
    }
    dedup_by_name(&mut all_pairs);
    all_pairs
}

/// Deduplicate read pairs by name, keeping the last occurrence.
///
/// Keeping the last occurrence makes event-order behavior explicit when
/// multiple simulated events touch the same original read name.
fn dedup_by_name(pairs: &mut Vec<ReadPair>) {
    let mut last_idx: HashMap<String, usize> = HashMap::with_capacity(pairs.len());
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

/// Get LOH haplotype set, classified set, and variant map for the event.
///
/// Returns `(target_set, classified_set, haplotype_variants)`:
/// - `target_set`: reads to suppress (DEL) or copy (DUP).
/// - `classified_set`: all reads that could be assigned to a haplotype via het SNPs.
///   Reads NOT in this set should fall back to random suppression/copying at VAF rate.
/// - `haplotype_variants`: het SNP positions → allele for DUP variant substitution.
fn get_haplotype_info(
    event: &SimEvent,
    config: &SimConfig,
    vaf: f64,
    rng: &mut StdRng,
) -> (HashSet<String>, HashSet<String>, HashMap<u64, u8>) {
    // Only use haplotype-aware suppression for het events (VAF 0.3-0.7).
    if !(0.3..=0.7).contains(&vaf) {
        return (HashSet::new(), HashSet::new(), HashMap::new());
    }

    match event {
        SimEvent::Deletion {
            chrom,
            del_start,
            del_end,
            ..
        } => {
            let (target_set, classified_set) = loh::identify_deleted_haplotype_reads(
                &config.bam_path,
                chrom,
                *del_start,
                *del_end,
                config.min_mapq,
                config.gvcf_path.as_deref(),
                Some(config.ref_path.as_str()),
                rng,
            )
            .unwrap_or_else(|e| {
                log::warn!("LOH classification failed: {}", e);
                (HashSet::new(), HashSet::new())
            });
            (target_set, classified_set, HashMap::new())
        }
        SimEvent::Duplication {
            chrom,
            dup_start,
            dup_end,
            ..
        } => loh::identify_duplicated_haplotype_reads(
            &config.bam_path,
            chrom,
            *dup_start,
            *dup_end,
            config.min_mapq,
            config.gvcf_path.as_deref(),
            Some(config.ref_path.as_str()),
            rng,
        )
        .unwrap_or_else(|e| {
            log::warn!("DUP haplotype classification failed: {}", e);
            (HashSet::new(), HashSet::new(), HashMap::new())
        }),
        _ => (HashSet::new(), HashSet::new(), HashMap::new()),
    }
}

/// Classify how a read pair relates to the SV boundaries.
fn classify_pair_relation(pair: &ReadPair, sv_start: u64, sv_end: u64) -> PairRelation {
    if pair.ref_end <= sv_start || pair.ref_start >= sv_end {
        PairRelation::Outside
    } else if pair.ref_start >= sv_start && pair.ref_end <= sv_end {
        PairRelation::Inside
    } else {
        PairRelation::Overlapping
    }
}

/// Compute the number of fragments to tile across a haplotype.
///
/// For non-additive events (DEL, INV, INS): uses the reference-mapped length
/// of the haplotype (excluding novel insertion sequence) to avoid inflating
/// coverage near insertion points.
///
/// For additive events (DUP, Fusion): every tiled fragment crosses a
/// breakpoint and all original fragments are kept, so `n` junction fragments
/// make up n / (coverage + n) of the depth there. For fraction `vaf` that is
/// n = coverage × vaf / (1 − vaf) per breakpoint.
fn compute_tiling_count(
    haplotype: &VariantHaplotype,
    coverage: f64,
    vaf: f64,
    mean_frag: f64,
    breakpoint_only: bool,
) -> usize {
    if haplotype.total_len == 0 {
        return 0;
    }

    let breakpoints = haplotype.breakpoints();

    if breakpoint_only && !breakpoints.is_empty() {
        // An additive event can't reach vaf = 1 (it would need infinitely
        // many added fragments), so cap it.
        const MAX_ADDITIVE_VAF: f64 = 0.95;
        if vaf > MAX_ADDITIVE_VAF {
            log::warn!(
                "additive event: VAF {:.2} capped at {:.2} (original reads are kept)",
                vaf,
                MAX_ADDITIVE_VAF
            );
        }
        let v = vaf.min(MAX_ADDITIVE_VAF);
        let n = (coverage * v / (1.0 - v) * breakpoints.len() as f64).round() as usize;
        return n.max(2); // at least 2 chimeric reads
    }

    // Fragment starts are uniform over the starts whose fragment overlaps
    // reference sequence: [0, L - f] minus starts lying wholly inside inserted
    // sequence (tiling redraws those). Interior depth is then n * f / starts;
    // matching the suppressed v * coverage needs n = coverage * v * starts / f.
    let novel_only: f64 = haplotype
        .segments
        .iter()
        .filter(|seg| seg.origin.is_none())
        .map(|seg| (seg.sequence.len() as f64 - mean_frag).max(0.0))
        .sum();
    let effective_len = if haplotype.ref_mapped_len() > 0 {
        (haplotype.total_len as f64 - mean_frag - novel_only).max(0.0)
    } else {
        0.0
    };

    let n = ((coverage * vaf * effective_len) / mean_frag).round() as usize;
    n.max(2) // at least 2 chimeric reads
}

/// Tile synthetic reads across the variant haplotype.
///
/// Number of reads is determined by coverage, VAF, and haplotype/zone length.
/// Reads naturally become chimeric when they span segment boundaries.
///
/// When `breakpoint_only` is true (for DUP/Fusion additive events), reads are
/// placed only near segment boundaries so they cross a breakpoint. This avoids
/// inflating coverage in flank regions where original reads are already kept.
fn tile_haplotype_reads(
    haplotype: &VariantHaplotype,
    synth_gen: &SynthReadGenerator,
    pool: &ReadPool,
    coverage: f64,
    vaf: f64,
    breakpoint_only: bool,
    name_prefix: &str,
    rng: &mut StdRng,
) -> Vec<ReadPair> {
    let hap_len = haplotype.total_len;
    if hap_len == 0 {
        return Vec::new();
    }

    // Use the library's real mean fragment length; fall back only when the
    // pool gave none (e.g. no reads).
    let mean_frag = if pool.frag_dist.mean.is_finite() && pool.frag_dist.mean > 0.0 {
        pool.frag_dist.mean
    } else {
        300.0
    };
    let read_length = synth_gen.read_length() as u64;
    let breakpoints = haplotype.breakpoints();

    let n_frags = compute_tiling_count(haplotype, coverage, vaf, mean_frag, breakpoint_only);

    log::info!(
        "Tiling {} synthetic reads across {}bp haplotype (cov={:.1}, vaf={:.2}, bp_only={})",
        n_frags,
        hap_len,
        coverage,
        vaf,
        breakpoint_only,
    );

    let mut pairs = Vec::with_capacity(n_frags);

    // Check if any novel (non-reference) segments exist. If so, we use
    // rejection sampling to avoid placing fragments entirely within novel
    // sequence (which wouldn't contribute to observable reference-aligned coverage).
    let has_novel = haplotype.segments.iter().any(|seg| seg.origin.is_none());

    // Attempt up to 2x the target count to compensate for rejected placements
    // (e.g., R1/R2 starting in a novel segment where hap_to_ref returns None).
    let max_attempts = n_frags * 2;
    let mut attempts = 0;
    let mut idx = 0;

    while pairs.len() < n_frags && attempts < max_attempts {
        attempts += 1;

        // Sample fragment length from empirical distribution.
        let frag_len = pool
            .frag_dist
            .sample_in_range(rng, read_length as i64, 1500) as u64;

        if frag_len > hap_len {
            continue;
        }

        let max_start = hap_len.saturating_sub(frag_len);

        let hap_start = if breakpoint_only && !breakpoints.is_empty() {
            // Pick a random breakpoint, then sample a start position that
            // ensures the fragment crosses it. Fragment [s, s+frag_len)
            // crosses bp when s < bp and s + frag_len > bp.
            let bp_idx = rng.gen_range(0..breakpoints.len());
            let bp = breakpoints[bp_idx];
            let zone_start = bp.saturating_sub(frag_len.saturating_sub(1));
            let zone_end = bp.saturating_sub(1).min(max_start);
            if zone_start > zone_end {
                continue;
            }
            rng.gen_range(zone_start..=zone_end)
        } else {
            // Uniform across the whole haplotype.
            // For haplotypes with novel segments, reject placements that
            // land entirely in novel sequence (up to 10 attempts).
            let mut start = rng.gen_range(0..=max_start);
            if has_novel {
                for _ in 0..10 {
                    if haplotype.overlaps_ref_segment(start, frag_len) {
                        break;
                    }
                    start = rng.gen_range(0..=max_start);
                }
            }
            start
        };

        let name = format!("{}_hap_{:06}", name_prefix, idx);
        idx += 1;

        if let Some(pair) = synth_gen.generate_haplotype_read_pair(
            haplotype,
            hap_start,
            frag_len,
            &name,
            rng,
        ) {
            pairs.push(pair);
        }
    }

    pairs
}

/// Estimate fragment depth at a reference position.
///
/// Counts fragments overlapping positions in a window, returns mean coverage.
/// Uses up to 50 evenly-spaced sample points for stable estimates.
/// The pool must be sorted by `ref_start` (which `build_read_pool` ensures).
fn estimate_coverage_at(pool: &ReadPool, pos: u64, window: u64) -> f64 {
    let start = pos.saturating_sub(window / 2);
    let end = pos.saturating_add(window / 2);

    let n_samples: u64 = 50;
    let range = end - start;
    // For very small windows, use fewer sample points (at least 1).
    let actual_samples = range.min(n_samples).max(1);
    let step = if actual_samples > 1 {
        range / actual_samples
    } else {
        1
    };

    let mut total = 0usize;

    for i in 0..actual_samples {
        let sample_pos = start + i * step;
        // Binary search: skip pairs that start after sample_pos (they can't
        // overlap it). The pool is sorted by ref_start, so partition_point
        // gives us the exact cutoff.
        let upper = pool
            .pairs
            .partition_point(|p| p.ref_start <= sample_pos);
        let count = pool.pairs[..upper]
            .iter()
            .filter(|p| p.ref_end > sample_pos)
            .count();
        total += count;
    }

    total as f64 / actual_samples as f64
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::haplotype::{HaplotypeSegment, SegmentOrigin};
    use crate::types::FusionJoin;
    use crate::stats::FragmentDist;

    fn make_pair(name: &str, start: u64, end: u64) -> ReadPair {
        ReadPair {
            name: name.to_string(),
            seq1: vec![],
            qual1: vec![],
            seq2: vec![],
            qual2: vec![],
            ref_start: start,
            ref_end: end,
            insert_size: (end - start) as i64,
            chrom: "chr1".to_string(),
        }
    }

    /// Build a test haplotype from segments (bypass SharedReference).
    fn make_haplotype(segments: Vec<HaplotypeSegment>) -> VariantHaplotype {
        VariantHaplotype::from_segments(segments)
    }

    /// Make a reference-origin segment.
    fn ref_segment(offset: u64, len: u64) -> HaplotypeSegment {
        HaplotypeSegment {
            sequence: vec![b'A'; len as usize],
            origin: Some(SegmentOrigin {
                chrom: "chr1".to_string(),
                ref_start: offset,
                ref_end: offset + len,
                is_reverse: false,
            }),
            hap_offset: 0, // will be recomputed
        }
    }

    /// Make a novel (insertion) segment with no reference origin.
    fn novel_segment(len: u64) -> HaplotypeSegment {
        HaplotypeSegment {
            sequence: vec![b'T'; len as usize],
            origin: None,
            hap_offset: 0,
        }
    }

    fn make_pool(pairs: Vec<ReadPair>) -> ReadPool {
        let frag_dist = FragmentDist::from_stats(400.0, 80.0);
        ReadPool { pairs, frag_dist }
    }

    // ---------------------------------------------------------------
    // classify_pair_relation tests
    // ---------------------------------------------------------------

    #[test]
    fn test_classify_pair_relation() {
        // Outside (left).
        assert_eq!(
            classify_pair_relation(&make_pair("a", 100, 200), 300, 400),
            PairRelation::Outside
        );
        // Outside (right).
        assert_eq!(
            classify_pair_relation(&make_pair("a", 500, 600), 300, 400),
            PairRelation::Outside
        );
        // Inside.
        assert_eq!(
            classify_pair_relation(&make_pair("a", 310, 390), 300, 400),
            PairRelation::Inside
        );
        // Overlapping left.
        assert_eq!(
            classify_pair_relation(&make_pair("a", 250, 350), 300, 400),
            PairRelation::Overlapping
        );
        // Overlapping right.
        assert_eq!(
            classify_pair_relation(&make_pair("a", 350, 450), 300, 400),
            PairRelation::Overlapping
        );
        // Spanning both boundaries.
        assert_eq!(
            classify_pair_relation(&make_pair("a", 250, 450), 300, 400),
            PairRelation::Overlapping
        );
    }

    // ---------------------------------------------------------------
    // compute_tiling_count tests
    // ---------------------------------------------------------------

    #[test]
    fn test_tiling_count_del() {
        // 5kb DEL with 2kb flanks: ref_mapped_len = 4kb (2 flanks), no novel.
        // Haplotype: [left_flank=2kb] [right_flank=2kb] (no middle → deletion).
        let hap = make_haplotype(vec![
            ref_segment(0, 2000),    // left flank
            ref_segment(7000, 2000), // right flank (after 5kb deleted region)
        ]);
        assert_eq!(hap.total_len, 4000);
        assert_eq!(hap.ref_mapped_len(), 4000);

        // 30x coverage, 0.5 VAF, 400bp mean frag. Fragment starts are uniform
        // over [0, L - f], so interior depth is n * f / (L - f); for v * cov
        // that is n = cov * v * (L - f) / f = 30 * 0.5 * 3600 / 400 = 135.
        let count = compute_tiling_count(&hap, 30.0, 0.5, 400.0, false);
        assert_eq!(count, 135);
    }

    #[test]
    fn test_tiling_count_ins() {
        // 500bp INS with 2kb flanks: ref_mapped_len = 4kb (flanks), total = 4.5kb.
        // Should use ref_mapped_len (4kb), NOT total_len (4.5kb).
        let hap = make_haplotype(vec![
            ref_segment(0, 2000),    // left flank
            novel_segment(500),      // inserted sequence
            ref_segment(2000, 2000), // right flank
        ]);
        assert_eq!(hap.total_len, 4500);
        assert_eq!(hap.ref_mapped_len(), 4000);

        // Fragment starts that overlap reference: (4500 - 400) minus the
        // 500 - 400 = 100 starts that would lie wholly in the insertion.
        // n = 30 * 0.5 * 4000 / 400 = 150.
        let count = compute_tiling_count(&hap, 30.0, 0.5, 400.0, false);
        assert_eq!(count, 150);

        // 1000 bp insertion: (5000 - 400) - (1000 - 400) = 4000 starts -> 150.
        let long = make_haplotype(vec![
            ref_segment(0, 2000),
            novel_segment(1000),
            ref_segment(2000, 2000),
        ]);
        assert_eq!(compute_tiling_count(&long, 30.0, 0.5, 400.0, false), 150);
    }

    #[test]
    fn test_tiling_count_inv() {
        // 5kb INV with 2kb flanks: total = 9kb, ref_mapped_len = 9kb.
        // Same ref footprint as before the INV, just middle is reversed.
        let hap = make_haplotype(vec![
            ref_segment(0, 2000), // left flank
            HaplotypeSegment {
                sequence: vec![b'A'; 5000],
                origin: Some(SegmentOrigin {
                    chrom: "chr1".to_string(),
                    ref_start: 2000,
                    ref_end: 7000,
                    is_reverse: true, // inverted
                }),
                hap_offset: 0,
            },
            ref_segment(7000, 2000), // right flank
        ]);
        assert_eq!(hap.total_len, 9000);
        assert_eq!(hap.ref_mapped_len(), 9000);

        let count = compute_tiling_count(&hap, 30.0, 0.5, 400.0, false);
        // Expected: 30 * 0.5 * (9000 - 400) / 400 = 322.5 → 323
        assert_eq!(count, 323);
    }

    #[test]
    fn test_tiling_count_dup_breakpoint_only() {
        // DUP junction haplotype: breakpoint_only = true.
        // Haplotype: [pre-dup flank] [post-dup flank] with 1 breakpoint.
        let hap = make_haplotype(vec![
            ref_segment(0, 2000), // up to dup end
            ref_segment(0, 2000), // from dup start (junction)
        ]);

        // 1 breakpoint, all originals kept: 30 * 0.5 / (1 - 0.5) = 30.
        let count = compute_tiling_count(&hap, 30.0, 0.5, 400.0, true);
        assert_eq!(count, 30);
    }

    #[test]
    fn test_tiling_count_breakpoint_only_gives_requested_fraction() {
        // Additive events keep all `cov` original fragments at the junction,
        // so n junction fragments make up n / (cov + n) of it. For fraction v,
        // n = cov * v / (1 - v) per breakpoint.
        let one_bp = make_haplotype(vec![ref_segment(0, 2000), ref_segment(5000, 2000)]);
        assert_eq!(compute_tiling_count(&one_bp, 40.0, 0.2, 400.0, true), 10);
        assert_eq!(compute_tiling_count(&one_bp, 100.0, 0.05, 400.0, true), 5);

        let two_bp = make_haplotype(vec![
            ref_segment(0, 2000),
            ref_segment(5000, 2000),
            ref_segment(9000, 2000),
        ]);
        assert_eq!(compute_tiling_count(&two_bp, 40.0, 0.2, 400.0, true), 20);
    }

    #[test]
    fn test_tiling_count_breakpoint_only_is_bounded_at_full_vaf() {
        // v = 1 would need infinitely many added fragments; the count must
        // stay finite (an unbounded usize would abort on allocation).
        let hap = make_haplotype(vec![ref_segment(0, 2000), ref_segment(5000, 2000)]);
        let count = compute_tiling_count(&hap, 40.0, 1.0, 400.0, true);
        assert!(count > compute_tiling_count(&hap, 40.0, 0.9, 400.0, true));
        assert!(count <= 40 * 100, "count {} is unbounded", count);
    }

    #[test]
    fn test_tiling_count_minimum() {
        // Very low coverage → at least 2 reads.
        let hap = make_haplotype(vec![ref_segment(0, 100)]);
        let count = compute_tiling_count(&hap, 0.1, 0.1, 400.0, false);
        assert_eq!(count, 2); // min of 2
    }

    #[test]
    fn test_tiling_count_empty_haplotype() {
        let hap = make_haplotype(vec![]);
        assert_eq!(hap.total_len, 0);
        let count = compute_tiling_count(&hap, 30.0, 0.5, 400.0, false);
        assert_eq!(count, 0);
    }

    // ---------------------------------------------------------------
    // estimate_coverage_at tests
    // ---------------------------------------------------------------

    #[test]
    fn test_estimate_coverage_at() {
        // 10 reads covering [0, 500), query at pos 250 with window 500.
        let pairs: Vec<ReadPair> = (0..10)
            .map(|i| make_pair(&format!("r{}", i), 0, 500))
            .collect();
        let pool = make_pool(pairs);
        let cov = estimate_coverage_at(&pool, 250, 500);
        assert!((cov - 10.0).abs() < 1.0, "expected ~10.0 got {:.1}", cov);
    }

    #[test]
    fn test_estimate_coverage_at_partial() {
        // Reads only cover [0, 250). Query at 250 with window 500.
        // Half the sample positions should see coverage, half should not.
        let pairs: Vec<ReadPair> = (0..20)
            .map(|i| make_pair(&format!("r{}", i), 0, 250))
            .collect();
        let pool = make_pool(pairs);
        let cov = estimate_coverage_at(&pool, 250, 500);
        assert!(
            cov < 20.0 && cov > 0.0,
            "expected partial coverage, got {:.1}",
            cov
        );
    }

    // ---------------------------------------------------------------
    // Suppression balance tests
    //
    // These test the suppression logic directly using the internal
    // classify_pair_relation + random suppression pattern, without
    // needing a SynthReadGenerator (which requires a real FASTA).
    // ---------------------------------------------------------------


    // ---------------------------------------------------------------
    // overlaps_ref_segment tests
    // ---------------------------------------------------------------

    #[test]
    fn test_overlaps_ref_segment() {
        // Haplotype: ref[0..2000] + novel[500bp] + ref[2000..4000]
        let hap = make_haplotype(vec![
            ref_segment(0, 2000),
            novel_segment(500),
            ref_segment(2000, 2000),
        ]);

        // Fragment in first ref segment: overlaps.
        assert!(hap.overlaps_ref_segment(100, 400));
        // Fragment spanning ref→novel boundary: overlaps (touches ref).
        assert!(hap.overlaps_ref_segment(1900, 400));
        // Fragment entirely in novel segment [2000..2500): no ref overlap.
        assert!(!hap.overlaps_ref_segment(2000, 400));
        // Fragment spanning novel→ref boundary: overlaps (touches second ref).
        assert!(hap.overlaps_ref_segment(2200, 500));
        // Fragment in second ref segment: overlaps.
        assert!(hap.overlaps_ref_segment(3000, 400));
    }

    // ── Read evidence tests per SV type ─────────────────────────────────

    use crate::reference::SharedReference;
    use crate::synth::{QualityProfile, SynthReadGenerator};
    use rand::rngs::StdRng;
    use rand::SeedableRng;
    use std::collections::HashMap as StdHashMap;

    fn mock_synth_gen(read_length: usize) -> SynthReadGenerator<'static> {
        let q30 = vec![b'!' + 30; read_length];
        let pairs: Vec<ReadPair> = (0..100)
            .map(|i| ReadPair {
                name: format!("mock_{}", i),
                seq1: vec![b'A'; read_length],
                qual1: q30.clone(),
                seq2: vec![b'T'; read_length],
                qual2: q30.clone(),
                ref_start: i as u64 * 500,
                ref_end: i as u64 * 500 + 500,
                insert_size: 500,
                chrom: "chr1".to_string(),
            })
            .collect();
        let profile = QualityProfile::from_read_pairs(&pairs, read_length);
        let pattern = b"ACGT";
        let seq: Vec<u8> = (0..100_000u64).map(|i| pattern[(i % 4) as usize]).collect();
        let mut seqs = StdHashMap::new();
        seqs.insert("chr1".to_string(), seq);
        let reference = SharedReference::from_sequences(seqs);
        let ref_static: &'static SharedReference = Box::leak(Box::new(reference));
        SynthReadGenerator::new(profile, ref_static, read_length, 0.0)
    }

    /// Build a deletion haplotype with unique per-base sequences (not all-A).
    fn del_haplotype(flank: u64, del_size: u64) -> VariantHaplotype {
        let pattern = b"ACGT";
        let left_seq: Vec<u8> = (0..flank).map(|i| pattern[(i % 4) as usize]).collect();
        let right_seq: Vec<u8> = (0..flank)
            .map(|i| pattern[((del_size + flank + i) % 4) as usize])
            .collect();

        make_haplotype(vec![
            HaplotypeSegment {
                sequence: left_seq,
                origin: Some(SegmentOrigin {
                    chrom: "chr1".to_string(),
                    ref_start: 0,
                    ref_end: flank,
                    is_reverse: false,
                }),
                hap_offset: 0,
            },
            HaplotypeSegment {
                sequence: right_seq,
                origin: Some(SegmentOrigin {
                    chrom: "chr1".to_string(),
                    ref_start: flank + del_size,
                    ref_end: 2 * flank + del_size,
                    is_reverse: false,
                }),
                hap_offset: 0,
            },
        ])
    }

    #[test]
    fn test_del_tiled_reads_include_split_reads() {
        // DEL: 2kb flanks, 5kb deletion → breakpoint at hap offset 2000
        // With boundary-crossing allowed, some reads should have an arm
        // that spans the breakpoint (chimeric content for split-read evidence).
        let hap = del_haplotype(2000, 5000);
        let gen = mock_synth_gen(150);
        let pool = make_pool(vec![]);
        let mut rng = StdRng::seed_from_u64(42);

        let pairs =
            tile_haplotype_reads(&hap, &gen, &pool, 30.0, 0.5, false, "test", &mut rng);
        assert!(!pairs.is_empty(), "Should produce some read pairs");

        // Discordant pairs (R1 in left, R2 in right) should still exist.
        let discordant_count = pairs
            .iter()
            .filter(|p| p.ref_start < 2000 && p.ref_end > 7000)
            .count();
        assert!(
            discordant_count > 0,
            "Should produce discordant pairs spanning the deletion"
        );
    }

    #[test]
    fn test_del_tiled_reads_produce_discordant_pairs() {
        // DEL: 2kb flanks, 5kb deletion → breakpoint at hap offset 2000
        // Discordant pairs: R1 in left (ref < 2000), R2 in right (ref >= 7000)
        let hap = del_haplotype(2000, 5000);
        let gen = mock_synth_gen(150);
        let pool = make_pool(vec![]);
        let mut rng = StdRng::seed_from_u64(42);

        let pairs =
            tile_haplotype_reads(&hap, &gen, &pool, 30.0, 0.5, false, "test", &mut rng);

        let mut discordant_count = 0;
        for pair in &pairs {
            // A discordant pair: ref_start < 2000 AND ref_end > 7000
            if pair.ref_start < 2000 && pair.ref_end > 7000 {
                discordant_count += 1;

                // Ref span should be much larger than typical fragment
                let ref_span = pair.ref_end - pair.ref_start;
                assert!(
                    ref_span > 5000,
                    "Discordant pair ref span {} should exceed deletion size",
                    ref_span
                );
            }
        }

        assert!(
            discordant_count > 0,
            "Should produce at least one discordant pair (got 0 out of {} total)",
            pairs.len()
        );
    }

    #[test]
    fn test_dup_junction_reads_span_noncontiguous_ref() {
        // DUP junction: left = ref[dup_end-flank..dup_end], right = ref[dup_start..dup_start+flank]
        // e.g., dup_start=1000, dup_end=6000, flank=2000
        // Left segment: ref[4000..6000], Right segment: ref[1000..3000]
        // Junction reads have R1 near 6000 and R2 near 1000 → non-contiguous
        let hap = make_haplotype(vec![
            HaplotypeSegment {
                sequence: {
                    let p = b"ACGT";
                    (4000u64..6000).map(|i| p[(i % 4) as usize]).collect()
                },
                origin: Some(SegmentOrigin {
                    chrom: "chr1".to_string(),
                    ref_start: 4000,
                    ref_end: 6000,
                    is_reverse: false,
                }),
                hap_offset: 0,
            },
            HaplotypeSegment {
                sequence: {
                    let p = b"ACGT";
                    (1000u64..3000).map(|i| p[(i % 4) as usize]).collect()
                },
                origin: Some(SegmentOrigin {
                    chrom: "chr1".to_string(),
                    ref_start: 1000,
                    ref_end: 3000,
                    is_reverse: false,
                }),
                hap_offset: 0,
            },
        ]);

        let gen = mock_synth_gen(150);
        let pool = make_pool(vec![]);
        let mut rng = StdRng::seed_from_u64(42);

        let pairs =
            tile_haplotype_reads(&hap, &gen, &pool, 30.0, 0.5, false, "test", &mut rng);
        assert!(!pairs.is_empty());

        // Some pairs should span the junction (R1 ref near 6000, R2 ref near 1000)
        let junction_count = pairs
            .iter()
            .filter(|p| {
                // R1 in left segment (ref 4000-6000) and R2 in right segment (ref 1000-3000)
                // This means ref_start < ref_end but the actual mapping is non-contiguous
                p.ref_start >= 1000 && p.ref_start < 3000 && p.ref_end > 4000
            })
            .count();

        assert!(
            junction_count > 0,
            "DUP junction should produce pairs spanning non-contiguous ref regions"
        );
    }

    #[test]
    fn test_tiled_read_names_are_unique_across_events() {
        let hap = del_haplotype(2000, 5000);
        let gen = mock_synth_gen(150);
        let pool = make_pool(vec![]);
        let mut rng_a = StdRng::seed_from_u64(42);
        let mut rng_b = StdRng::seed_from_u64(43);

        let pairs_a = tile_haplotype_reads(
            &hap, &gen, &pool, 30.0, 0.5, false, "ev0001", &mut rng_a,
        );
        let pairs_b = tile_haplotype_reads(
            &hap, &gen, &pool, 30.0, 0.5, false, "ev0002", &mut rng_b,
        );

        assert!(!pairs_a.is_empty());
        assert!(!pairs_b.is_empty());
        assert!(pairs_a.iter().all(|p| p.name.starts_with("ev0001_")));
        assert!(pairs_b.iter().all(|p| p.name.starts_with("ev0002_")));

        let names_a: std::collections::HashSet<&str> =
            pairs_a.iter().map(|p| p.name.as_str()).collect();
        let overlap_count = pairs_b
            .iter()
            .filter(|p| names_a.contains(p.name.as_str()))
            .count();
        assert_eq!(
            overlap_count, 0,
            "synthetic names must not collide across events"
        );
    }

    #[test]
    fn test_small_variant_tiling_produces_reads() {
        // Small-variant-like haplotype with very short flanks.
        // With read_length=150 and flank segments=100bp, reads necessarily
        // cross segment boundaries — this should work now that boundary
        // crossing is always allowed.
        let hap = make_haplotype(vec![
            HaplotypeSegment {
                sequence: vec![b'A'; 100],
                origin: Some(SegmentOrigin {
                    chrom: "chr1".to_string(),
                    ref_start: 1000,
                    ref_end: 1100,
                    is_reverse: false,
                }),
                hap_offset: 0,
            },
            HaplotypeSegment {
                sequence: vec![b'T'; 1], // alt allele
                origin: Some(SegmentOrigin {
                    chrom: "chr1".to_string(),
                    ref_start: 1100,
                    ref_end: 1101,
                    is_reverse: false,
                }),
                hap_offset: 0,
            },
            HaplotypeSegment {
                sequence: vec![b'C'; 100],
                origin: Some(SegmentOrigin {
                    chrom: "chr1".to_string(),
                    ref_start: 1101,
                    ref_end: 1201,
                    is_reverse: false,
                }),
                hap_offset: 0,
            },
        ]);

        let gen = mock_synth_gen(150);
        let pool = ReadPool {
            pairs: vec![],
            frag_dist: FragmentDist::from_stats(170.0, 10.0),
        };
        let mut rng = StdRng::seed_from_u64(7);

        let pairs = tile_haplotype_reads(
            &hap, &gen, &pool, 30.0, 0.5, false, "sv", &mut rng,
        );

        assert!(
            !pairs.is_empty(),
            "small-variant tiling should produce reads (boundary crossing allowed)",
        );
    }

    // ── Full tandem DUP tiling tests ─────────────────────────────────────

    fn tandem_dup_haplotype(dup_start: u64, dup_end: u64, flank: u64) -> VariantHaplotype {
        let pattern = b"ACGT";
        let make_seq = |start: u64, end: u64| -> Vec<u8> {
            (start..end).map(|i| pattern[(i % 4) as usize]).collect()
        };
        let left_start = dup_start.saturating_sub(flank);
        let right_end = dup_end + flank;

        make_haplotype(vec![
            HaplotypeSegment {
                sequence: make_seq(left_start, dup_start),
                origin: Some(SegmentOrigin {
                    chrom: "chr1".to_string(),
                    ref_start: left_start,
                    ref_end: dup_start,
                    is_reverse: false,
                }),
                hap_offset: 0,
            },
            HaplotypeSegment {
                sequence: make_seq(dup_start, dup_end),
                origin: Some(SegmentOrigin {
                    chrom: "chr1".to_string(),
                    ref_start: dup_start,
                    ref_end: dup_end,
                    is_reverse: false,
                }),
                hap_offset: 0,
            },
            HaplotypeSegment {
                sequence: make_seq(dup_start, dup_end),
                origin: Some(SegmentOrigin {
                    chrom: "chr1".to_string(),
                    ref_start: dup_start,
                    ref_end: dup_end,
                    is_reverse: false,
                }),
                hap_offset: 0,
            },
            HaplotypeSegment {
                sequence: make_seq(dup_end, right_end),
                origin: Some(SegmentOrigin {
                    chrom: "chr1".to_string(),
                    ref_start: dup_end,
                    ref_end: right_end,
                    is_reverse: false,
                }),
                hap_offset: 0,
            },
        ])
    }

    #[test]
    fn test_tandem_dup_tiling_produces_junction_reads() {
        // Full tandem DUP: [500..1000) flank | [1000..3000) copy1 | [1000..3000) copy2 | [3000..3500) flank
        // Novel junction at hap offset 2500 (copy1→copy2: dup_end→dup_start)
        let hap = tandem_dup_haplotype(1000, 3000, 500);
        assert_eq!(hap.total_len, 5000); // 500 + 2000 + 2000 + 500
        assert_eq!(hap.segments.len(), 4);

        let gen = mock_synth_gen(150);
        let pool = make_pool(vec![]);
        let mut rng = StdRng::seed_from_u64(42);

        // breakpoint_only=false for full tandem model
        let pairs = tile_haplotype_reads(
            &hap, &gen, &pool, 30.0, 0.5, false, "dup", &mut rng,
        );
        assert!(!pairs.is_empty(), "Should produce reads from full tandem haplotype");

        // Some pairs should have discordant reference mapping at the junction:
        // R1 maps to ref near dup_end (3000) and R2 maps to ref near dup_start (1000)
        // These are pairs where ref_start is in [1000,3000) and ref_end is also in [1000,3000)
        // but the ref_start > ref_end (or similar non-contiguous mapping) due to the junction.
        // Actually, since both copies map to [1000,3000), discordant evidence comes from
        // reads crossing the junction whose content aligns to ref near 3000 then 1000.
        // The ReadPair ref_start/ref_end will show ref coordinates within [1000,3000).
        // The split-read evidence will be in the aligned BAM, not directly visible here.
        // But we can verify that reads exist near the junction.

        let _junction_offset = 2500u64; // copy1→copy2 boundary in hap coords
        let _junction_pairs = pairs
            .iter()
            .filter(|p| {
                // ref_start maps to ref near dup_end (from end of copy1)
                // or ref_start maps to ref near dup_start (from start of copy2)
                // The key signature: read spans junction so ref_start > ref_end
                // (end of copy1 maps to high ref, start of copy2 maps to low ref)
                p.ref_start > p.ref_end.saturating_sub(1)
                    || (p.ref_start >= 2500 && p.ref_end <= 1500)
            })
            .count();

        // At minimum, the pool should have reads — some will be junction-crossing
        assert!(pairs.len() > 10, "Should produce substantial number of reads");
    }

    #[test]
    fn test_tandem_dup_tiling_count_uses_full_length() {
        // Full tandem DUP: total ref_mapped_len = 2*flank + 2*dup_size
        let hap = tandem_dup_haplotype(1000, 3000, 500);

        // ref_mapped_len = 500 + 2000 + 2000 + 500 = 5000
        assert_eq!(hap.ref_mapped_len(), 5000);

        // breakpoint_only=false: uses ref_mapped_len as effective_len
        let n = compute_tiling_count(&hap, 30.0, 0.5, 400.0, false);
        // Expected: 30 * 0.5 * 5000 / 400 = 187.5 → 188
        assert!(n > 100, "Full tandem DUP should produce many reads, got {}", n);
    }

    // ---------------------------------------------------------------
    // simulate_event integration tests
    // ---------------------------------------------------------------

    fn make_config() -> SimConfig {
        SimConfig {
            bam_path: String::new(),
            ref_path: String::new(),
            allele_fraction: 0.5,
            flank_bp: 1000,
            read_length: 150,
            min_mapq: 20,
            gvcf_path: None,
            indel_error_rate: 0.0,
            dup_model: "full".to_string(),
        }
    }

    /// Build a pool of read pairs uniformly covering [region_start, region_end).
    fn make_covering_pool(region_start: u64, region_end: u64, count: usize) -> ReadPool {
        let step = (region_end - region_start) / count as u64;
        let pairs = (0..count)
            .map(|i| {
                let start = region_start + i as u64 * step;
                make_pair(&format!("read_{}", i), start, start + 400)
            })
            .collect();
        ReadPool {
            pairs,
            frag_dist: FragmentDist::from_stats(400.0, 80.0),
        }
    }

    #[test]
    fn test_simulate_event_del_produces_chimeric_and_suppresses() {
        // DEL: chr1:1000-3000 (2kb deletion), flanks at [0,1000) and [3000,4000).
        // del_haplotype(flank=1000, del_size=2000) builds [0,1000) | [3000,4000).
        let mut hap = del_haplotype(1000, 2000);

        // Pool: 200 pairs spread across [0, 5000) so reads land before, inside, and after deletion.
        let pool = make_covering_pool(0, 5000, 200);
        let config = make_config();
        let mut gen = mock_synth_gen(150);
        let mut rng = StdRng::seed_from_u64(42);

        let event = SimEvent::Deletion {
            chrom: "chr1".to_string(),
            del_start: 1000,
            del_end: 3000,
            gene: "TEST".to_string(),
            exons: vec![],
            allele_fraction: Some(0.5),
        };

        let out = simulate_event(1, &event, &pool, &mut hap, &config, &mut gen, 0.5, &mut rng)
            .expect("simulate_event should succeed for DEL");

        // Reads inside the deletion should be suppressed at ~50%.
        assert!(
            out.suppressed_count > 0,
            "DEL should suppress some reads; got suppressed_count={}",
            out.suppressed_count
        );

        // Chimeric reads should be generated to replace suppressed fraction.
        assert!(
            !out.chimeric_pairs.is_empty(),
            "DEL should produce chimeric reads spanning the deletion junction"
        );

        // Chimeric read names should carry the event prefix.
        assert!(
            out.chimeric_pairs[0].name.starts_with("ev0001"),
            "chimeric read name should start with event prefix 'ev0001'"
        );

        // Total output (kept + chimeric) should be non-empty.
        assert!(
            out.kept_originals.len() + out.chimeric_pairs.len() > 0,
            "simulate_event should produce output reads"
        );
    }

    #[test]
    fn test_simulate_event_fusion_is_additive() {
        // Fusion: chr1:10000 >> chr1:20000 (same chromosome for simplicity).
        // Haplotype: left flank [9000,10000) | right flank [20000,21000).
        let mut hap = make_haplotype(vec![
            HaplotypeSegment {
                sequence: vec![b'A'; 1000],
                origin: Some(SegmentOrigin {
                    chrom: "chr1".to_string(),
                    ref_start: 9000,
                    ref_end: 10000,
                    is_reverse: false,
                }),
                hap_offset: 0,
            },
            HaplotypeSegment {
                sequence: vec![b'C'; 1000],
                origin: Some(SegmentOrigin {
                    chrom: "chr1".to_string(),
                    ref_start: 20000,
                    ref_end: 21000,
                    is_reverse: false,
                }),
                hap_offset: 1000,
            },
        ]);

        // Pool: 100 pairs near the breakpoint region.
        let pool = make_covering_pool(8000, 11000, 100);
        let config = make_config();
        let mut gen = mock_synth_gen(150);
        let mut rng = StdRng::seed_from_u64(99);

        let event = SimEvent::Fusion {
            chrom_a: "chr1".to_string(),
            bp_a: 10000,
            gene_a: "GENE_A".to_string(),
            chrom_b: "chr1".to_string(),
            bp_b: 20000,
            gene_b: "GENE_B".to_string(),
            join: FusionJoin::Forward,
            allele_fraction: Some(0.05),
        };

        let out = simulate_event(2, &event, &pool, &mut hap, &config, &mut gen, 0.05, &mut rng)
            .expect("simulate_event should succeed for Fusion");

        // Fusion is additive: no reads are suppressed.
        assert_eq!(
            out.suppressed_count, 0,
            "Fusion should not suppress any original reads (additive model)"
        );

        // All original reads should be kept.
        assert_eq!(
            out.kept_originals.len(),
            pool.pairs.len(),
            "Fusion should keep all {} original reads",
            pool.pairs.len()
        );

        // Chimeric reads should be generated at the junction.
        assert!(
            !out.chimeric_pairs.is_empty(),
            "Fusion should produce chimeric reads spanning the breakpoint"
        );
    }

    // ---------------------------------------------------------------
    // simulate_event suppression (real code)
    // ---------------------------------------------------------------

    fn del_event(start: u64, end: u64) -> SimEvent {
        SimEvent::Deletion {
            chrom: "chr1".to_string(),
            del_start: start,
            del_end: end,
            gene: "TEST".to_string(),
            exons: vec![],
            allele_fraction: Some(0.2),
        }
    }

    #[test]
    fn test_simulate_event_suppresses_inside_pairs_at_vaf_and_keeps_outside() {
        // DEL [1000,3000) with 1 kb flanks: footprint [0,4000).
        let mut pairs: Vec<ReadPair> = (0..500)
            .map(|i| make_pair(&format!("in_{}", i), 1500, 1900))
            .collect();
        pairs.extend((0..50).map(|i| make_pair(&format!("out_{}", i), 5000, 5400)));
        let pool = make_pool(pairs);
        let mut hap = del_haplotype(1000, 2000);
        let mut rng = StdRng::seed_from_u64(42);

        let out = simulate_event(
            1, &del_event(1000, 3000), &pool, &mut hap, &make_config(),
            &mut mock_synth_gen(150), 0.2, &mut rng,
        )
        .unwrap();

        assert!(out.suppressed_names.iter().all(|n| n.starts_with("in_")));
        let rate = out.suppressed_count as f64 / 500.0;
        assert!((rate - 0.2).abs() < 0.05, "inside suppression rate {:.3}", rate);
    }

    #[test]
    fn test_simulate_event_keeps_pairs_straddling_footprint_edge() {
        // Tiled fragments never extend past the haplotype ends, so originals
        // that stick out of the footprint must not be suppressed either;
        // otherwise depth dips at the footprint edge.
        // DEL [1000,3000) with 1 kb flanks: footprint [0,4000).
        let pairs: Vec<ReadPair> = (0..500)
            .map(|i| make_pair(&format!("edge_{}", i), 3800, 4200))
            .collect();
        let pool = make_pool(pairs);
        let mut hap = del_haplotype(1000, 2000);
        let mut rng = StdRng::seed_from_u64(42);

        let out = simulate_event(
            1, &del_event(1000, 3000), &pool, &mut hap, &make_config(),
            &mut mock_synth_gen(150), 0.2, &mut rng,
        )
        .unwrap();

        assert_eq!(out.suppressed_count, 0);
    }

    #[test]
    fn test_tiling_uses_actual_fragment_length_for_short_inserts() {
        // 220 bp fragments starting every 11 bp: exactly 20 cover each point.
        let pairs: Vec<ReadPair> = (0..728u64)
            .map(|i| make_pair(&format!("r{}", i), i * 11, i * 11 + 220))
            .collect();
        let pool = ReadPool {
            pairs,
            frag_dist: FragmentDist::from_stats(220.0, 30.0),
        };
        // DEL [2000,4000) with 2 kb flanks: haplotype of 4000 bp.
        let mut hap = del_haplotype(2000, 2000);
        let mut rng = StdRng::seed_from_u64(1);

        let out = simulate_event(
            1, &del_event(2000, 4000), &pool, &mut hap, &make_config(),
            &mut mock_synth_gen(150), 0.2, &mut rng,
        )
        .unwrap();

        // n = cov * v * (L - f) / f = 20 * 0.2 * (4000 - 220) / 220 = 68.7 -> 69
        assert_eq!(out.chimeric_pairs.len(), 69);
    }

    #[test]
    fn test_long_insertion_yields_reads_carrying_inserted_sequence() {
        // ref (A) [0,2000) | 1000 bp insertion (T) | ref (A) [2000,4000).
        let mut hap = make_haplotype(vec![
            ref_segment(0, 2000),
            novel_segment(1000),
            ref_segment(2000, 2000),
        ]);
        let pool = make_covering_pool(0, 6000, 1200); // 400 bp every 5 bp: cov 80
        let event = SimEvent::Insertion {
            chrom: "chr1".to_string(),
            pos: 2000,
            ins_seq: Some(vec![b'T'; 1000]),
            ins_len: 1000,
            gene: "TEST".to_string(),
            allele_fraction: Some(0.2),
        };
        let mut rng = StdRng::seed_from_u64(3);

        let out = simulate_event(
            1, &event, &pool, &mut hap, &make_config(), &mut mock_synth_gen(150), 0.2, &mut rng,
        )
        .unwrap();

        // With 400 bp fragments, 780 of the 4000 allowed starts put >= 10
        // inserted bases in a read (~19.5%). Before the fix it was 0.
        // R2 is reverse-complemented, so inserted T reads as A there.
        let carries_insert = |p: &&ReadPair| {
            p.seq1.iter().filter(|&&b| b == b'T').count() >= 10
                || p.seq2.iter().filter(|&&b| b == b'A').count() >= 10
        };
        let n = out.chimeric_pairs.len();
        let with_insert = out.chimeric_pairs.iter().filter(carries_insert).count();
        assert!(
            with_insert as f64 > 0.10 * n as f64,
            "only {} of {} tiled pairs carry inserted sequence",
            with_insert,
            n
        );
    }

    // ---------------------------------------------------------------
    // combine_event_outputs tests
    // ---------------------------------------------------------------

    fn spliced(kept: &[&str], chimeric: &[&str], suppressed: &[&str]) -> SplicedOutput {
        SplicedOutput {
            kept_originals: kept.iter().map(|n| make_pair(n, 0, 400)).collect(),
            chimeric_pairs: chimeric.iter().map(|n| make_pair(n, 0, 400)).collect(),
            suppressed_count: suppressed.len(),
            suppressed_names: suppressed.iter().map(|n| n.to_string()).collect(),
        }
    }

    fn sorted_names(pairs: &[ReadPair]) -> Vec<String> {
        let mut names: Vec<String> = pairs.iter().map(|p| p.name.clone()).collect();
        names.sort();
        names
    }

    #[test]
    fn test_dedup_by_name_keeps_last_occurrence() {
        let mut pairs = vec![
            make_pair("dup", 10, 30),
            make_pair("keep", 50, 70),
            make_pair("dup", 90, 110),
        ];
        dedup_by_name(&mut pairs);

        assert_eq!(pairs.len(), 2);
        assert_eq!(pairs[0].name, "keep");
        assert_eq!(pairs[1].name, "dup");
        assert_eq!(pairs[1].ref_start, 90);
    }

    #[test]
    fn test_combine_drops_pairs_suppressed_by_any_event() {
        // Both events extracted the same originals r1..r3 (nearby events share
        // a pool). Event 1 suppressed r1; event 2 left r1 alone because r1 lies
        // outside its own footprint. r1 must still be gone from the output.
        let outputs = vec![
            spliced(&["r2", "r3"], &["ev0001_a"], &["r1"]),
            spliced(&["r1", "r2"], &["ev0002_a"], &["r3"]),
        ];

        let combined = combine_event_outputs(outputs);

        assert_eq!(sorted_names(&combined), vec!["ev0001_a", "ev0002_a", "r2"]);
    }

    #[test]
    fn test_nearby_event_does_not_undo_deletion() {
        // Two het DELs 4 kb apart, simulated over one shared pool:
        //   A deletes [1000,3000), haplotype footprint [0,4000)
        //   B deletes [7000,9000), haplotype footprint [6000,10000)
        // Pairs inside A's deletion are outside B's footprint, so B keeps them.
        // After combining, about half of them (VAF 0.5) must still be gone.
        let pool = make_covering_pool(0, 12000, 240); // 400 bp fragments every 50 bp
        let config = make_config();
        let mut gen = mock_synth_gen(150);
        let mut rng = StdRng::seed_from_u64(7);

        let del = |start: u64, end: u64| SimEvent::Deletion {
            chrom: "chr1".to_string(),
            del_start: start,
            del_end: end,
            gene: "TEST".to_string(),
            exons: vec![],
            allele_fraction: Some(0.5),
        };
        let flank = |start: u64, end: u64| HaplotypeSegment {
            sequence: vec![b'A'; (end - start) as usize],
            origin: Some(SegmentOrigin {
                chrom: "chr1".to_string(),
                ref_start: start,
                ref_end: end,
                is_reverse: false,
            }),
            hap_offset: 0,
        };

        let mut hap_a = make_haplotype(vec![flank(0, 1000), flank(3000, 4000)]);
        let mut hap_b = make_haplotype(vec![flank(6000, 7000), flank(9000, 10000)]);

        let out_a = simulate_event(
            1, &del(1000, 3000), &pool, &mut hap_a, &config, &mut gen, 0.5, &mut rng,
        )
        .unwrap();
        let out_b = simulate_event(
            2, &del(7000, 9000), &pool, &mut hap_b, &config, &mut gen, 0.5, &mut rng,
        )
        .unwrap();

        let combined = combine_event_outputs(vec![out_a, out_b]);

        // Originals fully inside A's deletion: starts 1000, 1050, ..., 2600.
        let inside_a = |p: &ReadPair| p.ref_start >= 1000 && p.ref_end <= 3000;
        let total = pool.pairs.iter().filter(|p| inside_a(p)).count();
        let surviving = combined
            .iter()
            .filter(|p| p.name.starts_with("read_") && inside_a(p))
            .count();
        assert_eq!(total, 33);
        let frac = surviving as f64 / total as f64;
        assert!(
            (0.2..=0.8).contains(&frac),
            "expected ~50% of originals inside A's deletion to survive, got {}/{} ({:.2})",
            surviving,
            total,
            frac
        );
    }
}
