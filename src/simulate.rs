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
    event: &SimEvent,
    pool: &ReadPool,
    haplotype: &mut VariantHaplotype,
    config: &SimConfig,
    synth_gen: &mut SynthReadGenerator,
    vaf: f64,
    rng: &mut StdRng,
) -> Result<SplicedOutput> {
    // Get SV boundaries for read classification.
    let (sv_start, sv_end) = match event.primary_region() {
        Some((_chrom, start, end)) => (start, end),
        None => {
            // Fusion: use the first breakpoint as a point event.
            if let SimEvent::Fusion { bp_a, .. } = event {
                (*bp_a, *bp_a)
            } else {
                unreachable!()
            }
        }
    };

    // For DUPs and fusions, we keep all original reads and add chimeric on top.
    let is_additive = matches!(event, SimEvent::Duplication { .. } | SimEvent::Fusion { .. });

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
    let (hap_ref_start, hap_ref_end) = haplotype
        .ref_range()
        .unwrap_or((sv_start, sv_end));

    let mut kept = Vec::new();
    let mut suppressed = 0usize;

    for pair in &pool.pairs {
        if is_additive {
            kept.push(pair.clone());
            continue;
        }

        // Check if read overlaps the haplotype reference footprint.
        let in_haplotype = classify_pair_relation(pair, hap_ref_start, hap_ref_end);
        if matches!(in_haplotype, PairRelation::Outside) {
            // Entirely outside the haplotype range — keep, unaffected by the SV.
            kept.push(pair.clone());
            continue;
        }

        // Read is within (or overlapping) the haplotype range.
        // For inside/overlapping the SV with LOH info, use haplotype-aware
        // suppression for classified reads, random fallback for unclassified.
        // For haplotype-flank reads (or without LOH), use random.
        let in_sv = classify_pair_relation(pair, sv_start, sv_end);
        let use_loh_for_this = use_loh
            && matches!(
                in_sv,
                PairRelation::Inside | PairRelation::Overlapping
            );

        // For reads that only partially overlap the haplotype range, scale
        // suppression probability by the overlap fraction. This prevents
        // boundary artifacts where reads straddling the edge get fully
        // suppressed despite most of their coverage being outside.
        let overlap_frac = if matches!(in_haplotype, PairRelation::Overlapping) {
            let overlap_start = pair.ref_start.max(hap_ref_start);
            let overlap_end = pair.ref_end.min(hap_ref_end);
            let frag_len = pair.ref_end.saturating_sub(pair.ref_start).max(1) as f64;
            (overlap_end.saturating_sub(overlap_start)) as f64 / frag_len
        } else {
            1.0 // fully inside → full suppression rate
        };

        let effective_vaf = vaf * overlap_frac;

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
                    let hap_prob = (2.0 * vaf).min(1.0) * overlap_frac;
                    if rng.gen::<f64>() < hap_prob {
                        suppressed += 1;
                    } else {
                        kept.push(pair.clone());
                    }
                } else {
                    // Non-variant haplotype: only suppress when VAF > 0.5.
                    let other_prob = (2.0 * vaf - 1.0).max(0.0) * overlap_frac;
                    if rng.gen::<f64>() < other_prob {
                        suppressed += 1;
                    } else {
                        kept.push(pair.clone());
                    }
                }
            } else {
                // Read couldn't be classified (no het SNP overlap) —
                // fall back to random suppression at effective VAF rate.
                if rng.gen::<f64>() > effective_vaf {
                    kept.push(pair.clone());
                } else {
                    suppressed += 1;
                }
            }
        } else if rng.gen::<f64>() > effective_vaf {
            kept.push(pair.clone());
        } else {
            suppressed += 1;
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
        rng,
    );

    // For DUPs, also generate depth copies inside the region.
    let depth_copies = if let SimEvent::Duplication {
        dup_start, dup_end, ..
    } = event
    {
        synth_gen.generate_dup_depth_copies(
            pool,
            *dup_start,
            *dup_end,
            vaf,
            &loh_set,
            &classified_set,
            rng,
        )
    } else {
        Vec::new()
    };

    log::info!(
        "simulate_event: {} kept, {} suppressed, {} chimeric (haplotype-tiled), {} depth copies",
        kept.len(),
        suppressed,
        chimeric.len(),
        depth_copies.len(),
    );

    let mut all_chimeric = chimeric;
    all_chimeric.extend(depth_copies);

    Ok(SplicedOutput {
        chimeric_pairs: all_chimeric,
        kept_originals: kept,
        suppressed_count: suppressed,
    })
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
        } => {
            loh::identify_duplicated_haplotype_reads(
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
            })
        }
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
/// For additive events (DUP, Fusion): uses the aggregate breakpoint zone
/// length (2 × mean_frag per breakpoint) since reads are placed only near
/// segment boundaries.
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

    let effective_len = if breakpoint_only && !breakpoints.is_empty() {
        // Aggregate breakpoint zone length: 2 × mean_frag per breakpoint.
        let zone_per_bp = (2.0 * mean_frag) as u64;
        let total_zone = breakpoints.len() as u64 * zone_per_bp;
        total_zone.min(haplotype.total_len) as f64
    } else {
        // Use reference-mapped length (excludes novel insertion sequence).
        let ref_len = haplotype.ref_mapped_len();
        if ref_len > 0 { ref_len as f64 } else { haplotype.total_len as f64 }
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
    rng: &mut StdRng,
) -> Vec<ReadPair> {
    let hap_len = haplotype.total_len;
    if hap_len == 0 {
        return Vec::new();
    }

    let mean_frag = pool.frag_dist.mean.max(300.0);
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
    let has_novel = haplotype
        .segments
        .iter()
        .any(|seg| seg.origin.is_none());

    for i in 0..n_frags {
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

        let name = format!("sim_hap_{:06}", i);

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
/// Uses 50 evenly-spaced sample points in a 2000bp window for stable estimates.
fn estimate_coverage_at(pool: &ReadPool, pos: u64, window: u64) -> f64 {
    let start = pos.saturating_sub(window / 2);
    let end = pos + window / 2;

    let n_samples = 50;
    let step = (end - start) / n_samples as u64;
    if step == 0 {
        return 0.0;
    }
    let mut total = 0usize;

    for i in 0..n_samples {
        let sample_pos = start + i as u64 * step;
        let count = pool
            .pairs
            .iter()
            .filter(|p| p.ref_start <= sample_pos && p.ref_end > sample_pos)
            .count();
        total += count;
    }

    total as f64 / n_samples as f64
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::haplotype::{HaplotypeSegment, SegmentOrigin};
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
            ref_segment(0, 2000),     // left flank
            ref_segment(7000, 2000),  // right flank (after 5kb deleted region)
        ]);
        assert_eq!(hap.total_len, 4000);
        assert_eq!(hap.ref_mapped_len(), 4000);

        // 30x coverage, 0.5 VAF, 400bp mean frag.
        let count = compute_tiling_count(&hap, 30.0, 0.5, 400.0, false);
        // Expected: 30 * 0.5 * 4000 / 400 = 150
        assert_eq!(count, 150);
    }

    #[test]
    fn test_tiling_count_ins() {
        // 500bp INS with 2kb flanks: ref_mapped_len = 4kb (flanks), total = 4.5kb.
        // Should use ref_mapped_len (4kb), NOT total_len (4.5kb).
        let hap = make_haplotype(vec![
            ref_segment(0, 2000),     // left flank
            novel_segment(500),       // inserted sequence
            ref_segment(2000, 2000),  // right flank
        ]);
        assert_eq!(hap.total_len, 4500);
        assert_eq!(hap.ref_mapped_len(), 4000);

        let count = compute_tiling_count(&hap, 30.0, 0.5, 400.0, false);
        // Expected: 30 * 0.5 * 4000 / 400 = 150 (not 168.75 if using total_len).
        assert_eq!(count, 150);
    }

    #[test]
    fn test_tiling_count_inv() {
        // 5kb INV with 2kb flanks: total = 9kb, ref_mapped_len = 9kb.
        // Same ref footprint as before the INV, just middle is reversed.
        let hap = make_haplotype(vec![
            ref_segment(0, 2000),     // left flank
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
            ref_segment(7000, 2000),  // right flank
        ]);
        assert_eq!(hap.total_len, 9000);
        assert_eq!(hap.ref_mapped_len(), 9000);

        let count = compute_tiling_count(&hap, 30.0, 0.5, 400.0, false);
        // Expected: 30 * 0.5 * 9000 / 400 = 337.5 → 338
        assert_eq!(count, 338);
    }

    #[test]
    fn test_tiling_count_dup_breakpoint_only() {
        // DUP junction haplotype: breakpoint_only = true.
        // Haplotype: [pre-dup flank] [post-dup flank] with 1 breakpoint.
        let hap = make_haplotype(vec![
            ref_segment(0, 2000),     // up to dup end
            ref_segment(0, 2000),     // from dup start (junction)
        ]);

        // With mean_frag=400, zone_per_bp = 800bp per breakpoint.
        // 1 breakpoint → effective_len = 800.
        let count = compute_tiling_count(&hap, 30.0, 0.5, 400.0, true);
        // Expected: 30 * 0.5 * 800 / 400 = 30
        assert_eq!(count, 30);
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
        assert!(cov < 20.0 && cov > 0.0, "expected partial coverage, got {:.1}", cov);
    }

    // ---------------------------------------------------------------
    // Suppression balance tests
    //
    // These test the suppression logic directly using the internal
    // classify_pair_relation + random suppression pattern, without
    // needing a SynthReadGenerator (which requires a real FASTA).
    // ---------------------------------------------------------------

    /// Simulate the suppression loop from simulate_event for non-additive events.
    /// Returns (kept_count, suppressed_count).
    fn run_suppression(
        pairs: &[ReadPair],
        hap_ref_start: u64,
        hap_ref_end: u64,
        _sv_start: u64,
        _sv_end: u64,
        is_additive: bool,
        vaf: f64,
        rng: &mut StdRng,
    ) -> (usize, usize) {
        let mut kept = 0usize;
        let mut suppressed = 0usize;

        for pair in pairs {
            if is_additive {
                kept += 1;
                continue;
            }

            let in_haplotype = classify_pair_relation(pair, hap_ref_start, hap_ref_end);
            if matches!(in_haplotype, PairRelation::Outside) {
                kept += 1;
                continue;
            }

            // Scale suppression by overlap fraction for boundary-spanning reads.
            let overlap_frac = if matches!(in_haplotype, PairRelation::Overlapping) {
                let overlap_start = pair.ref_start.max(hap_ref_start);
                let overlap_end = pair.ref_end.min(hap_ref_end);
                let frag_len = pair.ref_end.saturating_sub(pair.ref_start).max(1) as f64;
                (overlap_end.saturating_sub(overlap_start)) as f64 / frag_len
            } else {
                1.0
            };

            // No LOH in these tests — always random suppression.
            if rng.gen::<f64>() > vaf * overlap_frac {
                kept += 1;
            } else {
                suppressed += 1;
            }
        }

        (kept, suppressed)
    }

    #[test]
    fn test_suppression_keeps_outside_reads() {
        use rand::SeedableRng;
        let mut rng = StdRng::seed_from_u64(42);

        // Reads at [100, 200) are outside haplotype range [200, 500).
        let outside_pairs: Vec<ReadPair> = (0..50)
            .map(|i| make_pair(&format!("outside_{}", i), 100, 200))
            .collect();

        let (kept, suppressed) = run_suppression(
            &outside_pairs,
            200, 500,   // hap_ref range
            300, 400,   // sv range
            false,      // not additive
            0.5,
            &mut rng,
        );

        assert_eq!(kept, 50);
        assert_eq!(suppressed, 0);
    }

    #[test]
    fn test_suppression_rate_matches_vaf() {
        use rand::SeedableRng;
        let mut rng = StdRng::seed_from_u64(42);

        // 500 reads inside the haplotype range — large N for statistical reliability.
        let inside_pairs: Vec<ReadPair> = (0..500)
            .map(|i| make_pair(&format!("in_{}", i), 200, 800))
            .collect();

        let (kept, suppressed) = run_suppression(
            &inside_pairs,
            0, 1000,    // hap_ref range
            100, 900,   // sv range
            false,
            0.5,
            &mut rng,
        );

        // Kept + suppressed = total.
        assert_eq!(kept + suppressed, 500);

        // Suppression rate ≈ 0.5 (within 10% for 500 reads).
        let rate = suppressed as f64 / 500.0;
        assert!(
            (rate - 0.5).abs() < 0.10,
            "suppression rate {:.3} too far from 0.5",
            rate
        );
    }

    #[test]
    fn test_suppression_additive_keeps_all() {
        use rand::SeedableRng;
        let mut rng = StdRng::seed_from_u64(42);

        // DUP: additive events keep all reads.
        let pairs: Vec<ReadPair> = (0..100)
            .map(|i| make_pair(&format!("r_{}", i), 400, 600))
            .collect();

        let (kept, suppressed) = run_suppression(
            &pairs,
            200, 800,   // hap_ref range
            300, 700,   // sv range
            true,       // additive (DUP)
            0.5,
            &mut rng,
        );

        assert_eq!(kept, 100);
        assert_eq!(suppressed, 0);
    }

    #[test]
    fn test_suppression_mixed_inside_outside() {
        use rand::SeedableRng;
        let mut rng = StdRng::seed_from_u64(123);

        // Mix of outside and inside reads.
        let mut pairs = Vec::new();
        for i in 0..50 {
            pairs.push(make_pair(&format!("out_{}", i), 10, 90)); // outside [100, 900)
        }
        for i in 0..200 {
            pairs.push(make_pair(&format!("in_{}", i), 300, 700)); // inside
        }

        let (kept, suppressed) = run_suppression(
            &pairs,
            100, 900,   // hap_ref range
            100, 900,   // sv range
            false,
            0.5,
            &mut rng,
        );

        // All 50 outside reads kept + roughly half of 200 inside reads kept.
        assert_eq!(kept + suppressed, 250);
        assert!(kept >= 50, "should keep at least the 50 outside reads, kept {}", kept);
        assert!(suppressed > 0, "should suppress some inside reads");

        // Outside reads are always kept, so suppressed must come from the 200 inside reads.
        assert!(suppressed <= 200);
        let inside_suppressed_rate = suppressed as f64 / 200.0;
        assert!(
            (inside_suppressed_rate - 0.5).abs() < 0.12,
            "inside suppression rate {:.3} too far from 0.5",
            inside_suppressed_rate
        );
    }

    #[test]
    fn test_suppression_overlap_fraction_reduces_rate() {
        use rand::SeedableRng;
        let mut rng = StdRng::seed_from_u64(42);

        // 500 reads that span the haplotype boundary: each read is [150, 550),
        // haplotype range is [200, 1000). Overlap = 350/400 = 87.5%.
        // Effective suppression rate = 0.5 * 0.875 = 0.4375.
        let overlapping_pairs: Vec<ReadPair> = (0..500)
            .map(|i| make_pair(&format!("ov_{}", i), 150, 550))
            .collect();

        let (kept, suppressed) = run_suppression(
            &overlapping_pairs,
            200, 1000,  // hap_ref range
            200, 1000,  // sv range
            false,
            0.5,
            &mut rng,
        );

        assert_eq!(kept + suppressed, 500);

        // Expected rate ≈ 0.4375 (less than the 0.5 VAF due to overlap fraction).
        let rate = suppressed as f64 / 500.0;
        assert!(
            rate < 0.5,
            "overlap fraction should reduce suppression rate below VAF, got {:.3}",
            rate
        );
        assert!(
            (rate - 0.4375).abs() < 0.08,
            "suppression rate {:.3} too far from expected 0.4375",
            rate
        );
    }

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
}
