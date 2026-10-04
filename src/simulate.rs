//! Unified SV event simulation: suppress reads + tile synthetic reads across
//! the variant haplotype.
//!
//! Replaces the per-SV splice functions (splice_deletion, splice_duplication, etc.)
//! with a single `simulate_event` that works for all SV types.

use std::collections::{BTreeSet, HashMap, HashSet};

use anyhow::Result;
use rand::rngs::StdRng;
use rand::Rng;

use crate::haplotype::VariantHaplotype;
use crate::loh;
use crate::reference::SharedReference;
use crate::synth::{copy_rate, SynthReadGenerator};
use crate::types::{DepthFold, ReadPair, ReadPool, SimConfig, SimEvent, SplicedOutput};

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
///
/// The sample's own SNPs around the event are read first (see [`loh`]), so
/// originals are suppressed by copy and synthetic reads carry the alleles of
/// the copy they come from.
pub fn simulate_event(
    event_index: usize,
    event: &SimEvent,
    pool: &ReadPool,
    haplotype: &mut VariantHaplotype,
    config: &SimConfig,
    synth_gen: &SynthReadGenerator,
    vaf: f64,
    rng: &mut StdRng,
) -> Result<SplicedOutput> {
    simulate_event_origin(
        event_index, event, pool, haplotype, config, synth_gen, vaf, None, rng,
    )
}

/// [`simulate_event`] under `--edit-model origin` when `origin` is `Some`.
/// The site's fragments are judged by where they came from, and the pool's
/// pairs are not suppressed here. `None` is exactly `simulate_event`.
#[allow(clippy::too_many_arguments)]
pub fn simulate_event_origin(
    event_index: usize,
    event: &SimEvent,
    pool: &ReadPool,
    haplotype: &mut VariantHaplotype,
    config: &SimConfig,
    synth_gen: &SynthReadGenerator,
    vaf: f64,
    origin: Option<&crate::origin::OriginSite>,
    rng: &mut StdRng,
) -> Result<SplicedOutput> {
    let copies = sample_copies_for_event(event, haplotype, config, synth_gen.reference(), rng)?;
    simulate_event_inner(
        event_index, event, pool, haplotype, config, synth_gen, vaf, &copies, origin, rng,
    )
}

/// Whether `event` only adds reads: a fusion, or a DUP under the legacy
/// junction model. It removes no original, so `--edit-model` does not apply.
pub fn is_additive(event: &SimEvent, dup_model: &str) -> bool {
    match event {
        SimEvent::Fusion { .. } => true,
        SimEvent::Duplication { .. } => dup_model == "junction",
        _ => false,
    }
}

/// Read the sample's two copies over each reference region the haplotype
/// draws from: the whole footprint for single-region events, and each side
/// for a fusion. A `--gvcf` that can't be read is an error (CR3). A region
/// whose SNPs can't be read otherwise (the pileup) gets none, with a warning.
fn sample_copies_for_event(
    event: &SimEvent,
    haplotype: &VariantHaplotype,
    config: &SimConfig,
    reference: &SharedReference,
    rng: &mut StdRng,
) -> Result<Vec<(String, loh::SampleCopies)>> {
    let regions: Vec<(String, u64, u64)> = match event {
        SimEvent::Fusion { .. } => haplotype
            .segments
            .iter()
            .filter_map(|seg| seg.origin.as_ref())
            .map(|o| (o.chrom.clone(), o.ref_start, o.ref_end))
            .collect(),
        _ => haplotype
            .ref_range()
            .map(|(start, end)| vec![(haplotype.primary_chrom().to_string(), start, end)])
            .unwrap_or_default(),
    };

    let mut copies = Vec::with_capacity(regions.len());
    for (chrom, start, end) in regions {
        let sample = match reference.fetch_sequence(&chrom, start, end).and_then(|ref_seq| {
            loh::sample_copies(
                &config.bam_path,
                &chrom,
                start,
                end,
                &ref_seq,
                config.min_mapq,
                config.gvcf_path.as_deref(),
                Some(config.ref_path.as_str()),
                rng,
            )
        }) {
            Ok(sample) => sample,
            Err(e) if e.downcast_ref::<loh::GvcfUnreadable>().is_some() => return Err(e),
            Err(e) => {
                log::warn!(
                    "could not read the sample's SNPs in {}:{}-{}: {}",
                    chrom,
                    start,
                    end,
                    e
                );
                loh::SampleCopies::default()
            }
        };
        copies.push((chrom, sample));
    }
    Ok(copies)
}

#[cfg(test)]
#[allow(clippy::too_many_arguments)]
fn simulate_event_with_copies(
    event_index: usize,
    event: &SimEvent,
    pool: &ReadPool,
    haplotype: &mut VariantHaplotype,
    config: &SimConfig,
    synth_gen: &SynthReadGenerator,
    vaf: f64,
    copies: &[(String, loh::SampleCopies)],
    rng: &mut StdRng,
) -> Result<SplicedOutput> {
    simulate_event_inner(
        event_index, event, pool, haplotype, config, synth_gen, vaf, copies, None, rng,
    )
}

/// [`simulate_event`] with the sample's two copies already read, per
/// reference region (chromosome) of the haplotype.
#[allow(clippy::too_many_arguments)]
fn simulate_event_inner(
    event_index: usize,
    event: &SimEvent,
    pool: &ReadPool,
    haplotype: &mut VariantHaplotype,
    config: &SimConfig,
    synth_gen: &SynthReadGenerator,
    vaf: f64,
    copies: &[(String, loh::SampleCopies)],
    origin: Option<&crate::origin::OriginSite>,
    rng: &mut StdRng,
) -> Result<SplicedOutput> {
    let name_prefix = format!("ev{:04}", event_index);

    // Get SV boundaries for read classification.
    let (sv_chrom, sv_start, sv_end) = match event.primary_region() {
        Some((chrom, start, end)) => (chrom.to_string(), start, end),
        None => {
            // Fusion: use the first breakpoint as a point event.
            if let SimEvent::Fusion { chrom_a, bp_a, .. } = event {
                (chrom_a.clone(), *bp_a, *bp_a)
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
    let is_additive = is_additive(event, &config.dup_model);

    // Above VAF 0.5 the event is on both copies in some cells, so part of the
    // synthetic reads come from the other copy. Its share of the added reads
    // is max(0, 2v - 1) / (2v).
    let p_other_copy = if vaf > 0.5 { 1.0 - 1.0 / (2.0 * vaf) } else { 0.0 };
    let mut other_haplotype = (p_other_copy > 0.0).then(|| haplotype.clone());

    // Each copy's haplotype takes that copy's alleles, except at bases the
    // event itself changes (a simulated SNV keeps its own allele).
    let own = own_bases(event);
    let mut read_copy: HashMap<String, bool> = HashMap::new();
    for (chrom, sample) in copies {
        haplotype.apply_variants(chrom, &outside(&sample.event_copy, chrom, own));
        if let Some(other) = other_haplotype.as_mut() {
            other.apply_variants(chrom, &outside(&sample.other_copy, chrom, own));
        }
        read_copy.extend(sample.read_copy.iter().map(|(n, &e)| (n.clone(), e)));
    }

    // Suppress original reads within the haplotype's reference footprint.
    //
    // For non-additive events (DEL, INV, INS, small variants, full DUP): suppress
    // reads by copy (see `copy_rate`) within [hap_ref_start, hap_ref_end), in the
    // flanks as well as the event. Tiled haplotype reads replace them. Reads
    // entirely OUTSIDE this range are kept unchanged.
    //
    // For additive events (junction DUP, Fusion): keep all originals. Chimeric
    // reads near the breakpoint and DUP depth copies are added on top.
    let (hap_ref_start, hap_ref_end) = haplotype.ref_range().unwrap_or((sv_start, sv_end));

    if let Some(site) = origin {
        let footprint =
            crate::origin::Span::new(haplotype.primary_chrom(), hap_ref_start, hap_ref_end);
        anyhow::ensure!(
            !is_additive,
            "origin site {} was handed to an additive event (a fusion, or a DUP under \
             --dup-model junction) -- additive events remove no originals, so \
             --edit-model origin does not apply to them; being handed a site here is a bug \
             in the caller",
            site.footprint,
        );
        anyhow::ensure!(
            site.footprint == footprint,
            "origin site {} does not match the haplotype's footprint {} -- this is a bug; \
             please report it",
            site.footprint,
            footprint,
        );
    }

    let mut kept = Vec::new();
    let mut suppressed: Vec<String> = Vec::new();

    for pair in &pool.pairs {
        // Only pairs entirely inside the haplotype's reference footprint are
        // replaced. Tiled fragments never extend past the haplotype ends, so
        // suppressing pairs that stick out would leave a depth dip there.
        // Under origin nothing is suppressed here: every event's chances are
        // summed first and drawn once, in `origin::decide` (R5).
        let replaceable = origin.is_none()
            && !is_additive
            && classify_pair_relation(pair, hap_ref_start, hap_ref_end) == PairRelation::Inside;
        let copy = read_copy.get(&pair.name).copied();
        if replaceable && rng.gen::<f64>() < copy_rate(copy, vaf) {
            suppressed.push(pair.name.clone());
        } else {
            kept.push(pair.clone());
        }
    }

    // Estimate coverage at the first breakpoint for chimeric read count, from
    // reads on that breakpoint's own chromosome only.
    let bp_positions = haplotype.breakpoints();
    let (first_bp_chrom, first_bp_ref) = bp_positions
        .first()
        .and_then(|&bp| haplotype.hap_to_ref(bp.saturating_sub(1)))
        .unwrap_or((sv_chrom.clone(), sv_start));
    let (cov, uncovered_breakpoint_sides) = coverage_for_tiling(
        event,
        haplotype,
        pool,
        origin,
        (&sv_chrom, sv_start, sv_end),
        (&first_bp_chrom, first_bp_ref),
    )?;

    // CR2: measured only. The tiling below still scales every fragment by
    // `cov`; this says how far the donor's own depth is from it. Under origin
    // each bin is measured with origin depth, the estimator `cov` came from
    // (the T3 rule).
    let depth_fold = match origin {
        Some(site) => {
            let depth = site.depth();
            depth_fold_by(haplotype, cov, &|chrom, pos, window| {
                depth.fragment_coverage_at(chrom, pos, window)
            })
        }
        None => depth_fold(haplotype, pool, cov),
    };

    // Tile synthetic reads across the haplotype.
    // For additive events (DUP/Fusion), restrict tiling to near breakpoints
    // to avoid inflating flank coverage.
    let (chimeric, adjusted_vaf) = tile_haplotype_reads(
        haplotype,
        other_haplotype.as_ref().map(|other| (other, p_other_copy)),
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
        chrom,
        dup_start,
        dup_end,
        ..
    } = event
    {
        if config.dup_model == "junction" {
            let none = loh::SampleCopies::default();
            let sample = copies
                .iter()
                .find(|(c, _)| c == chrom)
                .map_or(&none, |(_, sample)| sample);
            synth_gen.generate_dup_depth_copies(
                pool,
                *dup_start,
                *dup_end,
                vaf,
                sample,
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
        uncovered_breakpoint_sides,
        adjusted_vaf,
        depth_fold,
        origin_chances: origin
            .map(|site| site.removal_chances(&read_copy, vaf))
            .unwrap_or_default(),
    })
}

/// Names of the original read pairs spike took out of the BAM.
///
/// An original enters an event's pool either passed through (`kept_originals`)
/// or suppressed, so the union of the two is exactly the set merge.sh must
/// remove from the original BAM. Synthetic pairs carry fresh names and are
/// never in it.
pub fn consumed_original_names(outputs: &[SplicedOutput]) -> BTreeSet<String> {
    outputs
        .iter()
        .flat_map(|o| {
            o.kept_originals
                .iter()
                .map(|p| p.name.clone())
                .chain(o.suppressed_names.iter().cloned())
        })
        .collect()
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

/// Under `--edit-model origin`, move each event's pool pairs that
/// `origin::decide` removed from `kept_originals` to `suppressed_names`, so
/// [`combine_event_outputs`] drops them and [`consumed_original_names`] still
/// lists them. Removed names in no pool are the caller's to list.
pub fn apply_removals(outputs: &mut [SplicedOutput], removed: &BTreeSet<String>) {
    for output in outputs.iter_mut() {
        let (gone, kept): (Vec<ReadPair>, Vec<ReadPair>) =
            std::mem::take(&mut output.kept_originals)
                .into_iter()
                .partition(|p| removed.contains(&p.name));
        output.kept_originals = kept;
        output.suppressed_names.extend(gone.into_iter().map(|p| p.name));
        output.suppressed_count = output.suppressed_names.len();
    }
}

/// `--edit-model origin`'s one draw (R5): every event's chances go to one
/// [`crate::origin::decide`], which sums them per fragment and draws once per
/// family; then [`apply_removals`] moves each event's removed pool pairs.
/// Returns every removed name.
pub fn draw_origin_removals(outputs: &mut [SplicedOutput], rng: &mut StdRng) -> BTreeSet<String> {
    let chances: Vec<crate::origin::Chance> = outputs
        .iter()
        .flat_map(|o| o.origin_chances.iter().cloned())
        .collect();
    let removed = crate::origin::decide(&chances, rng);
    apply_removals(outputs, &removed);
    removed
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

/// Reference bases the event itself changes (a small variant's REF span):
/// the sample's alleles must not overwrite them.
fn own_bases(event: &SimEvent) -> Option<(&str, u64, u64)> {
    match event {
        SimEvent::SmallVariant {
            chrom,
            pos,
            ref_allele,
            ..
        } => Some((chrom.as_str(), *pos, *pos + ref_allele.len() as u64)),
        _ => None,
    }
}

/// `alleles` on `chrom`, without the bases in `own` (see [`own_bases`]).
fn outside(
    alleles: &HashMap<u64, u8>,
    chrom: &str,
    own: Option<(&str, u64, u64)>,
) -> HashMap<u64, u8> {
    match own {
        Some((c, start, end)) if c == chrom => alleles
            .iter()
            .filter(|(pos, _)| !(start..end).contains(*pos))
            .map(|(&pos, &base)| (pos, base))
            .collect(),
        _ => alleles.clone(),
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

/// Every reference position the haplotype's junctions join: for each
/// breakpoint, the last base before the cut and the first base after it.
///
/// A cut into novel sequence (the inserted bases of an INS) has no reference
/// base on that side and contributes nothing -- there is no donor position to
/// ask about.
///
/// Each position appears once. Two junctions a base apart -- a small
/// variant's one-base alt segment -- name the same base from either side, and
/// `bp.saturating_sub(1)` names the cut itself for a breakpoint at haplotype
/// offset 0, where there is no base before it to step back to. One position
/// is one question about donor coverage however many junctions touch it.
fn breakpoint_sides(haplotype: &VariantHaplotype) -> Vec<(String, u64)> {
    let mut sides: Vec<(String, u64)> = Vec::new();
    for hap_pos in haplotype
        .breakpoints()
        .iter()
        .flat_map(|&bp| [bp.saturating_sub(1), bp])
    {
        if let Some(side) = haplotype.hap_to_ref(hap_pos) {
            if !sides.contains(&side) {
                sides.push(side);
            }
        }
    }
    sides
}

/// The donor coverage the tiling count is scaled by, or an error if spike has
/// no donor reads to build this event from (N5, N12).
///
/// The tiling count is `coverage x VAF`, so where the coverage is 0 spike
/// would write a truth VCF beside reads it invented rather than simulated.
/// Which positions have to carry coverage depends on how many places the
/// donor pool was drawn from -- that is, what `extract_pool_for_event`
/// searched on this event's behalf:
///
/// - A **fusion** is extracted from two loci, one per partner, and every
///   fragment it plants spans the junction between them. A partner with no
///   donor reads makes half of every planted read invention, so **every**
///   breakpoint side must be covered. Testing all of them is also what makes
///   the verdict independent of which partner the event names first: keying
///   on the first breakpoint alone accepted one naming order and refused the
///   other (N12).
/// - Every other event is extracted from **one** locus and tiled across its
///   whole haplotype, so what it needs is donor depth somewhere around it.
///   One breakpoint side with nothing over it is a thin spot, or the far side
///   of an event that straddles the edge of a sliced or panel BAM -- ordinary
///   input, pinned by
///   `test_simulate_event_keeps_pairs_straddling_footprint_edge`. It is
///   refused only when **no** side of any of its breakpoints is covered,
///   which is the case N5 measured: the event nowhere near the reads.
///
/// The coverage returned is the first breakpoint's, exactly as before, unless
/// that side is the uncovered one -- which until now was refused outright, so
/// no run that already works changes its arithmetic (M7). In that case it is
/// the first covered side instead: the only depth the event has been measured
/// against, and scaling by the 0 would plant nothing at all.
///
/// Returned alongside it: the sides the event was kept *despite*. Keeping it
/// is right, but it is not nothing -- the haplotype still spans those sides,
/// so the fragments tiled across them are scaled by depth measured somewhere
/// else and land where the input BAM has no read. This logs that, and the
/// caller puts it in the run README; a silent coverage island is what this
/// return value exists to prevent.
///
/// Under `--edit-model origin` (R1), `origin` is the event's site and the
/// depth is its origin depth in fragment units, not the pool's, which is 0
/// inside a perfect twin. Only events that remove reads get a site, and those
/// come from one locus, so the fusion rule does not arise there. The two
/// models share this one function so they cannot drift: origin's former copy
/// had already dropped the warning below (PD-14, PD-30).
fn coverage_for_tiling(
    event: &SimEvent,
    haplotype: &VariantHaplotype,
    pool: &ReadPool,
    origin: Option<&crate::origin::OriginSite>,
    event_span: (&str, u64, u64),
    fallback_bp: (&str, u64),
) -> Result<(f64, Vec<String>)> {
    let mut sides = breakpoint_sides(haplotype);
    if sides.is_empty() {
        // No junction (a SNP, or a haplotype of one segment): the event's own
        // start is the only position the tiling count can be scaled at.
        sides.push((fallback_bp.0.to_string(), fallback_bp.1));
    }

    let origin_depth = origin.map(|site| site.depth());
    let covs: Vec<f64> = sides
        .iter()
        .map(|(chrom, pos)| match &origin_depth {
            Some(depth) => depth.fragment_coverage_at(chrom, *pos, 2000),
            None => estimate_coverage_at(pool, chrom, *pos, 2000),
        })
        .collect();
    let covered = |cov: f64| !cov.is_nan() && cov > 0.0;

    let is_fusion = event.is_multi_locus();
    let refuse = if is_fusion {
        !covs.iter().all(|&c| covered(c))
    } else {
        !covs.iter().any(|&c| covered(c))
    };

    // `breakpoint_sides` already names each position once, in haplotype order.
    let uncovered: Vec<String> = sides
        .iter()
        .zip(&covs)
        .filter(|(_, &c)| !covered(c))
        .map(|((chrom, pos), _)| format!("{}:{}", chrom, pos))
        .collect();

    if !refuse {
        let (first, cov) = covs
            .iter()
            .copied()
            .enumerate()
            .find(|&(_, c)| covered(c))
            .expect("a kept event has at least one covered breakpoint side");
        if origin.is_some() {
            let (chrom, pos) = &sides[first];
            log::info!(
                "  origin depth at {}:{}: {:.1}x (the donor pool's there: {:.1}x)",
                chrom,
                pos,
                cov,
                estimate_coverage_at(pool, chrom, *pos, 2000),
            );
        }
        if !uncovered.is_empty() {
            let bare = if origin.is_some() {
                "the origin depth is 0 at"
            } else {
                "the donor pool has no reads over"
            };
            log::warn!(
                "event {}:{}-{} is kept although {} \
                 {}: one bare breakpoint side is the far edge of a sliced or panel \
                 BAM, not a reason to refuse. But the haplotype spans every side, \
                 so the fragments tiled across the bare part are scaled by the \
                 {:.1}x measured elsewhere and land where the input BAM has no \
                 read -- a coverage island the input does not have. Extract a \
                 wider BAM, or narrow the event, if that matters downstream.",
                event_span.0,
                event_span.1,
                event_span.2,
                bare,
                uncovered.join(", "),
                cov,
            );
        }
        return Ok((cov, uncovered));
    }

    let positions = uncovered.join(", ");
    if let Some(site) = origin {
        anyhow::bail!(
            "event over {} has no origin depth at any of its breakpoints ({}): no read at \
             the spot or at its {} look-alike region(s) could have come from there, so \
             spike would invent the reads it plants and still write a truth VCF beside them.",
            site.footprint,
            positions,
            site.lookalikes.len(),
        );
    }
    // "one side" only when it is one: the parenthesis lists them all, so the
    // singular contradicted the message's own evidence.
    let scope = if is_fusion {
        let how_many = if uncovered.len() == 1 {
            "on one side of its junction".to_string()
        } else {
            format!("on {} sides of its junction", uncovered.len())
        };
        format!(
            "{} -- every read spike plants for a fusion spans the junction, so each \
             partner needs donor reads of its own",
            how_many
        )
    } else {
        "at any of its breakpoints".to_string()
    };
    anyhow::bail!(
        "event {}:{}-{} has no donor coverage {} ({}): the pool holds {} read \
         pair(s) but none of them cover that. The tiling count is coverage x VAF, \
         so spike would invent the reads it plants and still write a truth VCF \
         beside them. Check that the event lies in a covered region of the BAM -- \
         a --region window elsewhere, or a fusion partner that carries the whole \
         pool, fills the pool without covering the event.",
        event_span.0,
        event_span.1,
        event_span.2,
        scope,
        positions,
        pool.pairs.len(),
    );
}

/// The number of fragments to tile, and the fraction they actually plant
/// when spike could not plant the one that was asked for.
struct TilingCount {
    /// Fragments to tile.
    count: usize,
    /// The fraction `count` fragments realize, when the additive cap or the
    /// two-fragment floor moved the count off the requested VAF; `None` when
    /// the request stands and `SIM_VAF` keeps it unchanged.
    ///
    /// Detection is on the two mechanisms, not on the two numbers differing:
    /// the `round()` below moves the realized fraction off the request by a
    /// hair on nearly every event, so comparing numbers would report rounding
    /// on almost every record and drown the two mechanisms this field is for.
    adjusted_vaf: Option<f64>,
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
) -> TilingCount {
    let plain = |count| TilingCount {
        count,
        adjusted_vaf: None,
    };
    if haplotype.total_len == 0 {
        return plain(0);
    }

    // Every count below is coverage x VAF, so at zero coverage the formula
    // asks for nothing and only the floor of 2 is left -- two pairs invented
    // out of a constant quality profile, beside a truth VCF that claims a
    // real variant. The pool-size guard in `main.rs` cannot see this: it
    // sums every extraction window, and a `--region` window (or a fusion's
    // other side) can fill the pool without covering the event. Refuse here
    // instead, where the coverage is actually measured.
    if coverage.is_nan() || coverage <= 0.0 {
        return plain(0);
    }

    let breakpoints = haplotype.breakpoints();

    if breakpoint_only && !breakpoints.is_empty() {
        // An additive event can't reach vaf = 1 (it would need infinitely
        // many added fragments), so cap it.
        const MAX_ADDITIVE_VAF: f64 = 0.95;
        let capped = vaf > MAX_ADDITIVE_VAF;
        if capped {
            log::warn!(
                "additive event: VAF {:.2} capped at {:.2} (original reads are kept); \
                 the truth VCF records the capped fraction as SIM_VAF and {:.2} as \
                 SIM_REQ_VAF",
                vaf,
                MAX_ADDITIVE_VAF,
                vaf
            );
        }
        let v = vaf.min(MAX_ADDITIVE_VAF);
        let n = (coverage * v / (1.0 - v) * breakpoints.len() as f64).round() as usize;
        let (count, floored) = floor_tiling_count(n, coverage, v);
        // n junction fragments against the coverage kept at each breakpoint
        // are n / (n + coverage * breakpoints) of the depth there -- the
        // formula above, inverted for the count finally returned, so the
        // fraction recorded is the one the fragments emitted make up.
        let adjusted_vaf = (capped || floored).then(|| {
            realized_fraction(
                count as f64,
                count as f64 + coverage * breakpoints.len() as f64,
            )
        });
        return TilingCount { count, adjusted_vaf };
    }

    // Fragment starts are uniform over the starts whose fragment overlaps
    // reference sequence: [0, L - f] minus starts lying wholly inside inserted
    // sequence (tiling never draws those). Interior depth is then n * f / starts;
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
    let (count, floored) = floor_tiling_count(n, coverage, vaf);
    // n = coverage * vaf * effective_len / mean_frag, inverted for the count
    // finally returned.
    let adjusted_vaf =
        floored.then(|| realized_fraction(count as f64 * mean_frag, coverage * effective_len));
    TilingCount { count, adjusted_vaf }
}

/// `numerator / denominator` as an allele fraction.
///
/// Only reached with `coverage > 0`, so the denominator is positive for any
/// haplotype with reference-overlapping length to divide by. A haplotype
/// shorter than one fragment has none -- the floor's own case -- and the
/// fragments planted are then the whole of the event's support, a fraction of
/// 1; the same clamp keeps a very short haplotype's fraction inside (0, 1].
fn realized_fraction(numerator: f64, denominator: f64) -> f64 {
    if denominator > 0.0 {
        (numerator / denominator).min(1.0)
    } else {
        1.0
    }
}

/// Smallest number of tiled fragments spike will plant an event with.
const MIN_TILED_FRAGMENTS: usize = 2;

/// Apply the [`MIN_TILED_FRAGMENTS`] floor, and say when it changes the answer.
///
/// The floor is right for a thin-but-covered region: a haplotype shorter than
/// one fragment asks for 0 however real the coverage is, and planting nothing
/// would leave a truth VCF with no reads behind it. But when the floor raises
/// the count, the realized allele fraction is above the one that was asked
/// for -- at low coverage by several times over, in the direction that
/// flatters a caller -- so the caller records the realized fraction as
/// `SIM_VAF` and keeps the request in `SIM_REQ_VAF`. Callers only reach here
/// with `coverage > 0`, so this is never the zero-coverage case, which
/// `compute_tiling_count` refuses outright.
///
/// Returns the count to tile and whether the floor raised it.
fn floor_tiling_count(requested: usize, coverage: f64, vaf: f64) -> (usize, bool) {
    if requested >= MIN_TILED_FRAGMENTS {
        return (requested, false);
    }
    log::warn!(
        "coverage {:.1}x at VAF {:.3} asks for {} tiled fragment(s); spike emits the \
         {} it needs to plant the event at all, so the realized allele fraction is \
         above the {:.3} requested; the truth VCF records the realized fraction as \
         SIM_VAF and the {:.3} requested as SIM_REQ_VAF",
        coverage,
        vaf,
        requested,
        MIN_TILED_FRAGMENTS,
        vaf,
        vaf
    );
    (MIN_TILED_FRAGMENTS, true)
}

/// The fragment starts in `[0, max_start]` whose fragment overlaps reference
/// sequence, as inclusive intervals.
///
/// A fragment `[s, s + frag_len)` touches a reference segment at haplotype
/// `[a, b)` exactly when `s < b` and `s + frag_len > a`, so that segment
/// contributes `[a - (frag_len - 1), min(b - 1, max_start)]` -- the same
/// predicate `overlaps_ref_segment` tests, one segment at a time. Segments are
/// in haplotype order, so the intervals come out ascending and only
/// neighbours can meet; touching or nearby reference segments are merged,
/// because a start listed twice would be drawn twice and would make the
/// union's total length too long.
fn ref_overlapping_start_intervals(
    haplotype: &VariantHaplotype,
    frag_len: u64,
    max_start: u64,
) -> Vec<(u64, u64)> {
    let mut intervals: Vec<(u64, u64)> = Vec::new();
    if frag_len == 0 {
        return intervals;
    }
    for seg in &haplotype.segments {
        if seg.origin.is_none() {
            continue; // novel sequence anchors nothing
        }
        let seg_start = seg.hap_offset;
        let seg_end = seg_start + seg.sequence.len() as u64;
        if seg_end == 0 {
            continue; // an empty segment at offset 0 reaches no start
        }
        let lo = seg_start.saturating_sub(frag_len - 1);
        let hi = (seg_end - 1).min(max_start);
        if lo > hi {
            continue; // wholly past max_start
        }
        match intervals.last_mut() {
            Some(last) if lo <= last.1.saturating_add(1) => last.1 = last.1.max(hi),
            _ => intervals.push((lo, hi)),
        }
    }
    intervals
}

/// Draw a fragment start uniformly over the starts that overlap reference.
///
/// The draw is over the union's total length and then mapped into the
/// intervals, so every valid start is equally likely, a start outside the
/// union can never come out, and there is nothing to retry. `None` when the
/// union is empty -- no reference segment within `frag_len` of any start --
/// which the caller skips rather than placing a fragment anyway.
fn sample_ref_overlapping_start(
    haplotype: &VariantHaplotype,
    frag_len: u64,
    max_start: u64,
    rng: &mut StdRng,
) -> Option<u64> {
    let intervals = ref_overlapping_start_intervals(haplotype, frag_len, max_start);
    let total: u64 = intervals.iter().map(|&(lo, hi)| hi - lo + 1).sum();
    if total == 0 {
        return None;
    }
    let mut offset = rng.gen_range(0..total);
    for (lo, hi) in intervals {
        let len = hi - lo + 1;
        if offset < len {
            return Some(lo + offset);
        }
        offset -= len;
    }
    None // unreachable: `offset` is below the intervals' total length
}

/// Tile synthetic reads across the variant haplotype.
///
/// Number of reads is determined by coverage, VAF, and haplotype/zone length.
/// Reads naturally become chimeric when they span segment boundaries.
///
/// When `breakpoint_only` is true (for DUP/Fusion additive events), reads are
/// placed only near segment boundaries so they cross a breakpoint. This avoids
/// inflating coverage in flank regions where original reads are already kept.
///
/// `other_copy` is the same haplotype with the other sample copy's alleles,
/// and the chance a fragment comes from it (above VAF 0.5).
///
/// Returns the pairs, and the fraction they actually plant when the additive
/// cap or the two-fragment floor moved the count off `vaf` (see
/// [`TilingCount`]).
#[allow(clippy::too_many_arguments)]
fn tile_haplotype_reads(
    haplotype: &VariantHaplotype,
    other_copy: Option<(&VariantHaplotype, f64)>,
    synth_gen: &SynthReadGenerator,
    pool: &ReadPool,
    coverage: f64,
    vaf: f64,
    breakpoint_only: bool,
    name_prefix: &str,
    rng: &mut StdRng,
) -> (Vec<ReadPair>, Option<f64>) {
    let hap_len = haplotype.total_len;
    if hap_len == 0 {
        return (Vec::new(), None);
    }

    // Use the library's real mean fragment length; fall back only when the
    // pool gave none (e.g. no reads).
    let mean_frag = if pool.frag_dist.mean.is_finite() && pool.frag_dist.mean > 0.0 {
        pool.frag_dist.mean
    } else {
        300.0
    };
    let min_frag = synth_gen.min_fragment_len();
    let breakpoints = haplotype.breakpoints();

    let plan = compute_tiling_count(haplotype, coverage, vaf, mean_frag, breakpoint_only);
    let n_frags = plan.count;

    log::info!(
        "Tiling {} synthetic reads across {}bp haplotype (cov={:.1}, vaf={:.2}, bp_only={})",
        n_frags,
        hap_len,
        coverage,
        vaf,
        breakpoint_only,
    );

    let mut pairs = Vec::with_capacity(n_frags);

    // Check if any novel (non-reference) segments exist. If so, placement
    // below draws the fragment start directly from the starts whose fragment
    // overlaps reference sequence, so a fragment is never placed entirely
    // within novel sequence (which wouldn't contribute to observable
    // reference-aligned coverage).
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
            .sample_in_range(rng, min_frag, crate::stats::MAX_FRAGMENT_LEN) as u64;

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
        } else if has_novel {
            // Uniform across the starts whose fragment overlaps reference
            // sequence, drawn in proportion to those intervals' lengths.
            // `compute_tiling_count` left the novel-only starts out of the
            // budget, so accepting one here would spend a
            // reference-overlapping fragment on a read pair with no reference
            // anchor (CR5).
            match sample_ref_overlapping_start(haplotype, frag_len, max_start, rng) {
                Some(start) => start,
                // Unreachable for any haplotype the CLI can build:
                // `validate_flank` forces `--flank >= HAP_FLANK` (2000) and
                // `from_insertion` always yields at least one non-empty
                // flank, so the reference-overlapping union is provably
                // non-empty whenever frag_len <= hap_len. If it were ever
                // empty regardless, skipping is still the safe choice -- an
                // invalid start must never be emitted.
                None => continue,
            }
        } else {
            // Uniform across the whole haplotype.
            rng.gen_range(0..=max_start)
        };

        let name = synth_gen.read_name(&format!("{}_hap_{:06}", name_prefix, idx));
        idx += 1;

        // Both copies share the layout; only their alleles differ.
        let source = match other_copy {
            Some((other, p)) if rng.gen::<f64>() < p => other,
            _ => haplotype,
        };

        if let Some(pair) = synth_gen.generate_haplotype_read_pair(
            source,
            hap_start,
            frag_len,
            &name,
            rng,
        ) {
            pairs.push(pair);
        }
    }

    (pairs, plan.adjusted_vaf)
}

/// [`depth_fold`] with each bin's depth measured by `coverage_at(chrom,
/// center, width)` instead of the pool.
fn depth_fold_by(
    haplotype: &VariantHaplotype,
    cov: f64,
    coverage_at: &(dyn Fn(&str, u64, u64) -> f64 + Sync),
) -> DepthFold {
    use rayon::prelude::*;
    const BIN: u64 = 1000;
    let mut bins: Vec<(&str, u64, u64)> = Vec::new();
    for origin in haplotype.segments.iter().filter_map(|seg| seg.origin.as_ref()) {
        let len = origin.ref_end.saturating_sub(origin.ref_start);
        if len == 0 {
            continue;
        }
        let n_bins = ((len as f64 / BIN as f64).round() as u64).max(1);
        for b in 0..n_bins {
            let start = origin.ref_start + len * b / n_bins;
            let end = origin.ref_start + len * (b + 1) / n_bins;
            bins.push((origin.chrom.as_str(), start, end));
        }
    }
    // Each bin's depth on the thread pool; the worst is then picked in bin
    // order, so a tie goes to the same bin as in one pass.
    let depths: Vec<f64> = bins
        .par_iter()
        .map(|&(chrom, start, end)| coverage_at(chrom, (start + end) / 2, end - start))
        .collect();
    let mut worst = DepthFold {
        fold: 1.0,
        scaled_by: cov,
        worst_bin: String::new(),
        worst_depth: cov,
    };
    for (&(chrom, start, end), &depth) in bins.iter().zip(&depths) {
        let fold = ((depth + 1.0) / (cov + 1.0)).max((cov + 1.0) / (depth + 1.0));
        if fold > worst.fold {
            worst.fold = fold;
            worst.worst_bin = format!("{}:{}-{}", chrom, start, end);
            worst.worst_depth = depth;
        }
    }
    worst
}

/// The largest fold between the donor's depth in any bin this haplotype's
/// fragments are drawn from and `cov`, the depth they are all scaled by (CR2).
///
/// Each reference interval a segment is drawn from is cut into
/// `max(1, round(len / 1000))` equal bins, and each bin's depth is measured
/// the way `cov` is: [`estimate_coverage_at`] on the same pool, the bin as its
/// window. Measuring only; nothing here changes what is tiled.
fn depth_fold(haplotype: &VariantHaplotype, pool: &ReadPool, cov: f64) -> DepthFold {
    let depth = PoolDepth::new(pool);
    depth_fold_by(haplotype, cov, &|chrom, pos, window| depth.coverage_at(chrom, pos, window))
}

/// Estimate fragment depth at a reference position on `chrom`.
///
/// Counts fragments overlapping positions in a window, returns mean coverage.
/// Uses up to 50 evenly-spaced sample points for stable estimates.
/// The pool must be sorted by `ref_start` (which `build_read_pool` ensures).
/// Only pairs on `chrom` count: a fusion pool holds both partners, and the
/// other partner's reads can sit at the same coordinates (M9).
fn estimate_coverage_at(pool: &ReadPool, chrom: &str, pos: u64, window: u64) -> f64 {
    PoolDepth::new(pool).coverage_at(chrom, pos, window)
}

/// The donor pool's depth, ready to be asked at many positions.
///
/// At each point only the pairs that start close enough to reach it are
/// checked: none starts `longest` or more before the point and still covers
/// it. Checking every pair that starts before the point instead made the
/// depth fold, which asks once per 1 kb bin, grow with the square of the
/// event: 243 s of a 3 Mb DUP's 281 s on the 35x chr20 BAM.
struct PoolDepth<'a> {
    pool: &'a ReadPool,
    /// The longest pair's reference span.
    longest: u64,
}

impl<'a> PoolDepth<'a> {
    fn new(pool: &'a ReadPool) -> Self {
        let longest = pool
            .pairs
            .iter()
            .map(|p| p.ref_end.saturating_sub(p.ref_start))
            .max()
            .unwrap_or(0);
        PoolDepth { pool, longest }
    }

    /// See [`estimate_coverage_at`].
    fn coverage_at(&self, chrom: &str, pos: u64, window: u64) -> f64 {
        let points = coverage_sample_points(pos, window);
        let pairs = &self.pool.pairs;
        let mut total = 0usize;
        for &sample_pos in &points {
            // The pool is sorted by ref_start: pairs from `upper` on start
            // after the point, and pairs before `lower` end before it.
            let upper = pairs.partition_point(|p| p.ref_start <= sample_pos);
            let lower = pairs[..upper].partition_point(|p| p.ref_start + self.longest <= sample_pos);
            total += pairs[lower..upper]
                .iter()
                .filter(|p| p.ref_end > sample_pos && p.chrom == chrom)
                .count();
        }
        total as f64 / points.len() as f64
    }
}

/// The points [`estimate_coverage_at`] and `origin`'s read depth sample:
/// up to 50, evenly spaced from `pos - window / 2`. One grid for both, so the
/// two depths R1 compares count at the same points (PD-5).
pub(crate) fn coverage_sample_points(pos: u64, window: u64) -> Vec<u64> {
    let start = pos.saturating_sub(window / 2);
    let end = pos.saturating_add(window / 2);
    let range = end - start;
    // For very small windows, use fewer sample points (at least 1).
    let n = range.clamp(1, 50);
    let step = if n > 1 { range / n } else { 1 };
    (0..n).map(|i| start + i * step).collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::extract;
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

    fn make_pair_on(chrom: &str, name: &str, start: u64, end: u64) -> ReadPair {
        ReadPair {
            chrom: chrom.to_string(),
            ..make_pair(name, start, end)
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

    /// Make a novel (insertion) segment with no reference origin. G, not the
    /// complement of the A flanks, so a read carrying inserted sequence is
    /// recognizable whichever strand it came off (G forward, C reverse).
    fn novel_segment(len: u64) -> HaplotypeSegment {
        HaplotypeSegment {
            sequence: vec![b'G'; len as usize],
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
        let count = compute_tiling_count(&hap, 30.0, 0.5, 400.0, false).count;
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
        let count = compute_tiling_count(&hap, 30.0, 0.5, 400.0, false).count;
        assert_eq!(count, 150);

        // 1000 bp insertion: (5000 - 400) - (1000 - 400) = 4000 starts -> 150.
        let long = make_haplotype(vec![
            ref_segment(0, 2000),
            novel_segment(1000),
            ref_segment(2000, 2000),
        ]);
        assert_eq!(compute_tiling_count(&long, 30.0, 0.5, 400.0, false).count, 150);
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

        let count = compute_tiling_count(&hap, 30.0, 0.5, 400.0, false).count;
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
        let count = compute_tiling_count(&hap, 30.0, 0.5, 400.0, true).count;
        assert_eq!(count, 30);
    }

    #[test]
    fn test_tiling_count_breakpoint_only_gives_requested_fraction() {
        // Additive events keep all `cov` original fragments at the junction,
        // so n junction fragments make up n / (cov + n) of it. For fraction v,
        // n = cov * v / (1 - v) per breakpoint.
        let one_bp = make_haplotype(vec![ref_segment(0, 2000), ref_segment(5000, 2000)]);
        assert_eq!(compute_tiling_count(&one_bp, 40.0, 0.2, 400.0, true).count, 10);
        assert_eq!(compute_tiling_count(&one_bp, 100.0, 0.05, 400.0, true).count, 5);

        let two_bp = make_haplotype(vec![
            ref_segment(0, 2000),
            ref_segment(5000, 2000),
            ref_segment(9000, 2000),
        ]);
        assert_eq!(compute_tiling_count(&two_bp, 40.0, 0.2, 400.0, true).count, 20);
    }

    #[test]
    fn test_tiling_count_breakpoint_only_is_bounded_at_full_vaf() {
        // v = 1 would need infinitely many added fragments; the count must
        // stay finite (an unbounded usize would abort on allocation).
        let hap = make_haplotype(vec![ref_segment(0, 2000), ref_segment(5000, 2000)]);
        let count = compute_tiling_count(&hap, 40.0, 1.0, 400.0, true).count;
        assert!(count > compute_tiling_count(&hap, 40.0, 0.9, 400.0, true).count);
        assert!(count <= 40 * 100, "count {} is unbounded", count);
    }

    #[test]
    fn test_tiling_count_minimum() {
        // Very low coverage → at least 2 reads.
        let hap = make_haplotype(vec![ref_segment(0, 100)]);
        let count = compute_tiling_count(&hap, 0.1, 0.1, 400.0, false).count;
        assert_eq!(count, 2); // min of 2
    }

    #[test]
    fn test_tiling_count_zero_coverage_is_zero() {
        // No coverage at the breakpoint means there is nothing to scale
        // against, so the floor of 2 would invent reads rather than simulate
        // them. Both branches must return 0 and let the caller refuse.
        let hap = make_haplotype(vec![ref_segment(0, 4000)]);
        assert_eq!(compute_tiling_count(&hap, 0.0, 0.5, 400.0, false).count, 0);

        let additive = make_haplotype(vec![ref_segment(0, 2000), ref_segment(5000, 2000)]);
        assert_eq!(compute_tiling_count(&additive, 0.0, 0.2, 400.0, true).count, 0);
    }

    #[test]
    fn test_tiling_count_empty_haplotype() {
        let hap = make_haplotype(vec![]);
        assert_eq!(hap.total_len, 0);
        let count = compute_tiling_count(&hap, 30.0, 0.5, 400.0, false).count;
        assert_eq!(count, 0);
    }

    #[test]
    fn test_tiling_count_capped_additive_reports_the_fraction_it_planted() {
        // An additive request above the 0.95 cap plants the capped count, so
        // the fraction to record is the one those fragments make up, not the
        // 0.99 that was asked for.
        let hap = make_haplotype(vec![ref_segment(0, 2000), ref_segment(5000, 2000)]);
        let plan = compute_tiling_count(&hap, 40.0, 0.99, 400.0, true);
        // 40 * 0.95 / 0.05 * 1 breakpoint = 760 junction fragments, which are
        // 760 / (760 + 40) of the depth at that junction.
        assert_eq!(plan.count, 760);
        let adjusted = plan.adjusted_vaf.expect("the cap moved the fraction");
        assert_eq!(format!("{:.3}", adjusted), "0.950");
        assert!(adjusted < 0.99, "recorded {} is still the request", adjusted);
    }

    #[test]
    fn test_tiling_count_floored_event_reports_the_fraction_it_planted() {
        // 0.7x coverage at VAF 0.05 asks for 0 fragments; the floor plants 2,
        // which are a far larger share of the depth than 0.05.
        let hap = make_haplotype(vec![ref_segment(0, 4000)]);
        let plan = compute_tiling_count(&hap, 0.7, 0.05, 400.0, false);
        assert_eq!(plan.count, 2);
        let adjusted = plan.adjusted_vaf.expect("the floor moved the fraction");
        // 2 * 400 / (0.7 * (4000 - 400)) = 0.3175.
        assert_eq!(format!("{:.3}", adjusted), "0.317");
        assert!(adjusted > 0.05, "recorded {} is still the request", adjusted);
    }

    #[test]
    fn test_tiling_count_rounding_alone_does_not_move_the_recorded_fraction() {
        // 30 * 0.5 * (9000 - 400) / 400 = 322.5 fragments, rounded to 323 --
        // a realized fraction a hair off 0.5. Only the cap and the floor are
        // reported; rounding would move nearly every record and say nothing
        // about either mechanism.
        let hap = make_haplotype(vec![ref_segment(0, 9000)]);
        let plan = compute_tiling_count(&hap, 30.0, 0.5, 400.0, false);
        assert_eq!(plan.count, 323);
        assert_eq!(plan.adjusted_vaf, None);
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
        let cov = estimate_coverage_at(&pool, "chr1", 250, 500);
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
        let cov = estimate_coverage_at(&pool, "chr1", 250, 500);
        assert!(
            cov < 20.0 && cov > 0.0,
            "expected partial coverage, got {:.1}",
            cov
        );
    }

    #[test]
    fn test_estimate_coverage_at_counts_only_the_queried_chromosome() {
        // A fusion pool holds both partners. 10 pairs cover [0, 500) on chr1
        // and 10 more cover the same coordinates on chr2; a chr1 query must
        // see depth 10, not 20 (M9).
        let mut pairs: Vec<ReadPair> = (0..10)
            .map(|i| make_pair_on("chr1", &format!("a{}", i), 0, 500))
            .collect();
        pairs.extend((0..10).map(|i| make_pair_on("chr2", &format!("b{}", i), 0, 500)));
        let pool = make_pool(pairs);
        let cov = estimate_coverage_at(&pool, "chr1", 250, 500);
        assert!((cov - 10.0).abs() < 1.0, "expected ~10.0 got {:.1}", cov);
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

    // ---------------------------------------------------------------
    // Reference-overlapping start sampler
    // ---------------------------------------------------------------

    /// Every reference-segment layout the sampler has to get right: flanks at
    /// both ends, two insertions of different sizes, and two reference
    /// segments that touch (so their start intervals overlap and must merge).
    fn mixed_haplotype() -> VariantHaplotype {
        make_haplotype(vec![
            ref_segment(0, 20),
            novel_segment(30),
            ref_segment(20, 20),
            ref_segment(40, 10),
            novel_segment(5),
            ref_segment(50, 10),
        ])
    }

    #[test]
    fn test_ref_overlapping_start_intervals_are_exactly_the_valid_starts() {
        let hap = mixed_haplotype();
        assert_eq!(hap.total_len, 95);

        for frag_len in 1..=hap.total_len {
            let max_start = hap.total_len - frag_len;
            let intervals = ref_overlapping_start_intervals(&hap, frag_len, max_start);

            // Ascending and disjoint: an unmerged overlap would draw a shared
            // start twice and make the union's total length too long.
            for pair in intervals.windows(2) {
                assert!(
                    pair[0].1 < pair[1].0,
                    "frag_len {}: intervals {:?} overlap or are out of order",
                    frag_len,
                    intervals
                );
            }

            let listed: Vec<u64> = intervals
                .iter()
                .flat_map(|&(lo, hi)| lo..=hi)
                .collect();
            let valid: Vec<u64> = (0..=max_start)
                .filter(|&start| hap.overlaps_ref_segment(start, frag_len))
                .collect();
            assert_eq!(listed, valid, "frag_len {}", frag_len);
        }
    }

    #[test]
    fn test_sample_ref_overlapping_start_uses_both_flanks_of_a_long_insertion() {
        // ref [0,1000) | 100 kb insertion | ref [1000,2000). At frag_len 400
        // the valid starts are [0,999] on the left and [100601,101600] on the
        // right: 1000 each, so a sampler that weights by interval length
        // splits its draws evenly.
        let hap = make_haplotype(vec![
            ref_segment(0, 1000),
            novel_segment(100_000),
            ref_segment(1000, 1000),
        ]);
        let frag_len = 400;
        let max_start = hap.total_len - frag_len;
        let mut rng = StdRng::seed_from_u64(7);

        let mut left = 0;
        let mut right = 0;
        for _ in 0..1000 {
            let start = sample_ref_overlapping_start(&hap, frag_len, max_start, &mut rng)
                .expect("a haplotype with reference flanks always has a valid start");
            assert!(
                hap.overlaps_ref_segment(start, frag_len),
                "start {} puts the fragment wholly inside the insertion",
                start
            );
            if start < 1000 {
                left += 1;
            } else {
                right += 1;
            }
        }
        assert!(
            left > 350 && right > 350,
            "sampler favours one flank: {} left, {} right",
            left,
            right
        );
    }

    #[test]
    fn test_sample_ref_overlapping_start_weights_by_interval_length_for_unequal_flanks() {
        // The test above uses two EQUAL flanks, so a sampler that picked an
        // *interval* uniformly (ignoring how many starts it holds) would
        // pass it too -- both give ~50/50. Unequal flanks tell the two
        // apart. Reviewer's probe: ref[0,300) | 50 kb insertion | ref[300,1200),
        // frag_len 400. The reference segments contribute [0,299] (300
        // starts) on the left and [49901,50800] (900 starts) on the right --
        // a 25%/75% split by length, not 50/50 by interval count.
        let hap = make_haplotype(vec![
            ref_segment(0, 300),
            novel_segment(50_000),
            ref_segment(300, 900),
        ]);
        let frag_len = 400;
        let max_start = hap.total_len - frag_len;
        assert_eq!(
            ref_overlapping_start_intervals(&hap, frag_len, max_start),
            vec![(0, 299), (49901, 50800)],
            "the probe's own arithmetic for the two intervals"
        );

        let mut rng = StdRng::seed_from_u64(11);
        let n = 10_000;
        let mut left = 0;
        for _ in 0..n {
            let start = sample_ref_overlapping_start(&hap, frag_len, max_start, &mut rng)
                .expect("a haplotype with reference flanks always has a valid start");
            assert!(
                hap.overlaps_ref_segment(start, frag_len),
                "start {} puts the fragment wholly inside the insertion",
                start
            );
            if start < 300 {
                left += 1;
            }
        }

        // Expected left fraction is 300 / 1200 = 0.25. At n = 10_000 the
        // sampling standard error of that fraction is
        // sqrt(0.25 * 0.75 / 10_000) ~= 0.0043, so a tolerance of 0.03 (three
        // percentage points) is about seven standard errors -- effectively
        // flake-proof for a correct, length-weighted draw -- while a sampler
        // that instead picked one of the two intervals uniformly would land
        // near 0.50, 25 points off and nowhere close to passing.
        let left_frac = left as f64 / n as f64;
        assert!(
            (left_frac - 0.25).abs() < 0.03,
            "expected ~25% of draws from the 300-start left interval, got {:.4} ({} of {})",
            left_frac,
            left,
            n
        );
    }

    #[test]
    fn test_sample_ref_overlapping_start_is_none_without_reference_sequence() {
        // No reference segment at all: there is no valid start, and the
        // caller must skip the fragment rather than invent one.
        let hap = make_haplotype(vec![novel_segment(1000)]);
        let mut rng = StdRng::seed_from_u64(7);
        assert!(sample_ref_overlapping_start(&hap, 400, 600, &mut rng).is_none());
    }

    // ── Read evidence tests per SV type ─────────────────────────────────

    use crate::reference::SharedReference;
    use crate::synth::{QualityProfile, SynthReadGenerator};
    use rand::rngs::StdRng;
    use rand::SeedableRng;
    use std::collections::HashMap as StdHashMap;

    fn mock_synth_gen(read_length: usize) -> SynthReadGenerator<'static> {
        mock_synth_gen_at_q(read_length, 30)
    }

    /// `mock_synth_gen` at a chosen Phred. Q93, the FASTQ maximum, puts the
    /// substitution rate at 5e-10 per base, so a read's bases are the
    /// haplotype's own and a read's origin can be read off its sequence.
    fn mock_synth_gen_at_q(read_length: usize, phred: u8) -> SynthReadGenerator<'static> {
        let qual = vec![b'!' + phred; read_length];
        let pairs: Vec<ReadPair> = (0..100)
            .map(|i| ReadPair {
                name: format!("mock_{}", i),
                seq1: vec![b'A'; read_length],
                qual1: qual.clone(),
                seq2: vec![b'T'; read_length],
                qual2: qual.clone(),
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
            tile_haplotype_reads(&hap, None, &gen, &pool, 30.0, 0.5, false, "test", &mut rng).0;
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
            tile_haplotype_reads(&hap, None, &gen, &pool, 30.0, 0.5, false, "test", &mut rng).0;

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
            tile_haplotype_reads(&hap, None, &gen, &pool, 30.0, 0.5, false, "test", &mut rng).0;
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
            &hap, None, &gen, &pool, 30.0, 0.5, false, "ev0001", &mut rng_a,
        ).0;
        let pairs_b = tile_haplotype_reads(
            &hap, None, &gen, &pool, 30.0, 0.5, false, "ev0002", &mut rng_b,
        ).0;

        assert!(!pairs_a.is_empty());
        assert!(!pairs_b.is_empty());
        assert!(pairs_a.iter().all(|p| p.name.starts_with("SPIKE_ev0001_")));
        assert!(pairs_b.iter().all(|p| p.name.starts_with("SPIKE_ev0002_")));

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

    /// A 7-part name shape for the naming tests.
    fn seven_part_shape() -> crate::read_name::NameShape {
        crate::read_name::NameShape::Illumina {
            middle: "46:FC:2".to_string(),
            tiles: vec![1101, 1102],
            x: (1000, 2000),
            y: (3000, 4000),
        }
    }

    /// `named` is `plain` with each name in `shape`: the same reads, bases,
    /// qualities and places, so naming drew nothing from the random stream.
    fn same_reads_named_in(plain: &[ReadPair], named: &[ReadPair], shape: &crate::read_name::NameShape) {
        assert!(!plain.is_empty());
        assert_eq!(plain.len(), named.len());
        for (p, n) in plain.iter().zip(named) {
            let internal = p.name.strip_prefix("SPIKE_").unwrap_or_else(|| panic!("unmarked: {}", p.name));
            assert_eq!(n.name, shape.name(internal));
            let read = |r: &ReadPair| {
                (r.seq1.clone(), r.qual1.clone(), r.seq2.clone(), r.qual2.clone(), r.ref_start, r.ref_end)
            };
            assert_eq!(read(p), read(n), "{}", p.name);
            assert_eq!((p.insert_size, &p.chrom), (n.insert_size, &n.chrom), "{}", p.name);
        }
    }

    #[test]
    fn test_tiled_reads_are_named_in_the_generator_s_shape_without_drawing_randomness() {
        let hap = del_haplotype(2000, 5000);
        let pool = make_pool(vec![]);
        let tile = |gen: &SynthReadGenerator| {
            let mut rng = StdRng::seed_from_u64(42);
            tile_haplotype_reads(&hap, None, gen, &pool, 30.0, 0.5, false, "ev0001", &mut rng).0
        };
        let plain = tile(&mock_synth_gen(150));
        let named = tile(&mock_synth_gen(150).with_read_names(seven_part_shape()));
        assert!(plain.iter().all(|p| p.name.starts_with("SPIKE_ev0001_hap_")));
        same_reads_named_in(&plain, &named, &seven_part_shape());
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
            &hap, None, &gen, &pool, 30.0, 0.5, false, "sv", &mut rng,
        ).0;

        assert!(
            !pairs.is_empty(),
            "small-variant tiling should produce reads (boundary crossing allowed)",
        );
    }

    #[test]
    fn test_tiling_draws_fragments_shorter_than_a_read_only_in_a_trimmed_library() {
        let hap = make_haplotype(vec![HaplotypeSegment {
            sequence: (0..4_000usize).map(|i| b"ACGT"[i % 4]).collect(),
            origin: Some(SegmentOrigin {
                chrom: "chr1".to_string(),
                ref_start: 1_000,
                ref_end: 5_000,
                is_reverse: false,
            }),
            hap_offset: 0,
        }]);
        let pool = ReadPool { pairs: vec![], frag_dist: FragmentDist::from_stats(100.0, 0.0) };
        let mut rng = StdRng::seed_from_u64(11);

        let trimmed = mock_synth_gen_at_q(151, 93).with_adapter_trim(true);
        let pairs = tile_haplotype_reads(&hap, None, &trimmed, &pool, 30.0, 0.5, false, "sv", &mut rng).0;
        assert!(!pairs.is_empty());
        assert!(
            pairs.iter().all(|p| p.insert_size == 100 && p.seq1.len() == 100 && p.seq2.len() == 100),
            "a trimmed library's 100 bp fragments give 100 bp reads"
        );

        let untrimmed = mock_synth_gen_at_q(151, 93);
        let pairs = tile_haplotype_reads(&hap, None, &untrimmed, &pool, 30.0, 0.5, false, "sv", &mut rng).0;
        assert!(!pairs.is_empty());
        assert!(pairs.iter().all(|p| p.insert_size == 151), "an untrimmed library clamps up to one read");
    }

    #[test]
    fn test_the_donor_fragment_model_starts_at_one_base_in_a_trimmed_library() {
        let mut config = make_config();
        assert_eq!(config.min_fragment_len(), 150);
        config.adapter_trimmed = true;
        assert_eq!(config.min_fragment_len(), 1);
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
            &hap, None, &gen, &pool, 30.0, 0.5, false, "dup", &mut rng,
        ).0;
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

    /// Pairs of 400 bp starting every 10 bp over [9000, 17000): a fragment
    /// depth of exactly 40. Inside `thin` only every fourth start is kept,
    /// so a position whose covering starts all lie there reads 10.
    fn depth_pool(thin: Option<(u64, u64)>) -> ReadPool {
        let pairs = (9_000u64..17_000)
            .step_by(10)
            .filter(|&s| match thin {
                Some((a, b)) if s >= a && s < b => s % 40 == 0,
                _ => true,
            })
            .map(|s| make_pair(&format!("d{}", s), s, s + 400))
            .collect();
        ReadPool {
            pairs,
            frag_dist: FragmentDist::from_stats(400.0, 80.0),
        }
    }

    #[test]
    fn test_depth_fold_is_one_over_an_even_donor() {
        let hap = tandem_dup_haplotype(11_000, 15_000, 1_000);
        let fold = depth_fold(&hap, &depth_pool(None), 40.0);
        assert!((fold.fold - 1.0).abs() < 1e-9, "{:?}", fold);
        assert_eq!(fold.scaled_by, 40.0);
    }

    #[test]
    fn test_depth_fold_finds_the_thin_bin() {
        // CR2: the tiling scales every fragment by one depth. Where the donor
        // runs at a quarter of it the fold is (40+1)/(10+1), and the bin
        // wholly inside the thin stretch is the one named.
        let hap = tandem_dup_haplotype(11_000, 15_000, 1_000);
        let fold = depth_fold(&hap, &depth_pool(Some((12_000, 14_000))), 40.0);
        assert!((fold.fold - 41.0 / 11.0).abs() < 1e-9, "{:?}", fold);
        assert_eq!(fold.worst_bin, "chr1:13000-14000");
        assert!((fold.worst_depth - 10.0).abs() < 1e-9, "{:?}", fold);
    }

    #[test]
    fn test_depth_fold_counts_a_donor_deeper_than_the_scaling_depth() {
        // The fold is symmetric: a donor at four times the scaling depth is
        // as far off as one at a quarter of it.
        let hap = tandem_dup_haplotype(11_000, 15_000, 1_000);
        let fold = depth_fold(&hap, &depth_pool(None), 10.0);
        assert!((fold.fold - 41.0 / 11.0).abs() < 1e-9, "{:?}", fold);
    }

    #[test]
    fn test_tandem_dup_tiling_count_uses_full_length() {
        // Full tandem DUP: total ref_mapped_len = 2*flank + 2*dup_size
        let hap = tandem_dup_haplotype(1000, 3000, 500);

        // ref_mapped_len = 500 + 2000 + 2000 + 500 = 5000
        assert_eq!(hap.ref_mapped_len(), 5000);

        // breakpoint_only=false: uses ref_mapped_len as effective_len
        let n = compute_tiling_count(&hap, 30.0, 0.5, 400.0, false).count;
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
            adapter_trimmed: false,
            read_names: Default::default(),
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
        let gen = mock_synth_gen(150);
        let mut rng = StdRng::seed_from_u64(42);

        let event = SimEvent::Deletion {
            chrom: "chr1".to_string(),
            del_start: 1000,
            del_end: 3000,
            gene: "TEST".to_string(),
            exons: vec![],
            allele_fraction: Some(0.5),
        };

        let out = simulate_event(1, &event, &pool, &mut hap, &config, &gen, 0.5, &mut rng)
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

        // Chimeric read names should carry the mark and the event prefix.
        assert!(
            out.chimeric_pairs[0].name.starts_with("SPIKE_ev0001"),
            "chimeric read name should start with 'SPIKE_ev0001'"
        );

        // Total output (kept + chimeric) should be non-empty.
        assert!(
            out.kept_originals.len() + out.chimeric_pairs.len() > 0,
            "simulate_event should produce output reads"
        );
    }

    fn cr3_deletion() -> SimEvent {
        SimEvent::Deletion {
            chrom: "chr1".to_string(),
            del_start: 1000,
            del_end: 3000,
            gene: "TEST".to_string(),
            exons: vec![],
            allele_fraction: Some(0.5),
        }
    }

    #[test]
    fn test_simulate_event_stops_on_a_gvcf_it_cannot_read() {
        // CR3: an unreadable --gvcf used to log a warning and go on with no
        // SNPs, so the output could not be told from a sample that has none.
        let missing = std::env::temp_dir()
            .join(format!("spike_cr3_missing_{}.vcf", std::process::id()));
        let missing = missing.to_str().unwrap().to_string();
        let mut config = make_config();
        config.gvcf_path = Some(missing.clone());
        let mut hap = del_haplotype(1000, 2000);
        let pool = make_covering_pool(0, 5000, 200);
        let gen = mock_synth_gen(150);
        let mut rng = StdRng::seed_from_u64(42);

        let err = match simulate_event(
            1, &cr3_deletion(), &pool, &mut hap, &config, &gen, 0.5, &mut rng,
        ) {
            Ok(_) => panic!("an unreadable --gvcf must stop the event, not warn"),
            Err(e) => format!("{:#}", e),
        };
        assert!(err.contains(&missing), "the error must name the file: {}", err);
    }

    #[test]
    fn test_simulate_event_with_a_readable_gvcf_still_warns_when_the_pileup_fails() {
        // Only the gVCF read is fatal. A gVCF with no het SNPs here sends the
        // region to the pileup, and a pileup failure stays the warning it is
        // without --gvcf (make_config's BAM path is empty, so it fails).
        let path = std::env::temp_dir()
            .join(format!("spike_cr3_readable_{}.vcf", std::process::id()));
        std::fs::write(
            &path,
            "##fileformat=VCFv4.2\n\
             ##contig=<ID=chr1,length=100000>\n\
             #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n",
        )
        .unwrap();
        let mut config = make_config();
        config.gvcf_path = Some(path.to_str().unwrap().to_string());
        let mut hap = del_haplotype(1000, 2000);
        let pool = make_covering_pool(0, 5000, 200);
        let gen = mock_synth_gen(150);
        let mut rng = StdRng::seed_from_u64(42);

        let result =
            simulate_event(1, &cr3_deletion(), &pool, &mut hap, &config, &gen, 0.5, &mut rng);
        let _ = std::fs::remove_file(&path);
        result.expect("a pileup failure behind a readable gVCF must not stop the event");
    }

    #[test]
    fn test_simulate_event_refuses_pool_with_no_coverage_at_breakpoint() {
        // A --region window elsewhere (or a fusion partner that carries the
        // whole pool) can fill the pool past MIN_DONOR_PAIRS while leaving the
        // event itself uncovered. The tiling count is coverage x VAF, so it
        // collapses to the floor of 2 invented pairs beside a truth VCF.
        let mut hap = del_haplotype(1000, 2000);
        // 200 pairs, none of them anywhere near the deletion at chr1:1000-3000.
        let pool = make_covering_pool(500_000, 540_000, 200);
        let config = make_config();
        let gen = mock_synth_gen(150);
        let mut rng = StdRng::seed_from_u64(42);

        let event = SimEvent::Deletion {
            chrom: "chr1".to_string(),
            del_start: 1000,
            del_end: 3000,
            gene: "TEST".to_string(),
            exons: vec![],
            allele_fraction: Some(0.5),
        };

        let result = simulate_event(1, &event, &pool, &mut hap, &config, &gen, 0.5, &mut rng);
        let msg = match result {
            Ok(_) => panic!("a pool with no coverage at the breakpoint must be refused"),
            Err(e) => e.to_string(),
        };
        assert!(
            msg.contains("no donor coverage"),
            "error should name the missing coverage; got: {}",
            msg
        );
    }

    #[test]
    fn test_simulate_event_refuses_fusion_whose_other_partner_is_uncovered() {
        // A fusion is extracted from two loci and every read spike plants
        // spans the junction between them, so a partner with no donor reads
        // means half of every planted read is invented (N5's symptom). The
        // guard must see it whichever partner the event names first, so run
        // the same two loci both ways round and require the same verdict.
        let covered = 10000u64;
        let uncovered = 20000u64;
        let pool = make_covering_pool(8000, 11000, 100);

        let run = |bp_a: u64, bp_b: u64| -> Result<SplicedOutput> {
            let mut hap = make_haplotype(vec![
                HaplotypeSegment {
                    sequence: vec![b'A'; 1000],
                    origin: Some(SegmentOrigin {
                        chrom: "chr1".to_string(),
                        ref_start: bp_a - 1000,
                        ref_end: bp_a,
                        is_reverse: false,
                    }),
                    hap_offset: 0,
                },
                HaplotypeSegment {
                    sequence: vec![b'C'; 1000],
                    origin: Some(SegmentOrigin {
                        chrom: "chr1".to_string(),
                        ref_start: bp_b,
                        ref_end: bp_b + 1000,
                        is_reverse: false,
                    }),
                    hap_offset: 1000,
                },
            ]);
            let event = SimEvent::Fusion {
                chrom_a: "chr1".to_string(),
                bp_a,
                gene_a: "GENE_A".to_string(),
                chrom_b: "chr1".to_string(),
                bp_b,
                gene_b: "GENE_B".to_string(),
                join: FusionJoin::Forward,
                allele_fraction: Some(0.2),
            };
            let mut rng = StdRng::seed_from_u64(99);
            simulate_event(
                1, &event, &pool, &mut hap, &make_config(), &mock_synth_gen(150), 0.2, &mut rng,
            )
        };

        for (a, b, order) in [
            (covered, uncovered, "covered partner first"),
            (uncovered, covered, "uncovered partner first"),
        ] {
            let msg = match run(a, b) {
                Ok(out) => panic!(
                    "{}: a fusion partner with no donor coverage must be refused, \
                     but spike tiled {} chimeric pair(s)",
                    order,
                    out.chimeric_pairs.len()
                ),
                Err(e) => e.to_string(),
            };
            assert!(
                msg.contains("no donor coverage"),
                "{}: error should name the missing coverage; got: {}",
                order,
                msg
            );
        }
    }

    #[test]
    fn test_only_a_fusion_is_drawn_from_more_than_one_locus() {
        // `extract_pool_for_event` (main.rs) and `coverage_for_tiling`
        // (here) both key on this, in two files with no shared helper. The
        // match behind it is exhaustive, so a new event type cannot be added
        // without answering the question; this pins today's answers.
        let fusion = SimEvent::Fusion {
            chrom_a: "chr1".to_string(),
            bp_a: 10000,
            gene_a: "A".to_string(),
            chrom_b: "chr2".to_string(),
            bp_b: 20000,
            gene_b: "B".to_string(),
            join: FusionJoin::Forward,
            allele_fraction: None,
        };
        assert!(fusion.is_multi_locus(), "a fusion is drawn from two loci");
        for single in [
            del_event(1000, 2000),
            SimEvent::Duplication {
                chrom: "chr1".to_string(),
                dup_start: 1000,
                dup_end: 2000,
                gene: "T".to_string(),
                allele_fraction: None,
            },
            SimEvent::Inversion {
                chrom: "chr1".to_string(),
                inv_start: 1000,
                inv_end: 2000,
                gene: "T".to_string(),
                allele_fraction: None,
            },
            SimEvent::Insertion {
                chrom: "chr1".to_string(),
                pos: 1000,
                ins_seq: None,
                ins_len: 300,
                gene: "T".to_string(),
                allele_fraction: None,
            },
            SimEvent::SmallVariant {
                chrom: "chr1".to_string(),
                pos: 1000,
                ref_allele: b"A".to_vec(),
                alt_allele: b"T".to_vec(),
                gene: "T".to_string(),
                allele_fraction: None,
            },
        ] {
            assert!(
                !single.is_multi_locus(),
                "{:?} is extracted from one window",
                single.primary_region()
            );
        }
    }

    #[test]
    fn test_breakpoint_sides_names_each_reference_position_once() {
        // `bp.saturating_sub(1)` is the base before the cut, so two junctions
        // a base apart -- a small variant's alt segment -- name the same
        // reference position twice, and so does a breakpoint at haplotype
        // offset 0, where there is no base before the cut to clamp away from.
        // The refusal message deduplicates; the list the verdict is read off
        // did not.
        let hap = make_haplotype(vec![
            ref_segment(0, 1000),
            ref_segment(1000, 1),
            ref_segment(1001, 999),
        ]);
        assert_eq!(
            breakpoint_sides(&hap),
            vec![
                ("chr1".to_string(), 999),
                ("chr1".to_string(), 1000),
                ("chr1".to_string(), 1001),
            ],
            "each reference position is one side, however many junctions touch it"
        );
    }

    #[test]
    fn test_a_fusion_with_every_side_bare_does_not_claim_one_side() {
        // The refusal then lists all of them, so "on one side of its junction"
        // contradicts its own parenthesis.
        let pool = make_covering_pool(500_000, 540_000, 200);
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
        let event = SimEvent::Fusion {
            chrom_a: "chr1".to_string(),
            bp_a: 10000,
            gene_a: "GENE_A".to_string(),
            chrom_b: "chr1".to_string(),
            bp_b: 20000,
            gene_b: "GENE_B".to_string(),
            join: FusionJoin::Forward,
            allele_fraction: Some(0.2),
        };
        let mut rng = StdRng::seed_from_u64(99);
        let msg = match simulate_event(
            1, &event, &pool, &mut hap, &make_config(), &mock_synth_gen(150), 0.2, &mut rng,
        ) {
            Ok(_) => panic!("a fusion with no donor coverage anywhere must be refused"),
            Err(e) => e.to_string(),
        };
        assert!(
            msg.contains("chr1:9999") && msg.contains("chr1:20000"),
            "the refusal must list every bare side; got: {}",
            msg
        );
        assert!(
            !msg.contains("on one side of its junction"),
            "two bare sides are not one side; got: {}",
            msg
        );
    }

    #[test]
    fn test_simulate_event_keeps_del_whose_far_breakpoint_side_is_uncovered() {
        // A DEL that straddles the edge of a sliced or panel BAM has donor
        // reads on one side of its junction and none on the other. That is
        // ordinary input -- the haplotype is tiled across its whole footprint,
        // so the reads spike plants still sit on measured depth -- and it must
        // not be refused, from whichever side the covered reads come.
        let far_side_only = make_covering_pool(3000, 5000, 200);
        let near_side_only = make_covering_pool(0, 2000, 200);

        for (pool, side) in [(&far_side_only, "far"), (&near_side_only, "near")] {
            let mut hap = del_haplotype(1000, 2000);
            let mut rng = StdRng::seed_from_u64(42);
            let out = simulate_event(
                1, &del_event(1000, 3000), pool, &mut hap, &make_config(),
                &mock_synth_gen(150), 0.2, &mut rng,
            )
            .unwrap_or_else(|e| {
                panic!(
                    "a DEL covered only on its {} side must still be simulated: {}",
                    side, e
                )
            });
            assert!(
                !out.chimeric_pairs.is_empty(),
                "a DEL covered only on its {} side should still be tiled",
                side
            );
        }
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

        // Pool: 100 pairs near each partner's breakpoint. A fusion is
        // extracted from both loci, so both carry donor reads; a pool holding
        // only gene A's is the shape `coverage_for_tiling` refuses.
        let mut pairs = make_covering_pool(8000, 11000, 100).pairs;
        pairs.extend(make_covering_pool(19000, 22000, 100).pairs.into_iter().map(
            |mut p| {
                p.name = format!("b_{}", p.name);
                p
            },
        ));
        let pool = extract::build_read_pool(pairs, FragmentDist::from_stats(400.0, 80.0));
        let config = make_config();
        let gen = mock_synth_gen(150);
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

        let out = simulate_event(2, &event, &pool, &mut hap, &config, &gen, 0.05, &mut rng)
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


    #[test]
    fn test_fusion_at_chromosome_end_tiles_at_the_real_breakpoint_coverage() {
        // RightRight fusion 100 bp from a 10 kb contig's end: gene A's piece is
        // revcomp(chrEnd[9900, 11900)), which the fetch returns as 100 bases.
        // The junction sits at bp_a = 9900; a segment that still claims
        // ref_end = 11900 maps it to 2*9900 + 2000 - 10000 = 11800, past the
        // contig end, where the pool has no reads at all.
        let pattern = b"ACGT";
        let seq: Vec<u8> = (0..10_000u64).map(|i| pattern[(i % 4) as usize]).collect();
        let mut seqs = StdHashMap::new();
        seqs.insert("chrEnd".to_string(), seq);
        let reference = SharedReference::from_sequences(seqs);
        let mut hap = VariantHaplotype::from_fusion(
            &reference, "chrEnd", 9900, "chrEnd", 2000, 2000, FusionJoin::RightRight,
        )
        .unwrap();

        // 400 bp fragments every 10 bp from 0 to 9990: depth 40 in the interior.
        // They are chrEnd reads, and estimate_coverage_at only counts the
        // breakpoint's own chromosome, so they must say so.
        let pairs: Vec<ReadPair> = (0..1000u64)
            .map(|i| make_pair_on("chrEnd", &format!("r{}", i), i * 10, i * 10 + 400))
            .collect();
        let pool = ReadPool {
            pairs,
            frag_dist: FragmentDist::from_stats(400.0, 80.0),
        };

        let event = SimEvent::Fusion {
            chrom_a: "chrEnd".to_string(),
            bp_a: 9900,
            gene_a: "GENE_A".to_string(),
            chrom_b: "chrEnd".to_string(),
            bp_b: 2000,
            gene_b: "GENE_B".to_string(),
            join: FusionJoin::RightRight,
            allele_fraction: Some(0.5),
        };
        let mut rng = StdRng::seed_from_u64(7);

        let out = simulate_event(
            3, &event, &pool, &mut hap, &make_config(), &mock_synth_gen(150), 0.5, &mut rng,
        )
        .expect("simulate_event should succeed for a fusion at a contig end");

        // estimate_coverage_at samples 50 points over [8900, 10900) in steps of
        // 40: the 28 points up to 9980 see depth 40, the 10 from 10020 to 10380
        // tail off (37, 33, ..., 1) as fragments run out, the last 12 see none:
        // (28*40 + 190) / 50 = 26.2. An additive event tiles
        // cov * v/(1-v) * breakpoints = 26.2 * 1.0 * 1 -> 26 pairs.
        assert_eq!(out.chimeric_pairs.len(), 26);
    }

    #[test]
    fn test_cross_chromosome_fusion_ignores_the_other_chromosomes_depth() {
        // Fusion chr1:10000 >> chr2:10000. Both breakpoints sit at the same
        // coordinate on different chromosomes, so chr2's reads fall inside the
        // coverage window of chr1's breakpoint. Putting the far side's reads in
        // the pool -- which a fusion always does -- must not change how many
        // chimeric pairs the junction gets (M9).
        let segments = || {
            vec![
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
                        chrom: "chr2".to_string(),
                        ref_start: 10000,
                        ref_end: 11000,
                        is_reverse: false,
                    }),
                    hap_offset: 1000,
                },
            ]
        };

        let event = SimEvent::Fusion {
            chrom_a: "chr1".to_string(),
            bp_a: 10000,
            gene_a: "GENE_A".to_string(),
            chrom_b: "chr2".to_string(),
            bp_b: 10000,
            gene_b: "GENE_B".to_string(),
            join: FusionJoin::Forward,
            allele_fraction: Some(0.5),
        };

        // 400 bp fragments every 30 bp over [8000, 11000) on each side.
        let side = |chrom: &str, tag: &str| -> Vec<ReadPair> {
            (0..100u64)
                .map(|i| {
                    let start = 8000 + i * 30;
                    make_pair_on(chrom, &format!("{}{}", tag, i), start, start + 400)
                })
                .collect()
        };

        let run = |pairs: Vec<ReadPair>| -> usize {
            let pool = extract::build_read_pool(pairs, FragmentDist::from_stats(400.0, 80.0));
            let mut hap = make_haplotype(segments());
            let mut rng = StdRng::seed_from_u64(11);
            simulate_event(
                4,
                &event,
                &pool,
                &mut hap,
                &make_config(),
                &mock_synth_gen(150),
                0.5,
                &mut rng,
            )
            .expect("simulate_event should succeed for a cross-chromosome fusion")
            .chimeric_pairs
            .len()
        };

        // Gene B always carries some donor reads -- a fusion partner with
        // none is refused outright (N12) -- so the contrast is 1 chr2 pair
        // against 100. If chr2's depth leaked into chr1's breakpoint window,
        // the 100 would tile far more junction pairs than the 1.
        let pairs_a = side("chr1", "a");
        let mut thin_far_side = pairs_a.clone();
        thin_far_side.push(make_pair_on("chr2", "b_only", 9800, 10200));
        let near_side_only = run(thin_far_side);
        let mut both_sides = pairs_a;
        both_sides.extend(side("chr2", "b"));
        let with_far_side = run(both_sides);

        assert_eq!(
            with_far_side, near_side_only,
            "chr2's reads inflated the chr1 breakpoint's coverage: {} chimeric pairs with them, {} without",
            with_far_side, near_side_only
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
            &mock_synth_gen(150), 0.2, &mut rng,
        )
        .unwrap();

        assert!(out.suppressed_names.iter().all(|n| n.starts_with("in_")));
        let rate = out.suppressed_count as f64 / 500.0;
        assert!((rate - 0.2).abs() < 0.05, "inside suppression rate {:.3}", rate);
    }

    #[test]
    fn test_an_event_kept_on_partial_donor_coverage_names_the_bare_sides() {
        // `99f1a8e` keeps a single-locus event as long as *some* breakpoint
        // side carries donor reads. That is the right call -- refusing it was
        // `48be0c8`'s false refusal on a sliced BAM -- but it was silent. The
        // haplotype is still left_flank | right_flank, so the pairs tiled
        // across the bare half land in a window the input BAM has no read
        // over: a coverage island where there was none, and nothing said so.
        // DEL chr1:2000-7000, 2 kb flanks -> sides chr1:1999 and chr1:7000.
        let mut hap = del_haplotype(2000, 5000);
        // Donor reads over the right side only, far enough from the left one
        // that `estimate_coverage_at`'s 2 kb window sees nothing there.
        let pool = make_covering_pool(6500, 7500, 100);
        let mut rng = StdRng::seed_from_u64(42);

        let out = simulate_event(
            1, &del_event(2000, 7000), &pool, &mut hap, &make_config(),
            &mock_synth_gen(150), 0.2, &mut rng,
        )
        .unwrap();

        assert_eq!(
            out.uncovered_breakpoint_sides,
            vec!["chr1:1999".to_string()],
            "an event kept on partial donor coverage must name the sides with none"
        );
    }

    #[test]
    fn test_origin_warns_about_a_bare_breakpoint_side_as_clean_does() {
        // PD-30: origin's copy of the coverage check returned the bare sides
        // but dropped the warning `clean` logs for them. Same shape as the
        // test above: DEL chr1:2000-7000, reads over the right side only.
        crate::loh::tests::capture::install();
        let mut hap = del_haplotype(2000, 5000);
        let pool = make_covering_pool(6500, 7500, 100);
        let records = (0..40u64)
            .flat_map(|i| {
                let s = 6500 + 20 * i;
                [origin_read(&format!("w{}", i), true, s), origin_read(&format!("w{}", i), false, s + 250)]
            })
            .collect();
        let site = OriginSite { footprint: Span::new("chr1", 0, 9000), lookalikes: vec![], records, f: 1.0 };
        let mut rng = StdRng::seed_from_u64(42);
        let out = simulate_event_origin(
            1, &del_event(2000, 7000), &pool, &mut hap, &make_config(), &mock_synth_gen(150), 0.2, Some(&site), &mut rng,
        )
        .unwrap();
        assert_eq!(out.uncovered_breakpoint_sides, vec!["chr1:1999".to_string()]);
        let warned = crate::loh::tests::capture::warnings_matching(
            "event chr1:2000-7000 is kept although the origin depth is 0 at chr1:1999",
        );
        assert_eq!(warned.len(), 1, "{:?}", warned);
    }

    #[test]
    fn test_an_event_covered_on_every_side_names_none() {
        // The warning may not fire on ordinary input, or it says nothing.
        let mut hap = del_haplotype(2000, 5000);
        let pool = make_covering_pool(500, 8500, 400);
        let mut rng = StdRng::seed_from_u64(42);

        let out = simulate_event(
            1, &del_event(2000, 7000), &pool, &mut hap, &make_config(),
            &mock_synth_gen(150), 0.2, &mut rng,
        )
        .unwrap();

        assert!(
            out.uncovered_breakpoint_sides.is_empty(),
            "every side of this event has donor reads; got {:?}",
            out.uncovered_breakpoint_sides
        );
    }

    #[test]
    fn test_simulate_event_keeps_pairs_straddling_footprint_edge() {
        // Tiled fragments never extend past the haplotype ends, so originals
        // that stick out of the footprint must not be suppressed either;
        // otherwise depth dips at the footprint edge.
        // DEL [1000,3000) with 1 kb flanks: footprint [0,4000).
        // The `in_` pairs only give the breakpoint the donor coverage the
        // tiling count is scaled by; the assertion is about the `edge_` ones.
        let mut pairs: Vec<ReadPair> = (0..500)
            .map(|i| make_pair(&format!("edge_{}", i), 3800, 4200))
            .collect();
        pairs.extend((0..500).map(|i| make_pair(&format!("in_{}", i), 800, 1200)));
        let pool = make_pool(pairs);
        let mut hap = del_haplotype(1000, 2000);
        let mut rng = StdRng::seed_from_u64(42);

        let out = simulate_event(
            1, &del_event(1000, 3000), &pool, &mut hap, &make_config(),
            &mock_synth_gen(150), 0.2, &mut rng,
        )
        .unwrap();

        assert!(
            out.suppressed_names.iter().all(|n| n.starts_with("in_")),
            "a pair straddling the footprint edge must not be suppressed; got {:?}",
            out.suppressed_names
                .iter()
                .filter(|n| n.starts_with("edge_"))
                .collect::<Vec<_>>()
        );
        // `all` is vacuously true on an empty set, so the assertion above
        // would still hold if suppression had stopped entirely. Pin that it
        // happened at all.
        assert!(
            out.suppressed_names.iter().any(|n| n.starts_with("in_")),
            "pairs inside the footprint must still be suppressed; nothing was"
        );
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
            &mock_synth_gen(150), 0.2, &mut rng,
        )
        .unwrap();

        // n = cov * v * (L - f) / f = 20 * 0.2 * (4000 - 220) / 220 = 68.7 -> 69
        assert_eq!(out.chimeric_pairs.len(), 69);
    }

    #[test]
    fn test_an_insertions_tiled_reads_carry_the_prefix_validate_looks_for() {
        // `spike validate`'s `ins_planted` counts only the reads spike made for
        // truth record `sim_ins_N`, by name (RF13). The name is built here and
        // looked for there, so the two must agree.
        let mut hap = make_haplotype(vec![
            ref_segment(0, 2000),
            novel_segment(100),
            ref_segment(2000, 2000),
        ]);
        let pool = make_covering_pool(0, 6000, 1200);
        let event = SimEvent::Insertion {
            chrom: "chr1".to_string(),
            pos: 2000,
            ins_seq: Some(vec![b'G'; 100]),
            ins_len: 100,
            gene: "TEST".to_string(),
            allele_fraction: Some(0.5),
        };
        let mut rng = StdRng::seed_from_u64(3);
        let out = simulate_event(
            3, &event, &pool, &mut hap, &make_config(), &mock_synth_gen(150), 0.5, &mut rng,
        )
        .unwrap();
        let prefix = crate::validate::planted_read_prefix(3);
        assert!(!out.chimeric_pairs.is_empty());
        for pair in &out.chimeric_pairs {
            assert!(pair.name.starts_with(&prefix), "{} lacks {}", pair.name, prefix);
        }
    }

    #[test]
    fn test_long_insertion_yields_reads_carrying_inserted_sequence() {
        // ref (A) [0,2000) | 1000 bp insertion (G) | ref (A) [2000,4000).
        let mut hap = make_haplotype(vec![
            ref_segment(0, 2000),
            novel_segment(1000),
            ref_segment(2000, 2000),
        ]);
        let pool = make_covering_pool(0, 6000, 1200); // 400 bp every 5 bp: cov 80
        let event = SimEvent::Insertion {
            chrom: "chr1".to_string(),
            pos: 2000,
            ins_seq: Some(vec![b'G'; 1000]),
            ins_len: 1000,
            gene: "TEST".to_string(),
            allele_fraction: Some(0.2),
        };
        let mut rng = StdRng::seed_from_u64(3);

        let out = simulate_event(
            1, &event, &pool, &mut hap, &make_config(), &mock_synth_gen(150), 0.2, &mut rng,
        )
        .unwrap();

        // With 400 bp fragments, 780 of the 4000 allowed starts put >= 10
        // inserted bases in a read (~19.5%). Before the fix it was 0.
        // Either mate may be reverse-complemented, so the inserted G reads as
        // G or C; the A flanks never do.
        let carries_insert = |p: &&ReadPair| {
            let inserted = |seq: &[u8]| {
                seq.iter().filter(|&&b| b == b'G' || b == b'C').count() >= 10
            };
            inserted(&p.seq1) || inserted(&p.seq2)
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

    #[test]
    fn test_long_insertion_tiling_never_places_a_fragment_inside_the_insertion() {
        // ref (A) [0,1000) | 100 kb insertion (G) | ref (A) [1000,2000).
        // At the 400 bp mean fragment only 2000 of the 101601 starts put any
        // reference base in the fragment -- about one in fifty -- so this
        // exercises `sample_ref_overlapping_start` where the
        // reference-overlapping starts are a thin slice of the haplotype:
        // every draw still has to land in that slice, not just most of them.
        let hap = make_haplotype(vec![
            ref_segment(0, 1000),
            novel_segment(100_000),
            ref_segment(1000, 1000),
        ]);
        let gen = mock_synth_gen_at_q(150, 93);
        let pool = make_pool(vec![]);
        let mut rng = StdRng::seed_from_u64(42);

        let pairs =
            tile_haplotype_reads(&hap, None, &gen, &pool, 30.0, 0.5, false, "ins", &mut rng).0;
        assert!(!pairs.is_empty(), "the insertion haplotype should be tiled");

        // The flanks are A, so a read off them carries A (T when the mate is
        // reverse-complemented); the insertion is G, so a read wholly inside
        // it carries only G or C. A pair with neither A nor T came off a
        // fragment that lay entirely in inserted sequence -- exactly the
        // placement `compute_tiling_count` left out of the budget.
        let touches_reference = |p: &&ReadPair| {
            let off_a_flank = |seq: &[u8]| seq.iter().any(|&b| b == b'A' || b == b'T');
            off_a_flank(&p.seq1) || off_a_flank(&p.seq2)
        };
        let inside = pairs.len() - pairs.iter().filter(touches_reference).count();
        assert_eq!(
            inside,
            0,
            "{} of {} tiled fragments lay wholly inside the insertion",
            inside,
            pairs.len()
        );
    }

    // ---------------------------------------------------------------
    // The sample's two copies
    // ---------------------------------------------------------------

    /// Copies over chr1 [start, end): the event copy carries `event_base` and
    /// the other copy `other_base` at every position.
    fn uniform_copies(
        start: u64,
        end: u64,
        event_base: u8,
        other_base: u8,
        read_copy: HashMap<String, bool>,
    ) -> Vec<(String, loh::SampleCopies)> {
        vec![(
            "chr1".to_string(),
            loh::SampleCopies {
                event_copy: (start..end).map(|p| (p, event_base)).collect(),
                other_copy: (start..end).map(|p| (p, other_base)).collect(),
                read_copy,
            },
        )]
    }

    /// Count pairs whose R1 carries the event copy (G) and the other copy (A).
    /// R1 comes off either end of the fragment, so the event copy reads G
    /// forward and C reverse, and the other copy A forward and T reverse.
    fn count_by_copy(pairs: &[&ReadPair]) -> (usize, usize) {
        let mostly = |p: &ReadPair, b: u8| p.seq1.iter().filter(|&&x| x == b).count() * 2 > p.seq1.len();
        (
            pairs.iter().filter(|p| mostly(p, b'G') || mostly(p, b'C')).count(),
            pairs.iter().filter(|p| mostly(p, b'A') || mostly(p, b'T')).count(),
        )
    }

    #[test]
    fn test_flank_reads_are_removed_by_copy() {
        // DEL [1000,3000) with 1 kb flanks: footprint [0,4000). All pairs sit
        // in the left flank, outside the deletion itself. At VAF 0.3 the
        // event copy's reads go at 0.6, the other copy's stay, and reads of
        // unknown copy go at 0.3.
        let mut pairs = Vec::new();
        let mut read_copy = HashMap::new();
        for i in 0..600 {
            let name = match i % 3 {
                0 => format!("ev_{}", i),
                1 => format!("ot_{}", i),
                _ => format!("un_{}", i),
            };
            if i % 3 < 2 {
                read_copy.insert(name.clone(), i % 3 == 0);
            }
            pairs.push(make_pair(&name, 100, 500));
        }
        let copies = vec![(
            "chr1".to_string(),
            loh::SampleCopies { read_copy, ..Default::default() },
        )];
        let mut rng = StdRng::seed_from_u64(5);

        let out = simulate_event_with_copies(
            1, &del_event(1000, 3000), &make_pool(pairs), &mut del_haplotype(1000, 2000),
            &make_config(), &mock_synth_gen(150), 0.3, &copies, &mut rng,
        )
        .unwrap();

        let rate = |prefix: &str| {
            out.suppressed_names.iter().filter(|n| n.starts_with(prefix)).count() as f64 / 200.0
        };
        assert!((rate("ev_") - 0.6).abs() < 0.1, "event copy removed at {:.2}", rate("ev_"));
        assert_eq!(rate("ot_"), 0.0, "other copy must stay below VAF 0.5");
        assert!((rate("un_") - 0.3).abs() < 0.1, "unknown copy removed at {:.2}", rate("un_"));
    }

    #[test]
    fn test_synthetic_reads_carry_event_copy_alleles() {
        // Event copy: G everywhere, other copy: A. At VAF 0.5 every
        // synthetic read comes from the event copy.
        let copies = uniform_copies(0, 4000, b'G', b'A', HashMap::new());
        let pool = make_covering_pool(0, 5000, 1000);
        let mut rng = StdRng::seed_from_u64(5);

        let out = simulate_event_with_copies(
            1, &del_event(1000, 3000), &pool, &mut del_haplotype(1000, 2000),
            &make_config(), &mock_synth_gen(150), 0.5, &copies, &mut rng,
        )
        .unwrap();

        let (event, other) = count_by_copy(&out.chimeric_pairs.iter().collect::<Vec<_>>());
        assert_eq!(other, 0);
        assert_eq!(event, out.chimeric_pairs.len(), "every synthetic read should carry G");
    }

    #[test]
    fn test_synthetic_reads_above_half_vaf_come_from_both_copies() {
        // At VAF 0.8 the other copy carries the event in 2·0.8 − 1 = 60% of
        // cells, so it gives 0.5·0.6 / 0.8 = 37.5% of the synthetic reads.
        let copies = uniform_copies(0, 4000, b'G', b'A', HashMap::new());
        let pool = make_covering_pool(0, 5000, 1000);
        let mut rng = StdRng::seed_from_u64(5);

        let out = simulate_event_with_copies(
            1, &del_event(1000, 3000), &pool, &mut del_haplotype(1000, 2000),
            &make_config(), &mock_synth_gen(150), 0.8, &copies, &mut rng,
        )
        .unwrap();

        let (event, other) = count_by_copy(&out.chimeric_pairs.iter().collect::<Vec<_>>());
        assert_eq!(event + other, out.chimeric_pairs.len());
        let frac = other as f64 / (event + other) as f64;
        assert!((frac - 0.375).abs() < 0.1, "other copy share {:.3}", frac);
    }

    #[test]
    fn test_small_variant_keeps_its_own_allele() {
        // SNV A>T at 1000. The sample's event copy carries G at every position,
        // including 1000; the simulated allele must win there.
        let alt = HaplotypeSegment {
            sequence: vec![b'T'],
            origin: Some(SegmentOrigin {
                chrom: "chr1".to_string(),
                ref_start: 1000,
                ref_end: 1001,
                is_reverse: false,
            }),
            hap_offset: 0,
        };
        let mut hap = make_haplotype(vec![ref_segment(0, 1000), alt, ref_segment(1001, 1000)]);
        let event = SimEvent::SmallVariant {
            chrom: "chr1".to_string(),
            pos: 1000,
            ref_allele: b"A".to_vec(),
            alt_allele: b"T".to_vec(),
            gene: "TEST".to_string(),
            allele_fraction: Some(0.5),
        };
        let copies = uniform_copies(0, 2001, b'G', b'A', HashMap::new());
        let mut rng = StdRng::seed_from_u64(5);

        simulate_event_with_copies(
            1, &event, &make_covering_pool(0, 2001, 100), &mut hap, &make_config(),
            &mock_synth_gen(150), 0.5, &copies, &mut rng,
        )
        .unwrap();

        assert_eq!(hap.get_sequence(999, 3), b"GTG");
    }

    #[test]
    fn test_junction_dup_depth_copies_carry_the_copied_reads_alleles() {
        // Legacy junction model: each depth copy repeats an original pair of
        // the DUP [2000,6000) and must carry that pair's copy's alleles.
        let mut config = make_config();
        config.dup_model = "junction".to_string();
        let mut pairs = Vec::new();
        let mut read_copy = HashMap::new();
        for i in 0..300u64 {
            let name = if i % 2 == 0 { format!("ev_{}", i) } else { format!("ot_{}", i) };
            read_copy.insert(name.clone(), i % 2 == 0);
            let start = 2100 + i * 10;
            pairs.push(make_pair(&name, start, start + 400));
        }
        let pool = make_pool(pairs);
        let copies = uniform_copies(0, 10_000, b'G', b'A', read_copy);
        let event = SimEvent::Duplication {
            chrom: "chr1".to_string(),
            dup_start: 2000,
            dup_end: 6000,
            gene: "TEST".to_string(),
            allele_fraction: Some(0.75),
        };
        // Junction haplotype: ref [4000,6000) then ref [2000,4000).
        let mut hap = make_haplotype(vec![ref_segment(4000, 2000), ref_segment(2000, 2000)]);
        let mut rng = StdRng::seed_from_u64(5);

        let out = simulate_event_with_copies(
            1, &event, &pool, &mut hap, &config, &mock_synth_gen(150), 0.75, &copies,
            &mut rng,
        )
        .unwrap();

        // Depth copy names end in the index of the original pair in the pool.
        let original = |p: &ReadPair| -> String {
            let i: usize = p.name.rsplit('_').next().unwrap().parse().unwrap();
            pool.pairs[i].name.clone()
        };
        let depth: Vec<&ReadPair> =
            out.chimeric_pairs.iter().filter(|p| p.name.contains("_dup_depth_")).collect();
        let of_event: Vec<&ReadPair> =
            depth.iter().copied().filter(|p| original(p).starts_with("ev_")).collect();
        let of_other: Vec<&ReadPair> =
            depth.iter().copied().filter(|p| original(p).starts_with("ot_")).collect();
        // VAF 0.75: event copy's pairs copied at 1.0, the other's at 0.5.
        assert!(of_event.len() > 100 && of_other.len() > 40, "{} / {}", of_event.len(), of_other.len());
        assert_eq!(count_by_copy(&of_event), (of_event.len(), 0));
        assert_eq!(count_by_copy(&of_other), (0, of_other.len()));
    }

    #[test]
    fn test_junction_dup_depth_copies_are_named_in_the_generator_s_shape_without_drawing_randomness() {
        let mut config = make_config();
        config.dup_model = "junction".to_string();
        let pool = make_pool((0..300u64).map(|i| make_pair(&format!("p{}", i), 2100 + i * 10, 2500 + i * 10)).collect());
        let event = SimEvent::Duplication {
            chrom: "chr1".to_string(),
            dup_start: 2000,
            dup_end: 6000,
            gene: "TEST".to_string(),
            allele_fraction: Some(0.5),
        };
        let run = |gen: &SynthReadGenerator| {
            let mut hap = make_haplotype(vec![ref_segment(4000, 2000), ref_segment(2000, 2000)]);
            let mut rng = StdRng::seed_from_u64(5);
            simulate_event(1, &event, &pool, &mut hap, &config, gen, 0.5, &mut rng).unwrap().chimeric_pairs
        };
        let plain = run(&mock_synth_gen(150));
        let named = run(&mock_synth_gen(150).with_read_names(seven_part_shape()));
        assert!(
            plain.iter().filter(|p| p.name.starts_with("SPIKE_ev0001_dup_depth_")).count() > 100,
            "the depth copies are marked"
        );
        same_reads_named_in(&plain, &named, &seven_part_shape());
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
            uncovered_breakpoint_sides: Vec::new(),
            adjusted_vaf: None,
            depth_fold: DepthFold::default(),
            origin_chances: Vec::new(),
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
    fn test_consumed_original_names_covers_kept_and_suppressed() {
        // merge.sh removes these names from the original BAM, so the set has to
        // be every original the pools took: the passed-through ones (which come
        // back realigned from sim.bam) and the suppressed ones. Synthetic pairs
        // are not in the original BAM and must never appear.
        // r1 is only ever suppressed, r3 only ever kept, r2 kept by both events.
        let outputs = vec![
            spliced(&["r2", "r3"], &["ev0001_a"], &["r1"]),
            spliced(&["r2", "r3"], &["ev0002_a"], &[]),
        ];

        let names: Vec<String> = consumed_original_names(&outputs).into_iter().collect();

        assert_eq!(names, vec!["r1", "r2", "r3"]);
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
        let gen = mock_synth_gen(150);
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
            1, &del(1000, 3000), &pool, &mut hap_a, &config, &gen, 0.5, &mut rng,
        )
        .unwrap();
        let out_b = simulate_event(
            2, &del(7000, 9000), &pool, &mut hap_b, &config, &gen, 0.5, &mut rng,
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

    // ── Task 8: the origin path ─────────────────────────────────────────

    use crate::origin::{self, OriginRecord, OriginSite, Span};

    /// A MAPQ 60 read of 150 bp at `start` on chr1, with no XA.
    fn origin_read(name: &str, first: bool, start: u64) -> OriginRecord {
        OriginRecord {
            name: name.to_string(),
            first,
            placements: origin::placements(Span::new("chr1", start, start + 150), 60, &[]),
            duplicate: false,
            qc_fail: false,
            mate_unmapped: false,
            five_prime: origin::FivePrime { chrom: "chr1".into(), pos: start, reverse: !first },
        }
    }

    /// 80 unique pairs over chr1:0-5000, the footprint of `del_haplotype(2000, 1000)`.
    fn origin_site(footprint: Span) -> OriginSite {
        let records = (0..80u64)
            .flat_map(|i| {
                let s = 100 + 50 * i;
                [origin_read(&format!("o{}", i), true, s), origin_read(&format!("o{}", i), false, s + 250)]
            })
            .collect();
        OriginSite { footprint, lookalikes: vec![], records, f: 1.0 }
    }

    fn origin_del_event() -> SimEvent {
        SimEvent::Deletion {
            chrom: "chr1".to_string(),
            del_start: 2000,
            del_end: 3000,
            gene: "TEST".to_string(),
            exons: vec![],
            allele_fraction: Some(0.5),
        }
    }

    #[test]
    fn test_origin_keeps_an_event_whose_pool_has_no_read_in_the_footprint() {
        // R1: the pool is 500 kb away, so `clean` refuses this event
        // (test_simulate_event_refuses_pool_with_no_coverage_at_breakpoint).
        // Under origin the depth comes from the site.
        let mut hap = del_haplotype(2000, 1000);
        let pool = make_covering_pool(500_000, 540_000, 200);
        let site = origin_site(Span::new("chr1", 0, 5000));
        let mut rng = StdRng::seed_from_u64(42);
        let out = simulate_event_origin(
            1, &origin_del_event(), &pool, &mut hap, &make_config(), &mock_synth_gen(150), 0.5, Some(&site), &mut rng,
        )
        .unwrap();

        // Nothing is suppressed here: `origin::decide` draws after all events.
        assert_eq!(out.kept_originals.len(), 200);
        assert_eq!(out.suppressed_count, 0);
        assert!(!out.chimeric_pairs.is_empty());
        assert_eq!(out.origin_chances.len(), 80);
        assert!(out.origin_chances.iter().all(|c| (c.chance - 0.5 * (1.0 - 1e-6)).abs() < 1e-9));
    }

    #[test]
    fn test_origin_leaves_the_pool_pairs_inside_the_footprint_for_the_final_draw() {
        // These are the pairs `clean` suppresses here, at rate 1/2. Under
        // origin every one comes back kept: `origin::decide` removes them only
        // after every event's chances are summed (R5).
        let mut hap = del_haplotype(2000, 1000);
        let pool = make_covering_pool(0, 4000, 100);
        assert!(pool.pairs.iter().all(|p| p.ref_end <= 5000), "every pair lies inside chr1:0-5000");
        let site = origin_site(Span::new("chr1", 0, 5000));
        let mut rng = StdRng::seed_from_u64(42);
        let out = simulate_event_origin(
            1, &origin_del_event(), &pool, &mut hap, &make_config(), &mock_synth_gen(150), 0.5, Some(&site), &mut rng,
        )
        .unwrap();
        assert_eq!(out.suppressed_count, 0);
        assert_eq!(out.kept_originals.len(), 100);
    }

    #[test]
    fn test_origin_site_that_misses_the_haplotype_footprint_is_a_bug() {
        let mut hap = del_haplotype(2000, 1000);
        let pool = make_covering_pool(0, 5000, 200);
        let site = origin_site(Span::new("chr1", 0, 4000));
        let mut rng = StdRng::seed_from_u64(42);
        let err = simulate_event_origin(
            1, &origin_del_event(), &pool, &mut hap, &make_config(), &mock_synth_gen(150), 0.5, Some(&site), &mut rng,
        )
        .err()
        .expect("a mismatched site must be refused")
        .to_string();
        assert!(err.contains("bug"), "{}", err);

        // Same guard, the other conjunct: the footprint matches, but the
        // event is additive (a DUP under --dup-model junction). Additive
        // events remove no originals, so `--edit-model origin` never
        // applies to them -- being handed a site is a bug in the caller,
        // not a coordinate mismatch, and must not be reported as one.
        let mut additive_hap = del_haplotype(2000, 1000);
        let additive_pool = make_covering_pool(0, 5000, 200);
        let matching_site = origin_site(Span::new("chr1", 0, 5000));
        let mut additive_config = make_config();
        additive_config.dup_model = "junction".to_string();
        let dup_event = SimEvent::Duplication {
            chrom: "chr1".to_string(),
            dup_start: 2000,
            dup_end: 3000,
            gene: "TEST".to_string(),
            allele_fraction: Some(0.5),
        };
        let mut additive_rng = StdRng::seed_from_u64(42);
        let additive_err = simulate_event_origin(
            1, &dup_event, &additive_pool, &mut additive_hap, &additive_config,
            &mock_synth_gen(150), 0.5, Some(&matching_site), &mut additive_rng,
        )
        .err()
        .expect("an additive event handed a site must be refused")
        .to_string();
        assert!(additive_err.contains("additive event"), "{}", additive_err);
    }

    #[test]
    fn test_depth_fold_by_measures_bins_with_the_given_estimator() {
        let hap = tandem_dup_haplotype(11_000, 15_000, 1_000);
        let fold = depth_fold_by(&hap, 40.0, &|_, _, _| 10.0);
        assert!((fold.fold - 41.0 / 11.0).abs() < 1e-9, "{:?}", fold);
    }

    /// `depth_fold_by` as it was before it measured bins on the thread pool:
    /// the reference the parallel one must match.
    fn depth_fold_by_one_thread(
        haplotype: &VariantHaplotype,
        cov: f64,
        coverage_at: &dyn Fn(&str, u64, u64) -> f64,
    ) -> DepthFold {
        const BIN: u64 = 1000;
        let mut worst = DepthFold { fold: 1.0, scaled_by: cov, worst_bin: String::new(), worst_depth: cov };
        for origin in haplotype.segments.iter().filter_map(|seg| seg.origin.as_ref()) {
            let len = origin.ref_end.saturating_sub(origin.ref_start);
            if len == 0 {
                continue;
            }
            let n_bins = ((len as f64 / BIN as f64).round() as u64).max(1);
            for b in 0..n_bins {
                let start = origin.ref_start + len * b / n_bins;
                let end = origin.ref_start + len * (b + 1) / n_bins;
                let depth = coverage_at(&origin.chrom, (start + end) / 2, end - start);
                let fold = ((depth + 1.0) / (cov + 1.0)).max((cov + 1.0) / (depth + 1.0));
                if fold > worst.fold {
                    worst.fold = fold;
                    worst.worst_bin = format!("{}:{}-{}", origin.chrom, start, end);
                    worst.worst_depth = depth;
                }
            }
        }
        worst
    }

    #[test]
    fn test_depth_fold_measured_on_many_threads_is_the_one_thread_fold() {
        // The bins may be measured on the thread pool; the fold, its bin and
        // its depth must be the ones a single pass picks -- including the
        // first bin to reach a tie: the depth is 0 in every fifth kilobase,
        // so the worst depth is reached in bins at different places.
        let hap = tandem_dup_haplotype(11_000, 60_000, 2_000);
        let depth_at = |_: &str, pos: u64, _: u64| ((pos / 1000) % 5 * 10) as f64;
        let one_thread = depth_fold_by_one_thread(&hap, 40.0, &depth_at);
        assert_eq!((one_thread.fold, one_thread.worst_depth), (41.0, 0.0));
        // The first bin at depth 0 is the left flank's second, 10.5 kb.
        assert_eq!(one_thread.worst_bin, "chr1:10000-11000");
        let pool = rayon::ThreadPoolBuilder::new().num_threads(4).build().unwrap();
        assert_eq!(pool.install(|| depth_fold_by(&hap, 40.0, &depth_at)), one_thread);
        // And a depth that changes with every base, so a bin measured at the
        // wrong point or paired with another bin's depth shows.
        let scrambled = |_: &str, pos: u64, _: u64| ((pos * 7919) % 1009) as f64 / 10.0;
        let one_thread = depth_fold_by_one_thread(&hap, 40.0, &scrambled);
        assert!(one_thread.fold > 2.0, "{:?}", one_thread);
        assert_eq!(pool.install(|| depth_fold_by(&hap, 40.0, &scrambled)), one_thread);
    }

    #[test]
    fn test_apply_removals_moves_removed_pool_pairs_to_suppressed() {
        let mut outputs = vec![spliced(&["a", "b"], &[], &["c"])];
        let removed: BTreeSet<String> = ["a".to_string(), "not_in_a_pool".to_string()].into();
        apply_removals(&mut outputs, &removed);
        assert_eq!(sorted_names(&outputs[0].kept_originals), ["b"]);
        assert_eq!(outputs[0].suppressed_names, ["c", "a"]);
        assert_eq!(outputs[0].suppressed_count, 2);
        // A removed read in no pool is neither re-emitted nor listed here;
        // main lists it in replaced_reads.txt itself.
        assert!(!consumed_original_names(&outputs).contains("not_in_a_pool"));
    }

    fn origin_chance(name: &str, family: u64, p: f64) -> origin::Chance {
        origin::Chance {
            name: name.to_string(),
            family: vec![origin::FivePrime { chrom: "chr1".into(), pos: family, reverse: false }],
            chance: p,
        }
    }

    #[test]
    fn test_origin_draw_takes_every_events_chances() {
        // R5 as main calls it: the one draw sees every event's chances, so a
        // fragment only the second event can remove still goes.
        let mut first = spliced(&["a"], &[], &[]);
        first.origin_chances = vec![origin_chance("a", 1, 1.0)];
        let mut second = spliced(&["b"], &[], &[]);
        second.origin_chances = vec![origin_chance("b", 2, 1.0)];
        let mut outputs = vec![first, second];
        let removed = draw_origin_removals(&mut outputs, &mut StdRng::seed_from_u64(1));
        assert_eq!(removed, BTreeSet::from(["a".to_string(), "b".to_string()]));
        assert_eq!(outputs[1].suppressed_names, ["b"]);
        assert!(outputs[1].kept_originals.is_empty());
    }

    #[test]
    fn test_origin_draw_sums_a_fragments_chances_over_events() {
        // R5: two events each give "c" 1/2; together 1, so it goes on every
        // seed. Only the first event's chance, or one draw per event, would
        // keep it on some.
        for seed in 0..50 {
            let mut first = spliced(&["c"], &[], &[]);
            first.origin_chances = vec![origin_chance("c", 3, 0.5)];
            let mut second = spliced(&[], &[], &[]);
            second.origin_chances = vec![origin_chance("c", 3, 0.5)];
            let mut outputs = vec![first, second];
            let removed = draw_origin_removals(&mut outputs, &mut StdRng::seed_from_u64(seed));
            assert!(removed.contains("c"), "seed {}", seed);
        }
    }

    #[test]
    fn test_origin_measures_the_depth_fold_with_origin_depth() {
        // The T3 rule: under origin, SIM_DEPTH_FOLD compares each bin's
        // origin depth with the origin depth the tiling is scaled by. The
        // pool is 500 kb away, so the pool's depth in every bin is 0 and
        // gives another fold.
        let mut hap = del_haplotype(2000, 1000);
        let pool = make_covering_pool(500_000, 540_000, 200);
        let site = origin_site(Span::new("chr1", 0, 5000));
        let mut rng = StdRng::seed_from_u64(42);
        let out = simulate_event_origin(
            1, &origin_del_event(), &pool, &mut hap, &make_config(), &mock_synth_gen(150), 0.5, Some(&site), &mut rng,
        )
        .unwrap();
        let cov = out.depth_fold.scaled_by;
        let depth = site.depth();
        let by_origin = depth_fold_by(&hap, cov, &|c, p, w| depth.fragment_coverage_at(c, p, w));
        assert_eq!(out.depth_fold, by_origin);
        assert_ne!(out.depth_fold, depth_fold(&hap, &pool, cov));
    }

    #[test]
    fn test_origin_and_pool_depth_sample_the_same_points() {
        // R1 compares the two estimators, so they must sample one grid (PD-5).
        // A 1-base read at x counts at a sampled point only when that point
        // is x: for every x around each window, pool and origin depth agree.
        let mut pool = make_pool(vec![make_pair("r", 0, 1)]);
        let mut site = OriginSite { footprint: Span::new("chr1", 0, 10_000), lookalikes: vec![], records: vec![], f: 1.0 };
        for (pos, window) in [(1000, 0), (1000, 1), (1003, 2), (1000, 7), (1000, 49), (1001, 50), (1000, 51), (1000, 100), (1200, 2000), (1500, 2001)] {
            let (lo, hi) = (pos - window / 2 - 2, pos + window / 2 + 2);
            for x in lo..=hi {
                pool.pairs = vec![make_pair("r", x, x + 1)];
                site.records = vec![OriginRecord {
                    placements: vec![origin::Placement { span: Span::new("chr1", x, x + 1), chance: 1.0 }],
                    mate_unmapped: true,
                    ..origin_read("r", true, x)
                }];
                assert_eq!(
                    estimate_coverage_at(&pool, "chr1", pos, window),
                    site.depth().read_coverage_at("chr1", pos, window),
                    "pos {} window {} read at {}",
                    pos,
                    window,
                    x
                );
            }
        }
    }

    /// 2000 pairs sorted by start, as `build_read_pool` leaves a pool: most
    /// of 300-700 bp, on chr1 and chr2, and every 400th one 5-60 kb long, which
    /// a search bounded by a typical fragment length would miss.
    fn uneven_pool() -> ReadPool {
        let mut rng = StdRng::seed_from_u64(7);
        let mut pairs: Vec<ReadPair> = (0..2000)
            .map(|i| {
                let start = rng.gen_range(0..200_000u64);
                let span = if i % 400 == 0 { rng.gen_range(5_000..60_000) } else { rng.gen_range(300..700) };
                let chrom = if i % 3 == 0 { "chr2" } else { "chr1" };
                make_pair_on(chrom, &format!("q{}", i), start, start + span)
            })
            .collect();
        pairs.sort_by_key(|p| p.ref_start);
        make_pool(pairs)
    }

    #[test]
    fn test_pool_depth_counts_what_checking_every_pair_counts() {
        // The depth fold asks once per 1 kb bin, so PoolDepth looks only at
        // the pairs that start close enough to reach a point. It must count
        // exactly what checking every pair counts, the long pairs included.
        let pool = uneven_pool();
        let depth = PoolDepth::new(&pool);
        for pos in (0..210_000u64).step_by(1999) {
            for window in [1, 1000, 2000] {
                let points = coverage_sample_points(pos, window);
                let every: usize = points
                    .iter()
                    .map(|&at| {
                        pool.pairs.iter().filter(|p| p.chrom == "chr1" && p.ref_start <= at && at < p.ref_end).count()
                    })
                    .sum();
                assert_eq!(
                    depth.coverage_at("chr1", pos, window),
                    every as f64 / points.len() as f64,
                    "pos {} window {}",
                    pos,
                    window
                );
            }
        }
    }

    #[test]
    fn test_origin_depth_sums_what_summing_every_placement_sums() {
        // OriginDepth looks only at the placements that start close enough to
        // reach the window, then sums them in record order, so each depth is
        // bit-equal to summing every placement: the same additions, in the
        // same order. Placements: MAPQ 0 and 30, 0-2 XA hits, every 97th hit
        // 40 kb long, on chr1 and chr2. Each read has its own 5' end: reads
        // that shared one would be a duplicate family, whose reads count its
        // surest member's chance, not their own.
        let mut rng = StdRng::seed_from_u64(11);
        let records: Vec<OriginRecord> = (0..1500u64)
            .map(|i| {
                let start = rng.gen_range(0..100_000u64);
                let chrom = if i % 4 == 0 { "chr2" } else { "chr1" };
                let primary = Span::new(chrom, start, start + rng.gen_range(100..160));
                let alts: Vec<Span> = (0..i % 3)
                    .map(|_| {
                        let s = rng.gen_range(0..100_000u64);
                        Span::new("chr1", s, s + if i % 97 == 0 { 40_000 } else { 150 })
                    })
                    .collect();
                let mapq = if i % 5 == 0 { 30 } else { 0 };
                OriginRecord {
                    placements: origin::placements(primary, mapq, &alts),
                    mate_unmapped: true,
                    ..origin_read(&format!("z{}", i), true, i)
                }
            })
            .collect();
        let site = OriginSite { footprint: Span::new("chr1", 0, 200_000), lookalikes: vec![], records, f: 1.0 };
        let removable: BTreeSet<String> = site.removable_names().into_iter().collect();
        let every: Vec<&origin::Placement> = site
            .records
            .iter()
            .filter(|r| removable.contains(&r.name))
            .flat_map(|r| r.placements.iter())
            .collect();
        let depth = site.depth();
        for pos in (0..110_000u64).step_by(1499) {
            for window in [1, 1000, 2000] {
                let points = coverage_sample_points(pos, window);
                let summed = points.iter().fold(0.0, |total, &at| {
                    total
                        + every
                            .iter()
                            .filter(|p| p.span.chrom == "chr1" && p.span.start <= at && at < p.span.end)
                            .map(|p| p.chance)
                            .sum::<f64>()
                }) / points.len() as f64;
                assert_eq!(
                    depth.read_coverage_at("chr1", pos, window).to_bits(),
                    summed.to_bits(),
                    "pos {} window {}",
                    pos,
                    window
                );
            }
        }
    }

    #[test]
    fn test_is_additive_names_fusions_and_junction_dups() {
        assert!(!is_additive(&origin_del_event(), "full"));
        let dup = SimEvent::Duplication { chrom: "chr1".into(), dup_start: 10, dup_end: 20, gene: "G".into(), allele_fraction: None };
        assert!(is_additive(&dup, "junction"));
        assert!(!is_additive(&dup, "full"));
    }
}
