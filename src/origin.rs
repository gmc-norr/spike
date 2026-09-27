//! `--edit-model origin`: which original reads came from an event's edited
//! copy, wherever the aligner put them.
//!
//! The default (`clean`) edits only donor-pool pairs: both mates at
//! `--min-mapq` or above, a proper pair, no duplicate, secondary,
//! supplementary or QC-fail flag. Where the aligner cannot tell a spot from a
//! look-alike most reads fail that, so the event is barely planted (RF8).
//! `origin` gives every primary read a chance of having come from the
//! event's footprint, read from its MAPQ and its `XA` hits, and removes it
//! by that chance. The design is
//! `docs/superpowers/specs/2026-09-26-edit-model-origin-design.md`; comments
//! name its fixes R1-R5.

use std::collections::{BTreeMap, BTreeSet};

use anyhow::{bail, Result};
use noodles::sam::alignment::record::cigar::Op;

use crate::types::ReadPool;

/// A reference interval, 0-based half-open.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Span {
    pub chrom: String,
    pub start: u64,
    pub end: u64,
}

impl Span {
    pub fn new(chrom: &str, start: u64, end: u64) -> Self {
        Span {
            chrom: chrom.to_string(),
            start,
            end,
        }
    }

    /// Whether `other` lies wholly inside this span.
    pub fn holds(&self, other: &Span) -> bool {
        self.chrom == other.chrom && other.start >= self.start && other.end <= self.end
    }

    /// Whether `other` shares at least one base with this span.
    pub fn overlaps(&self, other: &Span) -> bool {
        self.chrom == other.chrom && other.start < self.end && self.start < other.end
    }
}

impl std::fmt::Display for Span {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}:{}-{}", self.chrom, self.start, self.end)
    }
}

/// One place a read may have come from, and the chance that it did.
#[derive(Debug, Clone, PartialEq)]
pub struct Placement {
    pub span: Span,
    pub chance: f64,
}

/// bwa-mem lists alternative hits in `XA` only when there are at most this
/// many (its `-h 5` default), so a MAPQ 0 read without `XA` has more.
const MAX_LISTED_HITS: usize = 5;

/// The chance a read came from where the aligner put it.
///
/// MAPQ 1 or more is the aligner's own estimate, `1 - 10^(-MAPQ/10)`. At
/// MAPQ 0 the hits are equally good: with `n_xa` listed alternatives the
/// primary is one of `1 + n_xa`; with none listed there are more than
/// [`MAX_LISTED_HITS`], so it is one of at least 6.
pub fn p_here(mapq: u8, n_xa: usize) -> f64 {
    if mapq > 0 {
        1.0 - 10f64.powf(-(mapq as f64) / 10.0)
    } else if n_xa > 0 {
        1.0 / (1 + n_xa) as f64
    } else {
        1.0 / (1 + MAX_LISTED_HITS) as f64
    }
}

/// The reference bases a CIGAR covers.
fn reference_length(ops: &[Op]) -> u64 {
    ops.iter()
        .filter(|op| op.kind().consumes_reference())
        .map(|op| op.len() as u64)
        .sum()
}

/// The reference bases a CIGAR string such as `100M1D51M` covers, or `None`
/// when it is malformed.
fn cigar_reference_length(cigar: &str) -> Option<u64> {
    let (mut len, mut n, mut digits) = (0u64, 0u64, false);
    for c in cigar.bytes() {
        if c.is_ascii_digit() {
            n = n * 10 + u64::from(c - b'0');
            digits = true;
            continue;
        }
        if !digits {
            return None;
        }
        match c {
            b'M' | b'D' | b'N' | b'=' | b'X' => len += n,
            b'I' | b'S' | b'H' | b'P' => {}
            _ => return None,
        }
        n = 0;
        digits = false;
    }
    (!digits).then_some(len)
}

/// The alternative hits in an `XA:Z` value: `chrom,±pos,CIGAR,NM;` each, with
/// `pos` 1-based and its sign the strand. Malformed entries are skipped.
pub fn parse_xa(xa: &str) -> Vec<Span> {
    xa.split(';')
        .filter(|entry| !entry.is_empty())
        .filter_map(|entry| {
            let mut fields = entry.split(',');
            let chrom = fields.next()?;
            let pos: i64 = fields.next()?.parse().ok()?;
            let len = cigar_reference_length(fields.next()?)?;
            let start = pos.unsigned_abs().checked_sub(1)?;
            Some(Span::new(chrom, start, start + len))
        })
        .collect()
}

/// A read's placements: its primary at [`p_here`], then each `XA` hit at an
/// equal share of the rest, `(1 - p_here) / k`.
pub fn placements(primary: Span, mapq: u8, alternatives: &[Span]) -> Vec<Placement> {
    let here = p_here(mapq, alternatives.len());
    let each = if alternatives.is_empty() {
        0.0
    } else {
        (1.0 - here) / alternatives.len() as f64
    };
    std::iter::once(Placement {
        span: primary,
        chance: here,
    })
    .chain(alternatives.iter().map(|span| Placement {
        span: span.clone(),
        chance: each,
    }))
    .collect()
}

/// The chance a read came from inside `footprint`: the sum over its
/// placements lying wholly inside it, the primary included (R2).
pub fn chance_within(placements: &[Placement], footprint: &Span) -> f64 {
    placements
        .iter()
        .filter(|p| footprint.holds(&p.span))
        .map(|p| p.chance)
        .sum()
}

use noodles::sam::alignment::record::cigar::op::Kind;

/// A read's unclipped 5' end: where its first sequenced base would align.
/// Duplicate marking keys a fragment on its two mates' 5' ends.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct FivePrime {
    pub chrom: String,
    pub pos: u64,
    pub reverse: bool,
}

/// Clipped bases at the start and at the end of a CIGAR.
fn clips(ops: &[Op]) -> (u64, u64) {
    let clipped = |op: &&Op| matches!(op.kind(), Kind::SoftClip | Kind::HardClip);
    let lead = ops.iter().take_while(clipped).map(|op| op.len() as u64).sum();
    let trail = ops.iter().rev().take_while(clipped).map(|op| op.len() as u64).sum();
    (lead, trail)
}

/// The unclipped 5' end of a read aligned at `start` (0-based) with `ops`.
pub fn five_prime(chrom: &str, start: u64, ops: &[Op], reverse: bool) -> FivePrime {
    let (lead, trail) = clips(ops);
    let pos = if reverse {
        start + reference_length(ops) + trail
    } else {
        start.saturating_sub(lead)
    };
    FivePrime {
        chrom: chrom.to_string(),
        pos,
        reverse,
    }
}

/// One primary record, as far as `origin` needs it.
#[derive(Debug, Clone)]
pub struct OriginRecord {
    pub name: String,
    /// Read 1 of its pair; the two records of one name differ here.
    pub first: bool,
    /// The primary placement first, then the `XA` hits.
    pub placements: Vec<Placement>,
    pub duplicate: bool,
    pub qc_fail: bool,
    pub mate_unmapped: bool,
    pub five_prime: FivePrime,
}

impl OriginRecord {
    /// Where the aligner put the read.
    pub fn primary(&self) -> &Placement {
        &self.placements[0]
    }
}

/// The records of one read name that `origin` read: one or both mates.
#[derive(Debug, Clone)]
pub struct Fragment<'a> {
    pub name: &'a str,
    pub mates: Vec<&'a OriginRecord>,
}

impl Fragment<'_> {
    /// The mate whose own placement is surest: the highest primary chance,
    /// ties to the higher chance of having come from `footprint`. One mate
    /// pinned uniquely somewhere pins the fragment there.
    fn surest(&self, footprint: &Span) -> &OriginRecord {
        self.mates
            .iter()
            .copied()
            .max_by(|a, b| {
                a.primary().chance.total_cmp(&b.primary().chance).then(
                    chance_within(&a.placements, footprint)
                        .total_cmp(&chance_within(&b.placements, footprint)),
                )
            })
            .expect("a fragment has at least one mate")
    }

    /// `p_origin`: the chance this fragment came from `footprint`, its
    /// surest mate's.
    pub fn chance(&self, footprint: &Span) -> f64 {
        chance_within(&self.surest(footprint).placements, footprint)
    }

    /// Whether spike's new reads can replace this fragment (R4): every mate
    /// has a placement wholly inside `footprint`, and a mate `origin` did not
    /// read is unmapped. The new fragments never reach past the footprint,
    /// so removing one that sticks out would leave a depth dip.
    pub fn removable(&self, footprint: &Span) -> bool {
        let every_mate_read = self.mates.len() == 2 || self.mates.iter().all(|m| m.mate_unmapped);
        every_mate_read
            && self
                .mates
                .iter()
                .all(|m| m.placements.iter().any(|p| footprint.holds(&p.span)))
    }

    /// Whether the aligner put a mate inside `footprint`, where the sample's
    /// phase call applies. Elsewhere the copy is unknown.
    pub fn at_spot(&self, footprint: &Span) -> bool {
        self.mates.iter().any(|m| footprint.holds(&m.primary().span))
    }

    /// The duplicate family: the mates' 5' ends, sorted. Duplicates of one
    /// molecule share it (R3).
    pub fn family(&self) -> Vec<FivePrime> {
        let mut ends: Vec<FivePrime> = self.mates.iter().map(|m| m.five_prime.clone()).collect();
        ends.sort();
        ends
    }
}

/// `XA` hits within this many bases of each other form one look-alike region.
const LOOKALIKE_GAP: u64 = 1000;
/// A region is a look-alike only when this many reads point into it.
const LOOKALIKE_MIN_READS: usize = 2;

/// The look-alike regions of `footprint`: the `XA` hits of the reads placed
/// over it that lie outside it, grouped when within [`LOOKALIKE_GAP`] of each
/// other and grown by `read_length` on each side. A region counts only when
/// at least [`LOOKALIKE_MIN_READS`] reads point into it.
pub fn lookalike_regions(spot: &[OriginRecord], footprint: &Span, read_length: u64) -> Vec<Span> {
    let mut hits: Vec<(&Span, (&str, bool))> = spot
        .iter()
        .flat_map(|r| {
            r.placements[1..]
                .iter()
                .map(move |p| (&p.span, (r.name.as_str(), r.first)))
        })
        .filter(|(span, _)| !footprint.overlaps(span))
        .collect();
    hits.sort();

    let mut regions = Vec::new();
    let mut i = 0;
    while i < hits.len() {
        let (first, _) = hits[i];
        let mut end = first.end;
        let mut reads: BTreeSet<(&str, bool)> = BTreeSet::new();
        while i < hits.len()
            && hits[i].0.chrom == first.chrom
            && hits[i].0.start <= end + LOOKALIKE_GAP
        {
            end = end.max(hits[i].0.end);
            reads.insert(hits[i].1);
            i += 1;
        }
        if reads.len() >= LOOKALIKE_MIN_READS {
            regions.push(Span::new(
                &first.chrom,
                first.start.saturating_sub(read_length),
                end + read_length,
            ));
        }
    }
    regions
}

/// Everything `origin` read for one event: the footprint the event's new
/// reads cover, the look-alike regions read besides it, every primary record
/// read in either, and the pool's fragment-to-read ratio.
#[derive(Debug, Clone)]
pub struct OriginSite {
    pub footprint: Span,
    pub lookalikes: Vec<Span>,
    pub records: Vec<OriginRecord>,
    /// See [`fragment_to_read_ratio`] (R1).
    pub f: f64,
}

impl OriginSite {
    /// The fragments with a placement inside the footprint, in name order.
    pub fn fragments(&self) -> Vec<Fragment<'_>> {
        let mut by_name: BTreeMap<&str, Vec<&OriginRecord>> = BTreeMap::new();
        for r in &self.records {
            by_name.entry(r.name.as_str()).or_default().push(r);
        }
        by_name
            .into_iter()
            .map(|(name, mates)| Fragment { name, mates })
            .filter(|f| {
                f.mates
                    .iter()
                    .any(|m| m.placements.iter().any(|p| self.footprint.holds(&p.span)))
            })
            .collect()
    }

    /// The names of the fragments spike may remove (R4).
    pub fn removable_names(&self) -> Vec<String> {
        self.fragments()
            .into_iter()
            .filter(|f| f.removable(&self.footprint))
            .map(|f| f.name.to_string())
            .collect()
    }

    /// Read depth that came from around `pos` on `chrom`. At up to 50 points
    /// over `window` (the points `simulate::estimate_coverage_at` samples),
    /// sum the chances of every placement covering the point, then average.
    /// Only fragments spike can remove count (R6): one kept by R4 stays in
    /// the BAM, so counting it would add new reads on top of it. Duplicate
    /// and QC-fail reads add nothing either, since spike's new reads are
    /// never flagged (R3).
    pub fn read_coverage_at(&self, chrom: &str, pos: u64, window: u64) -> f64 {
        let start = pos.saturating_sub(window / 2);
        let end = pos.saturating_add(window / 2);
        let range = end - start;
        let n = range.min(50).max(1);
        let step = if n > 1 { range / n } else { 1 };
        let removable: BTreeSet<String> = self.removable_names().into_iter().collect();
        let placed: Vec<&Placement> = self
            .records
            .iter()
            .filter(|r| !r.duplicate && !r.qc_fail && removable.contains(&r.name))
            .flat_map(|r| r.placements.iter())
            .filter(|p| p.span.chrom == chrom && p.span.start < end && start < p.span.end)
            .collect();
        let total: f64 = (0..n)
            .map(|i| {
                let at = start + i * step;
                placed
                    .iter()
                    .filter(|p| p.span.start <= at && at < p.span.end)
                    .map(|p| p.chance)
                    .sum::<f64>()
            })
            .sum();
        total / n as f64
    }

    /// [`read_coverage_at`](Self::read_coverage_at) in the fragment units the
    /// tiling count is in (R1).
    pub fn fragment_coverage_at(&self, chrom: &str, pos: u64, window: u64) -> f64 {
        self.read_coverage_at(chrom, pos, window) * self.f
    }
}

/// `f`: the pool's summed fragment spans over its summed read lengths (R1).
///
/// It depends on the library's fragment and read lengths, not on the spot,
/// so it is taken over the whole pool: an event inside a perfect twin has no
/// pool read in its footprint. Pool pairs keep no CIGAR, so read lengths are
/// whole reads. A soft-clipped read covers fewer bases than its length, so
/// where clips are common `f` runs a little low.
pub fn fragment_to_read_ratio(pool: &ReadPool) -> Result<f64> {
    let spans: u64 = pool.pairs.iter().map(|p| p.ref_end.saturating_sub(p.ref_start)).sum();
    let bases: u64 = pool.pairs.iter().map(|p| (p.seq1.len() + p.seq2.len()) as u64).sum();
    if bases == 0 {
        bail!(
            "the donor pool's {} read pair(s) hold no bases, so --edit-model origin cannot \
             convert read depth into fragment depth",
            pool.pairs.len()
        );
    }
    Ok(spans as f64 / bases as f64)
}

#[cfg(test)]
mod tests {
    use super::*;
    use noodles::sam::alignment::record::cigar::op::Kind;
    use crate::types::{ReadPair, ReadPool};

    fn close(a: f64, b: f64) -> bool {
        (a - b).abs() < 1e-9
    }

    /// A record placed at `start..start+100` with `alternatives` as its XA.
    pub(super) fn record(name: &str, first: bool, start: u64, mapq: u8, alternatives: &[Span]) -> OriginRecord {
        OriginRecord {
            name: name.to_string(),
            first,
            placements: placements(Span::new("chr1", start, start + 100), mapq, alternatives),
            duplicate: false,
            qc_fail: false,
            mate_unmapped: false,
            five_prime: FivePrime { chrom: "chr1".to_string(), pos: start, reverse: !first },
        }
    }

    fn fp() -> Span {
        Span::new("chr1", 0, 1000)
    }

    #[test]
    fn test_five_prime_steps_back_over_a_leading_clip_or_past_a_trailing_one() {
        let fwd = [Op::new(Kind::SoftClip, 5), Op::new(Kind::Match, 95)];
        assert_eq!(five_prime("chr1", 100, &fwd, false).pos, 95);
        let rev = [Op::new(Kind::Match, 95), Op::new(Kind::SoftClip, 5)];
        assert_eq!(five_prime("chr1", 300, &rev, true).pos, 400);
    }

    #[test]
    fn test_a_pair_takes_its_surest_mates_chance() {
        let a = record("p", true, 100, 60, &[]);
        let b = record("p", false, 300, 0, &[Span::new("chr1", 5000, 5100)]);
        let f = Fragment { name: "p", mates: vec![&a, &b] };
        assert!(close(f.chance(&fp()), 1.0 - 1e-6));
    }

    #[test]
    fn test_a_tie_on_the_primary_goes_to_the_higher_footprint_chance() {
        // Both MAPQ 0 with one hit (1/2 each). a's hit is outside (1/2 in the
        // footprint); b's is inside too (1).
        let a = record("p", true, 100, 0, &[Span::new("chr1", 5000, 5100)]);
        let b = record("p", false, 300, 0, &[Span::new("chr1", 600, 700)]);
        let f = Fragment { name: "p", mates: vec![&a, &b] };
        assert!(close(f.chance(&fp()), 1.0));
    }

    #[test]
    fn test_a_pair_is_removable_only_when_every_mate_could_come_from_the_footprint() {
        // R4.
        let a = record("p", true, 100, 60, &[]);
        let inside = record("p", false, 300, 60, &[]);
        let outside = record("p", false, 1200, 60, &[]);
        let outside_with_hit_inside = record("p", false, 1200, 0, &[Span::new("chr1", 700, 800)]);
        assert!(Fragment { name: "p", mates: vec![&a, &inside] }.removable(&fp()));
        assert!(!Fragment { name: "p", mates: vec![&a, &outside] }.removable(&fp()));
        assert!(Fragment { name: "p", mates: vec![&a, &outside_with_hit_inside] }.removable(&fp()));
    }

    #[test]
    fn test_an_unseen_mate_blocks_removal_unless_it_is_unmapped() {
        let alone = record("p", true, 100, 60, &[]);
        assert!(!Fragment { name: "p", mates: vec![&alone] }.removable(&fp()));
        let orphan = OriginRecord { mate_unmapped: true, ..record("p", true, 100, 60, &[]) };
        assert!(Fragment { name: "p", mates: vec![&orphan] }.removable(&fp()));
    }

    #[test]
    fn test_at_spot_means_a_mate_the_aligner_put_inside_the_footprint() {
        let here = record("p", true, 100, 0, &[Span::new("chr1", 5000, 5100)]);
        let there = record("q", true, 5000, 0, &[Span::new("chr1", 100, 200)]);
        assert!(Fragment { name: "p", mates: vec![&here] }.at_spot(&fp()));
        assert!(!Fragment { name: "q", mates: vec![&there] }.at_spot(&fp()));
    }

    #[test]
    fn test_duplicates_share_a_family_and_other_fragments_do_not() {
        let (a1, a2) = (record("a", true, 100, 60, &[]), record("a", false, 300, 60, &[]));
        let (d1, d2) = (record("d", true, 100, 60, &[]), record("d", false, 300, 60, &[]));
        let (o1, o2) = (record("o", true, 110, 60, &[]), record("o", false, 300, 60, &[]));
        let fam = |x: &OriginRecord, y: &OriginRecord| Fragment { name: "x", mates: vec![x, y] }.family();
        assert_eq!(fam(&a1, &a2), fam(&d2, &d1));
        assert_ne!(fam(&a1, &a2), fam(&o1, &o2));
    }

    #[test]
    fn test_p_here_follows_mapq_and_the_xa_count() {
        assert!(close(p_here(60, 0), 1.0 - 1e-6));
        assert!(close(p_here(3, 0), 1.0 - 10f64.powf(-0.3)));
        assert!(close(p_here(0, 1), 0.5));
        // No XA at MAPQ 0: bwa-mem lists at most 5, so there are more.
        assert!(close(p_here(0, 0), 1.0 / 6.0));
    }

    #[test]
    fn test_parse_xa_reads_chrom_start_strand_and_cigar_length() {
        let hits = parse_xa("cluster,+2002,151M,0;cluster,-1002,100M1D51M,1;");
        assert_eq!(
            hits,
            vec![Span::new("cluster", 2001, 2152), Span::new("cluster", 1001, 1153)]
        );
    }

    #[test]
    fn test_parse_xa_skips_a_malformed_entry() {
        assert_eq!(parse_xa("chr1,+10,5Q,0;chr1,+20,5M,0;"), vec![Span::new("chr1", 19, 24)]);
    }

    #[test]
    fn test_a_read_on_a_remote_copy_with_three_hits_in_the_footprint_gets_three_quarters() {
        // R2: bwa-mem put `cluster_1` on the remote copy at MAPQ 0 with all
        // three XA hits in one cluster (the review's probe).
        let alternatives = parse_xa("cluster,+2002,151M,0;cluster,+1002,151M,0;cluster,+1502,151M,0;");
        let placed = placements(Span::new("remote", 1001, 1152), 0, &alternatives);
        assert!(close(chance_within(&placed, &Span::new("cluster", 1000, 2400)), 0.75));
    }

    #[test]
    fn test_a_mapq0_read_whose_one_hit_is_also_in_the_footprint_gets_one() {
        let placed = placements(Span::new("chr1", 100, 250), 0, &[Span::new("chr1", 600, 750)]);
        assert!(close(chance_within(&placed, &Span::new("chr1", 0, 1000)), 1.0));
    }

    #[test]
    fn test_a_placement_sticking_out_of_the_footprint_does_not_count() {
        let placed = placements(Span::new("chr1", 950, 1100), 60, &[]);
        assert!(close(chance_within(&placed, &Span::new("chr1", 0, 1000)), 0.0));
    }

    #[test]
    fn test_hits_within_a_kilobase_form_one_region_grown_by_the_read_length() {
        let spot = [
            record("b", true, 200, 0, &[Span::new("chr1", 40_000, 40_100)]),
            record("b", false, 400, 0, &[Span::new("chr1", 40_900, 41_000)]),
        ];
        assert_eq!(
            lookalike_regions(&spot, &fp(), 150),
            vec![Span::new("chr1", 39_850, 41_150)]
        );
    }

    #[test]
    fn test_a_region_one_read_points_into_is_not_a_lookalike() {
        let spot = [
            record("b", true, 200, 0, &[Span::new("chr1", 40_000, 40_100)]),
            record("c", true, 300, 0, &[Span::new("chr1", 50_000, 50_100)]),
        ];
        assert!(lookalike_regions(&spot, &fp(), 150).is_empty());
    }

    #[test]
    fn test_hits_inside_the_footprint_are_not_lookalikes() {
        let spot = [
            record("b", true, 200, 0, &[Span::new("chr1", 600, 700)]),
            record("b", false, 400, 0, &[Span::new("chr1", 650, 750)]),
        ];
        assert!(lookalike_regions(&spot, &fp(), 150).is_empty());
    }

    #[test]
    fn test_hits_more_than_a_kilobase_apart_are_two_regions() {
        let spot = [
            record("b", true, 200, 0, &[Span::new("chr1", 40_000, 40_100), Span::new("chr1", 60_000, 60_100)]),
            record("c", true, 300, 0, &[Span::new("chr1", 40_050, 40_150), Span::new("chr1", 60_050, 60_150)]),
        ];
        assert_eq!(
            lookalike_regions(&spot, &fp(), 100),
            vec![Span::new("chr1", 39_900, 40_250), Span::new("chr1", 59_900, 60_250)]
        );
    }

    /// Twin spots L = chr1:2000-2150 and P = chr1:12000-12150: 20 reads at
    /// each, all MAPQ 0 with their one XA hit on the other. Each read's mate
    /// is unmapped, so a fragment is one record, and each has its own 5' end
    /// so no two share a duplicate family.
    pub(super) fn twin_site() -> OriginSite {
        let (l, p) = (Span::new("chr1", 2000, 2150), Span::new("chr1", 12_000, 12_150));
        let at = |name: String, end5: u64, here: &Span, there: &Span| OriginRecord {
            placements: placements(here.clone(), 0, std::slice::from_ref(there)),
            mate_unmapped: true,
            ..record(&name, true, end5, 0, &[])
        };
        let records = (0..20u64)
            .map(|i| at(format!("l{}", i), i, &l, &p))
            .chain((0..20u64).map(|i| at(format!("p{}", i), 100 + i, &p, &l)))
            .collect();
        OriginSite { footprint: Span::new("chr1", 0, 5000), lookalikes: vec![Span::new("chr1", 11_850, 12_300)], records, f: 1.5 }
    }

    #[test]
    fn test_origin_depth_at_a_twin_is_the_true_read_depth() {
        // 20 reads at L count 1/2 each and so do their 20 twins at P.
        let site = twin_site();
        assert!(close(site.read_coverage_at("chr1", 2075, 100), 20.0));
        assert!(close(site.fragment_coverage_at("chr1", 2075, 100), 30.0));
    }

    #[test]
    fn test_duplicate_and_qc_fail_reads_add_no_depth() {
        // R3.
        let mut site = twin_site();
        let copy = site.records[0].clone();
        site.records.push(OriginRecord { name: "dup".into(), duplicate: true, ..copy.clone() });
        site.records.push(OriginRecord { name: "qc".into(), qc_fail: true, ..copy });
        assert!(close(site.read_coverage_at("chr1", 2075, 100), 20.0));
    }

    #[test]
    fn test_reads_spike_cannot_remove_add_no_depth() {
        // R6: a read at L whose mapped mate spike never read is kept (R4),
        // so its depth must not pay for new reads on top of it.
        let mut site = twin_site();
        let copy = site.records[0].clone();
        site.records.push(OriginRecord { name: "stray".into(), mate_unmapped: false, ..copy });
        assert!(!site.removable_names().contains(&"stray".to_string()));
        assert!(close(site.read_coverage_at("chr1", 2075, 100), 20.0));
    }

    #[test]
    fn test_fragments_are_the_names_with_a_placement_in_the_footprint() {
        let mut site = twin_site();
        site.records.push(record("unique_at_p", true, 12_500, 60, &[]));
        let names: Vec<&str> = site.fragments().iter().map(|f| f.name).collect();
        assert_eq!(names.len(), 40);
        assert!(!names.contains(&"unique_at_p"));
    }

    fn pair(name: &str, start: u64, span: u64, bases: usize) -> ReadPair {
        ReadPair {
            name: name.to_string(),
            seq1: vec![b'A'; bases],
            qual1: vec![b'I'; bases],
            seq2: vec![b'A'; bases],
            qual2: vec![b'I'; bases],
            ref_start: start,
            ref_end: start + span,
            insert_size: span as i64,
            chrom: "chr1".to_string(),
        }
    }

    pub(super) fn pool(pairs: Vec<ReadPair>) -> ReadPool {
        ReadPool { pairs, frag_dist: crate::stats::FragmentDist::from_stats(400.0, 80.0) }
    }

    #[test]
    fn test_f_is_fragment_span_over_read_bases_across_the_whole_pool() {
        // R1: 400 bp fragments of two 150 bp reads, far from any footprint.
        let p = pool((0..30).map(|i| pair(&format!("r{}", i), 500_000 + 10 * i, 400, 150)).collect());
        assert!(close(fragment_to_read_ratio(&p).unwrap(), 400.0 / 300.0));
    }

    #[test]
    fn test_f_refuses_a_pool_without_bases() {
        let p = pool(vec![pair("r", 0, 400, 0)]);
        assert!(fragment_to_read_ratio(&p).is_err());
    }
}
