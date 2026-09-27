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
//! name its fixes R1-R8.

use std::collections::{BTreeMap, BTreeSet};
use std::collections::HashMap;

use anyhow::{bail, Context, Result};
use noodles::sam::alignment::record::cigar::Op;
use rand::rngs::StdRng;
use rand::Rng;

use crate::synth::copy_rate;
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
        let n = range.clamp(1, 50);
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

/// One event's chance of removing one fragment.
#[derive(Debug, Clone, PartialEq)]
pub struct Chance {
    pub name: String,
    pub family: Vec<FivePrime>,
    pub chance: f64,
}

impl OriginSite {
    /// This event's chance of removing each fragment it can remove:
    /// `p_origin x copy_rate(copy, vaf)`. At the spot `copy` is the sample's
    /// phase call for the fragment's duplicate family (R7): phasing skips
    /// duplicates (`src/loh.rs:628`), so a family takes the call of whichever
    /// member has one, and members that disagree get none. A fragment the
    /// aligner put only at a look-alike has no call, so its rate is `vaf`.
    pub fn removal_chances(&self, read_copy: &HashMap<String, bool>, vaf: f64) -> Vec<Chance> {
        let fragments: Vec<Fragment<'_>> = self
            .fragments()
            .into_iter()
            .filter(|f| f.removable(&self.footprint))
            .collect();
        let mut family_copy: BTreeMap<Vec<FivePrime>, Option<bool>> = BTreeMap::new();
        for f in &fragments {
            if let Some(&copy) = read_copy.get(f.name) {
                family_copy
                    .entry(f.family())
                    .and_modify(|call| {
                        if *call != Some(copy) {
                            *call = None;
                        }
                    })
                    .or_insert(Some(copy));
            }
        }
        fragments
            .into_iter()
            .filter_map(|f| {
                let family = f.family();
                let copy = if f.at_spot(&self.footprint) {
                    family_copy.get(&family).copied().flatten()
                } else {
                    None
                };
                let chance = f.chance(&self.footprint) * copy_rate(copy, vaf);
                (chance > 0.0).then(|| Chance {
                    name: f.name.to_string(),
                    family,
                    chance,
                })
            })
            .collect()
    }
}

/// Which fragments to remove, over every event's chances at once (R5).
///
/// A fragment came from one copy, so "it came from event 1's edited copy"
/// and "it came from event 2's" cannot both be true: its chances add,
/// capped at 1. One molecule has one origin, so a duplicate family is
/// removed or kept whole (R3, R7). There is one draw per family, in family
/// order, against its highest member's total.
pub fn decide(chances: &[Chance], rng: &mut StdRng) -> BTreeSet<String> {
    let mut totals: BTreeMap<&str, (&[FivePrime], f64)> = BTreeMap::new();
    for c in chances {
        totals
            .entry(c.name.as_str())
            .or_insert((c.family.as_slice(), 0.0))
            .1 += c.chance;
    }
    let mut families: BTreeMap<&[FivePrime], f64> = BTreeMap::new();
    for (family, total) in totals.values() {
        let highest = families.entry(*family).or_insert(0.0);
        *highest = highest.max(total.min(1.0));
    }
    let removed_families: BTreeSet<&[FivePrime]> = families
        .into_iter()
        .filter(|(_, highest)| rng.gen::<f64>() < *highest)
        .map(|(family, _)| family)
        .collect();
    totals
        .into_iter()
        .filter(|(_, (family, _))| removed_families.contains(family))
        .map(|(name, _)| name.to_string())
        .collect()
}

/// One [`OriginRecord`] from an alignment record, or `None` for a record
/// `origin` does not judge on its own. That is an unmapped, secondary or
/// supplementary record (those share their primary's name and go with it),
/// or one without a name. MAPQ 255, "unavailable", counts as 0.
fn origin_record(
    header: &noodles::sam::Header,
    buf: &noodles::sam::alignment::RecordBuf,
) -> Option<OriginRecord> {
    use noodles::sam::alignment::record::data::field::Tag;
    use noodles::sam::alignment::record_buf::data::field::Value;

    let flags = buf.flags();
    if flags.is_unmapped() || flags.is_secondary() || flags.is_supplementary() {
        return None;
    }
    let name = buf.name()?.to_string();
    let (chrom, _) = header
        .reference_sequences()
        .get_index(buf.reference_sequence_id()?)?;
    let chrom = chrom.to_string();
    let start = usize::from(buf.alignment_start()?) as u64 - 1;
    let ops: &[Op] = buf.cigar().as_ref();
    let xa = match buf.data().get(&Tag::new(b'X', b'A')) {
        Some(Value::String(s)) => Some(s.to_string()),
        _ => None,
    };
    let alternatives = xa.as_deref().map(parse_xa).unwrap_or_default();
    let mapq = buf.mapping_quality().map(|q| q.get()).unwrap_or(0);
    Some(OriginRecord {
        name,
        first: flags.is_first_segment(),
        placements: placements(
            Span::new(&chrom, start, start + reference_length(ops)),
            mapq,
            &alternatives,
        ),
        duplicate: flags.is_duplicate(),
        qc_fail: flags.is_qc_fail(),
        mate_unmapped: flags.is_mate_unmapped(),
        five_prime: five_prime(&chrom, start, ops, flags.is_reverse_complemented()),
    })
}

/// Every primary record overlapping `span`, from a BAM or a CRAM.
fn scan(bam_path: &str, ref_path: &str, span: &Span) -> Result<Vec<OriginRecord>> {
    let region = noodles::core::Region::new(
        span.chrom.as_str(),
        crate::extract::safe_noodles_position(span.start + 1)
            ..=crate::extract::safe_noodles_position(span.end),
    );
    let mut out = Vec::new();
    if crate::extract::is_cram(bam_path) {
        let repository = crate::extract::build_fasta_repository(ref_path)?;
        let (mut reader, header) =
            crate::extract::open_cram_reader_for_region(bam_path, &repository, &region)
                .context("failed to open CRAM for the origin scan")?;
        let queried = header
            .reference_sequences()
            .get_index_of(span.chrom.as_bytes());
        for result in reader.query(&header, &region)? {
            let buf = result?.try_into_alignment_record(&header)?;
            // A container holding several contigs is decoded whole (L2, N4).
            if !crate::extract::record_is_on_queried_reference(&buf, queried) {
                continue;
            }
            out.extend(origin_record(&header, &buf));
        }
    } else {
        let mut reader = noodles::bam::io::indexed_reader::Builder::default()
            .build_from_path(bam_path)
            .with_context(|| format!("failed to open BAM for the origin scan: {}", bam_path))?;
        let header = reader.read_header()?;
        for result in reader.query(&header, &region)? {
            let record = result?;
            let buf =
                noodles::sam::alignment::RecordBuf::try_from_alignment_record(&header, &record)?;
            out.extend(origin_record(&header, &buf));
        }
    }
    Ok(out)
}

/// How many records, from the start of the file, [`require_xa`] reads
/// looking for an `XA` tag. The 35x HG002 BAM's first `XA` is on record 64.
pub const XA_PROBE_RECORDS: usize = 100_000;

/// The number of records read up to and including the first that carries
/// `XA`, or `None` when none of the first `limit` does.
pub fn first_xa_record(bam_path: &str, ref_path: &str, limit: usize) -> Result<Option<usize>> {
    use noodles::sam::alignment::record::data::field::Tag;
    let xa = Tag::new(b'X', b'A');
    if crate::extract::is_cram(bam_path) {
        let repository = crate::extract::build_fasta_repository(ref_path)?;
        let mut reader = noodles::cram::io::reader::Builder::default()
            .set_reference_sequence_repository(repository)
            .build_from_path(bam_path)
            .with_context(|| format!("failed to open CRAM: {}", bam_path))?;
        let header = reader.read_header()?;
        for (i, result) in reader.records(&header).take(limit).enumerate() {
            let buf = result?.try_into_alignment_record(&header)?;
            if buf.data().get(&xa).is_some() {
                return Ok(Some(i + 1));
            }
        }
    } else {
        let mut reader = noodles::bam::io::reader::Builder
            .build_from_path(bam_path)
            .with_context(|| format!("failed to open BAM: {}", bam_path))?;
        reader.read_header()?;
        for (i, result) in reader.records().take(limit).enumerate() {
            if result?.data().get(&xa).is_some() {
                return Ok(Some(i + 1));
            }
        }
    }
    Ok(None)
}

/// Stop unless the file keeps the aligner's `XA` tags (R8); otherwise every
/// MAPQ 0 read would silently get 1/6 and no look-alike would be found.
///
/// The check is on the file, not the spot: a spot whose MAPQ 0 reads carry
/// no `XA` is normal, since under bwa-mem's `-h 5` rule their hits number
/// more than 5. Returns the record the first `XA` is on.
pub fn require_xa(bam_path: &str, ref_path: &str) -> Result<usize> {
    match first_xa_record(bam_path, ref_path, XA_PROBE_RECORDS)? {
        Some(n) => Ok(n),
        None => bail!(
            "--edit-model origin needs the aligner's XA tags (its alternative hits), but none \
             of the first {} records of {} carries one. bwa-mem and bwa-mem2 write XA by \
             default, and a later step can strip it. Re-align with one of them, or use \
             --edit-model clean.",
            XA_PROBE_RECORDS,
            bam_path,
        ),
    }
}

/// Read everything `origin` needs for one event: the footprint, its
/// look-alike regions and `f` (see [`OriginSite`]). Whether the file keeps
/// `XA` at all is [`require_xa`]'s check, made once per run.
pub fn gather(
    bam_path: &str,
    ref_path: &str,
    footprint: &Span,
    read_length: usize,
    pool: &ReadPool,
) -> Result<OriginSite> {
    let spot = scan(bam_path, ref_path, footprint)?;
    let lookalikes = lookalike_regions(&spot, footprint, read_length as u64);
    let mut seen: BTreeSet<(String, bool)> =
        spot.iter().map(|r| (r.name.clone(), r.first)).collect();
    let mut records = spot;
    for region in &lookalikes {
        for r in scan(bam_path, ref_path, region)? {
            if seen.insert((r.name.clone(), r.first)) {
                records.push(r);
            }
        }
    }
    Ok(OriginSite {
        footprint: footprint.clone(),
        lookalikes,
        records,
        f: fragment_to_read_ratio(pool)?,
    })
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

    use rand::rngs::StdRng;
    use rand::SeedableRng;
    use std::collections::HashMap;

    fn chance(name: &str, family: u64, p: f64) -> Chance {
        Chance { name: name.to_string(), family: vec![FivePrime { chrom: "chr1".into(), pos: family, reverse: false }], chance: p }
    }

    #[test]
    fn test_removal_chance_is_p_origin_times_the_copy_rate() {
        // At the spot the phase call applies; at a look-alike the rate is
        // vaf -- even when read_copy carries an entry for that fragment's own
        // name, since at a look-alike the code must pass None rather than
        // consult it (R7's at_spot gate).
        let mut site = twin_site();
        // R4/R6: a stray whose unseen mate blocks removal (mate_unmapped:
        // false) must be filtered out of removal_chances entirely.
        let copy = site.records[0].clone();
        site.records.push(OriginRecord { name: "stray".into(), mate_unmapped: false, ..copy });
        let read_copy: HashMap<String, bool> =
            [("l0".to_string(), true), ("p0".to_string(), true)].into();
        let chances = site.removal_chances(&read_copy, 0.5);
        let of = |n: &str| chances.iter().find(|c| c.name == n).unwrap().chance;
        assert!(close(of("l0"), 0.5 * 1.0));
        assert!(close(of("l1"), 0.5 * 0.5));
        assert!(close(of("p0"), 0.5 * 0.5));
        assert!(!chances.iter().any(|c| c.name == "stray"), "{:?}", chances);
    }

    #[test]
    fn test_two_events_at_both_twins_remove_half_not_seven_sixteenths() {
        // R5: each event gives 1/2 x 1/2 = 1/4; together 1/2. Separate draws
        // would give 1 - (3/4)^2 = 7/16.
        let chances: Vec<Chance> = (0..20_000u64)
            .flat_map(|i| [chance(&format!("f{}", i), i, 0.25), chance(&format!("f{}", i), i, 0.25)])
            .collect();
        let removed = decide(&chances, &mut StdRng::seed_from_u64(1));
        let share = removed.len() as f64 / 20_000.0;
        assert!((0.48..0.52).contains(&share), "removed share {}", share);
    }

    #[test]
    fn test_a_total_above_one_always_removes() {
        let chances = [chance("a", 1, 0.7), chance("a", 1, 0.7)];
        for seed in 0..20 {
            assert!(decide(&chances, &mut StdRng::seed_from_u64(seed)).contains("a"));
        }
    }

    #[test]
    fn test_a_duplicate_shares_its_originals_fate() {
        // R3: 1000 families of two, each at 1/2.
        let chances: Vec<Chance> = (0..1000u64)
            .flat_map(|i| [chance(&format!("o{}", i), i, 0.5), chance(&format!("d{}", i), i, 0.5)])
            .collect();
        let removed = decide(&chances, &mut StdRng::seed_from_u64(2));
        let split = (0..1000).filter(|i| removed.contains(&format!("o{}", i)) != removed.contains(&format!("d{}", i))).count();
        assert_eq!(split, 0);
        assert!((400..600).contains(&(removed.len() / 2)), "{} removed", removed.len());
    }

    #[test]
    fn test_a_duplicate_takes_its_familys_phase_call() {
        // R7: phasing skips duplicates (src/loh.rs:628), so only the original
        // is in read_copy. At VAF 0.5 it gets rate 1, and so must its duplicate.
        let orig = OriginRecord { mate_unmapped: true, ..record("orig", true, 100, 60, &[]) };
        let dup = OriginRecord { name: "dup".into(), duplicate: true, ..orig.clone() };
        let site = OriginSite { footprint: fp(), lookalikes: vec![], records: vec![orig, dup], f: 1.0 };
        let on_event: HashMap<String, bool> = [("orig".to_string(), true)].into();
        let chances = site.removal_chances(&on_event, 0.5);
        assert_eq!(chances.len(), 2);
        assert!(close(chances[0].chance, chances[1].chance), "{:?}", chances);
        // On the other copy the original's rate is 0, and so is its duplicate's.
        let on_other: HashMap<String, bool> = [("orig".to_string(), false)].into();
        assert!(site.removal_chances(&on_other, 0.5).is_empty());
    }

    #[test]
    fn test_a_family_shares_one_fate_even_when_its_chances_differ() {
        // R7: one draw per family, against its highest member's total. Member
        // by member, a draw between 0.5 and 0.99 would remove only "orig".
        let chances = [chance("orig", 7, 0.99), chance("dup", 7, 0.5)];
        let mut removed_count = 0u32;
        for seed in 0..200 {
            let removed = decide(&chances, &mut StdRng::seed_from_u64(seed));
            assert_eq!(removed.contains("orig"), removed.contains("dup"), "seed {}", seed);
            if removed.contains("orig") {
                removed_count += 1;
            }
        }
        // Measured over these same 200 seeds, in a scratch copy of this crate
        // (/home/parlar_ai/edit-model-run/scratch/t6mut), by mutating the
        // aggregation at src/origin.rs:492-496:
        //   max (real code)                                   -> removed_count = 199
        //   first-member (`or_insert(total.min(1.0))`, no `.max`) -> removed_count = 103
        //   mean (average the members' totals instead of max) -> removed_count = 155
        // 180 sits strictly between the real count and the highest wrong one.
        assert!(removed_count > 180, "removed_count {}", removed_count);
    }

    use noodles::sam::alignment::record::{Flags, MappingQuality};
    use noodles::sam::alignment::RecordBuf;

    fn bam_record(name: &str, flags: u16, pos: usize, mapq: u8, xa: Option<&str>, mate_pos: usize) -> RecordBuf {
        use noodles::sam::alignment::record::data::field::Tag;
        use noodles::sam::alignment::record_buf::data::field::Value;
        use noodles::sam::alignment::record_buf::{Data, QualityScores, Sequence};
        let data: Data = xa.map(|xa| (Tag::new(b'X', b'A'), Value::from(xa))).into_iter().collect();
        RecordBuf::builder()
            .set_name(name)
            .set_flags(Flags::from(flags))
            .set_reference_sequence_id(0)
            .set_alignment_start(noodles::core::Position::new(pos).unwrap())
            .set_mapping_quality(MappingQuality::new(mapq).unwrap())
            .set_cigar([Op::new(Kind::Match, 100)].into_iter().collect())
            .set_mate_reference_sequence_id(0)
            .set_mate_alignment_start(noodles::core::Position::new(mate_pos).unwrap())
            .set_sequence(Sequence::from(vec![b'A'; 100]))
            .set_quality_scores(QualityScores::from(vec![30u8; 100]))
            .set_data(data)
            .build()
    }

    fn scratch(label: &str) -> std::path::PathBuf {
        let dir = std::env::temp_dir().join(format!("spike_origin_{}_{}", label, std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        dir
    }

    fn bases_pool() -> ReadPool {
        pool((0..30).map(|i| pair(&format!("r{}", i), 500_000 + 10 * i, 400, 150)).collect())
    }

    /// chrT, 60 kb; footprint chrT:20000-25000, twin P at 42000.
    /// a: unique at the spot. b: MAPQ 0 at the spot, hits at P. c: MAPQ 0 at
    /// P, hits at the spot. d: unique at P. e: read 1 at the spot, read 2
    /// past the footprint's end.
    ///
    /// The flag rule needs more than one flag pattern to be exercised, so the
    /// spot also holds `dup` (0x400 on both mates), `qcfail` (0x200 on both),
    /// `orphan` (read 1 whose mate is unmapped, 0x8), `improper` (a pair
    /// without 0x2), a secondary record of `a` (0x100) and a supplementary
    /// record `chimera` (0x800) whose primary lies outside the scan. Only the
    /// last two are dropped; the rest are candidates, `duplicate` and
    /// `qc_fail` carried through.
    ///
    /// `t` is a near look-alike: both its mates sit at the footprint's right
    /// edge with XA hits just past it, so the region they make -- grown by the
    /// read length -- reaches back inside the footprint, and `t`'s read 2 is
    /// scanned twice (the `seen` dedup).
    ///
    /// Records are in position order, and the first XA is still on the third
    /// (`first_xa_record`'s tests count it).
    fn twin_bam(dir: &std::path::Path) -> String {
        let (r1, r2) = (0x63u16, 0x93u16);
        let records = [
            bam_record("a", r1, 21_001, 60, None, 21_201),
            bam_record("a", r2, 21_201, 60, None, 21_001),
            bam_record("b", r1, 22_001, 0, Some("chrT,+42001,100M,0;"), 22_201),
            bam_record("b", r2, 22_201, 0, Some("chrT,-42201,100M,0;"), 22_001),
            bam_record("dup", r1 | 0x400, 23_001, 60, None, 23_201),
            bam_record("dup", r2 | 0x400, 23_201, 60, None, 23_001),
            bam_record("qcfail", r1 | 0x200, 23_301, 60, None, 23_501),
            bam_record("qcfail", r2 | 0x200, 23_501, 60, None, 23_301),
            bam_record("orphan", 0x49, 23_601, 60, None, 23_601),
            bam_record("improper", 0x41, 23_701, 60, None, 24_101),
            bam_record("improper", 0x91, 24_101, 60, None, 23_701),
            bam_record("t", r1, 24_401, 0, Some("chrT,+25051,100M,0;"), 24_951),
            bam_record("a", r1 | 0x100, 24_501, 60, None, 21_201),
            bam_record("chimera", r1 | 0x800, 24_601, 60, None, 24_801),
            bam_record("e", r1, 24_801, 60, None, 25_101),
            bam_record("t", r2, 24_951, 0, Some("chrT,-25251,100M,0;"), 24_401),
            bam_record("e", r2, 25_101, 60, None, 24_801),
            bam_record("c", r1, 42_001, 0, Some("chrT,+22001,100M,0;"), 42_201),
            bam_record("c", r2, 42_201, 0, Some("chrT,-22201,100M,0;"), 42_001),
            bam_record("d", r1, 42_301, 60, None, 42_501),
            bam_record("d", r2, 42_501, 60, None, 42_301),
        ];
        crate::extract::test_fixtures::write_one_contig_bam(&dir.join("twin.bam"), "chrT", 60_000, &records)
    }

    #[test]
    fn test_gather_reads_the_spot_and_its_lookalike() {
        let dir = scratch("gather");
        let bam = twin_bam(&dir);
        let site = gather(&bam, "", &Span::new("chrT", 20_000, 25_000), 100, &bases_pool()).unwrap();

        assert_eq!(
            site.lookalikes,
            vec![
                Span::new("chrT", 24_950, 25_450),
                Span::new("chrT", 41_900, 42_400),
            ]
        );
        let names: Vec<&str> = site.fragments().iter().map(|f| f.name).collect();
        assert_eq!(
            names,
            ["a", "b", "c", "dup", "e", "improper", "orphan", "qcfail", "t"]
        );
        assert_eq!(
            site.removable_names(),
            ["a", "b", "c", "dup", "improper", "orphan", "qcfail"]
        );
        let c = site.fragments().into_iter().find(|f| f.name == "c").unwrap();
        assert!(!c.at_spot(&site.footprint));
        assert!(close(c.chance(&site.footprint), 0.5));
        assert!(close(site.f, 400.0 / 300.0));

        // The candidate rule is on the flags, and `gather` is the only place
        // that reads them off a file. A secondary and a supplementary record
        // are not judged on their own -- they share their primary's name and
        // go with it -- so `a` arrives with two records, not three, and
        // `chimera`, whose primary lies outside the scan, with none.
        assert_eq!(site.records.iter().filter(|r| r.name == "a").count(), 2);
        assert!(!site.records.iter().any(|r| r.name == "chimera"));
        // Every other primary is a candidate whatever its flags say, and
        // `duplicate` and `qc_fail` arrive as themselves (R3 reads them).
        let record_of = |name: &str| -> &OriginRecord {
            site.records
                .iter()
                .find(|r| r.name == name && r.first)
                .unwrap_or_else(|| panic!("{} is not a candidate", name))
        };
        assert!(record_of("dup").duplicate && !record_of("dup").qc_fail);
        assert!(record_of("qcfail").qc_fail && !record_of("qcfail").duplicate);
        assert!(record_of("orphan").mate_unmapped);
        assert!(!record_of("improper").duplicate && !record_of("improper").qc_fail);

        // A look-alike region is grown by the read length, so the near one
        // reaches back inside the footprint and `t`'s read 2 is scanned twice.
        // The `seen` dedup must take it once: twice would satisfy
        // `mates.len() == 2` from one read and make it look like a pair.
        assert!(site.lookalikes[0].overlaps(&site.footprint));
        assert_eq!(
            site.records.iter().filter(|r| r.name == "t" && !r.first).count(),
            1
        );
        let _ = std::fs::remove_dir_all(&dir);
    }

    /// Two MAPQ 0 reads at the spot, with no XA anywhere in the file.
    fn noxa_bam(dir: &std::path::Path) -> String {
        let records = [
            bam_record("z", 0x63, 21_001, 0, None, 21_201),
            bam_record("z", 0x93, 21_201, 0, None, 21_001),
        ];
        crate::extract::test_fixtures::write_one_contig_bam(&dir.join("noxa.bam"), "chrT", 60_000, &records)
    }

    #[test]
    fn test_gather_does_not_stop_at_a_spot_whose_mapq0_reads_lack_xa() {
        // R8: under bwa-mem's -h 5 rule such reads have more than 5 hits, so
        // each gets 1/6. Whether the file kept XA at all is `require_xa`'s question.
        let dir = scratch("spot_noxa");
        let bam = noxa_bam(&dir);
        let site = gather(&bam, "", &Span::new("chrT", 20_000, 25_000), 100, &bases_pool()).unwrap();
        let z = site.fragments().into_iter().find(|f| f.name == "z").unwrap();
        assert!(close(z.chance(&site.footprint), 1.0 / 6.0));
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn test_first_xa_record_counts_records_up_to_the_first_xa() {
        // In `twin_bam` the first record with XA is b's read 1, the third.
        let dir = scratch("first_xa");
        let bam = twin_bam(&dir);
        assert_eq!(first_xa_record(&bam, "", 100).unwrap(), Some(3));
        assert_eq!(first_xa_record(&bam, "", 2).unwrap(), None);
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn test_require_xa_stops_a_file_without_xa() {
        // R8.
        let dir = scratch("require_xa");
        let err = require_xa(&noxa_bam(&dir), "").unwrap_err().to_string();
        assert!(err.contains("XA") && err.contains("--edit-model clean"), "{}", err);
        assert_eq!(require_xa(&twin_bam(&dir), "").unwrap(), 3);
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn test_first_xa_record_reads_a_cram() {
        // The two-contig CRAM carries SA:Z on every record and no XA.
        let dir = scratch("cram_xa");
        let (fasta, cram) = crate::extract::test_fixtures::write_two_contig_cram(&dir);
        assert_eq!(first_xa_record(&cram, &fasta, 100).unwrap(), None);
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn test_gather_reads_only_the_queried_contig_of_a_cram() {
        // A container holding two contigs is decoded whole (L2, N4).
        let dir = scratch("cram");
        let (fasta, cram) = crate::extract::test_fixtures::write_two_contig_cram(&dir);
        let site = gather(&cram, &fasta, &Span::new("chrA", 0, 1000), 100, &bases_pool()).unwrap();
        // The guard is the only thing that decides this: chrB's records come
        // out of the query and must be dropped during the scan. `fragments()`
        // drops them on chromosome anyway, so its names cannot show the guard.
        assert_eq!(site.records.len(), 6);
        assert!(site
            .records
            .iter()
            .flat_map(|r| r.placements.iter())
            .all(|p| p.span.chrom == "chrA"));
        let names: Vec<&str> = site.fragments().iter().map(|f| f.name).collect();
        assert_eq!(names, ["chrA_pair0", "chrA_pair1", "chrA_pair2"]);
        let _ = std::fs::remove_dir_all(&dir);
    }
}
