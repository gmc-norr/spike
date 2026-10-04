//! Variant haplotype construction: builds a linear sequence representing the
//! variant allele of any SV type.
//!
//! The variant genome is described as an ordered list of **segments**, each from
//! a reference region (possibly reverse-complemented) or novel sequence. Reads
//! tiled uniformly across this linear sequence are automatically chimeric when
//! they span segment boundaries — no per-SV-type breakpoint logic needed.

use anyhow::Result;

use crate::extract::reverse_complement;
use crate::reference::SharedReference;
use crate::types::FusionJoin;

/// Origin of a haplotype segment in the reference genome.
#[derive(Debug, Clone)]
pub struct SegmentOrigin {
    pub chrom: String,
    pub ref_start: u64, // 0-based
    /// 0-based, exclusive. Always `ref_start + sequence.len()`: a fetch that
    /// runs off a chromosome end comes back short, and a segment must not
    /// claim more reference than its sequence covers.
    pub ref_end: u64,
    pub is_reverse: bool, // true for reverse-complemented segments (INV)
}

/// Where a planted pair shows its event, in haplotype coordinates.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Evidence {
    /// A small variant's changed ALT bases `[start, end)`; a pure deletion
    /// has `start == end`, the join. A read shows them when it covers them
    /// and one base on each side, the span the carried-allele check uses.
    Bases { start: u64, end: u64 },
    /// A structural event's junctions. A pair shows one when its fragment
    /// covers the base on each side: a split read or a discordant pair.
    Junctions(Vec<u64>),
}

impl Evidence {
    /// Whether the pair at `spans` shows the event.
    pub fn shown_by(&self, spans: &crate::synth::PairSpans) -> bool {
        // `[s, e)` covers bases `a` through `b` when `s <= a` and `e > b`.
        let covers = |(s, e): (u64, u64), a: u64, b: u64| s <= a && e > b;
        match self {
            Evidence::Bases { start, end } => spans
                .mates
                .iter()
                .any(|&read| *start > 0 && covers(read, start - 1, *end)),
            Evidence::Junctions(junctions) => junctions
                .iter()
                .any(|&j| j > 0 && covers(spans.fragment, j - 1, j)),
        }
    }
}

/// Whether the reference runs straight on from segment `a` into segment `b`:
/// both from the same chromosome and strand, and `b` starting where `a` ends.
fn continues(a: &HaplotypeSegment, b: &HaplotypeSegment) -> bool {
    match (&a.origin, &b.origin) {
        (Some(x), Some(y)) if x.chrom == y.chrom && x.is_reverse == y.is_reverse => {
            if x.is_reverse {
                y.ref_end == x.ref_start
            } else {
                x.ref_end == y.ref_start
            }
        }
        _ => false,
    }
}

/// A segment of the variant haplotype.
#[derive(Debug, Clone)]
pub struct HaplotypeSegment {
    /// Uppercase DNA sequence for this segment.
    pub sequence: Vec<u8>,
    /// Reference origin, or None for novel insertions.
    pub origin: Option<SegmentOrigin>,
    /// Offset of this segment's first base in the linear haplotype.
    pub hap_offset: u64,
}

/// Complete variant haplotype: linear sequence assembled from segments.
///
/// The haplotype includes flanking reference sequence on both sides of the SV
/// so that reads tiled near the edges form complete pairs. Reads landing fully
/// within a single segment are normal reference reads; reads crossing segment
/// boundaries are chimeric.
#[derive(Clone)]
pub struct VariantHaplotype {
    pub segments: Vec<HaplotypeSegment>,
    /// Total length of the linear haplotype (sum of all segment lengths).
    pub total_len: u64,
    /// Concatenated sequence for fast subsequence access.
    sequence: Vec<u8>,
}

impl VariantHaplotype {
    /// Build from a list of segments.
    pub fn from_segments(mut segments: Vec<HaplotypeSegment>) -> Self {
        let mut offset = 0u64;
        for seg in &mut segments {
            seg.hap_offset = offset;
            offset += seg.sequence.len() as u64;
        }

        let sequence: Vec<u8> = segments
            .iter()
            .flat_map(|s| s.sequence.iter().copied())
            .collect();
        let total_len = sequence.len() as u64;

        Self {
            segments,
            total_len,
            sequence,
        }
    }

    /// Build a deletion haplotype.
    ///
    /// `ref[start-flank..del_start] | ref[del_end..del_end+flank]`
    pub fn from_deletion(
        reference: &SharedReference,
        chrom: &str,
        del_start: u64,
        del_end: u64,
        flank: u64,
    ) -> Result<Self> {
        let left_start = del_start.saturating_sub(flank);
        let right_end = del_end.saturating_add(flank);

        let left_seq = fetch_upper(reference, chrom, left_start, del_start)?;
        let right_seq = fetch_upper(reference, chrom, del_end, right_end)?;
        // Near a chromosome end the fetch returns fewer bases than asked for,
        // so each origin ends where its sequence ends.
        let left_ref_end = left_start + left_seq.len() as u64;
        let right_ref_end = del_end + right_seq.len() as u64;

        Ok(Self::from_segments(vec![
            HaplotypeSegment {
                sequence: left_seq,
                origin: Some(SegmentOrigin {
                    chrom: chrom.to_string(),
                    ref_start: left_start,
                    ref_end: left_ref_end,
                    is_reverse: false,
                }),
                hap_offset: 0,
            },
            HaplotypeSegment {
                sequence: right_seq,
                origin: Some(SegmentOrigin {
                    chrom: chrom.to_string(),
                    ref_start: del_end,
                    ref_end: right_ref_end,
                    is_reverse: false,
                }),
                hap_offset: 0,
            },
        ]))
    }

    /// Build a tandem duplication junction haplotype.
    ///
    /// The DUP junction: `ref[dup_end-flank..dup_end] | ref[dup_start..dup_start+flank]`
    ///
    /// This covers only the junction reads (chimeric at dup_end→dup_start).
    /// Depth increase inside the DUP region is handled separately by `generate_dup_depth_copies`.
    pub fn from_duplication(
        reference: &SharedReference,
        chrom: &str,
        dup_start: u64,
        dup_end: u64,
        flank: u64,
    ) -> Result<Self> {
        let left_start = dup_end.saturating_sub(flank);
        let right_end = dup_start.saturating_add(flank);

        let left_seq = fetch_upper(reference, chrom, left_start, dup_end)?;
        let right_seq = fetch_upper(reference, chrom, dup_start, right_end)?;
        let left_ref_end = left_start + left_seq.len() as u64;
        let right_ref_end = dup_start + right_seq.len() as u64;

        Ok(Self::from_segments(vec![
            HaplotypeSegment {
                sequence: left_seq,
                origin: Some(SegmentOrigin {
                    chrom: chrom.to_string(),
                    ref_start: left_start,
                    ref_end: left_ref_end,
                    is_reverse: false,
                }),
                hap_offset: 0,
            },
            HaplotypeSegment {
                sequence: right_seq,
                origin: Some(SegmentOrigin {
                    chrom: chrom.to_string(),
                    ref_start: dup_start,
                    ref_end: right_ref_end,
                    is_reverse: false,
                }),
                hap_offset: 0,
            },
        ]))
    }

    /// Build a full tandem duplication haplotype.
    ///
    /// `ref[dup_start-flank..dup_start] | ref[dup_start..dup_end] | ref[dup_start..dup_end] | ref[dup_end..dup_end+flank]`
    ///
    /// The duplicated region appears twice in tandem. Reads tiled across this
    /// haplotype naturally produce both junction evidence (split reads and
    /// discordant pairs at the copy1→copy2 boundary) and correct depth increase
    /// (2x in the DUP region from two copies).
    pub fn from_tandem_duplication(
        reference: &SharedReference,
        chrom: &str,
        dup_start: u64,
        dup_end: u64,
        flank: u64,
    ) -> Result<Self> {
        let left_start = dup_start.saturating_sub(flank);
        let right_end = dup_end.saturating_add(flank);

        let left_seq = fetch_upper(reference, chrom, left_start, dup_start)?;
        let dup_seq_1 = fetch_upper(reference, chrom, dup_start, dup_end)?;
        let dup_seq_2 = fetch_upper(reference, chrom, dup_start, dup_end)?;
        let right_seq = fetch_upper(reference, chrom, dup_end, right_end)?;
        let left_ref_end = left_start + left_seq.len() as u64;
        let dup_ref_end = dup_start + dup_seq_1.len() as u64;
        let right_ref_end = dup_end + right_seq.len() as u64;

        Ok(Self::from_segments(vec![
            HaplotypeSegment {
                sequence: left_seq,
                origin: Some(SegmentOrigin {
                    chrom: chrom.to_string(),
                    ref_start: left_start,
                    ref_end: left_ref_end,
                    is_reverse: false,
                }),
                hap_offset: 0,
            },
            HaplotypeSegment {
                sequence: dup_seq_1,
                origin: Some(SegmentOrigin {
                    chrom: chrom.to_string(),
                    ref_start: dup_start,
                    ref_end: dup_ref_end,
                    is_reverse: false,
                }),
                hap_offset: 0,
            },
            HaplotypeSegment {
                sequence: dup_seq_2,
                origin: Some(SegmentOrigin {
                    chrom: chrom.to_string(),
                    ref_start: dup_start,
                    ref_end: dup_ref_end,
                    is_reverse: false,
                }),
                hap_offset: 0,
            },
            HaplotypeSegment {
                sequence: right_seq,
                origin: Some(SegmentOrigin {
                    chrom: chrom.to_string(),
                    ref_start: dup_end,
                    ref_end: right_ref_end,
                    is_reverse: false,
                }),
                hap_offset: 0,
            },
        ]))
    }

    /// Build an inversion haplotype.
    ///
    /// `ref[inv_start-flank..inv_start] | revcomp(ref[inv_start..inv_end]) | ref[inv_end..inv_end+flank]`
    pub fn from_inversion(
        reference: &SharedReference,
        chrom: &str,
        inv_start: u64,
        inv_end: u64,
        flank: u64,
    ) -> Result<Self> {
        let left_start = inv_start.saturating_sub(flank);
        let right_end = inv_end.saturating_add(flank);

        let left_seq = fetch_upper(reference, chrom, left_start, inv_start)?;
        let mut inv_seq = fetch_upper(reference, chrom, inv_start, inv_end)?;
        reverse_complement(&mut inv_seq);
        let right_seq = fetch_upper(reference, chrom, inv_end, right_end)?;
        let left_ref_end = left_start + left_seq.len() as u64;
        let inv_ref_end = inv_start + inv_seq.len() as u64;
        let right_ref_end = inv_end + right_seq.len() as u64;

        Ok(Self::from_segments(vec![
            HaplotypeSegment {
                sequence: left_seq,
                origin: Some(SegmentOrigin {
                    chrom: chrom.to_string(),
                    ref_start: left_start,
                    ref_end: left_ref_end,
                    is_reverse: false,
                }),
                hap_offset: 0,
            },
            HaplotypeSegment {
                sequence: inv_seq,
                origin: Some(SegmentOrigin {
                    chrom: chrom.to_string(),
                    ref_start: inv_start,
                    ref_end: inv_ref_end,
                    is_reverse: true,
                }),
                hap_offset: 0,
            },
            HaplotypeSegment {
                sequence: right_seq,
                origin: Some(SegmentOrigin {
                    chrom: chrom.to_string(),
                    ref_start: inv_end,
                    ref_end: right_ref_end,
                    is_reverse: false,
                }),
                hap_offset: 0,
            },
        ]))
    }

    /// Build an insertion haplotype.
    ///
    /// `ref[pos-flank..pos] | ins_seq | ref[pos..pos+flank]`
    pub fn from_insertion(
        reference: &SharedReference,
        chrom: &str,
        pos: u64,
        ins_seq: &[u8],
        flank: u64,
    ) -> Result<Self> {
        let left_start = pos.saturating_sub(flank);
        let right_end = pos.saturating_add(flank);

        let left_seq = fetch_upper(reference, chrom, left_start, pos)?;
        let right_seq = fetch_upper(reference, chrom, pos, right_end)?;
        let left_ref_end = left_start + left_seq.len() as u64;
        let right_ref_end = pos + right_seq.len() as u64;

        // Uppercase the insertion sequence for consistency.
        let ins_upper: Vec<u8> = ins_seq.iter().map(|b| b.to_ascii_uppercase()).collect();

        Ok(Self::from_segments(vec![
            HaplotypeSegment {
                sequence: left_seq,
                origin: Some(SegmentOrigin {
                    chrom: chrom.to_string(),
                    ref_start: left_start,
                    ref_end: left_ref_end,
                    is_reverse: false,
                }),
                hap_offset: 0,
            },
            HaplotypeSegment {
                sequence: ins_upper,
                origin: None, // novel insertion
                hap_offset: 0,
            },
            HaplotypeSegment {
                sequence: right_seq,
                origin: Some(SegmentOrigin {
                    chrom: chrom.to_string(),
                    ref_start: pos,
                    ref_end: right_ref_end,
                    is_reverse: false,
                }),
                hap_offset: 0,
            },
        ]))
    }

    /// Build a fusion haplotype.
    ///
    /// `bp_a` and `bp_b` are cuts between 0-based bases `bp - 1` and `bp`:
    /// - Forward:    `ref_A[bp_a-flank..bp_a] | ref_B[bp_b..bp_b+flank]`
    /// - LeftLeft:   `ref_A[bp_a-flank..bp_a] | revcomp(ref_B[bp_b-flank..bp_b])`
    /// - RightRight: `revcomp(ref_A[bp_a..bp_a+flank]) | ref_B[bp_b..bp_b+flank]`
    pub fn from_fusion(
        reference: &SharedReference,
        chrom_a: &str,
        bp_a: u64,
        chrom_b: &str,
        bp_b: u64,
        flank: u64,
        join: FusionJoin,
    ) -> Result<Self> {
        // One piece of each gene: the side of its cut that the join keeps.
        let piece = |chrom: &str, bp: u64, keep_left: bool, reverse: bool| -> Result<_> {
            let (start, end) = if keep_left {
                (bp.saturating_sub(flank), bp)
            } else {
                (bp, bp.saturating_add(flank))
            };
            let mut sequence = fetch_upper(reference, chrom, start, end)?;
            if reverse {
                reverse_complement(&mut sequence);
            }
            let ref_end = start + sequence.len() as u64;
            Ok(HaplotypeSegment {
                sequence,
                origin: Some(SegmentOrigin {
                    chrom: chrom.to_string(),
                    ref_start: start,
                    ref_end,
                    is_reverse: reverse,
                }),
                hap_offset: 0,
            })
        };

        let segments = match join {
            FusionJoin::Forward => vec![
                piece(chrom_a, bp_a, true, false)?,
                piece(chrom_b, bp_b, false, false)?,
            ],
            FusionJoin::LeftLeft => vec![
                piece(chrom_a, bp_a, true, false)?,
                piece(chrom_b, bp_b, true, true)?,
            ],
            FusionJoin::RightRight => vec![
                piece(chrom_a, bp_a, false, true)?,
                piece(chrom_b, bp_b, false, false)?,
            ],
        };
        Ok(Self::from_segments(segments))
    }

    /// Build a small variant haplotype (SNP, MNV, or small indel).
    ///
    /// `ref[pos-flank..pos] | alt_allele | ref[pos+len(ref_allele)..pos+len(ref_allele)+flank]`
    ///
    /// Works for all small variant types:
    /// - SNP (A→T): left flank + [T] + right flank, skipping 1 ref base
    /// - Small del (ACG→A): left flank + [A] + right flank, skipping 3 ref bases
    /// - Small ins (A→ACGT): left flank + [ACGT] + right flank, skipping 1 ref base
    /// - MNV (AC→TG): left flank + [TG] + right flank, skipping 2 ref bases
    pub fn from_small_variant(
        reference: &SharedReference,
        chrom: &str,
        pos: u64,
        ref_allele: &[u8],
        alt_allele: &[u8],
        flank: u64,
    ) -> Result<Self> {
        let left_start = pos.saturating_sub(flank);
        let ref_end_pos = pos + ref_allele.len() as u64;
        let right_end = ref_end_pos.saturating_add(flank);

        let left_seq = fetch_upper(reference, chrom, left_start, pos)?;
        let right_seq = fetch_upper(reference, chrom, ref_end_pos, right_end)?;
        let left_ref_end = left_start + left_seq.len() as u64;
        let right_ref_end = ref_end_pos + right_seq.len() as u64;
        let alt_upper: Vec<u8> = alt_allele.iter().map(|b| b.to_ascii_uppercase()).collect();

        // For SNPs/MNVs (equal length ref and alt), the alt segment has a 1:1
        // mapping to the reference. Give it a SegmentOrigin so hap_to_ref works
        // correctly through the variant site (needed for read coordinate mapping).
        // For indels (different lengths), there's no 1:1 mapping → origin = None.
        let alt_origin = if ref_allele.len() == alt_allele.len() {
            Some(SegmentOrigin {
                chrom: chrom.to_string(),
                ref_start: pos,
                ref_end: ref_end_pos,
                is_reverse: false,
            })
        } else {
            None
        };

        Ok(Self::from_segments(vec![
            HaplotypeSegment {
                sequence: left_seq,
                origin: Some(SegmentOrigin {
                    chrom: chrom.to_string(),
                    ref_start: left_start,
                    ref_end: left_ref_end,
                    is_reverse: false,
                }),
                hap_offset: 0,
            },
            HaplotypeSegment {
                sequence: alt_upper,
                origin: alt_origin,
                hap_offset: 0,
            },
            HaplotypeSegment {
                sequence: right_seq,
                origin: Some(SegmentOrigin {
                    chrom: chrom.to_string(),
                    ref_start: ref_end_pos,
                    ref_end: right_ref_end,
                    is_reverse: false,
                }),
                hap_offset: 0,
            },
        ]))
    }

    /// Get a subsequence from the linear haplotype.
    ///
    /// Returns a slice of the concatenated sequence at `[hap_start, hap_start+len)`.
    /// Clamps to haplotype bounds.
    pub fn get_sequence(&self, hap_start: u64, len: usize) -> &[u8] {
        let start = hap_start as usize;
        let end = (start + len).min(self.sequence.len());
        let start = start.min(self.sequence.len());
        &self.sequence[start..end]
    }

    /// Map a haplotype position back to reference coordinates.
    ///
    /// Returns `(chrom, ref_pos)` if the position falls within a reference-derived
    /// segment, or `None` if it falls in a novel insertion.
    pub fn hap_to_ref(&self, hap_pos: u64) -> Option<(String, u64)> {
        for seg in &self.segments {
            let seg_end = seg.hap_offset + seg.sequence.len() as u64;
            if hap_pos >= seg.hap_offset && hap_pos < seg_end {
                let origin = seg.origin.as_ref()?;
                let offset_in_seg = hap_pos - seg.hap_offset;
                let ref_pos = if origin.is_reverse {
                    // Reverse segment: haplotype offset 0 corresponds to ref_end-1.
                    origin
                        .ref_end
                        .saturating_sub(1)
                        .saturating_sub(offset_in_seg)
                } else {
                    origin.ref_start + offset_in_seg
                };
                return Some((origin.chrom.clone(), ref_pos));
            }
        }
        None
    }

    /// Like [`hap_to_ref`](Self::hap_to_ref), but a position in novel
    /// (inserted) sequence maps to the nearest reference-mapped base.
    /// Returns None only when no segment maps to the reference.
    pub fn hap_to_ref_nearest(&self, hap_pos: u64) -> Option<(String, u64)> {
        if let Some(mapped) = self.hap_to_ref(hap_pos) {
            return Some(mapped);
        }
        let nearest = self
            .segments
            .iter()
            .filter(|seg| seg.origin.is_some() && !seg.sequence.is_empty())
            .map(|seg| {
                let last = seg.hap_offset + seg.sequence.len() as u64 - 1;
                let pos = hap_pos.clamp(seg.hap_offset, last);
                (pos.abs_diff(hap_pos), pos)
            })
            .min()?;
        self.hap_to_ref(nearest.1)
    }

    /// Get the breakpoint positions in the haplotype (segment boundaries).
    ///
    /// Returns the haplotype offsets where one segment ends and the next begins.
    pub fn breakpoints(&self) -> Vec<u64> {
        let mut bps = Vec::new();
        for seg in &self.segments[..self.segments.len().saturating_sub(1)] {
            bps.push(seg.hap_offset + seg.sequence.len() as u64);
        }
        bps
    }

    /// The breakpoints where the reference does not simply continue: one side
    /// is novel sequence, or the chromosome or strand changes, or the two
    /// sides' reference positions do not meet. A tandem DUP's flanks run
    /// straight into its copies, so only its copy-to-copy boundary is one.
    pub fn junctions(&self) -> Vec<u64> {
        self.segments
            .windows(2)
            .filter(|pair| !continues(&pair[0], &pair[1]))
            .map(|pair| pair[1].hap_offset)
            .collect()
    }

    /// Get the reference chrom of the first segment (for ReadPair.chrom).
    pub fn primary_chrom(&self) -> &str {
        self.segments
            .first()
            .and_then(|s| s.origin.as_ref())
            .map(|o| o.chrom.as_str())
            .unwrap_or("?")
    }

    /// Total length of reference-mapped segments (excludes novel insertions).
    ///
    /// Used for tiling fragment count: the number of tiled reads should match
    /// the number of reads suppressed from the reference footprint. Novel
    /// sequence (insertions) adds haplotype length but doesn't add reference
    /// coverage, so it shouldn't inflate the tiling count.
    pub fn ref_mapped_len(&self) -> u64 {
        self.segments
            .iter()
            .filter_map(|seg| seg.origin.as_ref())
            .map(|o| o.ref_end.saturating_sub(o.ref_start))
            .sum()
    }

    /// Check if a haplotype range [start, start+len) crosses a segment boundary.
    ///
    /// Returns true if the range falls entirely within a single segment, false
    /// if it spans two or more segments.
    #[cfg(test)]
    pub fn is_within_single_segment(&self, start: u64, len: u64) -> bool {
        let end = start + len; // exclusive
        for seg in &self.segments {
            let seg_start = seg.hap_offset;
            let seg_end = seg.hap_offset + seg.sequence.len() as u64;
            if start >= seg_start && end <= seg_end {
                return true;
            }
        }
        false
    }

    /// Get the reference coordinate range covered by all segments.
    ///
    /// Returns `(min_ref_pos, max_ref_pos)` across all reference-derived segments.
    /// Used to scope read suppression: only reads within this range need
    /// suppression, since tiled haplotype reads only cover this footprint.
    pub fn ref_range(&self) -> Option<(u64, u64)> {
        let mut min_pos = u64::MAX;
        let mut max_pos = 0u64;
        let mut any = false;

        for seg in &self.segments {
            if let Some(origin) = &seg.origin {
                min_pos = min_pos.min(origin.ref_start);
                max_pos = max_pos.max(origin.ref_end);
                any = true;
            }
        }

        if any {
            Some((min_pos, max_pos))
        } else {
            None
        }
    }

    /// Apply the sample's SNP alleles on `chrom` to the haplotype sequence.
    ///
    /// `variants` maps reference position → allele base (uppercase). For each
    /// segment from `chrom`, every base whose reference position appears in the
    /// map is replaced with the allele (complemented in reversed segments). This
    /// makes reads tiled across the haplotype carry the sample's alleles instead
    /// of reference-only bases.
    pub fn apply_variants(&mut self, chrom: &str, variants: &std::collections::HashMap<u64, u8>) {
        if variants.is_empty() {
            return;
        }
        let mut applied = 0usize;
        for seg in &mut self.segments {
            let origin = match &seg.origin {
                Some(o) if o.chrom == chrom => o,
                _ => continue, // novel insertion, or another chromosome
            };
            let seg_len = seg.sequence.len() as u64;
            for offset in 0..seg_len {
                let ref_pos = if origin.is_reverse {
                    origin.ref_end.saturating_sub(1).saturating_sub(offset)
                } else {
                    origin.ref_start + offset
                };
                if let Some(&alt) = variants.get(&ref_pos) {
                    let base = if origin.is_reverse {
                        crate::extract::complement_base(alt)
                    } else {
                        alt
                    };
                    seg.sequence[offset as usize] = base;
                    // Also update the concatenated sequence.
                    let concat_idx = (seg.hap_offset + offset) as usize;
                    if concat_idx < self.sequence.len() {
                        self.sequence[concat_idx] = base;
                    }
                    applied += 1;
                }
            }
        }
        if applied > 0 {
            log::info!("Applied {} het SNP variants to haplotype sequence", applied,);
        }
    }

    /// Check if a haplotype range [start, start+len) overlaps any reference-mapped segment.
    ///
    /// Returns true if at least one base in the range comes from a reference segment
    /// (i.e., is not novel insertion sequence). Tiling placement no longer calls this
    /// directly -- `sample_ref_overlapping_start` in simulate.rs draws only from starts
    /// that satisfy this predicate -- but the tests use it as their oracle to check that
    /// sampler's output.
    #[cfg(test)]
    pub fn overlaps_ref_segment(&self, start: u64, len: u64) -> bool {
        let end = start + len;
        for seg in &self.segments {
            if seg.origin.is_none() {
                continue; // skip novel segments
            }
            let seg_end = seg.hap_offset + seg.sequence.len() as u64;
            // Overlap check: [start, end) ∩ [seg.hap_offset, seg_end) is non-empty.
            if start < seg_end && end > seg.hap_offset {
                return true;
            }
        }
        false
    }
}

/// Fetch reference sequence and uppercase it.
fn fetch_upper(reference: &SharedReference, chrom: &str, start: u64, end: u64) -> Result<Vec<u8>> {
    if start >= end {
        return Ok(Vec::new());
    }
    let seq = reference.fetch_sequence(chrom, start, end)?;
    Ok(seq.iter().map(|b| b.to_ascii_uppercase()).collect())
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Mock reference: returns a repeating ACGT pattern.
    struct MockRef;

    impl MockRef {
        fn fetch(&self, _chrom: &str, start: u64, end: u64) -> Vec<u8> {
            let pattern = b"ACGT";
            (start..end).map(|i| pattern[(i % 4) as usize]).collect()
        }
    }

    // Helper: build segments from mock reference since we can't use SharedReference
    // in unit tests without a real FASTA.
    fn mock_segments_deletion(del_start: u64, del_end: u64, flank: u64) -> VariantHaplotype {
        let mock = MockRef;
        let left_start = del_start.saturating_sub(flank);
        let right_end = del_end.saturating_add(flank);

        VariantHaplotype::from_segments(vec![
            HaplotypeSegment {
                sequence: mock.fetch("chr1", left_start, del_start),
                origin: Some(SegmentOrigin {
                    chrom: "chr1".to_string(),
                    ref_start: left_start,
                    ref_end: del_start,
                    is_reverse: false,
                }),
                hap_offset: 0,
            },
            HaplotypeSegment {
                sequence: mock.fetch("chr1", del_end, right_end),
                origin: Some(SegmentOrigin {
                    chrom: "chr1".to_string(),
                    ref_start: del_end,
                    ref_end: right_end,
                    is_reverse: false,
                }),
                hap_offset: 0,
            },
        ])
    }

    fn mock_segments_inversion(inv_start: u64, inv_end: u64, flank: u64) -> VariantHaplotype {
        let mock = MockRef;
        let left_start = inv_start.saturating_sub(flank);
        let right_end = inv_end.saturating_add(flank);

        let mut inv_seq = mock.fetch("chr1", inv_start, inv_end);
        reverse_complement(&mut inv_seq);

        VariantHaplotype::from_segments(vec![
            HaplotypeSegment {
                sequence: mock.fetch("chr1", left_start, inv_start),
                origin: Some(SegmentOrigin {
                    chrom: "chr1".to_string(),
                    ref_start: left_start,
                    ref_end: inv_start,
                    is_reverse: false,
                }),
                hap_offset: 0,
            },
            HaplotypeSegment {
                sequence: inv_seq,
                origin: Some(SegmentOrigin {
                    chrom: "chr1".to_string(),
                    ref_start: inv_start,
                    ref_end: inv_end,
                    is_reverse: true,
                }),
                hap_offset: 0,
            },
            HaplotypeSegment {
                sequence: mock.fetch("chr1", inv_end, right_end),
                origin: Some(SegmentOrigin {
                    chrom: "chr1".to_string(),
                    ref_start: inv_end,
                    ref_end: right_end,
                    is_reverse: false,
                }),
                hap_offset: 0,
            },
        ])
    }

    fn mock_segments_insertion(pos: u64, ins_seq: &[u8], flank: u64) -> VariantHaplotype {
        let mock = MockRef;
        let left_start = pos.saturating_sub(flank);
        let right_end = pos.saturating_add(flank);

        VariantHaplotype::from_segments(vec![
            HaplotypeSegment {
                sequence: mock.fetch("chr1", left_start, pos),
                origin: Some(SegmentOrigin {
                    chrom: "chr1".to_string(),
                    ref_start: left_start,
                    ref_end: pos,
                    is_reverse: false,
                }),
                hap_offset: 0,
            },
            HaplotypeSegment {
                sequence: ins_seq.to_vec(),
                origin: None,
                hap_offset: 0,
            },
            HaplotypeSegment {
                sequence: mock.fetch("chr1", pos, right_end),
                origin: Some(SegmentOrigin {
                    chrom: "chr1".to_string(),
                    ref_start: pos,
                    ref_end: right_end,
                    is_reverse: false,
                }),
                hap_offset: 0,
            },
        ])
    }

    #[test]
    fn test_deletion_haplotype_length() {
        let hap = mock_segments_deletion(1000, 2000, 500);
        // Left flank: [500, 1000) = 500bp. Right flank: [2000, 2500) = 500bp.
        assert_eq!(hap.total_len, 1000);
        assert_eq!(hap.segments.len(), 2);
    }

    #[test]
    fn test_deletion_breakpoint() {
        let hap = mock_segments_deletion(1000, 2000, 500);
        let bps = hap.breakpoints();
        assert_eq!(bps.len(), 1);
        assert_eq!(bps[0], 500); // breakpoint at offset 500 (end of left flank)
    }

    #[test]
    fn test_deletion_hap_to_ref() {
        let hap = mock_segments_deletion(1000, 2000, 500);

        // Position 0 in haplotype → ref position 500 (start of left flank).
        let (chrom, pos) = hap.hap_to_ref(0).unwrap();
        assert_eq!(chrom, "chr1");
        assert_eq!(pos, 500);

        // Position 499 → ref 999 (last base before deletion).
        let (_, pos) = hap.hap_to_ref(499).unwrap();
        assert_eq!(pos, 999);

        // Position 500 → ref 2000 (first base after deletion).
        let (_, pos) = hap.hap_to_ref(500).unwrap();
        assert_eq!(pos, 2000);

        // Position 999 → ref 2499.
        let (_, pos) = hap.hap_to_ref(999).unwrap();
        assert_eq!(pos, 2499);
    }

    #[test]
    fn test_deletion_get_sequence() {
        let hap = mock_segments_deletion(1000, 2000, 100);
        let seq = hap.get_sequence(0, 10);
        assert_eq!(seq.len(), 10);
        // MockRef pattern: ACGT... position 900 → A(900%4=0), C, G, T, A, C, G, T, A, C
        let expected: Vec<u8> = (900u64..910).map(|i| b"ACGT"[(i % 4) as usize]).collect();
        assert_eq!(seq, &expected[..]);
    }

    #[test]
    fn test_deletion_sequence_spans_breakpoint() {
        let hap = mock_segments_deletion(1000, 2000, 100);
        // Get sequence crossing the breakpoint at offset 100.
        let seq = hap.get_sequence(95, 10);
        assert_eq!(seq.len(), 10);
        // Positions 95-99: ref [995, 1000), positions 100-104: ref [2000, 2005)
        let mut expected = Vec::new();
        for i in 995u64..1000 {
            expected.push(b"ACGT"[(i % 4) as usize]);
        }
        for i in 2000u64..2005 {
            expected.push(b"ACGT"[(i % 4) as usize]);
        }
        assert_eq!(seq, &expected[..]);
    }

    #[test]
    fn test_inversion_haplotype() {
        let hap = mock_segments_inversion(100, 200, 50);
        // 3 segments: left [50,100), inverted [100,200), right [200,250)
        assert_eq!(hap.segments.len(), 3);
        assert_eq!(hap.total_len, 200); // 50 + 100 + 50

        let bps = hap.breakpoints();
        assert_eq!(bps.len(), 2);
        assert_eq!(bps[0], 50); // left→inverted boundary
        assert_eq!(bps[1], 150); // inverted→right boundary
    }

    #[test]
    fn test_inversion_coordinate_mapping() {
        let hap = mock_segments_inversion(100, 200, 50);

        // Left flank: hap pos 0 → ref 50.
        let (_, pos) = hap.hap_to_ref(0).unwrap();
        assert_eq!(pos, 50);

        // Inverted region: hap pos 50 → ref 199 (inv_end - 1, first base of inverted).
        let (_, pos) = hap.hap_to_ref(50).unwrap();
        assert_eq!(pos, 199);

        // Inverted region: hap pos 149 → ref 100 (inv_start, last base of inverted).
        let (_, pos) = hap.hap_to_ref(149).unwrap();
        assert_eq!(pos, 100);

        // Right flank: hap pos 150 → ref 200.
        let (_, pos) = hap.hap_to_ref(150).unwrap();
        assert_eq!(pos, 200);
    }

    #[test]
    fn test_insertion_haplotype() {
        let ins_seq = b"NNNNNNNNNN"; // 10bp insertion
        let hap = mock_segments_insertion(500, ins_seq, 100);

        // 3 segments: left [400,500), insert (10bp), right [500,600)
        assert_eq!(hap.segments.len(), 3);
        assert_eq!(hap.total_len, 210); // 100 + 10 + 100

        let bps = hap.breakpoints();
        assert_eq!(bps.len(), 2);
        assert_eq!(bps[0], 100); // left→insert boundary
        assert_eq!(bps[1], 110); // insert→right boundary
    }

    #[test]
    fn test_insertion_coordinate_mapping() {
        let ins_seq = b"NNNNNNNNNN";
        let hap = mock_segments_insertion(500, ins_seq, 100);

        // Left: hap 0 → ref 400.
        let (_, pos) = hap.hap_to_ref(0).unwrap();
        assert_eq!(pos, 400);

        // Last left: hap 99 → ref 499.
        let (_, pos) = hap.hap_to_ref(99).unwrap();
        assert_eq!(pos, 499);

        // Insert: hap 100-109 → None (novel sequence).
        assert!(hap.hap_to_ref(100).is_none());
        assert!(hap.hap_to_ref(109).is_none());

        // Right: hap 110 → ref 500.
        let (_, pos) = hap.hap_to_ref(110).unwrap();
        assert_eq!(pos, 500);
    }

    #[test]
    fn test_primary_chrom() {
        let hap = mock_segments_deletion(1000, 2000, 500);
        assert_eq!(hap.primary_chrom(), "chr1");
    }

    #[test]
    fn test_empty_flank() {
        // flank=0: haplotype is just the two sides of the deletion with no padding.
        let hap = mock_segments_deletion(100, 200, 0);
        assert_eq!(hap.total_len, 0);
    }

    #[test]
    fn test_get_sequence_out_of_bounds() {
        let hap = mock_segments_deletion(1000, 2000, 100);
        // Request beyond end.
        let seq = hap.get_sequence(195, 20);
        assert_eq!(seq.len(), 5); // only 5 bases left
    }

    // Helper for small variant haplotypes using MockRef.
    fn mock_segments_small_variant(
        pos: u64,
        ref_allele: &[u8],
        alt_allele: &[u8],
        flank: u64,
    ) -> VariantHaplotype {
        let mock = MockRef;
        let left_start = pos.saturating_sub(flank);
        let ref_end_pos = pos + ref_allele.len() as u64;
        let right_end = ref_end_pos.saturating_add(flank);

        // SNPs/MNVs (equal length) get a reference origin; indels do not.
        let alt_origin = if ref_allele.len() == alt_allele.len() {
            Some(SegmentOrigin {
                chrom: "chr1".to_string(),
                ref_start: pos,
                ref_end: ref_end_pos,
                is_reverse: false,
            })
        } else {
            None
        };

        VariantHaplotype::from_segments(vec![
            HaplotypeSegment {
                sequence: mock.fetch("chr1", left_start, pos),
                origin: Some(SegmentOrigin {
                    chrom: "chr1".to_string(),
                    ref_start: left_start,
                    ref_end: pos,
                    is_reverse: false,
                }),
                hap_offset: 0,
            },
            HaplotypeSegment {
                sequence: alt_allele.to_vec(),
                origin: alt_origin,
                hap_offset: 0,
            },
            HaplotypeSegment {
                sequence: mock.fetch("chr1", ref_end_pos, right_end),
                origin: Some(SegmentOrigin {
                    chrom: "chr1".to_string(),
                    ref_start: ref_end_pos,
                    ref_end: right_end,
                    is_reverse: false,
                }),
                hap_offset: 0,
            },
        ])
    }

    #[test]
    fn test_snp_haplotype() {
        // SNP at pos 500: ref=A, alt=T, flank=100
        let hap = mock_segments_small_variant(500, b"A", b"T", 100);
        // 3 segments: left [400,500) = 100bp, alt [T] = 1bp, right [501,601) = 100bp
        assert_eq!(hap.segments.len(), 3);
        assert_eq!(hap.total_len, 201); // 100 + 1 + 100

        let bps = hap.breakpoints();
        assert_eq!(bps.len(), 2);
        assert_eq!(bps[0], 100); // left→alt boundary
        assert_eq!(bps[1], 101); // alt→right boundary

        // The alt base should be T.
        let alt_seq = hap.get_sequence(100, 1);
        assert_eq!(alt_seq, b"T");
    }

    #[test]
    fn test_snp_coordinate_mapping() {
        let hap = mock_segments_small_variant(500, b"A", b"T", 100);

        // Left: hap 0 → ref 400.
        let (_, pos) = hap.hap_to_ref(0).unwrap();
        assert_eq!(pos, 400);

        // Last left: hap 99 → ref 499.
        let (_, pos) = hap.hap_to_ref(99).unwrap();
        assert_eq!(pos, 499);

        // SNP alt base at hap 100 → ref 500 (SNPs have 1:1 ref mapping).
        let (_, pos) = hap.hap_to_ref(100).unwrap();
        assert_eq!(pos, 500);

        // Right: hap 101 → ref 501.
        let (_, pos) = hap.hap_to_ref(101).unwrap();
        assert_eq!(pos, 501);
    }

    #[test]
    fn test_small_deletion_haplotype() {
        // Small del at pos 500: ref=ACG (3bp), alt=A (1bp) → delete CG
        let hap = mock_segments_small_variant(500, b"ACG", b"A", 100);
        // 3 segments: left [400,500)=100, alt [A]=1, right [503,603)=100
        assert_eq!(hap.segments.len(), 3);
        assert_eq!(hap.total_len, 201); // 100 + 1 + 100 (shorter than ref span of 203)

        // Right flank starts at ref 503 (pos + 3).
        let (_, pos) = hap.hap_to_ref(101).unwrap();
        assert_eq!(pos, 503);
    }

    #[test]
    fn test_small_insertion_haplotype() {
        // Small ins at pos 500: ref=A (1bp), alt=ACGT (4bp) → insert CGT
        let hap = mock_segments_small_variant(500, b"A", b"ACGT", 100);
        // 3 segments: left [400,500)=100, alt [ACGT]=4, right [501,601)=100
        assert_eq!(hap.segments.len(), 3);
        assert_eq!(hap.total_len, 204); // 100 + 4 + 100

        // Inserted sequence should be ACGT.
        let ins_seq = hap.get_sequence(100, 4);
        assert_eq!(ins_seq, b"ACGT");

        // Right flank starts at ref 501 (pos + 1).
        let (_, pos) = hap.hap_to_ref(104).unwrap();
        assert_eq!(pos, 501);
    }

    // ── is_within_single_segment tests ──────────────────────────────────

    #[test]
    fn test_within_single_segment_left() {
        // DEL: ref[500..1500) deleted, flanks 500bp each → hap [0..500) + [500..1000)
        let hap = mock_segments_deletion(500, 1500, 500);
        assert_eq!(hap.total_len, 1000);

        // Fully within left segment [0..500)
        assert!(hap.is_within_single_segment(0, 150));
        assert!(hap.is_within_single_segment(200, 150));
        assert!(hap.is_within_single_segment(350, 150)); // ends at 500 = boundary (exclusive, OK)
    }

    #[test]
    fn test_within_single_segment_right() {
        let hap = mock_segments_deletion(500, 1500, 500);

        // Fully within right segment [500..1000)
        assert!(hap.is_within_single_segment(500, 150));
        assert!(hap.is_within_single_segment(700, 150));
        assert!(hap.is_within_single_segment(850, 150)); // ends at 1000
    }

    #[test]
    fn test_within_single_segment_crossing_boundary() {
        let hap = mock_segments_deletion(500, 1500, 500);

        // Crosses the boundary at 500
        assert!(!hap.is_within_single_segment(400, 200)); // 400..600 crosses 500
        assert!(!hap.is_within_single_segment(490, 150)); // 490..640 crosses 500
        assert!(!hap.is_within_single_segment(351, 150)); // 351..501 crosses 500
    }

    #[test]
    fn test_within_single_segment_three_segments() {
        // INV: left flank + inverted middle + right flank
        let hap = mock_segments_inversion(1000, 2000, 500);
        // 3 segments: [0..500) [500..1500) [1500..2000)
        assert_eq!(hap.segments.len(), 3);
        assert_eq!(hap.total_len, 2000);

        // Within left segment
        assert!(hap.is_within_single_segment(0, 150));
        assert!(hap.is_within_single_segment(350, 150)); // ends at 500

        // Within inverted middle
        assert!(hap.is_within_single_segment(500, 150));
        assert!(hap.is_within_single_segment(1000, 150));
        assert!(hap.is_within_single_segment(1350, 150)); // ends at 1500

        // Within right segment
        assert!(hap.is_within_single_segment(1500, 150));
        assert!(hap.is_within_single_segment(1850, 150));

        // Crosses left→middle boundary at 500
        assert!(!hap.is_within_single_segment(450, 150));
        // Crosses middle→right boundary at 1500
        assert!(!hap.is_within_single_segment(1450, 150));
    }

    // ── tandem duplication tests ──────────────────────────────────────────

    fn mock_tandem_dup(dup_start: u64, dup_end: u64, flank: u64) -> VariantHaplotype {
        let pattern = b"ACGT";
        let left_start = dup_start.saturating_sub(flank);
        let right_end = dup_end + flank;

        let make_seq = |start: u64, end: u64| -> Vec<u8> {
            (start..end).map(|i| pattern[(i % 4) as usize]).collect()
        };

        VariantHaplotype::from_segments(vec![
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
    fn test_tandem_dup_haplotype_structure() {
        // DUP [1000, 2000) with 500bp flanks
        let hap = mock_tandem_dup(1000, 2000, 500);

        assert_eq!(hap.segments.len(), 4);
        // total_len = flank + dup + dup + flank = 500 + 1000 + 1000 + 500 = 3000
        assert_eq!(hap.total_len, 3000);

        // 3 breakpoints between 4 segments
        let bps = hap.breakpoints();
        assert_eq!(bps.len(), 3);
        assert_eq!(bps[0], 500);  // left_flank → copy1
        assert_eq!(bps[1], 1500); // copy1 → copy2 (the novel junction)
        assert_eq!(bps[2], 2500); // copy2 → right_flank
    }

    #[test]
    fn test_tandem_dup_ref_mapped_len() {
        let hap = mock_tandem_dup(1000, 2000, 500);
        // All 4 segments are reference-mapped
        // ref_mapped_len = 500 + 1000 + 1000 + 500 = 3000
        assert_eq!(hap.ref_mapped_len(), 3000);
    }

    #[test]
    fn test_tandem_dup_ref_range() {
        let hap = mock_tandem_dup(1000, 2000, 500);
        let (min_ref, max_ref) = hap.ref_range().unwrap();
        // min = dup_start - flank = 500
        assert_eq!(min_ref, 500);
        // max = dup_end + flank = 2500
        assert_eq!(max_ref, 2500);
    }

    #[test]
    fn test_tandem_dup_hap_to_ref_both_copies() {
        let hap = mock_tandem_dup(1000, 2000, 500);

        // Position in copy1 (offset 500 = start of copy1)
        let (chrom, ref_pos) = hap.hap_to_ref(500).unwrap();
        assert_eq!(chrom, "chr1");
        assert_eq!(ref_pos, 1000); // dup_start

        // Position in copy2 (offset 1500 = start of copy2)
        let (chrom, ref_pos) = hap.hap_to_ref(1500).unwrap();
        assert_eq!(chrom, "chr1");
        assert_eq!(ref_pos, 1000); // same ref coord as copy1 start

        // Mid-copy1 (offset 1000 = 500 into the dup region)
        let (_, ref_pos) = hap.hap_to_ref(1000).unwrap();
        assert_eq!(ref_pos, 1500);

        // Mid-copy2 (offset 2000 = 500 into the dup region in copy2)
        let (_, ref_pos) = hap.hap_to_ref(2000).unwrap();
        assert_eq!(ref_pos, 1500); // same ref coord
    }

    #[test]
    fn test_tandem_dup_junction_sequence() {
        // Verify the sequence at the copy1→copy2 junction is correct:
        // end of copy1 (ref near dup_end) followed by start of copy2 (ref near dup_start)
        let hap = mock_tandem_dup(1000, 2000, 500);
        let pattern = b"ACGT";

        // Last 10bp of copy1 = ref[1990..2000)
        let end_copy1 = hap.get_sequence(1490, 10);
        let expected_end: Vec<u8> = (1990u64..2000).map(|i| pattern[(i % 4) as usize]).collect();
        assert_eq!(end_copy1, &expected_end[..]);

        // First 10bp of copy2 = ref[1000..1010)
        let start_copy2 = hap.get_sequence(1500, 10);
        let expected_start: Vec<u8> = (1000u64..1010).map(|i| pattern[(i % 4) as usize]).collect();
        assert_eq!(start_copy2, &expected_start[..]);
    }

    // ---------------------------------------------------------------
    // from_fusion geometry, built from a real (in-memory) reference
    // ---------------------------------------------------------------

    /// chrA and chrB are distinct, non-palindromic 16 bp sequences.
    fn fusion_ref() -> SharedReference {
        let mut seqs = std::collections::HashMap::new();
        seqs.insert("chrA".to_string(), b"TTTTAAAACAGGTTTT".to_vec());
        seqs.insert("chrB".to_string(), b"GATTACACATGGCCAT".to_vec());
        SharedReference::from_sequences(seqs)
    }

    fn ref_pos(hap: &VariantHaplotype, hap_pos: u64) -> (String, u64) {
        hap.hap_to_ref(hap_pos).unwrap()
    }

    #[test]
    fn test_apply_variants_changes_only_the_named_chromosome() {
        // Two segments at the same positions on different chromosomes.
        let seg = |chrom: &str| HaplotypeSegment {
            sequence: b"AAAA".to_vec(),
            origin: Some(SegmentOrigin {
                chrom: chrom.to_string(),
                ref_start: 100,
                ref_end: 104,
                is_reverse: false,
            }),
            hap_offset: 0,
        };
        let mut hap = VariantHaplotype::from_segments(vec![seg("chrA"), seg("chrB")]);
        hap.apply_variants("chrB", &[(101, b'G')].into());
        assert_eq!(hap.get_sequence(0, 8), b"AAAAAGAA");
    }

    #[test]
    fn test_fusion_forward_joins_a_left_to_b_right() {
        // A[4..8) = AAAA, then B[8..12) = ATGG.
        let hap = VariantHaplotype::from_fusion(
            &fusion_ref(), "chrA", 8, "chrB", 8, 4, FusionJoin::Forward,
        )
        .unwrap();
        assert_eq!(hap.sequence, b"AAAAATGG");
        assert_eq!(ref_pos(&hap, 3), ("chrA".to_string(), 7));
        assert_eq!(ref_pos(&hap, 4), ("chrB".to_string(), 8));
    }

    #[test]
    fn test_fusion_left_left_joins_a_left_to_revcomp_of_b_left() {
        // A[4..8) = AAAA, then revcomp(B[4..8) = ACAC) = GTGT.
        // The base right after the junction is B's base just left of the cut.
        let hap = VariantHaplotype::from_fusion(
            &fusion_ref(), "chrA", 8, "chrB", 8, 4, FusionJoin::LeftLeft,
        )
        .unwrap();
        assert_eq!(hap.sequence, b"AAAAGTGT");
        assert_eq!(ref_pos(&hap, 3), ("chrA".to_string(), 7));
        assert_eq!(ref_pos(&hap, 4), ("chrB".to_string(), 7));
        assert_eq!(ref_pos(&hap, 7), ("chrB".to_string(), 4));
    }

    #[test]
    fn test_fusion_right_right_joins_revcomp_of_a_right_to_b_right() {
        // revcomp(A[8..12) = CAGG) = CCTG, then B[8..12) = ATGG.
        // The base right before the junction is A's base just right of the cut.
        let hap = VariantHaplotype::from_fusion(
            &fusion_ref(), "chrA", 8, "chrB", 8, 4, FusionJoin::RightRight,
        )
        .unwrap();
        assert_eq!(hap.sequence, b"CCTGATGG");
        assert_eq!(ref_pos(&hap, 0), ("chrA".to_string(), 11));
        assert_eq!(ref_pos(&hap, 3), ("chrA".to_string(), 8));
        assert_eq!(ref_pos(&hap, 4), ("chrB".to_string(), 8));
    }

    // ---------------------------------------------------------------
    // Chromosome-end clamping (M4), built from a real (in-memory) reference
    // ---------------------------------------------------------------

    /// A 10 kb contig: anything a constructor asks for past 10 000 comes back
    /// short, because `SharedReference::fetch_sequence` clamps.
    fn short_contig() -> SharedReference {
        let pattern = b"ACGT";
        let seq: Vec<u8> = (0..10_000u64).map(|i| pattern[(i % 4) as usize]).collect();
        let mut seqs = std::collections::HashMap::new();
        seqs.insert("chrEnd".to_string(), seq);
        SharedReference::from_sequences(seqs)
    }

    #[test]
    fn test_deletion_at_chromosome_end_does_not_overcount_ref_len() {
        // A DEL ending 100 bp from the chromosome end, 2 kb flanks:
        // left [1000, 3000) = 2000 bp, right [9900, 11900) → [9900, 10000) = 100 bp.
        let hap =
            VariantHaplotype::from_deletion(&short_contig(), "chrEnd", 3000, 9900, 2000).unwrap();
        // 2000 + 100 = 2100 bases of reference, not 2000 + 2000 = 4000.
        assert_eq!(hap.total_len, 2100);
        assert_eq!(hap.ref_mapped_len(), 2100);
        // The footprint ends at the contig end, not 1900 bp past it.
        assert_eq!(hap.ref_range(), Some((1000, 10_000)));
    }

    #[test]
    fn test_junction_duplication_at_chromosome_end_does_not_overcount_ref_len() {
        // DUP [9500, 9900) junction, 2 kb flanks: left [7900, 9900) = 2000 bp,
        // right [9500, 11500) → [9500, 10000) = 500 bp.
        let hap =
            VariantHaplotype::from_duplication(&short_contig(), "chrEnd", 9500, 9900, 2000).unwrap();
        // 2000 + 500 = 2500, not 2000 + 2000 = 4000.
        assert_eq!(hap.total_len, 2500);
        assert_eq!(hap.ref_mapped_len(), 2500);
        assert_eq!(hap.ref_range(), Some((7900, 10_000)));
    }

    #[test]
    fn test_tandem_duplication_at_chromosome_end_does_not_overcount_ref_len() {
        // DUP [3000, 9900) with 2 kb flanks: left 2000, two copies of the
        // 6900 bp region, right flank [9900, 11900) → 100 bp.
        let hap =
            VariantHaplotype::from_tandem_duplication(&short_contig(), "chrEnd", 3000, 9900, 2000)
                .unwrap();
        // 2000 + 6900 + 6900 + 100 = 15900, not ... + 2000 = 17800.
        assert_eq!(hap.total_len, 15_900);
        assert_eq!(hap.ref_mapped_len(), 15_900);
        assert_eq!(hap.ref_range(), Some((1000, 10_000)));
    }

    #[test]
    fn test_inversion_at_chromosome_end_does_not_overcount_ref_len() {
        // INV [3000, 9900) with 2 kb flanks: left 2000, inverted 6900,
        // right flank [9900, 11900) → 100 bp.
        let hap =
            VariantHaplotype::from_inversion(&short_contig(), "chrEnd", 3000, 9900, 2000).unwrap();
        // 2000 + 6900 + 100 = 9000, not ... + 2000 = 10900.
        assert_eq!(hap.total_len, 9000);
        assert_eq!(hap.ref_mapped_len(), 9000);
        assert_eq!(hap.ref_range(), Some((1000, 10_000)));
    }

    #[test]
    fn test_insertion_at_chromosome_end_does_not_overcount_ref_len() {
        // INS at 9900 with 2 kb flanks: left [7900, 9900) = 2000 bp, 50 bp of
        // novel sequence, right [9900, 11900) → 100 bp.
        let hap =
            VariantHaplotype::from_insertion(&short_contig(), "chrEnd", 9900, &[b'G'; 50], 2000)
                .unwrap();
        assert_eq!(hap.total_len, 2150);
        // Reference-mapped bases only: 2000 + 100 = 2100, not 2000 + 2000 = 4000.
        assert_eq!(hap.ref_mapped_len(), 2100);
        assert_eq!(hap.ref_range(), Some((7900, 10_000)));
    }

    #[test]
    fn test_small_variant_at_chromosome_end_does_not_overcount_ref_len() {
        // SNP G→T at 9990 (the contig repeats ACGT, so 9990 is a G) with 2 kb
        // flanks: left [7990, 9990) = 2000 bp, the 1 bp alt, right
        // [9991, 11991) → [9991, 10000) = 9 bp.
        let hap =
            VariantHaplotype::from_small_variant(&short_contig(), "chrEnd", 9990, b"G", b"T", 2000)
                .unwrap();
        // 2000 + 1 + 9 = 2010, not 2000 + 1 + 2000 = 4001.
        assert_eq!(hap.total_len, 2010);
        assert_eq!(hap.ref_mapped_len(), 2010);
        assert_eq!(hap.ref_range(), Some((7990, 10_000)));
    }

    #[test]
    fn test_fusion_reverse_piece_at_chromosome_end_maps_to_real_bases() {
        // RightRight keeps revcomp(A[bp_a, bp_a + flank)). With bp_a 100 bp
        // from the contig end the fetch returns 100 bases, so haplotype offset
        // 0 is the contig's last base (9999) and offset 99 is bp_a (9900).
        // An unclamped ref_end would map them 1900 bp past the contig end.
        let hap = VariantHaplotype::from_fusion(
            &short_contig(), "chrEnd", 9900, "chrEnd", 2000, 2000, FusionJoin::RightRight,
        )
        .unwrap();
        assert_eq!(hap.segments[0].sequence.len(), 100);
        assert_eq!(ref_pos(&hap, 0), ("chrEnd".to_string(), 9999));
        assert_eq!(ref_pos(&hap, 99), ("chrEnd".to_string(), 9900));
    }

    #[test]
    fn test_a_tandem_dups_only_junction_is_between_its_copies() {
        // Its left flank runs straight into copy 1, and copy 2 into its right
        // flank; only copy 1 -> copy 2 jumps back.
        let hap = mock_tandem_dup(1000, 2000, 500);
        assert_eq!(hap.segments.len(), 4);
        assert_eq!(hap.junctions(), vec![hap.segments[2].hap_offset]);
    }

    #[test]
    fn test_a_deletion_has_one_junction_an_inversion_two_and_an_insertion_two() {
        let del = mock_segments_deletion(1000, 2000, 500);
        assert_eq!(del.junctions(), vec![del.segments[1].hap_offset]);
        let inv = mock_segments_inversion(1000, 2000, 500);
        assert_eq!(inv.junctions(), vec![inv.segments[1].hap_offset, inv.segments[2].hap_offset]);
        let ins = mock_segments_insertion(1000, b"ACGTACGT", 500);
        assert_eq!(ins.junctions(), vec![ins.segments[1].hap_offset, ins.segments[2].hap_offset]);
    }

    use crate::synth::PairSpans;

    /// A pair over `fragment` whose two mates cover `mates`.
    fn spans(fragment: (u64, u64), mates: [(u64, u64); 2]) -> PairSpans {
        PairSpans { fragment, mates }
    }

    #[test]
    fn test_a_read_shows_a_snv_only_with_a_base_on_each_side() {
        let snv = Evidence::Bases { start: 100, end: 101 };
        let far = (300, 450);
        // A read ending on the SNV, or starting on it, lacks a side.
        assert!(!snv.shown_by(&spans((0, 450), [(0, 101), far])));
        assert!(!snv.shown_by(&spans((100, 450), [(100, 250), far])));
        assert!(snv.shown_by(&spans((0, 450), [(0, 102), far])));
        assert!(snv.shown_by(&spans((99, 450), [(99, 249), far])));
        // Either mate will do; the fragment alone will not.
        assert!(snv.shown_by(&spans((0, 450), [far, (0, 102)])));
        assert!(!snv.shown_by(&spans((0, 450), [(0, 90), (300, 450)])));
    }

    #[test]
    fn test_a_read_shows_a_pure_deletion_only_across_its_join() {
        let join = Evidence::Bases { start: 100, end: 100 };
        let far = (300, 450);
        assert!(!join.shown_by(&spans((0, 450), [(0, 100), far])));
        assert!(!join.shown_by(&spans((100, 450), [(100, 250), far])));
        assert!(join.shown_by(&spans((0, 450), [(0, 101), far])));
        assert!(join.shown_by(&spans((99, 450), [(99, 249), far])));
    }

    #[test]
    fn test_a_fragment_shows_a_junction_even_when_neither_read_crosses_it() {
        let junction = Evidence::Junctions(vec![1000]);
        assert!(junction.shown_by(&spans((800, 1200), [(800, 950), (1050, 1200)])));
        // A fragment that ends at the junction, or starts there, is all one side.
        assert!(!junction.shown_by(&spans((600, 1000), [(600, 750), (850, 1000)])));
        assert!(!junction.shown_by(&spans((1000, 1400), [(1000, 1150), (1250, 1400)])));
        assert!(junction.shown_by(&spans((999, 1400), [(999, 1150), (1250, 1400)])));
        // Any one junction will do.
        let two = Evidence::Junctions(vec![1000, 3000]);
        assert!(two.shown_by(&spans((2800, 3200), [(2800, 2950), (3050, 3200)])));
    }
}
