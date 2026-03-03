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

/// Origin of a haplotype segment in the reference genome.
#[derive(Debug, Clone)]
pub struct SegmentOrigin {
    pub chrom: String,
    pub ref_start: u64,   // 0-based
    pub ref_end: u64,     // 0-based, exclusive
    pub is_reverse: bool, // true for reverse-complemented segments (INV)
}

/// A segment of the variant haplotype.
#[derive(Debug)]
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

        Ok(Self::from_segments(vec![
            HaplotypeSegment {
                sequence: left_seq,
                origin: Some(SegmentOrigin {
                    chrom: chrom.to_string(),
                    ref_start: left_start,
                    ref_end: del_start,
                    is_reverse: false,
                }),
                hap_offset: 0,
            },
            HaplotypeSegment {
                sequence: right_seq,
                origin: Some(SegmentOrigin {
                    chrom: chrom.to_string(),
                    ref_start: del_end,
                    ref_end: right_end,
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

        Ok(Self::from_segments(vec![
            HaplotypeSegment {
                sequence: left_seq,
                origin: Some(SegmentOrigin {
                    chrom: chrom.to_string(),
                    ref_start: left_start,
                    ref_end: dup_end,
                    is_reverse: false,
                }),
                hap_offset: 0,
            },
            HaplotypeSegment {
                sequence: right_seq,
                origin: Some(SegmentOrigin {
                    chrom: chrom.to_string(),
                    ref_start: dup_start,
                    ref_end: right_end,
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

        Ok(Self::from_segments(vec![
            HaplotypeSegment {
                sequence: left_seq,
                origin: Some(SegmentOrigin {
                    chrom: chrom.to_string(),
                    ref_start: left_start,
                    ref_end: dup_start,
                    is_reverse: false,
                }),
                hap_offset: 0,
            },
            HaplotypeSegment {
                sequence: dup_seq_1,
                origin: Some(SegmentOrigin {
                    chrom: chrom.to_string(),
                    ref_start: dup_start,
                    ref_end: dup_end,
                    is_reverse: false,
                }),
                hap_offset: 0,
            },
            HaplotypeSegment {
                sequence: dup_seq_2,
                origin: Some(SegmentOrigin {
                    chrom: chrom.to_string(),
                    ref_start: dup_start,
                    ref_end: dup_end,
                    is_reverse: false,
                }),
                hap_offset: 0,
            },
            HaplotypeSegment {
                sequence: right_seq,
                origin: Some(SegmentOrigin {
                    chrom: chrom.to_string(),
                    ref_start: dup_end,
                    ref_end: right_end,
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

        Ok(Self::from_segments(vec![
            HaplotypeSegment {
                sequence: left_seq,
                origin: Some(SegmentOrigin {
                    chrom: chrom.to_string(),
                    ref_start: left_start,
                    ref_end: inv_start,
                    is_reverse: false,
                }),
                hap_offset: 0,
            },
            HaplotypeSegment {
                sequence: inv_seq,
                origin: Some(SegmentOrigin {
                    chrom: chrom.to_string(),
                    ref_start: inv_start,
                    ref_end: inv_end,
                    is_reverse: true,
                }),
                hap_offset: 0,
            },
            HaplotypeSegment {
                sequence: right_seq,
                origin: Some(SegmentOrigin {
                    chrom: chrom.to_string(),
                    ref_start: inv_end,
                    ref_end: right_end,
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

        // Uppercase the insertion sequence for consistency.
        let ins_upper: Vec<u8> = ins_seq.iter().map(|b| b.to_ascii_uppercase()).collect();

        Ok(Self::from_segments(vec![
            HaplotypeSegment {
                sequence: left_seq,
                origin: Some(SegmentOrigin {
                    chrom: chrom.to_string(),
                    ref_start: left_start,
                    ref_end: pos,
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
                    ref_end: right_end,
                    is_reverse: false,
                }),
                hap_offset: 0,
            },
        ]))
    }

    /// Build a fusion haplotype.
    ///
    /// Normal: `ref_A[bp_a-flank..bp_a] | ref_B[bp_b..bp_b+flank]`
    /// Inverted: `ref_A[bp_a-flank..bp_a] | revcomp(ref_B[bp_b..bp_b+flank])`
    pub fn from_fusion(
        reference: &SharedReference,
        chrom_a: &str,
        bp_a: u64,
        chrom_b: &str,
        bp_b: u64,
        flank: u64,
        inverted: bool,
    ) -> Result<Self> {
        let left_start = bp_a.saturating_sub(flank);
        let right_end = bp_b.saturating_add(flank);

        let left_seq = fetch_upper(reference, chrom_a, left_start, bp_a)?;
        let mut right_seq = fetch_upper(reference, chrom_b, bp_b, right_end)?;

        if inverted {
            reverse_complement(&mut right_seq);
        }

        Ok(Self::from_segments(vec![
            HaplotypeSegment {
                sequence: left_seq,
                origin: Some(SegmentOrigin {
                    chrom: chrom_a.to_string(),
                    ref_start: left_start,
                    ref_end: bp_a,
                    is_reverse: false,
                }),
                hap_offset: 0,
            },
            HaplotypeSegment {
                sequence: right_seq,
                origin: Some(SegmentOrigin {
                    chrom: chrom_b.to_string(),
                    ref_start: bp_b,
                    ref_end: right_end,
                    is_reverse: inverted,
                }),
                hap_offset: 0,
            },
        ]))
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
                    ref_end: pos,
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
                    ref_end: right_end,
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

    /// Apply het SNP variants to the haplotype sequence.
    ///
    /// `variants` maps reference position → alternate allele base (uppercase).
    /// For each reference-origin segment, every base whose reference position
    /// appears in the map is replaced with the alt allele. This ensures that
    /// reads tiled across the haplotype carry the correct het SNP alleles
    /// instead of reference-only bases.
    pub fn apply_variants(&mut self, variants: &std::collections::HashMap<u64, u8>) {
        if variants.is_empty() {
            return;
        }
        let mut applied = 0usize;
        for seg in &mut self.segments {
            let origin = match &seg.origin {
                Some(o) => o,
                None => continue, // novel insertion — no reference positions
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
    /// (i.e., is not novel insertion sequence). Used to reject tiling fragments that
    /// would land entirely in novel sequence and not contribute to observable coverage.
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
}
