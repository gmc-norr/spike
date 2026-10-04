//! A small variant the sample already carries (review finding 3).
//!
//! spike edits one copy of the sample and keeps the other copy's reads as
//! they were. When the sample already has another allele at the bases a small
//! variant changes, those reads still show it, so no edit gives the requested
//! fraction: a hom-alt site asked for at 0.5 came out all ALT, while its truth
//! record said 0.5. So the reads at each small variant's own bases are counted
//! first, and the event is refused when enough of them carry another allele.
//! The rule was measured on the hospital BAM before it was built
//! (`docs/superpowers/plans/2026-10-04-carried-allele.md`).

use std::io;

use anyhow::{Context, Result};
use noodles::sam::alignment::record::cigar::op::{Kind, Op};
use noodles::sam::alignment::record::Cigar as CigarTrait;
use noodles::sam::alignment::record::Flags;

/// Fewer reads than this spanning a site and it is not judged: the pileup's
/// own floor (`loh::call_snps`).
pub const MIN_READS: u32 = 10;

/// The share of the spanning reads that must carry another allele: the
/// pileup's lower bound for a het allele (`loh::call_snps`).
pub const MIN_SHARE: f64 = 0.2;

/// The reference bases a small variant changes, `[start, end)` 0-based: REF
/// and ALT without their common prefix, then their common suffix. A pure
/// insertion has `start == end`; its bases go between `start - 1` and `start`.
pub fn changed_span(pos: u64, ref_allele: &[u8], alt_allele: &[u8]) -> (u64, u64) {
    let prefix = ref_allele.iter().zip(alt_allele).take_while(|(r, a)| r == a).count();
    let (r, a) = (&ref_allele[prefix..], &alt_allele[prefix..]);
    let suffix = r.iter().rev().zip(a.iter().rev()).take_while(|(x, y)| x == y).count();
    let start = pos + prefix as u64;
    (start, start + (r.len() - suffix) as u64)
}

/// What one read shows at a site.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ReadAtSite {
    /// It does not span bases `start - 1` through `end`.
    Elsewhere,
    /// It spans them and shows the reference.
    Reference,
    /// It spans them and shows another allele: a base other than the
    /// reference inside the site, a deletion over it, or an insertion at one
    /// of its boundaries.
    Other,
}

/// What a read aligned at `align_start` (0-based) shows at `[start, end)`,
/// whose reference bases are `site_ref`.
pub fn read_at_site<I>(align_start: u64, seq: &[u8], cigar: I, start: u64, end: u64, site_ref: &[u8]) -> ReadAtSite
where
    I: IntoIterator<Item = io::Result<Op>>,
{
    let mut ref_pos = align_start;
    let mut seq_pos = 0usize;
    let mut other = false;
    for op in cigar {
        let Ok(op) = op else { break };
        let len = op.len();
        match op.kind() {
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                for p in ref_pos.max(start)..(ref_pos + len as u64).min(end) {
                    let base = seq.get(seq_pos + (p - ref_pos) as usize).map(u8::to_ascii_uppercase);
                    let want = site_ref.get((p - start) as usize).map(u8::to_ascii_uppercase);
                    if matches!(base, Some(b) if b != b'N' && Some(b) != want) {
                        other = true;
                    }
                }
                ref_pos += len as u64;
                seq_pos += len;
            }
            // Its bases go between `ref_pos - 1` and `ref_pos`.
            Kind::Insertion => {
                if start <= ref_pos && ref_pos <= end {
                    other = true;
                }
                seq_pos += len;
            }
            Kind::SoftClip => seq_pos += len,
            Kind::Deletion | Kind::Skip => {
                if ref_pos < end.max(start + 1) && ref_pos + len as u64 > start {
                    other = true;
                }
                ref_pos += len as u64;
            }
            Kind::HardClip | Kind::Pad => {}
        }
    }
    // The read spans [align_start, ref_pos); it must hold start - 1 through end.
    match (align_start < start && ref_pos > end, other) {
        (false, _) => ReadAtSite::Elsewhere,
        (true, false) => ReadAtSite::Reference,
        (true, true) => ReadAtSite::Other,
    }
}

/// The reads counted at one site.
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
pub struct SiteCount {
    /// Reads that span the site.
    pub spanning: u32,
    /// Of those, the reads that carry another allele.
    pub other: u32,
}

impl SiteCount {
    fn add(&mut self, read: ReadAtSite) {
        match read {
            ReadAtSite::Elsewhere => {}
            ReadAtSite::Reference => self.spanning += 1,
            ReadAtSite::Other => {
                self.spanning += 1;
                self.other += 1;
            }
        }
    }

    /// Whether the sample carries another allele here; `None` when fewer than
    /// [`MIN_READS`] reads span the site.
    pub fn carried(&self) -> Option<bool> {
        (self.spanning >= MIN_READS).then(|| self.other as f64 / self.spanning as f64 >= MIN_SHARE)
    }
}

/// The pileup's filter (`loh::count_alleles`).
fn counted(flags: Flags, mapq: Option<noodles::sam::alignment::record::MappingQuality>, min_mapq: u8) -> bool {
    !(flags.is_unmapped()
        || flags.is_secondary()
        || flags.is_supplementary()
        || flags.is_duplicate()
        || flags.is_qc_fail())
        && mapq.map_or(0, u8::from) >= min_mapq
}

/// Counts reads at sites of one BAM or CRAM, with the pileup's filter:
/// primary, not duplicate, not QC-fail, MAPQ at least `min_mapq`.
pub struct SiteReader {
    source: Source,
}

enum Source {
    Bam(Box<BamSource>),
    Cram {
        path: String,
        repository: noodles::fasta::Repository,
    },
}

struct BamSource {
    reader: noodles::bam::io::IndexedReader<noodles::bgzf::Reader<std::fs::File>>,
    header: noodles::sam::Header,
}

impl SiteReader {
    pub fn open(alignment_path: &str, ref_path: &str) -> Result<Self> {
        let source = if crate::extract::is_cram(alignment_path) {
            Source::Cram {
                path: alignment_path.to_string(),
                repository: crate::extract::build_fasta_repository(ref_path)?,
            }
        } else {
            let mut reader = noodles::bam::io::indexed_reader::Builder::default()
                .build_from_path(alignment_path)
                .with_context(|| format!("failed to open BAM for the carried-allele check: {}", alignment_path))?;
            let header = reader.read_header()?;
            Source::Bam(Box::new(BamSource { reader, header }))
        };
        Ok(SiteReader { source })
    }

    /// The reads at `[start, end)` on `chrom`, whose reference bases are `site_ref`.
    pub fn count(&mut self, chrom: &str, start: u64, end: u64, site_ref: &[u8], min_mapq: u8) -> Result<SiteCount> {
        let region = noodles::core::Region::new(
            chrom,
            crate::extract::safe_noodles_position(start.max(1))..=crate::extract::safe_noodles_position(end + 1),
        );
        let mut count = SiteCount::default();
        match &mut self.source {
            Source::Bam(bam) => {
                let BamSource { reader, header } = &mut **bam;
                for result in reader.query(header, &region)? {
                    let record = result?;
                    if !counted(record.flags(), record.mapping_quality(), min_mapq) {
                        continue;
                    }
                    let Some(Ok(p)) = record.alignment_start() else { continue };
                    let seq: Vec<u8> = record.sequence().iter().collect();
                    let cigar = record.cigar();
                    count.add(read_at_site(usize::from(p) as u64 - 1, &seq, cigar.iter(), start, end, site_ref));
                }
            }
            Source::Cram { path, repository } => {
                let (mut reader, header) = crate::extract::open_cram_reader_for_region(path, repository, &region)
                    .with_context(|| format!("failed to open CRAM for the carried-allele check: {}", path))?;
                let queried = header.reference_sequences().get_index_of(chrom.as_bytes());
                for result in reader.query(&header, &region)? {
                    let buf = result?.try_into_alignment_record(&header)?;
                    // A container holding several contigs is decoded whole (L2, N4).
                    if !crate::extract::record_is_on_queried_reference(&buf, queried) {
                        continue;
                    }
                    if !counted(buf.flags(), buf.mapping_quality(), min_mapq) {
                        continue;
                    }
                    let Some(p) = buf.alignment_start() else { continue };
                    let cigar = buf.cigar();
                    count.add(read_at_site(
                        usize::from(p) as u64 - 1,
                        buf.sequence().as_ref(),
                        CigarTrait::iter(&cigar),
                        start,
                        end,
                        site_ref,
                    ));
                }
            }
        }
        Ok(count)
    }
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;
    use noodles::sam::alignment::record::cigar::op::Kind;

    //                     0123456789012345678901234567890
    pub(crate) const REF: &[u8] = b"ACGTACGTACGTACGTACGTACGTACGTACGT";

    fn ops(cigar: &[(Kind, usize)]) -> Vec<io::Result<Op>> {
        cigar.iter().map(|&(k, n)| Ok(Op::new(k, n))).collect()
    }

    fn at(start: u64, seq: &[u8], cigar: &[(Kind, usize)], site: (u64, u64)) -> ReadAtSite {
        let site_ref = &REF[site.0 as usize..site.1 as usize];
        read_at_site(start, seq, ops(cigar), site.0, site.1, site_ref)
    }

    fn with(base: u8, at: usize) -> Vec<u8> {
        let mut seq = REF[..20].to_vec();
        seq[at] = base;
        seq
    }

    #[test]
    fn test_changed_span_keeps_only_the_changed_bases() {
        assert_eq!(changed_span(10, b"G", b"T"), (10, 11), "SNV");
        assert_eq!(changed_span(10, b"GT", b"CA"), (10, 12), "MNV");
        assert_eq!(changed_span(9, b"CGT", b"C"), (10, 12), "deletion: the anchor goes");
        assert_eq!(changed_span(9, b"C", b"CAA"), (10, 10), "insertion: between 9 and 10");
        assert_eq!(changed_span(9, b"CGTA", b"CTA"), (10, 11), "a shared suffix goes too");
    }

    #[test]
    fn test_a_read_showing_the_reference_spans_and_carries_nothing() {
        assert_eq!(at(0, &REF[..20], &[(Kind::Match, 20)], (10, 11)), ReadAtSite::Reference);
    }

    #[test]
    fn test_a_mismatch_counts_inside_the_site_only() {
        assert_ne!(REF[10], b'T');
        assert_eq!(at(0, &with(b'T', 10), &[(Kind::Match, 20)], (10, 11)), ReadAtSite::Other);
        assert_ne!(REF[12], b'C');
        assert_eq!(at(0, &with(b'C', 12), &[(Kind::Match, 20)], (10, 11)), ReadAtSite::Reference);
        assert_ne!(REF[9], b'A');
        assert_eq!(at(0, &with(b'A', 9), &[(Kind::Match, 20)], (10, 11)), ReadAtSite::Reference);
    }

    #[test]
    fn test_an_n_base_is_not_another_allele() {
        assert_eq!(at(0, &with(b'N', 10), &[(Kind::Match, 20)], (10, 11)), ReadAtSite::Reference);
    }

    #[test]
    fn test_a_deletion_over_the_site_is_another_allele() {
        let seq = [&REF[..10], &REF[12..20]].concat();
        assert_eq!(at(0, &seq, &[(Kind::Match, 10), (Kind::Deletion, 2), (Kind::Match, 8)], (10, 12)), ReadAtSite::Other);
        // One that ends before the site is not.
        let seq = [&REF[..7], &REF[9..20]].concat();
        assert_eq!(at(0, &seq, &[(Kind::Match, 7), (Kind::Deletion, 2), (Kind::Match, 11)], (10, 12)), ReadAtSite::Reference);
        // At an insertion site, a deletion of the base after it is another allele.
        let seq = [&REF[..10], &REF[11..20]].concat();
        assert_eq!(at(0, &seq, &[(Kind::Match, 10), (Kind::Deletion, 1), (Kind::Match, 9)], (10, 10)), ReadAtSite::Other);
    }

    #[test]
    fn test_an_insertion_is_another_allele_at_the_site_boundaries_only() {
        for (m, want) in [(10, ReadAtSite::Other), (11, ReadAtSite::Other), (12, ReadAtSite::Reference)] {
            let seq = [&REF[..m], b"AA", &REF[m..18]].concat();
            let cigar = [(Kind::Match, m), (Kind::Insertion, 2), (Kind::Match, 18 - m)];
            assert_eq!(at(0, &seq, &cigar, (10, 11)), want, "insertion at boundary {m}");
        }
    }

    #[test]
    fn test_a_read_that_stops_inside_the_site_does_not_span_it() {
        assert_eq!(at(0, &REF[..11], &[(Kind::Match, 11)], (10, 12)), ReadAtSite::Elsewhere);
        assert_eq!(at(11, &REF[11..21], &[(Kind::Match, 10)], (10, 12)), ReadAtSite::Elsewhere);
        // Bases s - 1 and e are needed too.
        assert_eq!(at(10, &REF[10..20], &[(Kind::Match, 10)], (10, 11)), ReadAtSite::Elsewhere);
        assert_eq!(at(0, &REF[..11], &[(Kind::Match, 11)], (10, 11)), ReadAtSite::Elsewhere);
        assert_eq!(at(0, &REF[..12], &[(Kind::Match, 12)], (10, 11)), ReadAtSite::Reference);
    }

    #[test]
    fn test_a_site_is_judged_from_ten_reads_and_a_fifth_of_them() {
        let count = |spanning, other| SiteCount { spanning, other }.carried();
        assert_eq!(count(9, 9), None);
        assert_eq!(count(10, 2), Some(true));
        assert_eq!(count(10, 1), Some(false));
        assert_eq!(count(20, 4), Some(true));
        assert_eq!(count(20, 3), Some(false));
    }

    /// chrT, 2 kb of `REF` repeated, with `total` reads of 100M starting at
    /// 1,450, 1,451 and so on (at most 48, so each spans 1,499-1,501). The
    /// first `alt_reads` carry T at 1,500, where the reference has none. Then
    /// two more carrying T there that the pileup does not count: a duplicate,
    /// and one at MAPQ 10.
    pub(crate) fn site_bam(dir: &std::path::Path, alt_reads: usize, total: usize) -> String {
        use noodles::sam::alignment::record::{Flags, MappingQuality};
        use noodles::sam::alignment::record_buf::Sequence;
        use noodles::sam::alignment::RecordBuf;
        let contig: Vec<u8> = (0..2_000).map(|i| REF[i % REF.len()]).collect();
        let mut records = Vec::new();
        for i in 0..total + 2 {
            let start = 1_450 + i;
            let mut seq = contig[start..start + 100].to_vec();
            if i < alt_reads || i >= total {
                seq[1_500 - start] = b'T';
            }
            let (flags, mapq) = match i.checked_sub(total) {
                Some(0) => (0x400, 60),
                Some(_) => (0, 10),
                None => (0, 60),
            };
            records.push(
                RecordBuf::builder()
                    .set_name(format!("r{i}"))
                    .set_flags(Flags::from(flags))
                    .set_reference_sequence_id(0)
                    .set_alignment_start(noodles::core::Position::new(start + 1).unwrap())
                    .set_mapping_quality(MappingQuality::new(mapq).unwrap())
                    .set_cigar([Op::new(Kind::Match, 100)].into_iter().collect())
                    .set_sequence(Sequence::from(seq))
                    .build(),
            );
        }
        crate::extract::test_fixtures::write_one_contig_bam(&dir.join(format!("site_{alt_reads}_{total}.bam")), "chrT", 2_000, &records)
    }

    pub(crate) fn site_dir(name: &str) -> std::path::PathBuf {
        let dir = std::env::temp_dir().join(format!("spike_carried_{}_{}", name, std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        dir
    }

    fn count_at_1500(alt_reads: usize, total: usize) -> SiteCount {
        let dir = site_dir(&format!("count_{alt_reads}_{total}"));
        let bam = site_bam(&dir, alt_reads, total);
        let site_ref = [REF[1_500 % REF.len()]];
        assert_ne!(site_ref[0], b'T');
        SiteReader::open(&bam, "unused.fa").unwrap().count("chrT", 1_500, 1_501, &site_ref, 20).unwrap()
    }

    #[test]
    fn test_a_bam_site_counts_the_pileup_s_reads_only() {
        // The duplicate and the MAPQ-10 read both carry T, and neither counts.
        assert_eq!(count_at_1500(0, 30), SiteCount { spanning: 30, other: 0 });
        assert_eq!(count_at_1500(30, 30), SiteCount { spanning: 30, other: 30 });
        assert_eq!(count_at_1500(15, 30), SiteCount { spanning: 30, other: 15 });
        assert_eq!(count_at_1500(9, 9).carried(), None);
    }
}
