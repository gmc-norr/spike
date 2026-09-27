//! CR4's resistant-read census: how much of the depth over an event spike
//! cannot edit.
//!
//! Only proper pairs whose reads both pass the MAPQ and flag filters enter an
//! event's donor pool, and the pool is everything the event can suppress. Any
//! other read over the event -- below `--min-mapq`, not a proper pair, a mate
//! unmapped or failing a filter -- stays in the merged BAM as it was. This
//! counts them. It changes nothing spike emits.

use anyhow::Result;
use std::collections::HashSet;

use crate::types::{DepthFold, SimEvent};

/// Above this resistant share spike warns: the reads it cannot touch are
/// more than a tenth of the event's depth, so what it realises is off the
/// request by more than a tenth. Locked in the CR4 plan before any share was
/// measured.
pub const WARN_ABOVE: f64 = 0.10;

/// The reads over one event: how many were counted, and how many of those
/// spike could not edit.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct Census {
    pub counted: usize,
    pub resistant: usize,
}

impl Census {
    /// The resistant share, 0 when no read was counted.
    pub fn fraction(&self) -> f64 {
        if self.counted == 0 {
            0.0
        } else {
            self.resistant as f64 / self.counted as f64
        }
    }
}

/// The reference spans whose reads the census counts: the span a DEL, DUP
/// or INV changes, the two bases around an insertion point, a small
/// variant's REF, and the two bases around each of a fusion's cuts.
pub fn census_spans(event: &SimEvent) -> Vec<(String, u64, u64)> {
    let around = |chrom: &str, point: u64| (chrom.to_string(), point.saturating_sub(1), point + 1);
    match event {
        SimEvent::Deletion {
            chrom,
            del_start: start,
            del_end: end,
            ..
        }
        | SimEvent::Duplication {
            chrom,
            dup_start: start,
            dup_end: end,
            ..
        }
        | SimEvent::Inversion {
            chrom,
            inv_start: start,
            inv_end: end,
            ..
        } => vec![(chrom.clone(), *start, *end)],
        SimEvent::Insertion { chrom, pos, .. } => vec![around(chrom, *pos)],
        SimEvent::SmallVariant {
            chrom,
            pos,
            ref_allele,
            ..
        } => vec![(chrom.clone(), *pos, *pos + ref_allele.len() as u64)],
        SimEvent::Fusion {
            chrom_a,
            bp_a,
            chrom_b,
            bp_b,
            ..
        } => vec![around(chrom_a, *bp_a), around(chrom_b, *bp_b)],
    }
}

/// Count the primary, mapped, non-duplicate, non-QC-fail records over
/// `spans`, at any MAPQ; a record whose read name is not in `editable` is
/// resistant. A BAM span is read in [`crate::extract::read_chunks`] chunks on
/// the thread pool.
pub fn count_resistant(
    alignment_path: &str,
    ref_path: &str,
    spans: &[(String, u64, u64)],
    editable: &HashSet<String>,
) -> Result<Census> {
    count_resistant_with(alignment_path, ref_path, spans, editable, &crate::extract::read_chunks)
}

/// [`count_resistant`] with every BAM span read in `n_chunks` chunks.
#[cfg(test)]
fn count_resistant_in(
    alignment_path: &str,
    ref_path: &str,
    spans: &[(String, u64, u64)],
    editable: &HashSet<String>,
    n_chunks: usize,
) -> Result<Census> {
    count_resistant_with(alignment_path, ref_path, spans, editable, &|_, _| n_chunks)
}

/// [`count_resistant`], with `chunks(start, end)` the number of chunks a BAM
/// span is read in (see [`crate::extract::fold_bam_region`]). Counts add, so
/// the chunks' sum is one query's count.
fn count_resistant_with(
    alignment_path: &str,
    ref_path: &str,
    spans: &[(String, u64, u64)],
    editable: &HashSet<String>,
    chunks: &(dyn Fn(u64, u64) -> usize + Sync),
) -> Result<Census> {
    let mut census = Census::default();
    for (chrom, start, end) in spans {
        if crate::extract::is_cram(alignment_path) {
            // MAPQ 0: the census counts the reads the pool's MAPQ filter left out.
            crate::validate::for_each_alignment(
                alignment_path,
                ref_path,
                chrom,
                *start,
                *end,
                0,
                &mut |name, _, _, _| {
                    census.counted += 1;
                    if !editable.contains(String::from_utf8_lossy(name).as_ref()) {
                        census.resistant += 1;
                    }
                },
            )?;
            continue;
        }
        let per_chunk = crate::extract::fold_bam_region(
            alignment_path,
            chrom,
            *start,
            *end,
            chunks(*start, *end),
            "failed to open BAM for a read scan",
            Census::default,
            |chunk: &mut Census, record| {
                // The records `for_each_alignment` hands on, at MAPQ 0.
                if !crate::validate::usable_alignment(record.flags(), record.mapping_quality().map(u8::from), 0) {
                    return Ok(());
                }
                let Some(name) = record.name() else {
                    return Ok(());
                };
                if !matches!(record.alignment_start(), Some(Ok(_))) {
                    return Ok(());
                }
                // It also parses each such record's CIGAR, and a bad one stops the scan.
                record.cigar().iter().collect::<std::io::Result<Vec<_>>>()?;
                let name: &[u8] = name.as_ref();
                chunk.counted += 1;
                if !editable.contains(String::from_utf8_lossy(name).as_ref()) {
                    chunk.resistant += 1;
                }
                Ok(())
            },
        )?;
        for chunk in per_chunk {
            census.counted += chunk.counted;
            census.resistant += chunk.resistant;
        }
    }
    Ok(census)
}

/// Above this depth fold spike warns (CR2): for a het DUP a bin at a third
/// less than the scaling depth comes out about a third too deep. Locked in the
/// CR2 plan before any fold was measured.
pub const DEPTH_FOLD_WARN_ABOVE: f64 = 1.5;

/// The warning for an event whose depth fold is above [`DEPTH_FOLD_WARN_ABOVE`].
pub fn depth_fold_warning(label: &str, fold: &DepthFold, edit_origin: bool) -> Option<String> {
    if fold.fold <= DEPTH_FOLD_WARN_ABOVE {
        return None;
    }
    Some(format!(
        "{}: the donor's depth over {} is {:.1}x, but every fragment this event tiles is \
         scaled by the {:.1}x measured at one of its breakpoints ({:.2}-fold). Where the \
         donor's depth differs from that, the event's depth there is wrong by about that \
         much; truth.vcf records the fold as SIM_DEPTH_FOLD (CR2).{}",
        label,
        fold.worst_bin,
        fold.worst_depth,
        fold.scaled_by,
        fold.fold,
        origin_suffix(edit_origin),
    ))
}

/// Above this resistant share spike refuses the event (RF8) unless told
/// `--allow-resistant`. Locked in the RF8 plan before any chr1 share was
/// measured.
pub const REFUSE_ABOVE: f64 = 0.5;

/// One line of [`refusal_message`] for an event whose resistant share is
/// above [`REFUSE_ABOVE`], unless `allow` is set.
///
/// The resistant reads are never removed or replaced, so the event that
/// reaches the reads is about `VAF x (1 - share)`: above one half they carry
/// less than half of what the truth record would claim (RF8).
pub fn refusal(label: &str, census: &Census, allow: bool) -> Option<String> {
    if allow || census.fraction() <= REFUSE_ABOVE {
        return None;
    }
    Some(format!(
        "{}: {} of {} reads over it ({:.1}%) are ones spike cannot edit",
        label,
        census.resistant,
        census.counted,
        census.fraction() * 100.0,
    ))
}

/// The error for a run whose events include any [`refusal`]: every refused
/// event at once, so a multi-event input can drop them all in one pass.
pub fn refusal_message(refusals: &[String], edit_origin: bool) -> Option<String> {
    if refusals.is_empty() {
        return None;
    }
    let events = if refusals.len() == 1 {
        "1 event".to_string()
    } else {
        format!("{} events", refusals.len())
    };
    Some(format!(
        "The reads over {} would carry less than half of what truth.vcf would claim, \
         so spike stops rather than write it (RF8):\n  {}\nThose reads are {}, and they stay in the merged \
         BAM as they are. If most of them are below --min-mapq, lowering it lets spike \
         edit them. Otherwise remove the events from the input, or pass \
         --allow-resistant to simulate them anyway; truth.vcf then records the share \
         as SIM_RESIST.{}",
        events,
        refusals.join("\n  "),
        if edit_origin {
            UNEDITED_UNDER_ORIGIN
        } else {
            "below --min-mapq, not a proper pair, or have a mate that fails a filter"
        },
        origin_suffix(edit_origin),
    ))
}

/// Said wherever the resistant share or the depth fold is reported under
/// `--edit-model origin` (PD-29, S6). Every threshold here was locked for the
/// `clean` measures -- the share of reads outside the donor pool, a fold
/// between pool depths -- and none has been measured for origin's, whose share
/// counts a read origin may remove as editable and whose fold compares origin
/// depths.
/// The reads spike cannot edit under `--edit-model origin`: R4 keeps a pair
/// with a mate that could not have come from the event's footprint.
const UNEDITED_UNDER_ORIGIN: &str = "neither in the donor pool nor in a pair origin may remove";

/// [`origin_thresholds_note`] after a space under origin, nothing under clean.
fn origin_suffix(edit_origin: bool) -> String {
    if edit_origin {
        format!(" {}", origin_thresholds_note())
    } else {
        String::new()
    }
}

pub fn origin_thresholds_note() -> String {
    format!(
        "Under --edit-model origin a read origin may remove counts as editable, and the \
         depth fold compares origin depths; the thresholds (warn above a {} share, refuse \
         above {}, warn above a {}-fold) were set for --edit-model clean and are not yet \
         checked for origin.",
        WARN_ABOVE, REFUSE_ABOVE, DEPTH_FOLD_WARN_ABOVE,
    )
}

/// The warning for an event whose resistant share is above [`WARN_ABOVE`].
pub fn warning(label: &str, census: &Census, edit_origin: bool) -> Option<String> {
    if census.fraction() <= WARN_ABOVE {
        return None;
    }
    Some(format!(
        "{}: {} of {} reads over it ({:.0}%) are ones spike cannot edit ({}). They stay in the \
         merged BAM as they are, so the event is weaker than requested; truth.vcf records \
         the share as SIM_RESIST (CR4).{}",
        label,
        census.resistant,
        census.counted,
        census.fraction() * 100.0,
        if edit_origin {
            UNEDITED_UNDER_ORIGIN
        } else {
            "below --min-mapq, not a proper pair, or a mate that fails a filter"
        },
        origin_suffix(edit_origin),
    ))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::FusionJoin;

    fn del(start: u64, end: u64) -> SimEvent {
        SimEvent::Deletion {
            chrom: "chr1".to_string(),
            del_start: start,
            del_end: end,
            gene: "G".to_string(),
            exons: vec![],
            allele_fraction: None,
        }
    }

    #[test]
    fn test_census_spans_cover_what_each_event_changes() {
        let span = |c: &str, s: u64, e: u64| (c.to_string(), s, e);
        assert_eq!(census_spans(&del(100, 200)), [span("chr1", 100, 200)]);
        let dup = SimEvent::Duplication {
            chrom: "chr2".to_string(),
            dup_start: 10,
            dup_end: 50,
            gene: "G".to_string(),
            allele_fraction: None,
        };
        assert_eq!(census_spans(&dup), [span("chr2", 10, 50)]);
        let inv = SimEvent::Inversion {
            chrom: "chr3".to_string(),
            inv_start: 5,
            inv_end: 9,
            gene: "G".to_string(),
            allele_fraction: None,
        };
        assert_eq!(census_spans(&inv), [span("chr3", 5, 9)]);
        let ins = SimEvent::Insertion {
            chrom: "chr1".to_string(),
            pos: 300,
            ins_seq: None,
            ins_len: 20,
            gene: "G".to_string(),
            allele_fraction: None,
        };
        assert_eq!(census_spans(&ins), [span("chr1", 299, 301)]);
        let snv = SimEvent::SmallVariant {
            chrom: "chr1".to_string(),
            pos: 40,
            ref_allele: b"ACG".to_vec(),
            alt_allele: b"A".to_vec(),
            gene: "G".to_string(),
            allele_fraction: None,
        };
        assert_eq!(census_spans(&snv), [span("chr1", 40, 43)]);
        let fusion = SimEvent::Fusion {
            chrom_a: "chr1".to_string(),
            bp_a: 1000,
            gene_a: "A".to_string(),
            chrom_b: "chr2".to_string(),
            bp_b: 5000,
            gene_b: "B".to_string(),
            allele_fraction: None,
            join: FusionJoin::Forward,
        };
        assert_eq!(
            census_spans(&fusion),
            [span("chr1", 999, 1001), span("chr2", 4999, 5001)]
        );
    }

    #[test]
    fn test_census_fraction_is_the_resistant_share() {
        assert_eq!(Census { counted: 0, resistant: 0 }.fraction(), 0.0);
        assert_eq!(Census { counted: 10, resistant: 5 }.fraction(), 0.5);
    }

    #[test]
    fn test_depth_fold_warns_only_above_the_threshold() {
        let at = |fold| {
            depth_fold_warning(
                "DUP chrT:10001-28000",
                &DepthFold {
                    fold,
                    scaled_by: 100.0,
                    worst_bin: "chrT:20000-21000".to_string(),
                    worst_depth: 25.0,
                },
                false,
            )
        };
        assert_eq!(at(1.5), None, "1.5 is not above 1.5");
        assert_eq!(at(1.0), None);
        let w = at(3.88).expect("3.88 is above 1.5");
        assert!(w.contains("DUP chrT:10001-28000"), "{}", w);
        assert!(w.contains("chrT:20000-21000"), "{}", w);
        assert!(w.contains("25.0x"), "{}", w);
        assert!(w.contains("100.0x"), "{}", w);
        assert!(w.contains("SIM_DEPTH_FOLD"), "{}", w);
    }

    #[test]
    fn test_census_warns_only_above_the_threshold() {
        let at = |counted, resistant| warning("DEL chr1:101-200", &Census { counted, resistant }, false);
        assert_eq!(at(100, 10), None, "0.10 is not above 0.10");
        assert_eq!(at(100, 0), None);
        let w = at(100, 11).expect("0.11 is above 0.10");
        assert!(w.contains("DEL chr1:101-200"), "{}", w);
        assert!(w.contains("11 of 100"), "{}", w);
        assert!(w.contains("--min-mapq"), "{}", w);
    }

    #[test]
    fn test_census_refuses_only_above_half() {
        let at = |counted, resistant| {
            refusal("DEL chr1:101-200", &Census { counted, resistant }, false)
        };
        assert_eq!(at(100, 50), None, "0.5 is not above 0.5");
        assert_eq!(at(100, 11), None, "above WARN_ABOVE is a warning, not a refusal");
        assert_eq!(at(0, 0), None);
        // CR4's lowmap probe: 1038 of 2074 is 0.5005.
        let r = at(2074, 1038).expect("0.5005 is above 0.5");
        assert!(r.contains("DEL chr1:101-200"), "{}", r);
        assert!(r.contains("1038 of 2074"), "{}", r);
    }

    #[test]
    fn test_allow_resistant_lifts_the_refusal() {
        // RF8's locus, del:chr20:27100000-27110000 on the 35x HG002 BAM.
        let census = Census { counted: 5749, resistant: 5731 };
        assert!(refusal("DEL chr20:27100001-27110000", &census, false).is_some());
        assert_eq!(refusal("DEL chr20:27100001-27110000", &census, true), None);
    }

    #[test]
    fn test_under_origin_the_messages_say_what_origin_counts_and_that_thresholds_are_unchecked() {
        // PD-29, S6: under origin a read in a pair origin may remove is
        // editable, so the clean definition in these messages is wrong there,
        // and every threshold was set for clean. Under clean nothing changes.
        let census = Census { counted: 100, resistant: 60 };
        let fold = DepthFold { fold: 3.88, scaled_by: 100.0, worst_bin: "chrT:20000-21000".to_string(), worst_depth: 25.0 };
        let refused = || vec![refusal("DEL chr1:101-200", &census, false).unwrap()];
        let note = origin_thresholds_note();
        for m in [
            warning("DEL chr1:101-200", &census, true).unwrap(),
            refusal_message(&refused(), true).unwrap(),
        ] {
            assert!(m.contains("neither in the donor pool nor in a pair origin may remove"), "{}", m);
            assert!(!m.contains("not a proper pair"), "{}", m);
            assert!(m.ends_with(&note), "{}", m);
        }
        assert!(depth_fold_warning("DUP chrT:10001-28000", &fold, true).unwrap().ends_with(&note));
        for m in [
            warning("DEL chr1:101-200", &census, false).unwrap(),
            refusal_message(&refused(), false).unwrap(),
            depth_fold_warning("DUP chrT:10001-28000", &fold, false).unwrap(),
        ] {
            assert!(!m.contains("--edit-model"), "{}", m);
        }
        assert!(warning("DEL chr1:101-200", &census, false)
            .unwrap()
            .contains("(below --min-mapq, not a proper pair, or a mate that fails a filter)"));
        for number in ["0.1", "0.5", "1.5"] {
            assert!(note.contains(number), "{}", note);
        }
    }

    #[test]
    fn test_refusal_message_names_every_refused_event() {
        assert_eq!(refusal_message(&[], false), None, "no refused event, no error");
        // Two of validate_pipeline.sh's events on its NA18488 background.
        let at = |label, counted, resistant| {
            refusal(label, &Census { counted, resistant }, false)
        };
        let first = at("DEL chr20:32723064-32724971", 407, 316).expect("0.776 is above 0.5");
        let second = at("DEL chr20:61943514-61945040", 200, 189).expect("0.945 is above 0.5");
        let m = refusal_message(&[first, second], false).expect("two refused events");
        assert!(m.contains("2 events"), "{}", m);
        assert!(m.contains("DEL chr20:32723064-32724971"), "{}", m);
        assert!(m.contains("316 of 407"), "{}", m);
        assert!(m.contains("DEL chr20:61943514-61945040"), "{}", m);
        assert!(m.contains("189 of 200"), "{}", m);
        assert!(m.contains("--allow-resistant"), "{}", m);
        assert!(m.contains("--min-mapq"), "{}", m);
        assert!(m.contains("SIM_RESIST"), "{}", m);
    }

    /// A one-contig BAM (`chrA`, 3000 bp) of 100 bp proper pairs, read 2
    /// 200 bp after read 1, each pair given as (name, 1-based start, MAPQ,
    /// extra flag bits), plus a hand-written `.bai` (see
    /// `extract::tests::write_three_pair_bam` for why one bin is enough).
    /// Pairs must be given in start order.
    fn write_pairs_bam(dir: &std::path::Path, pairs: &[(&str, usize, u8, u16)]) -> String {
        use noodles::sam::alignment::io::Write as _;
        use std::num::NonZeroUsize;

        const READ_LEN: usize = 100;
        let header = noodles::sam::Header::builder()
            .add_reference_sequence(
                "chrA",
                noodles::sam::header::record::value::Map::<
                    noodles::sam::header::record::value::map::ReferenceSequence,
                >::new(NonZeroUsize::try_from(3000).unwrap()),
            )
            .build();
        let record = |name: &str, start: usize, mapq: u8, extra: u16, first: bool| {
            let (pos, mate_pos) = if first { (start, start + 200) } else { (start + 200, start) };
            let span = 200 + READ_LEN;
            noodles::sam::alignment::RecordBuf::builder()
                .set_name(name)
                .set_flags(noodles::sam::alignment::record::Flags::from(
                    if first { 0x63u16 } else { 0x93u16 } | extra,
                ))
                .set_reference_sequence_id(0)
                .set_alignment_start(noodles::core::Position::new(pos).unwrap())
                .set_mapping_quality(
                    noodles::sam::alignment::record::MappingQuality::new(mapq).unwrap(),
                )
                .set_cigar(
                    [noodles::sam::alignment::record::cigar::Op::new(
                        noodles::sam::alignment::record::cigar::op::Kind::Match,
                        READ_LEN,
                    )]
                    .into_iter()
                    .collect(),
                )
                .set_mate_reference_sequence_id(0)
                .set_mate_alignment_start(noodles::core::Position::new(mate_pos).unwrap())
                .set_template_length(if first { span as i32 } else { -(span as i32) })
                .set_sequence(noodles::sam::alignment::record_buf::Sequence::from(
                    vec![b'A'; READ_LEN],
                ))
                .set_quality_scores(noodles::sam::alignment::record_buf::QualityScores::from(
                    vec![40u8; READ_LEN],
                ))
                .build()
        };
        // Every record, sorted by position as a queryable BAM must be.
        let mut records: Vec<(usize, noodles::sam::alignment::RecordBuf)> = Vec::new();
        for &(name, start, mapq, extra) in pairs {
            records.push((start, record(name, start, mapq, extra, true)));
            records.push((start + 200, record(name, start, mapq, extra, false)));
        }
        records.sort_by_key(|(pos, _)| *pos);

        let bam_path = dir.join("pairs.bam");
        {
            let mut writer = noodles::bam::io::writer::Builder
                .build_from_path(&bam_path)
                .unwrap();
            writer.write_header(&header).unwrap();
            for (_, r) in &records {
                writer.write_alignment_record(&header, r).unwrap();
            }
            writer.try_finish().unwrap();
        }
        let (first_record, end_of_file) = {
            let mut reader =
                noodles::bam::io::Reader::new(std::fs::File::open(&bam_path).unwrap());
            reader.read_header().unwrap();
            let first_record = reader.get_ref().virtual_position();
            let mut record = noodles::bam::Record::default();
            while reader.read_record(&mut record).unwrap() != 0 {}
            (first_record, reader.get_ref().virtual_position())
        };
        let mut bai: Vec<u8> = Vec::new();
        bai.extend_from_slice(b"BAI\x01");
        bai.extend_from_slice(&1u32.to_le_bytes());
        bai.extend_from_slice(&1u32.to_le_bytes());
        bai.extend_from_slice(&0u32.to_le_bytes());
        bai.extend_from_slice(&1u32.to_le_bytes());
        bai.extend_from_slice(&u64::from(first_record).to_le_bytes());
        bai.extend_from_slice(&u64::from(end_of_file).to_le_bytes());
        bai.extend_from_slice(&0u32.to_le_bytes());
        std::fs::write(format!("{}.bai", bam_path.display()), bai).unwrap();
        bam_path.to_str().unwrap().to_string()
    }

    #[test]
    fn test_count_resistant_counts_the_reads_the_pool_does_not_hold() {
        let dir = std::env::temp_dir().join(format!("spike_census_{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        // p_ok is editable. p_low (MAPQ 0) is resistant. p_dup is a
        // duplicate and p_sec secondary: neither is counted at all.
        let bam = write_pairs_bam(
            &dir,
            &[
                ("p_ok", 401, 60, 0),
                ("p_low", 601, 0, 0),
                ("p_dup", 801, 60, 0x400),
                ("p_sec", 1001, 60, 0x100),
            ],
        );
        let editable: HashSet<String> = ["p_ok".to_string()].into();
        let spans = [("chrA".to_string(), 300, 1400)];

        let census = count_resistant(&bam, "", &spans, &editable).unwrap();

        assert_eq!(census, Census { counted: 4, resistant: 2 });
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn test_the_census_read_in_chunks_counts_what_the_fixture_holds() {
        // The census may read a long span chunk by chunk on the thread pool.
        // Its counts must be what the fixture's own pair list says -- every
        // primary, mapped, non-duplicate, non-QC-fail read over the span, at
        // any MAPQ, resistant when its name is not editable -- in one query
        // and in five chunks on four threads.
        let dir = std::env::temp_dir().join(format!("spike_census_chunks_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let names: Vec<String> = (0..80).map(|i| format!("p{}", i)).collect();
        let pairs: Vec<(&str, usize, u8, u16)> = (0..80)
            .map(|i| {
                let extra = if i % 11 == 0 { 0x400 } else if i % 13 == 0 { 0x200 } else { 0 };
                (names[i].as_str(), 1 + 30 * i, (i % 7 * 10) as u8, extra)
            })
            .collect();
        let bam = write_pairs_bam(&dir, &pairs);
        let spans = [("chrA".to_string(), 100u64, 2_900u64)];
        let editable: HashSet<String> = names.iter().step_by(2).cloned().collect();
        let mut expected = Census::default();
        for &(name, start, _, extra) in &pairs {
            if extra & 0x600 != 0 {
                continue;
            }
            // 1-based reads of 100 bases at `start` and `start + 200`.
            for pos in [start as u64, start as u64 + 200] {
                if pos - 1 < 2_900 && pos + 99 > 100 {
                    expected.counted += 1;
                    expected.resistant += usize::from(!editable.contains(name));
                }
            }
        }
        assert!(expected.resistant > 0 && expected.resistant < expected.counted);
        assert_eq!(count_resistant_in(&bam, "", &spans, &editable, 1).unwrap(), expected);
        let chunked = rayon::ThreadPoolBuilder::new()
            .num_threads(4)
            .build()
            .unwrap()
            .install(|| count_resistant_in(&bam, "", &spans, &editable, 5).unwrap());
        assert_eq!(chunked, expected);
        let _ = std::fs::remove_dir_all(&dir);
    }
}
