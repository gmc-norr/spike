//! Truth VCF writer for simulated events.
//!
//! Follows the same coordinate conventions as sv_exon:
//! - VCF POS: 1-based preceding-base (numerically == 0-based start)
//! - VCF END: 1-based inclusive

use anyhow::{Context, Result};
use std::io::Write;

use crate::reference::SharedReference;
use crate::types::{FusionJoin, SimEvent};

/// POS and ALT of the two BND records of a fusion: gene A's record, then
/// gene B's (mate) record.
pub(crate) fn bnd_records(
    chrom_a: &str,
    bp_a: u64,
    chrom_b: &str,
    bp_b: u64,
    join: FusionJoin,
) -> [(u64, String); 2] {
    // POS is the 1-based base next to the junction on each side. For a cut
    // `bp` (between 0-based bases bp-1 and bp) that is the last kept base, bp,
    // when the left side is kept, and the first kept base, bp+1, when the
    // right side is kept.
    match join {
        // A-left then B-right: t[p[ on A, ]p]t on B.
        FusionJoin::Forward => [
            (bp_a, format!("N[{}:{}[", chrom_b, bp_b + 1)),
            (bp_b + 1, format!("]{}:{}]N", chrom_a, bp_a)),
        ],
        // A-left then B-left reversed: t]p] on both.
        FusionJoin::LeftLeft => [
            (bp_a, format!("N]{}:{}]", chrom_b, bp_b)),
            (bp_b, format!("N]{}:{}]", chrom_a, bp_a)),
        ],
        // A-right reversed then B-right: [p[t on both.
        FusionJoin::RightRight => [
            (bp_a + 1, format!("[{}:{}[N", chrom_b, bp_b + 1)),
            (bp_b + 1, format!("[{}:{}[N", chrom_a, bp_a + 1)),
        ],
    }
}

/// Put the reference base `t` into a BND ALT written with `N` for it
/// (`N[p[`, `N]p]`, `]p]N`, `[p[N`).
fn with_bnd_base(alt: &str, base: &str) -> String {
    if let Some(rest) = alt.strip_prefix('N') {
        format!("{}{}", base, rest)
    } else if let Some(rest) = alt.strip_suffix('N') {
        format!("{}{}", rest, base)
    } else {
        alt.to_string()
    }
}

/// Write a truth VCF describing the simulated events.
///
/// Each event may carry a per-event allele fraction; `default_af` is used as fallback.
/// That fraction is the *request*, and goes in `SIM_REQ_VAF`. `adjusted_afs`
/// holds its simulated counterpart, one entry per event in the same order:
/// `Some(v)` when the additive cap or the two-fragment floor moved the
/// fraction off the request, `None` when the request stands. `SIM_VAF` is
/// what was simulated, so it is `v` where there is one and the request
/// everywhere else.
pub fn write_truth_vcf(
    events: &[SimEvent],
    adjusted_afs: &[Option<f64>],
    default_af: f64,
    output_path: &str,
    ref_path: &str,
    reference: &SharedReference,
    contigs: &[(String, u64)],
) -> Result<()> {
    let mut f = std::fs::File::create(output_path)
        .with_context(|| format!("failed to create truth VCF: {}", output_path))?;

    // Header.
    writeln!(f, "##fileformat=VCFv4.3")?;
    writeln!(f, "##fileDate={}", chrono_date())?;
    writeln!(f, "##source=spike")?;
    writeln!(f, "##reference={}", ref_path)?;
    for (name, length) in contigs {
        writeln!(f, "##contig=<ID={},length={}>", name, length)?;
    }
    writeln!(f, "##ALT=<ID=DEL,Description=\"Deletion\">")?;
    writeln!(f, "##ALT=<ID=DUP,Description=\"Tandem duplication\">")?;
    writeln!(f, "##ALT=<ID=INV,Description=\"Inversion\">")?;
    writeln!(f, "##ALT=<ID=INS,Description=\"Insertion\">")?;
    writeln!(f, "##ALT=<ID=BND,Description=\"Translocation breakend\">")?;
    writeln!(
        f,
        "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">"
    )?;
    writeln!(
        f,
        "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position\">"
    )?;
    writeln!(
        f,
        "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"SV length\">"
    )?;
    writeln!(
        f,
        "##INFO=<ID=SIM_VAF,Number=1,Type=Float,Description=\"Allele fraction that was \
         simulated: the fraction of the depth the fragments spike planted make up. Below \
         SIM_REQ_VAF where the additive 0.95 cap applied and above it where the \
         two-fragment floor did. On an additive event (a fusion or a DUP under \
         --dup-model junction) it is the junction evidence -- the fraction the fragments \
         across the breakpoint make up; a junction DUP also plants interior depth copies \
         at SIM_REQ_VAF so a capped one's interior dosage is above this\">"
    )?;
    writeln!(
        f,
        "##INFO=<ID=SIM_REQ_VAF,Number=1,Type=Float,Description=\"Allele fraction requested \
         for this event; SIM_VAF is the fraction that was simulated\">"
    )?;
    writeln!(
        f,
        "##INFO=<ID=SIM_GENE,Number=1,Type=String,Description=\"Affected gene\">"
    )?;
    writeln!(
        f,
        "##INFO=<ID=SIM_EXONS,Number=.,Type=String,Description=\"Affected exons\">"
    )?;
    writeln!(
        f,
        "##INFO=<ID=MATEID,Number=1,Type=String,Description=\"ID of mate breakend\">"
    )?;
    writeln!(f, "##FILTER=<ID=PASS,Description=\"All filters passed\">")?;
    writeln!(
        f,
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
    )?;
    writeln!(
        f,
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE"
    )?;

    // Records as (CHROM, POS, rest of the line). Sorted before writing so the
    // file can be indexed.
    let mut records: Vec<(String, u64, String)> = Vec::new();
    // The reference base at a 1-based POS (N if it can't be read).
    let base_at = |chrom: &str, pos: u64| -> String {
        reference
            .fetch_sequence(chrom, pos.saturating_sub(1), pos)
            .ok()
            .filter(|b| b.len() == 1)
            .map(|b| (b[0] as char).to_string())
            .unwrap_or_else(|| "N".to_string())
    };

    for (i, event) in events.iter().enumerate() {
        let requested_af = event.allele_fraction().unwrap_or(default_af);
        // What was simulated: the request, unless the additive cap or the
        // two-fragment floor moved it.
        let event_af = adjusted_afs
            .get(i)
            .copied()
            .flatten()
            .unwrap_or(requested_af);
        // The genotype says which copies carry the event, which neither the
        // cap nor the floor changes, so it stays the request's.
        let gt = genotype_from_vaf(requested_af);

        match event {
            SimEvent::Deletion {
                chrom,
                del_start,
                del_end,
                gene,
                exons,
                ..
            } => {
                let sv_len = *del_end as i64 - *del_start as i64;
                let exons_str = if exons.is_empty() {
                    ".".to_string()
                } else {
                    exons.join(",")
                };
                // VCF POS: numerically == 0-based start; END: 0-based
                // half-open value == 1-based inclusive.
                records.push((
                    chrom.clone(),
                    *del_start,
                    format!(
                        "sim_del_{}\t{}\t<DEL>\t999\tPASS\tSVTYPE=DEL;END={};SVLEN=-{};SIM_VAF={:.3};SIM_REQ_VAF={:.3};SIM_GENE={};SIM_EXONS={}\tGT\t{}",
                        i + 1,
                        base_at(chrom, *del_start),
                        del_end,
                        sv_len,
                        event_af,
                        requested_af,
                        gene,
                        exons_str,
                        gt,
                    ),
                ));
            }
            SimEvent::Fusion {
                chrom_a,
                bp_a,
                gene_a,
                chrom_b,
                bp_b,
                gene_b,
                join,
                ..
            } => {
                let id_a = format!("sim_fus_{}", i + 1);
                let id_b = format!("sim_fus_{}_mate", i + 1);
                let [(pos_a, alt_a), (pos_b, alt_b)] =
                    bnd_records(chrom_a, *bp_a, chrom_b, *bp_b, *join);
                for (chrom, pos, id, alt, mate, gene) in [
                    (chrom_a, pos_a, &id_a, alt_a, &id_b, gene_a),
                    (chrom_b, pos_b, &id_b, alt_b, &id_a, gene_b),
                ] {
                    let base = base_at(chrom, pos);
                    records.push((
                        chrom.clone(),
                        pos,
                        format!(
                            "{}\t{}\t{}\t999\tPASS\tSVTYPE=BND;MATEID={};SIM_VAF={:.3};SIM_REQ_VAF={:.3};SIM_GENE={}\tGT\t{}",
                            id,
                            base,
                            with_bnd_base(&alt, &base),
                            mate,
                            event_af,
                            requested_af,
                            gene,
                            gt,
                        ),
                    ));
                }
            }
            SimEvent::Duplication {
                chrom,
                dup_start,
                dup_end,
                gene,
                ..
            } => {
                let sv_len = *dup_end as i64 - *dup_start as i64;
                records.push((
                    chrom.clone(),
                    *dup_start,
                    format!(
                        "sim_dup_{}\t{}\t<DUP>\t999\tPASS\tSVTYPE=DUP;END={};SVLEN={};SIM_VAF={:.3};SIM_REQ_VAF={:.3};SIM_GENE={}\tGT\t{}",
                        i + 1,
                        base_at(chrom, *dup_start),
                        dup_end,
                        sv_len,
                        event_af,
                        requested_af,
                        gene,
                        gt,
                    ),
                ));
            }
            SimEvent::Inversion {
                chrom,
                inv_start,
                inv_end,
                gene,
                ..
            } => {
                let sv_len = *inv_end as i64 - *inv_start as i64;
                records.push((
                    chrom.clone(),
                    *inv_start,
                    format!(
                        "sim_inv_{}\t{}\t<INV>\t999\tPASS\tSVTYPE=INV;END={};SVLEN={};SIM_VAF={:.3};SIM_REQ_VAF={:.3};SIM_GENE={}\tGT\t{}",
                        i + 1,
                        base_at(chrom, *inv_start),
                        inv_end,
                        sv_len,
                        event_af,
                        requested_af,
                        gene,
                        gt,
                    ),
                ));
            }
            SimEvent::Insertion {
                chrom,
                pos,
                ins_seq,
                ins_len,
                gene,
                ..
            } => {
                // Sequence-resolved ALT: the anchor base plus the inserted
                // bases, uppercased exactly as `from_insertion` uppercases
                // them, so the truth spells the sequence the reads carry.
                // An event whose sequence was never resolved has none to
                // write and keeps the symbolic ALT.
                let anchor = base_at(chrom, *pos);
                let alt = match ins_seq {
                    Some(seq) => {
                        let mut a = anchor.clone();
                        a.extend(seq.iter().map(|b| b.to_ascii_uppercase() as char));
                        a
                    }
                    None => "<INS>".to_string(),
                };
                records.push((
                    chrom.clone(),
                    *pos,
                    format!(
                        "sim_ins_{}\t{}\t{}\t999\tPASS\tSVTYPE=INS;SVLEN={};SIM_VAF={:.3};SIM_REQ_VAF={:.3};SIM_GENE={}\tGT\t{}",
                        i + 1,
                        anchor,
                        alt,
                        ins_len,
                        event_af,
                        requested_af,
                        gene,
                        gt,
                    ),
                ));
            }
            SimEvent::SmallVariant {
                chrom,
                pos,
                ref_allele,
                alt_allele,
                gene,
                ..
            } => {
                // VCF POS is 1-based.
                records.push((
                    chrom.clone(),
                    pos + 1,
                    format!(
                        "sim_var_{}\t{}\t{}\t999\tPASS\tSIM_VAF={:.3};SIM_REQ_VAF={:.3};SIM_GENE={}\tGT\t{}",
                        i + 1,
                        String::from_utf8_lossy(ref_allele),
                        String::from_utf8_lossy(alt_allele),
                        event_af,
                        requested_af,
                        gene,
                        gt,
                    ),
                ));
            }
        }
    }

    // Sort by contig order in the reference index, then position.
    let contig_rank: std::collections::HashMap<&str, usize> = contigs
        .iter()
        .enumerate()
        .map(|(rank, (name, _))| (name.as_str(), rank))
        .collect();
    records.sort_by_key(|(chrom, pos, _)| {
        (contig_rank.get(chrom.as_str()).copied().unwrap_or(usize::MAX), chrom.clone(), *pos)
    });
    for (chrom, pos, rest) in &records {
        writeln!(f, "{}\t{}\t{}", chrom, pos, rest)?;
    }

    log::info!(
        "Wrote truth VCF with {} events to {}",
        events.len(),
        output_path
    );
    Ok(())
}

fn genotype_from_vaf(vaf: f64) -> &'static str {
    if vaf >= 0.9 {
        "1/1"
    } else {
        "0/1"
    }
}

/// Simple date formatting without chrono dependency.
fn chrono_date() -> String {
    let now = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap_or(std::time::Duration::ZERO)
        .as_secs();
    let mut days = (now / 86400) as i64;

    let mut year = 1970i64;
    loop {
        let days_in_year = if year % 4 == 0 && (year % 100 != 0 || year % 400 == 0) {
            366
        } else {
            365
        };
        if days < days_in_year {
            break;
        }
        days -= days_in_year;
        year += 1;
    }

    let leap = year % 4 == 0 && (year % 100 != 0 || year % 400 == 0);
    let month_days: [i64; 12] = [
        31,
        if leap { 29 } else { 28 },
        31,
        30,
        31,
        30,
        31,
        31,
        30,
        31,
        30,
        31,
    ];
    let mut month = 0usize;
    for (i, &md) in month_days.iter().enumerate() {
        if days < md {
            month = i;
            break;
        }
        days -= md;
    }

    format!("{}{:02}{:02}", year, month + 1, days + 1)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_truth_vcf_is_sorted_with_contigs_and_reference_bases() {
        //            1-based: 1234567890123456
        let chr1 = b"GATTACAGATTACA".to_vec();
        let chr2 = b"CCCCGGGGTTTTAAAA".to_vec();
        let reference = SharedReference::from_sequences(
            [("chr1".to_string(), chr1), ("chr2".to_string(), chr2)].into(),
        );
        let contigs = vec![
            ("chr1".to_string(), 14),
            ("chr2".to_string(), 16),
            ("chr3".to_string(), 100),
        ];
        let del = |chrom: &str, start: u64, end: u64| SimEvent::Deletion {
            chrom: chrom.to_string(),
            del_start: start,
            del_end: end,
            gene: "G".to_string(),
            exons: vec![],
            allele_fraction: None,
        };
        // Given out of order on purpose.
        let events = vec![
            del("chr2", 6, 10),
            SimEvent::SmallVariant {
                chrom: "chr1".to_string(),
                pos: 9,
                ref_allele: b"T".to_vec(),
                alt_allele: b"C".to_vec(),
                gene: "G".to_string(),
                allele_fraction: None,
            },
            SimEvent::Fusion {
                chrom_a: "chr1".to_string(),
                bp_a: 4,
                gene_a: "A".to_string(),
                chrom_b: "chr2".to_string(),
                bp_b: 11,
                gene_b: "B".to_string(),
                allele_fraction: None,
                join: FusionJoin::Forward,
            },
            del("chr1", 2, 5),
        ];
        let path = std::env::temp_dir().join(format!("spike_truth_{}.vcf", std::process::id()));
        write_truth_vcf(&events, &[None; 4], 0.5, path.to_str().unwrap(), "ref.fa", &reference, &contigs)
            .unwrap();
        let text = std::fs::read_to_string(&path).unwrap();
        std::fs::remove_file(&path).ok();

        let contig_lines: Vec<&str> = text.lines().filter(|l| l.starts_with("##contig")).collect();
        assert_eq!(
            contig_lines,
            vec![
                "##contig=<ID=chr1,length=14>",
                "##contig=<ID=chr2,length=16>",
                "##contig=<ID=chr3,length=100>",
            ]
        );

        // CHROM, POS, ID, REF, ALT of each record, in file order.
        let records: Vec<Vec<&str>> = text
            .lines()
            .filter(|l| !l.starts_with('#'))
            .map(|l| l.split('\t').take(5).collect())
            .collect();
        assert_eq!(
            records,
            vec![
                vec!["chr1", "2", "sim_del_4", "A", "<DEL>"],
                vec!["chr1", "4", "sim_fus_3", "T", "T[chr2:12["],
                vec!["chr1", "10", "sim_var_2", "T", "C"],
                vec!["chr2", "6", "sim_del_1", "G", "<DEL>"],
                vec!["chr2", "12", "sim_fus_3_mate", "T", "]chr1:4]T"],
            ]
        );
    }

    fn records(bp_a: u64, bp_b: u64, join: FusionJoin) -> [(u64, String); 2] {
        bnd_records("chr1", bp_a, "chr2", bp_b, join)
    }

    #[test]
    fn test_bnd_records_forward() {
        // chr1 up to and including base 100, then chr2 from base 200.
        assert_eq!(
            records(100, 199, FusionJoin::Forward),
            [(100, "N[chr2:200[".to_string()), (200, "]chr1:100]N".to_string())]
        );
    }

    #[test]
    fn test_bnd_records_left_left() {
        // chr1 up to base 100, then chr2 up to base 200, reversed.
        assert_eq!(
            records(100, 200, FusionJoin::LeftLeft),
            [(100, "N]chr2:200]".to_string()), (200, "N]chr1:100]".to_string())]
        );
    }

    #[test]
    fn test_bnd_records_right_right() {
        // chr1 from base 100, reversed, then chr2 from base 200.
        assert_eq!(
            records(99, 199, FusionJoin::RightRight),
            [(100, "[chr2:200[N".to_string()), (200, "[chr1:100[N".to_string())]
        );
    }

    /// chr1 = GATTACAGATTACA; an insertion at 0-based 7 is anchored on the
    /// A at 0-based 6 (1-based POS 7).
    fn ins_reference() -> SharedReference {
        SharedReference::from_sequences(
            [("chr1".to_string(), b"GATTACAGATTACA".to_vec())].into(),
        )
    }

    /// The INS record's CHROM, POS, ID, REF, ALT and INFO. `tag` keeps
    /// concurrently running callers off each other's temporary file.
    fn ins_record(tag: &str, ins_seq: Option<Vec<u8>>, ins_len: u64) -> Vec<String> {
        let reference = ins_reference();
        let contigs = vec![("chr1".to_string(), 14)];
        let events = vec![SimEvent::Insertion {
            chrom: "chr1".to_string(),
            pos: 7,
            ins_seq,
            ins_len,
            gene: "G".to_string(),
            allele_fraction: None,
        }];
        let path = std::env::temp_dir()
            .join(format!("spike_truth_ins_{}_{}.vcf", std::process::id(), tag));
        write_truth_vcf(&events, &[None; 4], 0.5, path.to_str().unwrap(), "ref.fa", &reference, &contigs)
            .unwrap();
        let text = std::fs::read_to_string(&path).unwrap();
        std::fs::remove_file(&path).ok();
        let line = text
            .lines()
            .find(|l| !l.starts_with('#'))
            .expect("one INS record");
        let f: Vec<&str> = line.split('\t').collect();
        vec![
            f[0].to_string(),
            f[1].to_string(),
            f[2].to_string(),
            f[3].to_string(),
            f[4].to_string(),
            f[7].to_string(),
        ]
    }

    #[test]
    fn test_truth_ins_alt_is_anchor_plus_supplied_sequence() {
        let r = ins_record("supplied", Some(b"GGGG".to_vec()), 4);
        assert_eq!(&r[..5], ["chr1", "7", "sim_ins_1", "A", "AGGGG"]);
        assert!(
            r[5].starts_with("SVTYPE=INS;SVLEN=4;"),
            "SVTYPE and SVLEN unchanged, got {}",
            r[5]
        );
    }

    #[test]
    fn test_truth_ins_alt_uppercases_supplied_sequence() {
        // `from_insertion` uppercases before the reads are cut from it, so
        // truth must uppercase too or the ALT is not the reads' sequence.
        let r = ins_record("lowercase", Some(b"ggtt".to_vec()), 4);
        assert_eq!(r[4], "AGGTT");
    }

    /// A one-DEL truth VCF for a request of `af`, with `simulated` saying what
    /// the tiling actually planted (`None` when the request stands). `tag`
    /// keeps concurrently running callers off each other's temporary file.
    fn af_truth_text(tag: &str, af: f64, simulated: Option<f64>) -> String {
        let reference = SharedReference::from_sequences(
            [("chr1".to_string(), b"GATTACAGATTACA".to_vec())].into(),
        );
        let contigs = vec![("chr1".to_string(), 14)];
        let events = vec![SimEvent::Deletion {
            chrom: "chr1".to_string(),
            del_start: 2,
            del_end: 5,
            gene: "G".to_string(),
            exons: vec![],
            allele_fraction: Some(af),
        }];
        let path = std::env::temp_dir()
            .join(format!("spike_truth_af_{}_{}.vcf", std::process::id(), tag));
        write_truth_vcf(
            &events,
            &[simulated],
            0.5,
            path.to_str().unwrap(),
            "ref.fa",
            &reference,
            &contigs,
        )
        .unwrap();
        let text = std::fs::read_to_string(&path).unwrap();
        std::fs::remove_file(&path).ok();
        text
    }

    /// The INFO column of the single record in `text`.
    fn only_info(text: &str) -> String {
        let line = text
            .lines()
            .find(|l| !l.starts_with('#'))
            .expect("one record");
        line.split('\t').nth(7).expect("an INFO column").to_string()
    }

    #[test]
    fn test_truth_records_the_simulated_fraction_and_the_request() {
        // The two-fragment floor planted 0.317 where 0.05 was asked for.
        let text = af_truth_text("floored", 0.05, Some(0.31746031746031744));
        let info = only_info(&text);
        assert!(
            info.contains("SIM_VAF=0.317;SIM_REQ_VAF=0.050"),
            "SIM_VAF must be what was simulated and SIM_REQ_VAF the request, got {}",
            info
        );
    }

    #[test]
    fn test_truth_keeps_sim_vaf_at_the_request_when_nothing_moved_it() {
        // Neither cap nor floor applied, so SIM_VAF is exactly the text it
        // has always been, and the new field repeats it.
        let text = af_truth_text("unadjusted", 0.5, None);
        let info = only_info(&text);
        assert!(
            info.contains("SIM_VAF=0.500;SIM_REQ_VAF=0.500"),
            "an unadjusted event's SIM_VAF must not move, got {}",
            info
        );
    }

    #[test]
    fn test_truth_header_declares_the_requested_fraction() {
        let text = af_truth_text("header", 0.5, None);
        assert!(
            text.contains("##INFO=<ID=SIM_REQ_VAF,Number=1,Type=Float,"),
            "SIM_REQ_VAF must be declared, header was:\n{}",
            text
        );
        let sim_vaf = text
            .lines()
            .find(|l| l.starts_with("##INFO=<ID=SIM_VAF,"))
            .expect("a SIM_VAF header line");
        assert!(
            sim_vaf.contains("SIM_REQ_VAF"),
            "SIM_VAF's description must point at the request, got {}",
            sim_vaf
        );
        // "the fraction of the depth the fragments spike planted make up" is
        // not the whole story for an additive event: a junction DUP at a
        // capped af=0.99 plants 1900 junction fragments (0.950) *and* 249
        // interior depth copies scaled by the uncapped 0.99, so its interior
        // dosage realises about 0.996 while SIM_VAF records 0.950. The
        // recorded number is the junction evidence, and the description has
        // to say so.
        assert!(
            sim_vaf.contains("junction"),
            "SIM_VAF's description must say that on an additive event the \
             fraction is the junction evidence, got {}",
            sim_vaf
        );
    }
}
