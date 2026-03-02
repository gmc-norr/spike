//! VCF v4.3 input parser for spike.
//!
//! Reads SV records from a VCF file and converts them into SimEvents.
//! Supports DEL, INS, DUP, INV, and BND (fusion) SVTYPE records.

use anyhow::{bail, Context, Result};
use std::collections::HashSet;
use std::io::{BufRead, BufReader};

use crate::types::SimEvent;

/// Load simulation events from a VCF file.
///
/// Supports both plain `.vcf` and bgzip-compressed `.vcf.gz` files.
/// BND records are paired by MATEID to produce single Fusion events.
pub fn load_events_from_vcf(path: &str) -> Result<Vec<SimEvent>> {
    let file = std::fs::File::open(path)
        .with_context(|| format!("failed to open VCF: {}", path))?;

    let records = if path.ends_with(".gz") {
        let decoder = noodles::bgzf::Reader::new(file);
        let reader = BufReader::new(decoder);
        parse_vcf_records(reader)?
    } else {
        let reader = BufReader::new(file);
        parse_vcf_records(reader)?
    };

    let events = records_to_events(records)?;

    log::info!("Loaded {} events from VCF: {}", events.len(), path);
    Ok(events)
}

/// Raw parsed VCF record for SV processing.
#[derive(Debug)]
struct SvRecord {
    chrom: String,
    pos: u64,    // 0-based SV start (DEL/DUP/INV/INS) or 0-based breakpoint (BND)
    id: String,
    ref_allele: String,
    alt: String,
    info: String,
    sv_type: SvTypeTag,
}

#[derive(Debug, Clone, Copy, PartialEq)]
enum SvTypeTag {
    Del,
    Ins,
    Dup,
    Inv,
    Bnd,
    SmallVar,
}

fn parse_vcf_records<R: BufRead>(reader: R) -> Result<Vec<SvRecord>> {
    let mut records = Vec::new();

    for line in reader.lines() {
        let line = line?;
        if line.starts_with('#') {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 8 {
            continue;
        }

        let info = fields[7];
        let ref_col = fields[3];
        let alt_col = fields[4];

        let sv_type = match parse_info_field(info, "SVTYPE") {
            Some(t) => match t {
                "DEL" => SvTypeTag::Del,
                "INS" => SvTypeTag::Ins,
                "DUP" => SvTypeTag::Dup,
                "INV" => SvTypeTag::Inv,
                "BND" => SvTypeTag::Bnd,
                _ => continue,
            },
            None => {
                // No SVTYPE: check if this is a standard SNP/indel record.
                // Both REF and ALT must be pure DNA bases (no symbolic <...> alleles).
                if is_dna_allele(ref_col) && is_dna_allele(alt_col) {
                    SvTypeTag::SmallVar
                } else {
                    continue;
                }
            }
        };

        // VCF POS is 1-based. For symbolic SVs (DEL/DUP/INV/INS), POS is the
        // "preceding base" — its numeric value equals the 0-based SV start.
        // For BND, POS is the actual breakpoint position (1-based), so subtract 1.
        let raw_pos: u64 = match fields[1].parse::<u64>() {
            Ok(p) if p > 0 => p,
            _ => continue,
        };
        let pos = match sv_type {
            SvTypeTag::Bnd => raw_pos - 1, // BND: 1-based breakpoint → 0-based
            SvTypeTag::SmallVar => raw_pos - 1, // Small variant: 1-based → 0-based
            _ => raw_pos,                   // Others: 1-based preceding base == 0-based start
        };

        records.push(SvRecord {
            chrom: fields[0].to_string(),
            pos,
            id: fields[2].to_string(),
            ref_allele: ref_col.to_string(),
            alt: alt_col.to_string(),
            info: info.to_string(),
            sv_type,
        });
    }

    Ok(records)
}

/// Convert raw VCF records into SimEvents.
fn records_to_events(records: Vec<SvRecord>) -> Result<Vec<SimEvent>> {
    let mut events = Vec::new();
    let mut bnd_processed: HashSet<String> = HashSet::new();

    for record in &records {
        match record.sv_type {
            SvTypeTag::Del => {
                let end = parse_info_u64(&record.info, "END")
                    .or_else(|| {
                        parse_info_i64(&record.info, "SVLEN")
                            .map(|v| record.pos + v.unsigned_abs())
                    })
                    .unwrap_or(record.pos + 1);
                let gene = parse_info_field(&record.info, "SIM_GENE")
                    .unwrap_or("unknown")
                    .to_string();
                let af = extract_af(&record.info);

                events.push(SimEvent::Deletion {
                    chrom: record.chrom.clone(),
                    del_start: record.pos,
                    del_end: end,
                    gene,
                    exons: Vec::new(),
                    allele_fraction: af,
                });
            }
            SvTypeTag::Dup => {
                let end = parse_info_u64(&record.info, "END")
                    .or_else(|| {
                        parse_info_i64(&record.info, "SVLEN")
                            .map(|v| record.pos + v.unsigned_abs())
                    })
                    .unwrap_or(record.pos + 1);
                let gene = parse_info_field(&record.info, "SIM_GENE")
                    .unwrap_or("unknown")
                    .to_string();
                let af = extract_af(&record.info);

                events.push(SimEvent::Duplication {
                    chrom: record.chrom.clone(),
                    dup_start: record.pos,
                    dup_end: end,
                    gene,
                    allele_fraction: af,
                });
            }
            SvTypeTag::Inv => {
                let end = parse_info_u64(&record.info, "END")
                    .or_else(|| {
                        parse_info_i64(&record.info, "SVLEN")
                            .map(|v| record.pos + v.unsigned_abs())
                    })
                    .unwrap_or(record.pos + 1);
                let gene = parse_info_field(&record.info, "SIM_GENE")
                    .unwrap_or("unknown")
                    .to_string();
                let af = extract_af(&record.info);

                events.push(SimEvent::Inversion {
                    chrom: record.chrom.clone(),
                    inv_start: record.pos,
                    inv_end: end,
                    gene,
                    allele_fraction: af,
                });
            }
            SvTypeTag::Ins => {
                let ins_len = parse_info_i64(&record.info, "SVLEN")
                    .map(|v| v.unsigned_abs())
                    .unwrap_or(0);
                let gene = parse_info_field(&record.info, "SIM_GENE")
                    .unwrap_or("unknown")
                    .to_string();
                let af = extract_af(&record.info);

                // Check if ALT has explicit sequence (not symbolic <INS>).
                let ins_seq = if !record.alt.starts_with('<') && record.alt.len() > 1 {
                    // ALT contains the inserted sequence (first base = ref base).
                    let seq = record.alt.as_bytes()[1..].to_vec();
                    Some(seq)
                } else {
                    None
                };

                let effective_len = if let Some(ref seq) = ins_seq {
                    seq.len() as u64
                } else if ins_len > 0 {
                    ins_len
                } else {
                    log::warn!(
                        "INS record {} has no SVLEN and no explicit ALT sequence, skipping",
                        record.id
                    );
                    continue;
                };

                events.push(SimEvent::Insertion {
                    chrom: record.chrom.clone(),
                    pos: record.pos,
                    ins_seq,
                    ins_len: effective_len,
                    gene,
                    allele_fraction: af,
                });
            }
            SvTypeTag::SmallVar => {
                let af = extract_af(&record.info);
                let gene = parse_info_field(&record.info, "SIM_GENE")
                    .unwrap_or("unknown")
                    .to_string();

                events.push(SimEvent::SmallVariant {
                    chrom: record.chrom.clone(),
                    pos: record.pos,
                    ref_allele: record.ref_allele.as_bytes().to_vec(),
                    alt_allele: record.alt.as_bytes().to_vec(),
                    gene,
                    allele_fraction: af,
                });
            }
            SvTypeTag::Bnd => {
                // Skip if already processed as part of a MATEID pair.
                if bnd_processed.contains(&record.id) {
                    continue;
                }

                let (partner_chrom, partner_pos_1based, inverted) = parse_bnd_alt(&record.alt)
                    .with_context(|| {
                        format!(
                            "failed to parse BND ALT '{}' for record {}",
                            record.alt, record.id
                        )
                    })?;
                let partner_pos = partner_pos_1based.saturating_sub(1); // 1-based → 0-based

                // Mark both this record and its mate as processed.
                bnd_processed.insert(record.id.clone());
                if let Some(mate_id) = parse_info_field(&record.info, "MATEID") {
                    bnd_processed.insert(mate_id.to_string());
                }

                // Extract gene info.
                let gene_a = parse_info_field(&record.info, "GENE_A")
                    .unwrap_or("unknown")
                    .to_string();
                let gene_b = parse_info_field(&record.info, "GENE_B")
                    .unwrap_or("unknown")
                    .to_string();

                let af = extract_af(&record.info);

                events.push(SimEvent::Fusion {
                    chrom_a: record.chrom.clone(),
                    bp_a: record.pos,
                    gene_a,
                    chrom_b: partner_chrom,
                    bp_b: partner_pos,
                    gene_b,
                    allele_fraction: af,
                    inverted,
                });
            }
        }
    }

    Ok(events)
}

/// Parse BND ALT allele to extract partner chromosome, position, and orientation.
///
/// Handles all four VCF 4.3 BND orientations:
///   N[chr:pos[   — forward join (partner on + strand)
///   ]chr:pos]N   — forward join (reciprocal)
///   N]chr:pos]   — inverted join (partner on - strand)
///   [chr:pos[N   — inverted join (reciprocal)
///
/// Returns `(chrom, pos, inverted)`.
fn parse_bnd_alt(alt: &str) -> Result<(String, u64, bool)> {
    let open = alt.find('[').or_else(|| alt.find(']'));
    let close = alt.rfind('[').or_else(|| alt.rfind(']'));

    let (open_idx, close_idx) = match (open, close) {
        (Some(o), Some(c)) if o != c => (o, c),
        _ => bail!("no bracket pair found in BND ALT '{}'", alt),
    };

    let inner = &alt[open_idx + 1..close_idx];

    let (chrom, pos_str) = inner.split_once(':')
        .ok_or_else(|| anyhow::anyhow!("no chr:pos found in BND ALT '{}'", alt))?;

    let pos: u64 = pos_str.parse()
        .with_context(|| format!("invalid position in BND ALT '{}': '{}'", alt, pos_str))?;

    // Detect orientation from bracket characters:
    //   `[` brackets → forward join (partner + strand)
    //   `]` brackets → inverted join (partner - strand)
    let open_char = alt.as_bytes()[open_idx];
    let close_char = alt.as_bytes()[close_idx];
    let inverted = open_char == b']' && close_char == b']';

    Ok((chrom.to_string(), pos, inverted))
}

/// Extract allele fraction from INFO field.
/// Checks SIM_VAF, VAF, AF in order.
fn extract_af(info: &str) -> Option<f64> {
    for key in &["SIM_VAF", "VAF", "AF"] {
        if let Some(val) = parse_info_field(info, key) {
            if let Ok(v) = val.parse::<f64>() {
                if v > 0.0 && v <= 1.0 {
                    return Some(v);
                }
            }
        }
    }
    None
}

/// Extract a key=value from a VCF INFO field.
fn parse_info_field<'a>(info: &'a str, key: &str) -> Option<&'a str> {
    for field in info.split(';') {
        if let Some(value) = field.strip_prefix(key).and_then(|rest| rest.strip_prefix('=')) {
            return Some(value);
        }
    }
    None
}

/// Parse an integer INFO field value.
fn parse_info_u64(info: &str, key: &str) -> Option<u64> {
    parse_info_field(info, key)?.parse().ok()
}

/// Parse a signed integer INFO field value (e.g., SVLEN which can be negative).
fn parse_info_i64(info: &str, key: &str) -> Option<i64> {
    parse_info_field(info, key)?.parse().ok()
}

/// Check if a VCF allele string contains only valid DNA bases (A, C, G, T).
/// Returns false for symbolic alleles like `<DEL>`, empty strings, or alleles with non-DNA chars.
fn is_dna_allele(allele: &str) -> bool {
    !allele.is_empty()
        && !allele.starts_with('<')
        && allele
            .bytes()
            .all(|b| matches!(b.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T'))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_bnd_alt_n_bracket_open() {
        let (chrom, pos, inv) = parse_bnd_alt("N[chr9:130854070[").unwrap();
        assert_eq!(chrom, "chr9");
        assert_eq!(pos, 130854070);
        assert!(!inv, "N[chr:pos[ should be forward");
    }

    #[test]
    fn test_parse_bnd_alt_bracket_close_n() {
        let (chrom, pos, inv) = parse_bnd_alt("]chr22:23285370]N").unwrap();
        assert_eq!(chrom, "chr22");
        assert_eq!(pos, 23285370);
        assert!(inv, "]chr:pos]N should be inverted");
    }

    #[test]
    fn test_parse_bnd_alt_n_bracket_close() {
        let (chrom, pos, inv) = parse_bnd_alt("N]chr9:130854070]").unwrap();
        assert_eq!(chrom, "chr9");
        assert_eq!(pos, 130854070);
        assert!(inv, "N]chr:pos] should be inverted");
    }

    #[test]
    fn test_parse_bnd_alt_bracket_open_n() {
        let (chrom, pos, inv) = parse_bnd_alt("[chr22:23285370[N").unwrap();
        assert_eq!(chrom, "chr22");
        assert_eq!(pos, 23285370);
        assert!(!inv, "[chr:pos[N should be forward");
    }

    #[test]
    fn test_extract_af_sim_vaf() {
        assert_eq!(extract_af("SVTYPE=BND;SIM_VAF=0.050;GENE_A=BCR"), Some(0.05));
    }

    #[test]
    fn test_extract_af_vaf() {
        assert_eq!(extract_af("SVTYPE=DEL;VAF=0.121"), Some(0.121));
    }

    #[test]
    fn test_extract_af_none() {
        assert_eq!(extract_af("SVTYPE=DEL;END=100"), None);
    }

    #[test]
    fn test_parse_info_field() {
        assert_eq!(parse_info_field("SVTYPE=DEL;END=100;SVLEN=-500", "END"), Some("100"));
        assert_eq!(parse_info_field("SVTYPE=DEL;END=100;SVLEN=-500", "SVLEN"), Some("-500"));
        assert_eq!(parse_info_field("SVTYPE=DEL;END=100", "GENE"), None);
    }

    #[test]
    fn test_parse_snp_record() {
        // SNP: no SVTYPE, REF=A, ALT=T at POS=100 (1-based) → 0-based pos=99
        let vcf = "chr1\t100\ttest_snp\tA\tT\t.\t.\t.\n";
        let records = parse_vcf_records(vcf.as_bytes()).unwrap();
        let events = records_to_events(records).unwrap();
        assert_eq!(events.len(), 1);
        match &events[0] {
            SimEvent::SmallVariant { chrom, pos, ref_allele, alt_allele, .. } => {
                assert_eq!(chrom, "chr1");
                assert_eq!(*pos, 99);
                assert_eq!(ref_allele, b"A");
                assert_eq!(alt_allele, b"T");
            }
            _ => panic!("expected SmallVariant"),
        }
    }

    #[test]
    fn test_parse_small_deletion_record() {
        // Small del: REF=ACG, ALT=A at POS=100 (1-based) → 0-based pos=99
        let vcf = "chr1\t100\ttest_del\tACG\tA\t.\t.\t.\n";
        let records = parse_vcf_records(vcf.as_bytes()).unwrap();
        let events = records_to_events(records).unwrap();
        assert_eq!(events.len(), 1);
        match &events[0] {
            SimEvent::SmallVariant { pos, ref_allele, alt_allele, .. } => {
                assert_eq!(*pos, 99);
                assert_eq!(ref_allele, b"ACG");
                assert_eq!(alt_allele, b"A");
            }
            _ => panic!("expected SmallVariant"),
        }
    }

    #[test]
    fn test_parse_small_insertion_record() {
        // Small ins: REF=A, ALT=ACGT at POS=100 (1-based) → 0-based pos=99
        let vcf = "chr1\t100\ttest_ins\tA\tACGT\t.\t.\t.\n";
        let records = parse_vcf_records(vcf.as_bytes()).unwrap();
        let events = records_to_events(records).unwrap();
        assert_eq!(events.len(), 1);
        match &events[0] {
            SimEvent::SmallVariant { pos, ref_allele, alt_allele, .. } => {
                assert_eq!(*pos, 99);
                assert_eq!(ref_allele, b"A");
                assert_eq!(alt_allele, b"ACGT");
            }
            _ => panic!("expected SmallVariant"),
        }
    }

    #[test]
    fn test_parse_small_variant_with_af() {
        let vcf = "chr1\t100\tsnp1\tA\tT\t.\t.\tSIM_VAF=0.25\n";
        let records = parse_vcf_records(vcf.as_bytes()).unwrap();
        let events = records_to_events(records).unwrap();
        match &events[0] {
            SimEvent::SmallVariant { allele_fraction, .. } => {
                assert_eq!(*allele_fraction, Some(0.25));
            }
            _ => panic!("expected SmallVariant"),
        }
    }

    #[test]
    fn test_symbolic_alt_skipped() {
        // Symbolic ALT without SVTYPE should be skipped.
        let vcf = "chr1\t100\ttest\tA\t<DEL>\t.\t.\t.\n";
        let records = parse_vcf_records(vcf.as_bytes()).unwrap();
        assert!(records.is_empty());
    }

    #[test]
    fn test_is_dna_allele() {
        assert!(is_dna_allele("A"));
        assert!(is_dna_allele("ACGT"));
        assert!(is_dna_allele("acgt")); // lowercase ok
        assert!(!is_dna_allele("<DEL>"));
        assert!(!is_dna_allele(""));
        assert!(!is_dna_allele("N")); // N is not A/C/G/T
    }

    /// Test that VCF coordinates are parsed correctly for all SV types.
    /// VCF POS for DEL/DUP/INV/INS is 1-based preceding base (== 0-based SV start).
    /// VCF POS for BND is 1-based breakpoint (subtract 1 for 0-based).
    #[test]
    fn test_vcf_coordinate_parsing() {
        // DEL: POS=100 (preceding base), END=200 → 0-based [100, 200)
        let vcf = "chr1\t100\ttest_del\tN\t<DEL>\t.\t.\tSVTYPE=DEL;END=200\n";
        let records = parse_vcf_records(vcf.as_bytes()).unwrap();
        let events = records_to_events(records).unwrap();
        match &events[0] {
            SimEvent::Deletion { del_start, del_end, .. } => {
                assert_eq!(*del_start, 100);
                assert_eq!(*del_end, 200);
            }
            _ => panic!("expected Deletion"),
        }

        // DUP: POS=500, END=1000 → 0-based [500, 1000)
        let vcf = "chr1\t500\ttest_dup\tN\t<DUP>\t.\t.\tSVTYPE=DUP;END=1000\n";
        let records = parse_vcf_records(vcf.as_bytes()).unwrap();
        let events = records_to_events(records).unwrap();
        match &events[0] {
            SimEvent::Duplication { dup_start, dup_end, .. } => {
                assert_eq!(*dup_start, 500);
                assert_eq!(*dup_end, 1000);
            }
            _ => panic!("expected Duplication"),
        }

        // BND: POS=100 (1-based breakpoint) → 0-based bp = 99
        let vcf = "chr1\t100\ttest_bnd\tN\tN[chr2:200[\t.\t.\tSVTYPE=BND\n";
        let records = parse_vcf_records(vcf.as_bytes()).unwrap();
        let events = records_to_events(records).unwrap();
        match &events[0] {
            SimEvent::Fusion { bp_a, bp_b, inverted, .. } => {
                assert_eq!(*bp_a, 99);   // POS 100 → 0-based 99
                assert_eq!(*bp_b, 199);  // ALT pos 200 → 0-based 199
                assert!(!inverted, "N[chr:pos[ should be forward");
            }
            _ => panic!("expected Fusion"),
        }

        // INS: POS=300 (preceding base) → 0-based pos = 300
        let vcf = "chr1\t300\ttest_ins\tN\t<INS>\t.\t.\tSVTYPE=INS;SVLEN=50\n";
        let records = parse_vcf_records(vcf.as_bytes()).unwrap();
        let events = records_to_events(records).unwrap();
        match &events[0] {
            SimEvent::Insertion { pos, ins_len, .. } => {
                assert_eq!(*pos, 300);
                assert_eq!(*ins_len, 50);
            }
            _ => panic!("expected Insertion"),
        }
    }
}
