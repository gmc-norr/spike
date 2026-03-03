//! Exon BED parsing and event specification parsing.
//!
//! Replicated from sv_exon to avoid the fermi-lite C FFI build dependency.

use anyhow::{bail, Context, Result};
use std::collections::HashMap;

use crate::types::SimEvent;

/// Allele fraction specification from event syntax.
#[derive(Debug, Clone, PartialEq)]
pub enum AfSpec {
    /// Exact allele fraction value.
    Exact(f64),
    /// Germline heterozygous: sample from Beta(40,40) centered at 0.5.
    Het,
    /// Homozygous: fixed at 1.0.
    Hom,
}

/// Exon region definition.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct Exon {
    pub chrom: String,
    pub start: u64, // 0-based
    pub end: u64,   // 0-based half-open
    pub name: String,
    pub gene: String,
    pub number: u32, // 1-based exon ordinal
}

/// Gene target with exons.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct GeneTarget {
    pub gene: String,
    pub chrom: String,
    pub gene_start: u64, // min exon start (0-based)
    pub gene_end: u64,   // max exon end (0-based half-open)
    pub exons: Vec<Exon>,
}

/// Parse an exon BED file into gene targets.
///
/// Expected format: tab-separated, at least 4 columns:
///   chrom  start  end  name  [gene]
///
/// If 5th column (gene) is absent, gene is parsed from name (e.g. "LDLR_exon1" -> "LDLR").
pub fn parse_exon_bed(path: &str) -> Result<Vec<GeneTarget>> {
    let content =
        std::fs::read_to_string(path).with_context(|| format!("failed to read BED: {}", path))?;

    let mut gene_exons: HashMap<(String, String), Vec<Exon>> = HashMap::new();

    for (line_no, line) in content.lines().enumerate() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') || line.starts_with("track") {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 4 {
            bail!(
                "BED line {} has fewer than 4 columns: {}",
                line_no + 1,
                line
            );
        }

        let chrom = fields[0].to_string();
        let start: u64 = fields[1]
            .parse()
            .with_context(|| format!("invalid start at line {}", line_no + 1))?;
        let end: u64 = fields[2]
            .parse()
            .with_context(|| format!("invalid end at line {}", line_no + 1))?;
        let name = fields[3].to_string();

        let gene = if fields.len() >= 5 && !fields[4].is_empty() {
            fields[4].to_string()
        } else {
            name.split('_').next().unwrap_or(&name).to_string()
        };

        gene_exons
            .entry((gene.clone(), chrom.clone()))
            .or_default()
            .push(Exon {
                chrom,
                start,
                end,
                name,
                gene,
                number: 0,
            });
    }

    let mut targets = Vec::new();
    for ((gene, _chrom), mut exons) in gene_exons {
        exons.sort_by_key(|e| e.start);
        for (i, exon) in exons.iter_mut().enumerate() {
            exon.number = (i + 1) as u32;
        }
        let chrom = exons[0].chrom.clone();
        let gene_start = exons.iter().map(|e| e.start).min().unwrap();
        let gene_end = exons.iter().map(|e| e.end).max().unwrap();
        targets.push(GeneTarget {
            gene,
            chrom,
            gene_start,
            gene_end,
            exons,
        });
    }

    targets.sort_by(|a, b| a.chrom.cmp(&b.chrom).then(a.gene_start.cmp(&b.gene_start)));
    Ok(targets)
}

/// Parse an event specification string into a SimEvent and optional AF spec.
///
/// Supported formats:
///   "del:chr20:30000000-30005000"             — deletion by coordinates
///   "del:GENE:exon4-exon8"                   — deletion by exon range (requires --exon-bed)
///   "fusion:GENEA:exon14:GENEB:exon2"        — gene fusion by exon boundaries (requires --exon-bed)
///   "dup:chr20:30000000-30005000"             — tandem duplication
///   "inv:chr20:30000000-30005000"             — inversion
///   "ins:chr20:30000000:500"                  — insertion (500bp random seq)
///   "del:GENE:exon4-exon8;af=0.15"           — with explicit allele fraction
///   "fusion:GENEA:exon14:GENEB:exon2;af=het" — with heterozygous AF distribution
///   "del:GENE:exon4-exon8;af=hom"            — with homozygous AF (1.0)
pub fn parse_event_spec(spec: &str, genes: &[GeneTarget]) -> Result<(SimEvent, Option<AfSpec>)> {
    // Split on ';' to separate event spec from options.
    let (event_part, options_part) = match spec.split_once(';') {
        Some((ev, opts)) => (ev, Some(opts)),
        None => (spec, None),
    };

    let parts: Vec<&str> = event_part.split(':').collect();
    if parts.is_empty() {
        bail!("empty event specification");
    }

    let mut event = match parts[0].to_lowercase().as_str() {
        "del" => parse_region_spec(&parts[1..], "del", genes)?,
        "fusion" => parse_fusion_spec(&parts[1..], genes)?,
        "dup" => parse_region_spec(&parts[1..], "dup", genes)?,
        "inv" => parse_region_spec(&parts[1..], "inv", genes)?,
        "ins" => parse_insertion_spec(&parts[1..])?,
        "snp" => parse_snp_spec(&parts[1..])?,
        other => bail!(
            "unknown event type '{}', expected 'del', 'fusion', 'dup', 'inv', 'ins', or 'snp'",
            other
        ),
    };

    // Parse options (currently only af=...).
    let af_spec = if let Some(opts) = options_part {
        parse_event_options(opts, &mut event)?
    } else {
        None
    };

    Ok((event, af_spec))
}

/// Parse key=value options from the part after ';'.
fn parse_event_options(opts: &str, _event: &mut SimEvent) -> Result<Option<AfSpec>> {
    let mut af_spec = None;
    for kv in opts.split(';') {
        let kv = kv.trim();
        if kv.is_empty() {
            continue;
        }
        if let Some((key, value)) = kv.split_once('=') {
            match key.to_lowercase().as_str() {
                "af" => {
                    af_spec = Some(parse_af_value(value)?);
                }
                other => bail!("unknown event option '{}'", other),
            }
        } else {
            bail!("invalid option format '{}', expected key=value", kv);
        }
    }
    Ok(af_spec)
}

/// Parse an AF value: a number, "het", or "hom".
fn parse_af_value(value: &str) -> Result<AfSpec> {
    let value = value.trim().to_lowercase();
    match value.as_str() {
        "het" => Ok(AfSpec::Het),
        "hom" => Ok(AfSpec::Hom),
        _ => {
            let v: f64 = value.parse().with_context(|| {
                format!(
                    "invalid af value '{}', expected number, 'het', or 'hom'",
                    value
                )
            })?;
            if v <= 0.0 || v > 1.0 {
                bail!("af must be in (0.0, 1.0], got {}", v);
            }
            Ok(AfSpec::Exact(v))
        }
    }
}

/// Unified region parser for DEL, DUP, and INV events.
///
/// Supports both coordinate-based and exon-based specifications:
///   "del:chr20:30000000-30005000"     — coordinate-based
///   "del:GENE:exon4-exon8"           — exon-based (requires --exon-bed)
///   "dup:GENE:exon4-exon8"           — exon-based duplication
///   "inv:chr20:30000000-30005000"     — coordinate-based inversion
fn parse_region_spec(parts: &[&str], sv_type: &str, genes: &[GeneTarget]) -> Result<SimEvent> {
    if parts.len() < 2 {
        bail!(
            "{} spec requires at least 2 parts: '{}:CHR:START-END' or '{}:GENE:exonN-exonM'",
            sv_type,
            sv_type,
            sv_type
        );
    }

    let first = parts[0];
    let second = parts[1];

    // Try coordinate-based: second has a dash with numbers on both sides.
    if let Some((start_str, end_str)) = second.split_once('-') {
        if let (Ok(start), Ok(end)) = (start_str.parse::<u64>(), end_str.parse::<u64>()) {
            return make_region_event(
                sv_type,
                first.to_string(),
                start,
                end,
                "unknown".to_string(),
                Vec::new(),
            );
        }
    }

    // Exon-based: "sv_type:GENE:exonN-exonM"
    let gene_name = first;
    let gene = genes
        .iter()
        .find(|g| g.gene.eq_ignore_ascii_case(gene_name))
        .ok_or_else(|| {
            anyhow::anyhow!(
                "gene '{}' not found. Available: {}",
                gene_name,
                genes
                    .iter()
                    .map(|g| g.gene.as_str())
                    .collect::<Vec<_>>()
                    .join(", ")
            )
        })?;

    let exon_range = second;
    let (start_exon_str, end_exon_str) = exon_range.split_once('-').ok_or_else(|| {
        anyhow::anyhow!("expected exon range 'exonN-exonM', got '{}'", exon_range)
    })?;

    let start_num = parse_exon_number(start_exon_str)?;
    let end_num = parse_exon_number(end_exon_str)?;

    let (start_num, end_num) = if start_num <= end_num {
        (start_num, end_num)
    } else {
        (end_num, start_num)
    };

    let start_exon = gene
        .exons
        .iter()
        .find(|e| e.number == start_num)
        .ok_or_else(|| anyhow::anyhow!("exon {} not found in {}", start_num, gene_name))?;
    let end_exon = gene
        .exons
        .iter()
        .find(|e| e.number == end_num)
        .ok_or_else(|| anyhow::anyhow!("exon {} not found in {}", end_num, gene_name))?;

    let affected_exons: Vec<String> = gene
        .exons
        .iter()
        .filter(|e| e.number >= start_num && e.number <= end_num)
        .map(|e| e.name.clone())
        .collect();

    make_region_event(
        sv_type,
        gene.chrom.clone(),
        start_exon.start,
        end_exon.end,
        gene.gene.clone(),
        affected_exons,
    )
}

/// Create the appropriate SimEvent for a region-based SV type.
fn make_region_event(
    sv_type: &str,
    chrom: String,
    start: u64,
    end: u64,
    gene: String,
    exons: Vec<String>,
) -> Result<SimEvent> {
    match sv_type {
        "del" => Ok(SimEvent::Deletion {
            chrom,
            del_start: start,
            del_end: end,
            gene,
            exons,
            allele_fraction: None,
        }),
        "dup" => Ok(SimEvent::Duplication {
            chrom,
            dup_start: start,
            dup_end: end,
            gene,
            allele_fraction: None,
        }),
        "inv" => Ok(SimEvent::Inversion {
            chrom,
            inv_start: start,
            inv_end: end,
            gene,
            allele_fraction: None,
        }),
        _ => bail!("unexpected sv_type '{}' in make_region_event", sv_type),
    }
}

fn parse_fusion_spec(parts: &[&str], genes: &[GeneTarget]) -> Result<SimEvent> {
    // "fusion:GENEA:exonN:GENEB:exonM" or "fusion:GENEA:exonN:GENEB:exonM:inv"
    if parts.len() < 4 {
        bail!("fusion spec requires at least 4 parts: 'fusion:GENEA:exonN:GENEB:exonM[:inv]'");
    }

    let gene_a_name = parts[0];
    let exon_a_str = parts[1];
    let gene_b_name = parts[2];
    let exon_b_str = parts[3];

    // Check for optional ":inv" suffix.
    let inverted = parts.len() >= 5 && parts[4].eq_ignore_ascii_case("inv");

    let gene_a = genes
        .iter()
        .find(|g| g.gene.eq_ignore_ascii_case(gene_a_name))
        .ok_or_else(|| anyhow::anyhow!("gene '{}' not found", gene_a_name))?;
    let gene_b = genes
        .iter()
        .find(|g| g.gene.eq_ignore_ascii_case(gene_b_name))
        .ok_or_else(|| anyhow::anyhow!("gene '{}' not found", gene_b_name))?;

    let exon_a_num = parse_exon_number(exon_a_str)?;
    let exon_b_num = parse_exon_number(exon_b_str)?;

    let exon_a = gene_a
        .exons
        .iter()
        .find(|e| e.number == exon_a_num)
        .ok_or_else(|| anyhow::anyhow!("exon {} not found in {}", exon_a_num, gene_a_name))?;
    let exon_b = gene_b
        .exons
        .iter()
        .find(|e| e.number == exon_b_num)
        .ok_or_else(|| anyhow::anyhow!("exon {} not found in {}", exon_b_num, gene_b_name))?;

    // Breakpoint: end of exon A (intron boundary) joined to start of exon B.
    Ok(SimEvent::Fusion {
        chrom_a: gene_a.chrom.clone(),
        bp_a: exon_a.end, // first base NOT included from gene A
        gene_a: gene_a.gene.clone(),
        chrom_b: gene_b.chrom.clone(),
        bp_b: exon_b.start, // first base included from gene B
        gene_b: gene_b.gene.clone(),
        allele_fraction: None,
        inverted,
    })
}

/// Parse an insertion spec.
///
/// Formats:
///   "ins:chr20:30000000:500"          — random insertion of 500bp
///   "ins:chr20:30000000:ACGTACGT"     — explicit DNA insertion sequence
fn parse_insertion_spec(parts: &[&str]) -> Result<SimEvent> {
    if parts.len() < 3 {
        bail!("ins spec requires: 'ins:CHR:POS:LENGTH' or 'ins:CHR:POS:SEQUENCE'");
    }

    let chrom = parts[0].to_string();
    let pos: u64 = parts[1]
        .parse()
        .with_context(|| format!("invalid position: '{}'", parts[1]))?;

    // Try parsing third part as a length (u64). If that fails, treat as DNA sequence.
    if let Ok(ins_len) = parts[2].parse::<u64>() {
        if ins_len == 0 {
            bail!("insertion length must be > 0");
        }
        Ok(SimEvent::Insertion {
            chrom,
            pos,
            ins_seq: None,
            ins_len,
            gene: "unknown".to_string(),
            allele_fraction: None,
        })
    } else {
        // Validate as DNA sequence.
        let seq_str = parts[2];
        if seq_str.is_empty() {
            bail!("insertion sequence must be non-empty");
        }
        for &b in seq_str.as_bytes() {
            if !matches!(b.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T') {
                bail!(
                    "invalid base '{}' in insertion sequence, expected A/C/G/T",
                    b as char
                );
            }
        }
        let seq: Vec<u8> = seq_str.bytes().map(|b| b.to_ascii_uppercase()).collect();
        let len = seq.len() as u64;
        Ok(SimEvent::Insertion {
            chrom,
            pos,
            ins_seq: Some(seq),
            ins_len: len,
            gene: "unknown".to_string(),
            allele_fraction: None,
        })
    }
}

/// Parse a SNP/small variant spec.
///
/// Formats:
///   "snp:chr1:100:A:T"    — explicit REF and ALT bases
///   "snp:chr1:100:A>T"    — REF>ALT shorthand
///   "snp:chr1:100:ACG:A"  — small deletion (REF=ACG, ALT=A)
///   "snp:chr1:100:A:ACGT" — small insertion (REF=A, ALT=ACGT)
///
/// Position (`POS`) is 1-based in the CLI syntax and converted to internal
/// 0-based coordinates.
fn parse_snp_spec(parts: &[&str]) -> Result<SimEvent> {
    if parts.len() < 3 {
        bail!("snp spec requires: 'snp:CHR:POS:REF:ALT' or 'snp:CHR:POS:REF>ALT'");
    }

    let chrom = parts[0].to_string();
    let pos_1based: u64 = parts[1]
        .parse()
        .with_context(|| format!("invalid position: '{}'", parts[1]))?;
    if pos_1based == 0 {
        bail!("snp position must be >= 1 (1-based), got 0");
    }
    let pos = pos_1based - 1;

    // Parse REF and ALT: either "REF:ALT" (two parts) or "REF>ALT" (one part with >).
    let (ref_str, alt_str) = if parts.len() >= 4 {
        // "snp:chr:pos:REF:ALT"
        (parts[2], parts[3])
    } else if let Some((r, a)) = parts[2].split_once('>') {
        // "snp:chr:pos:REF>ALT"
        (r, a)
    } else {
        bail!("snp spec requires REF:ALT or REF>ALT, got '{}'", parts[2]);
    };

    let ref_allele: Vec<u8> = ref_str.as_bytes().to_vec();
    let alt_allele: Vec<u8> = alt_str.as_bytes().to_vec();

    // Validate: must be non-empty DNA bases.
    if ref_allele.is_empty() || alt_allele.is_empty() {
        bail!("REF and ALT alleles must be non-empty");
    }
    for &b in ref_allele.iter().chain(alt_allele.iter()) {
        if !matches!(b.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T') {
            bail!("invalid base '{}' in allele, expected A/C/G/T", b as char);
        }
    }

    // REF and ALT must differ.
    if ref_allele == alt_allele {
        bail!("REF and ALT alleles are identical: '{}'", ref_str);
    }

    Ok(SimEvent::SmallVariant {
        chrom,
        pos,
        ref_allele: ref_allele.iter().map(|b| b.to_ascii_uppercase()).collect(),
        alt_allele: alt_allele.iter().map(|b| b.to_ascii_uppercase()).collect(),
        gene: "unknown".to_string(),
        allele_fraction: None,
    })
}

/// Parse "exon4" or "exon14" or just "4" into exon number.
fn parse_exon_number(s: &str) -> Result<u32> {
    let s = s.trim();
    let num_str = if s.to_lowercase().starts_with("exon") {
        &s[4..]
    } else {
        s
    };
    num_str
        .parse::<u32>()
        .with_context(|| format!("invalid exon number: '{}'", s))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Synthetic gene fixtures for testing (no real coordinates).
    fn test_genes() -> Vec<GeneTarget> {
        // GENEA on chr1: 8 exons, 200bp each, 500bp spacing.
        let gene_a = GeneTarget {
            gene: "GENEA".to_string(),
            chrom: "chr1".to_string(),
            gene_start: 1000,
            gene_end: 1000 + 7 * 500 + 200, // 4700
            exons: (1..=8)
                .map(|i| {
                    let start = 1000 + (i as u64 - 1) * 500;
                    Exon {
                        chrom: "chr1".to_string(),
                        start,
                        end: start + 200,
                        name: format!("GENEA_exon{}", i),
                        gene: "GENEA".to_string(),
                        number: i,
                    }
                })
                .collect(),
        };
        // GENEB on chr2: 5 exons, 300bp each, 2000bp spacing.
        let gene_b = GeneTarget {
            gene: "GENEB".to_string(),
            chrom: "chr2".to_string(),
            gene_start: 10000,
            gene_end: 10000 + 4 * 2000 + 300, // 18300
            exons: (1..=5)
                .map(|i| {
                    let start = 10000 + (i as u64 - 1) * 2000;
                    Exon {
                        chrom: "chr2".to_string(),
                        start,
                        end: start + 300,
                        name: format!("GENEB_exon{}", i),
                        gene: "GENEB".to_string(),
                        number: i,
                    }
                })
                .collect(),
        };
        vec![gene_a, gene_b]
    }

    #[test]
    fn test_parse_coordinate_deletion() {
        let genes = test_genes();
        let (event, af) = parse_event_spec("del:chr20:30000000-30005000", &genes).unwrap();
        assert!(af.is_none());
        match event {
            SimEvent::Deletion {
                chrom,
                del_start,
                del_end,
                ..
            } => {
                assert_eq!(chrom, "chr20");
                assert_eq!(del_start, 30_000_000);
                assert_eq!(del_end, 30_005_000);
            }
            _ => panic!("expected Deletion"),
        }
    }

    #[test]
    fn test_parse_exon_deletion() {
        let genes = test_genes();
        // GENEA exon 4 start = 1000 + 3*500 = 2500, exon 8 end = 1000 + 7*500 + 200 = 4700
        let (event, af) = parse_event_spec("del:GENEA:exon4-exon8", &genes).unwrap();
        assert!(af.is_none());
        match event {
            SimEvent::Deletion {
                chrom,
                del_start,
                del_end,
                exons,
                ..
            } => {
                assert_eq!(chrom, "chr1");
                assert_eq!(del_start, 2500);
                assert_eq!(del_end, 4700);
                assert_eq!(exons.len(), 5);
            }
            _ => panic!("expected Deletion"),
        }
    }

    #[test]
    fn test_parse_fusion() {
        let genes = test_genes();
        // GENEA exon 3 end = 1000 + 2*500 + 200 = 2200
        // GENEB exon 2 start = 10000 + 1*2000 = 12000
        let (event, af) = parse_event_spec("fusion:GENEA:exon3:GENEB:exon2", &genes).unwrap();
        assert!(af.is_none());
        match event {
            SimEvent::Fusion {
                chrom_a,
                bp_a,
                chrom_b,
                bp_b,
                ..
            } => {
                assert_eq!(chrom_a, "chr1");
                assert_eq!(bp_a, 2200);
                assert_eq!(chrom_b, "chr2");
                assert_eq!(bp_b, 12000);
            }
            _ => panic!("expected Fusion"),
        }
    }

    #[test]
    fn test_parse_with_af_exact() {
        let genes = test_genes();
        let (_, af) = parse_event_spec("del:GENEA:exon4-exon8;af=0.15", &genes).unwrap();
        assert_eq!(af, Some(AfSpec::Exact(0.15)));
    }

    #[test]
    fn test_parse_with_af_het() {
        let genes = test_genes();
        let (_, af) = parse_event_spec("fusion:GENEA:exon3:GENEB:exon2;af=het", &genes).unwrap();
        assert_eq!(af, Some(AfSpec::Het));
    }

    #[test]
    fn test_parse_with_af_hom() {
        let genes = test_genes();
        let (_, af) = parse_event_spec("del:GENEA:exon4-exon8;af=hom", &genes).unwrap();
        assert_eq!(af, Some(AfSpec::Hom));
    }

    #[test]
    fn test_parse_af_invalid() {
        let genes = test_genes();
        assert!(parse_event_spec("del:GENEA:exon4-exon8;af=0.0", &genes).is_err());
        assert!(parse_event_spec("del:GENEA:exon4-exon8;af=1.5", &genes).is_err());
        assert!(parse_event_spec("del:GENEA:exon4-exon8;af=abc", &genes).is_err());
    }

    #[test]
    fn test_parse_snp_colon_format() {
        let genes = test_genes();
        let (event, af) = parse_event_spec("snp:chr1:100:A:T", &genes).unwrap();
        assert!(af.is_none());
        match event {
            SimEvent::SmallVariant {
                chrom,
                pos,
                ref_allele,
                alt_allele,
                ..
            } => {
                assert_eq!(chrom, "chr1");
                assert_eq!(pos, 99); // 1-based input -> 0-based internal
                assert_eq!(ref_allele, b"A");
                assert_eq!(alt_allele, b"T");
            }
            _ => panic!("expected SmallVariant"),
        }
    }

    #[test]
    fn test_parse_snp_arrow_format() {
        let genes = test_genes();
        let (event, _) = parse_event_spec("snp:chr1:100:A>T", &genes).unwrap();
        match event {
            SimEvent::SmallVariant {
                pos,
                ref_allele,
                alt_allele,
                ..
            } => {
                assert_eq!(pos, 99); // 1-based input -> 0-based internal
                assert_eq!(ref_allele, b"A");
                assert_eq!(alt_allele, b"T");
            }
            _ => panic!("expected SmallVariant"),
        }
    }

    #[test]
    fn test_parse_snp_with_af() {
        let genes = test_genes();
        let (_, af) = parse_event_spec("snp:chr1:100:A:T;af=0.3", &genes).unwrap();
        assert_eq!(af, Some(AfSpec::Exact(0.3)));
    }

    #[test]
    fn test_parse_snp_small_del() {
        let genes = test_genes();
        let (event, _) = parse_event_spec("snp:chr1:100:ACG:A", &genes).unwrap();
        match event {
            SimEvent::SmallVariant {
                pos,
                ref_allele,
                alt_allele,
                ..
            } => {
                assert_eq!(pos, 99); // 1-based input -> 0-based internal
                assert_eq!(ref_allele, b"ACG");
                assert_eq!(alt_allele, b"A");
            }
            _ => panic!("expected SmallVariant"),
        }
    }

    #[test]
    fn test_parse_snp_small_ins() {
        let genes = test_genes();
        let (event, _) = parse_event_spec("snp:chr1:100:A:ACGT", &genes).unwrap();
        match event {
            SimEvent::SmallVariant {
                pos,
                ref_allele,
                alt_allele,
                ..
            } => {
                assert_eq!(pos, 99); // 1-based input -> 0-based internal
                assert_eq!(ref_allele, b"A");
                assert_eq!(alt_allele, b"ACGT");
            }
            _ => panic!("expected SmallVariant"),
        }
    }

    #[test]
    fn test_parse_snp_invalid_base() {
        let genes = test_genes();
        assert!(parse_event_spec("snp:chr1:100:A:X", &genes).is_err());
    }

    #[test]
    fn test_parse_snp_same_alleles() {
        let genes = test_genes();
        assert!(parse_event_spec("snp:chr1:100:A:A", &genes).is_err());
    }

    #[test]
    fn test_parse_snp_zero_position_invalid() {
        let genes = test_genes();
        assert!(parse_event_spec("snp:chr1:0:A:T", &genes).is_err());
    }

    #[test]
    fn test_parse_exon_number() {
        assert_eq!(parse_exon_number("exon4").unwrap(), 4);
        assert_eq!(parse_exon_number("exon14").unwrap(), 14);
        assert_eq!(parse_exon_number("4").unwrap(), 4);
        assert!(parse_exon_number("exonXX").is_err());
    }

    #[test]
    fn test_parse_exon_duplication() {
        let genes = test_genes();
        // GENEA exon 4 start = 2500, exon 6 end = 1000 + 5*500 + 200 = 3700
        let (event, af) = parse_event_spec("dup:GENEA:exon4-exon6", &genes).unwrap();
        assert!(af.is_none());
        match event {
            SimEvent::Duplication {
                chrom,
                dup_start,
                dup_end,
                gene,
                ..
            } => {
                assert_eq!(chrom, "chr1");
                assert_eq!(dup_start, 2500);
                assert_eq!(dup_end, 3700);
                assert_eq!(gene, "GENEA");
            }
            _ => panic!("expected Duplication"),
        }
    }

    #[test]
    fn test_parse_exon_inversion() {
        let genes = test_genes();
        // GENEB exon 2 start = 12000, exon 4 end = 10000 + 3*2000 + 300 = 16300
        let (event, af) = parse_event_spec("inv:GENEB:exon2-exon4", &genes).unwrap();
        assert!(af.is_none());
        match event {
            SimEvent::Inversion {
                chrom,
                inv_start,
                inv_end,
                gene,
                ..
            } => {
                assert_eq!(chrom, "chr2");
                assert_eq!(inv_start, 12000);
                assert_eq!(inv_end, 16300);
                assert_eq!(gene, "GENEB");
            }
            _ => panic!("expected Inversion"),
        }
    }

    #[test]
    fn test_parse_coord_duplication() {
        let genes = test_genes();
        let (event, _) = parse_event_spec("dup:chr20:30000000-30005000", &genes).unwrap();
        match event {
            SimEvent::Duplication {
                chrom,
                dup_start,
                dup_end,
                ..
            } => {
                assert_eq!(chrom, "chr20");
                assert_eq!(dup_start, 30_000_000);
                assert_eq!(dup_end, 30_005_000);
            }
            _ => panic!("expected Duplication"),
        }
    }

    #[test]
    fn test_parse_coord_inversion() {
        let genes = test_genes();
        let (event, _) = parse_event_spec("inv:chr20:30000000-30005000", &genes).unwrap();
        match event {
            SimEvent::Inversion {
                chrom,
                inv_start,
                inv_end,
                ..
            } => {
                assert_eq!(chrom, "chr20");
                assert_eq!(inv_start, 30_000_000);
                assert_eq!(inv_end, 30_005_000);
            }
            _ => panic!("expected Inversion"),
        }
    }

    #[test]
    fn test_parse_explicit_insertion_sequence() {
        let genes = test_genes();
        let (event, _) = parse_event_spec("ins:chr20:30000000:ACGTACGT", &genes).unwrap();
        match event {
            SimEvent::Insertion {
                chrom,
                pos,
                ins_seq,
                ins_len,
                ..
            } => {
                assert_eq!(chrom, "chr20");
                assert_eq!(pos, 30_000_000);
                assert_eq!(ins_seq, Some(b"ACGTACGT".to_vec()));
                assert_eq!(ins_len, 8);
            }
            _ => panic!("expected Insertion"),
        }
    }

    #[test]
    fn test_parse_random_insertion_length() {
        let genes = test_genes();
        let (event, _) = parse_event_spec("ins:chr20:30000000:500", &genes).unwrap();
        match event {
            SimEvent::Insertion {
                ins_seq, ins_len, ..
            } => {
                assert!(ins_seq.is_none());
                assert_eq!(ins_len, 500);
            }
            _ => panic!("expected Insertion"),
        }
    }

    #[test]
    fn test_parse_insertion_invalid_sequence() {
        let genes = test_genes();
        assert!(parse_event_spec("ins:chr20:30000000:ACGXYZ", &genes).is_err());
    }

    #[test]
    fn test_parse_inverted_fusion() {
        let genes = test_genes();
        let (event, _) = parse_event_spec("fusion:GENEA:exon3:GENEB:exon2:inv", &genes).unwrap();
        match event {
            SimEvent::Fusion { inverted, .. } => {
                assert!(inverted, "fusion with :inv suffix should be inverted");
            }
            _ => panic!("expected Fusion"),
        }
    }

    #[test]
    fn test_parse_forward_fusion() {
        let genes = test_genes();
        let (event, _) = parse_event_spec("fusion:GENEA:exon3:GENEB:exon2", &genes).unwrap();
        match event {
            SimEvent::Fusion { inverted, .. } => {
                assert!(!inverted, "fusion without :inv suffix should be forward");
            }
            _ => panic!("expected Fusion"),
        }
    }
}
