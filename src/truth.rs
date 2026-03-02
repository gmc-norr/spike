//! Truth VCF writer for simulated events.
//!
//! Follows the same coordinate conventions as sv_exon:
//! - VCF POS: 1-based preceding-base (numerically == 0-based start)
//! - VCF END: 1-based inclusive

use anyhow::{Context, Result};
use std::io::Write;

use crate::types::SimEvent;

/// Write a truth VCF describing the simulated events.
///
/// Each event may carry a per-event allele fraction; `default_af` is used as fallback.
pub fn write_truth_vcf(
    events: &[SimEvent],
    default_af: f64,
    output_path: &str,
    ref_path: &str,
) -> Result<()> {
    let mut f = std::fs::File::create(output_path)
        .with_context(|| format!("failed to create truth VCF: {}", output_path))?;

    // Header.
    writeln!(f, "##fileformat=VCFv4.3")?;
    writeln!(f, "##fileDate={}", chrono_date())?;
    writeln!(f, "##source=spike")?;
    writeln!(f, "##reference={}", ref_path)?;
    writeln!(f, "##ALT=<ID=DEL,Description=\"Deletion\">")?;
    writeln!(f, "##ALT=<ID=DUP,Description=\"Tandem duplication\">")?;
    writeln!(f, "##ALT=<ID=INV,Description=\"Inversion\">")?;
    writeln!(f, "##ALT=<ID=INS,Description=\"Insertion\">")?;
    writeln!(
        f,
        "##ALT=<ID=BND,Description=\"Translocation breakend\">"
    )?;
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
        "##INFO=<ID=SIM_VAF,Number=1,Type=Float,Description=\"Simulated allele fraction\">"
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
    writeln!(
        f,
        "##FILTER=<ID=PASS,Description=\"All filters passed\">"
    )?;
    writeln!(
        f,
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
    )?;
    writeln!(
        f,
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE"
    )?;

    for (i, event) in events.iter().enumerate() {
        let event_af = event.allele_fraction().unwrap_or(default_af);
        let gt = genotype_from_vaf(event_af);

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

                writeln!(
                    f,
                    "{}\t{}\tsim_del_{}\tN\t<DEL>\t999\tPASS\tSVTYPE=DEL;END={};SVLEN=-{};SIM_VAF={:.3};SIM_GENE={};SIM_EXONS={}\tGT\t{}",
                    chrom,
                    del_start, // VCF POS: numerically == 0-based start
                    i + 1,
                    del_end, // VCF END: 0-based half-open value == 1-based inclusive
                    sv_len,
                    event_af,
                    gene,
                    exons_str,
                    gt,
                )?;
            }
            SimEvent::Fusion {
                chrom_a,
                bp_a,
                gene_a,
                chrom_b,
                bp_b,
                gene_b,
                inverted,
                ..
            } => {
                let id_a = format!("sim_fus_{}", i + 1);
                let id_b = format!("sim_fus_{}_mate", i + 1);

                // BND POS is 1-based breakpoint position (unlike DEL/DUP/INV/INS
                // where POS is 1-based preceding base == 0-based start).
                let bp_a_1based = bp_a + 1;
                let bp_b_1based = bp_b + 1;

                // BND bracket notation encodes orientation:
                //   Forward:  N[chr:pos[  and  ]chr:pos]N
                //   Inverted: N]chr:pos]  and  [chr:pos[N
                let (alt_a, alt_b) = if *inverted {
                    (
                        format!("N]{}:{}]", chrom_b, bp_b_1based),
                        format!("[{}:{}[N", chrom_a, bp_a_1based),
                    )
                } else {
                    (
                        format!("N[{}:{}[", chrom_b, bp_b_1based),
                        format!("]{}:{}]N", chrom_a, bp_a_1based),
                    )
                };
                writeln!(
                    f,
                    "{}\t{}\t{}\tN\t{}\t999\tPASS\tSVTYPE=BND;MATEID={};SIM_VAF={:.3};SIM_GENE={}\tGT\t{}",
                    chrom_a, bp_a_1based, id_a, alt_a, id_b, event_af, gene_a, gt,
                )?;

                writeln!(
                    f,
                    "{}\t{}\t{}\tN\t{}\t999\tPASS\tSVTYPE=BND;MATEID={};SIM_VAF={:.3};SIM_GENE={}\tGT\t{}",
                    chrom_b, bp_b_1based, id_b, alt_b, id_a, event_af, gene_b, gt,
                )?;
            }
            SimEvent::Duplication {
                chrom,
                dup_start,
                dup_end,
                gene,
                ..
            } => {
                let sv_len = *dup_end as i64 - *dup_start as i64;
                writeln!(
                    f,
                    "{}\t{}\tsim_dup_{}\tN\t<DUP>\t999\tPASS\tSVTYPE=DUP;END={};SVLEN={};SIM_VAF={:.3};SIM_GENE={}\tGT\t{}",
                    chrom, dup_start, i + 1, dup_end, sv_len, event_af, gene, gt,
                )?;
            }
            SimEvent::Inversion {
                chrom,
                inv_start,
                inv_end,
                gene,
                ..
            } => {
                let sv_len = *inv_end as i64 - *inv_start as i64;
                writeln!(
                    f,
                    "{}\t{}\tsim_inv_{}\tN\t<INV>\t999\tPASS\tSVTYPE=INV;END={};SVLEN={};SIM_VAF={:.3};SIM_GENE={}\tGT\t{}",
                    chrom, inv_start, i + 1, inv_end, sv_len, event_af, gene, gt,
                )?;
            }
            SimEvent::Insertion {
                chrom,
                pos,
                ins_len,
                gene,
                ..
            } => {
                writeln!(
                    f,
                    "{}\t{}\tsim_ins_{}\tN\t<INS>\t999\tPASS\tSVTYPE=INS;SVLEN={};SIM_VAF={:.3};SIM_GENE={}\tGT\t{}",
                    chrom, pos, i + 1, ins_len, event_af, gene, gt,
                )?;
            }
            SimEvent::SmallVariant {
                chrom,
                pos,
                ref_allele,
                alt_allele,
                gene,
                ..
            } => {
                let ref_str = String::from_utf8_lossy(ref_allele);
                let alt_str = String::from_utf8_lossy(alt_allele);
                // VCF POS is 1-based.
                writeln!(
                    f,
                    "{}\t{}\tsim_var_{}\t{}\t{}\t999\tPASS\tSIM_VAF={:.3};SIM_GENE={}\tGT\t{}",
                    chrom,
                    pos + 1, // 0-based → 1-based
                    i + 1,
                    ref_str,
                    alt_str,
                    event_af,
                    gene,
                    gt,
                )?;
            }
        }
    }

    log::info!("Wrote truth VCF with {} events to {}", events.len(), output_path);
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
        .unwrap()
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
        31, 30, 31, 30, 31, 31, 30, 31, 30, 31,
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
