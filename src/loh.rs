//! Haplotype-aware simulation for heterozygous SVs.
//!
//! **Deletions (LOH)**: one haplotype is removed. Het SNPs within the deletion
//! become homozygous (only the surviving allele remains). Reads from the deleted
//! haplotype are suppressed.
//!
//! **Duplications (allelic imbalance)**: one haplotype is duplicated. Het SNPs
//! shift from 50/50 to ~33/67 (2:1 ratio). Depth copies are drawn preferentially
//! from the duplicated haplotype.
//!
//! Two strategies for finding het SNP positions:
//! 1. **Pileup** (default): manual pileup to find positions with ~50/50 allele split
//! 2. **gVCF** (optional): extract het SNPs from a pre-called VCF (e.g., DeepVariant)

use std::cmp::Reverse;
use std::collections::{HashMap, HashSet};

use anyhow::{Context, Result};
use noodles::sam::alignment::record::cigar::op::Kind;
use noodles::sam::alignment::record::Cigar as CigarTrait;
use rand::rngs::StdRng;
use rand::Rng;

/// Return type for haplotype identification: (target_set, classified_set, haplotype_variants).
type HaplotypeResult = (HashSet<String>, HashSet<String>, HashMap<u64, u8>);

/// Return type for pileup: (het_snps, optional per-read alleles).
type PileupResult = (Vec<HetSnp>, Option<HashMap<String, Vec<(u64, u8)>>>);

/// A heterozygous SNP position with its two alleles.
#[derive(Debug)]
struct HetSnp {
    pos: u64,    // 0-based reference position
    allele1: u8, // major allele (A/C/G/T)
    allele2: u8, // minor allele
}

/// Identify reads from one haplotype in a genomic region.
///
/// Finds het SNPs in [region_start, region_end), randomly picks one allele
/// as the "target" at each position, then classifies reads by which allele
/// they carry.
///
/// Returns `(target_set, classified_set, haplotype_variants)`:
/// - `target_set`: reads predominantly carrying the target allele.
/// - `classified_set`: all reads that overlapped at least one het SNP.
/// - `haplotype_variants`: position → target allele for variant substitution.
///
/// If no het SNPs are found, all sets are empty and the caller should fall
/// back to random suppression.
#[allow(clippy::too_many_arguments)]
fn identify_haplotype_reads(
    alignment_path: &str,
    chrom: &str,
    region_start: u64,
    region_end: u64,
    min_mapq: u8,
    gvcf_path: Option<&str>,
    ref_path: Option<&str>,
    label: &str,
    rng: &mut StdRng,
) -> Result<HaplotypeResult> {
    // Strategy 1: gVCF provides het SNP positions → classify reads from BAM (1 BAM pass).
    // Strategy 2: pileup finds het SNPs AND records per-read alleles (1 BAM pass).
    let (het_snps, read_alleles) = if let Some(gvcf) = gvcf_path {
        let snps = load_het_snps_from_gvcf(gvcf, chrom, region_start, region_end)?;
        if snps.is_empty() {
            log::info!("{}: no het SNPs from gVCF, trying pileup fallback", label);
            pileup_and_collect(
                alignment_path,
                chrom,
                region_start,
                region_end,
                min_mapq,
                ref_path,
            )?
        } else {
            // gVCF gave us het SNPs; we still need one BAM/CRAM pass to classify reads.
            (snps, None)
        }
    } else {
        pileup_and_collect(
            alignment_path,
            chrom,
            region_start,
            region_end,
            min_mapq,
            ref_path,
        )?
    };

    if het_snps.is_empty() {
        log::info!(
            "{}: no het SNPs found in {}:{}-{}, falling back to random",
            label,
            chrom,
            region_start,
            region_end,
        );
        return Ok((HashSet::new(), HashSet::new(), HashMap::new()));
    }

    log::info!(
        "{}: found {} het SNPs in {}:{}-{}",
        label,
        het_snps.len(),
        chrom,
        region_start,
        region_end,
    );

    // Randomly pick one allele as the "target" at each het position.
    let mut target_allele: HashMap<u64, u8> = HashMap::new();
    for snp in &het_snps {
        let target = if rng.gen::<bool>() {
            snp.allele1
        } else {
            snp.allele2
        };
        target_allele.insert(snp.pos, target);
    }

    // Classify reads: use pre-collected alleles if available, otherwise BAM/CRAM pass.
    if let Some(read_alleles) = read_alleles {
        let (target_set, classified_set) =
            classify_from_collected(&target_allele, &read_alleles, label);
        Ok((target_set, classified_set, target_allele))
    } else {
        let (target_set, classified_set) = classify_reads_by_haplotype(
            alignment_path,
            chrom,
            region_start,
            region_end,
            min_mapq,
            &target_allele,
            ref_path,
            label,
        )?;
        Ok((target_set, classified_set, target_allele))
    }
}

/// Identify reads from the "deleted haplotype" in a deletion region.
///
/// Returns `(target_set, classified_set)`:
/// - `target_set`: read names that should be suppressed to simulate LOH.
/// - `classified_set`: all reads that overlapped at least one het SNP.
///   Reads NOT in `classified_set` couldn't be assigned and should fall back
///   to random suppression at the VAF rate.
///
/// If no het SNPs are found, both sets are empty.
#[allow(clippy::too_many_arguments)]
pub fn identify_deleted_haplotype_reads(
    alignment_path: &str,
    chrom: &str,
    del_start: u64,
    del_end: u64,
    min_mapq: u8,
    gvcf_path: Option<&str>,
    ref_path: Option<&str>,
    rng: &mut StdRng,
) -> Result<(HashSet<String>, HashSet<String>)> {
    let (target_set, classified_set, _variants) = identify_haplotype_reads(
        alignment_path,
        chrom,
        del_start,
        del_end,
        min_mapq,
        gvcf_path,
        ref_path,
        "LOH-DEL",
        rng,
    )?;
    Ok((target_set, classified_set))
}

/// Identify reads from the "duplicated haplotype" in a duplication region.
///
/// Returns `(target_set, classified_set, haplotype_variants)`:
/// - `target_set`: reads from the haplotype that should be preferentially
///   copied during depth increase, to create realistic allelic imbalance.
/// - `classified_set`: all reads that overlapped at least one het SNP.
///   Reads NOT in `classified_set` should fall back to random copying at VAF rate.
/// - `haplotype_variants`: position → allele for het SNPs on the duplicated
///   haplotype. Synthetic depth-copy reads should carry these alleles to
///   produce correct BAF (2:1 ratio at het sites).
#[allow(clippy::too_many_arguments)]
pub fn identify_duplicated_haplotype_reads(
    alignment_path: &str,
    chrom: &str,
    dup_start: u64,
    dup_end: u64,
    min_mapq: u8,
    gvcf_path: Option<&str>,
    ref_path: Option<&str>,
    rng: &mut StdRng,
) -> Result<HaplotypeResult> {
    identify_haplotype_reads(
        alignment_path,
        chrom,
        dup_start,
        dup_end,
        min_mapq,
        gvcf_path,
        ref_path,
        "AI-DUP",
        rng,
    )
}

/// Load het SNP positions from a VCF/gVCF file.
///
/// For `.vcf.gz` files, uses `bcftools view -H -r region` for efficient
/// indexed access. For plain `.vcf` files, reads and filters line by line.
fn load_het_snps_from_gvcf(
    gvcf_path: &str,
    chrom: &str,
    region_start: u64,
    region_end: u64,
) -> Result<Vec<HetSnp>> {
    let lines: Vec<String> = if gvcf_path.ends_with(".gz") {
        // Use bcftools for indexed access to bgzipped VCF.
        let region = format!("{}:{}-{}", chrom, region_start + 1, region_end);
        let output = std::process::Command::new("bcftools")
            .args(["view", "-H", "-r", &region, gvcf_path])
            .output()
            .context("failed to run bcftools for gVCF reading (is bcftools in PATH?)")?;

        if !output.status.success() {
            let stderr = String::from_utf8_lossy(&output.stderr);
            anyhow::bail!("bcftools failed: {}", stderr);
        }

        String::from_utf8_lossy(&output.stdout)
            .lines()
            .map(|l| l.to_string())
            .collect()
    } else {
        // Plain text VCF: read and filter.
        let content = std::fs::read_to_string(gvcf_path)
            .with_context(|| format!("failed to read gVCF: {}", gvcf_path))?;
        content.lines().map(|l| l.to_string()).collect()
    };

    let mut het_snps = Vec::new();

    for line in &lines {
        if line.starts_with('#') {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 10 {
            continue;
        }

        let line_chrom = fields[0];
        if line_chrom != chrom {
            continue;
        }

        // VCF POS is 1-based.
        let pos: u64 = match fields[1].parse::<u64>() {
            Ok(p) => p.saturating_sub(1), // convert to 0-based
            Err(_) => continue,
        };

        if pos < region_start || pos >= region_end {
            continue;
        }

        let ref_allele = fields[3].as_bytes();
        let alt_field = fields[4];

        // Handle multi-allelic: take the first ALT allele.
        let alt_allele = alt_field.split(',').next().unwrap_or(".").as_bytes();

        // Only consider SNPs (single-base REF and ALT).
        if ref_allele.len() != 1 || alt_allele.len() != 1 {
            continue;
        }
        if alt_allele == b"." || alt_allele == b"*" {
            continue;
        }

        // Check genotype for heterozygosity.
        let gt_field = fields[9].split(':').next().unwrap_or("");
        let is_het =
            gt_field == "0/1" || gt_field == "1/0" || gt_field == "0|1" || gt_field == "1|0";

        if is_het {
            het_snps.push(HetSnp {
                pos,
                allele1: ref_allele[0].to_ascii_uppercase(),
                allele2: alt_allele[0].to_ascii_uppercase(),
            });
        }
    }

    het_snps.sort_by_key(|s| s.pos);
    Ok(het_snps)
}

/// Single-pass pileup: find het SNPs AND collect per-read alleles.
///
/// Returns (het_snps, per_read_alleles). The per-read alleles map each read name
/// to its observed bases at every position in the region, which is later filtered
/// to het SNP positions during classification.
fn pileup_and_collect(
    alignment_path: &str,
    chrom: &str,
    region_start: u64,
    region_end: u64,
    min_mapq: u8,
    ref_path: Option<&str>,
) -> Result<PileupResult> {
    if crate::extract::is_cram(alignment_path) {
        let rp =
            ref_path.ok_or_else(|| anyhow::anyhow!("CRAM input requires a reference FASTA"))?;
        pileup_and_collect_cram(
            alignment_path,
            chrom,
            region_start,
            region_end,
            min_mapq,
            rp,
        )
    } else {
        pileup_and_collect_bam(alignment_path, chrom, region_start, region_end, min_mapq)
    }
}

/// BAM-specific pileup and allele collection.
fn pileup_and_collect_bam(
    bam_path: &str,
    chrom: &str,
    region_start: u64,
    region_end: u64,
    min_mapq: u8,
) -> Result<PileupResult> {
    let mut allele_counts: HashMap<u64, [u32; 4]> = HashMap::new();
    let mut read_alleles: HashMap<String, Vec<(u64, u8)>> = HashMap::new();

    let mut reader = noodles::bam::io::indexed_reader::Builder::default()
        .build_from_path(bam_path)
        .with_context(|| format!("failed to open BAM for pileup: {}", bam_path))?;
    let header = reader.read_header()?;

    let start_pos = crate::extract::safe_noodles_position(region_start + 1);
    let end_pos = crate::extract::safe_noodles_position(region_end);
    let region = noodles::core::Region::new(chrom, start_pos..=end_pos);
    let query = reader.query(&header, &region)?;

    for rec_result in query {
        let record = rec_result?;
        let flags = record.flags();

        if flags.is_unmapped()
            || flags.is_secondary()
            || flags.is_supplementary()
            || flags.is_duplicate()
            || flags.is_qc_fail()
        {
            continue;
        }

        let mq: u8 = match record.mapping_quality() {
            Some(q) => u8::from(q),
            None => 0,
        };
        if mq < min_mapq {
            continue;
        }

        let name = match record.name() {
            Some(n) => String::from_utf8_lossy(n.as_ref()).into_owned(),
            None => continue,
        };

        let align_start = match record.alignment_start() {
            Some(Ok(p)) => usize::from(p).saturating_sub(1) as u64,
            _ => continue,
        };

        let seq: Vec<u8> = record.sequence().iter().collect();
        let cigar = record.cigar();

        walk_cigar_pileup(
            &seq,
            Box::new(cigar.iter()),
            align_start,
            region_start,
            region_end,
            &name,
            &mut allele_counts,
            &mut read_alleles,
        );
    }

    Ok(find_het_snps(allele_counts, Some(read_alleles)))
}

/// CRAM-specific pileup and allele collection.
fn pileup_and_collect_cram(
    cram_path: &str,
    chrom: &str,
    region_start: u64,
    region_end: u64,
    min_mapq: u8,
    ref_path: &str,
) -> Result<PileupResult> {
    let mut allele_counts: HashMap<u64, [u32; 4]> = HashMap::new();
    let mut read_alleles: HashMap<String, Vec<(u64, u8)>> = HashMap::new();

    let repository = crate::extract::build_fasta_repository(ref_path)?;

    let mut reader = noodles::cram::io::indexed_reader::Builder::default()
        .set_reference_sequence_repository(repository)
        .build_from_path(cram_path)
        .with_context(|| format!("failed to open CRAM for pileup: {}", cram_path))?;
    let header = reader.read_header()?;

    let start_pos = crate::extract::safe_noodles_position(region_start + 1);
    let end_pos = crate::extract::safe_noodles_position(region_end);
    let region = noodles::core::Region::new(chrom, start_pos..=end_pos);
    let query = reader.query(&header, &region)?;

    for rec_result in query {
        let cram_record = rec_result?;
        let buf = cram_record
            .try_into_alignment_record(&header)
            .with_context(|| "failed to convert CRAM record")?;

        let flags = buf.flags();

        if flags.is_unmapped()
            || flags.is_secondary()
            || flags.is_supplementary()
            || flags.is_duplicate()
            || flags.is_qc_fail()
        {
            continue;
        }

        let mq: u8 = match buf.mapping_quality() {
            Some(q) => u8::from(q),
            None => 0,
        };
        if mq < min_mapq {
            continue;
        }

        let name = match buf.name() {
            Some(n) => String::from_utf8_lossy(n.as_ref()).into_owned(),
            None => continue,
        };

        let align_start = match buf.alignment_start() {
            Some(p) => usize::from(p).saturating_sub(1) as u64,
            None => continue,
        };

        let seq: Vec<u8> = buf.sequence().as_ref().to_vec();
        let cigar = buf.cigar();

        walk_cigar_pileup(
            &seq,
            CigarTrait::iter(&cigar),
            align_start,
            region_start,
            region_end,
            &name,
            &mut allele_counts,
            &mut read_alleles,
        );
    }

    Ok(find_het_snps(allele_counts, Some(read_alleles)))
}

/// Walk CIGAR operations and collect allele counts + per-read alleles.
/// Shared between BAM and CRAM pileup paths.
#[allow(clippy::too_many_arguments)]
fn walk_cigar_pileup(
    seq: &[u8],
    cigar_ops: Box<
        dyn Iterator<Item = std::io::Result<noodles::sam::alignment::record::cigar::Op>> + '_,
    >,
    align_start: u64,
    region_start: u64,
    region_end: u64,
    name: &str,
    allele_counts: &mut HashMap<u64, [u32; 4]>,
    read_alleles: &mut HashMap<String, Vec<(u64, u8)>>,
) {
    let mut ref_pos = align_start;
    let mut seq_pos = 0usize;

    for op_result in cigar_ops {
        let op = match op_result {
            Ok(o) => o,
            Err(_) => break,
        };
        let len = op.len();

        match op.kind() {
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                for i in 0..len {
                    let rp = ref_pos + i as u64;
                    if rp >= region_start && rp < region_end {
                        if let Some(&base) = seq.get(seq_pos + i) {
                            let base = base.to_ascii_uppercase();
                            let idx = match base {
                                b'A' => Some(0),
                                b'C' => Some(1),
                                b'G' => Some(2),
                                b'T' => Some(3),
                                _ => None,
                            };
                            if let Some(idx) = idx {
                                allele_counts.entry(rp).or_insert([0; 4])[idx] += 1;
                                read_alleles
                                    .entry(name.to_string())
                                    .or_default()
                                    .push((rp, base));
                            }
                        }
                    }
                }
                ref_pos += len as u64;
                seq_pos += len;
            }
            Kind::Insertion | Kind::SoftClip => {
                seq_pos += len;
            }
            Kind::Deletion | Kind::Skip => {
                ref_pos += len as u64;
            }
            Kind::HardClip | Kind::Pad => {}
        }
    }
}

/// Find het SNP positions from allele counts.
/// Returns (het_snps, read_alleles).
fn find_het_snps(
    allele_counts: HashMap<u64, [u32; 4]>,
    read_alleles: Option<HashMap<String, Vec<(u64, u8)>>>,
) -> PileupResult {
    let mut het_snps = Vec::new();

    for (&pos, counts) in &allele_counts {
        let total: u32 = counts.iter().sum();
        if total < 10 {
            continue;
        }

        let mut sorted: Vec<(u8, u32)> = [
            (b'A', counts[0]),
            (b'C', counts[1]),
            (b'G', counts[2]),
            (b'T', counts[3]),
        ]
        .iter()
        .filter(|(_, c)| *c > 0)
        .copied()
        .collect();

        sorted.sort_by_key(|&(_, c)| Reverse(c));

        if sorted.len() >= 2 {
            let f1 = sorted[0].1 as f64 / total as f64;
            let f2 = sorted[1].1 as f64 / total as f64;

            if (0.2..=0.8).contains(&f1) && (0.2..=0.8).contains(&f2) {
                het_snps.push(HetSnp {
                    pos,
                    allele1: sorted[0].0,
                    allele2: sorted[1].0,
                });
            }
        }
    }

    het_snps.sort_by_key(|s| s.pos);
    (het_snps, read_alleles)
}

/// Classify reads using pre-collected per-read alleles (no extra BAM pass).
///
/// Returns `(target_set, classified_set)`:
/// - `target_set`: reads that predominantly carry the target allele (to suppress/copy).
/// - `classified_set`: all reads that overlapped at least one het SNP (classifiable).
///   Reads NOT in `classified_set` couldn't be assigned to either haplotype and
///   should fall back to random suppression/copying at the VAF rate.
fn classify_from_collected(
    target_allele: &HashMap<u64, u8>,
    read_alleles: &HashMap<String, Vec<(u64, u8)>>,
    label: &str,
) -> (HashSet<String>, HashSet<String>) {
    let mut target_set = HashSet::new();
    let mut classified_set = HashSet::new();
    let mut n_ambiguous = 0u32;

    for (name, alleles) in read_alleles {
        let mut target_count = 0u32;
        let mut other_count = 0u32;

        for &(pos, base) in alleles {
            if let Some(&tgt) = target_allele.get(&pos) {
                if base == tgt {
                    target_count += 1;
                } else {
                    other_count += 1;
                }
            }
        }

        if target_count == 0 && other_count == 0 {
            continue;
        }
        classified_set.insert(name.clone());

        if target_count > other_count {
            target_set.insert(name.clone());
        } else if target_count == other_count && target_count > 0 {
            n_ambiguous += 1;
        }
    }

    let other = classified_set
        .len()
        .saturating_sub(target_set.len())
        .saturating_sub(n_ambiguous as usize);
    log::info!(
        "{}: classified {} fragments ({} target hap, {} other, {} ambiguous ties)",
        label,
        classified_set.len(),
        target_set.len(),
        other,
        n_ambiguous,
    );

    (target_set, classified_set)
}

/// Classify reads by which haplotype they belong to.
///
/// For each read overlapping the region, count how many het SNP positions
/// match the "target" allele vs the other allele. Reads with more target
/// allele matches are returned in the target set.
///
/// Returns `(target_set, classified_set)`:
/// - `target_set`: reads predominantly carrying the target allele.
/// - `classified_set`: all reads that overlapped at least one het SNP.
#[allow(clippy::too_many_arguments)]
fn classify_reads_by_haplotype(
    alignment_path: &str,
    chrom: &str,
    region_start: u64,
    region_end: u64,
    min_mapq: u8,
    target_allele: &HashMap<u64, u8>,
    ref_path: Option<&str>,
    label: &str,
) -> Result<(HashSet<String>, HashSet<String>)> {
    if crate::extract::is_cram(alignment_path) {
        let rp =
            ref_path.ok_or_else(|| anyhow::anyhow!("CRAM input requires a reference FASTA"))?;
        classify_reads_by_haplotype_cram(
            alignment_path,
            chrom,
            region_start,
            region_end,
            min_mapq,
            target_allele,
            rp,
            label,
        )
    } else {
        classify_reads_by_haplotype_bam(
            alignment_path,
            chrom,
            region_start,
            region_end,
            min_mapq,
            target_allele,
            label,
        )
    }
}

/// BAM-specific haplotype read classification.
#[allow(clippy::too_many_arguments)]
fn classify_reads_by_haplotype_bam(
    bam_path: &str,
    chrom: &str,
    region_start: u64,
    region_end: u64,
    min_mapq: u8,
    target_allele: &HashMap<u64, u8>,
    label: &str,
) -> Result<(HashSet<String>, HashSet<String>)> {
    let mut read_scores: HashMap<String, (u32, u32)> = HashMap::new();

    let mut reader = noodles::bam::io::indexed_reader::Builder::default()
        .build_from_path(bam_path)
        .with_context(|| {
            format!(
                "failed to open BAM for haplotype classification: {}",
                bam_path
            )
        })?;
    let header = reader.read_header()?;

    let start_pos = crate::extract::safe_noodles_position(region_start + 1);
    let end_pos = crate::extract::safe_noodles_position(region_end);
    let region = noodles::core::Region::new(chrom, start_pos..=end_pos);
    let query = reader.query(&header, &region)?;

    for rec_result in query {
        let record = rec_result?;
        let flags = record.flags();

        if flags.is_unmapped()
            || flags.is_secondary()
            || flags.is_supplementary()
            || flags.is_duplicate()
            || flags.is_qc_fail()
        {
            continue;
        }

        let mq: u8 = match record.mapping_quality() {
            Some(q) => u8::from(q),
            None => 0,
        };
        if mq < min_mapq {
            continue;
        }

        let name = match record.name() {
            Some(n) => String::from_utf8_lossy(n.as_ref()).into_owned(),
            None => continue,
        };

        let align_start = match record.alignment_start() {
            Some(Ok(p)) => usize::from(p).saturating_sub(1) as u64,
            _ => continue,
        };

        let seq: Vec<u8> = record.sequence().iter().collect();
        let cigar = record.cigar();

        let scores = read_scores.entry(name).or_insert((0, 0));
        walk_cigar_classify(
            &seq,
            Box::new(cigar.iter()),
            align_start,
            target_allele,
            scores,
        );
    }

    Ok(build_classification_sets(read_scores, label))
}

/// CRAM-specific haplotype read classification.
#[allow(clippy::too_many_arguments)]
fn classify_reads_by_haplotype_cram(
    cram_path: &str,
    chrom: &str,
    region_start: u64,
    region_end: u64,
    min_mapq: u8,
    target_allele: &HashMap<u64, u8>,
    ref_path: &str,
    label: &str,
) -> Result<(HashSet<String>, HashSet<String>)> {
    let mut read_scores: HashMap<String, (u32, u32)> = HashMap::new();

    let repository = crate::extract::build_fasta_repository(ref_path)?;

    let mut reader = noodles::cram::io::indexed_reader::Builder::default()
        .set_reference_sequence_repository(repository)
        .build_from_path(cram_path)
        .with_context(|| {
            format!(
                "failed to open CRAM for haplotype classification: {}",
                cram_path
            )
        })?;
    let header = reader.read_header()?;

    let start_pos = crate::extract::safe_noodles_position(region_start + 1);
    let end_pos = crate::extract::safe_noodles_position(region_end);
    let region = noodles::core::Region::new(chrom, start_pos..=end_pos);
    let query = reader.query(&header, &region)?;

    for rec_result in query {
        let cram_record = rec_result?;
        let buf = cram_record
            .try_into_alignment_record(&header)
            .with_context(|| "failed to convert CRAM record")?;

        let flags = buf.flags();

        if flags.is_unmapped()
            || flags.is_secondary()
            || flags.is_supplementary()
            || flags.is_duplicate()
            || flags.is_qc_fail()
        {
            continue;
        }

        let mq: u8 = match buf.mapping_quality() {
            Some(q) => u8::from(q),
            None => 0,
        };
        if mq < min_mapq {
            continue;
        }

        let name = match buf.name() {
            Some(n) => String::from_utf8_lossy(n.as_ref()).into_owned(),
            None => continue,
        };

        let align_start = match buf.alignment_start() {
            Some(p) => usize::from(p).saturating_sub(1) as u64,
            None => continue,
        };

        let seq: Vec<u8> = buf.sequence().as_ref().to_vec();
        let cigar = buf.cigar();

        let scores = read_scores.entry(name).or_insert((0, 0));
        walk_cigar_classify(
            &seq,
            CigarTrait::iter(&cigar),
            align_start,
            target_allele,
            scores,
        );
    }

    Ok(build_classification_sets(read_scores, label))
}

/// Walk CIGAR operations and score alleles against target haplotype.
/// Shared between BAM and CRAM classification paths.
fn walk_cigar_classify(
    seq: &[u8],
    cigar_ops: Box<
        dyn Iterator<Item = std::io::Result<noodles::sam::alignment::record::cigar::Op>> + '_,
    >,
    align_start: u64,
    target_allele: &HashMap<u64, u8>,
    scores: &mut (u32, u32),
) {
    let mut ref_pos = align_start;
    let mut seq_pos = 0usize;

    for op_result in cigar_ops {
        let op = match op_result {
            Ok(o) => o,
            Err(_) => break,
        };
        let len = op.len();

        match op.kind() {
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                for i in 0..len {
                    let rp = ref_pos + i as u64;
                    if let Some(&tgt_allele) = target_allele.get(&rp) {
                        if let Some(&base) = seq.get(seq_pos + i) {
                            let base = base.to_ascii_uppercase();
                            if base == tgt_allele {
                                scores.0 += 1;
                            } else if base != b'N' {
                                scores.1 += 1;
                            }
                        }
                    }
                }
                ref_pos += len as u64;
                seq_pos += len;
            }
            Kind::Insertion | Kind::SoftClip => {
                seq_pos += len;
            }
            Kind::Deletion | Kind::Skip => {
                ref_pos += len as u64;
            }
            Kind::HardClip | Kind::Pad => {}
        }
    }
}

/// Build target and classified sets from read scores.
/// Shared between BAM and CRAM classification paths.
fn build_classification_sets(
    read_scores: HashMap<String, (u32, u32)>,
    label: &str,
) -> (HashSet<String>, HashSet<String>) {
    let mut target_set = HashSet::new();
    let mut classified_set = HashSet::new();
    let mut n_ambiguous = 0u32;

    for (name, (target_count, other_count)) in &read_scores {
        if *target_count == 0 && *other_count == 0 {
            continue;
        }
        classified_set.insert(name.clone());

        if *target_count > *other_count {
            target_set.insert(name.clone());
        } else if *target_count == *other_count && *target_count > 0 {
            n_ambiguous += 1;
        }
    }

    let other = classified_set
        .len()
        .saturating_sub(target_set.len())
        .saturating_sub(n_ambiguous as usize);
    log::info!(
        "{}: classified {} fragments ({} target hap, {} other, {} ambiguous ties)",
        label,
        classified_set.len(),
        target_set.len(),
        other,
        n_ambiguous,
    );

    (target_set, classified_set)
}
