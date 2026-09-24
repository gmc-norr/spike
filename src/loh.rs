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
    /// Phase from a phased gVCF genotype: (phase set, the genotype's first
    /// haplotype carries allele2). The phase set is the PS value, or empty
    /// when the genotype is phased without PS (phased across the contig).
    phase: Option<(String, bool)>,
}

/// Identify reads from one haplotype in a genomic region.
///
/// Finds het SNPs in [region_start, region_end), phases them from reads that
/// cover two or more, picks one haplotype per phase block as the "target",
/// then classifies reads by which haplotype they carry.
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
            // gVCF gave us het SNPs; one BAM/CRAM pass collects read alleles.
            let positions: HashSet<u64> = snps.iter().map(|s| s.pos).collect();
            let alleles = collect_snp_alleles(
                alignment_path,
                chrom,
                region_start,
                region_end,
                min_mapq,
                &positions,
                ref_path,
            )?;
            (snps, Some(alleles))
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

    let read_alleles = read_alleles.unwrap_or_default();
    let target_allele = pick_target_alleles(&het_snps, &read_alleles, rng);
    let (target_set, classified_set) =
        classify_from_collected(&target_allele, &read_alleles, label);
    Ok((target_set, classified_set, target_allele))
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
    // Determine sample index from header (defaults to first sample, index 9).
    let sample_col: usize = 9;
    let mut het_snps = Vec::new();
    let mut first_other_chrom: Option<String> = None;

    if gvcf_path.ends_with(".gz") {
        // Use bcftools for indexed access to bgzipped VCF.
        // Stream stdout line-by-line via piped child process.
        let region = format!("{}:{}-{}", chrom, region_start + 1, region_end);
        let mut child = std::process::Command::new("bcftools")
            .args(["view", "-H", "-r", &region, gvcf_path])
            .stdout(std::process::Stdio::piped())
            .stderr(std::process::Stdio::piped())
            .spawn()
            .context("failed to run bcftools for gVCF reading (is bcftools in PATH?)")?;

        let stdout = child.stdout.take().unwrap();
        let reader = std::io::BufReader::new(stdout);
        use std::io::BufRead;
        for line_result in reader.lines() {
            let line = line_result.context("failed to read bcftools output")?;
            parse_gvcf_line(
                &line,
                chrom,
                region_start,
                region_end,
                sample_col,
                &mut het_snps,
                &mut first_other_chrom,
            );
        }

        let status = child.wait().context("failed to wait for bcftools")?;
        if !status.success() {
            anyhow::bail!("bcftools exited with status {}", status);
        }
    } else {
        // Plain text VCF: stream line-by-line to avoid loading entire file into memory.
        let file = std::fs::File::open(gvcf_path)
            .with_context(|| format!("failed to open gVCF: {}", gvcf_path))?;
        let reader = std::io::BufReader::new(file);
        use std::io::BufRead;
        for line_result in reader.lines() {
            let line = line_result
                .with_context(|| format!("failed to read gVCF: {}", gvcf_path))?;
            parse_gvcf_line(
                &line,
                chrom,
                region_start,
                region_end,
                sample_col,
                &mut het_snps,
                &mut first_other_chrom,
            );
        }
    }

    if het_snps.is_empty() {
        if let Some(other) = &first_other_chrom {
            log::warn!(
                "gVCF '{}': no records found for chromosome '{}', \
                 but found records for '{}'. \
                 Chromosome names may not match (e.g. 'chr1' vs '1'). \
                 Falling back to pileup-based het SNP detection.",
                gvcf_path,
                chrom,
                other,
            );
        }
    }

    het_snps.sort_by_key(|s| s.pos);
    Ok(het_snps)
}

/// Parse a single VCF/gVCF line and push any het SNP into `het_snps`.
fn parse_gvcf_line(
    line: &str,
    chrom: &str,
    region_start: u64,
    region_end: u64,
    sample_col: usize,
    het_snps: &mut Vec<HetSnp>,
    first_other_chrom: &mut Option<String>,
) {
    if line.starts_with('#') {
        return;
    }

    let fields: Vec<&str> = line.split('\t').collect();
    if fields.len() <= sample_col {
        return;
    }

    let line_chrom = fields[0];
    if line_chrom != chrom {
        // Record the first non-matching chromosome seen for a mismatch warning.
        if first_other_chrom.is_none() {
            *first_other_chrom = Some(line_chrom.to_string());
        }
        return;
    }

    // VCF POS is 1-based.
    let pos: u64 = match fields[1].parse::<u64>() {
        Ok(p) => p.saturating_sub(1), // convert to 0-based
        Err(_) => return,
    };

    if pos < region_start || pos >= region_end {
        return;
    }

    let ref_allele = fields[3].as_bytes();
    let alt_field = fields[4];

    // Handle multi-allelic: take the first ALT allele.
    let alt_allele = alt_field.split(',').next().unwrap_or(".").as_bytes();

    // Only consider SNPs (single-base REF and ALT).
    if ref_allele.len() != 1 || alt_allele.len() != 1 {
        return;
    }
    if alt_allele == b"." || alt_allele == b"*" {
        return;
    }

    // Check genotype for heterozygosity.
    let sample: Vec<&str> = fields[sample_col].split(':').collect();
    let gt_field = sample.first().copied().unwrap_or("");
    let is_het =
        gt_field == "0/1" || gt_field == "1/0" || gt_field == "0|1" || gt_field == "1|0";

    // A phased genotype also says which haplotype carries ALT; PS names the
    // phase set (absent: phased across the contig).
    let phase = match gt_field {
        "0|1" => Some(false),
        "1|0" => Some(true),
        _ => None,
    }
    .map(|first_is_alt| {
        let ps = fields[sample_col - 1]
            .split(':')
            .position(|key| key == "PS")
            .and_then(|i| sample.get(i))
            .filter(|v| **v != ".")
            .map(|v| v.to_string())
            .unwrap_or_default();
        (ps, first_is_alt)
    });

    if is_het {
        het_snps.push(HetSnp {
            pos,
            allele1: ref_allele[0].to_ascii_uppercase(),
            allele2: alt_allele[0].to_ascii_uppercase(),
            phase,
        });
    }
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
                    phase: None,
                });
            }
        }
    }

    het_snps.sort_by_key(|s| s.pos);
    (het_snps, read_alleles)
}

/// Phase het SNPs into blocks from fragments that cover two or more of them.
///
/// Returns SNP position → (block id, flip). Within a block, the haplotype
/// called "A" carries `allele1` at SNPs with flip = false and `allele2` at
/// SNPs with flip = true. SNPs no fragment links are blocks of their own.
fn phase_het_snps(
    het_snps: &[HetSnp],
    read_alleles: &HashMap<String, Vec<(u64, u8)>>,
) -> HashMap<u64, (usize, bool)> {
    let index: HashMap<u64, usize> =
        het_snps.iter().enumerate().map(|(i, s)| (s.pos, i)).collect();

    // Votes per SNP pair: +1 when a fragment carries allele1/allele2 of
    // both in the same way (cis), -1 when it carries them crossed (trans).
    let mut votes: HashMap<(usize, usize), i64> = HashMap::new();
    for alleles in read_alleles.values() {
        // (snp index, carries allele2?) at each SNP this fragment covers once.
        let mut seen: Vec<(usize, bool)> = Vec::new();
        let mut conflicted: HashSet<usize> = HashSet::new();
        for &(pos, base) in alleles {
            let Some(&i) = index.get(&pos) else { continue };
            let snp = &het_snps[i];
            let is_allele2 = if base == snp.allele1 {
                false
            } else if base == snp.allele2 {
                true
            } else {
                continue;
            };
            match seen.iter().find(|(j, _)| *j == i) {
                Some(&(_, prev)) if prev != is_allele2 => {
                    conflicted.insert(i); // mates disagree: drop this SNP
                }
                Some(_) => {}
                None => seen.push((i, is_allele2)),
            }
        }
        seen.retain(|(i, _)| !conflicted.contains(i));
        seen.sort_unstable();
        // Consecutive SNPs of the fragment are enough to link all of them.
        for pair in seen.windows(2) {
            let ((i, a), (j, b)) = (pair[0], pair[1]);
            *votes.entry((i, j)).or_insert(0) += if a == b { 1 } else { -1 };
        }
    }

    // A phased gVCF links SNPs of one phase set outright, overriding reads.
    const GVCF_LINK: i64 = 1_000_000;
    let mut last_in_set: HashMap<&str, (usize, bool)> = HashMap::new();
    for (j, snp) in het_snps.iter().enumerate() {
        let Some((set, first_is_allele2)) = &snp.phase else { continue };
        if let Some((i, prev)) = last_in_set.insert(set.as_str(), (j, *first_is_allele2)) {
            let vote = if prev == *first_is_allele2 { GVCF_LINK } else { -GVCF_LINK };
            *votes.entry((i, j)).or_insert(0) += vote;
        }
    }

    // Join SNPs along the strongest links first (union-find with parity).
    // A link that contradicts stronger ones already joined is ignored.
    let mut edges: Vec<((usize, usize), i64)> =
        votes.into_iter().filter(|(_, v)| *v != 0).collect();
    edges.sort_by_key(|&((i, j), v)| (Reverse(v.abs()), i, j));

    let n = het_snps.len();
    let mut parent: Vec<usize> = (0..n).collect();
    let mut parity: Vec<bool> = vec![false; n]; // flip relative to parent
    fn find(parent: &mut [usize], parity: &mut [bool], i: usize) -> (usize, bool) {
        if parent[i] == i {
            return (i, false);
        }
        let (root, p) = find(parent, parity, parent[i]);
        parity[i] ^= p;
        parent[i] = root;
        (root, parity[i])
    }
    for ((i, j), v) in edges {
        let (ri, pi) = find(&mut parent, &mut parity, i);
        let (rj, pj) = find(&mut parent, &mut parity, j);
        let differ = v < 0; // trans: flips must differ
        if ri != rj {
            parent[rj] = ri;
            parity[rj] = pi ^ pj ^ differ;
        }
    }

    het_snps
        .iter()
        .enumerate()
        .map(|(i, snp)| {
            let (root, flip) = find(&mut parent, &mut parity, i);
            (snp.pos, (root, flip))
        })
        .collect()
}

/// Pick the allele of one haplotype at every het SNP: phase the SNPs, then
/// flip one coin per phase block, so all SNPs of a block agree.
fn pick_target_alleles(
    het_snps: &[HetSnp],
    read_alleles: &HashMap<String, Vec<(u64, u8)>>,
    rng: &mut StdRng,
) -> HashMap<u64, u8> {
    let phase = phase_het_snps(het_snps, read_alleles);

    // One coin per block, flipped in order of the block's first SNP (SNPs
    // are sorted by position) so a seed always gives the same choice.
    let mut take_a: HashMap<usize, bool> = HashMap::new();
    let mut targets = HashMap::new();
    for snp in het_snps {
        let (block, flip) = phase[&snp.pos];
        let a = *take_a.entry(block).or_insert_with(|| rng.gen::<bool>());
        // Haplotype A carries allele2 where flip is set.
        let allele = if a != flip { snp.allele1 } else { snp.allele2 };
        targets.insert(snp.pos, allele);
    }
    let n_blocks = take_a.len();
    log::info!(
        "phased {} het SNPs into {} block(s); one haplotype picked per block",
        het_snps.len(),
        n_blocks
    );
    targets
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
        if target_count == other_count {
            // Can't tell the haplotype: leave unclassified (random suppression).
            n_ambiguous += 1;
            continue;
        }
        classified_set.insert(name.clone());
        if target_count > other_count {
            target_set.insert(name.clone());
        }
    }

    log::info!(
        "{}: classified {} fragments ({} target hap, {} other); {} ambiguous ties left unclassified",
        label,
        classified_set.len(),
        target_set.len(),
        classified_set.len() - target_set.len(),
        n_ambiguous,
    );

    (target_set, classified_set)
}

/// Collect, per fragment (mates pooled by name), the base each read shows
/// at the given SNP positions.
fn collect_snp_alleles(
    alignment_path: &str,
    chrom: &str,
    region_start: u64,
    region_end: u64,
    min_mapq: u8,
    positions: &HashSet<u64>,
    ref_path: Option<&str>,
) -> Result<HashMap<String, Vec<(u64, u8)>>> {
    if crate::extract::is_cram(alignment_path) {
        let rp =
            ref_path.ok_or_else(|| anyhow::anyhow!("CRAM input requires a reference FASTA"))?;
        collect_snp_alleles_cram(
            alignment_path,
            chrom,
            region_start,
            region_end,
            min_mapq,
            positions,
            rp,
        )
    } else {
        collect_snp_alleles_bam(
            alignment_path,
            chrom,
            region_start,
            region_end,
            min_mapq,
            positions,
        )
    }
}

/// BAM pass of [`collect_snp_alleles`].
fn collect_snp_alleles_bam(
    bam_path: &str,
    chrom: &str,
    region_start: u64,
    region_end: u64,
    min_mapq: u8,
    positions: &HashSet<u64>,
) -> Result<HashMap<String, Vec<(u64, u8)>>> {
    let mut read_alleles: HashMap<String, Vec<(u64, u8)>> = HashMap::new();

    let mut reader = noodles::bam::io::indexed_reader::Builder::default()
        .build_from_path(bam_path)
        .with_context(|| {
            format!(
                "failed to open BAM for SNP allele collection: {}",
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

        let alleles = read_alleles.entry(name).or_default();
        walk_cigar_collect(
            &seq,
            Box::new(cigar.iter()),
            align_start,
            positions,
            alleles,
        );
    }

    Ok(read_alleles)
}

/// CRAM pass of [`collect_snp_alleles`].
fn collect_snp_alleles_cram(
    cram_path: &str,
    chrom: &str,
    region_start: u64,
    region_end: u64,
    min_mapq: u8,
    positions: &HashSet<u64>,
    ref_path: &str,
) -> Result<HashMap<String, Vec<(u64, u8)>>> {
    let mut read_alleles: HashMap<String, Vec<(u64, u8)>> = HashMap::new();

    let repository = crate::extract::build_fasta_repository(ref_path)?;

    let mut reader = noodles::cram::io::indexed_reader::Builder::default()
        .set_reference_sequence_repository(repository)
        .build_from_path(cram_path)
        .with_context(|| {
            format!(
                "failed to open CRAM for SNP allele collection: {}",
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

        let alleles = read_alleles.entry(name).or_default();
        walk_cigar_collect(
            &seq,
            CigarTrait::iter(&cigar),
            align_start,
            positions,
            alleles,
        );
    }

    Ok(read_alleles)
}

/// Walk CIGAR operations and record the read's base at each SNP position.
/// Shared between BAM and CRAM passes.
fn walk_cigar_collect(
    seq: &[u8],
    cigar_ops: Box<
        dyn Iterator<Item = std::io::Result<noodles::sam::alignment::record::cigar::Op>> + '_,
    >,
    align_start: u64,
    positions: &HashSet<u64>,
    alleles: &mut Vec<(u64, u8)>,
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
                    if positions.contains(&rp) {
                        if let Some(&base) = seq.get(seq_pos + i) {
                            let base = base.to_ascii_uppercase();
                            if base != b'N' {
                                alleles.push((rp, base));
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

#[cfg(test)]
mod tests {
    use super::*;
    use rand::SeedableRng;

    fn snp(pos: u64) -> HetSnp {
        HetSnp { pos, allele1: b'A', allele2: b'G', phase: None }
    }

    fn phased(pos: u64, set: &str, first_is_allele2: bool) -> HetSnp {
        HetSnp { phase: Some((set.to_string(), first_is_allele2)), ..snp(pos) }
    }

    fn gvcf_snp(line: &str) -> Vec<HetSnp> {
        let mut snps = Vec::new();
        parse_gvcf_line(line, "chr1", 0, 1_000_000, 9, &mut snps, &mut None);
        snps
    }

    #[test]
    fn test_gvcf_phase_is_read_from_genotype_and_ps() {
        let with_ps = gvcf_snp("chr1\t101\t.\tA\tG\t50\tPASS\t.\tGT:GQ:PS\t1|0:50:77");
        assert_eq!(with_ps[0].phase, Some(("77".to_string(), true)));
        let no_ps = gvcf_snp("chr1\t101\t.\tA\tG\t50\tPASS\t.\tGT\t0|1");
        assert_eq!(no_ps[0].phase, Some((String::new(), false)));
        let unphased = gvcf_snp("chr1\t101\t.\tA\tG\t50\tPASS\t.\tGT\t0/1");
        assert_eq!(unphased[0].phase, None);
    }

    #[test]
    fn test_gvcf_phase_links_snps_no_read_covers() {
        let none: HashMap<String, Vec<(u64, u8)>> = HashMap::new();
        // Same phase set: one block; 1|0 and 0|1 get opposite flips.
        let snps = vec![phased(100, "7", true), phased(5000, "7", false), phased(9000, "7", true)];
        let phase = phase_het_snps(&snps, &none);
        assert_eq!(phase[&100].0, phase[&5000].0);
        assert_eq!(phase[&100].0, phase[&9000].0);
        assert_ne!(phase[&100].1, phase[&5000].1);
        assert_eq!(phase[&100].1, phase[&9000].1);
        // Different phase sets stay apart.
        let snps = vec![phased(100, "7", true), phased(5000, "8", true)];
        let phase = phase_het_snps(&snps, &none);
        assert_ne!(phase[&100].0, phase[&5000].0);
    }

    fn frag(alleles: &[(u64, u8)]) -> Vec<(u64, u8)> {
        alleles.to_vec()
    }

    /// Fragments linking 100-200 in trans (A with G) and 200-300 in cis.
    /// Haplotype 1: 100A 200G 300G; haplotype 2: 100G 200A 300A.
    fn linked_reads() -> HashMap<String, Vec<(u64, u8)>> {
        let mut reads = HashMap::new();
        for i in 0..5 {
            reads.insert(format!("h1a_{}", i), frag(&[(100, b'A'), (200, b'G')]));
            reads.insert(format!("h2a_{}", i), frag(&[(100, b'G'), (200, b'A')]));
            reads.insert(format!("h1b_{}", i), frag(&[(200, b'G'), (300, b'G')]));
            reads.insert(format!("h2b_{}", i), frag(&[(200, b'A'), (300, b'A')]));
        }
        reads
    }

    #[test]
    fn test_phase_links_snps_through_shared_fragments() {
        let snps = vec![snp(100), snp(200), snp(300)];
        let phase = phase_het_snps(&snps, &linked_reads());
        let (b100, f100) = phase[&100];
        let (b200, f200) = phase[&200];
        let (b300, f300) = phase[&300];
        assert!(b100 == b200 && b200 == b300, "one block: {:?}", phase);
        assert_ne!(f100, f200, "100A goes with 200G");
        assert_eq!(f200, f300, "200G goes with 300G");
    }

    #[test]
    fn test_phase_keeps_unlinked_snps_apart() {
        let snps = vec![snp(100), snp(200), snp(5000)];
        let phase = phase_het_snps(&snps, &linked_reads());
        assert_ne!(phase[&100].0, phase[&5000].0);
    }

    #[test]
    fn test_target_alleles_follow_one_haplotype_per_block() {
        // Whatever the coin flip, the targets are 100A/200G/300G or 100G/200A/300A.
        let snps = vec![snp(100), snp(200), snp(300)];
        let reads = linked_reads();
        for seed in 0..50 {
            let mut rng = StdRng::seed_from_u64(seed);
            let t = pick_target_alleles(&snps, &reads, &mut rng);
            let got = (t[&100], t[&200], t[&300]);
            assert!(
                got == (b'A', b'G', b'G') || got == (b'G', b'A', b'A'),
                "seed {}: mixed haplotype {:?}",
                seed,
                got
            );
        }
    }

    #[test]
    fn test_tied_fragment_is_unclassified() {
        // One target and one non-target allele: the fragment can't be
        // assigned, so it must fall back to random suppression.
        let target: HashMap<u64, u8> = [(100, b'A'), (200, b'G')].into();
        let reads: HashMap<String, Vec<(u64, u8)>> =
            [("tie".to_string(), frag(&[(100, b'A'), (200, b'A')]))].into();
        let (target_set, classified) = classify_from_collected(&target, &reads, "test");
        assert!(!target_set.contains("tie"));
        assert!(!classified.contains("tie"));
    }
}
