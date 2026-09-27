//! Haplotype-aware simulation: the sample has two copies of every region, and
//! an event goes on one of them.
//!
//! [`sample_copies`] finds the sample's SNPs around an event, phases the het
//! SNPs into blocks and picks, per block, the haplotype of the copy that gets
//! the event. The caller then:
//! - removes original reads by copy: the event copy's reads first, and the
//!   other copy's only above VAF 0.5, and
//! - gives synthetic reads the alleles of the copy they come from. Hom-alt
//!   SNPs are on both copies.
//!
//! So a het deletion turns het SNPs inside it homozygous (LOH), a het
//! duplication shifts them to 2:1, and SNPs in the flanks keep their balance.
//!
//! Two ways to find SNPs:
//! 1. **Pileup** (default): count the bases reads show at each position.
//! 2. **gVCF** (optional): het and hom-alt SNPs from a pre-called VCF
//!    (e.g., DeepVariant). Phased genotypes also link SNPs into blocks.

use std::cmp::Reverse;
use std::collections::{HashMap, HashSet};

use anyhow::{Context, Result};
use noodles::sam::alignment::record::cigar::op::Kind;
use noodles::sam::alignment::record::Cigar as CigarTrait;
use rand::rngs::StdRng;
use rand::Rng;

use crate::extract::record_is_on_queried_reference;

/// The sample's two copies of a region: the alleles each carries, and which
/// copy each original fragment came from.
#[derive(Debug, Default)]
pub struct SampleCopies {
    /// Alleles on the copy that gets the event: one haplotype per phase
    /// block at het SNPs, plus the hom-alt SNPs.
    pub event_copy: HashMap<u64, u8>,
    /// Alleles on the other copy: the other het alleles, plus the hom-alt SNPs.
    pub other_copy: HashMap<u64, u8>,
    /// Fragments assigned to a copy: name → true when on the event copy.
    /// Fragments that cover no het SNP, or tie, are left out.
    pub read_copy: HashMap<String, bool>,
}

/// SNPs found in a region.
#[derive(Debug, Default)]
struct RegionSnps {
    het: Vec<HetSnp>,
    /// Hom-alt SNPs: position → alt allele.
    hom_alt: HashMap<u64, u8>,
    /// Positions deleted on some haplotype: by a deletion allele it carries,
    /// or marked `*`.
    spanned: HashSet<u64>,
}

impl RegionSnps {
    /// A het SNP whose base another haplotype deletes has no copy carrying
    /// REF: both copies get ALT (a copy can't hold a deletion), so make it
    /// hom-alt.
    fn fold_spanning_deletions(&mut self) {
        let spanned = std::mem::take(&mut self.spanned);
        let (over_deletion, het): (Vec<HetSnp>, Vec<HetSnp>) = std::mem::take(&mut self.het)
            .into_iter()
            .partition(|s| spanned.contains(&s.pos));
        self.het = het;
        for s in over_deletion {
            self.hom_alt.insert(s.pos, s.allele2);
        }
    }
}

/// Call SNPs from pileup allele counts (A, C, G, T per position).
///
/// Het: the top two alleles each make up 20–80% of the reads. Hom-alt: one
/// allele makes up at least 90% and differs from the reference.
/// `ref_seq` holds the reference from `region_start`. Positions need 10 reads.
fn call_snps(counts: &HashMap<u64, [u32; 4]>, region_start: u64, ref_seq: &[u8]) -> RegionSnps {
    const BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];
    let mut snps = RegionSnps::default();

    for (&pos, c) in counts {
        let total: u32 = c.iter().sum();
        if total < 10 {
            continue;
        }
        let mut sorted: Vec<(u8, u32)> = BASES.iter().copied().zip(c.iter().copied()).collect();
        sorted.sort_by_key(|&(_, n)| Reverse(n));
        let f1 = sorted[0].1 as f64 / total as f64;
        let f2 = sorted[1].1 as f64 / total as f64;

        if (0.2..=0.8).contains(&f1) && (0.2..=0.8).contains(&f2) {
            snps.het.push(HetSnp {
                pos,
                allele1: sorted[0].0,
                allele2: sorted[1].0,
                phase: None,
            });
        } else if f1 >= 0.9 {
            let ref_base = pos
                .checked_sub(region_start)
                .and_then(|i| ref_seq.get(i as usize))
                .map(|b| b.to_ascii_uppercase());
            if matches!(ref_base, Some(r) if BASES.contains(&r) && r != sorted[0].0) {
                snps.hom_alt.insert(pos, sorted[0].0);
            }
        }
    }

    snps.het.sort_by_key(|s| s.pos);
    snps
}

/// Pick the event copy's haplotype (one coin per phase block), give the other
/// copy the other het alleles, add hom-alts to both, and assign fragments.
fn copies_from_snps(
    het_snps: &[HetSnp],
    hom_alt: &HashMap<u64, u8>,
    read_alleles: &HashMap<String, Vec<(u64, u8)>>,
    rng: &mut StdRng,
) -> SampleCopies {
    let mut event_copy = if het_snps.is_empty() {
        HashMap::new()
    } else {
        pick_target_alleles(het_snps, read_alleles, rng)
    };
    let read_copy = classify_from_collected(&event_copy, read_alleles);
    let mut other_copy: HashMap<u64, u8> = het_snps
        .iter()
        .map(|s| {
            let other = if event_copy[&s.pos] == s.allele1 { s.allele2 } else { s.allele1 };
            (s.pos, other)
        })
        .collect();
    for (&pos, &alt) in hom_alt {
        event_copy.insert(pos, alt);
        other_copy.insert(pos, alt);
    }
    SampleCopies {
        event_copy,
        other_copy,
        read_copy,
    }
}

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

/// Find the sample's SNPs in [region_start, region_end), pick the event
/// copy's haplotype (one coin per phase block) and assign fragments to copies.
///
/// `ref_seq` is the reference over the region; pileup needs it to tell
/// hom-alt from hom-ref. SNPs come from the gVCF when it has het SNPs here,
/// otherwise from pileup. With no het SNPs, no fragment is assigned and the
/// caller falls back to random suppression.
#[allow(clippy::too_many_arguments)]
pub fn sample_copies(
    alignment_path: &str,
    chrom: &str,
    region_start: u64,
    region_end: u64,
    ref_seq: &[u8],
    min_mapq: u8,
    gvcf_path: Option<&str>,
    ref_path: Option<&str>,
    rng: &mut StdRng,
) -> Result<SampleCopies> {
    let from_gvcf = match gvcf_path {
        Some(gvcf) => {
            let snps = load_snps_from_gvcf(gvcf, chrom, region_start, region_end)
                .with_context(|| GvcfUnreadable {
                    path: gvcf.to_string(),
                })?;
            if snps.het.is_empty() {
                log::info!("no het SNPs from gVCF, trying pileup fallback");
                None
            } else {
                Some(snps)
            }
        }
        None => None,
    };
    let snps = match from_gvcf {
        Some(snps) => snps,
        None => {
            let counts = count_alleles(
                alignment_path,
                chrom,
                region_start,
                region_end,
                min_mapq,
                ref_path,
            )?;
            call_snps(&counts, region_start, ref_seq)
        }
    };
    log::info!(
        "{}:{}-{}: {} het SNPs, {} hom-alt SNPs",
        chrom,
        region_start,
        region_end,
        snps.het.len(),
        snps.hom_alt.len(),
    );

    // A second pass reads each fragment's bases at the het SNPs only.
    let read_alleles = if snps.het.is_empty() {
        HashMap::new()
    } else {
        let positions: HashSet<u64> = snps.het.iter().map(|s| s.pos).collect();
        collect_snp_alleles(
            alignment_path,
            chrom,
            region_start,
            region_end,
            min_mapq,
            &positions,
            ref_path,
        )?
    };
    Ok(copies_from_snps(&snps.het, &snps.hom_alt, &read_alleles, rng))
}

/// The `--gvcf` could not be read. Unlike a region with no SNPs, this stops
/// the run: going on would simulate the region without the sample's SNPs,
/// and the output could not be told from a sample that has none (CR3).
/// The caller finds it with `downcast_ref` to tell it from a pileup failure.
#[derive(Debug)]
pub struct GvcfUnreadable {
    pub path: String,
}

impl std::fmt::Display for GvcfUnreadable {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "could not read the --gvcf '{}'; spike stops rather than simulate without \
             the sample's SNPs",
            self.path
        )
    }
}

/// What a gVCF's chromosome names say about the one being asked for.
#[derive(Debug, PartialEq)]
enum ContigMatch {
    /// It names this chromosome, or names nothing that contradicts it.
    Match,
    /// It names the other convention instead: 'chr20' vs '20'.
    Renamed(String),
    /// Its header names contigs, but no spelling of this one.
    Absent,
}

/// What spike does with a region the gVCF gave no SNPs for.
#[derive(Clone, Copy)]
enum NextStep {
    /// The gVCF was read and had none here: the pileup calls them instead.
    Pileup,
    /// The gVCF could not be read, so the run stops (CR3). It used to go on
    /// with no copies for the region, suppressing its reads at random.
    Stop,
}

impl NextStep {
    fn describe(self) -> &'static str {
        match self {
            NextStep::Pileup => "Falling back to pileup-based het SNP detection.",
            NextStep::Stop => "The run stops here.",
        }
    }
}

/// The `##contig=<ID=...>` names in a VCF header, plain or bgzipped.
///
/// Only the header is decompressed: the loop stops at the first record.
/// Best-effort: a file that cannot be read here reads as no contigs, and the
/// read below reports the failure.
fn header_contigs(gvcf_path: &str) -> Vec<String> {
    use std::io::BufRead;
    let Ok(file) = std::fs::File::open(gvcf_path) else {
        return Vec::new();
    };
    let reader: Box<dyn BufRead> = if gvcf_path.ends_with(".gz") {
        Box::new(std::io::BufReader::new(flate2::read::MultiGzDecoder::new(file)))
    } else {
        Box::new(std::io::BufReader::new(file))
    };

    let mut contigs = Vec::new();
    for line in reader.lines() {
        let Ok(line) = line else { break };
        if !line.starts_with('#') {
            break;
        }
        if let Some(fields) = line.strip_prefix("##contig=<") {
            if let Some(id) = fields.split(',').next().and_then(|f| f.strip_prefix("ID=")) {
                contigs.push(id.trim_end_matches('>').to_string());
            }
        }
    }
    contigs
}

/// The same chromosome under the other naming convention: 'chr20' <-> '20'.
fn other_spelling(chrom: &str) -> String {
    match chrom.strip_prefix("chr") {
        Some(bare) => bare.to_string(),
        None => format!("chr{}", chrom),
    }
}

/// Compare `chrom` against the names the gVCF uses. `seen_other` is a
/// chromosome seen in its records, for a header that lists no contigs.
fn match_contig(contigs: &[String], seen_other: Option<&str>, chrom: &str) -> ContigMatch {
    let other = other_spelling(chrom);
    if contigs.is_empty() {
        // Nothing to compare against. A record on another chromosome only
        // tells us something when it is this chromosome, spelled the other way.
        return match seen_other {
            Some(seen) if seen == other => ContigMatch::Renamed(seen.to_string()),
            _ => ContigMatch::Match,
        };
    }
    if contigs.iter().any(|c| c == chrom) {
        ContigMatch::Match
    } else if contigs.contains(&other) {
        ContigMatch::Renamed(other)
    } else {
        ContigMatch::Absent
    }
}

/// The warning for a gVCF that cannot hold SNPs for `chrom`, or None when
/// its names are fine. `next` says what spike does with the region instead.
fn contig_warning(
    matched: &ContigMatch,
    gvcf_path: &str,
    chrom: &str,
    next: NextStep,
) -> Option<String> {
    match matched {
        ContigMatch::Match => None,
        ContigMatch::Renamed(other) => Some(format!(
            "gVCF '{}' names chromosome '{}', not '{}': its chromosome naming does not \
             match the one asked for (e.g. 'chr1' vs '1'), so it has no SNPs here. {}",
            gvcf_path,
            other,
            chrom,
            next.describe(),
        )),
        ContigMatch::Absent => Some(format!(
            "gVCF '{}' has no chromosome '{}', so it has no SNPs here. {}",
            gvcf_path,
            chrom,
            next.describe(),
        )),
    }
}

/// Load het and hom-alt SNPs from a VCF/gVCF file.
///
/// For `.vcf.gz` files, uses `bcftools view -H -r region` for efficient
/// indexed access. For plain `.vcf` files, reads and filters line by line.
fn load_snps_from_gvcf(
    gvcf_path: &str,
    chrom: &str,
    region_start: u64,
    region_end: u64,
) -> Result<RegionSnps> {
    // Determine sample index from header (defaults to first sample, index 9).
    let sample_col: usize = 9;
    let mut snps = RegionSnps::default();
    let mut first_other_chrom: Option<String> = None;
    // The gVCF's own names tell whether it can hold SNPs here at all. An
    // empty result cannot: `bcftools view -r chr20:...` on a file naming
    // that chromosome '20' prints nothing and exits 0, exactly like a
    // region that genuinely has no SNPs.
    let contigs = header_contigs(gvcf_path);

    if gvcf_path.ends_with(".gz") {
        // Use bcftools for indexed access to bgzipped VCF.
        // Stream stdout line-by-line via piped child process.
        let region = format!("{}:{}-{}", chrom, region_start + 1, region_end);
        let mut child = std::process::Command::new("bcftools")
            .args(["view", "-H", "-r", &region, gvcf_path])
            .stdout(std::process::Stdio::piped())
            .stderr(std::process::Stdio::piped())
            .spawn()
            .with_context(|| {
                format!(
                    "failed to run bcftools for gVCF reading (is bcftools in PATH?). {}",
                    NextStep::Stop.describe()
                )
            })?;

        let stdout = child.stdout.take().unwrap();
        let mut stderr = child.stderr.take().unwrap();
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
                &mut snps,
                &mut first_other_chrom,
            );
        }

        // bcftools says why it failed (no index, not bgzipped): pass it on.
        let mut bcftools_error = String::new();
        use std::io::Read;
        let _ = stderr.read_to_string(&mut bcftools_error);

        let status = child.wait().context("failed to wait for bcftools")?;
        if !status.success() {
            // A failed read is not an empty region: the names may be wrong
            // as well, and neither is repaired by the pileup pass below.
            if let Some(warning) = contig_warning(
                &match_contig(&contigs, None, chrom),
                gvcf_path,
                chrom,
                NextStep::Stop,
            ) {
                log::warn!("{}", warning);
            }
            anyhow::bail!(
                "bcftools exited with status {} on gVCF '{}': {}. {}",
                status,
                gvcf_path,
                bcftools_error.trim(),
                NextStep::Stop.describe(),
            );
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
                &mut snps,
                &mut first_other_chrom,
            );
        }
    }

    // Warn only about a chromosome the gVCF cannot name, never about a
    // region that simply has no SNPs in it.
    if snps.het.is_empty() && snps.hom_alt.is_empty() {
        let matched = match_contig(&contigs, first_other_chrom.as_deref(), chrom);
        if let Some(warning) = contig_warning(&matched, gvcf_path, chrom, NextStep::Pileup) {
            log::warn!("{}", warning);
        }
    }

    snps.fold_spanning_deletions();
    snps.het.sort_by_key(|s| s.pos);
    Ok(snps)
}

/// Parse a single VCF/gVCF line and push any het SNP into `het_snps`.
fn parse_gvcf_line(
    line: &str,
    chrom: &str,
    region_start: u64,
    region_end: u64,
    sample_col: usize,
    snps: &mut RegionSnps,
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
    let sample: Vec<&str> = fields[sample_col].split(':').collect();
    let gt_field = sample.first().copied().unwrap_or("");

    // Bases deleted on a haplotype: by a deletion allele it carries, or
    // marked `*` (deleted by an upstream deletion).
    for (i, alt) in alt_field.split(',').enumerate() {
        let carried = gt_field.split(['/', '|']).any(|a| a == (i + 1).to_string());
        let kept = if alt == "*" { 0 } else { alt.len() };
        if carried && kept < ref_allele.len() {
            snps.spanned.extend(pos + kept as u64..pos + ref_allele.len() as u64);
        }
    }

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

    if gt_field == "1/1" || gt_field == "1|1" {
        snps.hom_alt.insert(pos, alt_allele[0].to_ascii_uppercase());
    }

    if is_het {
        snps.het.push(HetSnp {
            pos,
            allele1: ref_allele[0].to_ascii_uppercase(),
            allele2: alt_allele[0].to_ascii_uppercase(),
            phase,
        });
    }
}

/// Count the bases (A, C, G, T) reads show at each position of
/// [region_start, region_end).
fn count_alleles(
    alignment_path: &str,
    chrom: &str,
    region_start: u64,
    region_end: u64,
    min_mapq: u8,
    ref_path: Option<&str>,
) -> Result<HashMap<u64, [u32; 4]>> {
    if crate::extract::is_cram(alignment_path) {
        let rp =
            ref_path.ok_or_else(|| anyhow::anyhow!("CRAM input requires a reference FASTA"))?;
        count_alleles_cram(alignment_path, chrom, region_start, region_end, min_mapq, rp)
    } else {
        count_alleles_bam(alignment_path, chrom, region_start, region_end, min_mapq)
    }
}

/// BAM pass of [`count_alleles`], on the thread pool (see
/// [`crate::extract::fold_bam_region`]).
fn count_alleles_bam(
    bam_path: &str,
    chrom: &str,
    region_start: u64,
    region_end: u64,
    min_mapq: u8,
) -> Result<HashMap<u64, [u32; 4]>> {
    let n_chunks = crate::extract::read_chunks(region_start, region_end);
    count_alleles_bam_in(bam_path, chrom, region_start, region_end, min_mapq, n_chunks)
}

/// [`count_alleles_bam`] in `n_chunks` chunks. Counts add, so the chunks'
/// counts summed are the counts of one pass.
fn count_alleles_bam_in(
    bam_path: &str,
    chrom: &str,
    region_start: u64,
    region_end: u64,
    min_mapq: u8,
    n_chunks: usize,
) -> Result<HashMap<u64, [u32; 4]>> {
    let per_chunk = crate::extract::fold_bam_region(
        bam_path,
        chrom,
        region_start,
        region_end,
        n_chunks,
        "failed to open BAM for pileup:",
        HashMap::new,
        |allele_counts: &mut HashMap<u64, [u32; 4]>, record| {
            let flags = record.flags();

            if flags.is_unmapped()
                || flags.is_secondary()
                || flags.is_supplementary()
                || flags.is_duplicate()
                || flags.is_qc_fail()
            {
                return Ok(());
            }

            let mq: u8 = match record.mapping_quality() {
                Some(q) => u8::from(q),
                None => 0,
            };
            if mq < min_mapq {
                return Ok(());
            }

            let align_start = match record.alignment_start() {
                Some(Ok(p)) => usize::from(p).saturating_sub(1) as u64,
                _ => return Ok(()),
            };

            let seq: Vec<u8> = record.sequence().iter().collect();
            let cigar = record.cigar();

            walk_cigar_count(
                &seq,
                Box::new(cigar.iter()),
                align_start,
                region_start,
                region_end,
                allele_counts,
            );
            Ok(())
        },
    )?;

    let mut allele_counts: HashMap<u64, [u32; 4]> = HashMap::new();
    for chunk in per_chunk {
        for (pos, counts) in chunk {
            let total = allele_counts.entry(pos).or_insert([0; 4]);
            for (t, c) in total.iter_mut().zip(counts) {
                *t += c;
            }
        }
    }
    Ok(allele_counts)
}

/// CRAM pass of [`count_alleles`].
fn count_alleles_cram(
    cram_path: &str,
    chrom: &str,
    region_start: u64,
    region_end: u64,
    min_mapq: u8,
    ref_path: &str,
) -> Result<HashMap<u64, [u32; 4]>> {
    let mut allele_counts: HashMap<u64, [u32; 4]> = HashMap::new();

    let repository = crate::extract::build_fasta_repository(ref_path)?;

    let start_pos = crate::extract::safe_noodles_position(region_start + 1);
    let end_pos = crate::extract::safe_noodles_position(region_end);
    let region = noodles::core::Region::new(chrom, start_pos..=end_pos);
    let (mut reader, header) =
        crate::extract::open_cram_reader_for_region(cram_path, &repository, &region)
            .with_context(|| format!("failed to open CRAM for pileup: {}", cram_path))?;
    let query = reader.query(&header, &region)?;
    // `query` has already rejected an unknown contig, so this is `Some`.
    let queried_reference_sequence_id = header.reference_sequences().get_index_of(chrom.as_bytes());

    for rec_result in query {
        let cram_record = rec_result?;
        let buf = cram_record
            .try_into_alignment_record(&header)
            .with_context(|| "failed to convert CRAM record")?;

        // A container holding several contigs is decoded whole and `Query`
        // filters on coordinates alone, so another contig's bases would land
        // in this pileup and change which SNPs are called het (L2, N4).
        if !record_is_on_queried_reference(&buf, queried_reference_sequence_id) {
            continue;
        }

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

        let align_start = match buf.alignment_start() {
            Some(p) => usize::from(p).saturating_sub(1) as u64,
            None => continue,
        };

        let seq: Vec<u8> = buf.sequence().as_ref().to_vec();
        let cigar = buf.cigar();

        walk_cigar_count(
            &seq,
            CigarTrait::iter(&cigar),
            align_start,
            region_start,
            region_end,
            &mut allele_counts,
        );
    }

    Ok(allele_counts)
}

/// Walk CIGAR operations and count the bases in the region.
/// Shared between BAM and CRAM passes.
fn walk_cigar_count(
    seq: &[u8],
    cigar_ops: Box<
        dyn Iterator<Item = std::io::Result<noodles::sam::alignment::record::cigar::Op>> + '_,
    >,
    align_start: u64,
    region_start: u64,
    region_end: u64,
    allele_counts: &mut HashMap<u64, [u32; 4]>,
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
                        let idx = match seq.get(seq_pos + i).map(|b| b.to_ascii_uppercase()) {
                            Some(b'A') => 0,
                            Some(b'C') => 1,
                            Some(b'G') => 2,
                            Some(b'T') => 3,
                            _ => continue,
                        };
                        allele_counts.entry(rp).or_insert([0; 4])[idx] += 1;
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

/// Assign fragments to a copy by the het alleles they carry.
///
/// Returns fragment name → true when it carries more target (event copy)
/// alleles than other alleles. Fragments that cover no het SNP, or tie, are
/// left out: they fall back to random suppression/copying at the VAF rate.
fn classify_from_collected(
    target_allele: &HashMap<u64, u8>,
    read_alleles: &HashMap<String, Vec<(u64, u8)>>,
) -> HashMap<String, bool> {
    let mut read_copy = HashMap::new();
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
        read_copy.insert(name.clone(), target_count > other_count);
    }

    let n_target = read_copy.values().filter(|&&e| e).count();
    log::info!(
        "classified {} fragments ({} event copy, {} other); {} ambiguous ties left unclassified",
        read_copy.len(),
        n_target,
        read_copy.len() - n_target,
        n_ambiguous,
    );

    read_copy
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

/// BAM pass of [`collect_snp_alleles`], on the thread pool (see
/// [`crate::extract::fold_bam_region`]).
fn collect_snp_alleles_bam(
    bam_path: &str,
    chrom: &str,
    region_start: u64,
    region_end: u64,
    min_mapq: u8,
    positions: &HashSet<u64>,
) -> Result<HashMap<String, Vec<(u64, u8)>>> {
    let n_chunks = crate::extract::read_chunks(region_start, region_end);
    collect_snp_alleles_bam_in(bam_path, chrom, region_start, region_end, min_mapq, positions, n_chunks)
}

/// [`collect_snp_alleles_bam`] in `n_chunks` chunks. Each chunk holds its
/// records in file order and the chunks are joined in order, so every read
/// name's alleles come in the order one pass gives them.
fn collect_snp_alleles_bam_in(
    bam_path: &str,
    chrom: &str,
    region_start: u64,
    region_end: u64,
    min_mapq: u8,
    positions: &HashSet<u64>,
    n_chunks: usize,
) -> Result<HashMap<String, Vec<(u64, u8)>>> {
    let per_chunk = crate::extract::fold_bam_region(
        bam_path,
        chrom,
        region_start,
        region_end,
        n_chunks,
        "failed to open BAM for SNP allele collection:",
        HashMap::new,
        |read_alleles: &mut HashMap<String, Vec<(u64, u8)>>, record| {
            let flags = record.flags();

            if flags.is_unmapped()
                || flags.is_secondary()
                || flags.is_supplementary()
                || flags.is_duplicate()
                || flags.is_qc_fail()
            {
                return Ok(());
            }

            let mq: u8 = match record.mapping_quality() {
                Some(q) => u8::from(q),
                None => 0,
            };
            if mq < min_mapq {
                return Ok(());
            }

            let name = match record.name() {
                Some(n) => String::from_utf8_lossy(n.as_ref()).into_owned(),
                None => return Ok(()),
            };

            let align_start = match record.alignment_start() {
                Some(Ok(p)) => usize::from(p).saturating_sub(1) as u64,
                _ => return Ok(()),
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
            Ok(())
        },
    )?;

    let mut read_alleles: HashMap<String, Vec<(u64, u8)>> = HashMap::new();
    for chunk in per_chunk {
        for (name, alleles) in chunk {
            read_alleles.entry(name).or_default().extend(alleles);
        }
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

    let start_pos = crate::extract::safe_noodles_position(region_start + 1);
    let end_pos = crate::extract::safe_noodles_position(region_end);
    let region = noodles::core::Region::new(chrom, start_pos..=end_pos);
    let (mut reader, header) =
        crate::extract::open_cram_reader_for_region(cram_path, &repository, &region)
            .with_context(|| {
                format!(
                    "failed to open CRAM for SNP allele collection: {}",
                    cram_path
                )
            })?;
    let query = reader.query(&header, &region)?;
    // `query` has already rejected an unknown contig, so this is `Some`.
    let queried_reference_sequence_id = header.reference_sequences().get_index_of(chrom.as_bytes());

    for rec_result in query {
        let cram_record = rec_result?;
        let buf = cram_record
            .try_into_alignment_record(&header)
            .with_context(|| "failed to convert CRAM record")?;

        // A container holding several contigs is decoded whole and `Query`
        // filters on coordinates alone, so another contig's bases would land
        // in this pileup and change which SNPs are called het (L2, N4).
        if !record_is_on_queried_reference(&buf, queried_reference_sequence_id) {
            continue;
        }

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
pub(crate) mod tests {
    use super::*;

    /// chrT, 20 kb: 600 fragments whose mates overlap by 50 bases (read 1 at
    /// `p`, read 2 at `p + 50`, 100M each), bases varying by position and
    /// fragment, every seventh fragment at MAPQ 10.
    fn overlapping_mates_bam(dir: &std::path::Path) -> String {
        use noodles::sam::alignment::record::cigar::{op::Kind, Op};
        use noodles::sam::alignment::record::{Flags, MappingQuality};
        use noodles::sam::alignment::record_buf::Sequence;
        use noodles::sam::alignment::RecordBuf;
        let mut records: Vec<(usize, RecordBuf)> = Vec::new();
        for i in 0..600usize {
            let p = 2_000 + i * 17;
            let mapq = if i % 7 == 0 { 10 } else { 60 };
            for (flags, pos) in [(0x63u16, p), (0x93u16, p + 50)] {
                let seq: Vec<u8> = (0..100).map(|k| b"ACGT"[(pos + k + i % 3) % 4]).collect();
                let record = RecordBuf::builder()
                    .set_name(format!("f{}", i))
                    .set_flags(Flags::from(flags))
                    .set_reference_sequence_id(0)
                    .set_alignment_start(noodles::core::Position::new(pos).unwrap())
                    .set_mapping_quality(MappingQuality::new(mapq).unwrap())
                    .set_cigar([Op::new(Kind::Match, 100)].into_iter().collect())
                    .set_sequence(Sequence::from(seq))
                    .build();
                records.push((pos, record));
            }
        }
        records.sort_by_key(|(pos, _)| *pos);
        let records: Vec<RecordBuf> = records.into_iter().map(|(_, r)| r).collect();
        crate::extract::test_fixtures::write_one_contig_bam(&dir.join("mates.bam"), "chrT", 20_000, &records)
    }

    #[test]
    fn test_the_pileup_passes_read_in_chunks_give_what_one_query_gives() {
        // Both pileup passes may read their region chunk by chunk on the
        // thread pool. The counts must be the same, and each read name's
        // alleles must come in the same order -- mates overlap, so a name's
        // list interleaves its two records.
        let dir = std::env::temp_dir().join(format!("spike_pileup_chunks_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let bam = overlapping_mates_bam(&dir);
        let (start, end) = (3_000, 11_000);
        let positions: HashSet<u64> = (start..end).step_by(7).collect();
        let pool4 = rayon::ThreadPoolBuilder::new().num_threads(4).build().unwrap();

        let counts = count_alleles_bam_in(&bam, "chrT", start, end, 20, 1).unwrap();
        assert!(counts.len() > 7_000, "{} positions counted", counts.len());
        assert_eq!(pool4.install(|| count_alleles_bam_in(&bam, "chrT", start, end, 20, 6)).unwrap(), counts);

        let alleles = collect_snp_alleles_bam_in(&bam, "chrT", start, end, 20, &positions, 1).unwrap();
        assert!(alleles.values().any(|a| a.windows(2).any(|w| w[1].0 < w[0].0)), "no name interleaves its mates");
        assert_eq!(
            pool4.install(|| collect_snp_alleles_bam_in(&bam, "chrT", start, end, 20, &positions, 6)).unwrap(),
            alleles
        );
        let _ = std::fs::remove_dir_all(&dir);
    }
    use rand::SeedableRng;

    fn snp(pos: u64) -> HetSnp {
        HetSnp { pos, allele1: b'A', allele2: b'G', phase: None }
    }

    fn phased(pos: u64, set: &str, first_is_allele2: bool) -> HetSnp {
        HetSnp { phase: Some((set.to_string(), first_is_allele2)), ..snp(pos) }
    }

    fn gvcf_snps(line: &str) -> RegionSnps {
        let mut snps = RegionSnps::default();
        parse_gvcf_line(line, "chr1", 0, 1_000_000, 9, &mut snps, &mut None);
        snps
    }

    #[test]
    fn test_gvcf_phase_is_read_from_genotype_and_ps() {
        let with_ps = gvcf_snps("chr1\t101\t.\tA\tG\t50\tPASS\t.\tGT:GQ:PS\t1|0:50:77");
        assert_eq!(with_ps.het[0].phase, Some(("77".to_string(), true)));
        let no_ps = gvcf_snps("chr1\t101\t.\tA\tG\t50\tPASS\t.\tGT\t0|1");
        assert_eq!(no_ps.het[0].phase, Some((String::new(), false)));
        let unphased = gvcf_snps("chr1\t101\t.\tA\tG\t50\tPASS\t.\tGT\t0/1");
        assert_eq!(unphased.het[0].phase, None);
    }

    #[test]
    fn test_gvcf_hom_alt_snps_are_kept() {
        let unphased = gvcf_snps("chr1\t101\t.\tA\tG\t50\tPASS\t.\tGT\t1/1");
        assert!(unphased.het.is_empty());
        assert_eq!(unphased.hom_alt, [(100, b'G')].into());
        let phased = gvcf_snps("chr1\t101\t.\tA\tG\t50\tPASS\t.\tGT\t1|1");
        assert_eq!(phased.hom_alt, [(100, b'G')].into());
    }

    #[test]
    fn test_gvcf_het_snp_over_a_spanning_deletion_counts_as_hom_alt() {
        // One haplotype has the ALT base, the other has that base deleted.
        // No copy carries REF, so reads must not get it: the SNP is treated
        // as hom-alt. The deletion shows as ALT * (101) or as the deletion
        // record itself (300-301).
        let mut snps = RegionSnps::default();
        for line in [
            "chr1\t101\t.\tG\t*\t.\tPASS\t.\tGT\t1|0",
            "chr1\t101\t.\tG\tC\t.\tPASS\t.\tGT\t0|1",
            "chr1\t201\t.\tA\tT\t.\tPASS\t.\tGT\t0|1",
            "chr1\t300\t.\tCG\tC\t.\tPASS\t.\tGT\t1|0",
            "chr1\t301\t.\tG\tT\t.\tPASS\t.\tGT\t0|1",
        ] {
            parse_gvcf_line(line, "chr1", 0, 1_000_000, 9, &mut snps, &mut None);
        }
        snps.fold_spanning_deletions();
        assert_eq!(snps.het.iter().map(|s| s.pos).collect::<Vec<_>>(), vec![200]);
        assert_eq!(snps.hom_alt, [(100, b'C'), (300, b'T')].into());
    }

    /// Collects `log::warn!` messages so a test can assert which branch a
    /// call took. The logger is global, so tests filter by their own path.
    pub(crate) mod capture {
        use std::sync::{Mutex, OnceLock};

        static LINES: Mutex<Vec<String>> = Mutex::new(Vec::new());
        static INSTALLED: OnceLock<()> = OnceLock::new();

        struct Collector;

        impl log::Log for Collector {
            fn enabled(&self, metadata: &log::Metadata) -> bool {
                metadata.level() <= log::Level::Warn
            }
            fn log(&self, record: &log::Record) {
                if self.enabled(record.metadata()) {
                    LINES.lock().unwrap().push(record.args().to_string());
                }
            }
            fn flush(&self) {}
        }

        /// Install the collector once; every test may call this.
        pub fn install() {
            INSTALLED.get_or_init(|| {
                let _ = log::set_boxed_logger(Box::new(Collector));
                log::set_max_level(log::LevelFilter::Warn);
            });
        }

        /// The warnings logged so far that mention `needle`.
        pub fn warnings_matching(needle: &str) -> Vec<String> {
            LINES.lock().unwrap().iter().filter(|l| l.contains(needle)).cloned().collect()
        }
    }

    fn test_dir(name: &str) -> std::path::PathBuf {
        let dir = std::env::temp_dir().join(format!("spike_test_{}_{}", name, std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        dir
    }

    /// Gzipped, not bgzipped, and unindexed: bcftools cannot query it.
    fn write_gzip(path: &std::path::Path, text: &str) {
        let mut encoder = flate2::write::GzEncoder::new(
            std::fs::File::create(path).unwrap(),
            flate2::Compression::default(),
        );
        std::io::Write::write_all(&mut encoder, text.as_bytes()).unwrap();
        encoder.finish().unwrap();
    }

    /// Fail, by name, when the test's external dependency is missing.
    ///
    /// "Resolvable" here means what it means to the production path: that
    /// `Command::new("bcftools")` can spawn. `load_snps_from_gvcf` hands the
    /// bare name to the OS, which searches PATH itself, so spawning it is the
    /// only probe that cannot disagree with the real call -- a hand-rolled
    /// PATH walk can (a non-executable file of that name, a directory, a
    /// dangling symlink). The exit status says nothing about presence and is
    /// ignored; production likewise treats a spawn that succeeded as "present"
    /// and reports a bad status as its own, separate error.
    ///
    /// This diagnoses, it does not tolerate: a missing test dependency is a
    /// broken environment, so the test still fails. Only the message changes,
    /// so the reader installs bcftools instead of hunting a logic bug.
    #[track_caller]
    fn require_bcftools() {
        let spawned = std::process::Command::new("bcftools")
            .arg("--version")
            .stdout(std::process::Stdio::null())
            .stderr(std::process::Stdio::null())
            .status()
            .is_ok();
        assert!(
            spawned,
            "this test requires bcftools on PATH: it reads a .vcf.gz, which \
             load_snps_from_gvcf queries with `bcftools view`. Without bcftools \
             the read fails at the spawn and never reaches the behaviour under \
             test. Install bcftools and re-run."
        );
    }

    #[test]
    fn test_a_gvcf_naming_the_chromosome_warns_about_nothing() {
        // A region with no SNPs in it is normal: only the names decide.
        let contigs = ["chr20".to_string(), "chr21".to_string()];
        assert_eq!(match_contig(&contigs, Some("chr21"), "chr20"), ContigMatch::Match);
        assert_eq!(contig_warning(&ContigMatch::Match, "g.vcf", "chr20", NextStep::Pileup), None);
    }

    #[test]
    fn test_the_other_naming_convention_is_read_from_the_header() {
        let contigs = ["19".to_string(), "20".to_string()];
        assert_eq!(match_contig(&contigs, None, "chr20"), ContigMatch::Renamed("20".to_string()));
        let warning =
            contig_warning(&ContigMatch::Renamed("20".into()), "g.vcf.gz", "chr20", NextStep::Pileup)
                .unwrap();
        assert!(warning.contains("names chromosome '20', not 'chr20'"), "{}", warning);
        assert!(warning.contains("Falling back to pileup"), "{}", warning);
    }

    #[test]
    fn test_a_chromosome_the_gvcf_lacks_is_told_apart_from_a_renamed_one() {
        let contigs = ["chr20".to_string()];
        assert_eq!(match_contig(&contigs, None, "chrM"), ContigMatch::Absent);
        let warning =
            contig_warning(&ContigMatch::Absent, "g.vcf", "chrM", NextStep::Pileup).unwrap();
        assert!(warning.contains("has no chromosome 'chrM'"), "{}", warning);
    }

    #[test]
    fn test_a_header_without_contigs_warns_only_for_the_other_spelling() {
        // Nothing to compare against: a record on another chromosome says
        // nothing unless it is this chromosome under the other spelling.
        assert_eq!(match_contig(&[], Some("chr21"), "chr20"), ContigMatch::Match);
        assert_eq!(match_contig(&[], Some("20"), "chr20"), ContigMatch::Renamed("20".to_string()));
    }

    #[test]
    fn test_a_failed_gvcf_read_says_the_run_stops_not_pileup() {
        let warning = contig_warning(
            &ContigMatch::Renamed("20".into()),
            "g.vcf.gz",
            "chr20",
            NextStep::Stop,
        )
        .unwrap();
        assert!(warning.contains("The run stops here."), "{}", warning);
        assert!(!warning.contains("pileup"), "{}", warning);
    }

    #[test]
    fn test_contigs_are_read_from_plain_and_gzipped_headers() {
        // The .vcf.gz path must not need a bcftools query to learn the
        // names: a region query on a renamed file returns nothing, exit 0.
        let dir = test_dir("loh_contigs");
        let header = "##fileformat=VCFv4.2\n\
                      ##contig=<ID=20,length=64444167>\n\
                      ##contig=<ID=21,length=46709983>\n\
                      #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002\n\
                      20\t101\t.\tA\tG\t50\tPASS\t.\tGT\t0|1\n";

        let plain = dir.join("plain.vcf");
        std::fs::write(&plain, header).unwrap();
        assert_eq!(header_contigs(plain.to_str().unwrap()), ["20", "21"]);

        let gzipped = dir.join("gzipped.vcf.gz");
        write_gzip(&gzipped, header);
        assert_eq!(header_contigs(gzipped.to_str().unwrap()), ["20", "21"]);
    }

    #[test]
    fn test_a_renamed_plain_gvcf_warns() {
        capture::install();
        let dir = test_dir("loh_warn_renamed");

        // '20' naming: this file can hold no 'chr20' SNP at all.
        let renamed = dir.join("renamed.vcf");
        std::fs::write(
            &renamed,
            "##fileformat=VCFv4.2\n\
             ##contig=<ID=20,length=64444167>\n\
             #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002\n\
             20\t101\t.\tA\tG\t50\tPASS\t.\tGT\t0|1\n",
        )
        .unwrap();
        let snps = load_snps_from_gvcf(renamed.to_str().unwrap(), "chr20", 0, 1000).unwrap();
        assert!(snps.het.is_empty());
        let warnings = capture::warnings_matching(renamed.to_str().unwrap());
        assert_eq!(warnings.len(), 1, "{:?}", warnings);
        assert!(warnings[0].contains("names chromosome '20', not 'chr20'"), "{}", warnings[0]);
    }

    #[test]
    fn test_an_empty_region_of_a_matching_gvcf_stays_silent() {
        capture::install();
        let dir = test_dir("loh_warn_matching");

        // Matching names, no record in the region, records on another
        // chromosome: an ordinary quiet region.
        let matching = dir.join("matching.vcf");
        std::fs::write(
            &matching,
            "##fileformat=VCFv4.2\n\
             ##contig=<ID=chr20,length=64444167>\n\
             ##contig=<ID=chr21,length=46709983>\n\
             #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002\n\
             chr20\t5001\t.\tA\tG\t50\tPASS\t.\tGT\t0|1\n\
             chr21\t101\t.\tA\tG\t50\tPASS\t.\tGT\t0|1\n",
        )
        .unwrap();
        let snps = load_snps_from_gvcf(matching.to_str().unwrap(), "chr20", 0, 1000).unwrap();
        assert!(snps.het.is_empty());
        assert!(
            capture::warnings_matching(matching.to_str().unwrap()).is_empty(),
            "{:?}",
            capture::warnings_matching(matching.to_str().unwrap()),
        );
    }

    #[test]
    fn test_a_gvcf_read_that_fails_says_the_run_stops() {
        // Without bcftools the read fails at the spawn instead, and that
        // failure's own context ("failed to run bcftools ... (is bcftools in
        // PATH?)") carries the same NextStep::Stop sentence and no pileup
        // one -- so both assertions below hold and the test reports `ok` over
        // the wrong path, never reaching the `bcftools exited with status`
        // error it exists to pin. A false pass is worse than a failure.
        require_bcftools();
        // An error is not an empty result: there is no pileup pass after it.
        let dir = test_dir("loh_read_error");
        let path = dir.join("unindexed.vcf.gz");
        write_gzip(&path, "##fileformat=VCFv4.2\n##contig=<ID=chr20,length=64444167>\n");

        let err = load_snps_from_gvcf(path.to_str().unwrap(), "chr20", 0, 1000).unwrap_err();
        let message = format!("{}", err);
        assert!(message.contains("bcftools exited with status"), "{}", message);
        assert!(message.contains("The run stops here."), "{}", message);
        assert!(!message.contains("Falling back to pileup"), "{}", message);
    }

    #[test]
    fn test_a_renamed_gvcf_that_cannot_be_read_warns_about_the_stop_not_the_pileup() {
        // Without bcftools the read fails at the spawn, before the warning
        // this test is about: say so rather than report an empty vector.
        require_bcftools();
        capture::install();
        let dir = test_dir("loh_read_error_renamed");
        let path = dir.join("renamed_unindexed.vcf.gz");
        write_gzip(&path, "##fileformat=VCFv4.2\n##contig=<ID=20,length=64444167>\n");

        load_snps_from_gvcf(path.to_str().unwrap(), "chr20", 0, 1000).unwrap_err();
        let warnings = capture::warnings_matching(path.to_str().unwrap());
        assert_eq!(warnings.len(), 1, "{:?}", warnings);
        assert!(warnings[0].contains("names chromosome '20', not 'chr20'"), "{}", warnings[0]);
        assert!(warnings[0].contains("The run stops here."), "{}", warnings[0]);
        assert!(!warnings[0].contains("pileup"), "{}", warnings[0]);
    }

    #[test]
    fn test_pileup_calls_het_and_hom_alt_snps() {
        // Region starts at 1000; the reference there is all A.
        let ref_seq = vec![b'A'; 10];
        let counts: HashMap<u64, [u32; 4]> = [
            (1001, [10, 0, 10, 0]), // het A/G
            (1002, [0, 0, 0, 20]),  // hom-alt T
            (1003, [1, 0, 19, 0]),  // hom-alt G, one error read
            (1004, [20, 0, 0, 0]),  // hom-ref
            (1005, [0, 0, 0, 5]),   // too shallow to call
            (1006, [3, 0, 17, 0]),  // 85% G: neither het nor surely hom-alt
        ]
        .into();
        let snps = call_snps(&counts, 1000, &ref_seq);
        assert_eq!(snps.het.iter().map(|s| s.pos).collect::<Vec<_>>(), vec![1001]);
        assert_eq!(snps.hom_alt, [(1002, b'T'), (1003, b'G')].into());
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
        let read_copy = classify_from_collected(&target, &reads);
        assert!(!read_copy.contains_key("tie"));
    }

    #[test]
    fn test_copies_carry_one_haplotype_each_plus_hom_alt() {
        let snps = vec![snp(100), snp(200), snp(300)];
        let hom_alt: HashMap<u64, u8> = [(250, b'T')].into();
        let reads = linked_reads();
        let h1 = (b'A', b'G', b'G');
        let h2 = (b'G', b'A', b'A');
        let mut seen_h1 = false;
        for seed in 0..20 {
            let mut rng = StdRng::seed_from_u64(seed);
            let c = copies_from_snps(&snps, &hom_alt, &reads, &mut rng);
            let event = (c.event_copy[&100], c.event_copy[&200], c.event_copy[&300]);
            let other = (c.other_copy[&100], c.other_copy[&200], c.other_copy[&300]);
            assert!(
                (event, other) == (h1, h2) || (event, other) == (h2, h1),
                "seed {}: event {:?}, other {:?}",
                seed,
                event,
                other
            );
            seen_h1 |= event == h1;
            // Both copies carry the hom-alt allele.
            assert_eq!((c.event_copy[&250], c.other_copy[&250]), (b'T', b'T'));
            // A fragment is on the event copy when it carries its haplotype.
            assert_eq!(c.read_copy["h1a_0"], event == h1);
            assert_eq!(c.read_copy["h2b_3"], event == h2);
        }
        assert!(seen_h1, "the event copy should sometimes be haplotype 1");
    }

    // --- N4: the CRAM pileup must not read another contig's records ---

    /// The shared two-contig CRAM fixture in a scratch directory of its own.
    /// Returns `(dir, fasta_path, cram_path)`; the caller removes `dir`.
    fn two_contig_cram(tag: &str) -> (std::path::PathBuf, String, String) {
        let dir =
            std::env::temp_dir().join(format!("spike_test_loh_{}_{}", tag, std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        let (fasta, cram) = crate::extract::test_fixtures::write_two_contig_cram(&dir);
        (dir, fasta, cram)
    }

    /// Bases counted over the three intervals only chrB covers (0-based,
    /// half-open): all 400 of chrB's aligned bases land in them.
    fn foreign_bases(counts: &HashMap<u64, [u32; 4]>) -> u32 {
        (300..400u64)
            .chain(500..600)
            .chain(700..800)
            .filter_map(|pos| counts.get(&pos))
            .map(|c| c.iter().sum::<u32>())
            .sum()
    }

    #[test]
    fn test_count_alleles_cram_skips_another_contigs_reads() {
        // chrA and chrB cover disjoint intervals in the fixture, so any
        // count a chrA pileup reports over chrB's 301-400, 501-600 and
        // 701-800 is a base noodles-cram's `Query` handed over from the
        // shared container without comparing the record's reference id (L2).
        // This is the pileup het calls and phasing are built from, so a leak
        // here reaches the truth set: assert on the counts, not on what LOH
        // later decides.
        let (dir, fasta, cram) = two_contig_cram("count_alleles");

        let counts = count_alleles_cram(&cram, "chrA", 0, 2000, 20, &fasta).unwrap();

        assert_eq!(
            foreign_bases(&counts),
            0,
            "chrB's bases must not enter the chrA pileup"
        );
        // chrA's own six reads are still counted in full, so this is a filter
        // and not an empty result.
        let total: u32 = counts.values().map(|c| c.iter().sum::<u32>()).sum();
        assert_eq!(total, 600, "every chrA base must still be counted");

        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn test_collect_snp_alleles_cram_skips_another_contigs_reads() {
        // The second LOH pass reads each fragment's base at the het SNPs. A
        // chrB fragment that reaches it gets a haplotype verdict of its own
        // and then helps decide which chrA originals are suppressed.
        let (dir, fasta, cram) = two_contig_cram("collect_alleles");

        // 250 sits under chrA's own reads, 350 under chrB's only.
        let positions: HashSet<u64> = [250u64, 350].into_iter().collect();
        let alleles =
            collect_snp_alleles_cram(&cram, "chrA", 0, 2000, 20, &positions, &fasta).unwrap();

        let mut names: Vec<String> = alleles.keys().cloned().collect();
        names.sort();
        assert_eq!(
            names,
            vec!["chrA_pair0", "chrA_pair1", "chrA_pair2"],
            "only chrA fragments may be classified for a chrA region"
        );
        assert!(
            alleles
                .values()
                .all(|bases| bases.iter().all(|&(pos, _)| pos != 350)),
            "no fragment may report a base at a position only chrB covers"
        );

        let _ = std::fs::remove_dir_all(&dir);
    }
}
