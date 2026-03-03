//! Core types for spike.

/// A simulated event specification.
#[derive(Debug, Clone)]
pub enum SimEvent {
    /// Exon deletion: contiguous region removed from a single chromosome.
    Deletion {
        chrom: String,
        del_start: u64, // 0-based, inclusive
        del_end: u64,   // 0-based, exclusive (half-open)
        gene: String,
        exons: Vec<String>,           // affected exon names for annotation
        allele_fraction: Option<f64>, // per-event AF override (None = use global)
    },
    /// Gene fusion: two breakpoints on potentially different chromosomes joined.
    Fusion {
        chrom_a: String,
        bp_a: u64, // 0-based breakpoint in gene A (first base NOT included)
        gene_a: String,
        chrom_b: String,
        bp_b: u64, // 0-based breakpoint in gene B (first base included on right side)
        gene_b: String,
        allele_fraction: Option<f64>, // per-event AF override (None = use global)
        inverted: bool,               // true for reverse-complement join at gene B
    },
    /// Tandem duplication: region is duplicated in place.
    Duplication {
        chrom: String,
        dup_start: u64, // 0-based, inclusive
        dup_end: u64,   // 0-based, exclusive (half-open)
        gene: String,
        allele_fraction: Option<f64>,
    },
    /// Inversion: region is reversed in place.
    Inversion {
        chrom: String,
        inv_start: u64, // 0-based, inclusive
        inv_end: u64,   // 0-based, exclusive (half-open)
        gene: String,
        allele_fraction: Option<f64>,
    },
    /// Insertion: novel sequence inserted at a position.
    Insertion {
        chrom: String,
        pos: u64,                 // 0-based insertion point
        ins_seq: Option<Vec<u8>>, // explicit inserted sequence (None = generate random)
        ins_len: u64,             // length of insertion
        gene: String,
        allele_fraction: Option<f64>,
    },
    /// Small variant: SNP, MNV, or small indel represented by explicit REF/ALT alleles.
    SmallVariant {
        chrom: String,
        pos: u64,            // 0-based start position
        ref_allele: Vec<u8>, // reference allele bases (uppercase)
        alt_allele: Vec<u8>, // alternate allele bases (uppercase)
        gene: String,
        allele_fraction: Option<f64>,
    },
}

impl SimEvent {
    /// Get the per-event allele fraction override, if any.
    pub fn allele_fraction(&self) -> Option<f64> {
        match self {
            SimEvent::Deletion {
                allele_fraction, ..
            }
            | SimEvent::Fusion {
                allele_fraction, ..
            }
            | SimEvent::Duplication {
                allele_fraction, ..
            }
            | SimEvent::Inversion {
                allele_fraction, ..
            }
            | SimEvent::Insertion {
                allele_fraction, ..
            }
            | SimEvent::SmallVariant {
                allele_fraction, ..
            } => *allele_fraction,
        }
    }

    /// Returns (chrom, start, end) for single-region events.
    /// Returns None for Fusion (which has two regions).
    pub fn primary_region(&self) -> Option<(&str, u64, u64)> {
        match self {
            SimEvent::Deletion {
                chrom,
                del_start,
                del_end,
                ..
            } => Some((chrom, *del_start, *del_end)),
            SimEvent::Duplication {
                chrom,
                dup_start,
                dup_end,
                ..
            } => Some((chrom, *dup_start, *dup_end)),
            SimEvent::Inversion {
                chrom,
                inv_start,
                inv_end,
                ..
            } => Some((chrom, *inv_start, *inv_end)),
            SimEvent::Insertion { chrom, pos, .. } => Some((chrom, *pos, *pos)),
            SimEvent::SmallVariant {
                chrom,
                pos,
                ref_allele,
                ..
            } => Some((chrom, *pos, *pos + ref_allele.len() as u64)),
            SimEvent::Fusion { .. } => None,
        }
    }

    /// Set the per-event allele fraction.
    pub fn set_allele_fraction(&mut self, af: Option<f64>) {
        match self {
            SimEvent::Deletion {
                allele_fraction, ..
            }
            | SimEvent::Fusion {
                allele_fraction, ..
            }
            | SimEvent::Duplication {
                allele_fraction, ..
            }
            | SimEvent::Inversion {
                allele_fraction, ..
            }
            | SimEvent::Insertion {
                allele_fraction, ..
            }
            | SimEvent::SmallVariant {
                allele_fraction, ..
            } => *allele_fraction = af,
        }
    }
}

/// Configuration for the simulation.
#[derive(Debug, Clone)]
pub struct SimConfig {
    pub bam_path: String,
    pub ref_path: String,
    pub allele_fraction: f64, // 0.0 to 1.0
    pub flank_bp: u64,
    pub read_length: usize, // derived from BAM stats
    pub min_mapq: u8,
    /// Optional gVCF file for LOH simulation (het SNP positions).
    pub gvcf_path: Option<String>,
    /// Indel error rate: fraction of sequencing errors that are indels (vs substitutions).
    /// 0.0 = substitution-only errors (default), ~0.05 = typical Illumina.
    pub indel_error_rate: f64,
    /// Duplication model: "full" (full tandem haplotype) or "junction" (legacy junction-only).
    pub dup_model: String,
}

/// A read pair extracted from BAM, stored in FASTQ-ready form.
///
/// Sequences and qualities are in the ORIGINAL sequencing orientation (FASTQ order).
/// When extracted from BAM, reverse-complemented reads are flipped back so that
/// cycle position 0 = first base sequenced.
#[derive(Debug, Clone)]
pub struct ReadPair {
    pub name: String,
    pub seq1: Vec<u8>,  // read1 sequence (ASCII: A/C/G/T/N)
    pub qual1: Vec<u8>, // read1 phred+33 quality (ASCII)
    pub seq2: Vec<u8>,  // read2 sequence
    pub qual2: Vec<u8>, // read2 phred+33 quality
    /// Leftmost alignment position of fragment (0-based).
    pub ref_start: u64,
    /// Rightmost extent of fragment (0-based, exclusive).
    pub ref_end: u64,
    /// Template length from BAM.
    pub insert_size: i64,
    pub chrom: String,
}

/// Pool of real reads from a genomic region.
pub struct ReadPool {
    /// All read pairs, sorted by ref_start.
    pub pairs: Vec<ReadPair>,
    /// Fragment length distribution derived from these reads.
    pub frag_dist: crate::stats::FragmentDist,
}

/// Result of chimeric read generation for one event.
pub struct SplicedOutput {
    /// Chimeric read pairs to add.
    pub chimeric_pairs: Vec<ReadPair>,
    /// Original read pairs to keep (unmodified passthrough).
    pub kept_originals: Vec<ReadPair>,
    /// Count of original pairs suppressed (replaced or removed).
    pub suppressed_count: usize,
}
