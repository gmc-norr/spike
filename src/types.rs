//! Core types for spike.

/// How the two sides of a fusion are joined.
///
/// Each breakpoint `bp` is a cut between 0-based bases `bp - 1` and `bp`.
/// The join says which side of each cut is kept, and which piece is
/// reverse-complemented.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FusionJoin {
    /// A left of `bp_a`, then B right of `bp_b` (same strand). BND `t[p[`.
    Forward,
    /// A left of `bp_a`, then the reverse complement of B left of `bp_b`.
    /// BND `t]p]`.
    LeftLeft,
    /// The reverse complement of A right of `bp_a`, then B right of `bp_b`.
    /// BND `[p[t`.
    RightRight,
}

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
        bp_a: u64, // cut in gene A, between 0-based bases bp_a-1 and bp_a
        gene_a: String,
        chrom_b: String,
        bp_b: u64, // cut in gene B; `join` says which side of each cut is kept
        gene_b: String,
        allele_fraction: Option<f64>, // per-event AF override (None = use global)
        join: FusionJoin,
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
        ins_seq: Option<Vec<u8>>, // inserted sequence, uppercase (None = generate random)
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
    /// Whether this event's donor pool is drawn from more than one locus.
    ///
    /// `extract_pool_for_event` searches two windows for a multi-locus event
    /// and one for everything else, and `coverage_for_tiling` demands
    /// donor coverage at *every* breakpoint side of a multi-locus event but
    /// only *somewhere* around a single-locus one. Those two rules have to
    /// agree, and they used to be two independent `SimEvent::Fusion` patterns
    /// in two files with nothing linking them: a new multi-locus event type
    /// would silently take the permissive branch (an N5-class hole). The
    /// match below is exhaustive, so a new variant will not compile until
    /// someone answers this question for it.
    pub fn is_multi_locus(&self) -> bool {
        match self {
            SimEvent::Fusion { .. } => true,
            SimEvent::Deletion { .. }
            | SimEvent::Duplication { .. }
            | SimEvent::Inversion { .. }
            | SimEvent::Insertion { .. }
            | SimEvent::SmallVariant { .. } => false,
        }
    }

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
    pub read_length: usize, // cycles: the BAM's most common read length (bam_stats)
    pub min_mapq: u8,
    /// Optional gVCF file for LOH simulation (het SNP positions).
    pub gvcf_path: Option<String>,
    /// Indel error rate: fraction of sequencing errors that are indels (vs substitutions).
    /// 0.0 = substitution-only errors (default), ~0.05 = typical Illumina.
    pub indel_error_rate: f64,
    /// Duplication model: "full" (full tandem haplotype) or "junction" (legacy junction-only).
    pub dup_model: String,
    /// Whether the input library's reads were adapter-trimmed before
    /// alignment (`bam_stats`); see `SynthReadGenerator::with_adapter_trim`.
    pub adapter_trimmed: bool,
}

impl SimConfig {
    /// The shortest fragment the donor fragment model keeps: the one the
    /// generator draws down to (`SynthReadGenerator::min_fragment_len`).
    pub fn min_fragment_len(&self) -> usize {
        if self.adapter_trimmed {
            1
        } else {
            self.read_length
        }
    }
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
    /// Names of the original pairs suppressed by this event. Needed when
    /// combining events: another event's pool may contain the same pair and
    /// would otherwise pass it through as a kept original.
    pub suppressed_names: Vec<String>,
    /// The fraction this event's tiled fragments actually plant, when the
    /// additive cap or the two-fragment floor moved it off the requested VAF;
    /// `None` when the request stands. The truth VCF records it as `SIM_VAF`.
    pub adjusted_vaf: Option<f64>,
    /// Breakpoint sides of this event that the donor pool has no reads over,
    /// `chrom:pos`, deduplicated and in haplotype order. Empty for an event
    /// whose every side is covered. A single-locus event is kept when at
    /// least one side is covered, so this is the part of its footprint whose
    /// synthetic reads have no donor depth behind them.
    pub uncovered_breakpoint_sides: Vec<String>,
    /// How far the donor's depth, where this event's fragments are drawn,
    /// departs from the one depth they are all scaled by (CR2).
    pub depth_fold: DepthFold,
    /// Under `--edit-model origin`, this event's chance of removing each
    /// fragment it can remove. `origin::decide` sums them over every event
    /// and draws once per fragment. Empty under `clean`.
    pub origin_chances: Vec<crate::origin::Chance>,
}

/// The largest fold between the donor's depth in a bin an event's fragments
/// are drawn from and the one depth the tiling scales them all by (CR2). The
/// truth VCF records `fold` as `SIM_DEPTH_FOLD`.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct DepthFold {
    /// `max((D+1)/(C+1), (C+1)/(D+1))` over the bins; the +1 keeps an empty
    /// bin finite.
    pub fold: f64,
    /// `C`: the depth every fragment is scaled by.
    pub scaled_by: f64,
    /// The bin with the largest fold, `chrom:start-end` (0-based, half-open).
    pub worst_bin: String,
    /// `D` in that bin.
    pub worst_depth: f64,
}
