//! `--edit-model origin`: which original reads came from an event's edited
//! copy, wherever the aligner put them.
//!
//! The default (`clean`) edits only donor-pool pairs: both mates at
//! `--min-mapq` or above, a proper pair, no duplicate, secondary,
//! supplementary or QC-fail flag. Where the aligner cannot tell a spot from a
//! look-alike most reads fail that, so the event is barely planted (RF8).
//! `origin` gives every primary read a chance of having come from the
//! event's footprint, read from its MAPQ and its `XA` hits, and removes it
//! by that chance. The design is
//! `docs/superpowers/specs/2026-09-26-edit-model-origin-design.md`; comments
//! name its fixes R1-R5.

use noodles::sam::alignment::record::cigar::Op;

/// A reference interval, 0-based half-open.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Span {
    pub chrom: String,
    pub start: u64,
    pub end: u64,
}

impl Span {
    pub fn new(chrom: &str, start: u64, end: u64) -> Self {
        Span {
            chrom: chrom.to_string(),
            start,
            end,
        }
    }

    /// Whether `other` lies wholly inside this span.
    pub fn holds(&self, other: &Span) -> bool {
        self.chrom == other.chrom && other.start >= self.start && other.end <= self.end
    }

    /// Whether `other` shares at least one base with this span.
    pub fn overlaps(&self, other: &Span) -> bool {
        self.chrom == other.chrom && other.start < self.end && self.start < other.end
    }
}

impl std::fmt::Display for Span {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}:{}-{}", self.chrom, self.start, self.end)
    }
}

/// One place a read may have come from, and the chance that it did.
#[derive(Debug, Clone, PartialEq)]
pub struct Placement {
    pub span: Span,
    pub chance: f64,
}

/// bwa-mem lists alternative hits in `XA` only when there are at most this
/// many (its `-h 5` default), so a MAPQ 0 read without `XA` has more.
const MAX_LISTED_HITS: usize = 5;

/// The chance a read came from where the aligner put it.
///
/// MAPQ 1 or more is the aligner's own estimate, `1 - 10^(-MAPQ/10)`. At
/// MAPQ 0 the hits are equally good: with `n_xa` listed alternatives the
/// primary is one of `1 + n_xa`; with none listed there are more than
/// [`MAX_LISTED_HITS`], so it is one of at least 6.
pub fn p_here(mapq: u8, n_xa: usize) -> f64 {
    if mapq > 0 {
        1.0 - 10f64.powf(-(mapq as f64) / 10.0)
    } else if n_xa > 0 {
        1.0 / (1 + n_xa) as f64
    } else {
        1.0 / (1 + MAX_LISTED_HITS) as f64
    }
}

/// The reference bases a CIGAR covers.
fn reference_length(ops: &[Op]) -> u64 {
    ops.iter()
        .filter(|op| op.kind().consumes_reference())
        .map(|op| op.len() as u64)
        .sum()
}

/// The reference bases a CIGAR string such as `100M1D51M` covers, or `None`
/// when it is malformed.
fn cigar_reference_length(cigar: &str) -> Option<u64> {
    let (mut len, mut n, mut digits) = (0u64, 0u64, false);
    for c in cigar.bytes() {
        if c.is_ascii_digit() {
            n = n * 10 + u64::from(c - b'0');
            digits = true;
            continue;
        }
        if !digits {
            return None;
        }
        match c {
            b'M' | b'D' | b'N' | b'=' | b'X' => len += n,
            b'I' | b'S' | b'H' | b'P' => {}
            _ => return None,
        }
        n = 0;
        digits = false;
    }
    (!digits).then_some(len)
}

/// The alternative hits in an `XA:Z` value: `chrom,±pos,CIGAR,NM;` each, with
/// `pos` 1-based and its sign the strand. Malformed entries are skipped.
pub fn parse_xa(xa: &str) -> Vec<Span> {
    xa.split(';')
        .filter(|entry| !entry.is_empty())
        .filter_map(|entry| {
            let mut fields = entry.split(',');
            let chrom = fields.next()?;
            let pos: i64 = fields.next()?.parse().ok()?;
            let len = cigar_reference_length(fields.next()?)?;
            let start = pos.unsigned_abs().checked_sub(1)?;
            Some(Span::new(chrom, start, start + len))
        })
        .collect()
}

/// A read's placements: its primary at [`p_here`], then each `XA` hit at an
/// equal share of the rest, `(1 - p_here) / k`.
pub fn placements(primary: Span, mapq: u8, alternatives: &[Span]) -> Vec<Placement> {
    let here = p_here(mapq, alternatives.len());
    let each = if alternatives.is_empty() {
        0.0
    } else {
        (1.0 - here) / alternatives.len() as f64
    };
    std::iter::once(Placement {
        span: primary,
        chance: here,
    })
    .chain(alternatives.iter().map(|span| Placement {
        span: span.clone(),
        chance: each,
    }))
    .collect()
}

/// The chance a read came from inside `footprint`: the sum over its
/// placements lying wholly inside it, the primary included (R2).
pub fn chance_within(placements: &[Placement], footprint: &Span) -> f64 {
    placements
        .iter()
        .filter(|p| footprint.holds(&p.span))
        .map(|p| p.chance)
        .sum()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64) -> bool {
        (a - b).abs() < 1e-9
    }

    #[test]
    fn test_p_here_follows_mapq_and_the_xa_count() {
        assert!(close(p_here(60, 0), 1.0 - 1e-6));
        assert!(close(p_here(3, 0), 1.0 - 10f64.powf(-0.3)));
        assert!(close(p_here(0, 1), 0.5));
        // No XA at MAPQ 0: bwa-mem lists at most 5, so there are more.
        assert!(close(p_here(0, 0), 1.0 / 6.0));
    }

    #[test]
    fn test_parse_xa_reads_chrom_start_strand_and_cigar_length() {
        let hits = parse_xa("cluster,+2002,151M,0;cluster,-1002,100M1D51M,1;");
        assert_eq!(
            hits,
            vec![Span::new("cluster", 2001, 2152), Span::new("cluster", 1001, 1153)]
        );
    }

    #[test]
    fn test_parse_xa_skips_a_malformed_entry() {
        assert_eq!(parse_xa("chr1,+10,5Q,0;chr1,+20,5M,0;"), vec![Span::new("chr1", 19, 24)]);
    }

    #[test]
    fn test_a_read_on_a_remote_copy_with_three_hits_in_the_footprint_gets_three_quarters() {
        // R2: bwa-mem put `cluster_1` on the remote copy at MAPQ 0 with all
        // three XA hits in one cluster (the review's probe).
        let alternatives = parse_xa("cluster,+2002,151M,0;cluster,+1002,151M,0;cluster,+1502,151M,0;");
        let placed = placements(Span::new("remote", 1001, 1152), 0, &alternatives);
        assert!(close(chance_within(&placed, &Span::new("cluster", 1000, 2400)), 0.75));
    }

    #[test]
    fn test_a_mapq0_read_whose_one_hit_is_also_in_the_footprint_gets_one() {
        let placed = placements(Span::new("chr1", 100, 250), 0, &[Span::new("chr1", 600, 750)]);
        assert!(close(chance_within(&placed, &Span::new("chr1", 0, 1000)), 1.0));
    }

    #[test]
    fn test_a_placement_sticking_out_of_the_footprint_does_not_count() {
        let placed = placements(Span::new("chr1", 950, 1100), 60, &[]);
        assert!(close(chance_within(&placed, &Span::new("chr1", 0, 1000)), 0.0));
    }
}
