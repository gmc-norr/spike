//! VCF v4.3 input parser for spike.
//!
//! Reads SV records from a VCF file and converts them into SimEvents.
//! Supports DEL, INS, DUP, INV, and BND (fusion) SVTYPE records.

use anyhow::{bail, Context, Result};
use std::collections::HashSet;
use std::io::{BufRead, BufReader};

use crate::types::{FusionJoin, SimEvent};

/// Load simulation events from a VCF file.
///
/// Supports both plain `.vcf` and bgzip-compressed `.vcf.gz` files.
/// BND records are paired by MATEID to produce single Fusion events.
pub fn load_events_from_vcf(path: &str) -> Result<Vec<SimEvent>> {
    let file =
        std::fs::File::open(path).with_context(|| format!("failed to open VCF: {}", path))?;

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
    pos: u64, // 0-based SV start (DEL/DUP/INV/INS) or 0-based breakpoint (BND)
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
            _ => raw_pos,                  // Others: 1-based preceding base == 0-based start
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
                let Some((start, end)) = resolve_sv_span_or_warn(record, "DEL") else {
                    continue;
                };
                let gene = parse_info_field(&record.info, "SIM_GENE")
                    .unwrap_or("unknown")
                    .to_string();
                let af = extract_af(&record.info);

                events.push(SimEvent::Deletion {
                    chrom: record.chrom.clone(),
                    del_start: start,
                    del_end: end,
                    gene,
                    exons: Vec::new(),
                    allele_fraction: af,
                });
            }
            SvTypeTag::Dup => {
                let Some((start, end)) = resolve_sv_span_or_warn(record, "DUP") else {
                    continue;
                };
                let gene = parse_info_field(&record.info, "SIM_GENE")
                    .unwrap_or("unknown")
                    .to_string();
                let af = extract_af(&record.info);

                events.push(SimEvent::Duplication {
                    chrom: record.chrom.clone(),
                    dup_start: start,
                    dup_end: end,
                    gene,
                    allele_fraction: af,
                });
            }
            SvTypeTag::Inv => {
                let Some((start, end)) = resolve_sv_span_or_warn(record, "INV") else {
                    continue;
                };
                let gene = parse_info_field(&record.info, "SIM_GENE")
                    .unwrap_or("unknown")
                    .to_string();
                let af = extract_af(&record.info);

                events.push(SimEvent::Inversion {
                    chrom: record.chrom.clone(),
                    inv_start: start,
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

                // ALT with explicit sequence (not symbolic <INS>) carries the
                // inserted bases: the ones it adds to the flanks it shares
                // with REF. REF=AT ALT=ATGGG inserts GGG after the T; taking
                // ALT past its first base made it a 4 bp insertion of TGGG.
                let (ins_pos, ins_seq) = if alt_carries_sequence(&record.alt) {
                    let trimmed = trim_shared_flanks(&record.ref_allele, &record.alt);
                    if !trimmed.ref_rem.is_empty() || trimmed.alt_rem.is_empty() {
                        // REF keeps bases ALT drops: a complex record, not an
                        // insertion. Which bases are inserted and which are
                        // replaced is a guess, so drop it rather than write a
                        // plausible-looking wrong truth record.
                        log::warn!("{}", not_an_insertion_warning(record));
                        continue;
                    }
                    let start = event_start(record.pos, trimmed.prefix);
                    if start == 0 {
                        // No base 0 to anchor the insertion to, exactly as
                        // for the spans resolve_sv_span_or_warn rejects.
                        log::warn!("{}", before_first_base_warning(record, "INS"));
                        continue;
                    }
                    (start, Some(trimmed.alt_rem.to_vec()))
                } else {
                    (record.pos, None)
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
                    pos: ins_pos,
                    ins_seq,
                    ins_len: effective_len,
                    gene,
                    allele_fraction: af,
                });
            }
            SvTypeTag::SmallVar => {
                // REF and ALT must differ, case-insensitively: alleles that
                // are the same base but differ only in case (e.g. "A"/"a")
                // would otherwise become a no-op "variant" written straight
                // to the truth VCF. exon.rs's snp: spec already rejects
                // this; the VCF ingest path did not check at all.
                if record.ref_allele.eq_ignore_ascii_case(&record.alt) {
                    bail!(
                        "record {} has identical REF and ALT alleles: '{}'",
                        record.id,
                        record.ref_allele
                    );
                }

                let af = extract_af(&record.info);
                let gene = parse_info_field(&record.info, "SIM_GENE")
                    .unwrap_or("unknown")
                    .to_string();

                // Uppercase, as types.rs documents SmallVariant's alleles and
                // as exon.rs's snp: spec already stores them: synthesis
                // uppercases either way, so the stored case only ever reaches
                // the truth VCF text, and the same variant must not be written
                // one way via --vcf and another via --event.
                events.push(SimEvent::SmallVariant {
                    chrom: record.chrom.clone(),
                    pos: record.pos,
                    ref_allele: record.ref_allele.to_ascii_uppercase().into_bytes(),
                    alt_allele: record.alt.to_ascii_uppercase().into_bytes(),
                    gene,
                    allele_fraction: af,
                });
            }
            SvTypeTag::Bnd => {
                // Skip if already processed as part of a MATEID pair.
                // Never skip records with placeholder ID "." — they are not
                // unique identifiers and collapsing them would drop distinct BNDs.
                let has_real_id = record.id != ".";
                if has_real_id && bnd_processed.contains(&record.id) {
                    continue;
                }

                // record.pos is 0-based for BND; decode_bnd takes the VCF POS.
                let bnd = decode_bnd(&record.chrom, record.pos + 1, &record.alt).with_context(|| {
                    format!("failed to parse BND ALT '{}' for record {}", record.alt, record.id)
                })?;

                // Mark both this record and its mate as processed.
                if has_real_id {
                    bnd_processed.insert(record.id.clone());
                    if let Some(mate_id) = parse_info_field(&record.info, "MATEID") {
                        bnd_processed.insert(mate_id.to_string());
                    }
                }

                // Extract gene info.
                let gene_a = parse_info_field(&record.info, "GENE_A")
                    .unwrap_or("unknown")
                    .to_string();
                let gene_b = parse_info_field(&record.info, "GENE_B")
                    .unwrap_or("unknown")
                    .to_string();

                let (gene_a, gene_b) = if bnd.record_is_a {
                    (gene_a, gene_b)
                } else {
                    (gene_b, gene_a)
                };

                let af = extract_af(&record.info);

                events.push(SimEvent::Fusion {
                    chrom_a: bnd.chrom_a,
                    bp_a: bnd.bp_a,
                    gene_a,
                    chrom_b: bnd.chrom_b,
                    bp_b: bnd.bp_b,
                    gene_b,
                    allele_fraction: af,
                    join: bnd.join,
                });
            }
        }
    }

    Ok(events)
}

/// A BND record decoded into fusion coordinates (see [`FusionJoin`] for
/// what `bp_a`, `bp_b` and `join` mean).
#[derive(Debug, PartialEq)]
struct BndFusion {
    chrom_a: String,
    bp_a: u64,
    chrom_b: String,
    bp_b: u64,
    join: FusionJoin,
    /// False when the record's own side became gene B.
    record_is_a: bool,
}

/// Decode a BND record at `chrom`:`pos` (1-based VCF POS) with ALT `alt`.
///
/// The four VCF BND forms, with `t` the record's base and `p` the partner:
///   `t[p[` — t, then the partner from p rightwards          → Forward
///   `]p]t` — the partner up to p, then t rightwards          → Forward (record is B)
///   `t]p]` — t, then the partner up to p, reverse-complemented → LeftLeft
///   `[p[t` — the partner from p rightwards, reverse-complemented, then t
///            rightwards; the same molecule as the record's right side
///            reversed followed by the partner's right side     → RightRight
fn decode_bnd(chrom: &str, pos: u64, alt: &str) -> Result<BndFusion> {
    let bracket = alt
        .chars()
        .find(|c| *c == '[' || *c == ']')
        .ok_or_else(|| anyhow::anyhow!("no bracket in BND ALT '{}'", alt))?;
    let open = alt.find(bracket).unwrap();
    let close = alt.rfind(bracket).unwrap();
    if open == close {
        bail!("no bracket pair found in BND ALT '{}'", alt);
    }

    let (partner_chrom, partner_pos) = alt[open + 1..close]
        .rsplit_once(':')
        .ok_or_else(|| anyhow::anyhow!("no chr:pos found in BND ALT '{}'", alt))?;
    let partner_pos: u64 = partner_pos
        .parse()
        .with_context(|| format!("invalid position in BND ALT '{}'", alt))?;
    if partner_pos == 0 || pos == 0 {
        bail!("BND positions must be >= 1 in '{}'", alt);
    }

    // A cut between 0-based bases bp-1 and bp. The 1-based base x is 0-based
    // x-1, so keeping it as the last base of a left side means bp = x, and
    // keeping it as the first base of a right side means bp = x - 1.
    let (record, partner) = (chrom.to_string(), partner_chrom.to_string());
    let t_first = open > 0;
    let decoded = match (t_first, bracket) {
        (true, '[') => (record, pos, partner, partner_pos - 1, FusionJoin::Forward, true),
        (false, ']') => (partner, partner_pos, record, pos - 1, FusionJoin::Forward, false),
        (true, ']') => (record, pos, partner, partner_pos, FusionJoin::LeftLeft, true),
        (false, '[') => (record, pos - 1, partner, partner_pos - 1, FusionJoin::RightRight, true),
        _ => unreachable!("bracket is '[' or ']'"),
    };
    let (chrom_a, bp_a, chrom_b, bp_b, join, record_is_a) = decoded;
    Ok(BndFusion {
        chrom_a,
        bp_a,
        chrom_b,
        bp_b,
        join,
        record_is_a,
    })
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
        if let Some(value) = field
            .strip_prefix(key)
            .and_then(|rest| rest.strip_prefix('='))
        {
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

/// True when ALT carries sequence of its own, rather than being a symbolic
/// allele (`<DEL>`) or the bare anchor base. Same test the INS arm uses to
/// find an explicit inserted sequence.
fn alt_carries_sequence(alt: &str) -> bool {
    !alt.starts_with('<') && alt.len() > 1
}

/// REF and ALT with the flanks they share stripped off: the common prefix
/// first, then any common suffix the leftovers still share. `prefix` counts
/// the leading bases both alleles keep, so the last base the record leaves
/// untouched is 1-based POS + `prefix` - 1 and the event starts after it.
struct TrimmedAlleles<'a> {
    prefix: usize,
    ref_rem: &'a [u8],
    alt_rem: &'a [u8],
}

/// Strip the flanks REF and ALT share. Case-insensitive, so a soft-masked
/// REF is trimmed against an uppercase ALT like any other.
fn trim_shared_flanks<'a>(ref_allele: &'a str, alt: &'a str) -> TrimmedAlleles<'a> {
    let (r, a) = (ref_allele.as_bytes(), alt.as_bytes());
    let prefix = r
        .iter()
        .zip(a.iter())
        .take_while(|(x, y)| x.eq_ignore_ascii_case(y))
        .count();
    // Only bases past the shared prefix can also be a shared suffix:
    // REF=AT ALT=AGGGT is the 3 bp insertion GGG, not T -> GGGT.
    let (mut r_end, mut a_end) = (r.len(), a.len());
    while r_end > prefix && a_end > prefix && r[r_end - 1].eq_ignore_ascii_case(&a[a_end - 1]) {
        r_end -= 1;
        a_end -= 1;
    }
    TrimmedAlleles {
        prefix,
        ref_rem: &r[prefix..r_end],
        alt_rem: &a[prefix..a_end],
    }
}

/// The 0-based start of an event whose alleles share `prefix` leading bases:
/// the base just after the last shared one. `pos` is the record's VCF POS,
/// which the parser guarantees is >= 1, so alleles sharing nothing — an
/// equal-length substitution, which carries no padding base — start at POS
/// itself rather than at the base after it.
fn event_start(pos: u64, prefix: usize) -> u64 {
    pos + prefix as u64 - 1
}

/// True when `alt` is `ref_seq` reverse-complemented: the allele shape an
/// inversion written as an equal-length substitution has.
fn is_reverse_complement(ref_seq: &[u8], alt: &[u8]) -> bool {
    let mut rc = ref_seq.to_vec();
    crate::extract::reverse_complement(&mut rc);
    rc.eq_ignore_ascii_case(alt)
}

/// Resolve a DEL/DUP/INV record's span as 0-based `(start, end)`: prefer
/// INFO/END, then INFO/SVLEN, then the alleles. With END or SVLEN the start
/// is POS (the preceding base); read off the alleles it is wherever the
/// bases the alleles share end.
fn resolve_sv_span(record: &SvRecord) -> Option<(u64, u64)> {
    if let Some(end) = parse_info_u64(&record.info, "END") {
        return Some((record.pos, end));
    }
    if let Some(svlen) = parse_info_i64(&record.info, "SVLEN") {
        return Some((record.pos, record.pos + svlen.unsigned_abs()));
    }
    span_from_alleles(record)
}

/// Read a span off the alleles, for the shapes that state one unambiguously.
/// Anything else — a complex record, or one with no length anywhere —
/// returns `None` and the caller rejects it, rather than turning a malformed
/// record into a plausible-looking wrong truth record.
fn span_from_alleles(record: &SvRecord) -> Option<(u64, u64)> {
    if record.alt.starts_with('<') {
        // A symbolic ALT (<DEL>) carries no sequence to line REF up against,
        // so REF is the anchor plus the affected bases and the span is what
        // REF has past the anchor. A single base is *not* such an ALT: it is
        // the anchor only when it is REF's first base, and REF=ACGT ALT=T
        // deletes ACG at POS-1, not CGT at POS.
        return (record.ref_allele.len() > 1)
            .then(|| (record.pos, record.pos + record.ref_allele.len() as u64 - 1));
    }

    // An inversion spelled out base for base states its own span, so it is
    // read whole rather than trimmed; see [`inv_span`].
    if record.sv_type == SvTypeTag::Inv && record.ref_allele.len() == record.alt.len() {
        return inv_span(record);
    }

    let trimmed = trim_shared_flanks(&record.ref_allele, &record.alt);
    let start = event_start(record.pos, trimmed.prefix);
    let span = |len: usize| Some((start, start + len as u64));
    let (ref_len, alt_len) = (trimmed.ref_rem.len(), trimmed.alt_rem.len());

    if alt_len == 0 && ref_len > 0 {
        // REF keeps bases ALT drops. They are reference sequence at
        // [start, start + ref_len) however long the shared prefix is, so
        // REF=GTTAAAGTTTATCAGAAAATT ALT=GTTAAAG is the 14 bp event after
        // the 7 shared bases, not a 21 bp one at POS.
        return span(ref_len);
    }
    if ref_len == 0 && alt_len > 0 {
        // ALT is the anchor plus a duplicated copy. Only a bare single-base
        // REF says where that copy comes from; with sequence in REF the
        // stripped anchor would land on bases the record never spells out,
        // and the copy could equally be REF's own span.
        let ref_is_the_anchor = record.ref_allele.len() == 1 && trimmed.prefix == 1;
        if record.sv_type == SvTypeTag::Dup && ref_is_the_anchor {
            return span(alt_len);
        }
        return None;
    }
    None
}

/// Read an inversion's span off alleles that spell it out base for base.
/// An inversion's span is *stated* by the record, not derived from where its
/// alleles differ, so both spellings are matched against the whole allele:
/// an equal-length substitution, which carries no padding base and so starts
/// at POS itself, or a padding base followed by the inverted region.
/// Trimming shared flanks here narrows both spellings whenever the inverted
/// region's ends are their own complements, which about one equal-length INV
/// in four has: REF=AGTT ALT=AACT would come out as a 2 bp inversion of the
/// middle where the record spells out 4 bp, and the padded REF=TAGTT
/// ALT=TAACT as the same 2 bp instead of the 4 bp after its anchor.
fn inv_span(record: &SvRecord) -> Option<(u64, u64)> {
    let (r, a) = (record.ref_allele.as_bytes(), record.alt.as_bytes());
    if is_reverse_complement(r, a) {
        // No padding base: POS is the first inverted base, one past the
        // preceding base every other span here counts from. POS >= 1 is
        // guaranteed by the parser; a resulting start of 0 is rejected by
        // the caller, since no VCF record can name the base before the first.
        let start = record.pos - 1;
        return Some((start, start + r.len() as u64));
    }
    if r.len() > 1 && r[0].eq_ignore_ascii_case(&a[0]) && is_reverse_complement(&r[1..], &a[1..]) {
        // Padded: POS is the anchor, the inverted region follows it.
        return Some((record.pos, record.pos + r.len() as u64 - 1));
    }
    None
}

/// Resolve a DEL/DUP/INV span, warning when the record has to be dropped.
/// Shared by the three arms so a rejection reads the same whatever the type.
fn resolve_sv_span_or_warn(record: &SvRecord, sv_type: &str) -> Option<(u64, u64)> {
    match resolve_sv_span(record) {
        // Reading the span off the alleles can move the start left of POS,
        // and at POS=1 that is base 0 — not a position a VCF record can
        // name, so the truth record would come out as POS=0 with REF=N.
        Some((0, _)) => {
            log::warn!("{}", before_first_base_warning(record, sv_type));
            None
        }
        span @ Some(_) => span,
        None => {
            log::warn!("{}", no_length_warning(record, sv_type));
            None
        }
    }
}

/// Message logged when a DEL/DUP/INV record's span cannot be resolved.
/// Names chrom:pos as well as the ID, because SV records routinely carry
/// `ID=.`, which on its own does not say which record was dropped.
/// `record.pos` is the VCF POS column verbatim for these types.
fn no_length_warning(record: &SvRecord, sv_type: &str) -> String {
    format!(
        "{} record {} at {}:{} has no END and no SVLEN, and its REF/ALT alleles are not a shape \
         a length can be read from; skipping",
        sv_type, record.id, record.chrom, record.pos
    )
}

/// Message logged when an INS record's ALT is not its REF plus inserted
/// bases. Names chrom:pos as well as the ID, for the same reason
/// [`no_length_warning`] does: `ID=.` is common and identifies nothing.
fn not_an_insertion_warning(record: &SvRecord) -> String {
    format!(
        "INS record {} at {}:{} keeps REF bases its ALT drops, so the bases it inserts are \
         not a shape they can be read from; skipping",
        record.id, record.chrom, record.pos
    )
}

/// Message logged when a record's alleles place its event before the
/// chromosome's first base. Names chrom:pos for the same reason the others
/// do, so the dropped record can be found.
fn before_first_base_warning(record: &SvRecord, sv_type: &str) -> String {
    format!(
        "{} record {} at {}:{} has alleles that start the event before the chromosome's first \
         base, which no VCF record can name; skipping",
        sv_type, record.id, record.chrom, record.pos
    )
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

    /// Decode one ALT form for a record at chr1 POS 100, partner chr2:200.
    fn decode(alt: &str) -> BndFusion {
        decode_bnd("chr1", 100, alt).unwrap()
    }

    fn bnd(a: &str, bp_a: u64, b: &str, bp_b: u64, join: FusionJoin, record_is_a: bool) -> BndFusion {
        BndFusion {
            chrom_a: a.to_string(),
            bp_a,
            chrom_b: b.to_string(),
            bp_b,
            join,
            record_is_a,
        }
    }

    #[test]
    fn test_decode_bnd_t_open_p_open() {
        // t[p[: chr1 up to and including base 100 (cut after it), then chr2
        // from base 200 on (cut before it: 0-based 199).
        assert_eq!(
            decode("N[chr2:200["),
            bnd("chr1", 100, "chr2", 199, FusionJoin::Forward, true)
        );
    }

    #[test]
    fn test_decode_bnd_close_p_close_t() {
        // ]p]t: chr2 up to and including base 200, then chr1 from base 100 on.
        assert_eq!(
            decode("]chr2:200]N"),
            bnd("chr2", 200, "chr1", 99, FusionJoin::Forward, false)
        );
    }

    #[test]
    fn test_decode_bnd_t_close_p_close() {
        // t]p]: chr1 up to base 100, then chr2 up to base 200, reversed.
        assert_eq!(
            decode("N]chr2:200]"),
            bnd("chr1", 100, "chr2", 200, FusionJoin::LeftLeft, true)
        );
    }

    #[test]
    fn test_decode_bnd_open_p_open_t() {
        // [p[t: chr2 from base 200 on, reversed, then chr1 from base 100 on.
        // Same molecule as chr1-right reversed + chr2-right.
        assert_eq!(
            decode("[chr2:200[N"),
            bnd("chr1", 99, "chr2", 199, FusionJoin::RightRight, true)
        );
    }

    #[test]
    fn test_decode_bnd_mates_give_the_same_forward_fusion() {
        let a = decode_bnd("chr1", 100, "N[chr2:200[").unwrap();
        let b = decode_bnd("chr2", 200, "]chr1:100]N").unwrap();
        assert_eq!((a.chrom_a, a.bp_a, a.chrom_b, a.bp_b), (b.chrom_a, b.bp_a, b.chrom_b, b.bp_b));
    }

    #[test]
    fn test_truth_bnd_records_decode_to_the_written_fusion() {
        // Each record of a pair written by truth.rs must read back as the same
        // molecule. Symmetric joins may come back with A and B swapped.
        let cases = [
            (100, 199, FusionJoin::Forward),
            (100, 200, FusionJoin::LeftLeft),
            (99, 199, FusionJoin::RightRight),
        ];
        for (bp_a, bp_b, join) in cases {
            let written = ("chr1".to_string(), bp_a, "chr2".to_string(), bp_b);
            let swapped = ("chr2".to_string(), bp_b, "chr1".to_string(), bp_a);
            let [(pos_a, alt_a), (pos_b, alt_b)] =
                crate::truth::bnd_records("chr1", bp_a, "chr2", bp_b, join);
            for (chrom, pos, alt) in [("chr1", pos_a, alt_a), ("chr2", pos_b, alt_b)] {
                let d = decode_bnd(chrom, pos, &alt).unwrap();
                let got = (d.chrom_a, d.bp_a, d.chrom_b, d.bp_b);
                assert_eq!(d.join, join, "{:?}: {} {} {}", join, chrom, pos, alt);
                assert!(
                    got == written || (join != FusionJoin::Forward && got == swapped),
                    "{:?}: {} {} {} decoded to {:?}",
                    join, chrom, pos, alt, got
                );
            }
        }
    }

    #[test]
    fn test_extract_af_sim_vaf() {
        assert_eq!(
            extract_af("SVTYPE=BND;SIM_VAF=0.050;GENE_A=BCR"),
            Some(0.05)
        );
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
        assert_eq!(
            parse_info_field("SVTYPE=DEL;END=100;SVLEN=-500", "END"),
            Some("100")
        );
        assert_eq!(
            parse_info_field("SVTYPE=DEL;END=100;SVLEN=-500", "SVLEN"),
            Some("-500")
        );
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
            SimEvent::SmallVariant {
                chrom,
                pos,
                ref_allele,
                alt_allele,
                ..
            } => {
                assert_eq!(chrom, "chr1");
                assert_eq!(*pos, 99);
                assert_eq!(ref_allele, b"A");
                assert_eq!(alt_allele, b"T");
            }
            _ => panic!("expected SmallVariant"),
        }
    }

    #[test]
    fn test_parse_snp_record_rejects_identical_alleles() {
        // L9: exon.rs's snp: spec syntax already rejects REF == ALT, but the
        // VCF ingest path had no such check at all, so a REF=A ALT=A record
        // silently became a no-op "variant" in the truth VCF.
        let vcf = "chr1\t100\ttest_snp\tA\tA\t.\t.\t.\n";
        let records = parse_vcf_records(vcf.as_bytes()).unwrap();
        assert!(records_to_events(records).is_err());
    }

    #[test]
    fn test_parse_snp_record_rejects_identical_alleles_different_case() {
        // L9: REF and ALT that are the same base but differ only in case
        // (e.g. soft-masked casing) must be caught too, not just an exact
        // byte-for-byte match.
        let vcf = "chr1\t100\ttest_snp\tA\ta\t.\t.\t.\n";
        let records = parse_vcf_records(vcf.as_bytes()).unwrap();
        assert!(records_to_events(records).is_err());
    }

    #[test]
    fn test_parse_small_deletion_record() {
        // Small del: REF=ACG, ALT=A at POS=100 (1-based) → 0-based pos=99
        let vcf = "chr1\t100\ttest_del\tACG\tA\t.\t.\t.\n";
        let records = parse_vcf_records(vcf.as_bytes()).unwrap();
        let events = records_to_events(records).unwrap();
        assert_eq!(events.len(), 1);
        match &events[0] {
            SimEvent::SmallVariant {
                pos,
                ref_allele,
                alt_allele,
                ..
            } => {
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
            SimEvent::SmallVariant {
                pos,
                ref_allele,
                alt_allele,
                ..
            } => {
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
            SimEvent::SmallVariant {
                allele_fraction, ..
            } => {
                assert_eq!(*allele_fraction, Some(0.25));
            }
            _ => panic!("expected SmallVariant"),
        }
    }

    /// The alleles a SimEvent carries, for comparing two entry points.
    fn small_variant_alleles(event: &SimEvent) -> (Vec<u8>, Vec<u8>) {
        match event {
            SimEvent::SmallVariant {
                ref_allele,
                alt_allele,
                ..
            } => (ref_allele.clone(), alt_allele.clone()),
            other => panic!("expected SmallVariant, got {:?}", other),
        }
    }

    /// types.rs documents SmallVariant's alleles as uppercase and
    /// exon.rs's `snp:` spec upholds it, so a soft-masked lowercase record
    /// read from a VCF must not land in the truth VCF in a different case
    /// than the same variant typed on the command line.
    #[test]
    fn test_small_variant_alleles_are_uppercased_like_the_event_spec() {
        let vcf = "chr1\t100\tsnp1\ta\tc\t.\t.\t.\n";
        let records = parse_vcf_records(vcf.as_bytes()).unwrap();
        let events = records_to_events(records).unwrap();
        let (spec_event, _) = crate::exon::parse_event_spec("snp:chr1:100:a:c", &[]).unwrap();
        assert_eq!(
            small_variant_alleles(&events[0]),
            small_variant_alleles(&spec_event),
            "the same variant through --vcf and --event must store the same alleles"
        );
        assert_eq!(
            small_variant_alleles(&events[0]),
            (b"A".to_vec(), b"C".to_vec())
        );
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
            SimEvent::Deletion {
                del_start, del_end, ..
            } => {
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
            SimEvent::Duplication {
                dup_start, dup_end, ..
            } => {
                assert_eq!(*dup_start, 500);
                assert_eq!(*dup_end, 1000);
            }
            _ => panic!("expected Duplication"),
        }

        // BND: t[p[ at POS=100 keeps base 100 of chr1, then chr2 from base 200
        let vcf = "chr1\t100\ttest_bnd\tN\tN[chr2:200[\t.\t.\tSVTYPE=BND\n";
        let records = parse_vcf_records(vcf.as_bytes()).unwrap();
        let events = records_to_events(records).unwrap();
        match &events[0] {
            SimEvent::Fusion {
                bp_a,
                bp_b,
                join,
                ..
            } => {
                assert_eq!(*bp_a, 100); // t[p[ keeps base 100: cut after it
                assert_eq!(*bp_b, 199); // ALT pos 200 → 0-based 199
                assert_eq!(*join, FusionJoin::Forward, "N[chr:pos[ should be forward");
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

    /// DEL/DUP/INV with no END and no SVLEN must not silently become a 1 bp
    /// event. When REF is sequence-resolved (anchor + affected bases, ALT is
    /// the anchor alone), the length is derived from REF; when REF is a
    /// single base there is no length information at all, so the record is
    /// rejected loudly instead of guessed at.
    #[test]
    fn test_del_no_end_no_svlen_derives_length_from_ref() {
        // REF=ACGT, ALT=A: 3 deleted bases (ACGT minus the anchor A).
        let vcf = "chr1\t100\ttest_del\tACGT\tA\t.\t.\tSVTYPE=DEL\n";
        let records = parse_vcf_records(vcf.as_bytes()).unwrap();
        let events = records_to_events(records).unwrap();
        assert_eq!(events.len(), 1);
        match &events[0] {
            SimEvent::Deletion {
                del_start, del_end, ..
            } => {
                assert_eq!(*del_start, 100);
                assert_eq!(*del_end, 103);
            }
            _ => panic!("expected Deletion"),
        }
    }

    #[test]
    fn test_dup_no_end_no_svlen_derives_length_from_ref() {
        let vcf = "chr1\t100\ttest_dup\tACGT\tA\t.\t.\tSVTYPE=DUP\n";
        let records = parse_vcf_records(vcf.as_bytes()).unwrap();
        let events = records_to_events(records).unwrap();
        assert_eq!(events.len(), 1);
        match &events[0] {
            SimEvent::Duplication {
                dup_start, dup_end, ..
            } => {
                assert_eq!(*dup_start, 100);
                assert_eq!(*dup_end, 103);
            }
            _ => panic!("expected Duplication"),
        }
    }

    #[test]
    fn test_inv_no_end_no_svlen_derives_length_from_ref() {
        let vcf = "chr1\t100\ttest_inv\tACGT\tA\t.\t.\tSVTYPE=INV\n";
        let records = parse_vcf_records(vcf.as_bytes()).unwrap();
        let events = records_to_events(records).unwrap();
        assert_eq!(events.len(), 1);
        match &events[0] {
            SimEvent::Inversion {
                inv_start, inv_end, ..
            } => {
                assert_eq!(*inv_start, 100);
                assert_eq!(*inv_end, 103);
            }
            _ => panic!("expected Inversion"),
        }
    }

    #[test]
    fn test_del_no_end_no_svlen_single_base_ref_is_rejected() {
        // Symbolic ALT, single-base REF, no END, no SVLEN: no length
        // information exists. Must not silently become a 1 bp deletion.
        let vcf = "chr1\t100\ttest_del\tN\t<DEL>\t.\t.\tSVTYPE=DEL\n";
        let records = parse_vcf_records(vcf.as_bytes()).unwrap();
        let events = records_to_events(records).unwrap();
        assert!(events.is_empty());
    }

    #[test]
    fn test_dup_no_end_no_svlen_single_base_ref_is_rejected() {
        let vcf = "chr1\t100\ttest_dup\tN\t<DUP>\t.\t.\tSVTYPE=DUP\n";
        let records = parse_vcf_records(vcf.as_bytes()).unwrap();
        let events = records_to_events(records).unwrap();
        assert!(events.is_empty());
    }

    #[test]
    fn test_inv_no_end_no_svlen_single_base_ref_is_rejected() {
        let vcf = "chr1\t100\ttest_inv\tN\t<INV>\t.\t.\tSVTYPE=INV\n";
        let records = parse_vcf_records(vcf.as_bytes()).unwrap();
        let events = records_to_events(records).unwrap();
        assert!(events.is_empty());
    }

    /// Real chr20:38412500-38412520 reference sequence, used by the
    /// sequence-resolved cases below.
    const REF21: &str = "GTTAAAGTTTATCAGAAAATT";
    /// Its reverse complement.
    const RC21: &str = "AATTTTCTGATAAACTTTAAC";

    /// A DEL whose ALT carries sequence of its own is not "REF = anchor +
    /// deleted bases": here REF and ALT share a 7 bp prefix, so the real
    /// event is a 14 bp deletion at 38412507-38412520 anchored at 38412506.
    /// Deriving the span from REF alone gave END=38412520;SVLEN=-20 — a
    /// plausible-looking record starting 6 bases too early.
    #[test]
    fn test_del_sequence_resolved_multibase_alt_strips_the_shared_prefix() {
        let vcf = format!("chr20\t38412500\t.\t{}\tGTTAAAG\t.\t.\tSVTYPE=DEL\n", REF21);
        let records = parse_vcf_records(vcf.as_bytes()).unwrap();
        let events = records_to_events(records).unwrap();
        assert_eq!(events.len(), 1, "expected one DEL event, got {:?}", events);
        match &events[0] {
            SimEvent::Deletion {
                del_start, del_end, ..
            } => {
                // The 14 deleted bases are 1-based 38412507-38412520, so the
                // preceding base is 38412506 and END is 38412520.
                assert_eq!(*del_start, 38412506);
                assert_eq!(*del_end, 38412520);
            }
            _ => panic!("expected Deletion"),
        }

        // The same span must come back out of the truth VCF's END form.
        let round = "chr20\t38412506\t.\tG\t<DEL>\t999\tPASS\tSVTYPE=DEL;END=38412520;SVLEN=-14\n";
        let records = parse_vcf_records(round.as_bytes()).unwrap();
        let events = records_to_events(records).unwrap();
        match &events[0] {
            SimEvent::Deletion {
                del_start, del_end, ..
            } => {
                assert_eq!((*del_start, *del_end), (38412506, 38412520));
            }
            _ => panic!("expected Deletion"),
        }
    }

    /// An INV written as an equal-length substitution (REF = the inverted
    /// span, ALT = its reverse complement). VCF 4.3 uses no padding base
    /// when neither allele is empty, so POS is the first inverted base, not
    /// the preceding one: the span is the 21 bases 38412500-38412520.
    /// Deriving from REF as if POS were a padding base gave END=38412520
    /// over a 20 bp span starting at 38412501 — one base short and one base
    /// right of the truth.
    #[test]
    fn test_inv_equal_length_ref_and_alt_spans_ref_from_pos() {
        let vcf = format!(
            "chr20\t38412500\t.\t{}\t{}\t.\t.\tSVTYPE=INV\n",
            REF21, RC21
        );
        let records = parse_vcf_records(vcf.as_bytes()).unwrap();
        let events = records_to_events(records).unwrap();
        assert_eq!(events.len(), 1, "expected one INV event, got {:?}", events);
        match &events[0] {
            SimEvent::Inversion {
                inv_start, inv_end, ..
            } => {
                assert_eq!(*inv_start, 38412499);
                assert_eq!(*inv_end, 38412520);
                assert_eq!(*inv_end - *inv_start, REF21.len() as u64);
            }
            _ => panic!("expected Inversion"),
        }
    }

    /// An equal-length REF/ALT pair whose ALT is not REF reverse-complemented
    /// does not spell an inversion at all. Reading its length off REF anyway
    /// is exactly the plausible-looking-but-wrong record this arm must not
    /// emit, so it is rejected.
    #[test]
    fn test_inv_equal_length_alt_that_is_not_the_reverse_complement_is_rejected() {
        let not_rc = "ACGTACGTACGTACGTACGTA";
        assert_eq!(not_rc.len(), REF21.len());
        let vcf = format!(
            "chr20\t38412500\t.\t{}\t{}\t.\t.\tSVTYPE=INV\n",
            REF21, not_rc
        );
        let records = parse_vcf_records(vcf.as_bytes()).unwrap();
        let events = records_to_events(records).unwrap();
        assert!(
            events.is_empty(),
            "an equal-length pair that is not an inversion must not decode: {:?}",
            events
        );
    }

    /// The well-formed sequence-resolved DUP: single-base REF anchor, ALT =
    /// anchor + the duplicated copy. The length is in ALT, exactly as for a
    /// sequence INS, so it must be read from there rather than rejected.
    #[test]
    fn test_dup_no_end_no_svlen_derives_length_from_alt_sequence() {
        let vcf = format!("chr20\t38412499\t.\tT\tT{}\t.\t.\tSVTYPE=DUP\n", REF21);
        let records = parse_vcf_records(vcf.as_bytes()).unwrap();
        let events = records_to_events(records).unwrap();
        assert_eq!(events.len(), 1, "expected one DUP event, got {:?}", events);
        match &events[0] {
            SimEvent::Duplication {
                dup_start, dup_end, ..
            } => {
                // POS is the preceding base; the 21 duplicated bases are
                // 1-based 38412500-38412520 == 0-based [38412499, 38412520).
                assert_eq!(*dup_start, 38412499);
                assert_eq!(*dup_end, 38412520);
            }
            _ => panic!("expected Duplication"),
        }
    }

    /// A DUP with sequence in both REF and ALT is the same ambiguous shape
    /// as the DEL above (shared prefix), not the anchor+copy shape.
    #[test]
    fn test_dup_multibase_ref_and_alt_is_rejected() {
        let vcf = "chr20\t38412499\t.\tTGTT\tTGTTTGTT\t.\t.\tSVTYPE=DUP\n";
        let records = parse_vcf_records(vcf.as_bytes()).unwrap();
        let events = records_to_events(records).unwrap();
        assert!(
            events.is_empty(),
            "REF and ALT both carrying sequence must not be read as a length: {:?}",
            events
        );
    }

    /// A derived end must survive the round trip through the truth VCF's
    /// `SVTYPE=...;END=` form: END is 1-based inclusive and numerically
    /// equal to the 0-based half-open end, POS to the 0-based start.
    #[test]
    fn test_derived_sv_end_round_trips_through_end_info() {
        let vcf = format!("chr20\t38412499\t.\tT\tT{}\t.\t.\tSVTYPE=DUP\n", REF21);
        let records = parse_vcf_records(vcf.as_bytes()).unwrap();
        let events = records_to_events(records).unwrap();
        assert_eq!(events.len(), 1, "expected one DUP event, got {:?}", events);
        let (start, end) = match &events[0] {
            SimEvent::Duplication {
                dup_start, dup_end, ..
            } => (*dup_start, *dup_end),
            _ => panic!("expected Duplication"),
        };

        // The shape truth.rs writes for a Duplication.
        let round = format!(
            "chr20\t{}\t.\tT\t<DUP>\t999\tPASS\tSVTYPE=DUP;END={};SVLEN={}\n",
            start,
            end,
            end - start
        );
        let records = parse_vcf_records(round.as_bytes()).unwrap();
        let events = records_to_events(records).unwrap();
        match &events[0] {
            SimEvent::Duplication {
                dup_start, dup_end, ..
            } => {
                assert_eq!(*dup_start, start);
                assert_eq!(*dup_end, end);
            }
            _ => panic!("expected Duplication"),
        }
    }

    /// The inserted bases are what ALT adds, not everything past ALT's first
    /// base: REF=AT ALT=ATGGG inserts the 3 bases GGG after the T, and taking
    /// ALT[1..] made it a 4 bp insertion of TGGG.
    #[test]
    fn test_sequence_ins_strips_the_prefix_ref_and_alt_share() {
        let vcf = "chr1\t100\ttest_ins\tAT\tATGGG\t.\t.\tSVTYPE=INS\n";
        let records = parse_vcf_records(vcf.as_bytes()).unwrap();
        let events = records_to_events(records).unwrap();
        assert_eq!(events.len(), 1, "expected one INS event, got {:?}", events);
        match &events[0] {
            SimEvent::Insertion {
                pos,
                ins_seq,
                ins_len,
                ..
            } => {
                assert_eq!(*ins_len, 3);
                assert_eq!(ins_seq.as_deref(), Some(b"GGG".as_slice()));
                // The last base both alleles keep is 1-based 101 (the T), so
                // the novel bases go in after it.
                assert_eq!(*pos, 101);
            }
            _ => panic!("expected Insertion"),
        }
    }

    /// REF and ALT may share a suffix as well (REF=AT ALT=AGGGT is the same
    /// 3 bp insertion, written unnormalised). The anchor is then the shared
    /// prefix's last base, exactly as for a left-aligned record.
    #[test]
    fn test_sequence_ins_strips_a_shared_suffix_too() {
        let vcf = "chr1\t100\ttest_ins\tAT\tAGGGT\t.\t.\tSVTYPE=INS\n";
        let records = parse_vcf_records(vcf.as_bytes()).unwrap();
        let events = records_to_events(records).unwrap();
        assert_eq!(events.len(), 1, "expected one INS event, got {:?}", events);
        match &events[0] {
            SimEvent::Insertion {
                pos,
                ins_seq,
                ins_len,
                ..
            } => {
                assert_eq!(*ins_len, 3);
                assert_eq!(ins_seq.as_deref(), Some(b"GGG".as_slice()));
                assert_eq!(*pos, 100);
            }
            _ => panic!("expected Insertion"),
        }
    }

    /// The ordinary padded form — REF the anchor base alone — is unchanged.
    #[test]
    fn test_sequence_ins_with_single_base_ref_is_unchanged() {
        let vcf = "chr1\t100\ttest_ins\tA\tAGGG\t.\t.\tSVTYPE=INS\n";
        let records = parse_vcf_records(vcf.as_bytes()).unwrap();
        let events = records_to_events(records).unwrap();
        match &events[0] {
            SimEvent::Insertion {
                pos,
                ins_seq,
                ins_len,
                ..
            } => {
                assert_eq!(*ins_len, 3);
                assert_eq!(ins_seq.as_deref(), Some(b"GGG".as_slice()));
                assert_eq!(*pos, 100);
            }
            _ => panic!("expected Insertion"),
        }
    }

    /// REF bases that ALT drops make the record a complex indel, not an
    /// insertion: which bases are inserted and which replaced is a guess.
    /// Drop it loudly rather than invent an inserted sequence.
    #[test]
    fn test_sequence_ins_with_ref_bases_alt_drops_is_rejected() {
        let vcf = "chr1\t100\ttest_ins\tATT\tATGGG\t.\t.\tSVTYPE=INS\n";
        let records = parse_vcf_records(vcf.as_bytes()).unwrap();
        let events = records_to_events(records).unwrap();
        assert!(
            events.is_empty(),
            "a complex REF/ALT pair must not decode as an insertion: {:?}",
            events
        );
    }

    /// Like the DEL/DUP/INV skip warning, the INS one must locate the record
    /// it dropped — `ID=.` is common — so L11 can count what was lost.
    #[test]
    fn test_not_an_insertion_warning_identifies_record_by_chrom_and_pos() {
        let vcf = "chr20\t38412500\t.\tATT\tATGGG\t.\t.\tSVTYPE=INS\n";
        let records = parse_vcf_records(vcf.as_bytes()).unwrap();
        let msg = not_an_insertion_warning(&records[0]);
        assert!(
            msg.contains("chr20:38412500"),
            "warning must locate the record: {}",
            msg
        );
        assert!(msg.contains("INS"), "warning must name the type: {}", msg);
        assert!(
            msg.contains("skipping"),
            "warning must say the record is dropped: {}",
            msg
        );
    }

    /// SV records routinely carry `ID=.`, so the skip warning must name
    /// chrom:pos as well — otherwise the dropped record cannot be found.
    #[test]
    fn test_no_length_warning_identifies_record_by_chrom_and_pos() {
        let vcf = "chr20\t38412500\t.\tN\t<DEL>\t.\t.\tSVTYPE=DEL\n";
        let records = parse_vcf_records(vcf.as_bytes()).unwrap();
        let msg = no_length_warning(&records[0], "DEL");
        assert!(
            msg.contains("chr20:38412500"),
            "warning must locate the record: {}",
            msg
        );
        assert!(msg.contains("DEL"), "warning must name the type: {}", msg);
        assert!(
            msg.contains("skipping"),
            "warning must say the record is dropped: {}",
            msg
        );
    }

    /// A single-base ALT is not always REF's anchor base: `REF=ACGT ALT=T`
    /// shares its *suffix*, so the deleted bases are `ACG` at 1-based
    /// 100-102 and the preceding base is 99. Taking REF's length from POS
    /// claimed 100-103 — a record that deletes `CGT` and leaves `A`, where
    /// the alleles say the result is `T`.
    #[test]
    fn test_del_single_base_alt_that_is_not_the_anchor_strips_the_shared_suffix() {
        let vcf = "chr1\t100\ttest_del\tACGT\tT\t.\t.\tSVTYPE=DEL\n";
        let records = parse_vcf_records(vcf.as_bytes()).unwrap();
        let events = records_to_events(records).unwrap();
        assert_eq!(events.len(), 1, "expected one DEL event, got {:?}", events);
        match &events[0] {
            SimEvent::Deletion {
                del_start, del_end, ..
            } => {
                assert_eq!(*del_start, 99);
                assert_eq!(*del_end, 102);
            }
            _ => panic!("expected Deletion"),
        }

        // The normalised spelling of the same event still decodes the way it
        // always did: REF=ACGT ALT=A is the 3 bases after the anchor.
        let normalised = "chr1\t100\ttest_del\tACGT\tA\t.\t.\tSVTYPE=DEL\n";
        let records = parse_vcf_records(normalised.as_bytes()).unwrap();
        let events = records_to_events(records).unwrap();
        match &events[0] {
            SimEvent::Deletion {
                del_start, del_end, ..
            } => {
                assert_eq!((*del_start, *del_end), (100, 103));
            }
            _ => panic!("expected Deletion"),
        }
    }

    /// Stripping can move the start left of POS, and at POS=1 that is base
    /// 0 — a position no VCF record can name (truth.rs would write `POS=0`
    /// with `REF=N`). Reject those loudly instead.
    #[test]
    fn test_alleles_that_place_the_event_before_the_first_base_are_rejected() {
        // DEL: shared suffix only, so the preceding base would be 0.
        let del = "chr1\t1\ttest_del\tACGT\tT\t.\t.\tSVTYPE=DEL\n";
        let records = parse_vcf_records(del.as_bytes()).unwrap();
        let events = records_to_events(records).unwrap();
        assert!(events.is_empty(), "DEL before base 1 must not decode: {:?}", events);

        // INV written as an equal-length substitution: POS is the first
        // inverted base, so the preceding base would be 0.
        let inv = "chr1\t1\ttest_inv\tAGTT\tAACT\t.\t.\tSVTYPE=INV\n";
        let records = parse_vcf_records(inv.as_bytes()).unwrap();
        let events = records_to_events(records).unwrap();
        assert!(events.is_empty(), "INV before base 1 must not decode: {:?}", events);

        // INS whose alleles share only a suffix, same story.
        let ins = "chr1\t1\ttest_ins\tAT\tGGGAT\t.\t.\tSVTYPE=INS\n";
        let records = parse_vcf_records(ins.as_bytes()).unwrap();
        let events = records_to_events(records).unwrap();
        assert!(events.is_empty(), "INS before base 1 must not decode: {:?}", events);
    }

    /// The before-the-first-base rejection is a drop like any other, so it
    /// must name chrom:pos for L11 to count it.
    #[test]
    fn test_before_first_base_warning_identifies_record_by_chrom_and_pos() {
        let vcf = "chr20\t1\t.\tACGT\tT\t.\t.\tSVTYPE=DEL\n";
        let records = parse_vcf_records(vcf.as_bytes()).unwrap();
        let msg = before_first_base_warning(&records[0], "DEL");
        assert!(
            msg.contains("chr20:1"),
            "warning must locate the record: {}",
            msg
        );
        assert!(msg.contains("DEL"), "warning must name the type: {}", msg);
        assert!(
            msg.contains("skipping"),
            "warning must say the record is dropped: {}",
            msg
        );
    }

    /// An inversion's span is stated by the record, not derived from where
    /// its alleles happen to differ: `AGTT` and its reverse complement
    /// `AACT` begin and end with complementary bases, so trimming the flanks
    /// they share left a 2 bp inversion of the middle where the record spells
    /// out 4 bp. About one equal-length INV in four has such an end.
    #[test]
    fn test_inv_equal_length_self_complementary_ends_span_the_whole_record() {
        let vcf = "chr1\t100\ttest_inv\tAGTT\tAACT\t.\t.\tSVTYPE=INV\n";
        let records = parse_vcf_records(vcf.as_bytes()).unwrap();
        let events = records_to_events(records).unwrap();
        assert_eq!(events.len(), 1, "expected one INV event, got {:?}", events);
        match &events[0] {
            SimEvent::Inversion {
                inv_start, inv_end, ..
            } => {
                // No padding base, so POS is the first inverted base: the
                // 4 bases 1-based 100-103.
                assert_eq!(*inv_start, 99);
                assert_eq!(*inv_end, 103);
            }
            _ => panic!("expected Inversion"),
        }
    }

    /// The other spelling: a padding base, then the inverted region. The
    /// whole-allele pair is not a reverse-complement pair here, so it has to
    /// be checked past the anchor — but past the anchor only, never past the
    /// rest of the flanks the alleles share, which narrowed this record to
    /// the 2 bp its self-complementary ends leave over.
    #[test]
    fn test_inv_equal_length_with_a_padding_base_inverts_ref_past_the_anchor() {
        let vcf = "chr1\t100\ttest_inv\tTAGTT\tTAACT\t.\t.\tSVTYPE=INV\n";
        let records = parse_vcf_records(vcf.as_bytes()).unwrap();
        let events = records_to_events(records).unwrap();
        assert_eq!(events.len(), 1, "expected one INV event, got {:?}", events);
        match &events[0] {
            SimEvent::Inversion {
                inv_start, inv_end, ..
            } => {
                // POS is the anchor T; the 4 inverted bases are 101-104.
                assert_eq!(*inv_start, 100);
                assert_eq!(*inv_end, 104);
            }
            _ => panic!("expected Inversion"),
        }
    }

    /// A padded VCF ALT begins with the REF anchor. `REF=T ALT=GGGT` does
    /// not, so it does not spell "anchor + duplicated copy" and there is no
    /// saying where the extra bases come from. L7 took `ALT[1..]` (`GGT`);
    /// it is now rejected with the usual warning.
    #[test]
    fn test_dup_single_base_ref_that_is_not_alts_first_base_is_rejected() {
        let vcf = "chr1\t100\ttest_dup\tT\tGGGT\t.\t.\tSVTYPE=DUP\n";
        let records = parse_vcf_records(vcf.as_bytes()).unwrap();
        let events = records_to_events(records).unwrap();
        assert!(
            events.is_empty(),
            "an ALT that does not start with the REF anchor must not decode: {:?}",
            events
        );
    }

    /// The rejected DUP shape is the one whose ALT is *longer*
    /// (`TGTT`/`TGTTTGTT`), where the copy's source is ambiguous. An ALT
    /// that only drops REF bases states a span like any other, and is read
    /// as one for all three span types.
    #[test]
    fn test_dup_alt_shorter_than_ref_spans_the_bases_ref_keeps() {
        let vcf = "chr1\t100\ttest_dup\tTGTT\tTG\t.\t.\tSVTYPE=DUP\n";
        let records = parse_vcf_records(vcf.as_bytes()).unwrap();
        let events = records_to_events(records).unwrap();
        assert_eq!(events.len(), 1, "expected one DUP event, got {:?}", events);
        match &events[0] {
            SimEvent::Duplication {
                dup_start, dup_end, ..
            } => {
                assert_eq!((*dup_start, *dup_end), (101, 103));
            }
            _ => panic!("expected Duplication"),
        }
    }
}
