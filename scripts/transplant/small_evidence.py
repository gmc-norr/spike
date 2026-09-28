#!/usr/bin/env python3
"""Read evidence of SNVs, 1-49 bp indels and duplications, measured the same way for real and spiked ones.

Round 2 of the transplant test, locked in
docs/superpowers/plans/2026-09-28-transplant-round2.md. Reads as evidence.py
(round 1) reads them: primary, mapped, not duplicate, not QC-fail, any MAPQ;
a spiked BAM is the recipient BAM less the names spike replaced, plus sim.bam.

An event is (kind, chrom, POS, REF, ALT) as the truth VCF writes it (POS
1-based, REF and ALT with their anchor base); a duplication is the insertion
whose bases copy the reference next to it.

  SNV   n_carry  reads with ALT at POS         n_ref  reads with another base there
        A        n_carry / (n_carry + n_ref)
  indel repeat region [rs, re): the deleted or inserted bases as a unit,
        extended both ways while the reference repeats it
        n_carry  reads with an indel of the same type and length starting in
                 [rs - 1, re + 1] (an insertion: spelling the same sequence)
        n_ref    reads aligned over [rs - 10, re + 10) with no indel or clip there
        A        n_carry / (n_carry + n_ref)
        n_any    read pairs with a carrier, a same-type indel >= half the length
                 starting in [rs - 10, re + 10], or a soft clip >= 5 bp there
        E        n_any / flank depth (the 1 kb each side of [rs, re))
  DUP   segment [s, e): the copy the insertion repeats
        n_any    read pairs with an insertion >= half the segment within 20 bp
                 of it, a soft clip >= 10 bp within 10 bp of s or e, or a split
                 alignment ending near e with its other part starting near s
        J        n_any / flank depth (the 1 kb each side of [s, e))
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import evidence  # noqa: E402

FLANK = evidence.FLANK
CARRY_SLACK = 1
SPAN_SLACK = 10
MIN_SMALL_CLIP = 5
DUP_INS_SLACK = 20
COLUMNS = ["n_carry", "n_ref", "n_any", "flank_depth", "A", "E", "J"]
M_OPS = (0, 7, 8)  # M, =, X


def repeat_region(fetch, chrom, start, end, unit):
    """[rs, re): `unit`, laid from `start`, extended both ways while the reference repeats it."""
    n = len(unit)
    rs, re_ = start, end
    while True:
        chunk = fetch(chrom, re_, re_ + 1_000)
        i = 0
        while i < len(chunk) and chunk[i] == unit[(re_ + i - start) % n]:
            i += 1
        re_ += i
        if i < len(chunk) or not chunk:
            break
    while rs > 0:
        lo = max(0, rs - 1_000)
        chunk = fetch(chrom, lo, rs)
        i = len(chunk)
        while i > 0 and chunk[i - 1] == unit[(lo + i - 1 - start) % n]:
            i -= 1
        stop = i > 0 or lo == 0
        rs = lo + i
        if stop:
            break
    return rs, re_


def segment_of(fetch, chrom, pos, ref, alt):
    """The copy an insertion repeats, [s, e) 0-based, or None when it copies neither side."""
    ins = alt[len(ref):].upper()
    n = len(ins)
    if fetch(chrom, pos, pos + n) == ins:
        return (pos, pos + n)
    if fetch(chrom, pos - n, pos) == ins:
        return (pos - n, pos)
    return None


def cigar_ops(r):
    """(op, length, reference position, query position) at the start of every CIGAR op."""
    rpos, qpos = r.reference_start, 0
    for op, n in r.cigartuples:
        yield op, n, rpos, qpos
        if op in (0, 2, 3, 7, 8):
            rpos += n
        if op in (0, 1, 4, 7, 8):
            qpos += n


def soft_clips(r):
    """(clip point, length) of the soft clip at either end, if any."""
    ops = [c for c in r.cigartuples if c[0] != 5]
    out = []
    if ops and ops[0][0] == 4:
        out.append((r.reference_start, ops[0][1]))
    if len(ops) > 1 and ops[-1][0] == 4:
        out.append((r.reference_end, ops[-1][1]))
    return out


def flank_depth(blocks, lo, hi):
    return (evidence.depth(blocks, lo - FLANK, lo) + evidence.depth(blocks, hi, hi + FLANK)) / 2


def measure_snv(paths, chrom, pos, alt):
    p0 = pos - 1
    carry = other = 0
    for r in evidence.reads(paths, chrom, p0, p0 + 1):
        for op, n, rp, qp in cigar_ops(r):
            if rp <= p0 < rp + n and op in (0, 2, 3, 7, 8):
                if op in M_OPS:
                    if r.query_sequence[qp + p0 - rp].upper() == alt.upper():
                        carry += 1
                    else:
                        other += 1
                break
    return {"n_carry": carry, "n_ref": other, "A": carry / (carry + other) if carry + other else None}


def same_insertion(fetch, chrom, rs, re_, q, truth, p, seen):
    """Whether inserting `seen` at p spells what inserting `truth` at q does."""
    w0, w1 = min(rs, p, q) - 5, max(re_, p, q) + 5
    ref = fetch(chrom, w0, w1)
    return ref[:q - w0] + truth + ref[q - w0:] == ref[:p - w0] + seen.upper() + ref[p - w0:]


def measure_indel(paths, fetch, kind, chrom, pos, ref, alt):
    is_del = kind == "DEL"
    unit = (ref[1:] if is_del else alt[1:]).upper()
    n_len = len(unit)
    start, end = pos, pos + (n_len if is_del else 0)
    rs, re_ = repeat_region(fetch, chrom, start, end, unit)
    lo, hi = rs - SPAN_SLACK, re_ + SPAN_SLACK
    carry = refs = 0
    any_, blocks = set(), []
    for r in evidence.reads(paths, chrom, rs - FLANK - 200, re_ + FLANK + 200):
        blocks.append(r.get_blocks())
        carrier = indel_there = shows = False
        for op, n, rp, qp in cigar_ops(r):
            if op not in (1, 2):
                continue
            if (op == 1 and lo <= rp <= hi) or (op == 2 and rp <= hi and rp + n >= lo):
                indel_there = True
            if (op == 2) != is_del or not lo <= rp <= hi:
                continue
            if n >= 0.5 * n_len:
                shows = True
            if n == n_len and rs - CARRY_SLACK <= rp <= re_ + CARRY_SLACK and (
                    is_del or same_insertion(fetch, chrom, rs, re_, start, unit, rp, r.query_sequence[qp:qp + n])):
                carrier = True
        clips = [(p, n) for p, n in soft_clips(r) if lo <= p <= hi]
        if any(n >= MIN_SMALL_CLIP for _, n in clips):
            shows = True
        if carrier:
            carry += 1
            shows = True
        elif r.reference_start <= lo and r.reference_end >= hi and not indel_there and not clips:
            refs += 1
        if shows:
            any_.add(r.query_name)
    flank = flank_depth(blocks, rs, re_)
    return {"n_carry": carry, "n_ref": refs, "n_any": len(any_), "flank_depth": flank,
            "A": carry / (carry + refs) if carry + refs else None,
            "E": len(any_) / flank if flank > 0 else None}


def measure_dup(paths, fetch, chrom, pos, ref, alt):
    s, e = segment_of(fetch, chrom, pos, ref, alt)
    any_, blocks = set(), []
    for r in evidence.reads(paths, chrom, s - FLANK - 200, e + FLANK + 200):
        blocks.append(r.get_blocks())
        shows = any(op == 1 and n >= 0.5 * (e - s) and s - DUP_INS_SLACK <= rp <= e + DUP_INS_SLACK
                    for op, n, rp, _ in cigar_ops(r))
        # A tandem copy reads ...e][s...: the part before the join ends at e, the part after starts at s.
        if shows or evidence.is_clip(r, s, e) or evidence.is_split(r, chrom, e, s):
            any_.add(r.query_name)
    flank = flank_depth(blocks, s, e)
    return {"n_any": len(any_), "flank_depth": flank, "J": len(any_) / flank if flank > 0 else None}


def measure(bam, fetch, event, skip_names=frozenset(), extra_bams=()):
    """The evidence of `event` in `bam` (less `skip_names`) plus `extra_bams`; see the module doc."""
    kind, chrom, pos, ref, alt = event
    paths = [(bam, frozenset(skip_names))] + [(p, frozenset()) for p in extra_bams]
    row = dict.fromkeys(COLUMNS)
    if kind == "SNV":
        row.update(measure_snv(paths, chrom, pos, alt))
    elif kind in ("DEL", "INS"):
        row.update(measure_indel(paths, fetch, kind, chrom, pos, ref, alt))
    elif kind == "DUP":
        row.update(measure_dup(paths, fetch, chrom, pos, ref, alt))
    else:
        raise ValueError(f"unknown kind {kind!r}")
    return row


def recipient_count(row, kind):
    """What N1 counts in the recipient before spiking: carrying reads, or a duplication's pairs."""
    return row["n_any"] if kind == "DUP" else row["n_carry"]
