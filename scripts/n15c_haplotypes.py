#!/usr/bin/env python3
"""N15, third attempt: the haplotype rule with the truth set's nearby records.

Re-implements the rule docs/review/REVIEW.md's "N15 plan, third attempt" locks, for the
het indel sites in --sites, reading the other records from --truth (the file
`spike validate` is given):

- The window is N18's repeat region, the base on each side and 10 more.
- Every other truth record whose REF span lies inside the window is a
  candidate edit, one per ALT allele. Above 10 candidates the site uses none.
- The candidate haplotypes are the reference over the window with every
  subset of {the site} + candidates applied whose REF spans do not overlap.
- A read's bases run from the base aligned to the window's first position to
  the one aligned to its last, inserted bases included; a read with no
  aligned base at either end gives no vote.
- A read carries the site if every haplotype nearest its bases (Levenshtein)
  includes the site, spans it if none does, and gives no vote otherwise.
- One vote per fragment, mates that disagree give none (N16); MAPQ >= 20.

--no-neighbours drops the candidates, which is N15's second-attempt rule.

Writes one line per site: chrom:pos, carries, spans, candidate edits.

Usage:
    n15c_haplotypes.py --bam B --ref R --sites S.vcf --truth T.vcf --chrom C --out O.tsv
"""
import argparse
import itertools
import sys
from multiprocessing import Pool
from pathlib import Path

import pysam

sys.path.insert(0, str(Path(__file__).resolve().parent))
import n18_indel_flank as n18  # noqa: E402

FLANK = 10
MAX_EDITS = 10


def lev(a, b):
    if len(a) < len(b):
        a, b = b, a
    row = list(range(len(b) + 1))
    for i, x in enumerate(a):
        prev, row[0] = row[0], i + 1
        for j, y in enumerate(b):
            cur = row[j + 1]
            row[j + 1] = min(cur + 1, row[j] + 1, prev + (x != y))
            prev = cur
    return row[-1]


def bases_over(r, first, last):
    ref_pos, q = r.reference_start, 0
    qf = ql = None
    for op, n in r.cigartuples:
        if op in (0, 7, 8):
            if ref_pos <= first < ref_pos + n:
                qf = q + first - ref_pos
            if ref_pos <= last < ref_pos + n:
                ql = q + last - ref_pos
            ref_pos += n
            q += n
        elif op in (1, 4):
            q += n
        elif op in (2, 3):
            ref_pos += n
    if qf is None or ql is None:
        return None
    return r.query_sequence[qf:ql + 1].upper()


def apply(refseq, first, edits):
    out, at = [], first
    for pos, ref, alt in sorted(edits):
        out.append(refseq[at - first:pos - first])
        out.append(alt)
        at = pos + len(ref)
    out.append(refseq[at - first:])
    return "".join(out)


def overlaps(edits):
    spans = sorted((p, p + len(r)) for p, r, _ in edits)
    return any(b0 < a1 for (_, a1), (b0, _) in zip(spans, spans[1:]))


def work(args):
    chunk, bam_path, ref_path, truth, use_nb = args
    bam, fasta = pysam.AlignmentFile(bam_path), pysam.FastaFile(ref_path)
    out = []
    for chrom, vcf_pos, ref, alt in chunk:
        pos = vcf_pos - 1
        s, e = n18.repeat_region(fasta, chrom, pos, ref, alt)
        first, last = s - 1 - FLANK, e + FLANK
        refseq = fasta.fetch(chrom, first, last + 1).upper()
        target = (pos, ref, alt)
        cands = [x for x in truth.get(chrom, []) if use_nb and x != target and first <= x[0] and x[0] + len(x[1]) - 1 <= last]
        n_cands = len(cands)
        if n_cands > MAX_EDITS:
            cands = []
        haps = []
        for k in range(len(cands) + 1):
            for sub in itertools.combinations(cands, k):
                for with_t in (False, True):
                    edits = list(sub) + ([target] if with_t else [])
                    if not overlaps(edits):
                        haps.append((with_t, apply(refseq, first, edits)))
        cache, votes = {}, []
        for rd in bam.fetch(chrom, first, last + 1):
            if not n18.usable(rd):
                continue
            b = bases_over(rd, first, last)
            if b is None:
                continue
            if b not in cache:
                d = [(lev(b, h), w) for w, h in haps]
                m = min(x for x, _ in d)
                near = {w for x, w in d if x == m}
                cache[b] = "C" if near == {True} else "S" if near == {False} else None
            if cache[b] is not None:
                votes.append((rd.query_name, cache[b]))
        fv = n18.fragment_votes(votes)
        out.append((chrom, vcf_pos, fv.count("C"), fv.count("S"), n_cands))
    return out


def load_truth(vcf, chrom):
    """Every record on chrom in the truth file, one edit per ALT: (pos0, REF, ALT)."""
    truth = {}
    for line in open(vcf):
        if line.startswith("#"):
            continue
        c = line.split("\t")
        if c[0] == chrom:
            for a in c[4].split(","):
                truth.setdefault(chrom, []).append((int(c[1]) - 1, c[3].upper(), a.upper()))
    return truth


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--bam", required=True)
    ap.add_argument("--ref", required=True)
    ap.add_argument("--sites", required=True)
    ap.add_argument("--truth", required=True)
    ap.add_argument("--chrom", required=True)
    ap.add_argument("--no-neighbours", action="store_true")
    ap.add_argument("--out", required=True)
    ap.add_argument("--procs", type=int, default=32)
    a = ap.parse_args()
    sites = [s for s in n18.load_sites(a.sites) if s[0] == a.chrom]
    truth = {} if a.no_neighbours else load_truth(a.truth, a.chrom)
    chunks = [sites[i::a.procs * 2] for i in range(a.procs * 2)]
    with Pool(a.procs) as pool:
        res = [x for part in pool.map(work, [(c, a.bam, a.ref, truth, not a.no_neighbours) for c in chunks]) for x in part]
    with open(a.out, "w") as f:
        for chrom, vcf_pos, c, sp, k in sorted(res):
            f.write(f"{chrom}:{vcf_pos}\t{c}\t{sp}\t{k}\n")


if __name__ == "__main__":
    main()
