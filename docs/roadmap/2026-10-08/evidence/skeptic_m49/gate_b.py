"""Gate B of M49 on one existing DEL run: true origin (20-mer placement on the
DEL haplotype, as scripts/low_af/locate.py does) vs the aligner's placement in
sim.bam, primary records only.

Usage: gate_b.py SIM_BAM REFERENCE CHROM START0 END0 [FLANK]
START0/END0: the deleted interval, 0-based half-open as spike's haplotype joins
[START0-FLANK, START0) + [END0, END0+FLANK).
"""
import sys
from collections import Counter

import pysam

SEED = 20


def rc(s):
    return s.translate(str.maketrans('ACGTN', 'TGCAN'))[::-1]


def main(bam, ref, chrom, start, end, flank=2000):
    start, end, flank = int(start), int(end), int(flank)
    fa = pysam.FastaFile(ref)
    left = fa.fetch(chrom, start - flank, start).upper()
    right = fa.fetch(chrom, end, end + flank).upper()
    hap = left + right
    L = len(left)
    idx = {}
    for i in range(len(hap) - SEED + 1):
        idx.setdefault(hap[i:i + SEED], []).append(i)

    def to_ref(p):
        return start - flank + p if p < L else end + (p - L)

    def place(seq):
        votes = Counter()
        for strand, s in (('+', seq), ('-', rc(seq))):
            for o in range(0, len(s) - SEED + 1, SEED):
                hits = idx.get(s[o:o + SEED], [])
                if len(hits) == 1:
                    votes[(strand, hits[0] - o)] += 1
        return votes.most_common(1)[0][0] if votes else None

    rows = Counter()
    far = []
    for r in pysam.AlignmentFile(bam):
        if r.is_secondary or r.is_supplementary or not r.query_name.startswith('SPIKE_'):
            continue
        seq = r.query_sequence  # reference-forward orientation as aligned
        spot = place(seq)
        if spot is None:
            rows['unplaced_by_locator'] += 1
            continue
        strand, s = spot
        crosses = s < L < s + len(seq)
        kind = 'junction' if crosses else 'flank'
        if r.is_unmapped:
            rows[(kind, 'unmapped')] += 1
            continue
        if strand != '+':
            rows[(kind, 'wrong_strand')] += 1
            far.append((r.query_name, kind, 'strand', r.reference_name, r.reference_start, r.cigarstring))
            continue
        qs = r.query_alignment_start
        p = s + qs
        if not 0 <= p < len(hap):
            rows[(kind, 'off_hap')] += 1
            continue
        true_pos = to_ref(p)
        d = abs(r.reference_start - true_pos) if r.reference_name == chrom else 10**9
        ok = d <= 10
        clipped = any(op == 4 for op, _ in r.cigartuples)
        rows[(kind, 'within10' if ok else 'off', 'clipped' if clipped else 'unclipped',
              'mapq<20' if r.mapping_quality < 20 else 'mapq>=20')] += 1
        if not ok:
            far.append((r.query_name, kind, d, r.reference_name, r.reference_start, true_pos,
                        r.cigarstring, r.mapping_quality))
    for k, v in sorted(rows.items(), key=lambda kv: str(kv[0])):
        print(k, v)
    print('misplaced (>10 bp or wrong strand):', len(far))
    for f in far[:40]:
        print('  ', f)


if __name__ == '__main__':
    main(*sys.argv[1:])
