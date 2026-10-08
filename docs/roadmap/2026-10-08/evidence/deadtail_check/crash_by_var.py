"""K2b crash share (10+ of last 20 quals < Q15, sequencing order), real vs spike, split by whether a Q100
record (any; and indels only) lies within the read's reference span +-10 bp. Primary, MAPQ>=20, all CIGARs.
Usage: python3 -I crash_by_var.py BLOCKS VCF NAME=BAM ...
"""
import sys, bisect, collections, math
import numpy as np
import pysam

blocks_f, vcf_f = sys.argv[1:3]
sets = [a.split('=', 1) for a in sys.argv[3:]]
vcf = pysam.VariantFile(vcf_f)
blocks = []
for line in open(blocks_f):
    c, s, e = line.split(); s, e = int(s) - 2000, int(e) + 2000
    allp, indp = [], []
    for rec in vcf.fetch(c, s, e):
        allp.append(rec.pos - 1)
        if rec.alts and any(len(a) != len(rec.ref) for a in rec.alts):
            indp.append(rec.pos - 1)
    blocks.append((c, s, e, sorted(allp), sorted(indp)))


def fb(c, p):
    for b in blocks:
        if b[0] == c and b[1] + 200 <= p < b[2] - 400:
            return b


def hit(lst, a, z):
    k = bisect.bisect_left(lst, a)
    return k < len(lst) and lst[k] <= z


res = {}
for name, f in sets:
    t = collections.Counter()
    for r in pysam.AlignmentFile(f):
        if r.is_unmapped or r.is_secondary or r.is_supplementary or r.mapping_quality < 20:
            continue
        b = fb(r.reference_name, r.reference_start)
        if b is None:
            continue
        q = np.asarray(r.query_qualities)
        qs = q[::-1] if r.is_reverse else q
        cr = int((qs[-20:] < 15).sum() >= 10)
        a, z = r.reference_start - 10, r.reference_end + 10
        st = 'indel-under' if hit(b[4], a, z) else ('snv-under' if hit(b[3], a, z) else 'no-Q100')
        t[(st, 0)] += 1; t[(st, 1)] += cr
    res[name] = t


def z2(x1, n1, x2, n2):
    p = (x1 + x2) / (n1 + n2)
    se = math.sqrt(p * (1 - p) * (1 / n1 + 1 / n2)) or 1
    return (x2 / n2 - x1 / n1) / se


base = sets[0][0]
for st in ('no-Q100', 'snv-under', 'indel-under'):
    line = '%-12s' % st
    for name, _ in sets:
        t = res[name]
        line += '  %s %6d reads crash %.2f%%' % (name, t[(st, 0)], 100 * t[(st, 1)] / max(t[(st, 0)], 1))
        if name != base:
            line += ' (z %.2f)' % z2(res[base][(st, 1)], res[base][(st, 0)], t[(st, 1)], t[(st, 0)])
    print(line)
