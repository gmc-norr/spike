"""Overlapping mates: do both mates show the same non-reference base at the same position more
often than independent errors would? (library errors shared by the molecule: PCR, damage)
Usage: python3 mateconc.py BLOCKS REF VCF NAME=namegrouped.bam ..."""
import sys, collections, numpy as np, pysam
blocks_f, ref_f, vcf_f = sys.argv[1:4]
fa = pysam.FastaFile(ref_f); vcf = pysam.VariantFile(vcf_f)
blocks = []
for line in open(blocks_f):
    c, s, e = line.split(); s, e = int(s) - 2000, int(e) + 2000
    ref = fa.fetch(c, s, e).upper()
    mask = np.zeros(len(ref), dtype=bool)
    for rec in vcf.fetch(c, s, e):
        a = rec.pos - 1 - s; b = rec.pos - 1 + len(rec.ref) - s
        mask[max(a - 5, 0):max(b + 5, 0)] = True
    blocks.append((c, s, ref, mask))
def fb(c, p):
    for blk in blocks:
        if blk[0] == c and blk[1] + 200 <= p < blk[1] + len(blk[2]) - 400: return blk
def amap(r):
    d = {}
    for q, p in r.get_aligned_pairs(matches_only=True):
        d[p] = (r.query_sequence[q], r.query_qualities[q])
    return d
for a in sys.argv[4:]:
    name, f = a.split('=')
    det = []
    pend = {}; ov = 0; m1 = m2 = conc = disc_same_pos = 0
    for r in pysam.AlignmentFile(f):
        if r.is_unmapped or r.is_secondary or r.is_supplementary or r.mapping_quality < 20: continue
        if any(o not in (0, 4) for o, _ in r.cigartuples): continue
        o = pend.pop(r.query_name, None)
        if o is None:
            pend[r.query_name] = r; continue
        if o.reference_name != r.reference_name: continue
        blk = fb(r.reference_name, min(o.reference_start, r.reference_start))
        if blk is None: continue
        lo = max(o.reference_start, r.reference_start); hi = min(o.reference_end, r.reference_end)
        if hi <= lo: continue
        A = amap(o); B = amap(r)
        c, s, ref, mask = blk
        for p in range(lo, hi):
            i = p - s
            if i < 0 or i >= len(ref) or mask[i] or ref[i] == 'N' or p not in A or p not in B: continue
            (ba, qa), (bb, qb) = A[p], B[p]
            if qa < 10 or qb < 10 or ba == 'N' or bb == 'N': continue
            ov += 1
            ea = ba != ref[i]; eb = bb != ref[i]
            m1 += ea; m2 += eb
            if ea and eb:
                if ba == bb:
                    conc += 1
                    r1 = o if o.is_read1 else r
                    comp = str.maketrans('ACGT', 'TGCA')
                    rb, ab = (ref[i], ba) if not r1.is_reverse else (ref[i].translate(comp), ba.translate(comp))
                    fs = min(o.reference_start, r.reference_start); fe = max(o.reference_end, r.reference_end)
                    dist = min(p - fs, fe - 1 - p)
                    det.append((rb + '>' + ab, dist, min(qa, qb)))
                else: disc_same_pos += 1
    exp = m1 * m2 / max(ov, 1) / 3
    print('%s: overlapped bases %d; mate A errors %d, mate B errors %d; both wrong same base %d (independent expectation ~%.2f), both wrong different base %d' % (name, ov, m1, m2, conc, exp, disc_same_pos))

    import collections as C
    print('  R1-orientation change:', C.Counter(d[0] for d in det).most_common(8))
    print('  distance to nearest fragment end: <=5 %d, 6-20 %d, >20 %d' % (sum(d[1] <= 5 for d in det), sum(5 < d[1] <= 20 for d in det), sum(d[1] > 20 for d in det)))
    print('  min BQ of the two: <15 %d, 15-29 %d, >=30 %d' % (sum(d[2] < 15 for d in det), sum(15 <= d[2] < 30 for d in det), sum(d[2] >= 30 for d in det)))
