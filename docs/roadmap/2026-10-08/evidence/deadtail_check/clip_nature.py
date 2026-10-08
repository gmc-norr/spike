"""Real K2b reads: the clipped high-Q errors counted by spike's rule -- near a Q100 record (beyond the
+-5 bp mask), at a clip boundary shared by >=3 same-strand reads, or neither? Split crashed / not crashed.
Same filters as where_hiq.py. Usage: python3 -I clip_nature.py BLOCKS REF VCF BAM
"""
import sys, collections, bisect
import numpy as np
import pysam

blocks_f, ref_f, vcf_f, bam_f = sys.argv[1:5]
CODE = np.full(256, 4, dtype=np.int8)
for i, b in enumerate(b'ACGT'):
    CODE[b] = i
fa = pysam.FastaFile(ref_f); vcf = pysam.VariantFile(vcf_f)
blocks = []
for line in open(blocks_f):
    c, s, e = line.split(); s, e = int(s) - 2000, int(e) + 2000
    code = CODE[np.frombuffer(fa.fetch(c, s, e).upper().encode(), dtype=np.uint8)]
    mask = np.zeros(len(code), dtype=bool)
    recs = []
    for rec in vcf.fetch(c, s, e):
        a = rec.pos - 1 - s; b = rec.pos - 1 + len(rec.ref) - s
        mask[max(a - 5, 0):max(b + 5, 0)] = True
        recs.append((a, b, len(rec.ref) != len(rec.alts[0]) if rec.alts else False))
    starts = sorted(x[0] for x in recs)
    blocks.append(dict(chrom=c, start=s, code=code, mask=mask, starts=starts))


def fb(c, p):
    for b in blocks:
        if b['chrom'] == c and b['start'] + 200 <= p < b['start'] + len(b['code']) - 400:
            return b


def nearest(b, i):
    st = b['starts']
    k = bisect.bisect_left(st, i)
    d = 10 ** 9
    for j in (k - 1, k):
        if 0 <= j < len(st):
            d = min(d, abs(st[j] - i))
    return d


events = []  # (block id, boundary ref idx, strand, crashed, n hi errors in clip, dist to nearest Q100 from boundary)
bound_count = collections.Counter()
for r in pysam.AlignmentFile(bam_f):
    if r.is_unmapped or r.is_secondary or r.is_supplementary or r.mapping_quality < 20 or not r.is_proper_pair:
        continue
    b = fb(r.reference_name, r.reference_start)
    if b is None:
        continue
    cig = r.cigartuples
    if any(o not in (0, 4) for o, _ in cig):
        continue
    bid = blocks.index(b)
    q = np.asarray(r.query_qualities); L = len(q)
    sqo = q[::-1] if r.is_reverse else q
    crashed = int((sqo[-20:] < 15).sum() >= 10)
    lead = cig[0][1] if cig[0][0] == 4 else 0
    trail = cig[-1][1] if (cig[-1][0] == 4 and len(cig) > 1) else 0
    ri = r.reference_start - b['start'] - lead + np.arange(L)
    inb = (ri >= 0) & (ri < len(b['code'])); ric = np.clip(ri, 0, len(b['code']) - 1)
    rs = np.where(inb, b['code'][ric], 4); mk = np.where(inb, b['mask'][ric], True)
    sq = CODE[np.frombuffer(r.query_sequence.encode(), dtype=np.uint8)]
    valid = (rs < 4) & (sq < 4)
    for sl, bound in ((slice(0, lead), r.reference_start - b['start']), (slice(L - trail, L), r.reference_end - b['start'])):
        n = sl.stop - sl.start
        if n == 0:
            continue
        bound_count[(bid, bound, r.is_reverse)] += 1
        v = valid[sl]
        if n > 4:
            comp = int(v.sum()); same = int(((sq[sl] == rs[sl]) & v).sum())
            if not (comp > 0 and 2 * same >= comp):
                continue
        ok = v & (~mk[sl]) & (q[sl] >= 30)
        e = int((ok & (sq[sl] != rs[sl])).sum())
        if e:
            events.append((bid, bound, r.is_reverse, crashed, e, nearest(b, bound), n))

for cr in (0, 1):
    tot = collections.Counter(); reads = collections.Counter()
    for (bid, bound, rev, c, e, d, n) in events:
        if c != cr:
            continue
        shared = bound_count[(bid, bound, rev)] >= 3
        near = 'Q100<=20' if d <= 20 else ('Q100 21-150' if d <= 150 else 'Q100 >150')
        key = ('shared-spot' if shared else 'unshared') + ', ' + near
        tot[key] += e; reads[key] += 1
    T = sum(tot.values())
    print('%s: clipped hi errors %d in %d clipped ends' % ('crashed' if cr else 'not-crashed', T, sum(reads.values())))
    for k in sorted(tot):
        print('   %-28s errors %5d (%.1f%%)  ends %d' % (k, tot[k], 100 * tot[k] / max(T, 1), reads[k]))
