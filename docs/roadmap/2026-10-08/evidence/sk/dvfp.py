"""DeepVariant SNV calls in the 20 blocks, inside the Q100 bench BED: how many are not Q100 SNVs (FP)?
Usage: python3 dvfp.py BLOCKS BENCH_BED Q100_VCF CALL_VCF"""
import sys, collections, bisect
import pysam
blocks_f, bed_f, q_f, c_f = sys.argv[1:5]
blocks = [(l.split()[0], int(l.split()[1]), int(l.split()[2])) for l in open(blocks_f)]
bed = collections.defaultdict(list)
for line in open(bed_f):
    c, s, e = line.split()[:3]
    bed[c].append((int(s), int(e)))
for c in bed:
    bed[c].sort()
starts = {c: [x[0] for x in v] for c, v in bed.items()}
def inbed(c, p0):
    i = bisect.bisect_right(starts.get(c, []), p0) - 1
    return i >= 0 and bed[c][i][0] <= p0 < bed[c][i][1]
q = pysam.VariantFile(q_f); cv = pysam.VariantFile(c_f)
res = collections.Counter(); bp = 0; fps = []
for c, s, e in blocks:
    truth = {}
    near = set()
    for r in q.fetch(c, s, e):
        for a in r.alts or ():
            truth[(r.pos, a)] = 1
        for p in range(r.pos - 5, r.pos + len(r.ref) + 5):
            near.add(p)
    for r in cv.fetch(c, s, e):
        if len(r.ref) != 1 or not r.alts or any(len(a) != 1 for a in r.alts):
            continue
        if not inbed(c, r.pos - 1):
            continue
        gt = r.samples[0].get('GT')
        if gt is None or all(g in (0, None) for g in gt):
            res['ref_or_nocall'] += 1
            continue
        filt = ','.join(r.filter.keys()) or '.'
        tp = any((r.pos, a) in truth for a in r.alts)
        key = ('TP' if tp else ('FP_near_q100' if r.pos in near else 'FP')) + ':' + filt
        res[key] += 1
        if not tp and r.pos not in near:
            fps.append((c, r.pos, r.ref, r.alts, gt, filt, r.qual, r.samples[0].get('AD'), r.samples[0].get('DP')))
print(dict(res))
for f in fps:
    print(f)
