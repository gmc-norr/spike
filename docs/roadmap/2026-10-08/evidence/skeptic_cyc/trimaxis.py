"""Hospital blocks BAM: do reads cut short by their fragment end like mid-run or like end-of-run?
Low = Q<15 (spike LOW_Q). Sequencing orientation. Primary, MAPQ>=20.
Usage: python3 trimaxis.py BAM
"""
import sys, collections
import numpy as np
import pysam

bam = sys.argv[1]
CYC = 151
lowc = np.zeros(CYC + 1); totc = np.zeros(CYC + 1)        # full-length (151) reads, by physical cycle (1-based)
lowc150 = np.zeros(CYC + 1); totc150 = np.zeros(CYC + 1)
short_end = collections.defaultdict(lambda: [0, 0])          # length bin -> [low, total] over last 10 cycles
short_same = collections.defaultdict(lambda: [0.0, 0.0])     # length bin -> expected from full-length at same cycles
lens = collections.Counter()
last_q = collections.defaultdict(list)
for r in pysam.AlignmentFile(bam):
    if r.is_unmapped or r.is_secondary or r.is_supplementary or r.mapping_quality < 20:
        continue
    q = np.asarray(r.query_qualities)
    if r.is_reverse:
        q = q[::-1]
    L = len(q)
    lens[L] += 1
    lo = (q < 15)
    if L == 151:
        lowc[1:L + 1] += lo; totc[1:L + 1] += 1
        last_q['151'].append(q[-1])
    elif L == 150:
        lowc150[1:L + 1] += lo; totc150[1:L + 1] += 1
        last_q['150'].append(q[-1])
    elif 40 <= L < 140:
        b = (L // 20) * 20
        short_end[b][0] += int(lo[-10:].sum()); short_end[b][1] += 10
        last_q['short'].append(q[-1])
n = sum(lens.values())
print('reads', n, 'len151 %.3f len150 %.3f len<150 %.3f len<140 %.3f' % (lens[151] / n, lens[150] / n,
      sum(v for k, v in lens.items() if k < 150) / n, sum(v for k, v in lens.items() if k < 140) / n))
rate = lowc / np.maximum(totc, 1)
print('full-length 151: low share last 10 cycles %.3f%%, cycles 1-10 %.3f%%, 50-100 %.3f%%' % (
    100 * lowc[142:152].sum() / totc[142:152].sum(), 100 * lowc[1:11].sum() / totc[1:11].sum(), 100 * lowc[50:101].sum() / totc[50:101].sum()))
print('150-bp reads: low share last 10 cycles %.3f%%' % (100 * lowc150[141:151].sum() / totc150[141:151].sum()))
print('low share by cycle (151 reads): 141-145 %.3f 146-150 %.3f 151 %.3f' % (100 * lowc[141:146].sum() / totc[141:146].sum(), 100 * lowc[146:151].sum() / totc[146:151].sum(), 100 * rate[151]))
for b in sorted(k for k in short_end if not isinstance(k, tuple)):
    lw, t = short_end[b]
    # full-length reads' low share at the same physical cycles (b+1-10 .. b+20)
    a0, a1 = max(b - 9, 1), min(b + 20, 151)
    same = lowc[a0:a1 + 1].sum() / totc[a0:a1 + 1].sum()
    print('short reads len %d-%d: n=%d last-10 low share %.3f%%; full-length reads at cycles %d-%d %.3f%%' % (
        b, b + 19, t // 10, 100 * lw / t, a0, a1, 100 * same))
for k, v in last_q.items():
    v = np.asarray(v)
    print('last-base mean Q', k, '%.2f' % v.mean(), 'n=%d' % len(v))
