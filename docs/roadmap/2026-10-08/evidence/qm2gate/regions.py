"""Gate B (a): do quality statistics differ between genome regions? 20 blocks of 100 kb, one per
chromosome 1-20 at 30% of its length, MAPQ >= 20 primary non-duplicate reads: read-mean SD,
share of crashed reads (>= 10 of last 20 < Q15) and of perfect reads, per block.
Usage: regions.py BAM"""
import sys, pysam, numpy as np
bam = pysam.AlignmentFile(sys.argv[1])
lens = dict(zip(bam.references, bam.lengths))
rows = []
for c in [f"chr{i}" for i in range(1, 21)]:
    s = int(lens[c] * 0.3); qs = []
    for r in bam.fetch(c, s, s + 100_000):
        if r.flag & 0xF0C or r.mapping_quality < 20: continue
        q = np.array(r.query_qualities); qs.append(q[::-1] if r.is_reverse else q)
    if len(qs) < 1000: rows.append((c, len(qs), None, None, None)); continue
    top = max(q.max() for q in qs)
    m = np.array([q.mean() for q in qs])
    crash = np.mean([(q[-20:] < 15).sum() >= 10 for q in qs if len(q) >= 20])
    perf = np.mean([(q == top).all() for q in qs])
    rows.append((c, len(qs), m.std(ddof=1), 100 * crash, 100 * perf))
for r in rows:
    print(f"{r[0]:6s} reads {r[1]:6d}  " + ("too few" if r[2] is None else f"read-mean SD {r[2]:.2f}  crashed {r[3]:.2f}%  perfect {r[4]:.2f}%"))
ok = [r for r in rows if r[2] is not None]
for k, lab in ((2, "SD"), (3, "crashed %"), (4, "perfect %")):
    v = np.array([r[k] for r in ok]); print(f"{lab}: min {v.min():.2f} median {np.median(v):.2f} max {v.max():.2f}")
