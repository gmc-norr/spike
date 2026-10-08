"""Gate B (b): does a local read-class mix plus global per-class rates reproduce each region's own
crash share? Classes cut at the pooled 20 blocks' quantiles (as spike does). Per block:
predicted = sum over classes of the block's class share x the pooled crash rate of that class.
Also the same for perfect reads. Usage: mix.py BAM"""
import sys, pysam, numpy as np
bam = pysam.AlignmentFile(sys.argv[1]); lens = dict(zip(bam.references, bam.lengths))
blocks = {}
for c in [f"chr{i}" for i in range(1, 21)]:
    s = int(lens[c] * 0.3); qs = []
    for r in bam.fetch(c, s, s + 100_000):
        if r.flag & 0xF0C or r.mapping_quality < 20: continue
        q = np.array(r.query_qualities); qs.append(q[::-1] if r.is_reverse else q)
    if len(qs) >= 1000: blocks[c] = qs
allq = [q for v in blocks.values() for q in v]
top = max(q.max() for q in allq)
cuts = np.quantile([q.mean() for q in allq], [0.01, 0.03, 0.08, 0.2, 0.4, 0.6, 0.8])
cls = lambda q: int(np.searchsorted(cuts, q.mean(), side="right"))
crash = lambda q: len(q) >= 20 and (q[-20:] < 15).sum() >= 10
perf = lambda q: (q == top).all()
tot = np.zeros(8); cr = np.zeros(8); pf = np.zeros(8)
for q in allq:
    k = cls(q); tot[k] += 1; cr[k] += crash(q); pf[k] += perf(q)
rate_c, rate_p = cr / np.maximum(1, tot), pf / np.maximum(1, tot)
print("pooled crash rate by class:", np.round(100 * rate_c, 2))
err_mix, err_glob = [], []
for c, qs in blocks.items():
    n = np.zeros(8)
    for q in qs: n[cls(q)] += 1
    share = n / n.sum()
    obs_c = 100 * np.mean([crash(q) for q in qs]); pred_c = 100 * (share * rate_c).sum(); glob_c = 100 * cr.sum() / tot.sum()
    obs_p = 100 * np.mean([perf(q) for q in qs]); pred_p = 100 * (share * rate_p).sum()
    err_mix.append(abs(obs_c - pred_c)); err_glob.append(abs(obs_c - glob_c))
    print(f"{c:6s} crashed: observed {obs_c:5.2f}%  local mix {pred_c:5.2f}%  genome {glob_c:5.2f}%   perfect: observed {obs_p:5.2f}%  local mix {pred_p:5.2f}%")
print(f"mean |error| of the crash share: local class mix {np.mean(err_mix):.2f} points, genome-wide {np.mean(err_glob):.2f} points")
