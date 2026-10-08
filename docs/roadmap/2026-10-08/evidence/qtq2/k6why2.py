import pysam, numpy as np
def crash(r):
    q = np.array(r.query_qualities); q = q[::-1] if r.is_reverse else q
    return int((q[-20:] < 15).sum() >= 10)
s = pysam.AlignmentFile("k6/sim.bam")
t = {}
for r in s.fetch(until_eof=True):
    if r.flag & 0xF0C: continue
    k = ("SPIKE" if r.query_name.startswith("SPIKE_") else "orig", 1 if r.is_read1 else 2)
    a = t.setdefault(k, [0, 0]); a[0] += 1; a[1] += crash(r)
for k, (n, c) in sorted(t.items()): print(k, n, f"{100*c/n:.2f}%")
# spike reads by where they land
sp = [r for r in s.fetch(until_eof=True) if r.query_name.startswith("SPIKE_") and not r.flag & 0xF0C]
import collections
bins = collections.defaultdict(lambda: [0, 0])
for r in sp:
    b = "unmapped" if r.is_unmapped else f"{r.reference_name}:{r.reference_start//2000*2000}"
    bins[b][0] += 1; bins[b][1] += crash(r)
for b, (n, c) in sorted(bins.items(), key=lambda x: -x[1][0])[:15]: print(b, n, f"{100*c/n:.1f}%")
