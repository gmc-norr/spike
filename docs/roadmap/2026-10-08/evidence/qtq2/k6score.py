"""K6 (reported): spike's reads against the real reads in the two windows of the 31-value BAM run."""
import math, numpy as np, pysam
def z(k1, n1, k2, n2):
    p = (k1 + k2) / (n1 + n2); se = math.sqrt(p * (1 - p) * (1 / n1 + 1 / n2)); return (k1 / n1 - k2 / n2) / se if se else 0.0
bam = pysam.AlignmentFile("k6/merged.bam")
reads = {"spike": [], "real": []}
for c, a, b in [("chr20", 38547500, 38552500), ("chr20", 38895000, 38915000)]:
    for r in bam.fetch(c, a, b):
        if r.flag & 0xF0C: continue
        q = np.array(r.query_qualities); q = q[::-1] if r.is_reverse else q
        reads["spike" if r.query_name.startswith("SPIKE_") else "real"].append(q)
top = max(int(q.max()) for q in reads["real"])
vals = sorted({int(x) for q in reads["real"] for x in q}); svals = sorted({int(x) for q in reads["spike"] for x in q})
o = {}
for k, v in reads.items():
    m = np.array([q.mean() for q in v])
    o[k] = (len(v), m.std(ddof=1), sum(int((q == top).all()) for q in v), sum(int((q[-20:] < 15).sum() >= 10) for q in v))
s, r = o["spike"], o["real"]
print(f"K6 (31-value BAM): real values {len(vals)}, spike's values {len(svals)} (all from the real set: {set(svals) <= set(vals)})")
print(f"reads spike {s[0]}, real {r[0]}; read-mean SD spike {s[1]:.2f}, real {r[1]:.2f} (ratio {s[1]/r[1]:.3f}); perfect spike {100*s[2]/s[0]:.2f}% real {100*r[2]/r[0]:.2f}% (z {z(s[2],s[0],r[2],r[0]):+.2f}); crashed spike {100*s[3]/s[0]:.2f}% real {100*r[3]/r[0]:.2f}% (z {z(s[3],s[0],r[3],r[0]):+.2f})")
