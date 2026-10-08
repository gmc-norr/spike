"""Gate B (c): at ~200k pairs (20 blocks of 100 kb, hospital BAM), how many low-quality bases fall
in each crashed-tail error cell (class, Q, run bin 4-7 / 8+, end bin 1-10 / 11-30)? The table needs 200."""
import pysam, numpy as np, collections
B = "/home/parlar_ai/dev/cnv_validation/from_rv/raredisease_seracare_na12878_resample/D24-14230_Seq25-7600_30x_resample/raredisease_results/alignment/D24-14230_Seq25-7600_30x_sorted_md.bam"
bam = pysam.AlignmentFile(B); lens = dict(zip(bam.references, bam.lengths)); qs = []
for c in [f"chr{i}" for i in range(1, 21)]:
    s = int(lens[c] * 0.3)
    for r in bam.fetch(c, s, s + 100_000):
        if r.flag & 0xF0C or r.mapping_quality < 20: continue
        q = np.array(r.query_qualities); qs.append(q[::-1] if r.is_reverse else q)
cuts = np.quantile([q.mean() for q in qs], [0.01, 0.03, 0.08, 0.2, 0.4, 0.6, 0.8])
cells = collections.Counter()
for q in qs:
    k = int(np.searchsorted(cuts, q.mean(), side="right")); run = 0; n = len(q)
    for c, x in enumerate(q):
        run = run + 1 if x < 15 else 0
        if x < 15 and run >= 4 and n - c <= 30:
            cells[(k, int(x), "run 4-7" if run < 8 else "run 8+", "end 1-10" if n - c <= 10 else "end 11-30")] += 1
print(f"reads {len(qs)} (~{len(qs)//2} pairs)")
full = sum(1 for v in cells.values() if v >= 200)
print(f"crashed-tail cells with any base: {len(cells)}; with >= 200: {full}")
for k, v in sorted(cells.items(), key=lambda kv: -kv[1])[:12]: print(f"  {k}: {v}")
