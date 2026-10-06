"""Share of primary mapped reads with a soft clip in a BAM, with a Wilson 95% interval. Usage: count_clips.py BAM [name-prefix]"""
import sys, math, collections, pysam
def wilson(k, n, z=1.96):
    c = (k + z*z/2) / (n + z*z); h = z*math.sqrt(k*(n-k)/n + z*z/4) / (n + z*z); return 100*(c-h), 100*(c+h)
st = collections.Counter()
for r in pysam.AlignmentFile(sys.argv[1]):
    if r.is_secondary or r.is_supplementary or r.is_unmapped: continue
    m = r.query_name.split("_")[0]; st[(m, "n")] += 1; st[(m, "k")] += any(op == 4 for op, _ in r.cigartuples)
for m in sorted({m for m, _ in st}):
    k, n = st[(m, "k")], st[(m, "n")]; lo, hi = wilson(k, n)
    print(f"{m}: soft-clipped {100*k/n:.2f}% ({lo:.2f}-{hi:.2f}) n={k}/{n}")
