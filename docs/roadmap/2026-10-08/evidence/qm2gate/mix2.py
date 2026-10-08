"""Gate B for v2's per-event class mix: does the template's run mix explain each region's crash share as
well as the region's own read-class mix? 20 blocks of 100 kb (chr1-20 at 30%), primary MAPQ >= 20 reads.
Per block, predicted crash share = sum over bins of the block's share x the pooled rate, for: read class
(mix.py), the template's longest-run bin, and both together (class x run). Mean |error| over blocks.
Usage: mix2.py BAM REF"""
import sys, pysam, numpy as np, collections
bam = pysam.AlignmentFile(sys.argv[1]); fa = pysam.FastaFile(sys.argv[2]); lens = dict(zip(bam.references, bam.lengths))
comp = str.maketrans("ACGTN", "TGCAN")
def hp_bin(b, n):
    if b == "C": return 3 if n >= 7 else 2 if n >= 5 else 0
    return 3 if n >= 12 else 2 if n >= 9 else 1 if n >= 7 else 0
def read_hp(s):
    best, prev, n = 0, "", 0
    for b in s:
        n = n + 1 if b == prev and b != "N" else (b != "N"); prev = b; best = max(best, hp_bin(b, n))
    return best
blocks = {}
for c in [f"chr{i}" for i in range(1, 21)]:
    s = int(lens[c] * 0.3); rows = []
    ref = fa.fetch(c, s - 1000, s + 101_000).upper()
    for r in bam.fetch(c, s, s + 100_000):
        if r.flag & 0xF0C or r.mapping_quality < 20 or r.query_length != 151: continue
        q = np.array(r.query_qualities); q = q[::-1] if r.is_reverse else q
        cig = r.cigartuples; lead = cig[0][1] if cig[0][0] == 4 else 0
        a = (r.reference_start - lead) if not r.is_reverse else None
        if r.is_reverse:
            trail = cig[-1][1] if cig[-1][0] == 4 else 0; e = r.reference_end + trail
            t = ref[e - 151 - (s - 1000): e - (s - 1000)].translate(comp)[::-1]
        else:
            t = ref[a - (s - 1000): a + 151 - (s - 1000)]
        rows.append((q.mean(), (q[-20:] < 15).sum() >= 10, read_hp(t)))
    if len(rows) >= 1000: blocks[c] = rows
allr = [x for v in blocks.values() for x in v]
cuts = np.quantile([x[0] for x in allr], [0.01, 0.03, 0.08, 0.2, 0.4, 0.6, 0.8])
cls = lambda m: int(np.searchsorted(cuts, m, side="right"))
def rates(key):
    t = collections.Counter(); k = collections.Counter()
    for x in allr: t[key(x)] += 1; k[key(x)] += x[1]
    return {g: k[g] / t[g] for g in t}
keys = {"class": lambda x: cls(x[0]), "run": lambda x: x[2], "class x run": lambda x: (cls(x[0]), x[2])}
R = {n: rates(f) for n, f in keys.items()}
glob = np.mean([x[1] for x in allr])
err = collections.defaultdict(list)
for c, rows in blocks.items():
    obs = 100 * np.mean([x[1] for x in rows])
    line = f"{c:6s} observed {obs:5.2f}%"
    err["genome"].append(abs(obs - 100 * glob))
    for n, f in keys.items():
        pred = 100 * np.mean([R[n][f(x)] for x in rows]); err[n].append(abs(obs - pred)); line += f"  {n} {pred:5.2f}%"
    print(line)
print("mean |error| of the crash share, points: " + "  ".join(f"{n} {np.mean(v):.2f}" for n, v in err.items()))
