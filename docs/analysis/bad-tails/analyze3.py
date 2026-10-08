"""What happens in the bad tails? Step 4: does the DNA sequence set off the crash at fixed places?
For every crashed read with a 3' soft clip, the clip boundary (first clipped base's reference
position, in sequencing direction). If crashes were random, few would share a boundary on the same
strand. Also: the sample's true variants near those boundaries (GIAB Q100), and, at shared
boundaries, whether ALL reads on that strand (crashed or not) err more after the spot.
Usage: analyze3.py IN.pkl SLICE_BAM REF Q100_VCF"""
import sys, pickle, collections, numpy as np, pysam

R = pickle.load(open(sys.argv[1], "rb"))
bam = pysam.AlignmentFile(sys.argv[2]); fa = pysam.FastaFile(sys.argv[3]); vcf = pysam.VariantFile(sys.argv[4])
rng = np.random.default_rng(7)

cr = [r for r in R if r["crashed"] and r["clip3"] > 0]
print(f"crashed reads with a 3' clip: {len(cr)}; their causes:", collections.Counter(r["cause"] for r in cr).most_common())
def boundary(r):
    return int(r["pos"][151 - r["clip3"]])  # first clipped base, placed
site = collections.Counter()
for r in cr:
    site[(r["rev"], boundary(r))] += 1
def shared(r, tol=2):
    b = boundary(r)
    return sum(site[(r["rev"], b + d)] for d in range(-tol, tol + 1)) - 1
sh = np.array([shared(r) for r in cr])
print(f"\n9a. crashed 3'-clipped reads sharing their boundary (+-2 bp, same strand) with >= 1 other: {100*np.mean(sh>=1):.1f}%; with >= 2 others: {100*np.mean(sh>=2):.1f}%")
# the same with each read's boundary moved to a random cycle 100-150 of the same read (coverage kept)
for rep in range(3):
    s2 = collections.Counter(); bb = []
    for r in cr:
        c = int(rng.integers(100, 151)); p = int(r["pos"][c])
        bb.append((r["rev"], p)); s2[(r["rev"], p)] += 1
    sh2 = np.array([sum(s2[(rv, p + d)] for d in range(-2, 3)) - 1 for rv, p in bb])
    print(f"    shuffled boundaries, try {rep+1}: >= 1 other {100*np.mean(sh2>=1):.1f}%; >= 2 others {100*np.mean(sh2>=2):.1f}%")
opp = np.array([sum(site[(not r["rev"], boundary(r) + d)] for d in range(-2, 3)) for r in cr])
print(f"    on the OTHER strand at the same spot: >= 1 {100*np.mean(opp>=1):.1f}%")

# true variants near boundaries
lo = min(int(r["pos"][r["pos"] >= 0].min()) for r in R) ; hi = max(int(r["pos"].max()) for r in R)
var = sorted(v.start for v in vcf.fetch("chr20", lo, hi))
var = np.array(var)
def near(p, w=10):
    i = np.searchsorted(var, p)
    d = min(abs(var[j] - p) for j in (i - 1, i) if 0 <= j < len(var))
    return d <= w
print("\n9b. a true HG002 variant within 10 bp of the boundary:")
for lab, sel in (("shared boundary (>= 2 others)", sh >= 2), ("lone boundary (no other)", sh == 0)):
    rs = [r for r, ok in zip(cr, sel) if ok]
    print(f"    {lab:30s} reads {len(rs):5d}: {100*np.mean([near(boundary(r)) for r in rs]):.1f}%")
rp = rng.integers(lo, hi, 5000)
print(f"    {'random positions':30s}        : {100*np.mean([near(int(p)) for p in rp]):.1f}%")

# sequence just before the boundary, in sequencing direction
comp = str.maketrans("ACGT", "TGCA")
def before(rev, p, k=12):
    if not rev:
        return fa.fetch("chr20", p - k, p).upper()
    return fa.fetch("chr20", p + 1, p + 1 + k).upper().translate(comp)[::-1]
def after(rev, p, k=6):
    if not rev:
        return fa.fetch("chr20", p, p + k).upper()
    return fa.fetch("chr20", p + 1 - k, p + 1).upper().translate(comp)[::-1]
spots = sorted({(rv, b) for (rv, b), n in site.items() if n >= 3})
print(f"\n9c. spots where >= 3 crashed reads clip at one boundary on one strand: {len(spots)}")
def feats(seqs):
    out = collections.OrderedDict()
    out["G share"] = np.mean([s.count("G") / len(s) for s in seqs])
    out["has GGC"] = np.mean(["GGC" in s for s in seqs])
    out["has GG"] = np.mean(["GG" in s for s in seqs])
    out["homopolymer >= 4"] = np.mean([any(b * 4 in s for b in "ACGT") for s in seqs])
    out["G-run >= 3"] = np.mean(["GGG" in s for s in seqs])
    return out
bs = [before(rv, b) for rv, b in spots]
rnd = []
for _ in range(3000):
    p = int(rng.integers(lo + 50, hi - 50)); rnd.append(before(bool(rng.integers(0, 2)), p))
fs, fr = feats(bs), feats(rnd)
print("    12 bp before the spot, in sequencing direction:   " + "  ".join(f"{k}: {fs[k]:.2f} (random {fr[k]:.2f})" for k in fs))
print("    examples (12 bp before | 6 bp after):")
for rv, b in spots[:12]:
    print(f"      {'-' if rv else '+'} chr20:{b+1:<10d} {before(rv, b)} | {after(rv, b)}   crashed reads clipping here: {site[(rv, b)]}")

# at the spots: every read on that strand that covers the spot with 10+ bases on each side
print("\n9d. at those spots, all reads (crashed or not) on the strand that reads into the spot, against the other strand:")
acc = {"same": [0, 0, 0, 0], "other": [0, 0, 0, 0]}  # mismatches before, bases before, mismatches after, bases after
qacc = {"same": [0, 0, 0, 0], "other": [0, 0, 0, 0]}  # low-quality share before/after
for rv, b in spots:
    for r in bam.fetch("chr20", max(0, b - 1), b + 1):
        if r.flag & 0xF0C or r.mapping_quality < 20:
            continue
        pairs = {rp_: qp for qp, rp_ in r.get_aligned_pairs(matches_only=True)}
        if b not in pairs:
            continue
        refseq = None
        side = "same" if r.is_reverse == rv else "other"
        # "before" = the 10 bases read before the spot in that read's own sequencing direction
        step = 1 if not r.is_reverse else -1
        bef = [b - step * k for k in range(1, 11)]; aft = [b + step * k for k in range(0, 10)]
        seq = r.query_sequence; q = r.query_qualities
        for lst, (i_mm, i_n) in ((bef, (0, 1)), (aft, (2, 3))):
            for p in lst:
                qp = pairs.get(p)
                if qp is None:
                    continue
                rb = fa.fetch("chr20", p, p + 1).upper()
                acc[side][i_n] += 1; acc[side][i_mm] += seq[qp] != rb
                qacc[side][i_n] += 1; qacc[side][i_mm] += q[qp] < 15
for side in ("same", "other"):
    a, qa = acc[side], qacc[side]
    print(f"    {side:5s} strand: wrong before {a[0]/max(1,a[1]):.4f}  after {a[2]/max(1,a[3]):.4f}   "
          f"low quality before {qa[0]/max(1,qa[1]):.4f}  after {qa[2]/max(1,qa[3]):.4f}   ({a[1]} / {a[3]} bases)")
