"""Read-only reference arithmetic for the M13 skeptic check.
Deletion convention: 0-based half-open [S,E) removed.
demo VCF: POS = base before (1-based) -> S = POS; END = last deleted (1-based) -> E = END.
ClinVar list (cv_ldlr_dels.txt): columns start end len id; treated as S=start, E=end (BED-like)."""
import sys, random, bisect
import pysam

FA, RMSK, DEMO, CV = sys.argv[1:5]
fa = pysam.FastaFile(FA)
CH = "chr19"
LO, HI = 10900000, 11300000
seq = fa.fetch(CH, LO, HI).upper()


def ref(a, b):
    return seq[a - LO:b - LO]


def hom(s, e, cap=400):
    l = 0
    while l < cap and ref(s - l - 1, s - l) == ref(e - l - 1, e - l):
        l += 1
    r = 0
    while r < cap and ref(s + r, s + r + 1) == ref(e + r, e + r + 1):
        r += 1
    return l + r


alus = []
for line in open(RMSK):
    f = line.rstrip("\n").split("\t")
    if f[12] == "Alu":
        alus.append((int(f[6]), int(f[7]), f[9], f[10]))
alus.sort()
starts = [a[0] for a in alus]


def alu_at(p):
    i = bisect.bisect_right(starts, p) - 1
    while i >= 0 and i >= bisect.bisect_right(starts, p) - 5:
        a = alus[i]
        if a[0] <= p < a[1]:
            return a
        i -= 1
    return None


def alus_near(p, w):
    return [a for a in alus if a[1] > p - w and a[0] < p + w]


def same_or_pair(s, e, w):
    L = alus_near(s, w)
    R = alus_near(e, w)
    return any(a[2] == b[2] and a != b and a[0] < b[0] for a in L for b in R)


def lcs_exact(x, y):
    """longest common substring; returns (len, i_end_in_x, j_end_in_y)"""
    best = (0, 0, 0)
    prev = [0] * (len(y) + 1)
    for i in range(1, len(x) + 1):
        cur = [0] * (len(y) + 1)
        xi = x[i - 1]
        for j in range(1, len(y) + 1):
            if xi == y[j - 1]:
                v = prev[j - 1] + 1
                cur[j] = v
                if v > best[0]:
                    best = (v, i, j)
        prev = cur
    return best


def snap(s, e, w, minlen=10):
    """SIMPLEST: same-orientation Alu pair, A within w of s, B within w of e;
    junction at midpoint of the longest exact match >= minlen.
    Returns (s2, e2, matchlen, alu pair) for the pair whose snapped junction moves least."""
    best = None
    for a in alus_near(s, w):
        for b in alus_near(e, w):
            if a == b or a[2] != b[2] or a[1] > b[0]:
                continue
            x, y = ref(a[0], a[1]), ref(b[0], b[1])
            m, i, j = lcs_exact(x, y)
            if m < minlen:
                continue
            mid = m // 2
            s2 = a[0] + i - m + mid   # homologous position in A
            e2 = b[0] + j - m + mid   # same position in B
            move = abs(s2 - s) + abs(e2 - e)
            cand = (move, s2, e2, m, a[3], b[3])
            if best is None or cand < best:
                best = cand
    return best


def summary(label, ev):
    hs = sorted(hom(s, e) for s, e in ev)
    n = len(hs)
    print(f"{label}: n={n} median_hom={hs[n//2]} >=2bp={sum(h>=2 for h in hs)} blunt0={sum(h==0 for h in hs)} >=10bp={sum(h>=10 for h in hs)}")


demo = []
for line in open(DEMO):
    if line.startswith("#"):
        continue
    f = line.rstrip("\n").split("\t")
    end = int([x for x in f[7].split(";") if x.startswith("END=")][0][4:])
    demo.append((f[2], f[7].split(";")[0][7:], int(f[1]), end))
dels = [(s, e) for _, t, s, e in demo if t == "DEL"]
summary("demo DELs", dels)

cv = []
for line in open(CV):
    f = line.split()
    cv.append((int(f[0]), int(f[1]), f[3]))
summary("ClinVar list DELs (as S=start,E=end)", [(s, e) for s, e, _ in cv])
summary("ClinVar list DELs (as S=start-1,E=end-1)", [(s - 1, e - 1) for s, e, _ in cv])
summary("ClinVar list DELs (as S=start,E=end-1)", [(s, e - 1) for s, e, _ in cv])
summary("ClinVar list DELs (as S=start-1,E=end)", [(s - 1, e) for s, e, _ in cv])

# Alu context
def alu_stats(label, ev):
    n = len(ev)
    both = sum(1 for s, e in ev if alu_at(s) and alu_at(e - 1))
    one = sum(1 for s, e in ev if (alu_at(s) is None) != (alu_at(e - 1) is None))
    pair500 = sum(1 for s, e in ev if same_or_pair(s, e, 500))
    print(f"{label}: n={n} bothInAlu={both} oneInAlu={one} sameOrientPair500={pair500}")


alu_stats("demo all 22", [(s, e) for _, _, s, e in demo])
alu_stats("ClinVar list", [(s, e) for s, e, _ in cv])

random.seed(1)
lens = [e - s for _, _, s, e in demo]
bg = []
for _ in range(3000):
    L = random.choice(lens)
    s = random.randint(11075000, 11140000 - L)
    bg.append((s, s + L))
alu_stats("random same-length events in LDLR +-", bg)
summary("random same-length events", bg)

# Snap demo records
print("\nSNAP (W=1000, minlen 10) of demo records")
exons = {}
for line in open("/home/parlar_ai/dev/spike/data/ldlr_deletions/ldlr_exons_hg38.bed"):
    f = line.split()
    exons[int(f[3].split("exon")[1])] = (int(f[1]), int(f[2]))


def exons_hit(s, e):
    full = [k for k, (a, b) in sorted(exons.items()) if a >= s and b <= e]
    part = [k for k, (a, b) in sorted(exons.items()) if (a < e and b > s) and not (a >= s and b <= e)]
    return full, part


ok = 0
for name, t, s, e in demo:
    r = snap(s, e, 1000)
    f0, p0 = exons_hit(s, e)
    if r is None:
        print(f"{name}\t{t}\tno snap\tstated exons full={f0} part={p0}")
        continue
    ok += 1
    move, s2, e2, m, an, bn = r
    f1, p1 = exons_hit(s2, e2)
    h = hom(s2, e2)
    print(f"{name}\t{t}\tdS={s2-s}\tdE={e2-e}\tmatch={m}\thom_after={h}\t{an}/{bn}\texons {f0}{'+part'+str(p0) if p0 else ''} -> {f1}{'+part'+str(p1) if p1 else ''}\t{'SAME' if (f0,p0)==(f1,p1) else 'CHANGED'}")
print("snapped", ok, "of", len(demo))

# Accuracy: does snapping a rounded ClinVar event land near its stated (real) junction?
print("\nACCURACY on ClinVar list: round ends to nearest 100, snap W=1000, distance to stated")
for s, e, cid in cv:
    a_s, a_e = alu_at(s), alu_at(e - 1)
    rs, re_ = round(s, -2), round(e, -2)
    r = snap(rs, re_, 1000)
    h = hom(s, e)
    tag = f"{cid}\thom_stated={h}\tAlu_ends={'Y' if a_s else 'n'}{'Y' if a_e else 'n'}"
    if r is None:
        print(tag, "\tno snap")
    else:
        move, s2, e2, m, an, bn = r
        print(f"{tag}\tsnap dS={s2-s}\tdE={e2-e}\tmatch={m}\t{an}/{bn}")
