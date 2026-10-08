import sys, random, bisect
def load(p):
    d = {}
    for line in open(p):
        pos, ref, alt, gt = line.rstrip('\n').split('\t')
        if ',' in alt: continue
        g = gt.replace('|', '/')
        d[int(pos)] = ('het' if g in ('0/1', '1/0') else 'hom' if g == '1/1' else None, alt)
    return d
def bed(p, chrom):
    iv = []
    for line in open(p):
        f = line.split('\t')
        if f[0] == chrom: iv.append((int(f[1]), int(f[2])))
    return sorted(iv)
a = load(sys.argv[1]); b = load(sys.argv[2])  # a = sample (HG002), b = wrong gVCF (HG001)
ivs = bed(sys.argv[3], 'chr20'); ivs2 = bed(sys.argv[4], 'chr20')
def inside(iv, x, y):
    i = bisect.bisect_right([s for s, _ in iv], x) - 1
    return i >= 0 and iv[i][1] >= y
random.seed(1)
apos = sorted(a); bpos = sorted(b)
n = 0; stats = []
while n < 200:
    s0, e0 = random.choice(ivs)
    if e0 - s0 < 4001: continue
    c = random.randrange(s0 + 2000, e0 - 2000)
    if not inside(ivs2, c - 2000, c + 2001): continue
    lo, hi = c - 2000, c + 2001
    A = {p: a[p] for p in apos[bisect.bisect_left(apos, lo):bisect.bisect_left(apos, hi)]}
    B = {p: b[p] for p in bpos[bisect.bisect_left(bpos, lo):bisect.bisect_left(bpos, hi)]}
    b_only_hom = sum(1 for p, (g, al) in B.items() if g == 'hom' and p not in A)   # wrong person's hom-alt, sample hom-ref -> false alt on remade reads
    a_hom_lost = sum(1 for p, (g, al) in A.items() if g == 'hom' and p not in B)    # sample hom-alt missing from wrong gVCF -> ref written
    a_het = sum(1 for p, (g, al) in A.items() if g == 'het')
    a_het_shared = sum(1 for p, (g, al) in A.items() if g == 'het' and p in B and B[p][0] == 'het' and B[p][1] == al)
    b_het = sum(1 for p, (g, al) in B.items() if g == 'het')
    stats.append((b_only_hom, a_hom_lost, a_het, a_het_shared, b_het, len(A), len(B)))
    n += 1
import statistics as st
def frac(f): return sum(1 for s in stats if f(s)) / len(stats)
print("windows", len(stats))
print("windows with >=1 wrong-person hom-alt at a sample hom-ref site:", frac(lambda s: s[0] > 0))
print("windows with >=1 sample hom-alt absent from wrong gVCF:", frac(lambda s: s[1] > 0))
print("windows where wrong gVCF has >=1 het (so it replaces pileup):", frac(lambda s: s[4] > 0))
tot_ahet = sum(s[2] for s in stats); tot_shared = sum(s[3] for s in stats); tot_bhet = sum(s[4] for s in stats)
print("sample hets", tot_ahet, "of which het with same alt in wrong gVCF", tot_shared, "wrong-gVCF hets", tot_bhet)
print("mean wrong hom per window", st.mean(s[0] for s in stats), "mean lost hom", st.mean(s[1] for s in stats))
# self control: same file as both
