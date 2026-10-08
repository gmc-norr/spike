"""Clip test scoring (rule in README.txt, written before the run). For each BAM: primary, MAPQ >= 20,
proper pair, |TLEN| >= 151, 151 bp reads. A clip is bad-end when it is 1-4 bp or matches the reference
at >= 50% of its placed bases, and no Q100 variant lies within 10 bp of its boundary.
Reports per set: reads, bad-end share (95% Wilson CI) and z against the FIRST set, any clip, crashed,
mismatches per 100 aligned bases, and bad-end share by the template's longest-run bin.
Usage: clipscore.py REF Q100_VCF NAME=BAM [NAME=BAM ...]"""
import sys, math, collections, pysam, numpy as np

ref_path, vcf_path, *sets = sys.argv[1:]
fa = pysam.FastaFile(ref_path)
vcf = pysam.VariantFile(vcf_path)
comp = str.maketrans("ACGTN", "TGCAN")
var_cache = {}
def variants(chrom, pos):
    key = (chrom, pos // 50_000)
    if key not in var_cache:
        s = key[1] * 50_000
        var_cache[key] = np.array(sorted(v.start for v in vcf.fetch(chrom, max(0, s - 20), s + 50_020)), np.int64)
    return var_cache[key]
def near_variant(chrom, pos, w=10):
    v = variants(chrom, pos)
    i = np.searchsorted(v, pos)
    return any(0 <= j < len(v) and abs(int(v[j]) - pos) <= w for j in (i - 1, i))
def hp_bin(base, n):
    if base == "C":
        return 3 if n >= 7 else 2 if n >= 5 else 0
    return 3 if n >= 12 else 2 if n >= 9 else 1 if n >= 7 else 0
def read_hp(s):
    best, prev, n = 0, "", 0
    for b in s:
        n = n + 1 if b == prev and b != "N" else (b != "N")
        prev = b
        best = max(best, hp_bin(b, n))
    return best

def score(path):
    st = collections.Counter(); by_hp = collections.defaultdict(lambda: [0, 0])
    for r in pysam.AlignmentFile(path):
        if r.flag & 0x904 or not r.is_proper_pair or r.mapping_quality < 20 or abs(r.template_length) < 151 or r.query_length != 151:
            continue
        cig = r.cigartuples
        lead = cig[0][1] if cig[0][0] == 4 else 0
        trail = cig[-1][1] if cig[-1][0] == 4 else 0
        span = fa.fetch(r.reference_name, max(0, r.reference_start - lead), r.reference_end + trail).upper()
        seq = r.query_sequence
        bad = False
        for L, rs, s_, boundary in ((lead, r.reference_start - lead, seq[:lead], r.reference_start),
                                    (trail, r.reference_end, seq[len(seq) - trail:], r.reference_end)):
            if not L:
                continue
            refs = fa.fetch(r.reference_name, max(0, rs), rs + L).upper()
            ident = sum(a == b for a, b in zip(s_, refs)) / L
            if (L <= 4 or ident >= 0.5) and not near_variant(r.reference_name, boundary):
                bad = True
        q = np.array(r.query_qualities)
        qs = q[::-1] if r.is_reverse else q
        tmpl = span.translate(comp)[::-1] if r.is_reverse else span
        h = read_hp(tmpl)
        ins = sum(l for op, l in cig if op == 1); dl = sum(l for op, l in cig if op == 2); m = sum(l for op, l in cig if op in (0, 7, 8))
        st["reads"] += 1; st["bad"] += bad; st["clip"] += bool(lead or trail); st["crash"] += (qs[-20:] < 15).sum() >= 10
        st["mm"] += r.get_tag("NM") - ins - dl; st["aligned"] += m
        by_hp[h][0] += 1; by_hp[h][1] += bad
    return st, by_hp

def wilson(k, n):
    p = k / n; z = 1.96; d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d; h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return c - h, c + h
def ztest(k1, n1, k2, n2):
    p = (k1 + k2) / (n1 + n2); se = math.sqrt(p * (1 - p) * (1 / n1 + 1 / n2))
    return (k1 / n1 - k2 / n2) / se if se else 0.0

res = []
for item in sets:
    name, path = item.split("=", 1)
    res.append((name, *score(path)))
base = res[0][1]
print(f"{'set':6s} {'reads':>7s} {'bad-end clip':>13s} {'95% CI':>15s} {'z vs ' + res[0][0]:>10s} {'any clip':>9s} {'crashed':>8s} {'z':>6s} {'mm/100':>7s}")
for name, st, _ in res:
    n = st["reads"]; lo, hi = wilson(st["bad"], n)
    print(f"{name:6s} {n:7d} {100*st['bad']/n:12.2f}% {100*lo:6.2f}-{100*hi:5.2f}% {ztest(st['bad'], n, base['bad'], base['reads']):10.2f} "
          f"{100*st['clip']/n:8.2f}% {100*st['crash']/n:7.2f}% {ztest(st['crash'], n, base['crash'], base['reads']):6.2f} {100*st['mm']/st['aligned']:7.3f}")
print("bad-end clip share by the template's longest-run bin (0 / 1 / 2 / 3), reads in brackets:")
for name, _, by in res:
    print(f"  {name:6s} " + "  ".join(f"{100*by[h][1]/max(1,by[h][0]):6.2f}% ({by[h][0]})" for h in range(4)))
