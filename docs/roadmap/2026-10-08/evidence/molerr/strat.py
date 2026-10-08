"""Overlap concordance stratified by the two mates' BQ: per stratum, overlapped bases,
A-only errors, B-only errors, both-wrong same, both-wrong different. Also R1-orientation
spectrum of high-Q concordant sites and of high-Q discordant (one-mate) errors."""
import sys, collections as C, numpy as np, pysam
blocks_f, ref_f, vcf_f = sys.argv[1:4]
fa = pysam.FastaFile(ref_f); vcf = pysam.VariantFile(vcf_f)
blocks = []
for line in open(blocks_f):
    c, s, e = line.split(); s, e = int(s) - 2000, int(e) + 2000
    ref = fa.fetch(c, s, e).upper()
    mask = np.zeros(len(ref), dtype=bool)
    for rec in vcf.fetch(c, s, e):
        a = rec.pos - 1 - s; b = rec.pos - 1 + len(rec.ref) - s
        mask[max(a - 5, 0):max(b + 5, 0)] = True
    blocks.append((c, s, ref, mask))
def fb(c, p):
    for blk in blocks:
        if blk[0] == c and blk[1] + 200 <= p < blk[1] + len(blk[2]) - 400: return blk
def amap(r):
    d = {}
    for q, p in r.get_aligned_pairs(matches_only=True):
        d[p] = (r.query_sequence[q], r.query_qualities[q])
    return d
comp = str.maketrans('ACGT', 'TGCA')
def strat(qa, qb):
    m = min(qa, qb)
    return 'both>=30' if m >= 30 else ('min15-29' if m >= 15 else 'min<15')
for a in sys.argv[4:]:
    name, f = a.split('=')
    st = C.defaultdict(lambda: [0, 0, 0, 0, 0])
    spec_conc = C.Counter(); spec_one = C.Counter(); dups = 0
    pend = {}
    for r in pysam.AlignmentFile(f):
        if r.is_unmapped or r.is_secondary or r.is_supplementary or r.mapping_quality < 20: continue
        if r.is_duplicate: dups += 1
        if any(o not in (0, 4) for o, _ in r.cigartuples): continue
        o = pend.pop(r.query_name, None)
        if o is None: pend[r.query_name] = r; continue
        if o.reference_name != r.reference_name: continue
        blk = fb(r.reference_name, min(o.reference_start, r.reference_start))
        if blk is None: continue
        lo = max(o.reference_start, r.reference_start); hi = min(o.reference_end, r.reference_end)
        if hi <= lo: continue
        A = amap(o); B = amap(r); c, s, ref, mask = blk
        r1 = o if o.is_read1 else r
        for p in range(lo, hi):
            i = p - s
            if i < 0 or i >= len(ref) or mask[i] or ref[i] == 'N' or p not in A or p not in B: continue
            (ba, qa), (bb, qb) = A[p], B[p]
            if qa < 10 or qb < 10 or ba == 'N' or bb == 'N': continue
            k = strat(qa, qb); row = st[k]; row[0] += 1
            ea = ba != ref[i]; eb = bb != ref[i]
            def r1o(alt):
                return (ref[i] + '>' + alt) if not r1.is_reverse else (ref[i].translate(comp) + '>' + alt.translate(comp))
            if ea and eb:
                if ba == bb:
                    row[3] += 1
                    if k == 'both>=30': spec_conc[r1o(ba)] += 1
                else: row[4] += 1
            elif ea:
                row[1] += 1
                if k == 'both>=30': spec_one[r1o(ba)] += 1
            elif eb:
                row[2] += 1
                if k == 'both>=30': spec_one[r1o(bb)] += 1
    print(name, 'duplicate-flagged records:', dups)
    for k in ('both>=30', 'min15-29', 'min<15'):
        n, x, y, same, diff = st[k]
        exp = x * y / max(n, 1)  # expected both-wrong under independence (any base)
        print('  %-9s bases %8d  A-only %5d  B-only %5d  both-same %4d  both-diff %3d  (indep. both-wrong ~%.2f)  same/(same+diff) %.2f  concordant share of erroneous positions %.3f'
              % (k, n, x, y, same, diff, exp, same / max(same + diff, 1), same / max(same + diff + x + y, 1)))
    print('  hiQ concordant R1-orientation:', spec_conc.most_common(8))
    print('  hiQ one-mate R1-orientation  :', spec_one.most_common(8))
