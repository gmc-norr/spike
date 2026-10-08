"""Concordant same-base overlap errors (mateconc2 logic), then for each site count other reads in a
second (full) BAM that carry the same alt, split by min BQ of the pair. Usage: BLOCKS REF VCF PAIRBAM FULLBAM"""
import sys, collections as C, numpy as np, pysam
blocks_f, ref_f, vcf_f, pair_f, full_f = sys.argv[1:6]
fa = pysam.FastaFile(ref_f); vcf = pysam.VariantFile(vcf_f)
blocks = []
for line in open(blocks_f):
    c, s, e = line.split(); s, e = int(s) - 2000, int(e) + 2000
    ref = fa.fetch(c, s, e).upper(); mask = np.zeros(len(ref), dtype=bool)
    for rec in vcf.fetch(c, s, e):
        a = rec.pos - 1 - s; b = rec.pos - 1 + len(rec.ref) - s
        mask[max(a - 5, 0):max(b + 5, 0)] = True
    blocks.append((c, s, ref, mask))
def fb(c, p):
    for blk in blocks:
        if blk[0] == c and blk[1] + 200 <= p < blk[1] + len(blk[2]) - 400: return blk
def amap(r):
    return {p: (r.query_sequence[q], r.query_qualities[q]) for q, p in r.get_aligned_pairs(matches_only=True)}
pend = {}; sites = []; disc = []
for r in pysam.AlignmentFile(pair_f):
    if r.is_unmapped or r.is_secondary or r.is_supplementary or r.mapping_quality < 20: continue
    if any(o not in (0, 4) for o, _ in r.cigartuples): continue
    o = pend.pop(r.query_name, None)
    if o is None: pend[r.query_name] = r; continue
    if o.reference_name != r.reference_name: continue
    blk = fb(r.reference_name, min(o.reference_start, r.reference_start))
    if blk is None: continue
    lo = max(o.reference_start, r.reference_start); hi = min(o.reference_end, r.reference_end)
    if hi <= lo: continue
    A = amap(o); B = amap(r); c, s, ref, mask = blk
    for p in range(lo, hi):
        i = p - s
        if i < 0 or i >= len(ref) or mask[i] or ref[i] == 'N' or p not in A or p not in B: continue
        (ba, qa), (bb, qb) = A[p], B[p]
        if qa < 10 or qb < 10 or ba == 'N' or bb == 'N': continue
        if ba != ref[i] and bb != ref[i]:
            (sites if ba == bb else disc).append((c, p, ba, min(qa, qb), r.query_name))
full = pysam.AlignmentFile(full_f)
def others(c, p, alt, name):
    n_alt = n_tot = 0
    for col in full.pileup(c, p, p + 1, truncate=True, min_base_quality=10, stepper='nofilter', ignore_overlaps=False, max_depth=100000):
        for pr in col.pileups:
            if pr.is_del or pr.is_refskip or pr.alignment.is_secondary or pr.alignment.is_supplementary: continue
            if pr.alignment.query_name == name: continue
            n_tot += 1; n_alt += pr.alignment.query_sequence[pr.query_position] == alt
    return n_alt, n_tot
tab = C.Counter()
for (c, p, alt, mq, name) in sites:
    n_alt, n_tot = others(c, p, alt, name)
    tab[('Q>=30' if mq >= 30 else 'Q<30', min(n_alt, 3))] += 1
print('concordant sites', len(sites), 'discordant both-wrong', len(disc))
for k in sorted(tab): print(' ', k, tab[k])
tab2 = C.Counter()
for (c, p, alt, mq, name) in disc:
    n_alt, n_tot = others(c, p, alt, name)
    tab2[min(n_alt, 3)] += 1
print(' discordant: other reads with mate-A alt', dict(tab2))
