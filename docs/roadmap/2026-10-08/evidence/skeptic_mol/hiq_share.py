"""In overlaps: per-mate errors by the mate's own Q (>=30 or not), how many are concordant with the other mate;
background share of overlapped positions with both mates Q>=30. Usage: BLOCKS REF VCF NAME=BAM..."""
import sys, numpy as np, pysam
blocks_f, ref_f, vcf_f = sys.argv[1:4]
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
for a in sys.argv[4:]:
    name, f = a.split('=')
    pend = {}; ov = both_hi = 0; err_hi = err_lo = conc_hi = conc_lo = 0
    for r in pysam.AlignmentFile(f):
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
            ov += 1; both_hi += (qa >= 30 and qb >= 30)
            for (x, qx, y) in ((ba, qa, bb), (bb, qb, ba)):
                if x != ref[i]:
                    if qx >= 30: err_hi += 1; conc_hi += (y == x)
                    else: err_lo += 1; conc_lo += (y == x)
    print('%s: overlapped %d, both mates Q>=30 share %.3f; mate errors at own Q>=30: %d (concordant with mate %d, %.1f%%); at Q<30: %d (concordant %d, %.1f%%)' % (name, ov, both_hi / ov, err_hi, conc_hi, 100 * conc_hi / max(err_hi, 1), err_lo, conc_lo, 100 * conc_lo / max(err_lo, 1)))
