"""K2b BAMs: 1-4 bp soft clips by read end (sequencing orientation), and how many lie within 10 bp of a
Q100 record (the sample's own variants; spike's K2b reads are made from the reference, so carry none).
Also: for 3' 1-4 bp clips, how many clipped bases are high-Q (>=30) mismatches.
Usage: python3 clipnear.py REF VCF NAME=BAM ...
"""
import sys, collections
import pysam

ref_f, vcf_f = sys.argv[1:3]
fa = pysam.FastaFile(ref_f)
vcf = pysam.VariantFile(vcf_f)
cache = {}


def near_var(chrom, pos):
    key = (chrom, pos // 1000)
    if key not in cache:
        s = max(key[1] * 1000 - 50, 0)
        cache[key] = [(r.pos - 1, r.pos - 1 + len(r.ref)) for r in vcf.fetch(chrom, s, key[1] * 1000 + 1050)]
    return any(a - 10 <= pos <= b + 10 for a, b in cache[key])


for arg in sys.argv[3:]:
    name, f = arg.split('=')
    c = collections.Counter(); n = 0
    for r in pysam.AlignmentFile(f):
        if r.is_unmapped or r.is_secondary or r.is_supplementary or r.mapping_quality < 20:
            continue
        n += 1
        cig = r.cigartuples
        lc = cig[0][1] if cig[0][0] == 4 else 0
        rc = cig[-1][1] if cig[-1][0] == 4 else 0
        # boundary positions on the reference
        left_b = r.reference_start; right_b = r.reference_end
        for side, L, bpos in (('L', lc, left_b), ('R', rc, right_b)):
            if not (1 <= L <= 4):
                continue
            end = ('5p' if side == 'L' else '3p') if not r.is_reverse else ('3p' if side == 'L' else '5p')
            nv = near_var(r.reference_name, bpos)
            c[(end, 'all')] += 1
            c[(end, 'nearQ100')] += int(nv)
            # high-Q mismatches within the clip
            q = r.query_qualities; s = r.query_sequence
            if side == 'L':
                idx = range(0, L); rp = [left_b - L + k for k in range(L)]
            else:
                idx = range(len(s) - L, len(s)); rp = [right_b + k for k in range(L)]
            refs = fa.fetch(r.reference_name, min(rp), max(rp) + 1).upper()
            hq_mm = sum(1 for k, p in zip(idx, rp) if q[k] >= 30 and s[k] != refs[p - min(rp)] and s[k] != 'N')
            c[(end, 'hq_mm')] += hq_mm
            if not nv:
                c[(end, 'hq_mm_novar')] += hq_mm
    print(name, n, ' | '.join('%s: %d reads (%.3f%%), near Q100 %d (%.0f%%), hiQ mismatches in clip %d (not near Q100 %d)' % (
        e, c[(e, 'all')], 100 * c[(e, 'all')] / n, c[(e, 'nearQ100')], 100 * c[(e, 'nearQ100')] / max(c[(e, 'all')], 1),
        c[(e, 'hq_mm')], c[(e, 'hq_mm_novar')]) for e in ('5p', '3p')))
