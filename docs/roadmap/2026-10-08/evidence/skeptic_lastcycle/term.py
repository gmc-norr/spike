"""Terminal high-Q errors: per-cycle rate (aligned + clipped, placed), clumping vs independence, 5' clips near Q100.
Usage: python3 -I term.py BLOCKS REF VCF NAME=BAM ...  (primary, MAPQ>=20, M/S-only CIGAR, class = mean Q>=35)
"""
import sys, collections
import numpy as np
import pysam
blocks_f, ref_f, vcf_f = sys.argv[1:4]
sets = [a.split('=', 1) for a in sys.argv[4:]]
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
    for b in blocks:
        if b[0] == c and b[1] + 200 <= p < b[1] + len(b[2]) - 400:
            return b
for name, f in sets:
    K = 10
    hb = np.zeros(K + 1, dtype=np.int64); he = np.zeros(K + 1, dtype=np.int64); hea = np.zeros(K + 1, dtype=np.int64)
    sb = np.zeros(6, dtype=np.int64); se = np.zeros(6, dtype=np.int64)
    last4 = collections.Counter(); last4_clip = collections.Counter(); n_clean = 0
    c5 = 0; c5_near = 0; c3 = 0; c3_near = 0; n_all = 0
    c3_err2 = 0
    for r in pysam.AlignmentFile(f, check_sq=False):
        if r.is_unmapped or r.is_secondary or r.is_supplementary or r.mapping_quality < 20:
            continue
        b = fb(r.reference_name, r.reference_start)
        if b is None:
            continue
        cig = r.cigartuples
        if any(o not in (0, 4) for o, _ in cig):
            continue
        n_all += 1
        _, off, ref, mask = b
        seq = r.query_sequence; qual = r.query_qualities; L = len(seq)
        lc = cig[0][1] if cig[0][0] == 4 else 0
        rc = cig[-1][1] if cig[-1][0] == 4 else 0
        rev = r.is_reverse
        base0 = r.reference_start - off - lc   # ref index of query index 0
        def at(i):   # query index -> (refidx, mismatch, masked, clipped)
            ri = base0 + i
            if ri < 0 or ri >= len(ref):
                return None
            return ri, seq[i] != ref[ri] and seq[i] != 'N' and ref[ri] != 'N', bool(mask[ri]), (i < lc or i >= L - rc)
        five, three = (rc, lc) if rev else (lc, rc)
        # 5' / 3' short clips, near Q100?
        def near_end(qidx):
            return any(mask[base0 + i] for i in qidx if 0 <= base0 + i < len(ref))
        if 1 <= five <= 4:
            c5 += 1
            idx = range(L - 1, L - 15, -1) if rev else range(0, 14)
            c5_near += near_end(idx)
        if 1 <= three <= 4:
            c3 += 1
            idx = range(0, 14) if rev else range(L - 1, L - 15, -1)
            c3_near += near_end(idx)
        if np.mean(qual) < 35:
            continue
        n_clean += 1
        e4 = 0; ok4 = True
        for k in range(1, K + 1):   # cycles_left k
            i = (k - 1) if rev else (L - k)
            a = at(i)
            if a is None or a[2] or qual[i] < 30:
                if k <= 4: ok4 = False
                continue
            hb[k] += 1; he[k] += a[1]; hea[k] += a[1] and not a[3]
            if k <= 4: e4 += a[1]
        for c in range(1, 6):
            i = (L - c) if rev else (c - 1)
            a = at(i)
            if a is None or a[2] or qual[i] < 30: continue
            sb[c] += 1; se[c] += a[1]
        if ok4:
            last4[min(e4, 3)] += 1
            if 2 <= three <= 4: last4_clip[min(e4, 3)] += 1
    print('==', name, 'reads(M/S only)', n_all, 'clean(meanQ>=35)', n_clean)
    print('  5p 1-4bp clips %d (%.3f%%), with a Q100 record within ~14 bp of the 5p end: %d (%.0f%%)' % (c5, 100*c5/n_all, c5_near, 100*c5_near/max(c5,1)))
    print('  3p 1-4bp clips %d (%.3f%%), with a Q100 record within ~14 bp of the 3p end: %d (%.0f%%)' % (c3, 100*c3/n_all, c3_near, 100*c3_near/max(c3,1)))
    print('  clean reads, Q>=30 unmasked bases, by cycles_left: rate(placed incl clipped) / rate(aligned only) [errors]')
    for k in range(1, K + 1):
        print('    cl %2d  %.2e / %.2e  [%d / %d of %d]' % (k, he[k]/hb[k], hea[k]/hb[k], he[k], hea[k], hb[k]))
    print('  clean reads, cycles 1-5 from 5p: ' + ' '.join('c%d %.2e[%d]' % (c, se[c]/sb[c], se[c]) for c in range(1, 6)))
    p = he[1:5] / hb[1:5]
    exp2 = sum(p[i]*p[j] for i in range(4) for j in range(i+1, 4))
    n4 = sum(last4.values())
    print('  reads with all last 4 Q>=30 unmasked: %d; errors in last 4: 0:%d 1:%d 2:%d 3+:%d' % (n4, last4[0], last4[1], last4[2], last4[3]))
    print('  P(>=2 in last 4) observed %.2e (%d reads) vs independent from per-cycle rates %.2e (%.2f reads expected)' % (
        (last4[2]+last4[3])/n4, last4[2]+last4[3], exp2, exp2*n4))
    print('  of those >=2 reads, with a 2-4 bp 3p clip: %d' % (last4_clip[2]+last4_clip[3]))
