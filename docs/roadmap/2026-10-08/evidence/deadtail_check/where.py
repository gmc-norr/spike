"""Where do aligned Q>=30 mismatches sit, real vs spike (K2b same templates)?
Splits by: read crashed (>=10 of last 20 quals < 15), any aligned mismatch in the previous 30 cycles
(seq order), lows (<Q15) among the previous 16 cycles, and mate. MAPQ>=20 primary, no-indel reads, Q100 +-5 masked.
Usage: python3 -I where.py BLOCKS REF VCF NAME=BAM ...
"""
import sys, collections
import numpy as np
import pysam
blocks_f, ref_f, vcf_f = sys.argv[1:4]
sets = [a.split('=', 1) for a in sys.argv[4:]]
CODE = np.full(256, 4, dtype=np.int8)
for i, b in enumerate(b'ACGT'):
    CODE[b] = i
fa = pysam.FastaFile(ref_f); vcf = pysam.VariantFile(vcf_f)
blocks = []
for line in open(blocks_f):
    c, s, e = line.split(); s, e = int(s) - 2000, int(e) + 2000
    code = CODE[np.frombuffer(fa.fetch(c, s, e).upper().encode(), dtype=np.uint8)]
    mask = np.zeros(len(code), dtype=bool)
    for rec in vcf.fetch(c, s, e):
        a = rec.pos - 1 - s; b = rec.pos - 1 + len(rec.ref) - s
        mask[max(a - 5, 0):max(b + 5, 0)] = True
    blocks.append(dict(chrom=c, start=s, code=code, mask=mask))
def fb(c, p):
    for b in blocks:
        if b['chrom'] == c and b['start'] + 200 <= p < b['start'] + len(b['code']) - 400:
            return b
for name, f in sets:
    bases = collections.Counter(); errs = collections.Counter(); reads = collections.Counter()
    for r in pysam.AlignmentFile(f):
        if r.is_unmapped or r.is_secondary or r.is_supplementary or r.mapping_quality < 20:
            continue
        b = fb(r.reference_name, r.reference_start)
        if b is None:
            continue
        cig = r.cigartuples
        if any(o not in (0, 4) for o, _ in cig):
            continue
        q = np.asarray(r.query_qualities); L = len(q)
        q0 = cig[0][1] if cig[0][0] == 4 else 0
        ln = sum(l for o, l in cig if o == 0)
        rp = r.reference_start - b['start']
        seq = CODE[np.frombuffer(r.query_sequence.encode(), dtype=np.uint8)]
        # per query index: aligned?, mismatch?, ok(compare)
        al = np.zeros(L, dtype=bool); al[q0:q0 + ln] = True
        rs = np.full(L, 4, dtype=np.int8); rs[q0:q0 + ln] = b['code'][rp:rp + ln]
        mk = np.ones(L, dtype=bool); mk[q0:q0 + ln] = b['mask'][rp:rp + ln]
        ok = al & ~mk & (rs < 4) & (seq < 4)
        mm = ok & (rs != seq)
        if r.is_reverse:
            q = q[::-1]; ok = ok[::-1]; mm = mm[::-1]
        crashed = int((q[-20:] < 15).sum() >= 10)
        mate = 1 if r.is_read1 else 2
        mmc = np.concatenate([[0], np.cumsum(mm)])
        low = (q < 15).astype(int); lowc = np.concatenate([[0], np.cumsum(low)])
        nhi = 0
        for c in range(L):
            if not ok[c] or q[c] < 30:
                continue
            prev_e = mmc[c] - mmc[max(c - 30, 0)]
            prev_l = lowc[c] - lowc[max(c - 16, 0)]
            key = (crashed, min(int(prev_e), 2), int(prev_l > 0))
            bases[key] += 1; errs[key] += int(mm[c]); nhi += int(mm[c])
            bases[('mate', mate)] += 1; errs[('mate', mate)] += int(mm[c])
        reads[(crashed, min(nhi, 2))] += 1
    tb = sum(v for k, v in bases.items() if k[0] != 'mate'); te = sum(v for k, v in errs.items() if k[0] != 'mate')
    print('== %s: hi-Q aligned bases %d errors %d rate %.3e' % (name, tb, te, te / tb))
    for k in sorted(k for k in bases if k[0] != 'mate'):
        print('  crashed=%d prev_err30=%s lows16>0=%d  bases %9d errs %5d rate %.2e  share_of_errs %.3f' % (
            k[0], ['0', '1', '2+'][k[1]], k[2], bases[k], errs[k], errs[k] / max(bases[k], 1), errs[k] / te))
    for m in (1, 2):
        print('  mate %d: rate %.3e (%d/%d)' % (m, errs[('mate', m)] / bases[('mate', m)], errs[('mate', m)], bases[('mate', m)]))
    n = sum(reads.values())
    print('  reads: ' + ' '.join('crashed=%d hiQmm=%s %.4f' % (k[0], ['0', '1', '2+'][k[1]], v / n) for k, v in sorted(reads.items())))
