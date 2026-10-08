"""Where do counted high-Q (Q>=30) errors sit, real vs spike, on K2b's same templates?
Spike's learning rule (aligned bases + clipped ends that are 1-4 bp or >=50% reference), Q100 +-5 bp masked.
No-indel reads only (CIGAR of M and S), primary, MAPQ>=20, proper pair.
Splits: aligned / clip 1-4 bp / clip 5+ bp (bad end), crashed (10+ of last 20 quals < Q15) or not, read class.
Usage: python3 -I where_hiq.py BLOCKS REF VCF NAME=BAM ...
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
    bases = collections.Counter(); errs = collections.Counter()
    reads = collections.Counter()
    per_read_hi = collections.Counter()  # (crashed, n counted hi errors incl clips)
    for r in pysam.AlignmentFile(f):
        if r.is_unmapped or r.is_secondary or r.is_supplementary or r.mapping_quality < 20 or not r.is_proper_pair:
            continue
        b = fb(r.reference_name, r.reference_start)
        if b is None:
            continue
        cig = r.cigartuples
        if any(o not in (0, 4) for o, _ in cig):
            continue
        q = np.asarray(r.query_qualities)
        L = len(q)
        seqorder_q = q[::-1] if r.is_reverse else q
        crashed = int((seqorder_q[-20:] < 15).sum() >= 10)
        mq = q.mean(); cls = 0 if mq < 30 else (1 if mq < 35 else 2)
        lead = cig[0][1] if cig[0][0] == 4 else 0
        trail = cig[-1][1] if (cig[-1][0] == 4 and len(cig) > 1) else 0
        ri = r.reference_start - b['start'] - lead + np.arange(L)
        inb = (ri >= 0) & (ri < len(b['code']))
        ric = np.clip(ri, 0, len(b['code']) - 1)
        rs = np.where(inb, b['code'][ric], 4); mk = np.where(inb, b['mask'][ric], True)
        sq = CODE[np.frombuffer(r.query_sequence.encode(), dtype=np.uint8)]
        kind = np.zeros(L, dtype=np.int8)  # 0 aligned, 1 short clip (1-4), 2 long clip ok, 3 long clip foreign
        valid = (rs < 4) & (sq < 4)

        def judge(sl):
            n = sl.stop - sl.start
            if n == 0:
                return
            if n <= 4:
                kind[sl] = 1; return
            v = valid[sl]
            comp = int(v.sum()); same = int(((sq[sl] == rs[sl]) & v).sum())
            kind[sl] = 2 if (comp > 0 and 2 * same >= comp) else 3
        judge(slice(0, lead)); judge(slice(L - trail, L))
        ok = valid & (~mk) & (kind != 3) & (q >= 30)
        err = ok & (sq != rs)
        reads[(crashed, cls)] += 1
        n_hi = 0
        for k, lab in enumerate(['aligned', 'clip1-4', 'clip5+']):
            sel = ok & (kind == k)
            bases[(lab, crashed, cls)] += int(sel.sum())
            e = int((err & (kind == k)).sum())
            errs[(lab, crashed, cls)] += e
            n_hi += e
        per_read_hi[(crashed, min(n_hi, 6))] += 1
    tot_b = sum(bases.values()); tot_e = sum(errs.values())
    print('=== %s: reads %d (crashed %d = %.2f%%); counted hi bases %d, errors %d (%.3e/base)' % (
        name, sum(reads.values()), sum(v for (c, _), v in reads.items() if c), 100 * sum(v for (c, _), v in reads.items() if c) / sum(reads.values()),
        tot_b, tot_e, tot_e / tot_b))
    for lab in ['aligned', 'clip1-4', 'clip5+']:
        for cr in (0, 1):
            e = sum(errs[(lab, cr, c)] for c in range(3)); n = sum(bases[(lab, cr, c)] for c in range(3))
            print('  %-8s %-11s errors %5d (%.1f%% of all)  bases %9d  rate %.2e   by class <30/30-35/>=35: %s' % (
                lab, 'crashed' if cr else 'not-crashed', e, 100 * e / tot_e, n, e / max(n, 1),
                ' / '.join('%d' % errs[(lab, cr, c)] for c in range(3))))
    for cr in (0, 1):
        n = sum(v for (c, k), v in per_read_hi.items() if c == cr)
        print('  %s reads by counted hi errors: ' % ('crashed    ' if cr else 'not-crashed') +
              ' '.join('%d:%.4f' % (k, per_read_hi[(cr, k)] / max(n, 1)) for k in range(7)))
