"""Where do aligned high-quality (Q>=30) mismatches sit: by cycles to the 3' end and by read class.
Usage: python3 hiq.py BLOCKS REF VCF NAME=BAM ...  (no-indel reads, MAPQ>=20, Q100 +-5 bp masked)
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


EB = [(0, 10), (10, 30), (30, 60), (60, 200)]
for name, f in sets:
    tot = np.zeros((4, 3), dtype=np.int64); mm = np.zeros((4, 3), dtype=np.int64)
    per_read = collections.Counter(); reads = 0; clipped = 0
    for r in pysam.AlignmentFile(f):
        if r.is_unmapped or r.is_secondary or r.is_supplementary or r.mapping_quality < 20:
            continue
        b = fb(r.reference_name, r.reference_start)
        if b is None:
            continue
        cig = r.cigartuples
        if any(o not in (0, 4) for o, _ in cig):
            continue
        reads += 1
        if len(cig) > 1:
            clipped += 1
        quals = np.asarray(r.query_qualities)
        L = len(quals)
        mq = quals.mean()
        cls = 0 if mq < 30 else (1 if mq < 35 else 2)
        q0 = cig[0][1] if cig[0][0] == 4 else 0
        ln = sum(l for o, l in cig if o == 0)
        rp = r.reference_start - b['start']
        seq = CODE[np.frombuffer(r.query_sequence.encode(), dtype=np.uint8)[q0:q0 + ln]]
        qq = quals[q0:q0 + ln]
        rs = b['code'][rp:rp + ln]; mk = b['mask'][rp:rp + ln]
        cyc = np.arange(q0, q0 + ln)
        left = (L - cyc) if not r.is_reverse else (cyc + 1)  # cycles to the 3' end incl. this base
        ok = (~mk) & (rs < 4) & (seq < 4) & (qq >= 30)
        err = ok & (rs != seq)
        per_read[int(err.sum())] += 1
        for k, (lo, hi) in enumerate(EB):
            sel = (left > lo) & (left <= hi)
            tot[k, cls] += int((ok & sel).sum()); mm[k, cls] += int((err & sel).sum())
    print('%s: no-indel reads %d (clipped %.2f%%)' % (name, reads, 100 * clipped / reads))
    for k, (lo, hi) in enumerate(EB):
        print('  cycles to 3p end %3d-%3d: ' % (lo + 1, hi) + '  '.join(
            'class %s %.2e (%d/%d)' % (['meanQ<30', '30-35', '>=35'][c], mm[k, c] / max(tot[k, c], 1), mm[k, c], tot[k, c]) for c in range(3)))
    n = sum(per_read.values())
    print('  reads by aligned Q>=30 mismatches: ' + ' '.join('%d:%.4f' % (k, per_read[k] / n) for k in range(0, 6)) + ' 6+:%.5f' % (sum(v for k, v in per_read.items() if k >= 6) / n))
