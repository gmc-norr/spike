"""Are recurrent same-strand errors a property of the SITE? Split one BAM's reads in two halves
by name hash; find same-strand same-alt recurrent sites in half A; ask how often half B repeats
them, against the base rate at all sites.
Usage: python3 splithalf.py BAM REF VCF BENCH_BED chrom:start-end [...]
"""
import sys, zlib, collections
import numpy as np
import pysam

bam_f, ref_f, vcf_f, bed_f = sys.argv[1:5]
regions = sys.argv[5:]
CODE = np.full(256, 4, dtype=np.int8)
for i, b in enumerate(b'ACGT'):
    CODE[b] = i
fa = pysam.FastaFile(ref_f)
vcf = pysam.VariantFile(vcf_f)
bed = collections.defaultdict(list)
for line in open(bed_f):
    c, s, e = line.split()[:3]
    bed[c].append((int(s), int(e)))
bam = pysam.AlignmentFile(bam_f)
tot = collections.Counter()
for reg in regions:
    c, se = reg.split(':'); s, e = map(int, se.split('-'))
    code = CODE[np.frombuffer(fa.fetch(c, s, e).upper().encode(), dtype=np.uint8)]
    n = len(code)
    mask = np.zeros(n, dtype=bool)
    for rec in vcf.fetch(c, s, e):
        a = rec.pos - 1 - s; b = rec.pos - 1 + len(rec.ref) - s
        mask[max(a - 5, 0):max(b + 5, 0)] = True
    inbed = np.zeros(n, dtype=bool)
    for bs, be in bed.get(c, ()):
        if be > s and bs < e:
            inbed[max(bs - s, 0):min(be - s, n)] = True
    cnt = np.zeros((2, n, 2, 5), dtype=np.int32)  # half, pos, strand, base
    for r in bam.fetch(c, s, e):
        if r.flag & 0xF04 or r.mapping_quality < 20 or not r.is_proper_pair:
            continue
        cig = r.cigartuples
        if any(o not in (0, 4) for o, _ in cig):
            continue
        q0 = cig[0][1] if cig[0][0] == 4 else 0
        ln = sum(l for o, l in cig if o == 0)
        rp = r.reference_start - s
        if rp < 0 or rp + ln > n:
            continue
        seq = CODE[np.frombuffer(r.query_sequence.encode(), dtype=np.uint8)[q0:q0 + ln]]
        qual = np.asarray(r.query_qualities[q0:q0 + ln])
        ok = qual >= 10
        h = zlib.crc32(r.query_name.encode()) & 1
        np.add.at(cnt[h], (np.arange(rp, rp + ln)[ok], int(r.is_reverse), seq[ok]), 1)
    good = (~mask) & inbed & (code < 4)
    for p in np.nonzero(good)[0]:
        ref = code[p]
        for st in range(2):
            a_ = cnt[0, p, st, :4].copy(); b_ = cnt[1, p, st, :4].copy()
            if a_.sum() < 4 or b_.sum() < 4:
                continue
            a_[ref] = 0; b_[ref] = 0
            tot['sites'] += 1
            for alt in range(4):
                if alt == ref:
                    continue
                if b_[alt] >= 2:
                    tot['B>=2'] += 1
                if b_[alt] >= 1:
                    tot['B>=1'] += 1
                if a_[alt] >= 3:
                    tot['A>=3'] += 1
                    tot['A>=3 & B>=1'] += int(b_[alt] >= 1)
                    tot['A>=3 & B>=2'] += int(b_[alt] >= 2)
                    tot['A>=3 & other strand B>=1'] += int(cnt[1, p, 1 - st, alt] >= 1)
                if a_[alt] == 1:
                    tot['A==1'] += 1
                    tot['A==1 & B>=1'] += int(b_[alt] >= 1)
print(bam_f.split('/')[-1], dict(tot))
alts = 3 * tot['sites']
print('base rate per (site,strand,alt): P(B>=1) %.5f  P(B>=2) %.6f' % (tot['B>=1'] / alts, tot['B>=2'] / alts))
if tot['A>=3']:
    print('given A>=3 same-strand same-alt: P(B>=1) %.3f  P(B>=2) %.3f  (other strand in B >=1: %.3f); n=%d' % (
        tot['A>=3 & B>=1'] / tot['A>=3'], tot['A>=3 & B>=2'] / tot['A>=3'], tot['A>=3 & other strand B>=1'] / tot['A>=3'], tot['A>=3']))
if tot['A==1']:
    print('given A==1 (a single error): P(B>=1) %.4f; n=%d' % (tot['A==1 & B>=1'] / tot['A==1'], tot['A==1']))
