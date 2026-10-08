"""Vectorised re-count of recur.py's numbers, with a base-quality threshold argument.
Usage: python3 -I recur_fast.py BLOCKS REF VCF BED MINBQ NAME=BAM ...
Same filters as phys/recur.py: primary, MAPQ>=20, CIGAR only M/S, Q100 +-5 bp masked, bench BED, depth>=8.
Extra: same-strand ANY-alt (>=3 mismatches on one strand regardless of base) to separate the
substitution-matrix share from the site share.
"""
import sys, collections
import numpy as np
import pysam

blocks_f, ref_f, vcf_f, bed_f, minbq = sys.argv[1:6]
minbq = int(minbq)
sets = [a.split('=', 1) for a in sys.argv[6:]]
CODE = np.full(256, 4, dtype=np.int8)
for i, b in enumerate(b'ACGT'):
    CODE[b] = i
fa = pysam.FastaFile(ref_f)
vcf = pysam.VariantFile(vcf_f)
bed = collections.defaultdict(list)
for line in open(bed_f):
    c, s, e = line.split()[:3]
    bed[c].append((int(s), int(e)))
blocks = []
for line in open(blocks_f):
    c, s, e = line.split()
    s, e = int(s) - 2000, int(e) + 2000
    code = CODE[np.frombuffer(fa.fetch(c, s, e).upper().encode(), dtype=np.uint8)]
    mask = np.zeros(len(code), dtype=bool)
    for rec in vcf.fetch(c, s, e):
        a = rec.pos - 1 - s
        b = rec.pos - 1 + len(rec.ref) - s
        mask[max(a - 5, 0):max(b + 5, 0)] = True
    inbed = np.zeros(len(code), dtype=bool)
    for bs, be in bed.get(c, ()):
        if be > s and bs < e:
            inbed[max(bs - s, 0):min(be - s, len(code))] = True
    blocks.append(dict(chrom=c, start=s, code=code, mask=mask, inbed=inbed))
by_chrom = collections.defaultdict(list)
for b in blocks:
    by_chrom[b['chrom']].append(b)


def find_block(c, p):
    for b in by_chrom.get(c, ()):
        if b['start'] + 200 <= p < b['start'] + len(b['code']) - 400:
            return b
    return None


for name, bam_f in sets:
    counts = {id(b): np.zeros((len(b['code']), 2, 5), dtype=np.int32) for b in blocks}
    for r in pysam.AlignmentFile(bam_f):
        if r.is_unmapped or r.is_secondary or r.is_supplementary or r.mapping_quality < 20:
            continue
        b = find_block(r.reference_name, r.reference_start)
        if b is None:
            continue
        cig = r.cigartuples
        if any(o not in (0, 4) for o, _ in cig):
            continue
        q0 = cig[0][1] if cig[0][0] == 4 else 0
        ln = sum(l for o, l in cig if o == 0)
        rp = r.reference_start - b['start']
        seq = CODE[np.frombuffer(r.query_sequence.encode(), dtype=np.uint8)[q0:q0 + ln]]
        qual = np.asarray(r.query_qualities[q0:q0 + ln])
        ok = qual >= minbq
        idx = np.arange(rp, rp + ln)[ok]
        np.add.at(counts[id(b)], (idx, int(r.is_reverse), seq[ok]), 1)
    tot = collections.Counter()
    for b in blocks:
        c = counts[id(b)][:, :, :4].astype(np.int64)
        code = b['code']
        ok = (~b['mask']) & b['inbed'] & (code < 4)
        depth = c.sum(axis=(1, 2))
        sel = np.nonzero(ok & (depth >= 8))[0]
        cs = c[sel].copy()
        ref = code[sel].astype(np.int64)
        cs[np.arange(len(sel)), :, ref] = 0  # zero ref base on both strands
        alt_tot = cs.sum(axis=1)  # (n, 4)
        a = alt_tot.argmax(axis=1)
        n_a = alt_tot[np.arange(len(sel)), a]
        cand = (n_a >= 2) & (n_a / depth[sel] >= 0.12)
        one = cand & (np.minimum(cs[np.arange(len(sel)), 0, a], cs[np.arange(len(sel)), 1, a]) == 0)
        tot['scored'] += len(sel)
        tot['cand'] += int(cand.sum())
        tot['cand_one_strand'] += int(one.sum())
        mx = cs.max(axis=2)  # (n, 2) max same-alt per strand
        anyalt = cs.sum(axis=2)  # (n,2) mismatches per strand, any alt
        tot['ss_sa>=2'] += int((mx >= 2).sum())
        tot['ss_sa>=3'] += int((mx >= 3).sum())
        tot['ss_any>=3'] += int((anyalt >= 3).sum())
        tot['ss_any>=2'] += int((anyalt >= 2).sum())
        tot['mismatches'] += int(anyalt.sum())
        tot['bases'] += int(depth[sel].sum())
    print(name, 'minBQ', minbq, dict(tot), 'cand/Mb %.1f' % (1e6 * tot['cand'] / max(tot['scored'], 1)))
