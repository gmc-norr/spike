"""How much of the real-vs-spike recurrent-site gap is only the substitution matrix (which wrong base)?
Learn P(called | mate, strand, qbin, ref) from REAL mismatches; re-draw the called base of every spike
mismatch from that matrix (position and count of errors unchanged); recount recur.py's numbers.
Usage: python3 -I matrix_inject.py BLOCKS REF VCF BED REAL_BAM SPIKE_BAM SEED [SEED ...]
"""
import sys, collections
import numpy as np
import pysam

blocks_f, ref_f, vcf_f, bed_f, real_f, spike_f = sys.argv[1:7]
seeds = [int(x) for x in sys.argv[7:]]
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
for k, b in enumerate(blocks):
    b['k'] = k
    by_chrom[b['chrom']].append(b)


def find_block(c, p):
    for b in by_chrom.get(c, ()):
        if b['start'] + 200 <= p < b['start'] + len(b['code']) - 400:
            return b
    return None


def load(bam_f):
    """Return list of per-read arrays: (block k, positions, strand, mate, qbin, called code)."""
    out = []
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
        ok = qual >= 10
        pos = np.arange(rp, rp + ln)[ok]
        qq = qual[ok]
        qb = np.where(qq < 15, 0, np.where(qq < 30, 1, 2))
        out.append((b['k'], pos, int(r.is_reverse), 0 if r.is_read1 else 1, qb, seq[ok].copy()))
    return out


def count(reads):
    counts = [np.zeros((len(b['code']), 2, 5), dtype=np.int32) for b in blocks]
    for k, pos, st, m, qb, called in reads:
        np.add.at(counts[k], (pos, st, called), 1)
    tot = collections.Counter()
    for b in blocks:
        c = counts[b['k']][:, :, :4].astype(np.int64)
        code = b['code']
        ok = (~b['mask']) & b['inbed'] & (code < 4)
        depth = c.sum(axis=(1, 2))
        sel = np.nonzero(ok & (depth >= 8))[0]
        cs = c[sel].copy()
        ref = code[sel].astype(np.int64)
        cs[np.arange(len(sel)), :, ref] = 0
        alt_tot = cs.sum(axis=1)
        a = alt_tot.argmax(axis=1)
        n_a = alt_tot[np.arange(len(sel)), a]
        cand = (n_a >= 2) & (n_a / depth[sel] >= 0.12)
        tot['cand'] += int(cand.sum())
        mx = cs.max(axis=2)
        tot['ss_sa>=2'] += int((mx >= 2).sum())
        tot['ss_sa>=3'] += int((mx >= 3).sum())
        tot['ss_any>=3'] += int((cs.sum(axis=2) >= 3).sum())
    return dict(tot)


real = load(real_f)
spike = load(spike_f)
# matrix from real mismatches: cell (mate, strand, qbin, ref) -> counts of called
M = np.zeros((2, 2, 3, 4, 4), dtype=np.float64)
for k, pos, st, m, qb, called in real:
    ref = blocks[k]['code'][pos]
    mm = (ref < 4) & (called < 4) & (ref != called) & ~blocks[k]['mask'][pos]
    np.add.at(M, (m, st, qb[mm], ref[mm], called[mm]), 1)
print('real', count(real))
print('spike', count(spike))
for seed in seeds:
    rng = np.random.default_rng(seed)
    new = []
    for k, pos, st, m, qb, called in spike:
        ref = blocks[k]['code'][pos]
        mm = np.nonzero((ref < 4) & (called < 4) & (ref != called))[0]
        if len(mm):
            called = called.copy()
            for i in mm:
                w = M[m, st, qb[i], ref[i]].copy()
                w[ref[i]] = 0
                if w.sum() == 0:
                    continue
                called[i] = rng.choice(4, p=w / w.sum())
        new.append((k, pos, st, m, qb, called))
    print('spike+real-matrix seed', seed, count(new))
