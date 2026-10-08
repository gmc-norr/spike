"""Recurrent errors and caller-style candidates, real vs spike over the same templates (K2b).
Usage: python3 recur.py BLOCKS REF VCF BENCH_BED NAME=BAM ...
Per site (ref pos), counts bases with BQ>=10 from MAPQ>=20 primary reads with no indel in CIGAR
(aligned part only), per strand and called base. Sites within 5 bp of a Q100 record are masked;
only sites inside the Q100 v5.0q small-variant benchmark BED are scored.
"""
import sys, collections, bisect
import numpy as np
import pysam

blocks_f, ref_f, vcf_f, bed_f = sys.argv[1:5]
sets = [a.split('=', 1) for a in sys.argv[5:]]
CODE = np.full(256, 4, dtype=np.int8)
for i, b in enumerate(b'ACGT'):
    CODE[b] = i
fa = pysam.FastaFile(ref_f)
vcf = pysam.VariantFile(vcf_f)
bed = collections.defaultdict(list)
for line in open(bed_f):
    c, s, e = line.split()[:3]
    bed[c].append((int(s), int(e)))
for c in bed:
    bed[c].sort()

blocks = []
for line in open(blocks_f):
    c, s, e = line.split()
    s, e = int(s) - 2000, int(e) + 2000
    code = CODE[np.frombuffer(fa.fetch(c, s, e).upper().encode(), dtype=np.uint8)]
    mask = np.zeros(len(code), dtype=bool)
    near = np.zeros(len(code), dtype=bool)
    for rec in vcf.fetch(c, s, e):
        a = rec.pos - 1 - s
        b = rec.pos - 1 + len(rec.ref) - s
        mask[max(a - 5, 0):max(b + 5, 0)] = True
        near[max(a - 5, 0):max(b + 5, 0)] = True
    inbed = np.zeros(len(code), dtype=bool)
    for bs, be in bed.get(c, ()):
        if be > s and bs < e:
            inbed[max(bs - s, 0):min(be - s, len(code))] = True
    # longest one-letter run ending at each position (forward) and starting at each position (reverse)
    runf = np.zeros(len(code), dtype=np.int32)
    for i in range(len(code)):
        runf[i] = runf[i - 1] + 1 if i and code[i] == code[i - 1] and code[i] < 4 else 1
    runr = np.zeros(len(code), dtype=np.int32)
    for i in range(len(code) - 1, -1, -1):
        runr[i] = runr[i + 1] + 1 if i + 1 < len(code) and code[i] == code[i + 1] and code[i] < 4 else 1
    blocks.append(dict(chrom=c, start=s, code=code, mask=mask, near=near, inbed=inbed, runf=runf, runr=runr))
by_chrom = collections.defaultdict(list)
for b in blocks:
    by_chrom[b['chrom']].append(b)


def find_block(c, p):
    for k, b in enumerate(by_chrom.get(c, ())):
        if b['start'] + 200 <= p < b['start'] + len(b['code']) - 400:
            return b
    return None


def analyse(bam_f):
    # counts[block_id] : array (len, strand 2, base 5)
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
        ok = qual >= 10
        idx = np.arange(rp, rp + ln)[ok]
        np.add.at(counts[id(b)], (idx, int(r.is_reverse), seq[ok]), 1)
    return counts


res = {name: analyse(f) for name, f in sets}
for name, counts in res.items():
    cand = 0; cand_one_strand = 0; scored = 0; sites2 = 0; sites3 = 0; opp = 0
    learnmask_hq = [0, 0]
    after_run = collections.Counter(); base_after = collections.Counter()
    for b in blocks:
        c = counts[id(b)]
        code = b['code']
        ok = (~b['mask']) & b['inbed'] & (code < 4)
        depth = c[:, :, :4].sum(axis=(1, 2))
        for p in np.nonzero(ok & (depth >= 8))[0]:
            scored += 1
            ref = code[p]
            alt_tot = c[p, :, :4].sum(axis=0)
            alt_tot[ref] = 0
            a = int(alt_tot.argmax()); n = int(alt_tot[a])
            if n >= 2 and n / depth[p] >= 0.12:
                cand += 1
                if min(c[p, 0, a], c[p, 1, a]) == 0:
                    cand_one_strand += 1
            for s in range(2):
                ss = c[p, s, :4].copy(); ss[ref] = 0
                aa = int(ss.argmax())
                if ss[aa] >= 3:
                    sites3 += 1
                    if c[p, 1 - s, aa] > 0:
                        opp += 1
                    # longest run just before p in sequencing orientation, within 15 bp
                    if s == 0:
                        lo = max(p - 15, 0)
                        rl = int(b['runf'][lo:p].max()) if p > lo else 0
                    else:
                        hi = min(p + 16, len(code))
                        rl = int(b['runr'][p + 1:hi].max()) if hi > p + 1 else 0
                    after_run['>=7' if rl >= 7 else ('5-6' if rl >= 5 else '<5')] += 1
                if ss[aa] >= 2:
                    sites2 += 1
    print('%s: scored sites (depth>=8, in bench, unmasked) %d; candidates (alt>=2 & AF>=0.12) %d (%.1f per Mb), one-strand %d; '
          'same-strand same-alt >=2 reads %d, >=3 reads %d (opposite strand also has it: %d); >=3 sites by longest run within 15 bp before (seq orientation): %s'
          % (name, scored, cand, 1e6 * cand / max(scored, 1), cand_one_strand, sites2, sites3, opp, dict(after_run)))

# background: run-length distribution before random scored positions
tot = collections.Counter()
for b in blocks:
    code = b['code']
    ok = np.nonzero((~b['mask']) & b['inbed'] & (code < 4))[0]
    for p in ok[::50]:
        lo = max(p - 15, 0)
        rl = int(b['runf'][lo:p].max()) if p > lo else 0
        tot['>=7' if rl >= 7 else ('5-6' if rl >= 5 else '<5')] += 1
print('background positions by longest run within 15 bp before (forward):', dict(tot))
