"""Does spike's learning rule count the sample's own variants as high-quality errors?
Counts mismatches by quality on real reads the way quality::block_errors does (aligned bases plus
clipped ends judged bad ends: 1-4 bp, or >=50% reference), with spike's variant mask (>=5 aligned
reads and >=10% differ), then asks how many of the counted high-Q errors sit within 5 bp of a Q100 record.
Usage: python3 learnmask.py BLOCKS REF VCF BAM
"""
import sys, collections
import numpy as np
import pysam

blocks_f, ref_f, vcf_f, bam_f = sys.argv[1:5]
CODE = np.full(256, 4, dtype=np.int8)
for i, b in enumerate(b'ACGT'):
    CODE[b] = i
fa = pysam.FastaFile(ref_f)
vcf = pysam.VariantFile(vcf_f)
blocks = []
for line in open(blocks_f):
    c, s, e = line.split()
    s, e = int(s) - 2000, int(e) + 2000
    code = CODE[np.frombuffer(fa.fetch(c, s, e).upper().encode(), dtype=np.uint8)]
    near = np.zeros(len(code), dtype=bool)
    for rec in vcf.fetch(c, s, e):
        a = rec.pos - 1 - s; b = rec.pos - 1 + len(rec.ref) - s
        near[max(a - 5, 0):max(b + 5, 0)] = True
    blocks.append(dict(chrom=c, start=s, code=code, near=near, pile=np.zeros((len(code), 2), dtype=np.int32), reads=[]))


def find_block(c, p):
    for b in blocks:
        if b['chrom'] == c and b['start'] + 200 <= p < b['start'] + len(b['code']) - 400:
            return b
    return None


def placed(r, b):
    """(refidx, base, qual, kind) per query base; kind A aligned, C clipped, I inserted."""
    out = []
    rp = r.reference_start - b['start']
    seq = r.query_sequence; qual = r.query_qualities
    cig = r.cigartuples
    q = 0
    lead = cig[0][1] if cig[0][0] == 4 else 0
    for k in range(lead):
        out.append((rp - lead + k, seq[q], qual[q], 'C')); q += 1
    for op, ln in cig:
        if op == 4:
            continue
        if op in (0, 7, 8):
            for k in range(ln):
                out.append((rp, seq[q], qual[q], 'A')); rp += 1; q += 1
        elif op == 1:
            for k in range(ln):
                out.append((None, seq[q], qual[q], 'I')); q += 1
        elif op == 2:
            rp += ln
    while q < len(seq):
        out.append((rp, seq[q], qual[q], 'C')); rp += 1; q += 1
    return out


for r in pysam.AlignmentFile(bam_f):
    if r.is_unmapped or r.is_secondary or r.is_supplementary or r.mapping_quality < 20 or not r.is_proper_pair:
        continue
    b = find_block(r.reference_name, r.reference_start)
    if b is None:
        continue
    pl = placed(r, b)
    b['reads'].append(pl)
    for (i, base, q, kind) in pl:
        if kind == 'A' and i is not None and 0 <= i < len(b['code']) and base != 'N' and b['code'][i] < 4:
            b['pile'][i, 0] += 1
            b['pile'][i, 1] += int(CODE[ord(base)] != b['code'][i])

cnt = collections.Counter()
for b in blocks:
    code = b['code']; pile = b['pile']
    var = (pile[:, 0] >= 5) & (pile[:, 1] * 10 >= pile[:, 0])
    for pl in b['reads']:
        lead = 0
        while lead < len(pl) and pl[lead][3] == 'C':
            lead += 1
        trail = 0
        while trail < len(pl) and pl[len(pl) - 1 - trail][3] == 'C':
            trail += 1

        def clip_ok(rng):
            if len(rng) <= 4:
                return True
            same = comp = 0
            for k in rng:
                i, base, q, kind = pl[k]
                if 0 <= i < len(code) and code[i] < 4 and base != 'N':
                    comp += 1; same += int(CODE[ord(base)] == code[i])
            return comp > 0 and 2 * same >= comp
        lok = clip_ok(range(0, lead)); tok = clip_ok(range(len(pl) - trail, len(pl)))
        for k, (i, base, q, kind) in enumerate(pl):
            if kind == 'I' or i is None or not (0 <= i < len(code)) or code[i] == 4 or base == 'N':
                continue
            if kind == 'C' and not (lok if k < lead else tok):
                continue
            hq = 'hi' if q >= 30 else ('mid' if q >= 15 else 'lo')
            err = int(CODE[ord(base)] != code[i])
            if not var[i]:
                cnt[('spike-rule', hq, 'bases')] += 1
                cnt[('spike-rule', hq, 'err')] += err
                if err and b['near'][i]:
                    cnt[('spike-rule', hq, 'err_near_q100')] += 1
                if err and kind == "C":
                    cnt[(rule_clip := "spike-rule", hq, "clip_err_near")] += int(b["near"][i])
                    cnt[("cliplen", hq, (lead if k < lead else trail))] += 1
                    cnt[('spike-rule', hq, 'err_in_clip')] += 1
            if not b['near'][i] and kind == 'A':
                cnt[('truth-mask, aligned only', hq, 'bases')] += 1
                cnt[('truth-mask, aligned only', hq, 'err')] += err
for rule in ['spike-rule', 'truth-mask, aligned only']:
    for hq in ['lo', 'mid', 'hi']:
        n = cnt[(rule, hq, 'bases')]; e = cnt[(rule, hq, 'err')]
        extra = ''
        if rule == 'spike-rule':
            extra = '; errors within 5 bp of a Q100 record %d (%.1f%%), in counted clips %d' % (
                cnt[(rule, hq, 'err_near_q100')], 100 * cnt[(rule, hq, 'err_near_q100')] / max(e, 1), cnt[(rule, hq, 'err_in_clip')])
        print('%-26s %-3s bases %9d errors %7d rate %.3e%s' % (rule, hq, n, e, e / max(n, 1), extra))
for hq in ['lo', 'mid', 'hi']:
    print(hq, 'clip errors within 5 bp of a Q100 record:', cnt[('spike-rule', hq, 'clip_err_near')], 'by clip length:', sorted((k[2], v) for k, v in cnt.items() if k[0] == 'cliplen' and k[1] == hq)[:12])
