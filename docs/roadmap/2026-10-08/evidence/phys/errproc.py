"""Read-level error process, real vs spike over the same templates (K2b sets).
Usage: python3 -I errproc.py BLOCKS REF VCF NAME=BAM [NAME=BAM ...]
Counts in sequencing orientation, MAPQ>=20 primary reads, Q100 variant sites masked (+-5 bp).
"""
import sys, collections, bisect
import numpy as np
import pysam

blocks_f, ref_f, vcf_f = sys.argv[1:4]
sets = [a.split('=', 1) for a in sys.argv[4:]]
CODE = np.full(256, 4, dtype=np.int8)
for i, b in enumerate(b'ACGT'):
    CODE[b] = i
COMP = np.array([3, 2, 1, 0, 4], dtype=np.int8)
fa = pysam.FastaFile(ref_f)
vcf = pysam.VariantFile(vcf_f)

blocks = []
for line in open(blocks_f):
    c, s, e = line.split()
    s, e = int(s) - 2000, int(e) + 2000
    ref = fa.fetch(c, s, e).upper().encode()
    code = CODE[np.frombuffer(ref, dtype=np.uint8)]
    mask = np.zeros(len(ref), dtype=bool)
    for rec in vcf.fetch(c, s, e):
        a = rec.pos - 1 - 5 - s
        b = rec.pos - 1 + len(rec.ref) + 5 - s
        mask[max(a, 0):max(b, 0)] = True
    # homopolymer runs >= 4
    runs = []
    i = 0
    n = len(code)
    while i < n:
        j = i
        while j + 1 < n and code[j + 1] == code[i] and code[i] < 4:
            j += 1
        L = j - i + 1
        if L >= 4 and code[i] < 4 and not mask[max(i - 1, 0):j + 2].any():
            runs.append((i, j + 1, L, int(code[i])))
        i = j + 1
    runs.sort()
    blocks.append(dict(chrom=c, start=s, code=code, mask=mask, runs=runs, rstarts=[r[0] for r in runs]))
by_chrom = collections.defaultdict(list)
for b in blocks:
    by_chrom[b['chrom']].append(b)


def find_block(c, p):
    for b in by_chrom.get(c, ()):
        if b['start'] + 200 <= p < b['start'] + len(b['code']) - 400:
            return b
    return None


def hp_bin(L):
    for lo, name in [(15, '15+'), (12, '12-14'), (10, '10-11'), (8, '8-9'), (6, '6-7'), (4, '4-5')]:
        if L >= lo:
            return name
    return '<4'


BINS = ['4-5', '6-7', '8-9', '10-11', '12-14', '15+']


def analyse(bam_f):
    tot = np.zeros((2, 64), dtype=np.int64)     # mate x phred
    mm = np.zeros((2, 64), dtype=np.int64)
    spec = np.zeros((2, 3, 4, 4), dtype=np.int64)  # mate, qbin, true, called
    cyc_tot = np.zeros((2, 160), dtype=np.int64)
    cyc_mm = np.zeros((2, 160), dtype=np.int64)
    ctx_tot = np.zeros((2, 16), dtype=np.int64)  # high-Q (>=30) bases by preceding dinucleotide, seq orientation; [0]=all, [1]=mm
    recur = collections.Counter()
    q2_nonN = 0; q2_all = 0; n_bases = 0; n_N = 0
    hp_pass = collections.Counter(); hp_ind = collections.Counter()
    nonhp_bases = 0; nonhp_ind = 0
    reads = 0
    for r in pysam.AlignmentFile(bam_f):
        if r.is_unmapped or r.is_secondary or r.is_supplementary or r.mapping_quality < 20:
            continue
        b = find_block(r.reference_name, r.reference_start)
        if b is None:
            continue
        reads += 1
        mate = 0 if r.is_read1 else 1
        rev = r.is_reverse
        seq = np.frombuffer(r.query_sequence.encode(), dtype=np.uint8)
        qual = np.asarray(r.query_qualities, dtype=np.int64)
        L = len(seq)
        q2_all += int((qual == 2).sum()); q2_nonN += int(((qual == 2) & (seq != ord('N'))).sum())
        n_N += int((seq == ord('N')).sum()); n_bases += L
        off = b['start']
        code = b['code']; mask = b['mask']
        cig = r.cigartuples
        ops = [o for o, _ in cig]
        rpos = r.reference_start - off
        qpos = 0
        aligned_lo = rpos
        # indel events
        for k, (op, ln) in enumerate(cig):
            if op in (0, 7, 8):
                if all(o in (0, 4, 7, 8) for o in ops):  # substitutions from no-indel reads only
                    rs = code[rpos:rpos + ln]; qs = CODE[seq[qpos:qpos + ln]]; qq = qual[qpos:qpos + ln]
                    mk = mask[rpos:rpos + ln]
                    ok = (~mk) & (rs < 4) & (qs < 4)
                    cyc = np.arange(qpos, qpos + ln)
                    if rev:
                        cyc = L - 1 - cyc
                        t = COMP[rs]; c_ = COMP[qs]
                    else:
                        t = rs; c_ = qs
                    is_mm = ok & (rs != qs)
                    np.add.at(tot[mate], qq[ok], 1); np.add.at(mm[mate], qq[is_mm], 1)
                    np.add.at(cyc_tot[mate], cyc[ok], 1); np.add.at(cyc_mm[mate], cyc[is_mm], 1)
                    qb = np.where(qq < 15, 0, np.where(qq < 30, 1, 2))
                    idx = np.nonzero(is_mm)[0]
                    for i in idx:
                        spec[mate, qb[i], t[i], c_[i]] += 1
                        recur[(b['chrom'], off + rpos + int(i), rev, int(c_[i]))] += 1
                    # preceding dinucleotide in sequencing orientation, high-Q bases
                    p = np.arange(rpos, rpos + ln)
                    if rev:
                        a1 = code[np.minimum(p + 2, len(code) - 1)]; a2 = code[np.minimum(p + 1, len(code) - 1)]
                        a1 = COMP[a1]; a2 = COMP[a2]
                    else:
                        a1 = code[np.maximum(p - 2, 0)]; a2 = code[np.maximum(p - 1, 0)]
                    hq = ok & (qq >= 30) & (a1 < 4) & (a2 < 4)
                    cx = (a1.astype(np.int64) * 4 + a2)
                    np.add.at(ctx_tot[0], cx[hq], 1); np.add.at(ctx_tot[1], cx[hq & is_mm], 1)
                rpos += ln; qpos += ln
            elif op == 1:
                if 0 < k < len(cig) - 1 and ln <= 2 and not mask[max(rpos - 3, 0):rpos + 3].any():
                    # run touching the insertion point
                    j = bisect.bisect_right(b['rstarts'], rpos) - 1
                    hit = None
                    for jj in (j, j + 1):
                        if 0 <= jj < len(b['runs']):
                            s0, e0, Lr, base = b['runs'][jj]
                            if s0 <= rpos <= e0:
                                hit = Lr
                    if hit:
                        hp_ind[(hp_bin(hit), 'I')] += 1
                    else:
                        nonhp_ind += 1
                qpos += ln
            elif op == 2:
                if 0 < k < len(cig) - 1 and ln <= 2 and not mask[max(rpos - 3, 0):rpos + ln + 3].any():
                    j = bisect.bisect_right(b['rstarts'], rpos) - 1
                    hit = None
                    if 0 <= j < len(b['runs']):
                        s0, e0, Lr, base = b['runs'][j]
                        if s0 <= rpos < e0:
                            hit = Lr
                    if hit:
                        hp_ind[(hp_bin(hit), 'D')] += 1
                    else:
                        nonhp_ind += 1
                rpos += ln
            elif op == 4:
                qpos += ln
        aligned_hi = rpos
        # homopolymer passes: runs fully inside the aligned span with 3 bp margin
        lo = bisect.bisect_left(b['rstarts'], aligned_lo + 3)
        inside = 0
        for jj in range(lo, len(b['runs'])):
            s0, e0, Lr, base = b['runs'][jj]
            if s0 >= aligned_hi - 3:
                break
            if e0 <= aligned_hi - 3:
                hp_pass[hp_bin(Lr)] += 1
                inside += Lr
        nonhp_bases += (aligned_hi - aligned_lo) - inside
    return dict(tot=tot, mm=mm, spec=spec, cyc_tot=cyc_tot, cyc_mm=cyc_mm, ctx=ctx_tot, recur=recur,
                q2_all=q2_all, q2_nonN=q2_nonN, n_N=n_N, n_bases=n_bases, hp_pass=hp_pass, hp_ind=hp_ind,
                nonhp_bases=nonhp_bases, nonhp_ind=nonhp_ind, reads=reads)


res = {name: analyse(f) for name, f in sets}
B = 'ACGT'
for name, d in res.items():
    print('=' * 20, name, 'reads', d['reads'])
    print('N bases %d (%.2e/base); Q2 bases %d, Q2 on non-N %d' % (d['n_N'], d['n_N'] / d['n_bases'], d['q2_all'], d['q2_nonN']))
    for m in range(2):
        row = []
        for q in range(64):
            if d['tot'][m, q] > 1000:
                row.append('Q%d %.4f%% (n=%d)' % (q, 100 * d['mm'][m, q] / d['tot'][m, q], d['tot'][m, q]))
        print(' R%d mismatch by Q:' % (m + 1), '; '.join(row))
    for qb, qn in enumerate(['Q<15', 'Q15-29', 'Q>=30']):
        for m in range(2):
            s = d['spec'][m, qb]
            parts = []
            for t in range(4):
                tt = s[t].sum()
                if tt:
                    parts.append('%s->' % B[t] + '/'.join('%s%.2f' % (B[c], s[t, c] / tt) for c in range(4) if c != t) + '(%d)' % tt)
            mx = sum(s[t].max() for t in range(4)) / max(s.sum(), 1)
            print(' %s R%d top-alt share %.3f | %s' % (qn, m + 1, mx, ' '.join(parts)))
    for m in range(2):
        ct, cm = d['cyc_tot'][m], d['cyc_mm'][m]
        def rate(a, b_):
            return 100 * cm[a:b_].sum() / max(ct[a:b_].sum(), 1)
        print(' R%d mismatch%% cycles 1-5 %.3f, 6-10 %.3f, 50-100 %.3f, 141-145 %.3f, 146-150 %.3f, 151 %.3f' % (
            m + 1, rate(0, 5), rate(5, 10), rate(49, 100), rate(140, 145), rate(145, 150), rate(150, 151)))
    c0, c1 = d['ctx']
    overall = c1.sum() / c0.sum()
    parts = sorted(((c1[i] / c0[i] / overall, B[i // 4] + B[i % 4], c1[i], c0[i]) for i in range(16) if c0[i]), reverse=True)
    print(' high-Q (>=30) mismatch rate by preceding dinucleotide (seq orientation), ratio to overall %.2e:' % overall)
    print('   ' + ' '.join('%s %.2fx(%d)' % (k, v, e) for v, k, e, t in parts))
    rc = d['recur']
    nmm = sum(rc.values())
    k2 = sum(1 for v in rc.values() if v >= 2); k3 = sum(1 for v in rc.values() if v >= 3)
    e2 = sum(v for v in rc.values() if v >= 2)
    print(' recurrent same-strand same-alt error sites: keys %d, >=2 reads %d, >=3 reads %d; mismatches in >=2 sites %d of %d (%.2f%%)' % (
        len(rc), k2, k3, e2, nmm, 100 * e2 / max(nmm, 1)))
    print(' slips (1-2 bp indels) per read pass over a homopolymer run, by run length:')
    for bn in BINS:
        p = d['hp_pass'][bn]; i_ = d['hp_ind'][(bn, 'I')]; dd = d['hp_ind'][(bn, 'D')]
        print('   %6s passes %8d  ins %5d  del %5d  rate %.4f%%' % (bn, p, i_, dd, 100 * (i_ + dd) / max(p, 1)))
    print('   outside runs>=4: %d indels over %d bases (%.2e/base)' % (d['nonhp_ind'], d['nonhp_bases'], d['nonhp_ind'] / max(d['nonhp_bases'], 1)))
