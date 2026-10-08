"""Reported parts of K2 (v2 plan): the bad-end share by the template's longest-run bin, on plant.sh outputs
(k2.py's reads and bad-end rule).

In each event's flanks ([s-2000, s-300] and [e+300, e+2000]), primary non-duplicate reads at
MAPQ >= 20 lying wholly inside: spike's reads against the sample's own.
A clip is bad-end when it holds 1-4 bases, or >= 5 bases matching the reference at >= 50% of
placed positions, and fewer than 3 other reads in the event's window clip at the same boundary
(+-2 bp). K2 passes when the share of reads with a bad-end clip has two-proportion |z| < 3.
Reported: any soft clip, crash share, mismatches per 100 aligned bases, base dependence (K5),
and among spike pairs P(both crash) / (P(R1) P(R2)) (K3).
Usage: k2.py DIR EVENTS_TSV"""
import math, sys, statistics as st
from collections import Counter, defaultdict
import pysam
REF = pysam.FastaFile("/home/parlar_ai/dev/spike/data/giab_hg38/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta")
def crashes(q): return len(q) >= 20 and sum(x < 15 for x in q[-20:]) >= 10
def seq_order(r):
    q = list(r.query_qualities); s = r.query_sequence
    if r.is_reverse: q = q[::-1]; s = s.translate(str.maketrans('ACGTN', 'TGCAN'))[::-1]
    return s, q
def z(k1, n1, k2, n2):
    p = (k1 + k2) / (n1 + n2); se = math.sqrt(p * (1 - p) * (1 / n1 + 1 / n2)); return (k1 / n1 - k2 / n2) / se if se else 0.0
def clips(r):
    """(side, boundary, clipped bases, placed reference start) for each soft-clipped end."""
    out = []; cig = r.cigartuples; seq = r.query_sequence
    if cig[0][0] == 4: L = cig[0][1]; out.append(("L", r.reference_start, seq[:L], r.reference_start - L))
    if cig[-1][0] == 4: L = cig[-1][1]; out.append(("R", r.reference_end, seq[-L:], r.reference_end))
    return out
root, events = sys.argv[1], sys.argv[2]
rows = [l.rstrip('\n').split('\t') for l in open(events)][1:]
c = {True: Counter(), False: Counter()}; nm = {True: [], False: []}
base_n = {True: Counter(), False: Counter()}; base_low = {True: Counter(), False: Counter()}
mates = defaultdict(dict); n_events = 0
comp_t = str.maketrans('ACGTN', 'TGCAN')
def hp_bin(b, n):
    if b == 'C': return 3 if n >= 7 else 2 if n >= 5 else 0
    return 3 if n >= 12 else 2 if n >= 9 else 1 if n >= 7 else 0
def read_hp(t):
    best, prev, n = 0, '', 0
    for b in t:
        n = n + 1 if b == prev and b != 'N' else (b != 'N'); prev = b; best = max(best, hp_bin(b, n))
    return best
def template_bin(r, chrom):
    cig = r.cigartuples
    lead = cig[0][1] if cig[0][0] == 4 else 0; trail = cig[-1][1] if cig[-1][0] == 4 else 0
    t = REF.fetch(chrom, max(0, r.reference_start - lead), r.reference_end + trail).upper()
    return read_hp(t.translate(comp_t)[::-1] if r.is_reverse else t)
byrun = {True: defaultdict(lambda: [0, 0]), False: defaultdict(lambda: [0, 0])}
for ident, kind, chrom, s, e, cls in rows:
    s, e = int(s), int(e)
    try: b = pysam.AlignmentFile(f'{root}/{ident}/spiked.bam')
    except (FileNotFoundError, ValueError): print(f"  {ident}: no spiked.bam"); continue
    n_events += 1
    wlo, whi = s - 20000, e + 20000
    bound = Counter()
    for r in b.fetch(chrom, max(0, wlo), whi):
        if r.is_unmapped or r.is_secondary or r.is_supplementary or r.is_duplicate: continue
        for side, pos, _, _ in clips(r): bound[(side, pos)] += 1
    for lo, hi in [(s - 2000, s - 300), (e + 300, e + 2000)]:
        for r in b.fetch(chrom, lo, hi):
            if r.is_secondary or r.is_supplementary or r.is_duplicate or r.is_unmapped or r.mapping_quality < 20: continue
            if not (lo <= r.reference_start and r.reference_end <= hi): continue
            sp = r.query_name.startswith('SPIKE_')
            seq, q = seq_order(r)
            c[sp]['reads'] += 1
            cl = clips(r)
            c[sp]['clip'] += bool(cl)
            bad = False
            for side, pos, bases, pstart in cl:
                others = sum(bound[(side, pos + d)] for d in range(-2, 3)) - 1
                if others >= 3: continue
                if len(bases) <= 4: bad = True; continue
                ref = REF.fetch(chrom, max(0, pstart), pstart + len(bases)).upper()
                comp = [(x, y) for x, y in zip(bases.upper(), ref) if x != 'N' and y != 'N']
                if comp and sum(x == y for x, y in comp) * 2 >= len(comp): bad = True
            c[sp]['badend'] += bad
            h = template_bin(r, chrom); byrun[sp][h][0] += 1; byrun[sp][h][1] += bad
            c[sp]['crash'] += crashes(q)
            nm[sp].append(r.get_tag('NM') / r.query_alignment_length * 100)
            for base, x in zip(seq, q): base_n[sp][base] += 1; base_low[sp][base] += x < 15
    for r in b.fetch(chrom, max(0, s - 3000), e + 3000):
        if r.query_name.startswith('SPIKE_') and not r.is_secondary and not r.is_supplementary:
            mates[(ident, r.query_name)][1 if r.is_read1 else 2] = crashes(seq_order(r)[1])
sp, og = c[True], c[False]
print(f"events {n_events}; reads: spike {sp['reads']}, sample's own {og['reads']}")
zz = z(sp['badend'], sp['reads'], og['badend'], og['reads'])
print(f"K2 bad-end clip: spike {sp['badend']}/{sp['reads']} = {100*sp['badend']/sp['reads']:.2f}%, own {og['badend']}/{og['reads']} = {100*og['badend']/og['reads']:.2f}%, z = {zz:+.2f} -> {'PASS' if abs(zz) < 3 else 'FAIL'}")
for key, label in (('clip', 'any soft clip'), ('crash', 'crash')):
    print(f"reported {label}: spike {100*sp[key]/sp['reads']:.2f}%, own {100*og[key]/og['reads']:.2f}%, z = {z(sp[key], sp['reads'], og[key], og['reads']):+.2f}")
print(f"reported mismatches per 100 aligned bases: spike {st.mean(nm[True]):.3f}, own {st.mean(nm[False]):.3f}")
print('K5 share of bases < Q15 by base: ' + '; '.join(f"{k}: spike {base_low[True][k]/base_n[True][k]:.4f} own {base_low[False][k]/base_n[False][k]:.4f}" for k in 'ACGT'))
both = [(d[1], d[2]) for d in mates.values() if 1 in d and 2 in d]
n = len(both); c1 = sum(a for a, _ in both) / n; c2 = sum(b_ for _, b_ in both) / n; c12 = sum(a and b_ for a, b_ in both) / n
print(f"K3 spike pairs {n}: P(R1 crash) {c1:.4f} P(R2 crash) {c2:.4f} P(both) {c12:.5f} ratio {(c12 / (c1 * c2)) if c1 * c2 else float('nan'):.2f}")
for sp_, lab in ((True, 'spike'), (False, "own")):
    print(f"bad-end by run bin, {lab}: " + "  ".join(f"{h}: {100*byrun[sp_][h][1]/max(1,byrun[sp_][h][0]):.2f}% ({byrun[sp_][h][0]})" for h in range(4)))
