"""Re-score the existing hospital K2 data (qtq2/k2v2/*/spiked.bam) with k2runs.py's read set and
bad-end rule, split by: run bin, end of the bad clip in sequencing orientation (5'/3'), clip length,
read length (151 or shorter), |TLEN| >= 151, and whether a Q100 variant lies within 10 bp of the
boundary. Read-only on existing files."""
import sys, math
from collections import Counter, defaultdict
import pysam
REF = pysam.FastaFile("/home/parlar_ai/dev/spike/data/giab_hg38/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta")
Q100 = pysam.VariantFile("/home/parlar_ai/dev/spike/data/giab_hg38/HG002/GRCh38_HG2-T2TQ100-V1.1.vcf.gz")
root, events = sys.argv[1], sys.argv[2]
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
def clips(r):
    out = []; cig = r.cigartuples; seq = r.query_sequence
    if cig[0][0] == 4: L = cig[0][1]; out.append(("L", r.reference_start, seq[:L], r.reference_start - L))
    if cig[-1][0] == 4: L = cig[-1][1]; out.append(("R", r.reference_end, seq[-L:], r.reference_end))
    return out
def near_q100(chrom, pos, w=10):
    return any(True for _ in Q100.fetch(chrom, max(0, pos - w - 1), pos + w))
rows = [l.rstrip('\n').split('\t') for l in open(events)][1:]
tab = defaultdict(Counter)  # (sp, stratum) -> counts
for ident, kind, chrom, s, e, cls in rows:
    s, e = int(s), int(e)
    try: b = pysam.AlignmentFile(f'{root}/{ident}/spiked.bam')
    except (FileNotFoundError, ValueError): continue
    bound = Counter()
    for r in b.fetch(chrom, max(0, s - 20000), e + 20000):
        if r.is_unmapped or r.is_secondary or r.is_supplementary or r.is_duplicate: continue
        for side, pos, _, _ in clips(r): bound[(side, pos)] += 1
    for lo, hi in [(s - 2000, s - 300), (e + 300, e + 2000)]:
        for r in b.fetch(chrom, lo, hi):
            if r.is_secondary or r.is_supplementary or r.is_duplicate or r.is_unmapped or r.mapping_quality < 20: continue
            if not (lo <= r.reference_start and r.reference_end <= hi): continue
            sp = r.query_name.startswith('SPIKE_')
            h = template_bin(r, chrom)
            full = r.query_length == 151
            longt = abs(r.template_length) >= 151 and r.is_proper_pair
            k2ok = full and longt
            lc = 'L151' if r.query_length == 151 else 'L150' if r.query_length == 150 else 'L<150'
            pp = 'pp' if (r.is_proper_pair and abs(r.template_length) >= 151) else 'npp'
            keys = [f'bin{min(h,1)}_{lc}', f'bin{min(h,1)}_{lc}_{pp}', 'all_'+lc] if True else ['all', f'bin{h}', f'bin{h}_' + ('K2bset' if k2ok else 'notK2bset'), 'len151' if full else 'lenshort',
                    'tlen>=151' if longt else 'tlen<151']
            bad_ends = []
            for side, pos, bases, pstart in clips(r):
                others = sum(bound[(side, pos + d)] for d in range(-2, 3)) - 1
                if others >= 3: continue
                L = len(bases)
                isbad = L <= 4
                if not isbad:
                    ref = REF.fetch(chrom, max(0, pstart), pstart + L).upper()
                    cmp = [(x, y) for x, y in zip(bases.upper(), ref) if x != 'N' and y != 'N']
                    isbad = bool(cmp) and sum(x == y for x, y in cmp) * 2 >= len(cmp)
                if isbad:
                    end = ('3p' if side == 'L' else '5p') if r.is_reverse else ('5p' if side == 'L' else '3p')
                    bad_ends.append((end, '1-4' if L <= 4 else '5-20' if L <= 20 else '21+', near_q100(chrom, pos)))
            for k in keys:
                t = tab[(sp, k)]; t['reads'] += 1; t['bad'] += bool(bad_ends)
                t['bad_noq100'] += any(not q for _, _, q in bad_ends)
                for end, lb, q in bad_ends:
                    t[f'{end}_{lb}'] += 1
                    if not q: t[f'{end}_{lb}_noq100'] += 1
def z(k1, n1, k2, n2):
    p = (k1 + k2) / (n1 + n2); se = math.sqrt(p * (1 - p) * (1 / n1 + 1 / n2)) if 0 < p < 1 else 0; return (k1 / n1 - k2 / n2) / se if se else 0.0
strata = sorted({k for _, k in tab})
cols = ['bad', 'bad_noq100'] + [f'{e}_{l}' for e in ('5p', '3p') for l in ('1-4', '5-20', '21+')]
print('stratum              set    reads  ' + '  '.join(f'{c:>11s}' for c in cols) + '   z(bad)  z(bad_noq100)')
for k in strata:
    S, O = tab[(True, k)], tab[(False, k)]
    for name, t in (('spike', S), ('own', O)):
        n = max(1, t['reads'])
        extra = ''
        if name == 'spike' and O['reads'] and S['reads']:
            extra = f"   {z(S['bad'], S['reads'], O['bad'], O['reads']):+.2f}   {z(S['bad_noq100'], S['reads'], O['bad_noq100'], O['reads']):+.2f}"
        print(f'{k:20s} {name:5s} {t["reads"]:7d}  ' + '  '.join(f'{100*t[c]/n:6.2f}%({t[c]:3d})' for c in cols) + extra)
