"""Kill test K of docs/superpowers/plans/2026-10-04-carried-allele.md.

Applies the carried-allele rule to a BAM at sites the sample carries (K1: its
own calls) and sites it does not (K2: random, and another sample's variants),
and reports how often the rule refuses each group.

Usage: probe.py BAM REFERENCE CALLS_VCF OTHER_VCF OTHER_BED OUT_DIR
"""
import random
import sys
from collections import Counter, defaultdict

import pysam

CHROM = 'chr20'
MIN_MAPQ = 20
MIN_READS = 10
MIN_SHARE = 0.2
SPACING = 50
CALL_MARGIN = 100

# Pass bounds on checked sites (locked in the plan).
K1_MIN = {'het_snv': 0.95, 'hom_snv': 0.99, 'het_indel': 0.90, 'hom_indel': 0.95}
K2_MAX = 0.01


def trim(pos, ref, alt):
    """The changed reference bases [s, e) (0-based) of REF/ALT at 0-based pos:
    the common prefix first, then the common suffix."""
    i = 0
    while i < min(len(ref), len(alt)) and ref[i] == alt[i]:
        i += 1
    ref, alt, pos = ref[i:], alt[i:], pos + i
    j = 0
    while j < min(len(ref), len(alt)) and ref[-1 - j] == alt[-1 - j]:
        j += 1
    return pos, pos + len(ref) - j


def classify(read, s, e, ref):
    """(covers, carries) for one read; `ref` maps a 0-based position to its base."""
    if read.reference_start > s - 1 or read.reference_end < e + 1:
        return False, False
    seq = read.query_sequence
    ref_pos, q, carries = read.reference_start, 0, False
    for op, ln in read.cigartuples:
        if op in (0, 7, 8):  # M = X
            for k in range(ln):
                p = ref_pos + k
                if s <= p < e:
                    b = seq[q + k].upper()
                    if b != 'N' and b != ref(p):
                        carries = True
            ref_pos += ln
            q += ln
        elif op == 1:  # I: between ref_pos - 1 and ref_pos
            if s <= ref_pos <= e:
                carries = True
            q += ln
        elif op in (2, 3):  # D N
            if ref_pos < max(e, s + 1) and ref_pos + ln > s:
                carries = True
            ref_pos += ln
        elif op == 4:  # S
            q += ln
    return True, carries


def counted(read):
    return not (read.is_unmapped or read.is_secondary or read.is_supplementary or read.is_duplicate
                or read.is_qcfail) and read.mapping_quality >= MIN_MAPQ


def site_counts(bam, chrom, s, e, ref):
    n = k = 0
    for read in bam.fetch(chrom, max(s - 1, 0), e + 1):
        if not counted(read):
            continue
        covers, carries = classify(read, s, e, ref)
        n += covers
        k += carries
    return n, k


def decide(n, k):
    if n < MIN_READS:
        return 'unchecked'
    return 'refuse' if k / n >= MIN_SHARE else 'pass'


def baseline(bam, chrom, p, ref):
    """spike's pileup SNP rule (loh::call_snps) at p: 'refuse' when het or hom-alt."""
    counts = Counter()
    for col in bam.pileup(chrom, p, p + 1, truncate=True, stepper='nofilter', ignore_overlaps=False,
                          ignore_orphans=False, min_base_quality=0, min_mapping_quality=0):
        for pr in col.pileups:
            if pr.is_del or pr.is_refskip or not counted(pr.alignment):
                continue
            b = pr.alignment.query_sequence[pr.query_position].upper()
            if b in 'ACGT':
                counts[b] += 1
    total = sum(counts.values())
    if total < MIN_READS:
        return 'unchecked'
    top = counts.most_common(2) + [('N', 0)] * 2
    f1, f2 = top[0][1] / total, top[1][1] / total
    if 0.2 <= f1 <= 0.8 and 0.2 <= f2 <= 0.8:
        return 'refuse'
    if f1 >= 0.9 and top[0][0] != ref(p):
        return 'refuse'
    return 'pass'


def gt_kind(rec):
    gt = rec.samples[0]['GT']
    if gt is None or None in gt:
        return None
    if sorted(gt) == [0, 1]:
        return 'het'
    if tuple(gt) == (1, 1):
        return 'hom'
    return None


def var_kind(ref, alt):
    if len(ref) == 1 and len(alt) == 1:
        return 'snv'
    if len(ref) != len(alt) and max(len(ref), len(alt)) <= 51:
        return 'indel'
    return None


def read_bed(path, chrom):
    spans = []
    with open(path) as f:
        for line in f:
            c, a, b = line.split('\t')[:3]
            if c == chrom:
                spans.append((int(a), int(b)))
    return sorted(spans)


def inside(spans, a, b):
    lo, hi = 0, len(spans)
    while lo < hi:
        mid = (lo + hi) // 2
        if spans[mid][1] <= a:
            lo = mid + 1
        else:
            hi = mid
    return lo < len(spans) and spans[lo][0] <= a and b <= spans[lo][1]


def main(bam_path, ref_path, calls_path, other_path, other_bed, out_dir):
    rng = random.Random(1)
    fasta = pysam.FastaFile(ref_path)
    seq = fasta.fetch(CHROM).upper()
    ref = lambda p: seq[p]
    bam = pysam.AlignmentFile(bam_path)

    # The sample's calls: K1's draws, and the margin K2 keeps from them.
    near_call = bytearray(len(seq) + 1)
    k1 = defaultdict(list)
    for rec in pysam.VariantFile(calls_path).fetch(CHROM):
        a, b = max(rec.start - CALL_MARGIN, 0), min(rec.stop + CALL_MARGIN, len(seq))
        near_call[a:b] = b'\x01' * (b - a)
        if rec.alts is None or len(rec.alts) != 1:
            continue
        gq = rec.samples[0].get('GQ')
        g, v = gt_kind(rec), var_kind(rec.ref, rec.alts[0])
        if gq is not None and gq >= 20 and g and v:
            k1[f'{g}_{v}'].append((rec.start, rec.ref, rec.alts[0]))

    taken = []

    def free(p):
        return all(abs(p - t) >= SPACING for t in taken)

    def draw(pool, n):
        out = []
        for site in rng.sample(pool, len(pool)):
            if len(out) == n:
                break
            if free(site[0]):
                out.append(site)
                taken.append(site[0])
        return out

    groups = {}
    for name, n in [('het_snv', 300), ('hom_snv', 300), ('het_indel', 200), ('hom_indel', 200)]:
        groups['K1_' + name] = draw(k1[name], n)

    def clear(a, b):
        return all(c != 'N' for c in seq[a:b]) and not any(near_call[a:b])

    def random_sites(n, make):
        out = []
        while len(out) < n:
            p = rng.randrange(1000, len(seq) - 1000)
            site = make(p)
            if site and free(p) and clear(p - 1, p + len(site[1]) + 1):
                out.append(site)
                taken.append(p)
        return out

    groups['K2a_snv'] = random_sites(500, lambda p: (p, seq[p], rng.choice([b for b in 'ACGT' if b != seq[p]])))
    groups['K2b_del'] = random_sites(250, lambda p: (lambda L: (p, seq[p:p + 1 + L], seq[p]))(rng.randint(1, 10)))
    groups['K2b_ins'] = random_sites(
        250, lambda p: (p, seq[p], seq[p] + ''.join(rng.choice('ACGT') for _ in range(rng.randint(1, 10)))))

    bed = read_bed(other_bed, CHROM)
    other = defaultdict(list)
    for rec in pysam.VariantFile(other_path).fetch(CHROM):
        if rec.alts is None or len(rec.alts) != 1:
            continue
        v = var_kind(rec.ref, rec.alts[0])
        if v and inside(bed, rec.start, rec.stop) and not any(near_call[rec.start:rec.stop + 1]):
            other[v].append((rec.start, rec.ref, rec.alts[0]))
    groups['K2c_snv'] = draw(other['snv'], 300)
    groups['K2c_indel'] = draw(other['indel'], 200)

    rows = []
    for group, sites in groups.items():
        for pos, r, a in sites:
            s, e = trim(pos, r, a)
            n, k = site_counts(bam, CHROM, s, e, ref)
            base = baseline(bam, CHROM, s, ref) if len(r) == 1 and len(a) == 1 else '-'
            rows.append((group, CHROM, pos + 1, r, a, s, e, n, k, decide(n, k), base))

    with open(f'{out_dir}/sites.tsv', 'w') as f:
        f.write('group\tchrom\tpos\tref\talt\ts\te\tN\tK\tdecision\tbaseline\n')
        for row in rows:
            f.write('\t'.join(map(str, row)) + '\n')
    with open(f'{out_dir}/sites.vcf', 'w') as f:
        f.write('##fileformat=VCFv4.2\n')
        f.write(f'##contig=<ID={CHROM},length={len(seq)}>\n')
        f.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
        for group, chrom, pos, r, a, *_ in sorted(rows, key=lambda x: x[2]):
            f.write(f'{chrom}\t{pos}\t{group}\t{r}\t{a}\t.\tPASS\tAF=0.5\n')

    verdict = True
    print(f'{"group":12} {"sites":>5} {"checked":>7} {"refused":>8} {"share":>6} {"unchecked":>9} {"baseline refused":>16}')
    for group in groups:
        g = [row for row in rows if row[0] == group]
        checked = [row for row in g if row[9] != 'unchecked']
        refused = sum(row[9] == 'refuse' for row in checked)
        share = refused / len(checked) if checked else float('nan')
        b_checked = [row for row in g if row[10] not in ('-', 'unchecked')]
        b_ref = sum(row[10] == 'refuse' for row in b_checked)
        b_text = f'{b_ref}/{len(b_checked)}' if b_checked else '-'
        if group.startswith('K1_'):
            ok = share >= K1_MIN[group[3:]]
        else:
            ok = share <= K2_MAX
        verdict &= ok
        print(f'{group:12} {len(g):5} {len(checked):7} {refused:8} {share:6.3f} {len(g) - len(checked):9} '
              f'{b_text:>16}  {"PASS" if ok else "FAIL"}')
    print('\nK1 misses and K2 refusals (checked sites):')
    for row in rows:
        miss = row[0].startswith('K1_') and row[9] == 'pass'
        false = row[0].startswith('K2') and row[9] == 'refuse'
        if miss or false:
            print('  ' + '\t'.join(map(str, row)))
    print(f'\nK: {"PASS" if verdict else "FAIL"}')


if __name__ == '__main__':
    if len(sys.argv) != 7:
        sys.exit(__doc__)
    main(*sys.argv[1:])
