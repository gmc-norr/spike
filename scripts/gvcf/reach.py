"""Check R of docs/superpowers/plans/2026-10-04-gvcf-hom-alt.md (reported, not judged).

In random event-sized windows: how often the calls hold hom-alt SNVs and no het
SNV (so master drops the gVCF and uses the pileup alone), and how many of those
hom-alt calls spike's pileup rule (loh::call_snps) does not call hom-alt.

Usage: reach.py BAM REFERENCE CALLS_VCF OUT_TSV
"""
import random
import sys
from collections import Counter

import pysam

CHROM, WINDOWS, SIZE, MIN_MAPQ = 'chr20', 500, 4001, 20


def gt_kind(rec):
    gt = rec.samples[0]['GT']
    if gt is None or None in gt:
        return None
    return 'het' if sorted(gt) == [0, 1] else 'hom' if tuple(gt) == (1, 1) else None


def pileup_call(bam, p, ref_base):
    """spike's pileup call at p: ('hom', base), ('het', None), or ('none', None); with the depth."""
    counts = Counter()
    for col in bam.pileup(CHROM, p, p + 1, truncate=True, stepper='nofilter', ignore_overlaps=False,
                          ignore_orphans=False, min_base_quality=0, min_mapping_quality=0, max_depth=100000):
        for pr in col.pileups:
            r = pr.alignment
            if pr.is_del or pr.is_refskip or r.is_unmapped or r.is_secondary or r.is_supplementary \
                    or r.is_duplicate or r.is_qcfail or r.mapping_quality < MIN_MAPQ:
                continue
            b = r.query_sequence[pr.query_position].upper()
            if b in 'ACGT':
                counts[b] += 1
    total = sum(counts.values())
    if total < 10:
        return 'none', None, total
    top = counts.most_common(2) + [('N', 0)] * 2
    f1, f2 = top[0][1] / total, top[1][1] / total
    if 0.2 <= f1 <= 0.8 and 0.2 <= f2 <= 0.8:
        return 'het', None, total
    if f1 >= 0.9 and top[0][0] != ref_base:
        return 'hom', top[0][0], total
    return 'none', None, total


def main(bam_path, ref_path, calls_path, out_path):
    rng = random.Random(1)
    seq = pysam.FastaFile(ref_path).fetch(CHROM).upper()
    bam = pysam.AlignmentFile(bam_path)
    calls = pysam.VariantFile(calls_path)
    rows, affected = [], 0
    for w in range(WINDOWS):
        start = rng.randrange(1_000_000, 63_000_000)
        end = start + SIZE
        het, hom = 0, []
        for rec in calls.fetch(CHROM, start, end):
            alt = (rec.alts or ('.',))[0]
            if len(rec.ref) != 1 or len(alt) != 1 or alt in '.*' or not (start <= rec.start < end):
                continue
            kind = gt_kind(rec)
            het += kind == 'het'
            if kind == 'hom':
                hom.append((rec.start, alt.upper()))
        if het or not hom:
            rows.append((w, start, end, het, len(hom), '', '', '', ''))
            continue
        affected += 1
        for p, alt in hom:
            call, base, depth = pileup_call(bam, p, seq[p])
            kept = call == 'hom' and base == alt
            rows.append((w, start, end, het, len(hom), p + 1, alt, f'{call}:{depth}', 'kept' if kept else 'MISSED'))
    with open(out_path, 'w') as f:
        f.write('window\tstart\tend\thet_calls\thom_calls\thom_pos\talt\tpileup\tmaster\n')
        for row in rows:
            f.write('\t'.join(map(str, row)) + '\n')
    sites = [r for r in rows if r[8]]
    missed = [r for r in sites if r[8] == 'MISSED']
    print(f'windows {WINDOWS}; with a het call {len({r[0] for r in rows if r[3]})}; '
          f'hom-alt calls and no het call {affected} ({affected / WINDOWS:.3f}); no SNV call {len({r[0] for r in rows if not r[3] and not r[4]})}')
    print(f'hom-alt calls in those windows {len(sites)}; missed by the pileup {len(missed)} '
          f'({len(missed) / len(sites) if sites else float("nan"):.3f}); windows with a miss {len({r[0] for r in missed})}')
    for r in missed:
        print('  missed: window %s chr20:%s %s pileup %s' % (r[0], r[5], r[6], r[7]))


if __name__ == '__main__':
    if len(sys.argv) != 5:
        sys.exit(__doc__)
    main(*sys.argv[1:])
