"""Checks K1-K3 and B of docs/superpowers/plans/2026-10-04-low-af-evidence.md.

Runs spike at low allele fractions and compares the SIM_ALT_FRAGS it writes
with an independent count. That count rebuilds each event's haplotype from the
reference, places every spike read on it by 20-mer seeds (on either strand),
and applies the plan's rule:
  - a small variant shows in a pair when one read covers its changed bases and
    one base on each side;
  - a deletion shows when the pair's fragment spans the junction.

Usage: locate.py BAM REFERENCE MASTER_SPIKE NEW_SPIKE SITES_TSV OUT_DIR
"""
import gzip
import os
import re
import subprocess
import sys
from collections import Counter

import pysam

CHROM = 'chr20'
FLANK = 2000
SEED = 20
K1_AFS = (0.001, 0.01, 0.02, 0.05, 0.1)
K23_AFS = (0.02, 0.05)


def rc(s):
    return s.translate(str.maketrans('ACGTN', 'TGCAN'))[::-1]


def haplotype(fasta, kind, pos, ref, alt):
    """The event's haplotype and where it shows: ('bases', a, b) or ('junction', j).

    `pos` is the 1-based POS of a small variant, or the 0-based start of a DEL
    whose `ref` is its 0-based end."""
    if kind == 'del':
        start, end = pos, ref
        left = fasta.fetch(CHROM, max(0, start - FLANK), start).upper()
        right = fasta.fetch(CHROM, end, end + FLANK).upper()
        return left + right, ('junction', len(left))
    p0 = pos - 1
    left = fasta.fetch(CHROM, max(0, p0 - FLANK), p0).upper()
    right = fasta.fetch(CHROM, p0 + len(ref), p0 + len(ref) + FLANK).upper()
    prefix = 0
    while prefix < min(len(ref), len(alt)) and ref[prefix] == alt[prefix]:
        prefix += 1
    suffix = 0
    while (suffix < min(len(ref), len(alt)) - prefix
           and ref[len(ref) - 1 - suffix] == alt[len(alt) - 1 - suffix]):
        suffix += 1
    a, b = len(left) + prefix, len(left) + len(alt) - suffix
    return left + alt + right, ('bases', a, b)


def index(hap):
    idx = {}
    for i in range(len(hap) - SEED + 1):
        idx.setdefault(hap[i:i + SEED], []).append(i)
    return idx


def place(read, idx):
    """('+' or '-', start on the haplotype) by majority over the seeds, or None."""
    votes = Counter()
    for strand, seq in (('+', read), ('-', rc(read))):
        for o in range(0, len(seq) - SEED + 1, SEED):
            hits = idx.get(seq[o:o + SEED], [])
            if len(hits) == 1:
                votes[(strand, hits[0] - o)] += 1
    if not votes:
        return None
    return votes.most_common(1)[0][0]


def spike_pairs(out):
    reads = {}
    for mate in ('R1', 'R2'):
        with gzip.open(f'{out}/{mate}.fq.gz', 'rt') as fq:
            lines = fq.read().splitlines()
        for i in range(0, len(lines), 4):
            name = re.sub(r'/[12]$', '', lines[i][1:].split()[0])
            if name.startswith('SPIKE_'):
                reads.setdefault(name, []).append(lines[i + 1])
    return reads


def count(out, hap, shows):
    """(pairs showing, pairs placed, spike pairs) by the locator."""
    idx = index(hap)
    shown = placed = 0
    pairs = spike_pairs(out)
    for name, seqs in pairs.items():
        spots = [place(s, idx) for s in seqs]
        if len(seqs) != 2 or None in spots:
            continue
        mates = [(start, start + len(s), strand) for s, (strand, start) in zip(seqs, spots)]
        fwd = [m for m in mates if m[2] == '+']
        rev = [m for m in mates if m[2] == '-']
        if len(fwd) != 1 or len(rev) != 1:
            continue
        placed += 1
        if shows[0] == 'bases':
            a, b = shows[1], shows[2]
            hit = any(s <= a - 1 and e >= b + 1 for s, e, _ in mates)
        else:
            j = shows[1]
            hit = fwd[0][0] <= j - 1 and rev[0][1] >= j + 1
        shown += hit
    return shown, placed, len(pairs)


def alt_kmer_reads(out, fasta, pos, ref, alt):
    left = fasta.fetch(CHROM, pos - 11, pos - 1).upper()
    right = fasta.fetch(CHROM, pos - 1 + len(ref), pos - 1 + len(ref) + 10).upper()
    k = left + alt + right
    n = 0
    for seqs in spike_pairs(out).values():
        n += sum(k in s or rc(k) in s for s in seqs)
    return n


def run(spike, bam, ref_path, event, out):
    """{'ok', 'frags' (SIM_ALT_FRAGS, None where the binary writes none), 'log', 'error'}."""
    p = subprocess.run([spike, '--bam', bam, '--reference', ref_path, '--event', event,
                        '--seed', '1', '--threads', '16', '-o', out], capture_output=True, text=True)
    log = p.stdout + p.stderr
    open(f'{out}.log', 'w').write(log)
    if p.returncode != 0:
        errors = [line for line in log.splitlines() if line.startswith('Error:')]
        if not errors:
            raise SystemExit(f'spike exited {p.returncode} without an Error line for {event}:\n{log[-2000:]}')
        return {'ok': False, 'frags': None, 'log': log, 'error': errors[0][:150]}
    info = [line for line in open(f'{out}/truth.vcf') if not line.startswith('#')][0].split('\t')[7]
    m = re.search(r'SIM_ALT_FRAGS=(\d+)', info)
    return {'ok': True, 'frags': int(m.group(1)) if m else None, 'log': log, 'error': None}


def stripped_truth(path):
    keep = []
    for line in open(path):
        if line.startswith('##INFO=<ID=SIM_ALT_FRAGS,') or line.startswith('##INFO=<ID=SIM_VAF,'):
            continue
        keep.append(re.sub(r';SIM_ALT_FRAGS=[0-9.]+', '', line))
    return keep


def same_reads(a, b):
    files = ('R1.fq.gz', 'R2.fq.gz', 'replaced_reads.txt', 'fastq_removed_reads.txt')
    same = all(open(f'{a}/{f}', 'rb').read() == open(f'{b}/{f}', 'rb').read() for f in files)
    return same and stripped_truth(f'{a}/truth.vcf') == stripped_truth(f'{b}/truth.vcf')


def main(bam, ref_path, master, new, sites_tsv, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    fasta = pysam.FastaFile(ref_path)
    sites = [line.split('\t') for line in open(sites_tsv).read().splitlines()[1:]]
    rows = []
    for pos, ref, alt in sites:
        pos = int(pos)
        del5 = fasta.fetch(CHROM, pos - 1, pos + 4).upper()
        cases = [('snv', af, f'snp:{CHROM}:{pos}:{ref}:{alt};af={af}', (pos, ref, alt)) for af in K1_AFS]
        cases += [('small_del', af, f'snp:{CHROM}:{pos}:{del5}:{del5[0]};af={af}', (pos, del5, del5[0]))
                  for af in K23_AFS]
        cases += [('del', af, f'del:{CHROM}:{pos}-{pos + 1000};af={af}', (pos, pos + 1000, None))
                  for af in K23_AFS]
        for kind, af, event, (p, r, a) in cases:
            tag = f'{kind}_{pos}_{af}'
            new_out, master_out = f'{out_dir}/new/{tag}', f'{out_dir}/master/{tag}'
            os.makedirs(os.path.dirname(new_out), exist_ok=True)
            os.makedirs(os.path.dirname(master_out), exist_ok=True)
            n = run(new, bam, ref_path, event, new_out)
            m = run(master, bam, ref_path, event, master_out)
            if not n['ok']:
                rows.append((kind, af, pos, 'refused', n['error'], 'master too' if not m['ok'] else 'MASTER RAN'))
                print(rows[-1], flush=True)
                continue
            if n['frags'] is None:
                raise SystemExit(f'the new binary wrote no SIM_ALT_FRAGS for {event}')
            hap, shows = haplotype(fasta, 'del' if kind == 'del' else 'small', p, r, a)
            located, placed, total = count(new_out, hap, shows)
            kmer = alt_kmer_reads(new_out, fasta, p, r, a) if kind == 'snv' else ''
            warned = 'none of the' in n['log'] and 'shows it' in n['log']
            same = same_reads(new_out, master_out) if m['ok'] else 'MASTER REFUSED'
            rows.append((kind, af, pos, n['frags'], located, placed, total, kmer, warned, same))
            print(rows[-1], flush=True)

    with open(f'{out_dir}/k.tsv', 'w') as f:
        f.write('kind\taf\tpos\tSIM_ALT_FRAGS\tlocator\tplaced_pairs\tspike_pairs\talt_kmer_reads\twarned\treads_same_as_master\n')
        for row in rows:
            f.write('\t'.join(map(str, row)) + '\n')

    print()
    for kind, need in (('snv', 98), ('small_del', 39), ('del', 39)):
        ran = [r for r in rows if r[0] == kind and r[3] != 'refused']
        refused = [r for r in rows if r[0] == kind and r[3] == 'refused']
        equal = sum(r[3] == r[4] for r in ran)
        worst = max((abs(r[3] - r[4]) for r in ran), default=0)
        print(f'{kind}: ran {len(ran)}, refused {len(refused)}; equal {equal}, largest difference {worst}; '
              f'unplaced pairs {sum(r[6] - r[5] for r in ran)}')
        zeros = [r for r in ran if r[3] == 0]
        print(f'  SIM_ALT_FRAGS=0 in {len(zeros)}: locator 0 in {sum(r[4] == 0 for r in zeros)}, '
              f'ALT 21-mer 0 in {sum(r[7] == 0 for r in zeros if r[7] != "")}, warned in {sum(r[8] for r in zeros)}; '
              f'warned with SIM_ALT_FRAGS>0: {sum(r[8] for r in ran if r[3] > 0)}')
        by_af = Counter((r[1], r[3] == 0) for r in ran)
        print('  share with SIM_ALT_FRAGS=0 by AF: ' + ', '.join(
            f'{af}: {by_af[(af, True)]}/{by_af[(af, True)] + by_af[(af, False)]}'
            for af in sorted({r[1] for r in ran})))
    k1 = [r for r in rows if r[0] == 'snv' and r[3] != 'refused']
    print(f'B (K1 runs): reads and stripped truth identical to master in '
          f'{sum(r[9] is True for r in k1)} of {len(k1)}')


if __name__ == '__main__':
    if len(sys.argv) != 7:
        sys.exit(__doc__)
    main(*sys.argv[1:])
