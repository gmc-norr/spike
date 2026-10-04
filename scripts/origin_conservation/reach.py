"""Check R of docs/superpowers/plans/2026-10-04-origin-conservation.md (reported, not judged).

Runs `--edit-model origin` at random chr20 SNV sites under two spike binaries
(master and the fix) and reports how far the fix moves the origin depth and
the number of tiled reads.

Usage: reach.py BAM REFERENCE MASTER_SPIKE NEW_SPIKE OUT_DIR [SITES]
"""
import os
import random
import re
import shutil
import subprocess
import sys

import pysam

CHROM = 'chr20'
EDGE = 100_000
FLANK = 2000
DEPTH = re.compile(r'origin depth at \S+: ([0-9.]+)x \(the donor pool\'s there: ([0-9.]+)x\)')
TILING = re.compile(r'Tiling (\d+) synthetic reads')


def draw(ref_path, n):
    """`n` SNV sites, seed 1: (1-based pos, REF, ALT)."""
    rng = random.Random(1)
    fasta = pysam.FastaFile(ref_path)
    length = fasta.get_reference_length(CHROM)
    sites = []
    while len(sites) < n:
        pos = rng.randrange(EDGE, length - EDGE)
        ref = fasta.fetch(CHROM, pos, pos + 1).upper()
        if ref in 'ACGT':
            sites.append((pos + 1, ref, 'ACGT'[('ACGT'.index(ref) + 1) % 4]))
    return sites


def run(spike, bam, ref, site, out):
    pos, r, a = site
    cmd = [spike, '--bam', bam, '--reference', ref, '--event', f'snp:{CHROM}:{pos}:{r}:{a}',
           '--edit-model', 'origin', '--allow-resistant', '--seed', '1', '--threads', '16', '-o', out]
    p = subprocess.run(cmd, capture_output=True, text=True)
    log = p.stdout + p.stderr
    shutil.rmtree(out, ignore_errors=True)
    if p.returncode != 0:
        errors = [line for line in log.splitlines() if line.startswith('Error:')]
        if not errors:
            raise SystemExit(f'spike exited {p.returncode} without an Error line at {site}:\n{log[-2000:]}')
        return {'refused': errors[0][:160]}
    depth, tiling = DEPTH.search(log), TILING.search(log)
    if depth is None or tiling is None:
        raise SystemExit(f'no origin depth or tiling line at {site}:\n{log[-2000:]}')
    return {'depth': float(depth.group(1)), 'pool': float(depth.group(2)), 'tiled': int(tiling.group(1))}


def mapq0_share(bam, pos):
    n = z = 0
    for read in bam.fetch(CHROM, pos - 1 - FLANK, pos + FLANK):
        if read.is_secondary or read.is_supplementary or read.is_unmapped:
            continue
        n += 1
        z += read.mapping_quality == 0
    return z / n if n else float('nan')


def main(bam_path, ref_path, master, new, out_dir, n='150'):
    os.makedirs(out_dir, exist_ok=True)
    bam = pysam.AlignmentFile(bam_path)
    rows = []
    for i, site in enumerate(draw(ref_path, int(n))):
        old = run(master, bam_path, ref_path, site, f'{out_dir}/run_master')
        fix = run(new, bam_path, ref_path, site, f'{out_dir}/run_new')
        rows.append((i, site, old, fix, mapq0_share(bam, site[0])))
        print(i, site, old, fix, flush=True)

    with open(f'{out_dir}/r.tsv', 'w') as f:
        f.write('i\tpos\tref\talt\tmapq0_share\tmaster\tnew\n')
        for i, (pos, r, a), old, fix, z in rows:
            f.write(f'{i}\t{pos}\t{r}\t{a}\t{z:.3f}\t{old}\t{fix}\n')

    ran = [row for row in rows if 'tiled' in row[2] and 'tiled' in row[3]]
    refused_master = [row[0] for row in rows if 'refused' in row[2]]
    refused_new = [row[0] for row in rows if 'refused' in row[3]]
    print(f'\nsites {len(rows)}; refused by master {len(refused_master)}, by new {len(refused_new)}; '
          f'same sites refused: {refused_master == refused_new}')
    changed = [row for row in ran if row[2] != row[3]]
    print(f'ran under both {len(ran)}; depth or tiled count changed at {len(changed)}')

    def summary(name, values):
        values = sorted(values)
        print(f'  {name}: min {values[0]:.4f}, median {values[len(values) // 2]:.4f}, max {values[-1]:.4f}')

    summary('new / master origin depth', [row[3]['depth'] / row[2]['depth'] for row in ran if row[2]['depth'] > 0])
    summary('new / master tiled reads', [row[3]['tiled'] / row[2]['tiled'] for row in ran if row[2]['tiled'] > 0])
    print('five largest changes in tiled reads:')
    biggest = sorted(ran, key=lambda row: -abs(row[3]['tiled'] - row[2]['tiled']))[:5]
    for i, (pos, r, a), old, fix, z in biggest:
        print(f'  {CHROM}:{pos} mapq0 {z:.3f}: tiled {old["tiled"]} -> {fix["tiled"]}, '
              f'depth {old["depth"]} -> {fix["depth"]} (pool {old["pool"]})')


if __name__ == '__main__':
    if len(sys.argv) not in (6, 7):
        sys.exit(__doc__)
    main(*sys.argv[1:])
