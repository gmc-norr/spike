#!/bin/bash
# README "Run time", clean rows: old (spike-dp) and new (spike-doc) binaries back to back, the
# table's method: --seed 1 --allow-resistant, 3 runs at 1 and 8 threads, 1 at 4 and 16.
set -uo pipefail
S="$(cd "$(dirname "$0")" && pwd)"
BAM=/home/parlar_ai/dev/sv_caller/data/giab_hg38/HG002/HG002.novaseq.pcr-free.35x.bwamem2.dedup.grch38_no_alt.bam
REF=/home/parlar_ai/dev/spike/data/giab_hg38/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta
mkdir -p "$S/rt"
for ev in "del:chr20:7119236-7120236" "del:chr20:14550000-17550000" "dup:chr20:14550000-17550000"; do
  for t in 1 4 8 16; do
    n=1; [ "$t" = 1 ] || [ "$t" = 8 ] && n=3
    for i in $(seq 1 $n); do
      for b in spike-dp spike-doc; do
        d="$S/rt/out_$b"; rm -rf "$d"
        /usr/bin/time -f "%e %M" -o "$S/rt/time.txt" "$S/$b" --bam "$BAM" --reference "$REF" --event "$ev" \
            --seed 1 --allow-resistant --threads "$t" -o "$d" > "$S/rt/log_$b.txt" 2>&1
        rc=$?
        echo "$ev $b t$t run$i exit $rc $(cat "$S/rt/time.txt")"
      done
    done
  done
done
