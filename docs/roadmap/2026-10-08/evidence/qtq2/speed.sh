#!/usr/bin/env bash
# Speed check of the v2 plan: dup:chr20:14550000-15550000 on the 35x BAM, --seed 1 --threads 8,
# master and v2 alternating, 3 runs each; wall seconds.
cd /tmp/claude-1066/-home-parlar-ai-dev-spike/23a35c29-9600-4da2-a2c9-618425a1e0d3/scratchpad/qtq2
REF=/home/parlar_ai/dev/spike/data/giab_hg38/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta
BAM=/home/parlar_ai/dev/spike/data/giab_hg38/HG002/HG002.novaseq.pcr-free.35x.bwamem2.dedup.grch38_no_alt.bam
for i in 1 2 3; do
  for b in master v2; do
    bin=../qtq/spike-master; [ $b = v2 ] && bin=./spike-v2
    rm -rf speed_$b
    t0=$(date +%s.%N)
    $bin --bam $BAM --reference $REF --event dup:chr20:14550000-15550000 --seed 1 --threads 8 -o speed_$b > speed_$b.log 2>&1
    t1=$(date +%s.%N)
    echo "$b $i $(echo "$t1 - $t0" | bc) exit $?"
  done
done
