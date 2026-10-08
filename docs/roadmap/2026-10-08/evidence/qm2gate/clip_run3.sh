#!/usr/bin/env bash
# Seed check of the clip test: newA and newAB with seeds 22 and 23 (same model training, other draws).
set -euo pipefail
export SPIKE_N7_BAM=/home/parlar_ai/dev/spike/data/giab_hg38/HG002/HG002.novaseq.pcr-free.35x.bwamem2.dedup.grch38_no_alt.bam
export SPIKE_REF=/home/parlar_ai/dev/spike/data/giab_hg38/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta
export CARGO_TARGET_DIR=$PWD/../qm-target SPIKE_GATE_OUT=$PWD/../qm2gate/clips
VCF=/home/parlar_ai/dev/spike/data/giab_hg38/HG002/GRCh38_HG2-T2TQ100-V1.1.vcf.gz
args=()
for seed in 22 23; do
  for spec in newA:63 newAB:191; do
    set_=${spec%%:*}_s$seed; qx=${spec##*:}
    if [ "$qx" = 191 ]; then export SPIKE_GATE_MASK=$PWD/../qm2gate/q100_blocks.tsv; else unset SPIKE_GATE_MASK; fi
    SPIKE_GATE_SEED=$seed SPIKE_GATE_SET=$set_ SPIKE_QX=$qx cargo test --release -- --ignored zz_gate_clips --nocapture 2>&1 | grep -E "^set|panicked"
    bwa-mem2 mem -t 32 "$SPIKE_REF" $SPIKE_GATE_OUT/${set_}_R1.fq $SPIKE_GATE_OUT/${set_}_R2.fq 2> $SPIKE_GATE_OUT/${set_}.align.log | samtools view -b -o $SPIKE_GATE_OUT/${set_}.bam -
    args+=("$set_=$SPIKE_GATE_OUT/$set_.bam")
  done
done
python3 ../qm2gate/clipscore.py "$SPIKE_REF" $VCF real=$SPIKE_GATE_OUT/real.bam newA=$SPIKE_GATE_OUT/newA.bam newAB=$SPIKE_GATE_OUT/newAB.bam "${args[@]}"
