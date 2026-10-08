#!/usr/bin/env bash
# Clip test runs (README.txt). From the qm worktree with hp_all.patch + hp_apply5.py (8 classes) applied.
set -euo pipefail
export SPIKE_N7_BAM=/home/parlar_ai/dev/spike/data/giab_hg38/HG002/HG002.novaseq.pcr-free.35x.bwamem2.dedup.grch38_no_alt.bam
export SPIKE_REF=/home/parlar_ai/dev/spike/data/giab_hg38/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta
export CARGO_TARGET_DIR=$PWD/../qm-target SPIKE_GATE_OUT=$PWD/../qm2gate/clips
VCF=/home/parlar_ai/dev/spike/data/giab_hg38/HG002/GRCh38_HG2-T2TQ100-V1.1.vcf.gz
for spec in real:0 v1:0 new:63; do
  set_=${spec%%:*}; qx=${spec##*:}
  SPIKE_GATE_SET=$set_ SPIKE_QX=$qx cargo test --release -- --ignored zz_gate_clips --nocapture 2>&1 | grep -E "^set|panicked"
  bwa-mem2 mem -t 32 "$SPIKE_REF" $SPIKE_GATE_OUT/${set_}_R1.fq $SPIKE_GATE_OUT/${set_}_R2.fq 2> $SPIKE_GATE_OUT/${set_}.align.log | samtools view -b -o $SPIKE_GATE_OUT/${set_}.bam -
done
python3 ../qm2gate/clipscore.py "$SPIKE_REF" $VCF real=$SPIKE_GATE_OUT/real.bam v1=$SPIKE_GATE_OUT/v1.bam new=$SPIKE_GATE_OUT/new.bam
