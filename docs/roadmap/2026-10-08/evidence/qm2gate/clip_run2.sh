#!/usr/bin/env bash
# Clip test rerun (README.txt): new-A (QX 63) and new-AB (QX 191, slips with the Q100 mask), scored with
# the same real.bam. From the qm worktree with clip_all.patch + hp_apply6.py applied.
set -euo pipefail
export SPIKE_N7_BAM=/home/parlar_ai/dev/spike/data/giab_hg38/HG002/HG002.novaseq.pcr-free.35x.bwamem2.dedup.grch38_no_alt.bam
export SPIKE_REF=/home/parlar_ai/dev/spike/data/giab_hg38/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta
export CARGO_TARGET_DIR=$PWD/../qm-target SPIKE_GATE_OUT=$PWD/../qm2gate/clips
VCF=/home/parlar_ai/dev/spike/data/giab_hg38/HG002/GRCh38_HG2-T2TQ100-V1.1.vcf.gz
for spec in newA:63 newAB:191; do
  set_=${spec%%:*}; qx=${spec##*:}
  if [ "$qx" = 191 ]; then export SPIKE_GATE_MASK=$PWD/../qm2gate/q100_blocks.tsv; fi
  SPIKE_GATE_SET=$set_ SPIKE_QX=$qx cargo test --release -- --ignored zz_gate_clips --nocapture 2>&1 | grep -E "^set|^slips|^  [ACGT]:|panicked"
  bwa-mem2 mem -t 32 "$SPIKE_REF" $SPIKE_GATE_OUT/${set_}_R1.fq $SPIKE_GATE_OUT/${set_}_R2.fq 2> $SPIKE_GATE_OUT/${set_}.align.log | samtools view -b -o $SPIKE_GATE_OUT/${set_}.bam -
done
python3 ../qm2gate/clipscore.py "$SPIKE_REF" $VCF real=$SPIKE_GATE_OUT/real.bam v1=$SPIKE_GATE_OUT/v1.bam new=$SPIKE_GATE_OUT/new.bam newA=$SPIKE_GATE_OUT/newA.bam newAB=$SPIKE_GATE_OUT/newAB.bam
python3 ../qm2gate/clipdiag.py "$SPIKE_REF" real=$SPIKE_GATE_OUT/real.bam newA=$SPIKE_GATE_OUT/newA.bam newAB=$SPIKE_GATE_OUT/newAB.bam
