#!/usr/bin/env bash
# Throwaway Gate B runs for hp_apply.py (run from the qm worktree with the change applied).
export SPIKE_N7_BAM=/home/parlar_ai/dev/spike/data/giab_hg38/HG002/HG002.novaseq.pcr-free.35x.bwamem2.dedup.grch38_no_alt.bam
export SPIKE_REF=/home/parlar_ai/dev/spike/data/giab_hg38/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta
export CARGO_TARGET_DIR=$PWD/../qm-target
for qx in 0 1 2 3; do
  SPIKE_QX=$qx cargo test --release -- --ignored zz_gate_large_sample --nocapture 2>&1 | grep -E "^train|^QX|^  reads|panicked"
done
for qx in 0 4 8 12; do
  SPIKE_QX=$qx cargo test --release -- --ignored zz_gate_errors --nocapture 2>&1 | grep -E "^QX|^  |panicked"
done
