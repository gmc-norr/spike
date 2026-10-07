#!/usr/bin/env bash
# K2b of docs/superpowers/plans/2026-10-07-quality-model-v2.md: the big held-out clip test.
# The ignored test synth::tests::measure_clip_share samples BAM as a run does, learns from half the
# pairs and writes the other half as sequenced (real) and as spike makes them (spike, one set per
# seed); each set is aligned with bwa-mem2 and scored with clipscore.py.
#   BAM      an indexed BAM (the plan used HG002 35x, novaseq pcr-free, bwa-mem2)
#   REF      GRCh38 no-alt FASTA, bwa-mem2 indexed
#   Q100_VCF GIAB HG002 T2T-Q100 v1.1 benchmark VCF (bad-end clips near its variants are left out)
#   OUT      an empty directory for the FASTQ, BAM and scores
# Run from the repository root.
set -euo pipefail
: "${BAM:?}" "${REF:?}" "${Q100_VCF:?}" "${OUT:?}"
here=$(cd "$(dirname "$0")" && pwd)
make_set() {  # make_set NAME SEED
  SPIKE_K2B_BAM=$BAM SPIKE_K2B_REF=$REF SPIKE_K2B_OUT=$OUT SPIKE_K2B_SET=$1 SPIKE_K2B_SEED=$2 \
    cargo test --release -- --ignored measure_clip_share --nocapture 2>&1 | grep -E "^set |panicked"
}
align() { bwa-mem2 mem -t 32 "$REF" "$OUT/$1_R1.fq" "$OUT/$1_R2.fq" 2> "$OUT/$1.align.log" | samtools view -b -o "$OUT/$1.bam" -; }
make_set real 21 && align real
sets=(real="$OUT/real.bam")
for seed in 21 22 23; do
  make_set spike "$seed" && mv "$OUT/spike_R1.fq" "$OUT/spike$seed"_R1.fq && mv "$OUT/spike_R2.fq" "$OUT/spike$seed"_R2.fq
  align "spike$seed"; sets+=("spike$seed=$OUT/spike$seed.bam")
done
python3 "$here/clipscore.py" "$REF" "$Q100_VCF" "${sets[@]}" | tee "$OUT/k2b.txt"
