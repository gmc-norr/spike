#!/usr/bin/env bash
# Offline quality-model tests (Python, outside spike). Real reads from the SNV run's 25 windows.
#   REF      GRCh38 no-alt FASTA, bwa-mem2 indexed
#   SLICE    the chr20 slice BAM (training data for e5b, e7, e8, e9, e9m)
#   SNV_BAM  merged.bam of the 25-SNV run (docs/presentation/scripts/runs.sh)
# e4m.py and e9m.py need ../clips/nonerr_mask.pkl (run ../clips/run.sh first).
set -euo pipefail
cd "$(dirname "$0")"
: "${REF:?}" "${SLICE:?}" "${SNV_BAM:?}"
align() { bwa-mem2 mem -t 32 "$REF" "$1" "$2" 2> "$3.align.log" | samtools view -b -o "$3.bam" -; }
python3 load.py                     # reads.pkl: qualities, mismatches (clipped bases placed), variant sites masked
python3 e12.py                      # quality by base; error rate by Q and read class
python3 e3.py                       # chain vs copied strings vs 4-class chain, held out
python3 e4.py && align R1.fq.gz R2.fq.gz e4 && python3 score_e4.py      # soft clips: real, chain, copyQ, copyQE
python3 e5.py && python3 e5b.py     # bits per quality: how much the bases add
python3 e6.py                       # fqzcomp-style generator, first version
python3 e7.py && align e7_R1.fq.gz e7_R2.fq.gz e7 && python3 count_clips.py e7.bam
python3 e8.py                       # bigger contexts with backoff, bits per quality
python3 e9.py && align e9_R1.fq.gz e9_R2.fq.gz e9 && python3 count_clips.py e9.bam
python3 clumping.py                 # errors in the last 20 cycles, all clips counted as errors (superseded)
python3 e4m.py && align e4m_R1.fq.gz e4m_R2.fq.gz e4m
python3 e9m.py && align e9m_R1.fq.gz e9m_R2.fq.gz e9m
python3 score_masked.py             # clips and clumping with non-error tails masked
