#!/usr/bin/env bash
# Every result run but speed, in parallel groups.
cd /tmp/claude-1066/-home-parlar-ai-dev-spike/23a35c29-9600-4da2-a2c9-618425a1e0d3/scratchpad/qtq2
( ./runs.sh ./spike-v2 v2_s42 42 > runs_v2.out 2>&1; echo "v2 runs done" ) &
( for s in 42 2 3 4 5; do ./runs.sh ../qtq/spike-master master_s$s $s; done > runs_master.out 2>&1; echo "master runs done" ) &
( ./k2_all.sh > k2_all.out 2>&1; echo "k2 done" ) &
( bash k6.sh > k6.out 2>&1; echo "k6 done" ) &
( cd ../qm && BAM=/home/parlar_ai/dev/spike/data/giab_hg38/HG002/HG002.novaseq.pcr-free.35x.bwamem2.dedup.grch38_no_alt.bam \
   REF=/home/parlar_ai/dev/spike/data/giab_hg38/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta \
   Q100_VCF=/home/parlar_ai/dev/spike/data/giab_hg38/HG002/GRCh38_HG2-T2TQ100-V1.1.vcf.gz OUT=$PWD/../qtq2/k2b \
   CARGO_TARGET_DIR=$PWD/../qm-target bash -c 'mkdir -p $OUT && docs/analysis/quality-model-v2/k2b.sh' > ../qtq2/k2b.out 2>&1; echo "k2b done" ) &
wait
echo "ALL DONE"
