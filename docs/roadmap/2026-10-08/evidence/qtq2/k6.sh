#!/usr/bin/env bash
# K6: one SNV and the 10 kb DEL on the novoalign HG002 chr20 BAM (31 quality values), new binary.
set -u
cd /tmp/claude-1066/-home-parlar-ai-dev-spike/23a35c29-9600-4da2-a2c9-618425a1e0d3/scratchpad/qtq2
REF=/home/parlar_ai/dev/spike/data/giab_hg38/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta
BAM=/home/parlar_ai/dev/spike/data/giab_hg38/HG002/HG002.GRCh38.chr20.bam
rm -rf k6
/tmp/claude-1066/-home-parlar-ai-dev-spike/23a35c29-9600-4da2-a2c9-618425a1e0d3/scratchpad/qtq2/spike-v2 --bam $BAM --reference "$REF" --event "snp:chr20:38550000:T:C;af=0.1" --event del:chr20:38900000-38910000 --threads 8 -o k6 > k6.spike.log 2>&1; echo "spike exit $?"
bash k6/align.sh > k6.align.log 2>&1; echo "align exit $?"
bash k6/merge.sh > k6.merge.log 2>&1; echo "merge exit $?"
