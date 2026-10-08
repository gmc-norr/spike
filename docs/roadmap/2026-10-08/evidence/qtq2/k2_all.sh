#!/bin/bash
cd /tmp/claude-1066/-home-parlar-ai-dev-spike/23a35c29-9600-4da2-a2c9-618425a1e0d3/scratchpad/qtq2
tail -n +2 /tmp/claude-1066/-home-parlar-ai-dev-spike/23a35c29-9600-4da2-a2c9-618425a1e0d3/scratchpad/inspect/events.tsv | while IFS=$'\t' read id set chrom start end cls; do echo "$id $chrom $start $end"; done | xargs -P 6 -L 1 bash plant.sh > k2_plant.log 2>&1
echo "k2 planting done"
