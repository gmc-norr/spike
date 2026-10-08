#!/bin/bash
# usage: run.sh OUTDIR
W=/tmp/claude-1066/-home-parlar-ai-dev-spike/23a35c29-9600-4da2-a2c9-618425a1e0d3/scratchpad/m26
G=/home/parlar_ai/dev/spike/data/giab_hg38
EV=$(sed 's/^/--event /' $W/events.txt | tr '\n' ' ')
/usr/bin/time -f "WALL %e USER %U SYS %S MAXRSS_KB %M" /home/parlar_ai/dev/spike/target/release/spike --bam $G/HG002/HG002.novaseq.pcr-free.35x.bwamem2.dedup.grch38_no_alt.bam --reference $G/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta $EV --seed 1 --threads 1 --allow-resistant -o $1 > $1.log 2>&1 &
TP=$!
sleep 0.3
PID=$(pgrep -P $TP spike)
: > $1.rss
while kill -0 $TP 2>/dev/null; do r=$(awk '/VmRSS/{print $2}' /proc/$PID/status 2>/dev/null); [ -n "$r" ] && echo "$(date +%s.%N | cut -c1-14) $r" >> $1.rss; sleep 0.25; done
wait $TP; echo exit $?
