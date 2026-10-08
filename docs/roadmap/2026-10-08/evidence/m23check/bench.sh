#!/bin/bash
set -u
S=$1; B=$2
T=$(mktemp -d -p "$B")
mkfifo "$T/raw2"
pigz -dc "$S/fq/raw2.fq.gz" > "$T/raw2" &
start=$(date +%s.%N)
pigz -dc "$S/fq/raw1.fq.gz" | SPIKE_REMOVED="$S/run/fastq_removed_reads.txt" SPIKE_RAW2="$T/raw2" SPIKE_OUT2=/dev/null \
  SPIKE_STYLE1="$T/s1" SPIKE_STYLE2="$T/s2" SPIKE_WHY="$T/why" SPIKE_COUNT="$T/count" LC_ALL=C awk -f "$B/pair.awk" > /dev/null
wait
end=$(date +%s.%N); echo "pair-awk no-compress: $(echo "$end - $start" | bc) s; count: $(cat $T/count)"
# per-mate in parallel, with names to FIFOs compared by cmp
mkfifo "$T/n1" "$T/n2"
start=$(date +%s.%N)
cmp "$T/n1" "$T/n2" > "$T/cmp.out" 2>&1 & CMP=$!
( pigz -dc "$S/fq/raw1.fq.gz" | L="$S/run/fastq_removed_reads.txt" NAMES=/dev/null LC_ALL=C awk 'BEGIN{while((getline n<ENVIRON["L"])>0)r[n]=1} NR%4==1{name=$1;sub(/^@/,"",name);sub(/\/[12]$/,"",name);print name > "'"$T/n1"'";skip=(name in r);if(skip)d++} !skip{print} END{print d+0 > "/dev/stderr"}' > /dev/null ) & A1=$!
( pigz -dc "$S/fq/raw2.fq.gz" | L="$S/run/fastq_removed_reads.txt" LC_ALL=C awk 'BEGIN{while((getline n<ENVIRON["L"])>0)r[n]=1} NR%4==1{name=$1;sub(/^@/,"",name);sub(/\/[12]$/,"",name);print name > "'"$T/n2"'";skip=(name in r);if(skip)d++} !skip{print} END{print d+0 > "/dev/stderr"}' > /dev/null ) & A2=$!
wait $A1 $A2; wait $CMP; rc=$?
end=$(date +%s.%N); echo "per-mate + cmp no-compress: $(echo "$end - $start" | bc) s; cmp rc=$rc $(cat $T/cmp.out)"
rm -rf "$T"
