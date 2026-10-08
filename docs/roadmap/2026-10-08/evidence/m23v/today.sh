#!/bin/bash
set -euo pipefail
cd "$(dirname "$0")"
T=$(mktemp -d -p .)
mkfifo $T/raw2 $T/out2
cat r2.fq > $T/raw2 &
cat < $T/out2 > o2.fq &
SPIKE_REMOVED=removed.txt SPIKE_RAW2=$T/raw2 SPIKE_OUT2=$T/out2 SPIKE_STYLE1=$T/s1 SPIKE_STYLE2=$T/s2 SPIKE_WHY=$T/why SPIKE_COUNT=$T/count LC_ALL=C gawk -f pair.awk < r1.fq > o1.fq
wait
cat $T/count
rm -rf $T
