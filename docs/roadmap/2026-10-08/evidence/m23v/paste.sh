#!/bin/bash
set -euo pipefail
cd "$(dirname "$0")"
T=$(mktemp -d -p .)
mkfifo $T/out2
cat < $T/out2 > o2.fq &
paste <(paste - - - - < r1.fq) <(paste - - - - < r2.fq) | SPIKE_REMOVED=removed.txt SPIKE_OUT2=$T/out2 SPIKE_WHY=$T/why SPIKE_COUNT=$T/count LC_ALL=C gawk -f paste.awk > o1.fq
wait
cat $T/count
rm -rf $T
