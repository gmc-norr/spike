#!/bin/bash
set -euo pipefail
cd "$(dirname "$0")"
T=$(mktemp -d -p .)
mkfifo $T/n1 $T/n2
cmp $T/n1 $T/n2 > $T/cmp 2>&1 &
C=$!
SPIKE_REMOVED=removed.txt NAMES=$T/n1 WHY=$T/w1 COUNT=$T/c1 LC_ALL=C gawk -f mate.awk < r1.fq > o1.fq &
A1=$!
SPIKE_REMOVED=removed.txt NAMES=$T/n2 WHY=$T/w2 COUNT=$T/c2 LC_ALL=C gawk -f mate.awk < r2.fq > o2.fq &
A2=$!
wait $A1; wait $A2; wait $C
cat $T/c1 $T/c2
rm -rf $T
