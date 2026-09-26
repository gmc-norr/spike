#!/bin/bash
# RF8 criteria C1, C2 and C3 on real data: each event run three ways.
# Usage: rf8_compare.sh OUT_DIR MASTER_SPIKE NEW_SPIKE BAM REFERENCE EVENTS_FILE
# EVENTS_FILE holds one --event spec per line. For event <n> it writes
# OUT_DIR/<n>/{master,new,allow}/{exit,log,run/}: master's binary, the new
# binary, and the new binary with --allow-resistant, all at --seed 1. Then it
# prints one line per event: the three exits and whether each output file
# matches master's (truth.vcf compared without its ##fileDate line).
set -euo pipefail
out="$1"; master="$2"; new="$3"; bam="$4"; ref="$5"; events="$6"
rm -rf "$out"; mkdir -p "$out"
one() {
    local n="$1" spec="$2" tag bin extra d
    for tag in master new allow; do
        bin="$new"; extra=""
        [ "$tag" = master ] && bin="$master"
        [ "$tag" = allow ] && extra="--allow-resistant"
        d="$out/$n/$tag"; mkdir -p "$d"
        set +e
        "$bin" --bam "$bam" --reference "$ref" --event "$spec" --seed 1 $extra \
            -o "$d/run" > "$d/log" 2>&1
        echo $? > "$d/exit"
        set -e
    done
}
export -f one; export out master new bam ref
nl -nln -w1 -s' ' "$events" | xargs -P 4 -L 1 bash -c 'one "$0" "$1"'

digest() {  # md5 of one output file, or "absent"
    local f="$1"
    [ -e "$f" ] || { echo absent; return; }
    case "$f" in
        *truth.vcf) grep -v '^##fileDate' "$f" | md5sum | cut -d' ' -f1 ;;
        *) md5sum < "$f" | cut -d' ' -f1 ;;
    esac
}
verdict() {  # "same" only when both files exist and match: two absent files are not a match
    local a="$1" b="$2"
    if [ "$a" = absent ] || [ "$b" = absent ]; then echo absent
    elif [ "$a" = "$b" ]; then echo same
    else echo DIFF
    fi
}
n=0
while read -r spec; do
    n=$((n + 1))
    line="$n $spec exit m/n/a=$(cat "$out/$n/master/exit")/$(cat "$out/$n/new/exit")/$(cat "$out/$n/allow/exit")"
    for f in R1.fq.gz R2.fq.gz replaced_reads.txt truth.vcf events.bed; do
        m=$(digest "$out/$n/master/run/$f")
        line="$line $f:new=$(verdict "$m" "$(digest "$out/$n/new/run/$f")")"
        line="$line,allow=$(verdict "$m" "$(digest "$out/$n/allow/run/$f")")"
    done
    echo "$line"
done < "$events"
