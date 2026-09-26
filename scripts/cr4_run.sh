#!/bin/bash
# CR4 criterion C4: run each of cr4_placements.py's 40 deletions on its own.
# Usage: cr4_run.sh OUT_DIR SPIKE BAM REFERENCE BENCHMARK_BED [CHROM]
# Writes OUT_DIR/cr4-c4-events.txt and OUT_DIR/cr4-c4/<n>/{exit,log,run/};
# score with cr4_score.py OUT_DIR. CHROM defaults to chr20 (RF8 uses chr1).
set -euo pipefail
out="$1"; spike="$2"; bam="$3"; ref="$4"; bed="$5"; chrom="${6:-chr20}"
here="$(cd "$(dirname "$0")" && pwd)"
python3 "$here/cr4_placements.py" "$bed" "$chrom" > "$out/cr4-c4-events.txt"
rm -rf "$out/cr4-c4"; mkdir -p "$out/cr4-c4"
one() {
    local d="$out/cr4-c4/$1"
    mkdir -p "$d"
    set +e
    "$spike" --bam "$bam" --reference "$ref" --event "$2" --seed 1 -o "$d/run" > "$d/log" 2>&1
    echo $? > "$d/exit"
}
export -f one; export out spike bam ref
nl -nln -w1 -s' ' "$out/cr4-c4-events.txt" | xargs -P 4 -L 1 bash -c 'one "$0" "$1"'
