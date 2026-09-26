#!/bin/bash
# RF6 criterion C2: the junction row must go red when the spike-in is wrong.
# Usage: rf6_red.sh SPIKE REFERENCE RUN_DIR...
# Each RUN_DIR is a slice_loop.sh output whose run/merged.bam was kept. Two
# validate runs each, printing the exit and the junction_sequence row:
#   shifted: run/truth.vcf with END moved 50 bp right, on run/merged.bam
#   donor:   run/truth.vcf as written, on slice.bam (never spiked)
set -euo pipefail
spike="$1"; ref="$2"; shift 2
for d in "$@"; do
    [ -s "$d/run/merged.bam" ] || { echo "$d: no run/merged.bam" >&2; exit 2; }
    LC_ALL=C awk 'BEGIN { FS = OFS = "\t" }
        /^#/ { print; next }
        { n = split($8, a, ";")
          for (i = 1; i <= n; i++) if (a[i] ~ /^END=/) a[i] = "END=" (substr(a[i], 5) + 50)
          s = a[1]; for (i = 2; i <= n; i++) s = s ";" a[i]
          $8 = s; print }' "$d/run/truth.vcf" > "$d/truth.shifted.vcf"
    for case in shifted donor; do
        if [ "$case" = shifted ]; then bam="$d/run/merged.bam"; truth="$d/truth.shifted.vcf"
        else bam="$d/slice.bam"; truth="$d/run/truth.vcf"; fi
        set +e
        "$spike" validate --bam "$bam" --truth "$truth" --reference "$ref" --json \
            > "$d/red.$case.json" 2> "$d/red.$case.log"
        rc=$?
        set -e
        row=$(python3 -c '
import json, sys
rows = [r for r in json.load(open(sys.argv[1]))["checks"] if r["check"] == "junction_sequence"]
print(" | ".join(f"{r[\"observed\"]} {\"PASS\" if r[\"pass\"] else \"FAIL\"}" for r in rows) or "absent")
' "$d/red.$case.json")
        echo "$(basename "$d") $case exit=$rc junction_sequence=[$row] END=$(grep -v '^#' "$truth" | grep -o 'END=[0-9]*')"
    done
done
