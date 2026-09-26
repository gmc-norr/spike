#!/bin/bash
# rf11_k1.sh -- RF11 plan, K1: correct insertions of 20-49 bp through the real
# spike -> align -> merge -> validate loop, one slice_loop.sh run each.
#
# Usage: rf11_k1.sh OUT_DIR SPIKE BAM REFERENCE SITES_TSV [THREADS]
#
# Runs every `k1` site of SITES_TSV (rf11_sites.py) at each length in LENGTHS,
# one after another, into OUT_DIR/ins_<pos>_<len>/. rf11_score.py scores them.
set -uo pipefail
out="$1"; spike="$2"; bam="$3"; ref="$4"; sites="$5"; threads="${6:-12}"
here="$(cd "$(dirname "$0")" && pwd)"
LENGTHS="20 30 40 45 49"

mkdir -p "$out"
awk -F'\t' '$1 == "k1" {print $2, $3}' "$sites" | while read -r chrom pos; do
  for len in $LENGTHS; do
    bash "$here/slice_loop.sh" "$out/ins_${pos}_${len}" "$spike" "$bam" "$ref" \
      "ins:${chrom}:${pos}:${len}" "$threads" < /dev/null   # keep the site list on stdin
    echo "ins:${chrom}:${pos}:${len} spike=$(cat "$out/ins_${pos}_${len}/spike.exit" 2>/dev/null) validate=$(cat "$out/ins_${pos}_${len}/validate.json.exit" 2>/dev/null)"
  done
done
