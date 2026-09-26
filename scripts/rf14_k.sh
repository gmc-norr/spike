#!/bin/bash
# rf14_k.sh -- RF14, K: correct deletions through the real
# spike -> align -> merge -> validate loop, one slice_loop.sh run each.
#
# Usage: rf14_k.sh OUT_DIR SPIKE BAM REFERENCE SITES_TSV [THREADS]
#
# From SITES_TSV (rf14_sites.py):
#   rand sites: deletions of 50, 300, 1000 and 10000 bp starting there (VAF 0.5),
#               and 1000 bp at VAF 0.1       -> OUT_DIR/r_<pos>_<len>[_af01]/
#   hg sites:   HG002's own deletion, at its own POS and END (VAF 0.5)
#                                            -> OUT_DIR/h_<pos>/
# Every run passes --allow-resistant to spike, as validate_pipeline.sh does (RF8).
# rf14_planted.py scores them.
set -uo pipefail
out="$1"; spike="$2"; bam="$3"; ref="$4"; sites="$5"; threads="${6:-12}"
here="$(cd "$(dirname "$0")" && pwd)"
export SPIKE_ARGS="--allow-resistant"

run() {  # dir event
  bash "$here/slice_loop.sh" "$out/$1" "$spike" "$bam" "$ref" "$2" "$threads" < /dev/null
  echo "$1 spike=$(cat "$out/$1/spike.exit" 2>/dev/null) validate=$(cat "$out/$1/validate.json.exit" 2>/dev/null)"
}

mkdir -p "$out"
awk -F'\t' '$1 == "rand" {print $2, $3}' "$sites" | while read -r chrom pos; do
  for len in 50 300 1000 10000; do
    run "r_${pos}_${len}" "del:${chrom}:${pos}-$((pos + len))"
  done
  run "r_${pos}_1000_af01" "del:${chrom}:${pos}-$((pos + 1000));af=0.1"
done
awk -F'\t' '$1 == "hg" {print $2, $3, $4}' "$sites" | while read -r chrom pos end; do
  run "h_${pos}" "del:${chrom}:${pos}-${end}"
done
