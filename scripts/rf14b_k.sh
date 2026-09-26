#!/bin/bash
# rf14b_k.sh -- RF14 second attempt, K: correct deletions through the real
# spike -> align -> merge -> validate loop, one slice_loop.sh run each.
#
# Usage: rf14b_k.sh OUT_DIR SPIKE BAM REFERENCE SITES_TSV [THREADS]
#
# From SITES_TSV (rf14b_sites.py):
#   rand sites: deletions of 50, 300, 1000 and 10000 bp starting there (VAF 0.5)
#                                            -> OUT_DIR/r_<pos>_<len>/
#   hg sites:   HG002's own deletion, at its own POS and END (VAF 0.5)
#                                            -> OUT_DIR/h_<pos>/
#   low sites:  a 1000 bp deletion at VAF 0.1 -> OUT_DIR/l_<pos>/
# Every run passes --allow-resistant to spike, as validate_pipeline.sh does (RF8);
# rf14b_score.py judges the runs above SIM_RESIST 0.5 apart, as the plan locks.
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
done
awk -F'\t' '$1 == "hg" {print $2, $3, $4}' "$sites" | while read -r chrom pos end; do
  run "h_${pos}" "del:${chrom}:${pos}-${end}"
done
awk -F'\t' '$1 == "low" {print $2, $3}' "$sites" | while read -r chrom pos; do
  run "l_${pos}" "del:${chrom}:${pos}-$((pos + 1000));af=0.1"
done
