#!/bin/bash
# rf13b_k.sh -- RF13 second attempt, K: correct insertions through the real
# spike -> align -> merge -> validate loop, one slice_loop.sh run each.
#
# Usage: rf13_k.sh OUT_DIR SPIKE BAM REFERENCE SITES_TSV [THREADS]
#
# From SITES_TSV (rf13b_sites.py):
#   rand sites: random bases at 1, 2, 4, 15, 45 and 300 bp (VAF 0.5), and
#               45 bp at VAF 0.1       -> OUT_DIR/r_<pos>_<len>[_af01]/
#   hg sites:   HG002's own inserted bases at its own POS (VAF 0.5)
#                                      -> OUT_DIR/h_<pos>/
# rf13b_planted.py scores them.
set -uo pipefail
out="$1"; spike="$2"; bam="$3"; ref="$4"; sites="$5"; threads="${6:-12}"
here="$(cd "$(dirname "$0")" && pwd)"

run() {  # dir event
  bash "$here/slice_loop.sh" "$out/$1" "$spike" "$bam" "$ref" "$2" "$threads" < /dev/null
  echo "$1 spike=$(cat "$out/$1/spike.exit" 2>/dev/null) validate=$(cat "$out/$1/validate.json.exit" 2>/dev/null)"
}

mkdir -p "$out"
awk -F'\t' '$1 == "rand" {print $2, $3}' "$sites" | while read -r chrom pos; do
  for len in 1 2 4 15 45 300; do
    run "r_${pos}_${len}" "ins:${chrom}:${pos}:${len}"
  done
  run "r_${pos}_45_af01" "ins:${chrom}:${pos}:45;af=0.1"
done
awk -F'\t' '$1 == "hg" {print $2, $3, $4}' "$sites" | while read -r chrom pos seq; do
  run "h_${pos}" "ins:${chrom}:${pos}:${seq}"
done
