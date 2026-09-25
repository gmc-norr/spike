#!/bin/bash
# probe_loop.sh -- one review-probe event through the whole loop, on a merged BAM.
#
# The Codex review's probes (`scripts/review_sv_model.py`) build a synthetic donor
# BAM on a 40 kb `chrT`, and its `depths` helper reconstructs merge.sh in Python
# rather than producing a BAM. `spike validate`'s coverage and split-read checks
# need a real merged BAM, so this does the same loop `slice_loop.sh` does -- spike,
# align.sh, merge.sh, spike validate -- on a probe BAM, with no slicing.
#
# The probe reference must be bwa-mem2-indexed; a 40 kb reference indexes in ~1 s.
#
# Usage: probe_loop.sh OUT_DIR SPIKE DONOR_BAM INDEXED_REF EVENT [THREADS] [EXTRA...]
# Same output layout as slice_loop.sh, minus slice.bam.
set -uo pipefail
out="$1"; spike="$2"; bam="$3"; ref="$4"; event="$5"; threads="${6:-4}"
shift 6 2>/dev/null || shift 5

mkdir -p "$out"
rm -rf "$out/run"
# --flank 2000 is what the review's own `Probe.run` passes.
"$spike" --bam "$bam" --reference "$ref" --event "$event" --flank 2000 --seed 17 \
  -o "$out/run" > "$out/spike.log" 2>&1
echo $? > "$out/spike.exit"
[ "$(cat "$out/spike.exit")" = 0 ] || exit 0

bash "$out/run/align.sh" "$ref" "$threads" > "$out/align.log" 2>&1 || exit 0
bash "$out/run/merge.sh" "$bam" "$ref" "$threads" > "$out/merge.log" 2>&1 || exit 0

"$spike" validate --bam "$out/run/merged.bam" --truth "$out/run/truth.vcf" \
  --reference "$ref" --flank 2000 "$@" > "$out/validate.txt" 2> "$out/validate.log"
echo $? > "$out/validate.exit"
"$spike" validate --bam "$out/run/merged.bam" --truth "$out/run/truth.vcf" \
  --reference "$ref" --flank 2000 --json "$@" > "$out/validate.json" 2> "$out/validate.json.log"
echo $? > "$out/validate.json.exit"
