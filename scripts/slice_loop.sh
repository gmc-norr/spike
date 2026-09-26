#!/bin/bash
# slice_loop.sh -- one real-data spike -> align -> merge -> validate run on a slice.
#
# merge.sh over the whole-genome BAM is slow and large, so each event is run on a
# slice of it: the event's span grown by PAD on each side, indexed. That keeps the
# merge to seconds and the output to a few megabytes.
#
# Usage: slice_loop.sh OUT_DIR SPIKE BAM REFERENCE EVENT [THREADS] [EXTRA...]
#
# Writes, under OUT_DIR:
#   slice.bam(.bai)      the sliced donor
#   slice.exit           samtools view's exit status for the slice itself
#   run/                 spike's output, plus sim.bam and merged.bam
#   spike.exit           spike's exit status        spike.log
#   align.log merge.log
#   validate.txt         `spike validate`'s text table (stdout)
#   validate.log         its log (stderr)
#   validate.exit        its exit status
# The caller scores those; nothing here judges anything.
#
# `slice.exit` and `spike.exit` are separate files on purpose: a failed slice --
# a wrong BAM path, a contig name this BAM does not carry -- is a failure of this
# harness, and writing it into `spike.exit` made the scorer report it as "spike
# declined this event". A run whose slice failed leaves `slice.exit` non-zero and
# **no** `spike.exit` at all, which is a different state from spike refusing.
#
# Every report this run could write is deleted before anything runs: `align.sh`
# or `merge.sh` failing leaves the previous run's `validate.json` in place
# otherwise, and the scorer would score the old report as this run's.
#
# EXTRA... is passed to `spike validate` (e.g. --strict, --json).
set -uo pipefail
out="$1"; spike="$2"; bam="$3"; ref="$4"; event="$5"; threads="${6:-4}"
shift 6 2>/dev/null || shift 5
PAD=100000

# The event spec is `<type>:<chrom>:<start>-<end>` or `ins:<chrom>:<pos>:<len>`.
chrom=$(echo "$event" | cut -d: -f2)
coords=$(echo "$event" | cut -d: -f3)
start=$(echo "$coords" | cut -d- -f1)
end=$(echo "$coords" | cut -d- -f2)
[ "$end" = "$start" ] && end=$start
lo=$(( start > PAD ? start - PAD : 1 ))
hi=$(( end + PAD ))

mkdir -p "$out"
# A rerun must not be scorable on the last run's reports (see the note above).
rm -f "$out/spike.exit" "$out/slice.exit" \
      "$out/validate.txt" "$out/validate.log" "$out/validate.exit" \
      "$out/validate.json" "$out/validate.json.log" "$out/validate.json.exit" \
      "$out/validate.before.json" "$out/validate.before.log" "$out/validate.before.exit"
rm -rf "$out/run"

samtools view -@ "$threads" -b -o "$out/slice.bam" "$bam" "$chrom:$lo-$hi" \
  > "$out/slice.log" 2>&1
slice_rc=$?
echo "$slice_rc" > "$out/slice.exit"
# No `spike.exit` is written here: a slice that failed is this harness failing,
# not spike declining the event, and the scorer tells the two apart by which
# file it finds.
[ "$slice_rc" = 0 ] || exit 1
samtools index -@ "$threads" "$out/slice.bam" >> "$out/slice.log" 2>&1

"$spike" --bam "$out/slice.bam" --reference "$ref" --event "$event" --seed 1 \
  -o "$out/run" > "$out/spike.log" 2>&1
echo $? > "$out/spike.exit"
[ "$(cat "$out/spike.exit")" = 0 ] || exit 0   # a refused event is the caller's to report

bash "$out/run/align.sh" "$ref" "$threads" > "$out/align.log" 2>&1 || exit 0
bash "$out/run/merge.sh" "$out/slice.bam" "$ref" "$threads" > "$out/merge.log" 2>&1 || exit 0

"$spike" validate --bam "$out/run/merged.bam" --truth "$out/run/truth.vcf" \
  --reference "$ref" "$@" > "$out/validate.txt" 2> "$out/validate.log"
echo $? > "$out/validate.exit"

# The same report as JSON. The text table is space-padded, so a scorer that
# parses it has to guess where the event label ends; `--json` names every field.
# Both come from the same merged BAM, which is deleted right after.
"$spike" validate --bam "$out/run/merged.bam" --truth "$out/run/truth.vcf" \
  --reference "$ref" --json "$@" > "$out/validate.json" 2> "$out/validate.json.log"
echo $? > "$out/validate.json.exit"

# A second binary, given as $BEFORE_SPIKE in the environment, is validated against
# the SAME merged BAM. That is what a "the default did not move" check needs: two
# binaries on one alignment, not two alignments.
if [ -n "${BEFORE_SPIKE:-}" ]; then
  "$BEFORE_SPIKE" validate --bam "$out/run/merged.bam" --truth "$out/run/truth.vcf" \
    --reference "$ref" --json "$@" > "$out/validate.before.json" 2> "$out/validate.before.log"
  echo $? > "$out/validate.before.exit"
fi

rm -f "$out/run/merged_tmp.bam" "$out/run/outside.bam"
