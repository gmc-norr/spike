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
#   run/                 spike's output, plus sim.bam and merged.bam
#   spike.exit           spike's exit status        spike.log
#   align.log merge.log
#   validate.txt         `spike validate`'s text table (stdout)
#   validate.log         its log (stderr)
#   validate.exit        its exit status
# The caller scores those; nothing here judges anything.
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
samtools view -@ "$threads" -b -o "$out/slice.bam" "$bam" "$chrom:$lo-$hi" \
  > "$out/slice.log" 2>&1 || { echo "slice failed" > "$out/spike.exit"; exit 1; }
samtools index -@ "$threads" "$out/slice.bam" >> "$out/slice.log" 2>&1

rm -rf "$out/run"
"$spike" --bam "$out/slice.bam" --reference "$ref" --event "$event" --seed 1 \
  -o "$out/run" > "$out/spike.log" 2>&1
echo $? > "$out/spike.exit"
[ "$(cat "$out/spike.exit")" = 0 ] || exit 0   # a refused event is the caller's to report

bash "$out/run/align.sh" "$ref" "$threads" > "$out/align.log" 2>&1 || exit 0
bash "$out/run/merge.sh" "$out/slice.bam" "$ref" "$threads" > "$out/merge.log" 2>&1 || exit 0

"$spike" validate --bam "$out/run/merged.bam" --truth "$out/run/truth.vcf" \
  --reference "$ref" "$@" > "$out/validate.txt" 2> "$out/validate.log"
echo $? > "$out/validate.exit"

# Keep the disk use down: the merged BAM is the only large file and it has been
# scored by the time we get here. `--keep-merged` in EXTRA is not a validate flag,
# so a caller that needs it keeps it by copying before calling again.
rm -f "$out/run/merged_tmp.bam" "$out/run/outside.bam"
