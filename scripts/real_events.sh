#!/bin/bash
# real_events.sh -- run a list of events through slice_loop.sh, N at a time.
#
# T2's C4, T5's and T6's false-failure rates all need the same thing: many
# correct real-data events, each scored by `spike validate`. This runs a file of
# `--event` specs one per slice and leaves each run's directory for a scorer.
#
# Usage: real_events.sh OUT_DIR SPIKE BAM REFERENCE EVENTS_FILE [PARALLEL] [THREADS]
# Writes OUT_DIR/<n>/ per event (slice_loop.sh's layout) and OUT_DIR/events.txt.
# Keep PARALLEL*THREADS at or under 16: other sessions share the machine.
set -uo pipefail
out="$1"; spike="$2"; bam="$3"; ref="$4"; events="$5"
parallel="${6:-4}"; threads="${7:-4}"
here="$(cd "$(dirname "$0")" && pwd)"

mkdir -p "$out"
cp "$events" "$out/events.txt"
# One blank line and every later result is misattributed: this numbered events
# with `nl`, which numbers non-empty lines only, while real_events_score.py
# numbers with `enumerate` over every line. Both now refuse a blank line, and the
# numbering below is `awk NR` -- every line, blank ones included -- so the two
# cannot disagree even if this check is ever relaxed.
if grep -n '^[[:space:]]*$' "$out/events.txt" > /dev/null; then
    echo "real_events.sh: $events has a blank line, which would misnumber every" \
         "event after it (real_events_score.py numbers every line):" >&2
    grep -n '^[[:space:]]*$' "$out/events.txt" >&2
    exit 2
fi
one() { bash "$HERE/slice_loop.sh" "$OUT/$1" "$SPIKE" "$BAM" "$REF" "$2" "$THREADS"; }
export -f one
export HERE="$here" OUT="$out" SPIKE="$spike" BAM="$bam" REF="$ref" THREADS="$threads"
# BEFORE_SPIKE, if set, is validated against each run's own merged BAM too.
export BEFORE_SPIKE="${BEFORE_SPIKE:-}"
# SPIKE_ARGS, if set, is passed to spike itself (see slice_loop.sh).
export SPIKE_ARGS="${SPIKE_ARGS:-}"
awk '{printf "%d %s\n", NR, $0}' "$out/events.txt" \
  | xargs -P "$parallel" -L 1 bash -c 'one "$0" "$1"'

# The merged BAMs are the only large files and each has been scored by now.
find "$out" -name 'merged.bam*' -delete
