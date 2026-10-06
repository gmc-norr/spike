#!/usr/bin/env bash
# The spike runs behind the deck's figures, on a chr20 slice (38.5-40.2 Mb) of the GIAB HG002
# NovaSeq PCR-free 35x BAM (bwa-mem2 2.2.1). The deck used spike at master 97aa3b4.
#   REF    GRCh38 no-alt analysis set FASTA (bwa-mem2 indexed)
#   SLICE  the chr20 slice BAM
#   SPIKE  spike binary (default: spike on PATH)
# Writes het/, hom/ and snv/ next to this script, each with merged.bam and the validate output.
set -u
cd "$(dirname "$0")"
REF="${REF:?set REF to the GRCh38 no-alt FASTA}"
SLICE="${SLICE:?set SLICE to the chr20 slice BAM}"
S="${SPIKE:-spike}"
run() {
  local name="$1"; shift
  rm -rf "$name"
  $S --bam "$SLICE" --reference "$REF" "$@" --threads 8 -o "$name" > "$name.spike.log" 2>&1 || { echo "$name spike exit $?"; tail -5 "$name.spike.log"; return; }
  bash "$name/align.sh" > "$name.align.log" 2>&1 || { echo "$name align exit $?"; return; }
  bash "$name/merge.sh" > "$name.merge.log" 2>&1 || { echo "$name merge exit $?"; return; }
  $S validate --bam "$name/merged.bam" --truth "$name/truth.vcf" --reference "$REF" > "$name.validate.txt" 2> "$name.validate.log"
  echo "$name validate exit $?"
}
run het --event del:chr20:38900000-38910000 --event ins:chr20:39000000:300
run hom --event "del:chr20:38900000-38910000;af=1.0"
args=(); while read -r e; do args+=(--event "$e"); done < snv_events.txt
run snv "${args[@]}"
