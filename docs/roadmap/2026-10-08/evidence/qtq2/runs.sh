#!/usr/bin/env bash
# K1 and K7: the round trips of docs/presentation on the chr20 slice. runs.sh BIN OUTDIR SEED
set -u
BIN=$(readlink -f "$1"); OUT=$2; SEED=$3
mkdir -p "$OUT"; cd "$OUT"
REF=/home/parlar_ai/dev/spike/data/giab_hg38/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta
SLICE=/tmp/claude-1066/-home-parlar-ai-dev-spike/23a35c29-9600-4da2-a2c9-618425a1e0d3/scratchpad/ddp/slice.bam
run() {
  local name="$1"; shift
  rm -rf "$name"
  $BIN --bam $SLICE --reference "$REF" "$@" --seed $SEED --threads 8 -o "$name" > "$name.spike.log" 2>&1 || { echo "$OUT/$name spike exit $?"; tail -5 "$name.spike.log"; return; }
  bash "$name/align.sh" > "$name.align.log" 2>&1 || { echo "$OUT/$name align exit $?"; return; }
  bash "$name/merge.sh" > "$name.merge.log" 2>&1 || { echo "$OUT/$name merge exit $?"; return; }
  $BIN validate --bam "$name/merged.bam" --truth "$name/truth.vcf" --reference "$REF" > "$name.validate.txt" 2> "$name.validate.log"
  echo "$OUT/$name validate exit $?"
}
run het --event del:chr20:38900000-38910000 --event ins:chr20:39000000:300
run m4 --event snp:chr20:39000000:TGG:T --event snp:chr20:39100000:T:TCCGG --event snp:chr20:39200000:AT:GC
run hom --event "del:chr20:38900000-38910000;af=1.0"
args=(); while read -r e; do args+=(--event "$e"); done < /tmp/claude-1066/-home-parlar-ai-dev-spike/23a35c29-9600-4da2-a2c9-618425a1e0d3/scratchpad/pres/snv_events.txt
run snv "${args[@]}"
