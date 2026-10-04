#!/bin/bash
# The duplicates check (docs/superpowers/plans/2026-10-04-duplicates.md): spike a full FASTQ with
# --into-fastq, then align and mark duplicates the way raredisease does, beside the unspiked stand-in
# and C1's stand-in. Every step's stderr goes to a log next to its output.
#
# Usage: run.sh SPIKE SOURCE_BAM REFERENCE CALLS_VCF STANDIN_DIR PICARD_JAR OUT
#   STANDIN_DIR holds raw_R1.fastq.gz and raw_R2.fastq.gz (the full-FASTQ check's stand-in).
set -euo pipefail
SPIKE=$1 SRC=$2 REF=$3 CALLS=$4 STANDIN=$5 PICARD=$6 OUT=$7
HERE="$(cd "$(dirname "$0")" && pwd)"
REGION=chr20:9500000-12500000
SOURCE_REGION=chr20:9000000-13000000
RG='@RG\tID:sim\tPL:ILLUMINA\tSM:sim'
mkdir -p "$OUT/tmp"

align() {  # R1 R2 OUT_BAM: raredisease's bwa-mem2 command, then sorted
    bwa-mem2 mem -M -K 100000000 -t 16 -R "$RG" "$REF" "$1" "$2" 2> "$3.bwa.err" \
        | samtools sort -@ 4 -o "$3" - 2> "$3.sort.err"
    samtools index "$3"
}

mark() {  # IN OUT [EXTRA...]: Picard 3.3.0 MarkDuplicates with raredisease's options
    local in=$1 out=$2
    shift 2
    java -Xmx8g -jar "$PICARD" MarkDuplicates --INPUT "$in" --OUTPUT "$out" \
        --METRICS_FILE "${out%.bam}.metrics.txt" --TMP_DIR "$OUT/tmp" --REFERENCE_SEQUENCE "$REF" \
        --MAX_SEQUENCES_FOR_DISK_READ_ENDS_MAP 50000 "$@" 2> "${out%.bam}.picard.err"
    samtools index "$out"
}

records() { echo $(( $(pigz -dc "$1" | wc -l) / 4 )); }

python3 -B "$HERE/draw_events.py" --bam "$SRC" --reference "$REF" --calls "$CALLS" --region "$REGION" \
    --seed 7 --out "$OUT/events.vcf" > "$OUT/events.txt"
echo "events: $(grep -vc '^#' "$OUT/events.vcf")"

mkdir -p "$OUT/spike"
(cd "$OUT/spike" && "$SPIKE" --bam "$SRC" --reference "$REF" --vcf "$OUT/events.vcf" --seed 1 --threads 16 \
    --aligner "bwa-mem2 mem -M -K 100000000 -t 16 -R '$RG'" --align \
    --into-fastq "$STANDIN/raw_R1.fastq.gz" "$STANDIN/raw_R2.fastq.gz" -o out > spike.log 2>&1)
echo "spike: done"

align "$STANDIN/raw_R1.fastq.gz" "$STANDIN/raw_R2.fastq.gz" "$OUT/baseline.sorted.bam"
mark "$OUT/baseline.sorted.bam" "$OUT/baseline.bam"
align "$OUT/spike/out/spiked_R1.fastq.gz" "$OUT/spike/out/spiked_R2.fastq.gz" "$OUT/spiked.sorted.bam"
mark "$OUT/spiked.sorted.bam" "$OUT/spiked.bam"
echo "baseline and spiked: aligned and marked"

# C1: drop the kept pair of 100 duplicate sets, then mark again, tagging on.
mark "$OUT/baseline.sorted.bam" "$OUT/baseline.tagged.bam" --TAG_DUPLICATE_SET_MEMBERS true
python3 -B "$HERE/control.py" pick --bam "$OUT/baseline.tagged.bam" --region "$SOURCE_REGION" --seed 7 \
    --sets "$OUT/c1_sets.tsv" --drop "$OUT/c1_drop.txt"
for m in 1 2; do
    pigz -dc "$STANDIN/raw_R$m.fastq.gz" \
        | DROP="$OUT/c1_drop.txt" LC_ALL=C awk 'BEGIN { while ((getline n < ENVIRON["DROP"]) > 0) drop[n] = 1 }
            NR % 4 == 1 { name = substr($1, 2); sub(/\/[12]$/, "", name); skip = (name in drop) }
            !skip' \
        | pigz -p 8 > "$OUT/c1_R$m.fastq.gz"
    before=$(records "$STANDIN/raw_R$m.fastq.gz") after=$(records "$OUT/c1_R$m.fastq.gz")
    echo "C1 R$m: $before records, $after after dropping $(wc -l < "$OUT/c1_drop.txt")"
    [ $(( before - after )) -eq "$(wc -l < "$OUT/c1_drop.txt")" ]
done
align "$OUT/c1_R1.fastq.gz" "$OUT/c1_R2.fastq.gz" "$OUT/c1.sorted.bam"
mark "$OUT/c1.sorted.bam" "$OUT/c1.tagged.bam" --TAG_DUPLICATE_SET_MEMBERS true
python3 -B "$HERE/control.py" count --bam "$OUT/c1.tagged.bam" --region "$SOURCE_REGION" \
    --sets "$OUT/c1_sets.tsv" > "$OUT/c1.txt"
cat "$OUT/c1.txt"

python3 -B "$HERE/measure.py" --spiked "$OUT/spiked.bam" --baseline "$OUT/baseline.bam" --source "$SRC" \
    --sim "$OUT/spike/out/sim.bam" --reference "$REF" --calls "$CALLS" --events "$OUT/events.vcf" \
    --replaced "$OUT/spike/out/replaced_reads.txt" --source-region "$SOURCE_REGION" --region "$REGION" \
    --spiked-metrics "$OUT/spiked.metrics.txt" --baseline-metrics "$OUT/baseline.metrics.txt" \
    --c1 "$OUT/c1.txt" | tee "$OUT/report.txt"
