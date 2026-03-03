#!/bin/bash
set -euo pipefail
# Align simulated reads and sort.
# Usage: bash align.sh [REF] [THREADS]
REF="${1:-data/giab_hg38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta}"
THREADS="${2:-8}"
SAMTOOLS="samtools"
DIR="$(cd "$(dirname "$0")" && pwd)"

echo "Aligning $DIR/R1.fq.gz + R2.fq.gz (bwa-mem2, $THREADS threads)..."
bwa-mem2 mem -t "$THREADS" \
    -R '@RG\tID:sim\tSM:SIM\tPL:ILLUMINA' \
    "$REF" \
    "$DIR/R1.fq.gz" "$DIR/R2.fq.gz" \
    2>"$DIR/align.log" | \
    "$SAMTOOLS" sort -@ "$THREADS" -o "$DIR/sim.bam" -

"$SAMTOOLS" index "$DIR/sim.bam"

TOTAL=$("$SAMTOOLS" view -c "$DIR/sim.bam" 2>/dev/null || echo "?")
SA_COUNT=$("$SAMTOOLS" view "$DIR/sim.bam" | grep -c "SA:Z:" 2>/dev/null || echo "0")
echo "Done: $DIR/sim.bam ($TOTAL reads, $SA_COUNT with SA tags)"
