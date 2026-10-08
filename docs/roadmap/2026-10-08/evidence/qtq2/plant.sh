#!/bin/bash
# plant.sh ID CHROM START END -- one event: spike (current master) into HG002, the real
# hospital aligner command, and three small BAMs over the event +/- 20 kb:
#   spiked.bam (HG002 minus spike's replaced reads, plus sim.bam), base.bam (HG002 as is),
#   real.bam (HG001 at the same place; for 'forward' events HG001 carries it for real).
set -uo pipefail
id=$1; chrom=$2; start=$3; end=$4
D=/tmp/claude-1066/-home-parlar-ai-dev-spike/23a35c29-9600-4da2-a2c9-618425a1e0d3/scratchpad/qtq2/k2v2/$id
H2=/home/parlar_ai/dev/cnv_validation/from_rv/raredisease_seracare_na12878_resample/D24-14230_Seq25-7600_30x_resample/raredisease_results/alignment/D24-14230_Seq25-7600_30x_sorted_md.bam
H1=/home/parlar_ai/dev/cnv_validation/from_rv/raredisease_seracare_na12878_resample/D25-7403_Seq25-7598_30x_resample/raredisease_results/alignment/D25-7403_Seq25-7598_30x_sorted_md.bam
REF=/home/parlar_ai/dev/spike/data/giab_hg38/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta
SPIKE=/tmp/claude-1066/-home-parlar-ai-dev-spike/23a35c29-9600-4da2-a2c9-618425a1e0d3/scratchpad/qtq2/spike-v2
PX="pixi run --manifest-path /home/parlar_ai/sv_slice_test/wt/slice_safety/pixi.toml"
ALIGNER="bwa-mem2 mem -M -K 100000000 -t 4 -R '@RG\tID:sim\tPL:ILLUMINA\tSM:sim'"
rm -rf "$D"; mkdir -p "$D"
win="$chrom:$((start - 20000))-$((end + 20000))"
echo "$win" > "$D/window"
$PX "$SPIKE" --bam "$H2" --reference "$REF" --event "del:$chrom:$start-$end;af=het" --seed 1 --threads 4 \
    --aligner "$ALIGNER" -o "$D/run" > "$D/spike.log" 2>&1 || { echo "$id refused: $(grep -m1 Error "$D/spike.log")"; exit 0; }
$PX bash "$D/run/align.sh" "$REF" 4 > "$D/align.log" 2>&1 || { echo "$id align failed"; exit 1; }
$PX samtools view -b -o "$D/base.bam" "$H2" "$win" && $PX samtools index "$D/base.bam"
$PX samtools view -b -o "$D/real.bam" "$H1" "$win" && $PX samtools index "$D/real.bam"
$PX samtools view -b -N "$D/run/replaced_reads.txt" -U "$D/kept.bam" -o /dev/null "$D/base.bam"
$PX samtools merge -f -o "$D/merged.tmp.bam" "$D/kept.bam" "$D/run/sim.bam"
$PX samtools sort -o "$D/spiked.bam" "$D/merged.tmp.bam" && $PX samtools index "$D/spiked.bam"
rm -f "$D/merged.tmp.bam" "$D/kept.bam"
echo "$id ok $(grep -v '^#' "$D/run/truth.vcf" | cut -f8 | grep -o 'SIM_ALT_FRAGS=[0-9]*;SIM_RESIST=[0-9.]*')"
