#!/bin/bash
# Check H of the --into-fastq plan (docs/superpowers/plans/2026-10-04-full-fastq.md).
set -uo pipefail
S="$(cd "$(dirname "$0")" && pwd)"
SRC=~/dev/cnv_validation/from_rv/raredisease_seracare_na12878_resample/D24-14230_Seq25-7600_30x_resample/raredisease_results/alignment/D24-14230_Seq25-7600_30x_sorted_md.bam
REF=~/dev/spike/data/giab_hg38/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta
STANDIN="$S/../fastq/standin/main"
RG='@RG\tID:sim\tPL:ILLUMINA\tSM:sim'
OLD="$S/run-fix/spike/out"
run() {  # DIR ARGS...: the after-check's spike step with ARGS in place of --raw-fastq
    local d=$1; shift
    rm -rf "$S/h/$d"; mkdir -p "$S/h/$d"
    (cd "$S/h/$d" && "$S/spike-if" --bam "$SRC" --reference "$REF" --vcf "$S/run-fix/events.vcf" --seed 1 --threads 16 \
        --aligner "bwa-mem2 mem -M -K 100000000 -t 16 -R '$RG'" --align "$@" -o out > spike.log 2>&1)
    echo "$d: exit $?"
}
same() { if cmp -s <(pigz -dc "$1") <(pigz -dc "$2"); then echo "PASS  $3 identical (decompressed)"; else echo "FAIL  $3 differs"; fi; }
samefile() { if cmp -s "$1" "$2"; then echo "PASS  $3 identical"; else echo "FAIL  $3 differs"; fi; }

echo "== H1"
run h1 --into-fastq "$STANDIN/raw_R1.fastq.gz" "$STANDIN/raw_R2.fastq.gz"
for m in 1 2; do same "$S/h/h1/out/spiked_R$m.fastq.gz" "$OLD/spiked_R$m.fastq.gz" "spiked_R$m.fastq.gz"; done
for m in 1 2; do same "$S/h/h1/out/R$m.fq.gz" "$OLD/R$m.fq.gz" "R$m.fq.gz"; done
for f in replaced_reads.txt fastq_removed_reads.txt truth.vcf; do samefile "$S/h/h1/out/$f" "$OLD/$f" "$f"; done
echo "control (spiked_R1 against the raw stand-in R1, must differ):"; same "$S/h/h1/out/spiked_R1.fastq.gz" "$STANDIN/raw_R1.fastq.gz" "spiked_R1 vs raw R1"

echo "== H2"
run h2 --raw-fastq "$STANDIN/raw_R1.fastq.gz" "$STANDIN/raw_R2.fastq.gz"
[ -e "$S/h/h2/out" ] && echo "FAIL  an output folder was made" || echo "PASS  no output folder"
grep -c "Extracting" "$S/h/h2/spike.log" | sed 's/^/Extracting lines (want 0): /'
echo "clap says:"; cat "$S/h/h2/spike.log"

echo "== H3"
echo "last log line:"; grep -v "^\s*$" "$S/h/h1/spike.log" | tail -1
grep -n "R1.fq.gz\|spiked_R1" "$S/h/h1/out/README.md"

echo "== H4"
run h4 --into-fastq "$STANDIN/raw_R1.fastq.gz" "$STANDIN/raw_R2.fastq.gz" --fastq-prefix S1
for m in 1 2; do same "$S/h/h4/out/S1_R$m.fastq.gz" "$OLD/spiked_R$m.fastq.gz" "S1_R$m.fastq.gz against the after-check's spiked_R$m"; done
ls "$S/h/h4/out" | grep -c "^spiked_" | sed 's/^/spiked_ files (want 0): /'
echo "last log line:"; grep -v "^\s*$" "$S/h/h4/spike.log" | tail -1
grep -n "S1_R1\|spiked_R1" "$S/h/h4/out/README.md"
run h4slash --into-fastq "$STANDIN/raw_R1.fastq.gz" "$STANDIN/raw_R2.fastq.gz" --fastq-prefix a/b
[ -e "$S/h/h4slash/out" ] && echo "FAIL  an output folder was made" || echo "PASS  no output folder"; grep -c "Extracting" "$S/h/h4slash/spike.log" | sed 's/^/Extracting lines (want 0): /'; tail -1 "$S/h/h4slash/spike.log"
run h4alone --fastq-prefix S1
[ -e "$S/h/h4alone/out" ] && echo "FAIL  an output folder was made" || echo "PASS  no output folder"; grep -c "Extracting" "$S/h/h4alone/spike.log" | sed 's/^/Extracting lines (want 0): /'; cat "$S/h/h4alone/spike.log"
