#!/bin/bash
# K3 of the duplicates fix: --edit-model origin, old and new binaries, same events and seed.
set -uo pipefail
S="$(cd "$(dirname "$0")" && pwd)"
SRC=~/dev/cnv_validation/from_rv/raredisease_seracare_na12878_resample/D24-14230_Seq25-7600_30x_resample/raredisease_results/alignment/D24-14230_Seq25-7600_30x_sorted_md.bam
REF=~/dev/spike/data/giab_hg38/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta
for b in spike-dp spike-fix; do
    rm -rf "$S/k3/$b"; mkdir -p "$S/k3/$b"
    (cd "$S/k3/$b" && /usr/bin/time -f "wall %e s" "$S/$b" --bam "$SRC" --reference "$REF" --vcf "$S/run/events.vcf" \
        --seed 1 --threads 16 --edit-model origin -o out > spike.log 2>&1)
    echo "$b: exit $? ($(tail -1 "$S/k3/$b/spike.log"))"
done
ok=1
for f in replaced_reads.txt fastq_removed_reads.txt truth.vcf; do
    if cmp -s "$S/k3/spike-dp/out/$f" "$S/k3/spike-fix/out/$f"; then echo "PASS  $f identical ($(wc -l < "$S/k3/spike-dp/out/$f") lines)"; else echo "FAIL  $f differs"; ok=0; fi
done
for f in R1.fq.gz R2.fq.gz; do
    if cmp -s <(pigz -dc "$S/k3/spike-dp/out/$f") <(pigz -dc "$S/k3/spike-fix/out/$f"); then echo "PASS  $f identical decompressed ($(( $(pigz -dc "$S/k3/spike-dp/out/$f" | wc -l) / 4 )) records)"; else echo "FAIL  $f differs"; ok=0; fi
done
grep -c "Duplicates of the removed originals" "$S/k3/spike-fix/spike.log" | sed 's/^/clean-mode duplicate log lines in the origin run (want 0): /'
[ "$ok" = 1 ] && echo "K3: PASS" || echo "K3: FAIL"
