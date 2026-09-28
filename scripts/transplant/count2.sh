#!/usr/bin/env bash
# count2.sh: the inputs of transplant round 2 (SNVs, 1-49 bp indels, 50-299 bp duplications).
#   small.bed      Platinum smallvar region AND Q100 v5.0q smvar region, autosomes, minus masks
#   h1.all.tsv     every HG001 small-variant record, split and left-normalised: CHROM POS REF ALT GT
#   h2.all.tsv     the same for HG002
#   sv.bed         round 1's SV region (Q100 stvar AND Platinum svs, autosomes, minus masks)
#   insbench/      truvari bench of the two samples' insertions >= 50 bp inside sv.bed
# Every tool's stderr goes to a log beside its output; the last line says the script finished.
set -euo pipefail
D=${1:?usage: count2.sh OUT_DIR}; mkdir -p "$D"
CV=${CNV_VALIDATION:-$HOME/dev/cnv_validation}
TRUVARI=$CV/.pixi/envs/default/bin/truvari
BEDTOOLS=$CV/.pixi/envs/default/bin/bedtools
REF=$HOME/dev/spike/data/giab_hg38/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta
H1=$CV/platinum_pedigree_truthset_v1.2/NA12878_hq_v1.2.1.latest-smallvar.vcf.gz
H1BED=$CV/platinum_pedigree_truthset_v1.2/NA12878_hq_v1.2.smallvar.bed.gz
H1SV=$CV/platinum_pedigree_truthset_v1.2/NA12878_hq_v1.2.svs.vcf.gz
H1SVBED=$CV/platinum_pedigree_truthset_v1.2/NA12878_hq_v1.2.svs.bed.gz
H2=$HOME/dev/spike/data/giab_hg38/HG002/GRCh38_HG2-T2TQ100-V1.1.vcf.gz
H2BED=$HOME/dev/spike/data/giab_hg38/HG002/HG002_GRCh38_v5.0q_smvar.benchmark.bed
H2SVBED=$HOME/dev/spike/data/giab_hg38/HG002/GRCh38_HG2-T2TQ100-V1.1_stvar.benchmark.bed
AUTO=$(seq -s, -f 'chr%g' 1 22)
cd "$D"

# Masks: the 7 SeraCare genes +-2 kb, and +-2 kb around both ends of the cell line's MSH2-intron <-> IGL junctions.
zcat $CV/resources/genes.bed.gz | LC_ALL=C awk -F'\t' -v OFS='\t' '$4 ~ /^(BRCA1|BRCA2|MSH2|MSH6|MLH1|PMS2|CDKN2A)$/ {s=$2-2000; if (s<0) s=0; print $1, s, $3+2000}' > mask.raw
printf 'chr2\t47580850\t47584865\nchr22\t23851678\t23858364\n' >> mask.raw
sort -k1,1 -k2,2n mask.raw | $BEDTOOLS merge > mask.bed

region() {  # region A B OUT: A AND B, autosomes, minus the masks
  $BEDTOOLS intersect -a "$1" -b "$2" | LC_ALL=C awk '$1 ~ /^chr([1-9]|1[0-9]|2[0-2])$/' | sort -k1,1 -k2,2n \
    | $BEDTOOLS merge | $BEDTOOLS subtract -a - -b mask.bed > "$3"
}
cut -f1-3 $H2BED | sort -k1,1 -k2,2n > h2.small.bed
zcat $H1BED | cut -f1-3 | sort -k1,1 -k2,2n > h1.small.bed
region h2.small.bed h1.small.bed small.bed
cut -f1-3 $H2SVBED | sort -k1,1 -k2,2n > h2.sv.bed
zcat $H1SVBED | cut -f1-3 | sort -k1,1 -k2,2n > h1.sv.bed
region h2.sv.bed h1.sv.bed sv.bed

# Small-variant records. -c x drops a record whose REF does not match the reference (Q100 has one, at chr9:70701801).
for S in h1:$H1 h2:$H2; do N=${S%%:*}; V=${S#*:}
  bcftools view -r $AUTO $V -Ou 2> $N.view.log | bcftools norm -m-any -c x -f $REF -Ou 2> $N.norm.log \
    | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT]\n' > $N.all.tsv 2> $N.query.log
  test "$(cut -f1 $N.all.tsv | sort -u | wc -l)" -eq 22 || { echo "$N.all.tsv does not hold 22 autosomes" >&2; exit 1; }
done

# Insertions >= 50 bp, biallelic, on autosomes, matched between the samples by truvari.
for S in h1:$H1SV h2:$H2; do N=${S%%:*}; V=${S#*:}
  bcftools view -m2 -M2 -i '(strlen(ALT)-strlen(REF))>=50' -r $AUTO $V -Oz -o $N.ins.vcf.gz 2> $N.ins.log
  bcftools index -t -f $N.ins.vcf.gz
done
rm -rf insbench
$TRUVARI bench -b h1.ins.vcf.gz -c h2.ins.vcf.gz --includebed sv.bed --sizemin 50 --sizemax 100000000 --passonly \
  -o insbench > insbench.log 2>&1 || { tail -20 insbench.log >&2; exit 1; }

for F in small.bed sv.bed; do echo "$F: $(LC_ALL=C awk '{s+=$3-$2} END {printf "%.0f Mb", s/1e6}' $F)"; done
echo "h1.all.tsv: $(wc -l < h1.all.tsv) records; h2.all.tsv: $(wc -l < h2.all.tsv) records"
echo "count2.sh: done"
