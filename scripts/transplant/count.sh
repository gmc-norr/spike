#!/usr/bin/env bash
# count.sh: deletions >= 50 bp usable for an HG001 <-> HG002 transplant test.
set -euo pipefail
D=${1:?usage: count.sh OUT_DIR}; mkdir -p "$D"
CV=${CNV_VALIDATION:-$HOME/dev/cnv_validation}
TRUVARI=$CV/.pixi/envs/default/bin/truvari
BEDTOOLS=$CV/.pixi/envs/default/bin/bedtools
H2=$HOME/dev/spike/data/giab_hg38/HG002/GRCh38_HG2-T2TQ100-V1.1.vcf.gz
H2BED=$HOME/dev/spike/data/giab_hg38/HG002/GRCh38_HG2-T2TQ100-V1.1_stvar.benchmark.bed
H1=$CV/platinum_pedigree_truthset_v1.2/NA12878_hq_v1.2.svs.vcf.gz
H1BED=$CV/platinum_pedigree_truthset_v1.2/NA12878_hq_v1.2.svs.bed.gz
FAI=$HOME/dev/spike/data/giab_hg38/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta.fai
cd $D

# Masks: the 7 SeraCare genes +-2 kb, and +-2 kb around both ends of the cell line's MSH2-intron <-> IGL junctions.
zcat $CV/resources/genes.bed.gz | LC_ALL=C awk -F'\t' -v OFS='\t' '$4 ~ /^(BRCA1|BRCA2|MSH2|MSH6|MLH1|PMS2|CDKN2A)$/ {s=$2-2000; if (s<0) s=0; print $1, s, $3+2000}' > mask.raw
printf 'chr2\t47580850\t47584865\nchr22\t23851678\t23858364\n' >> mask.raw
sort -k1,1 -k2,2n mask.raw | $BEDTOOLS merge > mask.bed

# Region both truth sets vouch for, autosomes only, minus the masks.
( if [[ $H2BED == *.gz ]]; then zcat $H2BED; else cat $H2BED; fi ) | cut -f1-3 | sort -k1,1 -k2,2n > h2.bed
zcat $H1BED | cut -f1-3 | sort -k1,1 -k2,2n > h1.bed
$BEDTOOLS intersect -a h2.bed -b h1.bed | LC_ALL=C awk '$1 ~ /^chr([1-9]|1[0-9]|2[0-2])$/' | sort -k1,1 -k2,2n | $BEDTOOLS merge | $BEDTOOLS subtract -a - -b mask.bed > common.bed

# Deletions >= 50 bp (REF longer than ALT by >= 50), biallelic, on autosomes.
for S in h1:$H1 h2:$H2; do N=${S%%:*}; V=${S#*:}
  bcftools view -m2 -M2 -i '(strlen(REF)-strlen(ALT))>=50' -r $(seq -s, -f 'chr%g' 1 22) $V -Oz -o $N.del.vcf.gz 2>/dev/null || \
  bcftools view -m2 -M2 -i '(strlen(REF)-strlen(ALT))>=50' $V -Oz -o $N.del.vcf.gz
  bcftools index -t -f $N.del.vcf.gz
done

rm -rf bench
$TRUVARI bench -b h1.del.vcf.gz -c h2.del.vcf.gz --includebed common.bed --sizemin 50 --sizemax 100000000 --passonly -o bench > bench.log 2>&1 || { tail -20 bench.log; exit 1; }

# Shared: HG001 deletions matched in HG002 (tp-base). HG001-only: fn. HG002-only: fp.
bin() { LC_ALL=C awk -F'\t' '{ l=length($4)-length($5); split($10,g,":"); gt=g[1]; gsub(/\|/,"/",gt); h=(gt=="0/1"||gt=="1/0")?"het":((gt=="1/1")?"hom":"other"); b=(l<300)?"50-299":(l<1000)?"300-999":(l<10000)?"1k-10k":"10k+"; c[b" "h]++ } END { for (k in c) print k, c[k] }' | sort; }
for F in tp-base fn fp; do echo "== $F"; bcftools view -H bench/$F.vcf.gz | bin; done
echo "common region: $(LC_ALL=C awk '{s+=$3-$2} END {printf "%.0f Mb", s/1e6}' common.bed)"
