#!/usr/bin/env bash
# From the runs of runs.sh to the files in ../data, then to the figures in ../figures.
# Run runs.sh first (same REF, SLICE and SPIKE). Needs samtools and Python 3 with numpy, matplotlib and pysam;
# the charts use the Noto Sans font.
set -euo pipefail
cd "$(dirname "$0")"
REF="${REF:?set REF to the GRCh38 no-alt FASTA}"
SLICE="${SLICE:?set SLICE to the chr20 slice BAM}"
S="${SPIKE:-spike}"

# Figure 8: depth across the 10 kb deletion, 250 bp bins
R=chr20:38885001-38925000
samtools depth -a -r $R "$SLICE" > depth_orig.txt
samtools depth -a -r $R het/merged.bam > depth_het.txt
samtools depth -a -r $R hom/merged.bam > depth_hom.txt
python3 depth_bins.py

# Reads in the 25 SNV windows (5 kb each): Figures 5, 6, 10-13 (Figure 5 by runbins.py, which needs pysam)
regs=$(for i in $(seq 0 24); do p=$((38550000 + i*60000)); echo "chr20:$((p-2500))-$((p+2500))"; done | tr '\n' ' ')
samtools view -F 0xF0C snv/merged.bam $regs > snv_reads.sam
python3 reads_insert_quality.py      # insert.csv, qual_cycle.csv
python3 qual_stats.py                # qual_stats.json: transitions, low-quality tails, mismatch by Q
python3 qvar.py                      # qual_var.csv, qual_var.json: quality spread
python3 qvar_ctrl.py                 # controls for the spread (printed)

# Figures 7 and 9: allele fractions, recounted independently of spike validate
$S validate --bam snv/merged.bam --truth snv/truth.vcf --reference "$REF" --json > snv.validate.json 2>/dev/null || true
python3 snv_recount.py               # snv_af_counts.csv
python3 stats2.py                    # snv_accept.csv; Wilson CIs and the KS test (printed)

python3 charts.py && python3 charts2.py && python3 charts3.py
python3 runbins.py                   # crash_runbin.csv, chart_runbin.png: Figure 5
mkdir -p ../data ../figures
mv -f snv.validate.json depth_bins.csv insert.csv qual_cycle.csv qual_stats.json qual_var.csv qual_var.json snv_af_counts.csv snv_accept.csv crash_runbin.csv ../data/
mv -f chart_*.png ../figures/
