DUP Validation Test - 2026-03-03
=================================

Goal: Validate spike's DUP simulation by comparing synthetic tandem duplications
(spiked into HG002) against real tandem duplications (observed in HG001/NA12878).

Source of DUP variants
----------------------
- File: data/giab_hg38/HG001/HG001_GRCh38.pbsv.vcf.gz
  (PacBio pbsv callset from GIAB FTP: ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/NA12878/analysis/PacBio_pbsv_05212019/)
- GIAB benchmarks contain no DUP truth variants. However, pbsv encodes tandem
  duplications as INS with SVANN=TANDEM annotation (17,750 of 28,090 INS calls).
- We filtered for PASS, SVANN=TANDEM, size 200-5000bp, overlapping HG001 v4.2.1
  high-confidence BED regions, yielding 18 candidates.
- Selected 5 across different chromosomes and sizes.

Conversion: tandem INS -> DUP
-----------------------------
For each tandem INS at VCF POS with SVLEN, the duplicated region is [POS, POS+SVLEN).
The ALT sequence of these tandem insertions matches the reference immediately after
the insertion point, confirming they are tandem duplications of the downstream region.
Converted to SVTYPE=DUP with END=POS+SVLEN for spike input.

Input VCF: data/giab_hg38/HG001/hg001_tandem_dups_for_spike.vcf

Selected DUP events
--------------------
  ID       Chr    POS          SVLEN   DUP region
  dup_01   chr3   184754689    1472    chr3:184754689-184756161
  dup_02   chr4   26209775     2709    chr4:26209775-26212484
  dup_03   chr6   28712877     1876    chr6:28712877-28714753
  dup_04   chr8   992901       1890    chr8:992901-994791
  dup_05   chr20  819293       3287    chr20:819293-822580

All are heterozygous (0/1) in HG001.

Spike command
-------------
  target/release/spike \
    --bam data/giab_hg38/HG002/HG002.novaseq.pcr-free.35x.bwamem2.dedup.grch38_no_alt.bam \
    --reference data/giab_hg38/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta \
    --vcf data/giab_hg38/HG001/hg001_tandem_dups_for_spike.vcf \
    --allele-fraction 0.5 \
    --output data/validation/dup_test \
    --align

Output files
------------
  sim.bam       - Aligned BAM with spiked DUP reads (original + synthetic)
  sim.bam.bai   - BAM index
  truth.vcf     - Truth VCF with 5 DUP events
  R1.fq.gz      - Forward reads (spiked FASTQ)
  R2.fq.gz      - Reverse reads (spiked FASTQ)
  align.sh      - Alignment script used
  align.log     - bwa-mem2 alignment log

IGV comparison
--------------
To validate realism, compare side-by-side in IGV at the DUP loci above:
  Track 1: sim.bam (synthetic DUPs in HG002 background)
  Track 2: HG001 BAM (real tandem DUPs) - use NA12878.final.cram from
           1000 Genomes (ftp.sra.ebi.ac.uk/vol1/run/ERR323/ERR3239334/)

What to look for:
  - Depth increase (~1.5x for het DUP) within duplicated region
  - Junction/chimeric reads at DUP boundaries
  - Allelic imbalance at het SNPs within the DUP region (~33/67 ratio)
