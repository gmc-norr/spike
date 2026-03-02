# spike

Haplotype-based read spike-in simulator for genomic variants.

spike takes a real BAM/CRAM file and a set of variant specifications, then produces paired FASTQ files containing reads that carry the specified variants at a controlled allele fraction. Unlike simple read-cloning approaches, spike generates independent synthetic reads with realistic quality profiles that survive deduplication tools.

Designed for **validating bioinformatic pipelines** — variant callers, structural variant detection, fusion detection, copy number analysis, and any downstream tool that operates on aligned sequencing data.

## How it works

```
Real BAM/CRAM + Reference FASTA + Variant specs
                  |
                  v
    +-----------------------------+
    |  1. Extract read pairs      |  Real reads from the event region
    |  2. Learn quality           |  Per-cycle Q score distributions from real data
    |  3. Build haplotype         |  Linear variant sequence from ordered segments
    |  4. Suppress reads          |  Remove reads at VAF rate within SV boundaries
    |  5. Tile synthetic          |  New reads across haplotype with learned Q profile
    |  6. LOH simulation          |  Het SNPs -> hom for DELs; BAF shift for DUPs
    +-----------------------------+
                  |
                  v
     Paired FASTQ + Truth VCF + Alignment script
```

**Variant haplotype model**: The variant allele is represented as an ordered list of segments, each from a reference region (possibly reverse-complemented) or novel sequence. Reads tiled uniformly across this linear sequence are automatically chimeric when they span segment boundaries — no per-SV-type breakpoint logic needed.

**Quality-aware synthesis**: Instead of cloning real reads (which produces exact duplicates flagged by dedup tools), spike learns per-cycle, base-conditioned quality score distributions from the donor reads and generates independent synthetic reads with correlated sequencing errors.

**LOH and allelic imbalance**: For heterozygous deletions, het SNPs within the deleted region become homozygous (only the surviving haplotype allele remains). For duplications, het SNPs shift from ~50/50 to ~33/67 allele balance. Het SNP positions are found via pileup (default) or from a pre-called gVCF.

## Supported variant types

| Type | Event spec | VCF SVTYPE | Description |
|------|-----------|------------|-------------|
| Deletion | `del:chr:start-end` or `del:GENE:exon4-exon8` | DEL | Region removed from one haplotype |
| Duplication | `dup:chr:start-end` or `dup:GENE:exon4-exon8` | DUP | Tandem duplication in place |
| Inversion | `inv:chr:start-end` or `inv:GENE:exon4-exon8` | INV | Region reversed in place |
| Insertion | `ins:chr:pos:length` or `ins:chr:pos:ACGT` | INS | Novel sequence inserted at position |
| Fusion | `fusion:GENEA:exonN:GENEB:exonM` | BND | Two breakpoints joined across genes/chromosomes |
| Inverted Fusion | `fusion:GENEA:exonN:GENEB:exonM:inv` | BND | Fusion with reverse-complement join |
| SNP/Indel | `snp:chr:pos:REF:ALT` | Standard REF/ALT | SNPs, MNVs, small insertions/deletions |

All types support per-event allele fraction control.

## Installation

### Prerequisites

- **Rust** (edition 2021 or later)
- **Reference FASTA** with `.fai` index (e.g., from `samtools faidx`)
- **BAM/CRAM** input must be coordinate-sorted and indexed (`.bai` / `.crai`)
- **bcftools** (only needed if using `--gvcf` with `.vcf.gz` files)
- An **aligner** for the optional `--align` step (default: `bwa-mem2`; also supports `minimap2`, `bowtie2`, or any custom aligner)

### Build

```bash
cargo build --release
# Binary at target/release/spike
```

## Quick start

```bash
# Simulate a 10kb heterozygous deletion on chr17
spike \
  --bam sample.bam \
  --reference GRCh38.fasta \
  --event "del:chr17:43045000-43055000" \
  -o output/

# Output: output/R1.fq.gz, output/R2.fq.gz, output/truth.vcf, output/align.sh
```

## Usage examples

### Coordinate-based structural variants

Specify events directly with genomic coordinates. For coordinate-based `del`, `dup`, and `inv` specs, the coordinates follow the VCF convention where the start value is equivalent to the 0-based SV start:

```bash
# 10kb deletion
spike --bam sample.bam --reference GRCh38.fasta \
  --event "del:chr17:43045000-43055000" \
  -o output/

# 5kb tandem duplication
spike --bam sample.bam --reference GRCh38.fasta \
  --event "dup:chr20:30000000-30005000" \
  -o output/

# 3kb inversion
spike --bam sample.bam --reference GRCh38.fasta \
  --event "inv:chr9:21970000-21973000" \
  -o output/

# 500bp insertion (random sequence)
spike --bam sample.bam --reference GRCh38.fasta \
  --event "ins:chr20:30000000:500" \
  -o output/

# Insertion with explicit sequence
spike --bam sample.bam --reference GRCh38.fasta \
  --event "ins:chr20:30000000:ACGTACGTACGT" \
  -o output/
```

### Gene/exon-based events

When an exon BED file is provided, events can be specified using gene names and exon ranges. This is useful for simulating clinically relevant variants like multi-exon deletions:

```bash
# Delete exons 4-8 of LDLR
spike --bam sample.bam --reference GRCh38.fasta \
  --exon-bed gene_exons.bed \
  --event "del:LDLR:exon4-exon8" \
  -o output/

# Duplicate exons 2-5 of BRCA1
spike --bam sample.bam --reference GRCh38.fasta \
  --exon-bed gene_exons.bed \
  --event "dup:BRCA1:exon2-exon5" \
  -o output/

# Invert exons 3-6 of a gene
spike --bam sample.bam --reference GRCh38.fasta \
  --exon-bed gene_exons.bed \
  --event "inv:TP53:exon3-exon6" \
  -o output/
```

The exon BED file should be tab-separated with at least 4 columns: `chrom start end name [gene]`. If the 5th column (gene) is absent, the gene is parsed from the name (e.g., `LDLR_exon1` -> `LDLR`).

### Gene fusions

Fusions join two breakpoints, potentially across different chromosomes. The `:inv` suffix creates an inverted fusion (reverse-complement join at gene B):

```bash
# BCR-ABL1 fusion (forward)
spike --bam sample.bam --reference GRCh38.fasta \
  --exon-bed gene_exons.bed \
  --event "fusion:BCR:exon14:ABL1:exon2" \
  -o output/

# EML4-ALK inverted fusion
spike --bam sample.bam --reference GRCh38.fasta \
  --exon-bed gene_exons.bed \
  --event "fusion:EML4:exon13:ALK:exon20:inv" \
  -o output/

# Fusion at low somatic VAF
spike --bam sample.bam --reference GRCh38.fasta \
  --exon-bed gene_exons.bed \
  --event "fusion:BCR:exon14:ABL1:exon2;af=0.05" \
  -o output/
```

### SNPs and small indels

Small variants use the `snp:` prefix (applies to SNPs, MNVs, small deletions, and small insertions). Positions are 1-based:

```bash
# Single nucleotide variant
spike --bam sample.bam --reference GRCh38.fasta \
  --event "snp:chr17:7577120:C:T" \
  -o output/

# Alternative syntax with '>'
spike --bam sample.bam --reference GRCh38.fasta \
  --event "snp:chr17:7577120:C>T" \
  -o output/

# Small deletion (3bp -> 1bp)
spike --bam sample.bam --reference GRCh38.fasta \
  --event "snp:chr17:7577530:ACG:A" \
  -o output/

# Small insertion (1bp -> 4bp)
spike --bam sample.bam --reference GRCh38.fasta \
  --event "snp:chr1:100000:A:ACGT" \
  -o output/
```

### Multiple events with per-event allele fractions

Combine multiple events in a single run. Each event can have its own allele fraction:

```bash
spike --bam sample.bam --reference GRCh38.fasta \
  --exon-bed gene_exons.bed \
  --event "del:chr17:43045000-43055000;af=0.3" \
  --event "dup:chr20:30000000-30005000;af=het" \
  --event "snp:chr17:7577120:C:T;af=0.15" \
  --event "fusion:BCR:exon14:ABL1:exon2;af=0.05" \
  -o output/
```

AF specifiers:
- `af=0.15` — exact allele fraction (any value in (0, 1])
- `af=het` — sample from Beta(40,40), centered at ~0.5 (simulates germline heterozygous)
- `af=hom` — fixed at 1.0 (homozygous)
- *(omitted)* — uses the global `--allele-fraction` (default: 0.5)

### VCF input

Load variant specifications from a VCF file instead of (or in addition to) `--event` flags:

```bash
# From VCF only
spike --bam sample.bam --reference GRCh38.fasta \
  --vcf variants.vcf \
  -o output/

# Combine VCF and manual events
spike --bam sample.bam --reference GRCh38.fasta \
  --vcf structural_variants.vcf \
  --event "snp:chr17:7577120:C:T;af=0.3" \
  -o output/
```

Supported VCF records:
- **DEL, DUP, INV, INS** — standard SVTYPE records with END or SVLEN
- **BND** — breakend notation, paired by MATEID into Fusion events (detects inverted orientation from bracket pattern)
- **SNP/indel** — standard REF/ALT records without SVTYPE
- **AF from INFO** — reads `SIM_VAF`, `VAF`, or `AF` fields (checked in that order)

Both plain `.vcf` and bgzip-compressed `.vcf.gz` files are supported.

### CRAM input

spike supports CRAM files transparently — just pass a `.cram` file instead of `.bam`:

```bash
spike --bam sample.cram --reference GRCh38.fasta \
  --event "del:chr17:43045000-43055000" \
  -o output/
```

The `--reference` FASTA is required for CRAM decoding (it is also required for haplotype construction, so there is no extra burden). Both `.crai` and `.cram.crai` index conventions are supported.

### LOH simulation with a gVCF

For more accurate haplotype-aware read suppression, provide a pre-called VCF (e.g., from DeepVariant) with het SNP genotypes:

```bash
spike --bam sample.bam --reference GRCh38.fasta \
  --gvcf deepvariant_calls.g.vcf.gz \
  --event "del:chr19:11090000-11133000" \
  -o output/
```

Without `--gvcf`, spike uses an automatic pileup approach to discover het SNPs. The gVCF approach is more accurate when calls are available.

### Configurable aligner

The generated `align.sh` script uses `bwa-mem2` by default. Override with `--aligner`:

```bash
# Use minimap2
spike --bam sample.bam --reference GRCh38.fasta \
  --event "del:chr17:43045000-43055000" \
  --aligner minimap2 \
  -o output/

# Use bowtie2
spike --bam sample.bam --reference GRCh38.fasta \
  --event "del:chr17:43045000-43055000" \
  --aligner bowtie2 \
  -o output/

# Custom aligner (must accept <ref> <r1.fq.gz> <r2.fq.gz>, produce SAM on stdout)
spike --bam sample.bam --reference GRCh38.fasta \
  --event "del:chr17:43045000-43055000" \
  --aligner "my_aligner -t 8" \
  -o output/

# Use a different samtools binary
spike --bam sample.bam --reference GRCh38.fasta \
  --event "del:chr17:43045000-43055000" \
  --samtools /opt/samtools/bin/samtools \
  -o output/
```

### Auto-align after simulation

Run alignment automatically after FASTQ generation:

```bash
spike --bam sample.bam --reference GRCh38.fasta \
  --event "del:chr20:40000000-40005000" \
  --align \
  -o output/
# Produces output/sim.bam (sorted and indexed)
```

Or run the generated script manually:

```bash
bash output/align.sh                    # Uses defaults from spike run
bash output/align.sh /path/to/ref 8    # Override reference and thread count
```

### Controlling the read extraction region

By default, spike extracts reads from a region around each event (event +/- `--flank`). Override this for full-gene coverage:

```bash
spike --bam sample.bam --reference GRCh38.fasta \
  --event "del:LDLR:exon4-exon8" \
  --exon-bed gene_exons.bed \
  --region "chr19:11080000-11140000" \
  -o output/
```

### Indel error model

By default, synthetic read errors are substitution-only. To include realistic indel errors:

```bash
spike --bam sample.bam --reference GRCh38.fasta \
  --event "del:chr17:43045000-43055000" \
  --indel-error-rate 0.05 \
  -o output/
```

The `--indel-error-rate` specifies the fraction of sequencing errors that are indels (vs substitutions). Typical Illumina values are 0.0-0.05.

## Output files

| File | Description |
|------|-------------|
| `R1.fq.gz` | Forward reads (gzipped FASTQ) |
| `R2.fq.gz` | Reverse reads (gzipped FASTQ) |
| `truth.vcf` | VCF with simulated variant records and AF annotations |
| `align.sh` | Shell script for alignment + samtools sort/index |
| `sim.bam` | Aligned BAM (only if `--align` used) |

### Truth VCF

The truth VCF contains one record per simulated event with:
- Standard VCF fields (CHROM, POS, REF, ALT)
- `SVTYPE` and `END` / `SVLEN` for structural variants
- `SIM_VAF` in the INFO field with the actual allele fraction used
- `SIM_GENE` with the associated gene name
- BND records for fusions (with `]`/`[` notation reflecting orientation)

## CLI reference

```
spike --help

Options:
  -b, --bam <BAM>                  Input BAM/CRAM file (coordinate-sorted, indexed)
  -r, --reference <FASTA>          Reference FASTA (with .fai index)
  -e, --event <SPEC>               Event specification(s), can be repeated
      --vcf <VCF>                  Input VCF file with variant records
      --exon-bed <BED>             Exon BED file (required for gene-based events)
      --allele-fraction <AF>       Global target allele fraction [default: 0.5]
  -o, --output <DIR>               Output directory [default: /tmp/spike]
      --seed <SEED>                Random seed for reproducibility [default: 42]
  -t, --threads <N>                Threads for BAM reading [default: 4]
      --region <chr:start-end>     Read extraction region (overrides event +/- flank)
      --flank <BP>                 Flanking region around events [default: 10000]
      --min-mapq <MAPQ>            Minimum mapping quality [default: 20]
      --aligner <CMD>              Aligner for align script [default: bwa-mem2]
                                   Presets: bwa-mem2, minimap2, bowtie2, or custom
      --samtools <PATH>            Path to samtools binary [default: samtools]
      --align                      Run alignment after FASTQ generation
      --indel-error-rate <RATE>    Indel error fraction [default: 0.0]
      --gvcf <VCF>                 gVCF/VCF with het SNP calls for LOH simulation
```

## Architecture

```
main.rs          CLI, event parsing, orchestration
types.rs         Core types: SimEvent, ReadPair, ReadPool, SimConfig
haplotype.rs     Variant haplotype construction (segment-based)
simulate.rs      Read suppression + synthetic read tiling
synth.rs         Quality-profiled synthetic read generation
extract.rs       BAM/CRAM read pair extraction
stats.rs         Fragment length distribution
loh.rs           Loss of heterozygosity / allelic imbalance
exon.rs          Exon BED and event spec parsing
vcf_input.rs     VCF input parser (DEL/INS/DUP/INV/BND/SNP/indel)
truth.rs         Truth VCF output
fastq.rs         Gzipped paired FASTQ writer
reference.rs     Indexed FASTA reading + in-memory sequence store
bam_stats.rs     BAM/CRAM insert size and read length statistics
```

## Simulation model details

### Read suppression

Within the SV boundaries, original reads are suppressed at the target VAF rate. Reads fully inside the SV region are suppressed at `P = VAF`. Reads spanning SV boundaries are suppressed at `P = VAF * overlap_fraction`, preventing boundary depth artifacts. Reads fully outside the SV region are always kept.

For additive events (DUP, Fusion), no reads are suppressed — only depth copies or chimeric reads are added.

### Synthetic read tiling

Synthetic reads are tiled uniformly across the variant haplotype at a rate proportional to `coverage * VAF`. Fragment lengths are sampled from the empirical distribution of the donor reads. For insertions, rejection sampling ensures fragments are not placed entirely within novel sequence.

For additive events (DUP, Fusion), reads are placed only near breakpoints to avoid inflating flank coverage.

### Depth copies (DUP)

For duplications, reads overlapping the duplicated region by >= 50% of fragment length are copied to simulate the increased coverage of the extra copy. When haplotype information is available, copies are drawn preferentially from the duplicated haplotype for correct allelic imbalance.

### LOH simulation

For heterozygous deletions, het SNP positions within the deleted region are identified (via pileup or gVCF), and reads from the deleted haplotype are preferentially suppressed. Synthetic reads carry only the surviving haplotype's allele at het SNP positions, producing realistic loss-of-heterozygosity signal.

For heterozygous duplications, the extra copy's reads carry one haplotype's alleles, shifting het SNPs from ~50/50 to ~33/67 allele balance.

### Quality profile

The quality model is learned per-cycle and per-base from the donor reads. Two levels of conditioning:
1. **Base-conditioned**: `(read_number, cycle, sequenced_base)` — captures base-specific effects like the Illumina GG quality dip
2. **Cycle-only fallback**: `(read_number, cycle)` — used when a base-conditioned bin has too few observations

Error rates are derived from the sampled quality scores: `P(error) = 10^(-Q/10)`.

### Indel error model

When `--indel-error-rate` is set above 0, a fraction of sequencing errors are modeled as insertions or deletions (50/50 split) rather than substitutions. This maintains fixed read length: insertions consume an output position without advancing the reference, and deletions skip a reference base without consuming an output position.

## Pipeline validation examples

### Validate a deletion caller

```bash
# 1. Simulate a 10kb het deletion in BRCA1
spike --bam NA12878.bam --reference GRCh38.fasta \
  --event "del:chr17:43045000-43055000;af=het" \
  --align \
  -o sim_brca1_del/

# 2. Run your variant caller on sim_brca1_del/sim.bam
my_sv_caller --input sim_brca1_del/sim.bam --output calls.vcf

# 3. Compare calls.vcf against sim_brca1_del/truth.vcf
```

### Validate a fusion caller

```bash
# Simulate BCR-ABL1 at 5% VAF (somatic-like)
spike --bam tumor.bam --reference GRCh38.fasta \
  --exon-bed gene_exons.bed \
  --event "fusion:BCR:exon14:ABL1:exon2;af=0.05" \
  --align \
  -o sim_bcr_abl/
```

### Sensitivity titration

```bash
# Test caller sensitivity across a range of VAFs
for vaf in 0.01 0.05 0.10 0.25 0.50; do
  spike --bam sample.bam --reference GRCh38.fasta \
    --event "del:chr17:43045000-43055000;af=${vaf}" \
    --align \
    -o "sim_vaf_${vaf}/"
done
```

### Multi-event truth set

```bash
# Simulate a panel of variants from a VCF truth set
spike --bam sample.bam --reference GRCh38.fasta \
  --vcf panel_variants.vcf \
  --align \
  -o sim_panel/
```
