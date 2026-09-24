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
    |  2. Learn quality           |  Markov chain Q score model from real data
    |  3. Build haplotype         |  Linear variant sequence from ordered segments
    |  4. Read sample's SNPs      |  Het + hom-alt SNPs, phased; pick the event copy
    |  5. Suppress reads          |  Remove reads by copy within the haplotype footprint
    |  6. Tile synthetic          |  New reads carrying their copy's alleles, learned Q
    +-----------------------------+
                  |
                  v
     Paired FASTQ + Truth VCF + Alignment script
```

### Core concepts

**Variant haplotype model**: Every variant — from a single SNP to a multi-kilobase structural rearrangement — is represented as an ordered list of *segments*, each drawn from a reference region (possibly reverse-complemented) or from novel sequence. These segments are concatenated into a single linear haplotype sequence. Reads tiled uniformly across this linear sequence become automatically chimeric when they span a segment boundary. This single mechanism handles all SV types without any per-type breakpoint logic.

**Read suppression and replacement**: For non-additive events (DEL, INV, INS, SNP, full-model DUP), original reads within the haplotype's reference footprint are suppressed at the target VAF rate, and new synthetic reads tiled across the variant haplotype replace the removed fraction. For additive events (Fusion, junction-model DUP), all original reads are kept and synthetic reads are added on top.

**Quality-aware synthesis**: Instead of cloning real reads (which produces exact duplicates flagged by dedup tools), spike learns a first-order Markov chain quality model from the donor reads — capturing both per-cycle quality degradation and the inter-position correlation of quality scores — and generates independent synthetic reads with realistic quality profiles and correlated sequencing errors.

**The sample's two copies**: The event goes on one of the sample's two copies of the region. spike reads the sample's own SNPs around every event (het and hom-alt, via pileup or a pre-called gVCF), phases the het SNPs, and picks the event copy's haplotype. Original reads are removed by copy, and synthetic reads carry the alleles of the copy they come from. So SNPs in the flanks keep their allele balance and hom-alt SNPs stay hom-alt; a het deletion turns het SNPs inside it homozygous (LOH), and a het duplication shifts them to ~33/67. See [The sample's SNPs](#the-samples-snps).

## Supported variant types

| Type | Event spec | VCF SVTYPE | Description |
|------|-----------|------------|-------------|
| Deletion | `del:chr:start-end` or `del:GENE:exon4-exon8` | DEL | Region removed from one haplotype |
| Duplication | `dup:chr:start-end` or `dup:GENE:exon4-exon8` | DUP | Tandem duplication in place |
| Inversion | `inv:chr:start-end` or `inv:GENE:exon4-exon8` | INV | Region reversed in place |
| Insertion | `ins:chr:pos:length` or `ins:chr:pos:ACGT` | INS | Novel sequence inserted at position |
| Fusion | `fusion:GENEA:exonN:GENEB:exonM` | BND | Two breakpoints joined across genes/chromosomes; orientation follows the gene strands |
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

# Output: output/R1.fq.gz, output/R2.fq.gz, output/truth.vcf, output/align.sh,
#         output/events.bed, output/replaced_reads.txt, output/merge.sh,
#         output/README.md
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

Exon numbers are read from the names (`TP53_exon1` → exon 1), so they follow transcript order on both strands: on a minus-strand gene, exon 1 has the highest coordinates. The strand is inferred from this numbering. If no exon name of a gene carries a number, exons are numbered by genomic position (with a warning), which is backwards for minus-strand genes. Keep one transcript per gene: a gene whose exon numbers repeat, or where only some names carry a number, is rejected.

### Gene fusions

Fusions join two breakpoints, potentially across different chromosomes. Gene A keeps exon N and everything upstream; gene B keeps exon M and everything downstream. The join orientation follows the gene strands from the exon BED, so genes on opposite strands (e.g. EML4-ALK) get a reverse-complement join automatically. The old `:inv` suffix is no longer accepted:

```bash
# BCR-ABL1 fusion (forward)
spike --bam sample.bam --reference GRCh38.fasta \
  --exon-bed gene_exons.bed \
  --event "fusion:BCR:exon14:ABL1:exon2" \
  -o output/

# EML4-ALK (genes on opposite strands)
spike --bam sample.bam --reference GRCh38.fasta \
  --exon-bed gene_exons.bed \
  --event "fusion:EML4:exon13:ALK:exon20" \
  -o output/

# Fusion at low somatic VAF
spike --bam sample.bam --reference GRCh38.fasta \
  --exon-bed gene_exons.bed \
  --event "fusion:BCR:exon14:ABL1:exon2;af=0.05" \
  -o output/
```

### SNPs and small indels

Small variants use the `snp:` prefix (applies to SNPs, MNVs, small deletions, and small insertions). Positions are 1-based in the CLI/API syntax (converted to 0-based internally):

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

By default, overlapping events on the same chromosome are rejected to keep effects independent. Use `--allow-overlap` to override (overlaps are then simulated independently and merged, which is approximate in overlap zones).

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
- **BND** — breakend notation, paired by MATEID into Fusion events. All four forms are read (`t[p[`, `]p]t`, `t]p]`, `[p[t`); either record of a mate pair gives the same fusion
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

### The sample's SNPs from a gVCF

For more accurate haplotype-aware simulation, provide a pre-called VCF (e.g., from DeepVariant) with SNP genotypes:

```bash
spike --bam sample.bam --reference GRCh38.fasta \
  --gvcf deepvariant_calls.g.vcf.gz \
  --event "del:chr19:11090000-11133000" \
  -o output/
```

Without `--gvcf`, spike uses an automatic pileup approach to discover het and hom-alt SNPs. The gVCF approach is more accurate when calls are available, and a phased VCF also links SNPs that no read spans.

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

### Merging into the original BAM

After aligning, `merge.sh` substitutes the spiked reads back into the original BAM:

```bash
bash output/align.sh                              # Step 1: produce sim.bam
bash output/merge.sh                              # Step 2: produce merged.bam (full genome)
bash output/merge.sh /moved/original.bam ref.fa 8 # Override paths/threads (still the SAME BAM spike ran on)
```

The originals to replace are named, not located. spike writes every read name it extracted to `replaced_reads.txt`, and `merge.sh` drops exactly those records (`samtools view -N`, so samtools >= 1.13 is required) before merging `sim.bam` in. Removing by name rather than by event region matters in both directions:

- Records inside the event regions that spike never extracted — PCR duplicates, non-proper pairs, reads below `--min-mapq`, reads whose mate was filtered — are **kept**. Removing by region deleted them without putting anything back, which cost real depth.
- Mates that lie outside the event regions but whose pair spike did extract are **removed**, because `sim.bam` carries their replacement. Removing by region left them in, so they appeared twice.
- Pairs spike extracted but could not use — a mate with no stored base qualities, see [Missing or unusable donor base qualities](#missing-or-unusable-donor-base-qualities) — are **removed without replacement**. That costs depth; left in, they would dilute the realised VAF of every event they overlap.

`ORIGINAL_BAM` (the first argument to `merge.sh`) must be the exact BAM spike was run on: `replaced_reads.txt` names reads spike found in that BAM, and `-N` only removes names it can find. A different BAM shares essentially no read names, so nothing would be removed and `sim.bam` would be merged onto full, un-thinned original depth — `merge.sh` guards against this by comparing how many of the expected names it actually matched against how many `replaced_reads.txt` lists, and aborts with an error instead of silently producing a wrong `merged.bam`.

Records spike extracted and then suppressed (the deleted copy of a heterozygous deletion, for instance) stay gone: that absence *is* the simulated variant.

The records kept because spike never extracted them (PCR duplicates, non-proper pairs, low-MAPQ or orphaned-mate reads) are real original reads that now sit inside an event's footprint, so an event's residual depth/allele fraction in `merged.bam` is no longer exactly the simulated value. Measured on an HG002 chr20 run: `validate.rs` only skips secondary/supplementary/duplicate/QC-fail and low-MAPQ reads — it has no proper-pair or mate-unmapped filter — so of the 1,436 records recovered by this change, the 87 non-proper-pair and 39 orphaned-mate records (126 total, 1.2% of the 10,501 in-BED records) reached an AF or depth measurement in `spike validate`; a consumer that counts duplicates rather than skipping them could see the residual shift by up to the full recovered fraction (13.7%).

`align.sh` tags the simulated reads `@RG ID:sim SM:<sample>`, where `<sample>` is the `SM` of the original BAM's first `@RG` line, so `merged.bam` stays single-sample. If the original BAM's read groups carry different `SM` values it is already multi-sample; the first one still wins and spike logs a warning. A BAM with no `@RG SM` at all falls back to `SM:SIM`. Characters outside `[A-Za-z0-9._+@:-]` are replaced with `_` so the name is safe inside the generated `@RG` line.

`merged.bam` is appropriate for end-to-end testing where the caller needs to see the full genome (e.g., tools that estimate background noise from off-target regions). `sim.bam` is sufficient for targeted callers or focused benchmarking.

### Validating the spike-in

Use `spike validate` to automatically verify that the spike-in reads in a simulated BAM look realistic — correct depth, allele fraction, and read-level signals at each event:

```bash
spike validate \
  --bam output/sim.bam \
  --truth output/truth.vcf \
  --reference GRCh38.fasta

# JSON output for automated pipelines
spike validate \
  --bam output/sim.bam \
  --truth output/truth.vcf \
  --reference GRCh38.fasta \
  --json
```

Options:

```
spike validate --bam <BAM> --truth <VCF> --reference <FASTA> [OPTIONS]

  --bam, -b        Simulated BAM/CRAM file (required)
  --truth, -t      Truth VCF from spike (required)
  --reference, -r  Reference FASTA with .fai index (required)
  --min-mapq       Minimum MAPQ for counting reads (default: 20)
  --flank          Flanking bp for coverage comparison (default: 5000)
  --json           Output JSON instead of text table
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

`--region` does not replace the event window. For each event on the same
chromosome, spike extracts the region and the event +/- `--flank`, merged into
one query when they overlap or touch and kept as two queries when they do not.
The gap between a region and a distant event on the same chromosome is never
read, so a fusion partner megabases away costs one extra event-sized window
rather than every read in between. A region on another chromosome than the
event is ignored for that event. Reads shared by two windows enter the donor
pool once.

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
| `events.bed` | Extraction regions (event ± flank) used to build the spike-in |
| `replaced_reads.txt` | Names of the originals spike extracted, including pairs dropped for unusable quality; `merge.sh` removes exactly these |
| `align.sh` | Aligns R1/R2 → `sim.bam` (event regions only) |
| `merge.sh` | Merges `sim.bam` into the original BAM → `merged.bam` (full genome) |
| `README.md` | Run log: command, events table, read counts, next-step instructions |
| `sim.bam` | Aligned BAM covering event regions (produced by `align.sh`) |
| `merged.bam` | Original BAM with spiked reads substituted (produced by `merge.sh`) |

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

Haplotype-based read spike-in simulator for genomic variants

Usage: spike [OPTIONS] --bam <BAM> --reference <REFERENCE>

Options:
  -b, --bam <BAM>
          Input BAM file (coordinate-sorted, indexed)

  -r, --reference <REFERENCE>
          Reference FASTA (with .fai index)

  -e, --event <EVENT>
          Event specification(s). Can be repeated. Formats: --event "del:chr20:30000000-30005000"           (coordinate-based) --event "del:GENE:exon4-exon8"                  (gene-based, requires --exon-bed) --event "dup:GENE:exon4-exon8"                  (gene-based duplication) --event "inv:GENE:exon4-exon8"                  (gene-based inversion) --event "fusion:GENEA:exon14:GENEB:exon2"       (fusion, requires --exon-bed) --event "dup:chr20:30000000-30005000" --event "inv:chr20:30000000-30005000" --event "ins:chr20:30000000:500"                (random insertion sequence) --event "ins:chr20:30000000:ACGTACGT"           (explicit insertion sequence) --event "snp:chr20:30000000:A:T"                (SNP/small variant, POS is 1-based) --event "snp:chr20:30000000:A>T"                (alternate syntax, POS is 1-based) --event "snp:chr20:30000000:ACG:A"              (small deletion, POS is 1-based) --event "snp:chr20:30000000:A:ACGT"             (small insertion, POS is 1-based) Per-event AF (appended with ;): --event "del:GENE:exon4-exon8;af=0.15" --event "fusion:GENEA:exon14:GENEB:exon2;af=het" af=<number>: exact AF, af=het: Beta(40,40)~0.5, af=hom: 1.0

      --vcf <VCF>
          Input VCF file with variant records. Supports DEL, INS, DUP, INV, BND, and standard SNP/indel records (no SVTYPE, explicit REF/ALT alleles). Can be combined with --event. At least one of --event or --vcf required

      --exon-bed <EXON_BED>
          Exon BED file. Required when using gene-based --event specs (e.g. "del:GENE:exon4-exon8")

      --allele-fraction <ALLELE_FRACTION>
          Target allele fraction (0.0-1.0)
          
          [default: 0.5]

  -o, --output <OUTPUT>
          Output directory
          
          [default: /tmp/spike]

      --seed <SEED>
          Random seed for reproducibility
          
          [default: 42]

  -t, --threads <THREADS>
          Number of threads for BAM reading
          
          [default: 4]

      --region <REGION>
          Extra read extraction region (e.g. "chr19:11080000-11140000"), on top of event ± flank -- it does not replace the event window. Use this to ensure the output BAM covers the full gene/region of interest. For an event on the same chromosome, the region and the event ± flank are merged into one query when they overlap or touch, and kept as two queries when they do not, so a distant fusion partner costs one extra event-sized window rather than every read in between. A region on another chromosome than the event is ignored for that event

      --flank <FLANK>
          Flanking region (bp) to include around events. Always defines the event window; when --region is also set, the region adds another window alongside it (see --region) rather than replacing this one
          
          [default: 10000]

      --min-mapq <MIN_MAPQ>
          Minimum mapping quality for donor reads
          
          [default: 20]

      --aligner <ALIGNER>
          Aligner for alignment script. Presets: "bwa-mem2" (default), "minimap2", "bowtie2", or a custom command that accepts <ref> <r1.fq.gz> <r2.fq.gz> and produces SAM on stdout
          
          [default: bwa-mem2]

      --samtools <SAMTOOLS>
          Path to samtools binary (used for sort/index in align script)
          
          [default: samtools]

      --align
          Automatically run alignment after FASTQ generation

      --indel-error-rate <INDEL_ERROR_RATE>
          Indel error rate per base in synthetic reads (fraction of total error that is indel rather than substitution). Default 0.0 means substitution-only. Typical Illumina: 0.0 to 0.05
          
          [default: 0]

      --gvcf <GVCF>
          Optional gVCF/VCF with SNP calls for LOH simulation. When provided, het SNP positions are extracted from this file to determine which reads belong to the deleted haplotype (more accurate than the default pileup-based approach). Supports .vcf and .vcf.gz (requires bcftools in PATH for .vcf.gz)

      --allow-overlap
          Allow overlapping events on the same chromosome.
          
          By default, overlapping events are rejected to keep event effects independent and truth interpretation unambiguous.

      --dup-model <DUP_MODEL>
          Duplication model: "full" (default) builds a full tandem haplotype with duplicated region appearing twice, producing both junction reads and correct depth increase from a single tiling pass. "junction" uses the legacy junction-only haplotype with separate depth copies
          
          [default: full]

  -h, --help
          Print help (see a summary with '-h')

  -V, --version
          Print version
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
loh.rs           The sample's SNPs: two copies, phasing, read assignment
exon.rs          Exon BED and event spec parsing
vcf_input.rs     VCF input parser (DEL/INS/DUP/INV/BND/SNP/indel)
truth.rs         Truth VCF output
fastq.rs         Gzipped paired FASTQ writer
reference.rs     Indexed FASTA reading + in-memory sequence store
bam_stats.rs     BAM/CRAM insert size and read length statistics
validate.rs      `spike validate` subcommand: automated spike-in quality checks
```

## Simulation model details

All SV types share a common simulation framework: (1) build a linear variant haplotype from segments, (2) suppress original reads in the affected region at the VAF rate, (3) tile synthetic reads across the haplotype to replace the suppressed fraction. The details of what the haplotype looks like, how reads are suppressed, and what observable signals result differ by variant type.

A segment whose flank runs off the end of a contig is truncated to what the reference actually holds: the haplotype is shorter there, its reference footprint stops at the contig end, and the count of synthetic reads follows the shorter haplotype.

### Deletion (DEL)

**Haplotype structure** (2 segments):
```
[left_flank]            [right_flank]
ref[start-F .. start]   ref[end .. end+F]
```
where `F` = flank size (default 10kb). The deleted region `[start, end)` is absent from the haplotype; the left and right flanking segments are placed adjacent.

**Simulation**: Original read pairs lying entirely within the haplotype footprint `[start-F, end+F)` are suppressed at rate `P = VAF`; pairs that stick out of it are kept. Synthetic reads are tiled uniformly across the two-segment haplotype. Reads that span the junction between left_flank and right_flank are chimeric — when re-aligned to the reference, they produce split reads and discordant pairs that span the deletion breakpoint.

**Observable signals in the output BAM**:
- Reduced depth in `[start, end)` proportional to VAF (e.g., ~0.5x for het)
- Split reads / soft-clipped reads at both breakpoints
- Discordant read pairs spanning the deleted region (larger insert size than expected)
- LOH at het SNP positions within the deletion (see LOH section)

**Haplotype-aware suppression (LOH)**: Reads of the event copy are suppressed first (see [Read suppression details](#read-suppression-details)), so at VAF 0.5 every read of the deleted copy is gone inside the deletion: het SNPs there become homozygous in the output. In the flanks, the synthetic reads that replace the event copy's reads carry its alleles, so flank SNPs keep their balance.

### Tandem duplication (DUP)

spike supports two models for simulating tandem duplications, selected via `--dup-model`:

#### Full tandem model (default, `--dup-model full`)

**Haplotype structure** (4 segments):
```
[left_flank]   [dup_copy_1]       [dup_copy_2]       [right_flank]
ref[S-F .. S]  ref[S .. E]        ref[S .. E]        ref[E .. E+F]
```
where `S` = dup_start, `E` = dup_end, `F` = flank size. The duplicated region appears twice in tandem.

**Simulation**: This model uses the same suppress-and-replace approach as deletions. Original reads within `[S-F, E+F)` are suppressed at the VAF rate, and synthetic reads are tiled *uniformly* across the full 4-segment haplotype. Because the haplotype contains two copies of the duplicated region, tiling naturally produces:
- ~2x as many reads mapping to `[S, E)` per unit length compared to the flanks
- Chimeric reads at the copy1-to-copy2 junction (the internal breakpoint at position E in haplotype space, which maps to the E→S join)

This model is preferred because a single tiling pass produces both the correct depth profile and the junction evidence, without needing separate depth-copy logic.

**Observable signals**:
- Increased depth in `[S, E)` proportional to `1 + VAF` (e.g., ~1.5x for het)
- Split reads / soft-clipped reads at the tandem junction
- Discordant read pairs with reduced insert size spanning the junction
- Allelic imbalance at het SNPs within the DUP (see below)

#### Junction model (legacy, `--dup-model junction`)

**Haplotype structure** (2 segments):
```
[left_of_junction]    [right_of_junction]
ref[E-F .. E]         ref[S .. S+F]
```
This covers only the junction region where the end of the duplicated region meets the start of the duplicate copy (the E→S chimeric breakpoint).

**Simulation**: This model is additive — all original reads are kept. Two separate sets of synthetic reads are generated:

1. **Junction reads**: Tiled across the 2-segment haplotype near the breakpoint. These produce split reads and discordant pairs at the E→S junction.
2. **Depth copies**: For each real read pair overlapping the DUP region by >= 50% of its fragment length, a synthetic copy is generated with new quality scores and error profile, at a rate set by the pair's copy:
   - Event copy reads: copied at rate `min(1, 2*VAF)`
   - Other copy reads: copied at rate `max(0, 2*VAF - 1)`
   - Reads of unknown copy: copied at rate `VAF`

   Each depth copy carries the alleles of the copy it repeats.

This model is kept for backward compatibility but produces less naturally integrated results than the full model.

**Allelic imbalance for DUPs**: For heterozygous duplications, one haplotype has 2 copies while the other has 1. Het SNPs shift from ~50/50 allele balance to ~33/67 (2:1 ratio), because the synthetic reads carry the duplicated copy's alleles. In the full model the alleles are written into the haplotype sequence before tiling; in the junction model, into each depth copy.

### Inversion (INV)

**Haplotype structure** (3 segments):
```
[left_flank]          [inverted_region]            [right_flank]
ref[S-F .. S]         revcomp(ref[S .. E])         ref[E .. E+F]
```
The middle segment is the reverse complement of the original reference sequence.

**Simulation**: Same suppress-and-replace approach as deletions. Reads tiled across the haplotype that span the left_flank→inverted boundary or the inverted→right_flank boundary are chimeric and produce split reads at both inversion breakpoints. Reads landing entirely within the inverted segment align to the reference in the opposite orientation.

**Observable signals**:
- No depth change (the region is the same length, just reversed)
- Split reads at both breakpoints (left and right)
- Read pairs with unexpected orientation near breakpoints (FR→RF or FF/RR)
- Reads in the inverted region may show soft-clipping at the boundaries

### Insertion (INS)

**Haplotype structure** (3 segments):
```
[left_flank]          [novel_sequence]       [right_flank]
ref[P-F .. P]         <inserted bases>       ref[P .. P+F]
```
The inserted sequence is either user-specified or randomly generated. It has no reference origin (novel sequence).

**Simulation**: Suppress-and-replace. Reads spanning left_flank→novel or novel→right_flank produce chimeric reads at the insertion point. Rejection sampling during tiling ensures fragments are not placed entirely within the novel sequence (such reads wouldn't align to the reference at all).

**Observable signals**:
- No depth change in flanking regions
- Soft-clipped reads at the insertion point (with clipped bases matching the inserted sequence)
- Discordant insert sizes for pairs where one read is in the insertion and the other is in flanking reference

### Gene fusion (BND)

**Haplotype structure** (2 segments; `bpA`, `bpB` are cuts between two bases):
```
Forward     (t[p[):  ref_A[bpA-F .. bpA]           |  ref_B[bpB .. bpB+F]
Left-left   (t]p]):  ref_A[bpA-F .. bpA]           |  revcomp(ref_B[bpB-F .. bpB])
Right-right ([p[t):  revcomp(ref_A[bpA .. bpA+F])  |  ref_B[bpB .. bpB+F]
```
Two breakpoints on potentially different chromosomes are joined. Left-left and right-right are the two junctions of an inversion-type rearrangement; each keeps the side of both cuts named, with one piece reverse-complemented. The truth VCF writes the matching BND form for each, with POS on the base next to the junction.

**Simulation**: Fusions are additive — all original reads are kept. Synthetic reads are tiled near the breakpoint (breakpoint-only mode) so that they cross the junction. This avoids inflating coverage in the flanking regions where real reads already provide normal depth.

**Observable signals**:
- Split reads at the fusion junction (one side mapping to gene A, the other to gene B)
- Discordant read pairs with mates on different chromosomes (or unexpectedly far apart on the same chromosome)
- No depth change in flanking regions (additive only near the breakpoint)

### SNP / small indel (SmallVariant)

**Haplotype structure** (3 segments):
```
[left_flank]          [alt_allele]          [right_flank]
ref[P-F .. P]         <ALT bases>           ref[P+len(REF) .. P+len(REF)+F]
```
Works for all small variant types:
- **SNP** (A→T): alt segment = `[T]`, skips 1 ref base
- **MNV** (AC→TG): alt segment = `[TG]`, skips 2 ref bases
- **Small deletion** (ACG→A): alt segment = `[A]`, skips 3 ref bases
- **Small insertion** (A→ACGT): alt segment = `[ACGT]`, skips 1 ref base

For equal-length substitutions (SNPs/MNVs), the alt segment retains a reference origin mapping for correct coordinate translation. For indels, the alt segment has no reference origin.

**Simulation**: Suppress-and-replace. Reads crossing the variant position carry the ALT allele; elsewhere they carry the sample's own alleles. The sample's alleles never overwrite the variant's own bases.

**Observable signals**:
- ALT allele at the expected frequency in pileup
- For small indels, soft-clipped reads near the variant position

## Read suppression details

The suppress-and-replace model classifies each original read pair relative to the haplotype's reference footprint:

| Relation | Behavior |
|---|---|
| **Outside** (both reads entirely outside footprint) | Always kept |
| **Inside** (both reads entirely within footprint) | Suppressed by copy (below) |
| **Overlapping** (fragment straddles footprint boundary) | Always kept |

Only pairs entirely inside the footprint are replaced, because synthetic fragments never extend past the haplotype ends either. Near an edge, both the suppressed and the synthetic depth taper off the same way, so total depth stays flat.

Every event reports the names it suppressed, and a read suppressed by any event stays out of the output, even if a nearby event's pool also contains it.

Inside pairs are suppressed by the copy they come from (see [Read classification](#read-classification)), for every event type and VAF, in the flanks as well as the event:
- **Event copy reads**: suppressed at `P = min(1, 2*VAF)`
- **Other copy reads**: suppressed at `P = max(0, 2*VAF - 1)`
- **Reads of unknown copy** (no het SNP, or a tie): suppressed at `P = VAF`

Averaged over both copies this is `VAF`. At VAF 0.5 it removes every read of the event copy and none of the other — correct LOH. Above 0.5 the event is on both copies in some cells.

## Synthetic read tiling

The number of synthetic reads to tile is:

- **Non-additive events**: `n = round(coverage * VAF * starts / mean_fragment_length)`, where `starts` is the number of fragment start positions tiling can use: `haplotype_length - mean_fragment_length`, minus starts that would lie wholly inside inserted sequence. Starts are uniform, so the flanks get `VAF * coverage` synthetic depth, replacing what was suppressed.
- **Additive events** (breakpoint-only tiling): every original read is kept, so `n = round(coverage * VAF / (1 - VAF))` per breakpoint makes junction fragments a `VAF` fraction of the depth there (VAF capped at 0.95).

`mean_fragment_length` is the library's own mean. `coverage` is the donor pool's mean fragment depth in a 2 kb window around the first breakpoint, counted only on that breakpoint's own chromosome: a fusion's pool holds both partners, and reads from the far side would otherwise be added to the near side's depth.

Fragment lengths are sampled from the empirical distribution of the donor reads. Each fragment is placed at a random position on the haplotype and a read pair is synthesized with quality scores from the learned Markov model. The fragment's left end is always read forward and its right end reverse (FR), and a coin flip out of the same seeded stream decides which of the two is R1: about half the pairs come out F1R2 and half F2R1, as in a real library, so read-orientation filters (Mutect2's, for one) see a balanced strand mix. Whichever mate is R1 is sampled from the R1 quality model.

Each fragment comes from one of the sample's copies and carries its alleles. Up to VAF 0.5 all fragments come from the event copy. Above 0.5 the other copy gives `max(0, 2*VAF - 1) / (2*VAF)` of them, matching what was suppressed from it.

For additive events, fragment placement is restricted to positions that cross a segment boundary (breakpoint). For non-additive events, placement is uniform across the haplotype, with rejection sampling to avoid placing fragments entirely within novel (non-reference) sequence.

## The sample's SNPs

For every event, spike reads the sample's SNPs over the haplotype's whole reference footprint (event ± 2 kb; each side of a fusion separately), at every VAF. Het SNPs are phased and assign reads to copies; het and hom-alt alleles go into the synthetic reads. Two strategies are supported:

### Pileup (default)

A pileup counts the bases (A/C/G/T) at each position (MAPQ >= threshold, excluding secondary/supplementary/duplicate/QC-fail). Positions with >= 10 reads are called:
- **het**: two alleles each at 20%-80%,
- **hom-alt**: one allele at >= 90% that differs from the reference.

A second pass records each fragment's bases at the het SNPs only, so memory grows with the number of het SNPs, not with the region's size.

### gVCF (optional, `--gvcf`)

SNPs are loaded from a pre-called VCF (e.g., DeepVariant gVCF): het (`0/1`, `1/0`) and hom-alt (`1/1`) biallelic SNPs. A het SNP whose base a deletion on the other haplotype removes (a carried deletion allele, or ALT `*`) has no copy carrying REF, so it counts as hom-alt. A BAM pass then records which allele each read carries at the het SNPs. Phased genotypes (`0|1`, `1|0`, with an optional `PS` phase set) are used for phasing (below); a phased VCF such as a GIAB/T2T benchmark gives the most realistic result. If the gVCF has no het SNPs in the region, spike falls back to pileup.

### Phasing

A deletion or duplication affects one whole haplotype, so every het SNP in it must lose (or gain) the allele of the same haplotype. spike phases the het SNPs into blocks:

- SNPs in the same phase set of a phased gVCF are linked outright.
- Otherwise, SNPs are linked by fragments (mates pooled) that cover two or more of them; links are joined strongest first, and a link that contradicts stronger ones is ignored.

One coin flip per block picks the event copy's haplotype. SNPs that no read or phase set links (usually more than a fragment length apart) form separate blocks and get their own coin flip; with short reads alone their relative phase is unknown.

### Read classification

Each fragment is scored by how many het SNPs show the event copy's vs. the other copy's allele. More event-copy matches → event copy; more other matches → other copy. Ties, and fragments covering no het SNP, stay unclassified and fall back to random handling at the VAF rate.

## Quality profile

The quality model uses a first-order Markov chain learned from the donor reads: each position's quality score depends on the previous position's quality, capturing the autocorrelation seen in real Illumina data (runs of low quality tend to cluster together).

Quality scores are sampled using a 4-level fallback hierarchy, from most specific to least:

1. **Markov + base**: `P(Q_i | cycle, base, prev_q_bin)` — full model with base-specific effects (e.g., Illumina GG quality dip) and inter-position correlation
2. **Markov + cycle**: `P(Q_i | cycle, prev_q_bin)` — drops base conditioning when base-specific bins are sparse
3. **Base-only**: `P(Q_i | cycle, base)` — no Markov (used at cycle 0, or when Markov bins have too few observations)
4. **Cycle-only**: `P(Q_i | cycle)` — final fallback

The previous quality is quantized into 4 bins (Q0-9, Q10-19, Q20-29, Q30+) to keep transition tables tractable. Each level requires at least 30 observations before it is used; otherwise sampling falls through to the next level.

Error rates are derived from the sampled quality scores: `P(error) = 10^(-Q/10)`. When an error occurs, a random incorrect base is substituted.

### Missing or unusable donor base qualities

SAM's QUAL field is all-or-nothing per record: a read either has a full quality string or none at all (`*`). The two containers encode `*` differently and spike sees both — BAM keeps a per-base array with every byte `0xFF`, while noodles normalises a CRAM record's all-`0xFF` buffer to an *empty* one before spike ever sees it. A record with either shape, or carrying a raw quality above Q93 (the SAM maximum — a malformed record), is dropped during extraction: it is not included in the learned quality profile and does not contribute a read pair to the output. Dropped records are counted and logged as a warning (spike logs to stderr), e.g.:

```
WARN spike::extract] 457 record(s) considered for the chr20:38402500-38432500 donor pool had no quality scores (SAM '*'); their pairs were dropped from the pool and merge.sh removes them from the merged BAM
```

**A dropped pair is removed from the merged BAM even though nothing replaces it.** spike has no quality string to write for it, so it cannot come back through the FASTQ — but it is still listed in `replaced_reads.txt`, so `merge.sh` drops it. That costs real depth. Leaving it in costs more: the pair would sit inside every event it overlaps as reference support that no allele fraction can suppress and no synthetic read can replace, so the realised VAF would come out diluted by the unusable-quality fraction while the truth VCF still claimed the full one. Because the same pairs are dropped across the whole extraction window, not just inside the event, the loss cancels in a flank-normalised ratio; the dilution would not. Measured on an HG002 chr20 slice with 5% of donor pairs' quality stripped to `*`, over a 10 kb deletion: 147 of 2868 eligible records inside the event (5.13%) carried `*` quality and were un-suppressible; removing them takes that to 0.00%, at a cost of those same 147 records of depth.

Without the drop, a missing quality used to decode to an invalid FASTQ quality byte (space) for kept reads and poisoned the learned quality model, so synthetic reads sampled from it came out mostly-Q0 with effectively random bases. On the CRAM path it produced a FASTQ record with a full-length SEQ line next to a zero-length QUAL line, which `samtools import` rejects outright (`truncated file`).

As a last line of defense, `write_paired_fastq` validates every pair *before* it creates either output file — so a refusal leaves no half-written `.fq.gz` in `--output` — and returns an error if a quality string's length does not match its SEQ, or if any byte falls outside the printable Phred+33 range `!`-`~` (33-126).

### Indel error model

When `--indel-error-rate` is set above 0, a fraction of sequencing errors are modeled as insertions or deletions (50/50 split) rather than substitutions. This maintains fixed read length: insertions consume an output position without advancing the reference, and deletions skip a reference base without consuming an output position.

## Pipeline validation examples

### Validate a deletion caller

```bash
# 1. Simulate a 10kb het deletion in BRCA1
spike --bam NA12878.bam --reference GRCh38.fasta \
  --event "del:chr17:43045000-43055000;af=het" \
  -o sim_brca1_del/

# 2. Align the simulated reads
bash sim_brca1_del/align.sh

# 3a. Run caller on sim.bam (event regions only — good for targeted callers)
my_sv_caller --input sim_brca1_del/sim.bam --output calls.vcf

# 3b. Or merge into original BAM for full-genome callers
bash sim_brca1_del/merge.sh
my_sv_caller --input sim_brca1_del/merged.bam --output calls.vcf

# 4. Optionally verify the spike-in looks correct before running the caller
spike validate --bam sim_brca1_del/sim.bam \
  --truth sim_brca1_del/truth.vcf \
  --reference GRCh38.fasta

# 5. Compare calls against truth
# (compare calls.vcf against sim_brca1_del/truth.vcf)
```

### Validate a fusion caller

```bash
# Simulate BCR-ABL1 at 5% VAF (somatic-like)
spike --bam tumor.bam --reference GRCh38.fasta \
  --exon-bed gene_exons.bed \
  --event "fusion:BCR:exon14:ABL1:exon2;af=0.05" \
  -o sim_bcr_abl/

bash sim_bcr_abl/align.sh
spike validate --bam sim_bcr_abl/sim.bam \
  --truth sim_bcr_abl/truth.vcf \
  --reference GRCh38.fasta
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
