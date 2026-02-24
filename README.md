# spike

Haplotype-based read spike-in simulator for genomic variants.

spike takes a real BAM file and a set of variant specifications, then produces paired FASTQ files containing reads that carry the specified variants at a controlled allele fraction. Unlike simple read-cloning approaches, spike generates independent synthetic reads with realistic quality profiles that survive deduplication tools.

## How it works

```
Real BAM + Reference FASTA + Variant specs
                  │
                  ▼
    ┌─────────────────────────┐
    │  1. Extract read pairs  │  Real reads from the event region
    │  2. Learn quality       │  Per-cycle Q score distributions from real data
    │  3. Build haplotype     │  Linear variant sequence from ordered segments
    │  4. Suppress reads      │  Remove reads at VAF rate within SV boundaries
    │  5. Tile synthetic      │  New reads across haplotype with learned Q profile
    │  6. LOH simulation      │  Het SNPs → hom for DELs; BAF shift for DUPs
    └─────────────────────────┘
                  │
                  ▼
     Paired FASTQ + Truth VCF + Alignment script
```

**Variant haplotype model**: The variant allele is represented as an ordered list of segments, each from a reference region (possibly reverse-complemented) or novel sequence. Reads tiled uniformly across this linear sequence are automatically chimeric when they span segment boundaries — no per-SV-type breakpoint logic needed.

**Quality-aware synthesis**: Instead of cloning real reads (which produces exact duplicates flagged by dedup tools), spike learns per-cycle, base-conditioned quality score distributions from the donor reads and generates independent synthetic reads with correlated sequencing errors.

**LOH and allelic imbalance**: For heterozygous deletions, het SNPs within the deleted region become homozygous (only the surviving haplotype allele remains). For duplications, het SNPs shift from ~50/50 to ~33/67 allele balance. Het SNP positions are found via pileup (default) or from a pre-called gVCF.

## Supported variant types

| Type | Event spec | VCF SVTYPE | Description |
|------|-----------|------------|-------------|
| Deletion | `del:chr:start-end` | DEL | Region removed from one haplotype |
| Duplication | `dup:chr:start-end` | DUP | Tandem duplication in place |
| Inversion | `inv:chr:start-end` | INV | Region reversed in place |
| Insertion | `ins:chr:pos:length` | INS | Novel sequence inserted at position |
| Fusion | `fusion:GENEA:exonN:GENEB:exonM` | BND | Two breakpoints joined across genes/chromosomes |

All types support per-event allele fraction control.

## Installation

```bash
cargo build --release
# Binary at target/release/spike
```

Requires Rust edition 2021.

## Usage

### Basic: coordinate-based events

```bash
spike \
  --bam sample.bam \
  --reference GRCh38.fasta \
  --event "del:chr17:43045000-43055000" \
  -o /tmp/spike_output
```

### Multiple events with per-event allele fractions

```bash
spike \
  --bam sample.bam \
  --reference GRCh38.fasta \
  --event "del:chr17:43045000-43055000;af=0.3" \
  --event "dup:chr20:30000000-30005000;af=het" \
  -o /tmp/spike_output
```

AF specifiers: `af=0.15` (exact), `af=het` (Beta(40,40) ~ 0.5), `af=hom` (1.0).

### Gene-based events (requires exon BED)

```bash
spike \
  --bam sample.bam \
  --reference GRCh38.fasta \
  --exon-bed ldlr_exons_hg38.bed \
  --event "del:LDLR:exon4-exon8" \
  -o /tmp/spike_output
```

### Fusion events

```bash
spike \
  --bam sample.bam \
  --reference GRCh38.fasta \
  --exon-bed gene_exons.bed \
  --event "fusion:BCR:exon14:ABL1:exon2;af=0.15" \
  -o /tmp/spike_output
```

### VCF input

```bash
spike \
  --bam sample.bam \
  --reference GRCh38.fasta \
  --vcf variants.vcf \
  -o /tmp/spike_output
```

Supports DEL, INS, DUP, INV, and BND (paired by MATEID) from VCF 4.3 files. Plain `.vcf` and bgzip-compressed `.vcf.gz` are both supported. Per-event AF is read from the `VAF` or `SIM_VAF` INFO field if present.

### With LOH from a gVCF

```bash
spike \
  --bam sample.bam \
  --reference GRCh38.fasta \
  --gvcf deepvariant_calls.g.vcf.gz \
  --event "del:chr19:11090000-11133000" \
  -o /tmp/spike_output
```

### Auto-align after simulation

```bash
spike \
  --bam sample.bam \
  --reference GRCh38.fasta \
  --event "del:chr20:40000000-40005000" \
  --align \
  -o /tmp/spike_output
```

Uses bwa-mem2 via pixi. A convenience `align.sh` script is always generated in the output directory regardless of `--align`.

## Output files

| File | Description |
|------|-------------|
| `R1.fq.gz` | Forward reads (gzipped FASTQ) |
| `R2.fq.gz` | Reverse reads (gzipped FASTQ) |
| `truth.vcf` | VCF with simulated variant records and AF annotations |
| `align.sh` | Shell script for bwa-mem2 alignment + samtools sort/index |
| `sim.bam` | Aligned BAM (only if `--align` used) |

## CLI reference

```
spike --help

Options:
  -b, --bam <BAM>                  Input BAM file (coordinate-sorted, indexed)
  -r, --reference <FASTA>          Reference FASTA (with .fai index)
  -e, --event <SPEC>               Event specification(s), can be repeated
      --vcf <VCF>                  Input VCF file with SV records
      --exon-bed <BED>             Exon BED file (required for gene-based events)
      --allele-fraction <AF>       Target allele fraction [default: 0.5]
  -o, --output <DIR>               Output directory [default: /tmp/spike]
      --seed <SEED>                Random seed [default: 42]
  -t, --threads <N>                Threads for BAM reading [default: 4]
      --region <chr:start-end>     Read extraction region (overrides event ± flank)
      --flank <BP>                 Flanking region around events [default: 10000]
      --min-mapq <MAPQ>            Minimum mapping quality [default: 20]
      --align                      Run bwa-mem2 alignment after FASTQ generation
      --gvcf <VCF>                 gVCF/VCF with SNP calls for LOH simulation
```

## Architecture

```
main.rs          CLI, event parsing, orchestration
types.rs         Core types: SimEvent, ReadPair, ReadPool, SimConfig
haplotype.rs     Variant haplotype construction (segment-based)
simulate.rs      Read suppression + synthetic read tiling
synth.rs         Quality-profiled synthetic read generation
extract.rs       BAM read pair extraction
stats.rs         Fragment length distribution
loh.rs           Loss of heterozygosity / allelic imbalance
exon.rs          Exon BED and event spec parsing
vcf_input.rs     VCF input parser (DEL/INS/DUP/INV/BND)
truth.rs         Truth VCF output
fastq.rs         Gzipped paired FASTQ writer
reference.rs     Indexed FASTA reading + in-memory sequence store
bam_stats.rs     BAM insert size and read length statistics
```

## Simulation model details

### Read suppression

Within the SV boundaries, original reads are suppressed at the target VAF rate. Reads fully inside the SV region are suppressed at `P = VAF`. Reads spanning SV boundaries are suppressed at `P = VAF * overlap_fraction`, preventing boundary depth artifacts. Reads fully outside the SV region are always kept.

For additive events (DUP, Fusion), no reads are suppressed — only depth copies or chimeric reads are added.

### Synthetic read tiling

Synthetic reads are tiled uniformly across the variant haplotype at a rate proportional to `coverage * VAF`. Fragment lengths are sampled from the empirical distribution of the donor reads. For insertions, rejection sampling ensures fragments are not placed entirely within novel sequence.

### Depth copies (DUP)

For duplications, reads overlapping the duplicated region by >= 50% of fragment length are copied to simulate the increased coverage of the extra copy.

### LOH simulation

For heterozygous deletions, het SNP positions within the deleted region are identified (via pileup or gVCF), and reads from the deleted haplotype are preferentially suppressed. Synthetic reads carry only the surviving haplotype's allele at het SNP positions, producing realistic loss-of-heterozygosity signal.

For heterozygous duplications, the extra copy's reads carry one haplotype's alleles, shifting het SNPs from ~50/50 to ~33/67 allele balance.
