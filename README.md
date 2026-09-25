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
- **Reference FASTA** with `.fai` index (e.g., from `samtools faidx`). May be
  plain or bgzip-compressed (`.gz`/`.bgz`, detected by extension); a
  bgzipped FASTA also needs the matching `.gzi` index that `samtools faidx`
  produces alongside the `.fai` for a bgzipped input.
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

The exon BED file should be tab-separated with at least 4 columns. spike needs two things from it: a **gene symbol** and an **exon number**.

| Column | Standard BED | What spike does with it |
| --- | --- | --- |
| 1-3 | `chrom start end` | the exon's interval (0-based, half-open). Fewer than 4 columns on any line is a hard error -- spike does not fall back to BED3 |
| 4 | `name` | the exon name. The exon number is read from it (`LDLR_exon1` -> 1), and it is also where the gene symbol comes from when column 5 does not give one: the part before the first `_` (`LDLR_exon1` -> `LDLR`) |
| 5 | `score` | read as the gene symbol, **unless** it is empty, or a BED score -- `.` or an integer 0-1000 -- in which case it is ignored and the gene comes from the name |
| 6+ | `strand`, ... | ignored; the strand is inferred from the exon numbering (see below) |

So both layouts work: spike's own `chrom start end LDLR_exon1 LDLR` (as in `data/ldlr_deletions/ldlr_exons_hg38.bed`), and a standard BED6 `chrom start end LDLR_exon1 0 +`, which before carried its score into the gene symbol and named every gene `0`. There is no format flag and no guessing beyond that one rule: a gene symbol is never `.` and never a bare number 0-1000.

The rule is one-directional, and the other direction is where it can be wrong. `0-1000` is what the BED spec allows a score to be, but tools do not all honour it -- bedtools, MACS and others routinely emit scores above 1000, and UCSC tolerates them -- so a file scoring its exons `5000` has that read as the gene symbol and every exon lands under a gene named `5000`. The outcome is loud rather than silent: the event then fails with `gene 'LDLR' not found. Available: 5000`, or earlier still with `gene 5000: exon number 1 appears more than once` once several genes collapse into it.

**To disambiguate a file whose gene symbols look like scores, put the symbol before the first `_` of column 4**: a gene literally named `100` is read as a score in column 5 and lost, but `100_exon1` in column 4 with a score in column 5 resolves to `100` (measured). The precedence runs the other way, though, so this does not rescue an out-of-spec score: a non-empty column 5 that is not a score always wins, and `LDLR_exon1` in column 4 beside `5000` in column 5 still gives a gene named `5000` (also measured). For that file, fix column 5 -- put the score back inside 0-1000, set it to `.`, or cut the file down to its first 4 columns.

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
- **DEL, DUP, INV, INS** — standard SVTYPE records with END or SVLEN. A DEL/DUP/INV with neither falls back to the alleles, but only where they give an unambiguous span. Except for an inversion (below), the flanks `REF` and `ALT` share are stripped first (the common prefix, then any common suffix left over), so the event starts after the last base both alleles keep, not at `POS`: `POS=38412500 REF=GTTAAAGTTTATCAGAAAATT ALT=GTTAAAG SVTYPE=DEL` is the 14 bp deletion 38412507-38412520, anchored at 38412506. A single-base `ALT` is stripped like any other — it is the anchor base only when it really is `REF`'s first base, so `REF=ACGT ALT=T` deletes `ACG` and is anchored at `POS`-1, not `CGT` anchored at `POS`. What is left of the alleles then has to spell the event — `REF` keeping bases `ALT` drops, in which case the span is those bases, for all three types (`REF=ACGT ALT=A` and `REF=TGTT ALT=TG` both state one); or, for DUP, a single-base `REF` anchor whose `ALT` is that anchor plus the duplicated copy (`REF=G ALT=GACGT...`, the same form a sequence-resolved INS uses). An INV written out base for base is read whole instead, because an inversion's span is *stated* by the record rather than derived from where its alleles differ: `REF` is either the inverted region itself with `ALT` its reverse complement (an equal-length substitution, which carries no padding base and so starts at `POS` itself) or a padding base followed by that pair. Stripping shared flanks here would read `REF=AGTT ALT=AACT`, whose ends are their own complements, as a 2 bp inversion of the middle. Anything else — a complex pair, an `ALT` that is not the reverse complement where an INV needs one, a DUP whose `ALT` is longer and carries sequence in both alleles (its copy could equally be `REF`'s own span) or whose `ALT` does not begin with the `REF` anchor, a record with no length anywhere, or alleles that would start the event before the chromosome's first base (`POS=1` with no shared prefix) — is rejected with a warning naming its type and `chrom:pos` (logged to stderr) instead of being silently treated as a 1 bp event. A sequence-resolved INS is read the same way: `REF=AT ALT=ATGGG` inserts the 3 bases `GGG` after the `T`, and a record whose `REF` keeps bases its `ALT` drops is not an insertion at all and is skipped with the same kind of warning
- **BND** — breakend notation, paired by MATEID into Fusion events. All four forms are read (`t[p[`, `]p]t`, `t]p]`, `[p[t`); either record of a mate pair gives the same fusion
- **SNP/indel** — standard REF/ALT records without SVTYPE; a record whose REF and ALT are identical once case is normalized (e.g. `A`/`a`) is rejected rather than turned into a no-op "variant" in the truth VCF
- **AF from INFO** — reads `SIM_VAF` then `VAF`. Plain `AF` is **not** read by default: in a population VCF (gnomAD, 1000G) `AF` is the allele frequency in the population, not the fraction of this sample's reads that should carry the allele, so using it silently produces a truth set at the wrong VAF. A record that carries `AF` but has no usable `SIM_VAF` or `VAF` value falls back to `--allele-fraction`, and the run reports how many did — this fires for `AF=0.001` alone just as it does for `SIM_VAF=nan;AF=0.001`, not only when `AF` is the sole VAF key present. Pass `--vcf-info-af` to read `AF` as the VAF, for a VCF that really does state one there. A VAF key whose value is not a fraction in (0, 1] (`SIM_VAF=nan`, an unparseable number) is likewise reported and, once every key it has is exhausted, falls back to `--allele-fraction` — but a later key that does hold a usable value is used instead, so `SIM_VAF=nan;VAF=0.3` is not counted as a fallback at all

Both plain `.vcf` and bgzip-compressed `.vcf.gz` files are supported.

**Every record spike does not simulate is counted and reported on stderr**, per reason, so a truth set that is short of records says why rather than leaving it to be noticed downstream:

```
WARN skipped 3 VCF record(s): 1 short line (fewer than 8 columns); 1 multi-allelic ALT; 1 SVTYPE spike does not simulate
```

The reasons are: a short line (fewer than the 8 mandatory columns; a blank line is not a record at all and is not counted here), a `POS` that is not a positive integer, a multi-allelic `ALT` (spike simulates one allele per record, and taking the first of `A,T` would silently simulate half of it), an `SVTYPE` spike has no model for (`CNV`), a record with no `SVTYPE` whose alleles are not plain DNA to fall back on, and the three allele-shape rejections described above (no length or span, an `INS` whose `ALT` is not its `REF` plus inserted bases, and an event that would start before the chromosome's first base).

Every ingest also logs one line unconditionally, whether or not anything was skipped, so a clean run and a run whose counting silently failed are not both silence:

```
INFO variants.vcf: read 42 VCF record(s), simulated 42, skipped 0
```

`SVTYPE` subtypes are read: VCF v4.3 writes them with a colon. For DEL and INS the subtype names what was deleted or inserted without changing the event's shape, so any suffix is accepted — `DEL:ME:ALU` is simulated as a deletion, `INS:ME:L1` as an insertion, identically to the plain type. DUP is narrower: only the bare type and `DUP:TANDEM` are simulated, as the same local duplication; `DUP:DISPERSED` and `DUP:INT` place the copy elsewhere in the genome, not adjacent to the original, so simulating them as tandem would be a wrong span rather than a dropped record — they, and any other DUP subtype, are counted as an unsimulated SVTYPE instead. INV and BND have no subtypes in the spec, so a colon after either is likewise unrecognised.

spike acts on neither `FILTER` nor `GT`, so a non-PASS or `0/0` record is simulated like any other. It counts them and says so, so that is visible rather than assumed:

```
INFO VCF ingest ignores FILTER and GT: 1 record(s) it read are not PASS and 0 are homozygous reference; all are simulated like any other
```

### CRAM input

spike supports CRAM files transparently — just pass a `.cram` file instead of `.bam`:

```bash
spike --bam sample.cram --reference GRCh38.fasta \
  --event "del:chr17:43045000-43055000" \
  -o output/
```

The `--reference` FASTA is required for CRAM decoding (it is also required for haplotype construction, so there is no extra burden). The index must be the input's name with `.crai` appended (`sample.cram` → `sample.cram.crai`, what `samtools index` writes); the `sample.crai` spelling is not looked for.

Every region query on a CRAM — read extraction, the LOH pileup that phases heterozygous SNPs, and the `spike validate` checks — seeks only the containers that can overlap the requested region. noodles otherwise decodes every container on the chromosome: on a 10 Mb chr20 CRAM, a 30 kb deletion took 50.4 s and now takes 3.1 s, against 1.9 s from the equivalent BAM.

A CRAM container can hold several contigs at once (`samtools view -C --output-fmt-option multi_seq_per_slice=1` writes them that way, and htslib does it by itself on files with many short contigs). Such a container is decoded as a whole, and noodles then returns every record in it whose *position* falls in the queried window, whatever contig the record is on. All three query paths above now drop those foreign records. On a chr20+chr21 CRAM merging the full HG002 reads for both chromosomes over the same window, a 30 kb chr20 window yielded 7875 read pairs before and yields 4132 now, and the LOH pileup called 10,510 heterozygous SNPs in it and now calls 15. On a second CRAM built to put the chr21 reads only inside the deletion (so the leak's effect on coverage shows up clearly), `spike validate`'s `coverage_ratio` check (event depth ÷ flank depth) read 1.77x for a deletion that should show *reduced* coverage, and now reads 0.91x — matching the same check run on a chr20-only control. Each "now" figure is exactly what the identical reads give from a chr20-only CRAM. A read on the queried contig whose *mate* lies on another one is not foreign and is kept, the same as from a BAM; read extraction additionally requires both mates on the queried contig, because it builds pairs.

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

`--align` passes `--reference` straight through to the aligner as its index
prefix, and `align.sh` does not build that index itself — it expects one
already present at that exact path. A bgzipped `--reference` now loads for
simulation, but `--align` still fails unless a bwa-mem2 index was built
under that same bgzipped name (`bwa-mem2 index ref.fa.bgz`); an index built
from a differently-named, e.g. uncompressed, copy of the reference will not
be found.

Or run the generated script manually:

```bash
bash output/align.sh                    # Uses defaults from spike run
bash output/align.sh /path/to/ref 8    # Override reference and thread count
```

The defaults baked into `align.sh` and `merge.sh` are absolute paths, resolved
when spike writes the scripts, so a relative `--reference` or `--bam` still
works when the script is run from another directory. They are also
single-quoted, so a path holding a space, `}`, `"`, `$`, a backtick or a single
quote is passed through unchanged rather than being expanded or breaking the
script. `--samtools` is resolved the same way only when it contains a `/`; a
bare command name is left alone for `$PATH` to find. A custom `--aligner` is a
command line, not a path, so it is emitted verbatim -- quote it yourself if it
contains anything the shell would act on.

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

`align.sh` tags the simulated reads `@RG ID:sim SM:<sample>`, where `<sample>` is the `SM` of the original BAM's first `@RG` line, so `merged.bam` stays single-sample. If the original BAM's read groups carry different `SM` values it is already multi-sample; the first one still wins and spike logs a warning. A BAM with no `@RG SM` at all falls back to `SM:SIM`. The generated scripts quote the sample name, so one holding a space or an apostrophe (`SM:Patient 123`) reaches the aligner intact and keeps matching the original read groups; control characters and a backslash are replaced with `_`, because a tab ends the `SM` field and a newline ends the `@RG` line whatever the quoting, and bwa-mem2/minimap2 unescape `\t`/`\n` inside the `-R` string themselves -- a shell cannot quote against that.

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

The same table of checks-by-event-type is printed by `spike validate --help`.

The three `[global]` checks — `insert_size`, `dup_rate` and `mean_mapq` — are
computed from one sample of the reads in the truth events' own windows (event
+/- `--flank`), read with the same indexed query the per-event checks use, not
from the head of the file: the first 100k records of a whole-genome BAM are
chr1's telomere, mean MAPQ 10.0, which fails a check the rest of the file
passes. The sample is capped at 200,000 records and **every** event gets an
equal share of that cap (200,000 / number of events, at least one record each),
so no event is left out however many the truth VCF holds. `spike validate`
therefore needs the BAM's `.bai` / the CRAM's `.crai`, which its per-event
checks already required.

Each region's contribution and the size of the whole sample are logged at
`INFO` (`[global] sampled 117 records from chr20:37500081-37500320`,
`[global] sample: 117 records over 1 of 1 event regions`): a check computed
over 117 records prints exactly like one computed over 200,000, so read the log
before trusting a global number from a narrow `--flank`. Within one window the
sample is that window's first records in coordinate order, not a spread over
it, so an event longer than its share covers (above roughly 1 Mb at 35x with
one event) is represented by its start. A region that cannot be queried — a
truth event on a contig the alignment file does not have — is logged as a `WARN`
and skipped; the sample fails only if no region could be read at all.

A global check whose sample cannot answer it **fails** rather than passing on a
default. `dup_rate` distinguishes three cases when no sampled record carries the
duplicate flag: if a `@PG` record names a duplicate marker (`samtools markdup`,
Picard/GATK `MarkDuplicates`, `sambamba markdup`, `bammarkduplicates`,
`umi_tools dedup`) the rate is a measured `0.0%` and passes — the pipeline
decided those reads are not duplicates; if fewer than 1,000 records were sampled
it reports `too few reads`, because at a 1% duplicate rate a 41-record window
holds no duplicate about two times in three; otherwise it reports `no dup flags`
and asks for duplicates to be marked. A window holding no properly-paired record
reports `insert_size` as `no pairs`. Each also logs a `WARN` naming the reason.
A truth VCF with no events leaves no window to sample, so all three global
checks fail rather than reporting whole-file statistics for a run that validated
nothing.

A truth event **no check applies to** is reported as a failed `event_checked`
result rather than left out: without this, a truth VCF of such events scored
`3/3 PASS` on the three global checks alone, having verified nothing about the
events it was given. Which check covers which type:

| Truth event | Per-event checks |
| --- | --- |
| DEL, DUP | `coverage_ratio`, `split_reads` |
| INV, BND | `split_reads` |
| INS | `ins_reads` |
| SNP (single-base REF and ALT) | `allele_freq` |
| anything else (e.g. `SVTYPE=CNV`) | none -- `event_checked` FAIL |

`ins_reads` counts reads whose alignment **leaves the reference at POS**: an
insertion has no second breakpoint and no reference span, so neither
`coverage_ratio` nor `split_reads` can see it, but an aligner still has to put
the inserted bases somewhere -- an `I` CIGAR operation when the insertion fits
inside a read that anchors on both sides, a soft clip at the insertion point
when it does not. Either counts, if it is at least `min(SVLEN, 50)` bases long
and its reference boundary is within 100 bp of POS, and the check passes at
two such reads (the same threshold `split_reads` uses: one clipped read is
background anywhere, two at the same point are not). An INS record with no
usable `SVLEN` has no length to look for and is a failed check.

Measured on a DEL+INS run on the HG002 chr20 slice, aligned with `align.sh`
and merged with `merge.sh`: **21** reads carry the planted 300 bp insertion at
chr20:39000000, against **0, 0, 0, 0 and 1** at five control positions in the
same BAM where nothing was planted. Before this check existed, spike's own
round trip could not succeed for insertions -- the same run scored `5/6 PASS`
and exited **1** on the `event_checked` row, and now scores `6/6 PASS` and
exits 0.

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

### Contig names containing `:`

Some references name contigs with colons in them -- GRCh38's HLA alleles
(`HLA-A*01:01:01:01`) are the common case -- so the first `:` in an `--event`
or `--region` argument is not reliably the end of the contig name. spike
resolves the name against the reference's `.fai` contigs first, taking the
longest one that matches, and splits on `:` only when none does:

```bash
spike --bam sample.bam --reference GRCh38_full_analysis_set.fasta \
  --event "del:HLA-A*01:01:01:01:1000-2000" \
  --region "HLA-A*01:01:01:01:500-3000" \
  -o output/
```

A name the reference does not list is still split at its first `:`, so a typo
in a contig name fails the way it always did, and a gene name is never shadowed
by a contig. BND ALT strings from `--vcf` (`N[HLA-A*01:01:01:01:200[`) need no
contig list: their last `:` is always the one before the position.

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

      --vcf-info-af
          Read INFO/AF from --vcf records as the allele fraction to simulate.
          
          Off by default: in a population VCF (gnomAD, 1000G) AF is the allele frequency in the population, not the fraction of this sample's reads that should carry the allele, so using it silently produces a truth set at the wrong VAF. Without this flag only SIM_VAF and VAF are read, and a record that carries AF but has no usable SIM_VAF or VAF falls back to --allele-fraction, with a count of how many did on stderr.

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

spike generates fixed-length paired-end reads and needs a fragment at least as long as one read; it caps the fragment lengths it draws at 1500bp, so it does not support long-read (PacBio/ONT) libraries. If the input BAM/CRAM's mean read length exceeds 1500bp, spike exits with an error before doing any work rather than generating reads that don't match the library.

For the same reason spike needs a **paired-end** library: every donor read comes from an extracted pair, and the quality profile is trained on R1 and R2 separately. A single-end BAM/CRAM yields no donor material, so spike exits with an error naming the file rather than emitting a handful of synthetic pairs built on an empty quality profile. Single-end is detected from SAM flag `0x1`, which a paired library sets on every read whether or not the pair aligned properly — a paired BAM with no proper pairs in it is still accepted, and falls back to a 350 ± 50 bp insert size. Reaching 50,000 primary records with no `0x1` on any of them **is** how the scan establishes that, so it stops there instead of reading the file to the end. The message says whether the cap bound: at end of file the count it reports is every primary record in the file and the verdict is certain, while at the cap it is a window, and the message says so — a file that does hold pairs but puts 50,000 consecutive primary records without `0x1` at its head would be reported the same way. That is bounded rather than impossible: spike's input is a coordinate-sorted indexed BAM/CRAM, where a mixed library's single- and paired-end reads interleave at every locus, so no supported input is known to reach it.

### Deletion (DEL)

**Haplotype structure** (2 segments):
```
[left_flank]            [right_flank]
ref[start-F .. start]   ref[end .. end+F]
```
where `F` = the haplotype flank (`HAP_FLANK` in `src/main.rs`, 2000bp) — a fixed internal constant, not the user-configurable `--flank` option (which controls how wide a window of *original* reads is extracted around the event, and defaults to 10kb). The deleted region `[start, end)` is absent from the haplotype; the left and right flanking segments are placed adjacent.

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

Whether the gVCF can hold SNPs for the chromosome at all is read from its own `##contig` names, not guessed from an empty result: `bcftools view -r chr20:...` on a file that calls that chromosome `20` prints nothing and exits 0, exactly like a region that genuinely has no SNPs. spike reads the header itself (plain or bgzipped, no index needed) and warns only when the names cannot match — the other convention (`20` for `chr20`) or no such chromosome at all. A region that simply has no SNPs in it is not warned about. Each message says what happens next: a fallback to pileup when the file was read, and — when it could not be read at all, e.g. an unindexed `.vcf.gz` that bcftools refuses, with bcftools' own reason quoted — that LOH is skipped for that region and its original reads are suppressed at random instead.

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

Both mates are generated in **sequencing order**, 5'→3' along the read. The reverse mate's template is complemented and walked right to left along the reference before generation, rather than being generated along the reference and reverse-complemented afterwards, so its Markov chain runs with the cycle counter like the forward mate's instead of against it.

Error rates are derived from the sampled quality scores: `P(error) = 10^(-Q/10)`. When an error occurs, a random incorrect base is substituted.

Every `N` a synthetic read emits is reported at **Q2**, whatever put it there — an `N` in the reference, padding for the part of a read that runs past a contig end, or padding after a deletion sequencing error exhausted the template. An `N` is a no-call, and a real Illumina no-call is always Q2 (all 682 `N` bases in 173,338 HG002 NovaSeq reads on chr20 are Q2); the learned profile knows nothing about `N` and would otherwise hand one an ordinary score, usually Q37. The Q2 is what the Markov chain carries forward, so the base after an `N` is sampled from the after-a-Q2 transition bin.

The profile is still sampled for an `N` and the result discarded, so the *quality* draw stays one per template base — but an `N` skips the `P(error)` draw a called base makes, so it consumes strictly fewer random numbers than a called base. The stream is not unchanged: two runs that differ only in whether one template base is `N` diverge from that base onwards. (Level 4 above has one further last resort, a fixed Q20 for a cycle past the end of the profile, which makes no draw at all.)

Carrying Q2 forward costs the bases *after* an `N` nothing in practice. Measured over 40 seeds at a 1 bp reference gap on chr2 (an `N`-containing read there is 98.3% real sequence), the non-`N` bases of `N`-containing reads average **Q35.621** when the chain carries Q2 and **Q35.615** when it does not — a difference of +0.006 Q against a 0.022 standard error. The reason is that real Illumina Q2 is rare enough (6 of 389,429 donor bases in that window) that the after-Q2 transition bin never reaches the 30-observation threshold, so sampling falls straight through to the non-Markov levels.

### Too few donor reads

Every simulated read is built from the donor pool extracted for its event -- the quality profile above, the fragment-length distribution and the coverage the tiling count is scaled by all come from it. A pool holding fewer than **30 read pairs** after deduplication is refused: the run exits non-zero, naming the event and the windows it searched, and writes nothing.

```
Error: event DEL  chr20:30000001-30010000 (10000bp) has no usable donor reads: 0 read pair(s)
extracted from chr20:29990000-30020000, fewer than the 30 spike needs. ...
```

Without that check a starved window is silent. spike logs `Built read pool: 0 pairs`, then falls through to every substitute in turn -- the constant Q20 last resort above for every base, the default 400 +/- 80 fragment distribution, coverage 0 with the 2-read tiling floor -- and exits **0** with a truth VCF and two invented read pairs beside it. An event in a zero-coverage region, an off-target panel BAM and a mistyped `--region` all reach it.

30 is the observation count the quality model itself requires before it will sample from a bin, and a pool of *n* pairs puts *n* observations in each cycle-only bin (level 4 above), so it is the smallest pool at which any level of the model is trained to its own threshold. It is a floor on "measured from this library at all", not a coverage requirement: the 30 kb window of `del:chr20:38412500-38422500` yields 4559 pairs on the 35x HG002 BAM, so it would have to fall to roughly 0.2x before 30 pairs bound. If a real event does sit in a region that thin, widen `--flank` or `--region`, lower `--min-mapq`, or use a BAM that covers it.

### No donor coverage at the breakpoint

The 30-pair floor above counts the pool as a whole, summed over every window
the event was extracted from. That is not the same question the tiling count
asks: the number of synthetic fragments is `coverage x VAF`, and `coverage` is
measured in a 2 kb window around the *first breakpoint* only. A pool can clear
30 pairs and still measure coverage 0 there -- `--region` pointing somewhere
the event is not, or a fusion whose other partner carries the whole pool. The
tiling count then collapses to its floor of 2, and spike used to write those 2
invented pairs and a truth VCF beside them and exit **0**.

A breakpoint with no donor coverage is now refused the same way a starved pool
is: the run exits non-zero and writes nothing.

```
Error: event chr20:30000000-30010000 has no donor coverage at its first breakpoint
chr20:29999999: the pool holds 6117 read pair(s) but none of them cover that
position. ...
```

When the coverage is real but thin, the floor still applies -- a haplotype
shorter than one fragment asks for 0 fragments however good the coverage is,
and planting nothing would leave a truth VCF with no reads behind it. But the
floor then plants *more* support than `--allele-fraction` asked for, while
`SIM_VAF` in the truth VCF still records the request, so spike warns with both
numbers:

```
WARN spike::simulate] coverage 0.7x at VAF 0.050 asks for 0 tiled fragment(s); spike
emits the 2 it needs to plant the event at all, so the realized allele fraction will
be above the 0.050 recorded as SIM_VAF in the truth VCF
```

### Missing or unusable donor base qualities

SAM's QUAL field is all-or-nothing per record: a read either has a full quality string or none at all (`*`). The two containers encode `*` differently and spike sees both — BAM keeps a per-base array with every byte `0xFF`, while noodles normalises a CRAM record's all-`0xFF` buffer to an *empty* one before spike ever sees it. A record with either shape, or carrying a raw quality above Q93 (the SAM maximum — a malformed record), is dropped during extraction: it is not included in the learned quality profile and does not contribute a read pair to the output. Dropped records are counted and logged as a warning (spike logs to stderr), e.g.:

```
WARN spike::extract] 457 record(s) considered for the chr20:38402500-38432500 donor pool had no quality scores (SAM '*'); their pairs were dropped from the pool and merge.sh removes them from the merged BAM
```

**A dropped pair is removed from the merged BAM even though nothing replaces it.** spike has no quality string to write for it, so it cannot come back through the FASTQ — but it is still listed in `replaced_reads.txt`, so `merge.sh` drops it. That costs real depth. Leaving it in costs more: the pair would sit inside every event it overlaps as reference support that no allele fraction can suppress and no synthetic read can replace, so the realised VAF would come out diluted by the unusable-quality fraction while the truth VCF still claimed the full one. Because the same pairs are dropped across the whole extraction window, not just inside the event, the loss cancels in a flank-normalised ratio; the dilution would not. Measured on an HG002 chr20 slice with 5% of donor pairs' quality stripped to `*`, over a 10 kb deletion: 147 of 2868 eligible records inside the event (5.13%) carried `*` quality and were un-suppressible; removing them takes that to 0.00%, at a cost of those same 147 records of depth.

Without the drop, a missing quality used to decode to an invalid FASTQ quality byte (space) for kept reads and poisoned the learned quality model, so synthetic reads sampled from it came out mostly-Q0 with effectively random bases. On the CRAM path it produced a FASTQ record with a full-length SEQ line next to a zero-length QUAL line, which `samtools import` rejects outright (`truncated file`).

As a last line of defense, `write_paired_fastq` validates every pair *before* it creates either output file — so a refusal leaves no half-written `.fq.gz` in `--output` — and returns an error if a quality string's length does not match its SEQ, or if any byte falls outside the printable Phred+33 range `!`-`~` (33-126).

`write_paired_fastq` also flushes each gzip stream's returned writer explicitly after `finish()`, rather than only letting it drop. `finish()` empties flate2's own internal buffer into the `BufWriter` it hands back, but that write can land entirely inside the `BufWriter`'s own buffer without ever reaching the file — a write failure (e.g. a full disk) then surfaces only when the `BufWriter` is dropped, where `Drop`'s own flush swallows any error. Calling `write_paired_fastq` with `R1.fq.gz` symlinked to `/dev/full` and a single small read pair returned `Ok(...)` before this fix and now returns `Err("No space left on device (os error 28)")`. Both streams' `finish()`+`flush()` are computed unconditionally before either error is propagated, so an R1 failure can't short-circuit past R2's — leaving R2's `GzEncoder` to be dropped unfinished would rediscover the exact same swallowed-error bug for R2 alone; when both streams fail, the reported error now names both.

### Indel error model

When `--indel-error-rate` is set above 0, a fraction of sequencing errors are modeled as insertions or deletions (50/50 split) rather than substitutions. This maintains fixed read length: insertions consume an output position without advancing the reference, and deletions skip a reference base without consuming an output position.

Because a deletion error consumes a template base without emitting one, the read runs past the end of the window it would otherwise need. Generation therefore fetches 10 extra template bases past each read's 3' end — to the right of a forward read, to the left of a reverse one — so those bases are real sequence. `N` padding is now only a last resort at a contig or haplotype end.

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

### End-to-end validation harness

`scripts/validate_pipeline.sh` runs the whole loop on real data: it takes het
DELs from the GIAB HG002 T2TQ100 truth set, spikes them into a clean
1000 Genomes background BAM at several VAFs, aligns and merges with the
generated `align.sh` / `merge.sh`, runs `spike validate` and Delly, scores the
calls with Truvari, and runs the same caller on the *unspiked* background so the
verdict can be attributed to the spike-in rather than to the background sample.

```bash
# Full chr20 run
bash scripts/validate_pipeline.sh \
  --giab-dir /path/to/giab_hg38 \
  --background-bam NA18488.chr20.bam

# Small, fast run: one 2.3 Mb window, two allele fractions
bash scripts/validate_pipeline.sh \
  --giab-dir /path/to/giab_hg38 \
  --background-bam NA18488.chr20.bam \
  --region chr20:61900000-64200000 \
  --vafs "0.5 0.25" \
  --outdir /scratch/spike_validation
```

Output goes to `--outdir` (or `$OUTDIR`), which defaults to
`<repo>/validation_run` (untracked). A relative `--outdir` is resolved against
the current directory. The script refuses to run when it names a directory
whose contents git tracks, so it cannot overwrite committed fixtures.

Requires `samtools`, `bwa-mem2`, `bcftools`, `bgzip`, `tabix`, `delly`,
`truvari` and `python3`. Each is taken from `PATH` and can be overridden with an
environment variable of the same name (`SAMTOOLS=...`, `DELLY=...`, ...); the
GIAB paths default to `<repo>/data/giab_hg38` and can be moved with
`--giab-dir` or with `GIAB_DIR` / `REFERENCE` / `TRUTH_VCF` / `BENCH_BED`
(inside a git worktree `data/giab_hg38` is a dangling symlink, so pass
`--giab-dir` there). The output directory can also be set with `OUTDIR`.

Options: `--region chr:beg-end` restricts the run to one window (it slices the
background BAM, so the whole run stays small), `--max-events N` caps the number
of truth DELs, `--min-events N` is the floor below which the run aborts
(default 5), `--vafs "0.5 0.25 0.1"` sets the allele fractions, `--min-recall F`
fails the run when any VAF recalls less than `F`, and `--min-gain N` is how many
truth events the highest VAF must recover *beyond the background control*
(default 1). `--min-gain 0` is the weakest setting, not an off switch: the
highest VAF must still at least match the control, and a run with no control to
compare against still fails. `--skip-to N` resumes an existing `--outdir` at
step N (1-8); a value outside that range is refused rather than silently
skipping the whole pipeline.

The background BAM must be aligned to the same reference: Delly refuses a BAM
whose header names contigs the FASTA does not have (an `_alt` background next
to a `no_alt` reference, say), so step 1 checks that up front and stops.

Read the recall column against the background, not against zero: the truth DELs
are common HG002 variants, and a 1000 Genomes background often carries the same
ones, so Delly recovers some of them with no spike-in at all. On a 2.3 Mb chr20
window (8 truth DELs) Delly finds 3 of them in the *unspiked* background — three
of the four recovered at VAF 0.5. Step 7b measures that floor for you: it runs
the same `delly call` and the same `truvari bench` on the background BAM alone
and writes the result as the `background` row of the summary table. Its calls
are cached in `<outdir>/background_control/`, keyed on the background BAM, the
reference, `--region` and the truth VCF, so re-running the same outdir with a
different window or truth set recomputes the floor instead of scoring a fresh
spiked number against a stale one.

The script **exits non-zero** when it did not validate anything: a missing tool
or data file, fewer than `--min-events` truth events, a `spike` or `merge.sh`
error, an unparseable `spike validate` report, a missing Truvari summary, a
missing background control, or a highest-VAF result that does not recover at
least `--min-gain` more truth events than that control. The verdict is reached
whichever steps ran: `--skip-to 8` over an outdir whose highest VAF has no
Truvari output fails there, rather than printing a row of `N/A` and passing. That last gate is the
one that makes the verdict mean something: a plain "recovered more than zero"
test passes on the background alone, so a run in which the spike-in contributed
nothing would still print `VALIDATION PASSED`. A summary table is written to
`<outdir>/validation_summary.tsv`, and the control's own calls to
`<outdir>/background_control/`.

The gate counts recovered truth events, so it catches a run that contributed
nothing at all; on a single replicate it cannot separate a very weak spike-in
from re-alignment noise. Measured on the 2.3 Mb window: at VAF 0.5 the spiked
run recovers `sim_del_6` (chr20:63636171) on top of the control's three — the
one truth DEL the background does *not* carry, which is the right answer; at
`--vafs 0.001` it also recovers one more than the control, but that one is
`sim_del_8`, a DEL the background carries and that re-alignment flipped from FN
to TP. Separating those needs replicates or a VAF titration, not a single run.

Note that `spike validate`'s per-event checks are *not* part of the verdict
beyond "the report parsed and at least one check passed". On a cross-sample
spike-in most of them compare against a background that already carries the
event, so they fail for reasons that have nothing to do with the injection —
see `N1` in `REVIEW.md`.
