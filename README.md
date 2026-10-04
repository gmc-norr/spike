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

**Read suppression and replacement**: For non-additive events (DEL, INV, INS, SNP, full-model DUP), original reads within the haplotype's reference footprint are suppressed at the target VAF rate, and new synthetic reads tiled across the variant haplotype replace the removed fraction. For additive events (Fusion, junction-model DUP), all original reads are kept and synthetic reads are added on top. By default a non-additive event removes only reads from its donor pool; the experimental `--edit-model origin` removes reads by their chance of having come from the event, at the event and at its look-alikes elsewhere in the genome.

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

### Test-time prerequisites

What `cargo test` needs on PATH is not the same list as what a spike *run*
needs:

- **bcftools** — required to run the suite, unconditionally, even though a
  run needs it only for `--gvcf` with a `.vcf.gz`. Two tests in `src/loh.rs`,
  `loh::tests::test_a_renamed_gvcf_that_cannot_be_read_warns_about_the_stop_not_the_pileup`
  and `loh::tests::test_a_gvcf_read_that_fails_says_the_run_stops`, read a
  `.vcf.gz` the way spike does, through `bcftools view`. Without bcftools both
  fail, by name and on purpose: a missing test dependency is a broken
  environment, not a test to skip.
- **bash**, for the tests that run the generated `align.sh` / `merge.sh` and
  `scripts/validate_pipeline.sh`. Those tests write their own stub `samtools`
  and stub aligner and put them on the script's PATH, so a *real* `samtools`,
  aligner, `bgzip`, `tabix`, `delly` or `truvari` is **not** needed — measured:
  with all of them off PATH the suite is `667 passed; 2 failed; 1 ignored`, the
  two failures being the bcftools tests above.

No reference FASTA, BAM or CRAM is needed for `cargo test`: the tests build
the tiny inputs they need. The one `#[ignore]`d measurement
(`synth::tests::measure_n7_quality_drift`) does need an HG002 BAM, named by
`SPIKE_N7_BAM`, and is not run by default.

### Build

```bash
cargo build --release
# Binary at target/release/spike
```

### Measuring a change: one target dir per commit

A before/after comparison is only as good as the two binaries behind it. In
this repo's own review work, three different source trees built with one
shared `CARGO_TARGET_DIR` gave three **byte-identical** binaries, and cargo
printed `Finished ... in 0.06s` each time. The comparison would have shown no
effect for a change that has one, with no error anywhere. (A milder form of
the same trap: measuring with a `target/release/spike` that was not rebuilt
after the last edit.)

So give each commit its own target dir, and check the binaries differ before
trusting the numbers:

```bash
CARGO_TARGET_DIR=/tmp/spike-before cargo build --release   # at the "before" commit
CARGO_TARGET_DIR=/tmp/spike-after  cargo build --release   # at the "after" commit
md5sum /tmp/spike-before/release/spike /tmp/spike-after/release/spike
```

Two equal checksums for two commits that change code mean a stale build, not
a zero effect. See `NF4` in `docs/review/CLINICAL_SV_NEW_FINDINGS.md`.

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
#         output/fastq_removed_reads.txt, output/fastq.sh, output/README.md
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

**A site the sample already carries is refused.** spike edits one copy of the sample and keeps the other copy's reads as they were. If the sample already has another allele at the bases a small variant changes, those reads keep it, and the requested fraction cannot be reached. Asked for at 0.5, a site where the sample is hom-alt for the same change came out 25 ALT and 0 REF reads, while `truth.vcf` said 0.5 (review finding 3).

So before any other work, spike counts the reads at each small variant's own bases. It uses the pileup's filter: primary, not duplicate, not QC-fail, MAPQ at least `--min-mapq`. A read that spans the site carries another allele when it has, inside the changed bases, a base other than the reference, a deletion, or an insertion at either edge.
- spike refuses the event when at least 10 reads span the site and at least a fifth of them carry another allele. These are the pileup's own floor and its het bound.
- Every refused event is listed in one error, and nothing is written.
- With fewer than 10 spanning reads, spike warns that it could not check, and goes on.
- Structural events (`del:`, `ins:`, `dup:`, `inv:`, fusions) are not checked.

On the hospital HG002 30x BAM (chr20), the rule refused:
- 99.3% of the het SNVs, 100% of the hom SNVs, 97.5% of the het indels and 98.5% of the hom indels that the sample's own DeepVariant calls hold;
- 0.2% of random SNVs and 0% of random 1-10 bp indels where it has no call;
- 0% of HG001's SNVs and 0.5% of its indels.

The plan and the numbers are in `docs/superpowers/plans/2026-10-04-carried-allele.md`.

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
- `af=het` — exactly 0.5: one copy of two (germline heterozygous), the same as
  `af=0.5`. Earlier versions drew it from Beta(40,40), which moved the event
  itself, not just what the reads show: above 0.5 the other copy lost reads
  too (CR7(c) in `docs/review/CLINICAL_SV_DESIGN_NOTES.md`).
- `af=hom` — fixed at 1.0 (homozygous)
- *(omitted)* — uses the global `--allele-fraction` (default: 0.5)

By default, events on the same chromosome that land on each other are rejected to keep effects independent. Use `--allow-overlap` to override (they are then simulated independently and merged, which is approximate where they meet).

What counts as landing on each other is the **replacement footprint**, not the event span. Each event is simulated on its own against the original BAM, and it replaces reads — removes the originals and tiles synthetic ones — across its span *plus 3500 bp on each side*: 2000 bp of haplotype flank (the reference spike carries either side of the edit, `HAP_FLANK`) plus 1500 bp for the longest fragment it can emit (`MAX_FRAGMENT_LEN`). So an event's footprint is `start - 3500 .. end + 3500`, and **two events must leave at least 7000 bp between the end of one span and the start of the next.** Closer than that and each event's synthetic flank writes plain reference over the other event's edit, undoing it: two 1 kb homozygous deletions 1 kb apart recovered 70x and 74x of a 75x baseline instead of 0x. The error names both events, both spans as you gave them, and both footprints:

```
Error: overlapping events detected (default is to reject overlaps).
Each event replaces reads across its span grown by 3500bp on each side (2000bp of haplotype flank + 1500bp of fragment), and two such replacement footprints may not intersect: two spans on one chromosome must be at least 7000bp apart.
Use --allow-overlap to override.
  - events 1 and 2 have intersecting replacement footprints on chrT: spans 10000-11000 and 12000-13000 (footprints 6500-14500 and 8500-16500)
```

`--allow-overlap` still accepts them with a warning, and the composition is still approximate there — it does not make nearby events compose correctly.

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

- Records inside the event regions that spike never extracted — duplicates of the pairs it kept, non-proper pairs, reads below `--min-mapq`, reads whose mate was filtered — are **kept**. Removing by region deleted them without putting anything back, which cost real depth.
- **Duplicates of a pair spike removed are removed with it.** A duplicate is another read of the same molecule, so it goes where that molecule goes. spike never takes a duplicate as donor material. After the run, it reads the event windows again (± 1,500 bp) and lists every duplicate pair that shares both 5' ends with a pair it removed. Left behind, they would do no harm in `merged.bam`, where they keep their flag. But a pipeline that marks duplicates again from scratch, as one starting from the [full FASTQ](#full-fastq) does, finds nothing for them to copy and counts one per family as a read of the old allele. Measured before this fix: 3.4% of the reads at spiked hom SNVs, against 0.17% at real ones, and 2-3 reads inside hom deletions (`docs/superpowers/plans/2026-10-04-duplicates.md`). Duplicates of a pair spike kept stay, with their flag.
- Mates that lie outside the event regions but whose pair spike did extract are **removed**, because `sim.bam` carries their replacement. Removing by region left them in, so they appeared twice.
- Pairs spike extracted but could not use — a mate with no stored base qualities, see [Missing or unusable donor base qualities](#missing-or-unusable-donor-base-qualities) — are **removed without replacement**. That costs depth; left in, they would dilute the realised VAF of every event they overlap.

`ORIGINAL_BAM` (the first argument to `merge.sh`) must be the exact BAM spike was run on: `replaced_reads.txt` names reads spike found in that BAM, and `-N` only removes names it can find. A different BAM shares essentially no read names, so nothing would be removed and `sim.bam` would be merged onto full, un-thinned original depth — `merge.sh` guards against this by comparing how many of the expected names it actually matched against how many `replaced_reads.txt` lists, and aborts with an error instead of silently producing a wrong `merged.bam`.

Records spike extracted and then suppressed (the deleted copy of a heterozygous deletion, for instance) stay gone: that absence *is* the simulated variant.

Under `--edit-model origin`, `replaced_reads.txt` also names every read `origin` removed, so `merge.sh` can remove records `clean` keeps: reads below `--min-mapq`, and reads at a look-alike elsewhere in the genome. Duplicates of a removed pair go under either model. See [Editing hard spots](#editing-hard-spots---edit-model-origin-experimental).

The records kept because spike never extracted them (duplicates of the pairs it kept, non-proper pairs, low-MAPQ or orphaned-mate reads) are real original reads that now sit inside an event's footprint, so an event's residual depth/allele fraction in `merged.bam` is no longer exactly the simulated value. Measured on an HG002 chr20 run: `validate.rs` only skips secondary/supplementary/duplicate/QC-fail and low-MAPQ reads — it has no proper-pair or mate-unmapped filter — so of the 1,436 records recovered by this change, the 87 non-proper-pair and 39 orphaned-mate records (126 total, 1.2% of the 10,501 in-BED records) reached an AF or depth measurement in `spike validate`; a consumer that counts duplicates rather than skipping them could see the residual shift by up to the full recovered fraction (13.7%).

`align.sh` tags the simulated reads `@RG ID:sim SM:<sample>`, where `<sample>` is the `SM` of the original BAM's first `@RG` line, so `merged.bam` stays single-sample. If the original BAM's read groups carry different `SM` values it is already multi-sample; the first one still wins and spike logs a warning. A BAM with no `@RG SM` at all falls back to `SM:SIM`. The generated scripts quote the sample name, so one holding a space or an apostrophe (`SM:Patient 123`) reaches the aligner intact and keeps matching the original read groups; control characters and a backslash are replaced with `_`, because a tab ends the `SM` field and a newline ends the `@RG` line whatever the quoting, and bwa-mem2/minimap2 unescape `\t`/`\n` inside the `-R` string themselves -- a shell cannot quote against that.

`merged.bam` is appropriate for end-to-end testing where the caller needs to see the full genome (e.g., tools that estimate background noise from off-target regions). `sim.bam` is sufficient for targeted callers or focused benchmarking.

#### Reads spike cannot edit

spike counts those kept reads for every event and reports them (CR4). Which reads it edits does not change; it only says how many it could not. Per event, the log prints `Reads over the event spike cannot edit: <resistant> of <counted>`, the run README gets a **Resistant reads** column, and `truth.vcf` records the share as `SIM_RESIST`, which `spike validate` reads back as its advisory `resistant` row. Above **0.10** spike also warns, and the run README names the event:

```
DEL  chrT:10001-14000 (4000bp): 1038 of 2074 reads over it (50%) are ones spike cannot edit (below --min-mapq, not a proper pair, or a mate that fails a filter). They stay in the merged BAM as they are, so the event is weaker than requested; truth.vcf records the share as SIM_RESIST (CR4).
```

That is the review's `lowmap` probe, half of whose pairs are at MAPQ 0: a deletion asked for at AF=1 keeps half its depth. (Its share is 1038 of 2074, just over one half, so since RF8 it takes `--allow-resistant` to run; see below.) On ordinary loci the share is small. Measured on 40 seeded 10 kb deletions inside the HG002 SV benchmark on chr20 (35x, `--min-mapq` 20): median **0.010**, largest **0.083**, and none above 0.10. The same draw on chr1: median **0.011**, largest **0.286**, and 3 of 40 above 0.10. Duplicates are not counted, because callers usually skip them; a consumer that counts them sees more surviving depth than `SIM_RESIST` says.

**Above 0.5 spike refuses (RF8).** The event that reaches the reads is about `VAF x (1 - SIM_RESIST)`, because the reads spike cannot edit are never removed or replaced. Above one half, the reads carry less than half of what `truth.vcf` would claim. So spike exits non-zero, writes nothing, and lists every such event at once (the message is in [What spike refuses](#what-spike-refuses)). `--allow-resistant` runs them anyway, byte-identical to what spike wrote before the refusal existed.

- **Random spots rarely hit it.** 0 of the 40 chr1 deletions above are over 0.5.
- **Real SV sites often do, because they sit in repeats.** `scripts/validate_pipeline.sh` spikes 20 of HG002's real chr20 deletions into NA18488. On that background **6 of the 20** are over 0.5 (0.611 to 0.945). At `chr20:61943514-61945040`, 188 of the 200 reads are below MAPQ 20. The script therefore passes `--allow-resistant` and logs how many of its events are over 0.5. How to score those events is not yet decided.
- **Lowering `--min-mapq` often fixes it, because the pool then holds those reads.** RF8's locus, `del:chr20:27100000-27110000` on the 35x HG002 BAM:
  - at `--min-mapq` 20, 5731 of 5749 reads are uneditable, and spike refuses;
  - at `--min-mapq 0`, only 58 of 5749 (0.010) are, and spike accepts it.

  On the pipeline's six, `--min-mapq 0` brings five under 0.5. The sixth, `chr20:63093346-63094243`, stays at 11 of 18. Those 11 are all at MAPQ 20 or above, but none is a proper pair (one has an unmapped mate), so MAPQ is not what keeps them out.
- **`--edit-model origin` is the other way** (next section). At `del:chr20:7119236-7120236` on the 35x HG002 BAM it takes the share from 0.940 to 0.016. Under it a read counts as editable when it is in the donor pool or in a pair `origin` may remove. The warning, the refusal and the run README say so, and add that the 0.10 and 0.5 thresholds were set for `clean` and are not yet checked for `origin`.

#### Editing hard spots: `--edit-model origin` (experimental)

`--min-mapq 0` makes the MAPQ 0 reads editable, but it removes each one as if it certainly came from the event. It also never touches the look-alike copy the aligner split those reads with. `--edit-model origin` models what a real variant does there instead:

- Every primary read over the event's footprint (the event ± 2 kb) is a candidate; at a look-alike, only a read whose `XA` names a position inside the footprint is. The chance of having come from the footprint comes from its MAPQ and its `XA` alternative hits.
- It is removed by that chance times its copy's rate, but only if spike's new reads can replace it: every mapped mate could have come from inside the footprint.
- Chances from several events add up. A duplicate shares its original's fate.
- The number of new reads comes from where reads came from (the origin depth). It does not come from the donor pool, which can be empty inside a perfect twin.

It needs the aligner's `XA` tags. bwa-mem and bwa-mem2 write them by default, and spike stops if none of the first 100,000 records of the BAM carries one. A spot whose MAPQ 0 reads lack `XA` is normal: under bwa-mem's `-h 5` rule their hits number more than 5, so each gets a chance of 1/6. At `chr20:7117236-7121236` in the 35x HG002 BAM, 822 reads are MAPQ 0 and 1 carries `XA`, yet the file's first 100,000 records hold 15,255 with it.

**What changes in the output:**
- **`replaced_reads.txt` also names every read `origin` removed**, so `merge.sh` takes them out of the whole BAM, including reads `clean` never touches. Measured on `del:chr20:14530000-14531000` over a 100 kb slice of the 35x HG002 BAM (chr20:14,500,000-14,600,000), again on 2026-10-04 after `clean` began removing duplicates too: `origin` lists 2,748 names and `clean` 2,736. 25 are `origin`'s alone: 24 duplicates, and 1 pair that is not proper and is below MAPQ 20. 13 are `clean`'s alone, all duplicates. Both models remove the duplicates of the pairs they remove, and they remove different pairs: each of `clean`'s 13 goes with a pair only `clean` removed, and 23 of `origin`'s 24 with a pair only `origin` removed. On the same slice the binary from just before that change listed 48 names more under `origin` (96 records), 47 of them duplicates. The figure here before, from master 55b95c4, was 41 names and 82 records, 76 of them duplicates.
- **The new reads are scaled by the origin depth.** Around each point, that is the summed chance of every read `origin` may remove, turned into fragment depth. The log prints it beside the pool's depth, e.g. `origin depth at chr20:7119235: 33.7x (the donor pool's there: 8.5x)`. With no origin depth at any breakpoint, spike stops (see [What spike refuses](#what-spike-refuses)).
- **The census means something else.** `SIM_RESIST` counts a read as editable when it is in the donor pool or in a pair `origin` may remove, and `SIM_DEPTH_FOLD` compares origin depths. Their thresholds (warn above 0.10, refuse above 0.5, warn above a 1.5-fold) were set for `clean` and are not yet checked for `origin`. The warnings, the refusal, the run README and a `##spike_edit_model=origin` header line in `truth.vcf` all say so. `spike validate` does not read that line, so its `resistant` and `depth_fold` rows use the same two thresholds under either model.

**At RF14's site**, `del:chr20:7119236-7120236` on the 35x HG002 BAM with default flags:
- `clean` refuses: 359 of 382 reads over it (94%) are ones it cannot edit.
- `origin` runs: 6 of 382 (`SIM_RESIST=0.016`), and `SIM_DEPTH_FOLD=1.30`. It removed 219 fragments at the spot and **0** at its 14 look-alikes.

**The look-alike half does little on real data so far.** Those 14 look-alikes all hold reads (3518 primary reads), so the BAM is not what is missing. Of the 52 fragments there with a read whose `XA` hit lies in the footprint, 25 have a mate that could not have come from the footprint, so the second bullet's rule keeps them. 18 more have a mate outside the regions spike reads. The other 9 pass, yet none was removed in that run. (The count is awk over `samtools view` using `XA` start positions only; see `.claude/judgment-gate-cases.md`.)

**It needs a whole-genome BAM.** A look-alike can be anywhere in the genome, and a BAM cut to one chromosome or region does not hold its reads. spike then warns, once per event. This is `del:chr20:14530000-14531000 --edit-model origin` on the 100 kb slice above:

```
WARN  spike::origin] origin: 2 of 2 look-alike region(s) of chr20:14528000-14533000 hold
no read in this BAM (e.g. chr13:59757249-59757716). A look-alike exists because reads at
the event list it as an equally good place, so a whole-genome BAM holds reads there; a BAM
cut to one chromosome or region does not. origin removes nothing from those regions, and
the event's origin depth misses the share their reads would add. Use a whole-genome BAM so
those regions can be read.
```

A look-alike on a contig the BAM's header does not list is skipped, with a warning of its own. That happens to a BAM aligned with alt contigs and then re-headered to the no_alt set. Over chr20:3-53 Mb of the pipeline's NA18488 background, the reads' `XA` hits name 1636 contigs, and 1462 of them are not among the header's 195.

`origin` takes longer than `clean` on a big event: see [Run time](#run-time).

On a made-up 80 kb genome built with an exact twin, `origin` reproduced the depth a real het deletion leaves at both the event and its look-alike copy, where `--min-mapq 0` did not and `clean` refused for want of donor coverage — a physics test whose rule was locked before the run returned SUPPORTED (`docs/review/REVIEW.md`, "`--edit-model origin`: physics test"). The default stays `clean` until `origin` is tested against real data. The design is in `docs/superpowers/specs/2026-09-26-edit-model-origin-design.md`.

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
  --strict         Count the advisory checks in the exit status
  --help, -h       Show this help
```

`spike validate --help` prints the same table of checks-by-event-type, with the
real checks alone: the advisory rows below are described here and in the
report itself, not in the usage text.

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
| DEL | `coverage_ratio`, `del_planted`, and the advisory `coverage_any_mapq`, `split_reads` and `split_reads_each_end` |
| DUP | `coverage_ratio`, `split_reads`, and the advisory `coverage_any_mapq` and `split_reads_each_end` |
| INV, BND | `split_reads`, and the advisory `split_reads_each_end` |
| INS | `ins_planted`, and the advisory `ins_reads` and `ins_sequence` |
| SNP, small indel and MNV (explicit REF and ALT) | `allele_freq` |
| anything else (e.g. `SVTYPE=CNV`) | none -- `event_checked` FAIL |

Six rows are **advisory**: printed and counted with the rest, but left out of
the exit status unless `--strict` is given (see below). `split_reads` is advisory
for a DEL as well, and only for a DEL (RF14).

`coverage_any_mapq` is `coverage_ratio` recomputed with **no MAPQ floor**. It
covers the same event types, DEL and DUP, and it is the same code over the same
`--flank` window, judged against the same expectation (`1 - VAF` for a DEL,
`1 + VAF` for a DUP) with the same 0.30 tolerance. The one difference is which
reads are counted: every primary, non-duplicate, non-QC-fail read whatever its
MAPQ, where `coverage_ratio` counts only reads at `--min-mapq` (default 20) or
above. Both rows are always printed -- including at `--min-mapq 0`, where they
are identical by construction -- because a report whose shape depended on a
flag would hide the comparison the row exists to make. On the one-contig chrA
CRAM this repo's own unit tests write -- `mapq_hidden_deletion_cram` in
`src/validate.rs`, not one of the review probes -- whose deletion window holds
nothing but MAPQ 0 reads at half the flank depth, the two rows read:

```
Event                               Check              Expected                  Observed        Status
DEL chrA:10000-11000 (hidden)       coverage_ratio     0.00                      0.00            PASS
DEL chrA:10000-11000 (hidden)       coverage_any_mapq  0.00                      0.50            FAIL (advisory)
```

`split_reads_each_end` is `split_reads` with its minimum required at **each**
breakpoint instead of pooled over both. It covers the same event types, DEL,
DUP, INV and BND, and it is built from the **same two queries** as the pooled
row -- the same 500 bp windows, the same `SA:Z` test, the same minimum of two
distinct read names -- so the two rows can never disagree about what evidence
exists. The pooled row counts the union of the two windows' read names and
passes on two *in total*; this row takes each window's own count and passes only
if each reaches two. Its `Expected` reads `>=2 at each end` and its `Observed`
is the two counts, `<at the event's own breakpoint>/<at the partner's>`. On the
one-contig chrA CRAM this repo's own unit tests write --
`one_sided_split_read_cram` in `src/validate.rs` -- whose middle junction has
three reads joining the left breakpoint to the right and nothing joining back,
the two rows read:

```
Event                               Check              Expected                  Observed        Status
DEL chrA:10001-12000 (one_sided)    split_reads        >=2 joining chrA:12001    3               PASS
DEL chrA:10001-12000 (one_sided)    split_reads_each_end >=2 at each end           3/0             FAIL (advisory)
```

That second line is two characters wider than it looks: `split_reads_each_end`
is 20 characters and the Check column is 18, so **this one row's own later
columns sit two characters right of every other row's**. The alternative was
renaming the row or widening the column, and widening it would move the spacing
of every non-advisory row, so the raggedness is deliberate. It affects no other
row's *alignment*: the column pads short names and only overflows long ones. It
does affect *parsing*. Every other column is truncated one character below its
width, so before this row every column boundary in every row was at least two
spaces wide and a reader splitting on runs of two or more spaces worked; the
check name is the one field passed through untruncated, and at 20 characters
`split_reads_each_end` is the first content ever to overflow its column. On that
one row exactly **one** space separates Check from Expected, so `awk -F'  +'` or
`re.split(r"\s{2,}")` merges those two fields and shifts every later field along
-- measured on the fixture run those two lines come from, where every each-end
row splits into four fields and every other row into five. Nothing in this repo
parses the table that way, so this costs nothing today. `--json` is the
parseable form: use it rather than the text table.

`ins_sequence` is `ins_reads` with the inserted bases actually read. It covers
INS events, the same type `ins_reads` covers, and it is the first check here
that verifies inserted **sequence** rather than a CIGAR: the inserted bases come
out of the truth record's own ALT (the anchor base and then the insertion, which
spike writes for a random sequence as well as an explicit one), and the row
looks for them in the reads' own bases. Two probe k-mers are taken from the
insertion, its **first** and **last** `k = min(inserted length, 31)` bases, and
a read supports the insertion if its bases hold either one in either
orientation. The reads counted are the distinct read names within 150 bp of POS
that pass the same filter `ins_reads` applies -- primary, mapped, not a
duplicate, not QC-fail, at `--min-mapq` or above -- and the minimum is the same
two. Its `Expected` reads `>=2 with a <k>bp alt kmer` and its `Observed` is the
supporting count. On the one-contig chrA CRAM this repo's own unit tests write
-- `inserted_sequence_cram` in `src/validate.rs` -- graded once against the
truth ALT the reads carry and once against a **different** sequence of the same
length, the two rows read:

```
Event                               Check              Expected                  Observed        Status
INS chrA:5000 (carried)             ins_reads          >=2 reads with >=50bp...  3               PASS
INS chrA:5000 (carried)             ins_sequence       >=2 with a 31bp alt kmer  3               PASS (advisory)
INS chrA:5000 (wrong_bases)         ins_reads          >=2 reads with >=50bp...  3               PASS
INS chrA:5000 (wrong_bases)         ins_sequence       >=2 with a 31bp alt kmer  0               FAIL (advisory)
```

(Measured before RF13. `ins_reads` has printed `(advisory)` since then, and an
`ins_planted` row comes first.)

Two cases the row cannot answer, each reported as a **failed** row rather than a
silent pass, and each costing nothing in the exit status because the row is
advisory. An insertion of fewer than 12 bases gives a probe k-mer too short to
mean anything inside a read -- a given random 12-mer is expected in about one
151 bp read in 10^5, an 8-mer in about one in 400 -- so the row reads
`Expected >=2 with a >=12bp kmer` and `Observed alt is <n>bp`. And if either
probe k-mer is already in the **reference** within 1 kb of POS, an unedited read
would match it and the count would mean nothing, so the row reads
`Observed kmer in ref`. Both log a `WARN` naming the reason. A random insertion
makes the second vanishingly unlikely; an explicit one copied from nearby
sequence does not.

A record whose ALT recorded **no sequence** gets **no row at all**, and that is
not a failure: a symbolic `<INS>`, or a bare anchor base, is an older spike's
truth VCF, written before the inserted sequence was put in the ALT, and there is
nothing for this row to look for. The same rule the census rows below apply to a
truth record carrying no census field.

The other two advisory rows are the **census spike recorded**, and a record gets
each one whenever it carries that row's field, whatever the event's type. `resistant` reports the record's
`SIM_RESIST` and `depth_fold` its `SIM_DEPTH_FOLD`, each against the threshold
spike already warns at: `<=0.100` for the resistant share, `<=1.50` for the
depth fold. A record carrying neither field -- an older spike's truth VCF --
gets neither row, and that is not a failure; a field that is there and cannot be
read is an advisory FAIL (the row's Observed reads `bad: <value>`), because a
census that cannot be read has not been checked. One value counts as absent
rather than unreadable: `.`, VCF's own missing value, which spike writes itself
whenever it has no number for that event -- grading it would make spike
advisory-FAIL its own output. An advisory row is never a check *of* the event
either: an event no check covers still reports `event_checked` FAIL beside its
two advisory rows.

Such a row prints as `PASS (advisory)` or `FAIL (advisory)` in the Status
column, counts in the `Result:` line with every other row, and is summarised on
a line of its own below it:

```
Result: 6/8 PASS
Advisory: 2 checks, 1 PASS, 1 FAIL (not in the exit status; --strict includes them)
```

The advisory rows are **out of the exit status** unless `--strict` is given: the
error that run exits with counts the six real checks alone (`1/6 validation
checks failed`), so a run that failed one check before these rows existed still
exits with exactly that. `--strict` counts every row instead (`2/8 validation
checks failed`), and the summary line then reads `(in the exit status:
--strict)`. A truth VCF carrying no census prints no `resistant` and no
`depth_fold` row, and its non-advisory rows and its exit status are exactly what
they were before either row existed -- an older spike's truth VCF is not a FAIL.
The advisory rows that do not read the census still print beside them, so such a
report does still carry an `Advisory:` line: measured on a one-DEL truth VCF with
both census INFO fields and both census header lines stripped, the report prints
`coverage_any_mapq` and `split_reads_each_end` and the line reads `Advisory: 2
checks, 2 PASS, 0 FAIL`. In `--json` every check object carries `"advisory":
true` or `"advisory": false` beside its `"pass"`, and the `summary` object
carries `counted_total`, `counted_pass`, `counted_fail` and `strict` beside
`total`, `pass` and `fail`: the three original keys are over every row printed,
while the `counted_*` trio is over the rows the exit status is computed from --
the non-advisory rows, or every row under `--strict` -- so a consumer can tell
which rows that status counted. `scripts/validate_pipeline.sh` is such a
consumer: its step-5 guard insists on at least one *non-advisory* check passing,
because `resistant` and `depth_fold` are read back from the truth VCF and pass
whatever the BAM holds.

#### What each check establishes, and what it does not

Each check above is evidence at one locus, and each is narrower than the claim
it is easy to read into it. What each one settles, and what it leaves open:

`coverage_ratio` establishes that the event's **mean** depth over its whole
span, divided by the mean of its two flanks (`--flank` bp on each side, default
5000), is within 0.30 of `1 - VAF` for a DEL or `1 + VAF` for a DUP. It does
not establish that the depth is right anywhere in particular: one average over
the whole span hides a local error the rest of the span cancels. Measured by
the review's reproduction script (`scripts/review_sv_model.py`) on a synthetic
probe, a het DUP of `chrT:10000-28000` over a donor whose interior section is
18.75x PASSed at an observed **1.32** against an expected **1.50** while that
interior read **81.09x** -- 4.32x, where a locally proportional CN2->CN3
predicts 28.13x (`CR2` in `docs/review/REVIEW.md` for the depth, `CR9` for the check passing). spike itself now measures that mismatch at simulation time, as `SIM_DEPTH_FOLD` in the truth VCF (3.88 on that probe; see [One depth for the whole event](#one-depth-for-the-whole-event)), and `spike validate` now reports it as the advisory `depth_fold` row beside this one. Nor does the check see a read its own
MAPQ filter rejects: on the same probes, with half the donor pairs at MAPQ 0, a
deletion requested at AF=1 kept **37.5x** of its reads inside the deletion and
still PASSed at an observed **0.00**, because `--min-mapq` (default 20) hides
exactly the reads that survived (`CR4`). spike itself now counts those reads
at simulation time, as `SIM_RESIST` in the truth VCF (0.500 on that probe; see
[Reads spike cannot edit](#reads-spike-cannot-edit)), and `spike validate` now
reports it as the advisory `resistant` row beside this one -- and recounts this
very ratio with no floor at all, as `coverage_any_mapq`.

`coverage_any_mapq` establishes that the same ratio, counted over **every
primary, non-duplicate, non-QC-fail read regardless of MAPQ**, is within 0.30 of
the same expectation. That is strictly more than `coverage_ratio` sees: on the
unit-test fixture above, a deletion window the floored row reads as **0.00**
reads **0.50** here, because the reads that survived the edit are the reads
`--min-mapq` throws away. What the row cannot do is tell **"spike could not
edit these reads"** from **"this locus is hard"**. A real locus of low
mappability reads thin in the donor pool and thick at any MAPQ whether or not
the library is uneven, and this row sees only the second half of that: it counts
the depth that is there, never why it is there. So a failing row is a reason to
read the locus's mappability and the run's own `SIM_RESIST` census, not by
itself evidence that the spike-in is wrong -- which is why the row is advisory
and out of the exit status unless `--strict` is given. It also inherits
`coverage_ratio`'s own limit: one mean over the whole span, hiding a local error
the rest of the span cancels.

`split_reads` establishes that at least **two** distinct read names, pooled
over the two breakpoints, sit within 500 bp of one breakpoint and carry an
`SA:Z` entry naming the partner's contig at a position within 500 bp of the
other. Pooled is the operative word: the two reads may both sit at the *same*
breakpoint, with nothing seen at the other, and the check still passes -- its
own `expected` string (`>=2 joining chr:pos`) reads as though each end must
contribute, and it does not. The advisory `split_reads_each_end` row beside it
is the same evidence with the minimum required at each end, and it is the row
that says so. That is the whole of it. The SA parser reads an entry's first two
fields -- contig and position -- and stops: the strand, CIGAR, MAPQ and NM
fields of the same entry are never looked at. So the check does not verify that
the two segments are on opposite strands, which for an INV is the one thing
separating its junctions from any other pair of splits; it does not locate the
breakpoint the CIGAR implies; and it reads no sequence at all. It is a count
with no denominator, so it says nothing about allele fraction -- two split
reads pass at any depth. In the other direction, a short deletion an aligner
represents as a `D` operation rather than as a supplementary alignment carries
no SA tag at all, so a spike-in whose sequence is right can still fail:
measured on the harness's chr20 window, a correctly planted 505 bp DEL emitted
188 ALT fragment pairs of which **0** carried an `SA:Z` tag, and `split_reads`
scored 0 -- a FAIL (`N1` in `docs/review/REVIEW.md`).

`split_reads_each_end` establishes the same thing at **each** breakpoint rather
than pooled over both: two distinct read names within 500 bp of the event's own
breakpoint carrying an `SA:Z` entry that names the partner, **and** two within
500 bp of the partner carrying one back -- for any event whose two breakpoints
are more than 500 bp apart. (For a shorter one the same two records answer both
halves; that blind zone is spelt out below.) That is what the pooled row's own
`expected` string reads as though it asked for, and it is what catches the
one-sided case above -- three reads at one end and nothing at the other reads
`3/0` here and FAILs, on the same input the pooled row PASSes. What it does
**not** add is any new kind of evidence: it is built from the same two queries
and reads the same two fields of an `SA:Z` entry, the contig and the position,
so it checks neither strand -- for an INV, still the one thing separating its
junctions from any other pair of splits -- nor the breakpoint the CIGAR implies,
nor a single base of sequence, nor any fraction. Both rows remain counts without
denominators.

And it does **not** tell one-sided evidence from two-sided for a short event.
The row inherits the pooled row's two 500 bp windows and the pooled row's 500 bp
tolerance on where an `SA:Z` entry may land, so when an event's own breakpoint
and its partner are within 500 bp of each other the two windows overlap *and*
each window's SA test reaches the other breakpoint: one read sitting only at the
left breakpoint, whose `SA:Z` names the right one, satisfies both halves at
once. The two counts become the same number and the row degenerates into the
pooled row, so for that whole size class it cannot distinguish one-sided
evidence from two-sided. Measured on this repo's own fixture, with three reads
at the left breakpoint naming the right one and **nothing at the right
breakpoint at all**: at `END - POS == 500` the row reads `3/3` and PASSes, and
at `END - POS == 501` the same layout reads `3/0` and FAILs. That boundary is
not academic -- `scripts/validate_pipeline.sh` sets `MIN_DEL_SIZE=500` and spike
itself puts no floor on SV size, so the pipeline's smallest admissible deletion
sits exactly on it. Reusing the pooled row's window was what this row was
specified to do; the blind zone is the cost, and it is recorded here rather than
worked around.

And because it is the pooled row's own evidence at a strictly stronger bar, it
can only fail **at least as often**: every event the pooled row fails, this row
fails too. That includes the case above where the aligner represented a correct
deletion with no `SA:Z` entry at all, which is why the row is advisory and out
of the exit status unless `--strict` is given. For an event whose breakpoints
are more than 500 bp apart it also adds the events whose evidence is one-sided;
at 500 bp or less it adds nothing, for the reason just given.

How much more often it fails on correct real data, measured: on 40 seeded 10 kb
deletions inside the HG002 T2T-Q100 SV benchmark on chr20, each spiked into the
35x HG002 BAM at the default AF and run through `align.sh`, `merge.sh` and
`spike validate`, the pooled `split_reads` row failed on **4 of 40** and
`split_reads_each_end` on **6 of 40** -- so the per-breakpoint row's own
contribution is **2 of 40**. Both of those two carried **7** joining reads at one
breakpoint and **1** at the other, which the pooled row passes at `observed 8`:
they are the one-sided evidence this row exists to name. The other four are the
pooled row's own failures, where no `SA:Z` entry *names the partner* at either end.
The entries exist: in those runs' `sim.bam` the junction reads are split, and 10,
3+6, 11 and 1+7 of them carry `SA:Z` at the two ends. But 35 of those 38 entries put
the supplementary piece on another chromosome, 32 of them at MAPQ 0. The sequence just
past the far breakpoint is repeated elsewhere, so the aligner cannot place the short
piece. A real deletion there would align the same way (RF6 in `docs/review/REVIEW.md`).
That is why, for a DEL, `split_reads` is advisory since RF14 and `del_planted`
decides.

`del_planted` establishes that at least **one** of the reads spike made for a
deletion is in the BAM within 500 bp of START or END (the windows `split_reads`
reads), carrying the deletion's join. It is the DEL row that decides (RF14).
- **Whose reads.** spike's reads are known by name: truth record `sim_del_N` goes
  with reads named `SPIKE_evNNNN_hap_...` (see [Read names](#read-names)).
- **Carrying.** A read carries the deletion when its bases hold the event
  haplotype's 31 bases across the join, `ref[START-15, START) + ref[END, END+16)`,
  in either orientation, with at most 2 bases differing. Those are the sample's
  own SNPs, which spike writes onto the event copy, and sequencing errors.
- **Which records.** Every record counts but a secondary or supplementary one:
  MAPQ, duplicate, QC-fail and unmapped flags are the aligner's verdict, not the
  question. No read of the sample's own can count, so one read is enough.

What it does **not** establish:
- **How the aligner split the reads.** That is what `split_reads` and the realism
  probe are about.
- **That the bases were removed, where the join is already reference
  sequence.** A deletion of repeat units, or one between two copies of a repeat,
  has a join that unedited reads spell too. There the row shows spike's reads are
  present, and `coverage_ratio`, still counted, is what measures the removal. At
  RF14's fresh real HG002 sites this was 9 of 12 and 10 of 12.
- **A truth VCF from another tool.** Without a `sim_del_N` ID the row is a failed
  not-evaluable row.

Measured on RF14's fresh sites in the 35x HG002 BAM:
- **At VAF 0.5**, every correct deletion spike accepts by default carries: 33 of 33
  (8 to 35 reads). These were 50 bp to 10 kb, 9 of them real HG002 deletions. The
  other 3 real ones are events spike refuses by default (`SIM_RESIST` above
  0.5); forced through with `--allow-resistant`, they carry too.
- **At VAF 0.1**, 2 of 30 correct 1 kb deletions had no read of spike's across
  the join, by chance (95% interval 0.018-0.213). The row fails those. The old
  `split_reads` failed 14 of the same 30.
- **The same truth, made wrong** (END + 50, START - 50, moved 1 kb, the ID
  renumbered, or the unspiked BAM): 0 carriers in 330 tries.
- **200 empty sites:** 0 pass.

`ins_planted` establishes that at least **one** of the reads spike made for the
event is in the BAM within 150 bp of POS, carrying its inserted bases across a
junction. It is the INS row that decides (RF13). spike's reads are known by
name: truth record `sim_ins_N` goes with reads named `SPIKE_evNNNN_hap_...` (see [Read names](#read-names)). A read
carries the insertion when its bases hold one of two 31-base junction probes
from the event's own haplotype (15 reference bases, then 16 past the junction,
at each end of the insertion), with every **inserted** base matching exactly and
at most 2 **reference** bases differing (the sample's own SNPs, which spike
writes onto the event copy, and sequencing errors). Every such record counts but
a secondary or supplementary one: its MAPQ and its duplicate, QC-fail and
unmapped flags are the aligner's verdict, not the question. No read of the
sample's own can count, so there is no background to rise above, and one read
is enough. It does **not** establish how the aligner wrote the insertion down --
that is what `ins_reads` and the realism probe are about -- and it cannot judge
a truth VCF from another tool: without a `sim_ins_N` ID, or without inserted
bases in the ALT (a symbolic `<INS>`, which spike wrote before CR7), the row is a
failed not-evaluable row.

`ins_reads` establishes that at least **two** reads leave the reference within
100 bp of POS, by an `I` operation of at least `min(SVLEN, 50)` bases or by a
soft clip of that length once the threshold has reached 50. It works from the
CIGAR alone -- the record's own sequence is never read -- so it does not
establish that the inserted bases are the bases the truth VCF names. An
insertion of roughly the right length in roughly the right place passes,
whatever it spells. The advisory `ins_sequence` row beside it is the row that
reads the bases. Like `split_reads` it is a count, not a fraction.

`ins_sequence` establishes that at least **two** distinct read names within
150 bp of POS carry, somewhere in their own bases, one of the two probe k-mers
taken from the truth record's ALT -- its first and last `k = min(inserted
length, 31)` bases -- in one orientation or the other. It is the first check
here to read a base of inserted sequence rather than a CIGAR, and what it
catches is the case `ins_reads` cannot see at all: an insertion of the right
length in the right place spelling something else entirely, which fails this row
and passes that one on the same input.

What it does **not** establish is anything about *where* those bases are. It
asks only that they are **present**: not that they are at the right offset, not
that the insertion is in the orientation the truth record names, not that they
occur the right number of times, and not at what allele fraction -- a read
carrying the k-mer anywhere in its 151 bases counts, and like the two rows above
it is a count without a denominator. Nor does it read the whole insertion. Two
limits stack there, and only the second is the row's own. First, only sequence
near the insertion's two **ends** is reachable at all: a read anchored beside
the insertion point reads into it by soft-clipping, and by at most about its own
length, so of a 2000 bp insertion roughly the first and last 150 bases can sit
inside a read and the middle cannot -- that insertion's middle k-mer sits 1000
bases in, far past the reach of a 151 bp read, which is why the probes are taken
from the ends and not from the middle. Second, of what is reachable the row
probes `k = 31` bases at each end and no more, not the ~150 a read could show
it. So for a **2000 bp insertion it verifies 62 of 2000 bases** and says nothing
whatever about the other 1938. And it believes the file it reads: the inserted
bases come from the truth VCF's ALT, so a truth VCF edited between the run that
wrote it and the validation that reads it is taken at its word -- which is also
what makes the row's own separation measurable, by validating one run's BAM
against another run's truth VCF.

Two inputs it declines rather than grades, both reported as failed rows and both
harmless in the exit status because the row is advisory: an insertion of fewer
than 12 bases, whose k-mer would be too short to be specific inside a read, and
an ALT one of whose probe k-mers the reference already holds within 1 kb of POS,
where an unedited read would match. A symbolic `<INS>` ALT recorded no sequence
and gets no row at all.

How often the row fails on correct real data, measured: on 24 insertions -- four
at each of 50, 100, 250, 500, 1000 and 2000 bp, at the first 24 starts of the
same seeded list, each on its own +/-100 kb slice of the 35x HG002 BAM through
`align.sh`, `merge.sh` and `spike validate` -- it failed on **0 of 24**, and on
**0 of 4** in every one of the six length classes. Supporting reads: min 10,
median 23, max 34. The count does not fall off with length, because the probes
are the insertion's first and last `k` bases and a read reaches those from either
side whatever sits between them.

`allele_freq` is the one per-event check that measures a fraction rather than
counting evidence; the table below says what it can and cannot read off a truth
record.

`resistant` and `depth_fold` establish nothing at all about the BAM. Both
numbers are **read back from the truth VCF**, not recomputed from the reads:
`resistant` is the `SIM_RESIST` share spike counted while it was editing the
donor, `depth_fold` the `SIM_DEPTH_FOLD` fold it measured while it was scaling
the event's fragments. So the two rows report what spike measured at simulation
time, and they inherit its blind spots -- under `clean`, `SIM_DEPTH_FOLD`'s
donor pool holds only reads at `--min-mapq` or above, so a bin thin only in
mappable reads still counts as a fold. They also believe the file they read: a truth VCF edited
between the run that wrote it and the validation that reads it is taken at its
word, because nothing here re-derives either number. Nor do they read which
edit model wrote it: a truth VCF from `--edit-model origin` is graded against
the same two thresholds, which were set for `clean`. What a failing row does
establish is that spike's own census of the run went past the threshold spike
warns at, which is a reason to read that run's log and README rather than this
BAM.

The three `[global]` checks are **fixed library heuristics, not comparisons
against the donor**. `insert_size` passes when the sampled mean is 50-1000 bp
and its SD is 5-300 bp, `dup_rate` when the duplicate rate is under 50%, and
`mean_mapq` when the mean MAPQ is above 20. `spike validate` never reads the
donor BAM's own measured fragment distribution, so none of the three asks
whether the simulated reads resemble the library they came from. **A nonzero
exit therefore need not mean an event is wrong:** any failed check fails the
run, and a global one can fail on a property of the file that has nothing to do
with the spike-in. On the review's synthetic probes every fragment is exactly
400 bp, so `insert_size` observed `400+/-0` and failed on the SD floor of 5,
and nothing had marked duplicates, so `dup_rate` reported `no dup flags` and
failed; all four probes exited **1** with both of those rows failing. Read the
failing rows, not the exit code alone.

**A check that cannot measure its answer is a failed check, at event level as
well as globally.** `allele_freq` reads a truth record's REF/ALT pair, picks
the counting rule that fits its shape, and fails outright when none does:

| Truth record | Counted as | Observed |
| --- | --- | --- |
| `A` > `T` (a substitution) | the alt base against the pileup depth at POS | a fraction |
| `ACG` > `A` (a small deletion) | reads whose bases over the site are nearer the reference with the deletion made, against reads whose bases are nearer the reference | a fraction |
| `A` > `ACCGG` (a small insertion) | reads whose bases over the site are nearer the reference with the insertion made, against reads whose bases are nearer the reference | a fraction |
| `AT` > `GC` (an MNV) | reads whose bases are the *whole* alt run, against reads whose bases are the whole ref run | a fraction |
| `AC` > `GTT` (a complex allele) | nothing -- no single operation or allele run to count | `N/A (complex allele)` FAIL |
| depth below 5 | nothing | `low depth (n)` FAIL |
| an allele that is not A/C/G/T | nothing | `unknown alt base` / `unknown allele base` FAIL |

An indel is read from each read's own bases rather than the pileup, because
an indel is not a column in one. `validate` writes out short sequences over the
site, from 10 bp before the indel's repeat region to 10 bp past it (below): the
reference with the truth record's edit made, and without it. It takes the bases
a read shows between those two ends, inserted bases included, and counts the
single-base edits that turn them into each sequence (the Levenshtein distance).
A read nearer a sequence with the edit carries the allele, a read nearer one
without it spans it, and a read equally near both kinds does not vote. So it
does not matter where the aligner put the gap: along a repeat it may write one
deletion anywhere, and the read's bases are the same wherever it goes. A gap of
the same size somewhere else is not taken for the truth's either: outside a
repeat it removes other bases, and three or more bases off it is nearer the
reference.

**Other truth records near the indel go into those sequences too.** A truth
set often writes one local change as several records a few bases apart, and
the reads carry the whole change: `ACG` > `A` next to `T` > `TCG` is `CGT`
become `TCG`, which an aligner writes as three mismatches. Against the
reference plus one of the two records, such a read is two edits from each, so
it could not vote. So every other record in the truth VCF whose REF lies
inside the stretch is a candidate. The sequences are the reference with every
combination of the site and the candidates applied (never two that overlap),
and the read's own bases pick the nearest. No genotype is needed, so an
unphased truth set works. With more than 10 candidates the site is compared
with its own record alone. A record is only seen if it is in the truth VCF
`validate` is given. An MNV's bases are read **jointly**, one read at a time: a
fraction per base would answer a different question at each offset, and a read
carrying only one of the two substitutions is not this variant.

Every count in the table is of **read pairs**, not reads. Where the two mates
of a pair overlap the variant they read the same DNA molecule, so they give
one vote, and a pair whose mates disagree gives none. "Depth below 5" is five
pairs. On HG002 at 35x (fragments 418 ± 178 bp, 151 bp reads) mates overlap
often enough that counting them twice changed the fraction at 75% of 1,959
real het sites.

Unless the truth VCF lists it, the rule cannot tell a *different* insertion
of the same size at the same place from the truth's: it is as many edits from
the edited sequence as from the reference, or fewer, so it counts as carrying
or not at all, never as spanning. When the other allele is a truth record of
its own, its reads are nearest it and count as spanning.

**Only reads that reach well past an indel vote on it, either way.** Near
its end a read's indel is written as mismatches or a clip rather than a gap,
and a read that stops inside the repeat the indel sits in cannot show an
extra or missing unit at all. Either way it aligns as the reference, whatever
it carries. So `validate` first finds the indel's repeat region: the deleted
or inserted unit, extended along the reference for as long as it repeats. A
read then votes only if it has a base aligned to the reference at each end of
the stretch above -- the base on each side of that region, and 10 more beyond
it. A read that stops short, or has a gap or a clip there, enters neither
count. The test is the same for a carrier and for a reference read.

Measured on HG002 35x, graded against GIAB het indels at 0.5:

| | before | after |
| --- | --- | --- |
| chr20 (6,663 indels): out of range | 9.3% | 3.7% |
| chr20: mean fraction | 0.41 | 0.48 |
| chr21 + chr22 (8,008 indels, not used to choose 10): out of range | 9.3% | 4.1% |
| fragments kept | | 81% |

Comparing each read with every combination of the truth records near an
indel (above) was then measured the same way, on chromosomes not used to
design it (chr17-chr19, 25,903 het indels, with GIAB's other records within
500 bp in the truth file):

| out of range | before | after |
| --- | --- | --- |
| indels with another GIAB record within 25 bp (4,482) | 17.6% | 8.7% |
| isolated indels (21,421) | 1.4% | 1.4% |
| all | 4.2% | 2.6% |

SNVs fail at about 0.6-0.8% on the same BAM.

The three rows that measure nothing still push a result row, so the event
counts as covered, and the row says out loud that nothing was measured. Only a
**complex** allele -- one that changes length *and* rewrites the anchor base,
such as `AC` > `GTT` or `A` > `CG` -- is left unmeasured.

An **unrecognised `SVTYPE`** (`CNV`, `DEL:ME`, …) keeps its own type. It used
to fall through to the small-variant arm, where a symbolic ALT such as `<CNV>`
is longer than one base and took the indel exit above — one silent PASS per
unknown type. It is now reported as `event_checked FAIL`, with a `WARN` naming
the type.

A truth record whose **`END` is at or before its own `POS`** is refused when
the truth VCF is read, naming the record; the run exits non-zero without
grading anything. The event region is empty, `count_depth_in_region` answers a
zero-length region with `0.0`, and a DEL at `SIM_VAF=0.9` therefore PASSed
`coverage_ratio` (expected 0.10, observed 0.00) over a region no query had
read.

`ins_reads` counts reads whose alignment **leaves the reference at POS**: an
insertion has no second breakpoint and no reference span, so neither
`coverage_ratio` nor `split_reads` can see it, but an aligner still has to put
the inserted bases somewhere -- an `I` CIGAR operation when the insertion fits
inside a read that anchors on both sides, a soft clip at the insertion point
when it does not. Either counts, if it is at least `min(SVLEN, 50)` bases long
and its reference boundary is within 100 bp of POS, and the check passes at
two such reads (the same threshold `split_reads` uses). The two need not be at
the same point, and at 35x two such clips turn up by chance more often than
that suggests (see below). An INS record with no
usable `SVLEN` has no length to look for and is a failed check. Nothing in this
check reads a base; the advisory `ins_sequence` row beside it does, by looking
for the truth record's own inserted bases in the reads (see
[Validating the spike-in](#validating-the-spike-in)).

A **soft clip only counts once the insertion is 50 bp or longer** -- the same
cap `min(SVLEN, 50)` applies. A read that anchors both sides of a short
insertion writes an `I` operation, so below the cap a clip is background.
Measured on the merged HG002 chr20 slice, a truth record of
`SVTYPE=INS;SVLEN=3` at five positions where nothing was planted found **1, 3,
0, 0 and 1** clipped reads, and 3 is over the two-read pass mark -- one of the
five PASSed on an insertion that was never there. Counting only `I` operations
below the cap, the same five positions give **0, 0, 0, 0, 0**, while genuinely
planted 3 bp and 12 bp insertions still find 12 reads each and PASS.

**What `ins_reads` gets wrong, measured (RF11 in `docs/review/REVIEW.md`).** Both errors
come from counting CIGAR marks near POS rather than reading bases, and they pull
in opposite directions, so no clip threshold fixes both:

- **Correct 40-49 bp insertions can FAIL.** Under bwa-mem2's default penalties
  an `I` costs `6 + SVLEN` and a clip costs 5, so a read with fewer than about
  SVLEN bases on one side of the insertion scores better clipped. Near 50 bp
  that is almost every read, and below 50 a clip does not count. (At a 40 bp
  insertion, all 11 reads carrying it were clipped.) On the 35x HG002 BAM, correct spike-ins at 8 sites
  failed `ins_reads` 0 of 8 times at 20 and 30 bp, 1 of 8 at 40 bp, 5 of 8 at
  45 bp and **8 of 8 at 49 bp**. `ins_sequence` passed all 40.
- **An insertion that is not there can PASS at 50 bp and up.** With nothing
  planted, 12 of 200 random chr20 sites and 46 of 200 sites in simple repeats
  already have two reads clipped by 50 or more bases within 100 bp (measured
  with `spike validate` itself). At the random sites most (9 of 12) are
  scattered clips, not a breakpoint.
- Counting clips from 20, 30 or 40 bp would fix the first and worsen the
  second: 22, 18 or 15 random sites and 64, 58 or 50 repeat sites.

So a FAIL on `ins_reads` for a 40-49 bp insertion usually means bwa-mem2 clipped
the reads, not that the insertion is missing: look at `ins_sequence`. And a PASS
at 50 bp or more is weaker evidence than two reads suggests.

**Real variants fail these checks too** ("Realism probe" in `docs/review/REVIEW.md`). Run on
HG002's own BAM against GIAB's truth:
- `split_reads` fails 12 of the pipeline's 20 real deletions.
- `ins_reads` fails 47 of 112 real het insertions of 20-39 bp. spike's own
  insertions of that size, at random unique sites, fail 0 of 16.

So a FAIL on either row does not show that a spike-in is wrong. And spike's short
insertions are easier to align than real ones (RF15 in
`docs/review/CLINICAL_SV_NEW_FINDINGS.md`).

**So `ins_reads` is advisory, and `ins_planted` decides an insertion** (RF13 in
`docs/review/REVIEW.md`). On the 35x HG002 BAM:
- **54 correct insertions.** These were 1, 2, 4, 15, 45 and 300 bp random ones at
  6 sites, 45 bp at VAF 0.1, and 12 of HG002's own insertions at their own
  positions with their own bases. `ins_planted` found spike's reads carrying
  every one: 2 to 37 reads.
- **Master's exit status** was 0 on 42 of the 54 and is 0 on all 54 now. The 12
  that moved are the 45 bp ones `ins_reads` failed.
- **The same truth, made wrong.** Every inserted base swapped, POS moved by
  1 kb, the ID renumbered, or the unspiked BAM: 0 carriers in 216 tries.
- **200 empty sites.** 0 pass `ins_planted`; `ins_reads` still passes 10, now
  only as information.

**And for a DEL, `split_reads` is advisory and `del_planted` decides** (RF14 in
`docs/review/REVIEW.md`; the measured rates are in the `del_planted` paragraph above). Master's
and the new `validate` were run on the same merged BAMs, 178 saved runs:
- **No run master passed now fails.**
- **36 runs flipped from fail to pass** (109 → 145 exit 0), including RF6's four
  chr20 deletions.
- **The pipeline's 20 real deletions on NA18488:** 5 → 11 exit 0. The 9 still
  failing fail `coverage_ratio` alone, at sites where the background already
  carries the deletion (RF12).

Measured on a DEL+INS run on the HG002 chr20 slice, aligned with `align.sh`
and merged with `merge.sh`: **14** reads carry the planted 300 bp insertion at
chr20:39000000. Before this check existed, spike's own round trip could not
succeed for insertions -- the same recipe scored `5/6 PASS` and exited **1** on
the `event_checked` row -- and this run scores `13/13 PASS` and exits 0. The
total rose because seven of those thirteen rows are the advisory rows the two
events now carry beside their own checks: the same run prints `Advisory: 7
checks, 7 PASS, 0 FAIL`. Its control is from an **earlier** run of the same
recipe, which read **21** reads at chr20:39000000 against **0, 0, 0, 0 and 1**
at five control positions in the same BAM where nothing was planted; the
re-measured run above was not re-run against those five positions, so the 14 and
the `13/13` are one run, and the 21 and those five controls are the other.

The same round trip works for small indels and MNVs. A run of
`snp:chr20:39000000:TGG:T` (a 2 bp deletion), `snp:chr20:39100000:T:TCCGG` (a
4 bp insertion) and `snp:chr20:39200000:AT:GC` (an MNV) on the same slice,
aligned and merged the same way, scored **3/6 PASS and exited 1** with all
three rows reading `N/A (indel or MNV)` -- a verdict reached before the BAM was
opened -- and now scores **12/12 PASS, exit 0** at **0.50, 0.53 and 0.42**
against a `SIM_VAF` of 0.50, the total again rising because six of the twelve
rows are advisory (`Advisory: 6 checks, 6 PASS, 0 FAIL`). Those three fractions
and that total are one run, measured together. The controls beside them are from
an **earlier** run of the same three records, which read 0.40, 0.41 and 0.46
off **17, 19 and 13** reads carrying the variant at the three planted sites,
against **0, 0, 0, 0 and 0** for each of them at five positions where nothing
was planted (38600000, 38900000, 39500000, 39750000, 40100000), where all
fifteen checks read 0.00 and FAIL, and where the same three truth records
against the **unspiked** BAM read 0.00, 0.00 and 0.00 and all FAIL. The
re-measured run above was not re-run against those five positions or against the
unspiked BAM, so do not read its fractions and those controls as one
measurement.

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

### Run time

`--threads` (default 4) sets spike's own threads as well as the ones `align.sh` and `merge.sh` give the aligner and samtools. spike splits only the work that draws no random numbers (reading regions of the BAM, the quality profile, `origin`'s look-alike reads and removal chances, writing R1 and R2), so a given `--seed` gives the same files at any thread count. Tiling reads, phasing and drawing which reads to remove stay on one thread.

Measured on the 35x HG002 BAM with `--seed 1 --allow-resistant`, one run at a time on a 72-core machine: wall time and peak memory from `/usr/bin/time`, the mean of three runs at 1 and 8 threads and one run at 4 and 16.

| Event | Model | 1 thread | 4 (default) | 8 | 16 | Peak memory at 1 / 4 / 8 |
| --- | --- | --- | --- | --- | --- | --- |
| `del:chr20:7119236-7120236` (1 kb) | clean | 0.55 s | 0.57 s | 0.56 s | 0.56 s | 194 / 194 / 194 MB |
| | origin | 0.63 s | 0.75 s | 0.74 s | 0.82 s | 188 / 222 / 338 MB |
| `del:chr20:14550000-17550000` (3 Mb) | clean | 24.0 s | 12.5 s | 10.7 s | 9.9 s | 1231 / 1380 / 1711 MB |
| | origin | 49.2 s | 23.7 s | 18.5 s | 16.3 s | 1819 / 1908 / 2168 MB |
| `dup:chr20:14550000-17550000` (3 Mb) | clean | 44.3 s | 29.8 s | 28.2 s | 26.5 s | 1493 / 1787 / 2028 MB |
| | origin | 66.6 s | 39.4 s | 34.4 s | 31.4 s | 2126 / 2215 / 2500 MB |

A small event under `origin` takes slightly longer on more than one thread than on one.

The `clean` rows were measured again on 2026-10-04, after `clean` began reading each event's windows a second time to find the duplicates of the pairs it removes. The binary from just before that change, run alternately with the new one, took 19.2 s and 38.5 s at one thread on the 3 Mb deletion and duplication, and 10.1 s and 26.8 s at four. So the second read costs 5-6 s at one thread on a 3 Mb event, 2.5-3 s at four, and under 0.1 s on the 1 kb one. Peak memory moved by under 5% either way. The `origin` rows are from the earlier measurement; `origin` does not take this step.

## Output files

| File | Description |
|------|-------------|
| `R1.fq.gz` | spike's reads around the events only, read 1 (gzipped FASTQ): the originals it kept and its own new reads, **not the whole sample** |
| `R2.fq.gz` | The same, read 2 |
| `truth.vcf` | VCF with simulated variant records and AF annotations |
| `events.bed` | Extraction regions (event ± flank) used to build the spike-in |
| `replaced_reads.txt` | Names of the originals spike extracted, including pairs dropped for unusable quality, the duplicates of every pair it removed, and under `--edit-model origin` every read `origin` removed; `merge.sh` removes exactly these |
| `align.sh` | Aligns R1/R2 → `sim.bam` (event regions only) |
| `merge.sh` | Merges `sim.bam` into the original BAM → `merged.bam` (full genome) |
| `README.md` | Run log: command, events table (a **Requested VAF** and a **Simulated VAF** column per event -- the same pair `truth.vcf` records as `SIM_REQ_VAF` and `SIM_VAF` -- then kept, chimeric and suppressed reads, pairs dropped for unusable quality, **Resistant reads** and **Depth fold**), the pairs dropped in total, read counts, next-step instructions |
| `sim.bam` | Aligned BAM covering event regions (produced by `align.sh`) |
| `merged.bam` | Original BAM with spiked reads substituted (produced by `merge.sh`) |
| `fastq_removed_reads.txt` | The originals spike removed and does not write back, their duplicates among them: `replaced_reads.txt` minus the pairs in R1/R2. `fastq.sh` removes exactly these from the raw FASTQ |
| `fastq.sh` | Builds the full spiked FASTQ pair from the sample's raw FASTQ pair (see [Full FASTQ](#full-fastq)) |
| `NAME_R1.fastq.gz`, `NAME_R2.fastq.gz` | With `--into-fastq`: **the whole sample, spiked**, the pair for a pipeline that starts from FASTQ. `NAME` is `--fastq-prefix`, default `spiked` (see [Full FASTQ](#full-fastq)) |

### Full FASTQ

A pipeline that starts from raw FASTQ, such as nf-core/raredisease, needs the whole sample as a FASTQ pair. `--into-fastq` spikes the events into the sample's own raw FASTQ:

```bash
spike --bam sample.bam --reference ref.fa --vcf variants.vcf -o output/S1 \
      --into-fastq RAW_R1.fastq.gz RAW_R2.fastq.gz --fastq-prefix S1
```

This writes `output/S1/S1_R1.fastq.gz` and `output/S1/S1_R2.fastq.gz`: the whole sample, spiked. These are the pair to give the pipeline.
- **The prefix.** Without `--fastq-prefix` the pair is `spiked_R1.fastq.gz` and `spiked_R2.fastq.gz`. With many samples, a prefix per sample keeps the pairs apart. It is a plain name, with no `/`, because `-o` picks the folder.
- **Not `R1.fq.gz` and `R2.fq.gz`.** Those hold only spike's reads around the events (the originals it kept and its own new ones), not the whole sample.
- **Checks.** Before it starts, spike checks that both raw files exist and are two different files, that neither is a file the run writes (such as `S1_R1.fastq.gz` or `R1.fq.gz` in `-o`), and that the prefix is a plain name.
- **At the end** it runs the output's `fastq.sh`, and fails if `fastq.sh` refuses the raw files.

`fastq.sh` can also be run by hand on a finished run:

```bash
bash output/S1/fastq.sh RAW_R1.fastq.gz RAW_R2.fastq.gz S1_R1.fastq.gz S1_R2.fastq.gz [THREADS]
```

What `fastq.sh` does:

- **It keeps every raw read exactly as it was** (bases, qualities, order and header), except the originals listed in `fastq_removed_reads.txt`. Those are the pairs spike removed and did not write back, and their duplicates. The originals spike kept stay as the raw reads they are, not as their copies from the BAM. The BAM's copies have been through the pipeline's own trimming and correction (fastp, in raredisease), and the raw reads the pipeline later filtered out are still there too.
- **It then adds spike's own reads** (`SPIKE_...`) at the end of each file. Their headers take the raw file's style, copied from its first record: for example ` 1:N:0:ACGTACGT+TGCATGCA` after the name, or `/1`.
- **A read is matched by name:** its header's first word, without `@` and a trailing `/1` or `/2`. That is the name the aligner gave it in the BAM spike was run on.
- **It reads both mates together,** record by record. They must hold the same names in the same order, because an aligner pairs R1 and R2 by position.
- **It stops, and changes nothing,** when:
  - an output is the same file as an input (a raw file, or spike's own `R1.fq.gz`, `R2.fq.gz` or `fastq_removed_reads.txt`), or both outputs are one file, or both raw files are one file. Writing an output over an input would destroy it;
  - the mates differ in name at some record, or one has more records;
  - a record is broken: a header without `@`, a third line without `+`, a sequence and quality of different lengths, or a file that ends inside a record;
  - a listed original is missing from the raw pair, or is in it more than once. Missing means the raw FASTQ is not the full FASTQ of that BAM's run. A sample sequenced on several lanes is concatenated first (`cat` joins gzip files), in the same lane order for R1 and R2.
- **Outputs appear only when finished.** Each is written to a hidden temporary file beside it and moved into place once both mates are complete and checked. A run that stops leaves files already at those names as they were.
- It needs no `align.sh` or `sim.bam`. It uses `pigz` when installed, and `gzip` otherwise.

The checks are in `docs/superpowers/plans/2026-10-04-full-fastq.md` and `docs/superpowers/plans/2026-10-04-fastq-safety.md`.

### Read names

Every read spike makes is named `SPIKE_<event>_<kind>_<counter>`: for example `SPIKE_ev0001_hap_000123` (tiled reads), or `SPIKE_ev0001_dup_depth_000123` (depth copies under `--dup-model junction`). `samtools view merged.bam | grep SPIKE_` finds every one of them. spike never edits a real read in place, so every other read in `merged.bam` is an original.

Tools also read a flowcell position out of a read's name. Picard MarkDuplicates, which raredisease runs, splits the name on `:`. When there are 5 or 7 fields, it takes the last three as tile, x and y, and uses them to find optical duplicates. Any other name gets no position, and Picard warns once that its `READ_NAME_REGEX` did not match.

Real names differ by machine: for example `A00744:46:HV3C3DSXX:2:1221:8775:9361` (NovaSeq 6000) and `D00360:96:H2YLYBCXX:1:2105:5916:51581` (HiSeq 2500). So spike learns the shape from the input's first 50,000 primary records and names its own reads the same way. The shape is the one most of those names have:

- **7 fields** (`MACHINE:RUN:FLOWCELL:LANE:TILE:X:Y`) or **5 fields** (`MACHINE:LANE:TILE:X:Y`), with the last three numbers. spike's reads keep the input's most common `RUN:FLOWCELL:LANE` (or `LANE`), and get a tile seen with it and an x and y within the range seen. From the first input above, that is `SPIKE_ev0001_hap_000123:46:HV3C3DSXX:2:<tile>:<x>:<y>`. The machine field is the only one changed, and Picard does not read it.
- **Anything else** (or a tie): spike's reads are named `SPIKE_ev0001_hap_000123`, with no position. Picard cannot read a position from the real reads either.

The tile, x and y come from a hash of the read's name, not from the run's random stream, so the name changes nothing else about a run. The log line `Read names: ...` says which shape was learned. The checks are in `docs/superpowers/plans/2026-10-04-read-names.md`.

### Truth VCF

The truth VCF contains one record per simulated event with:
- Standard VCF fields (CHROM, POS, REF, ALT)
- `SVTYPE` and `END` / `SVLEN` for structural variants
- A sequence-resolved `ALT` for an insertion (the anchor base at `POS` plus the inserted bases, not a symbolic `<INS>`), so the file grows by roughly one byte per inserted base
- `SIM_VAF` in the INFO field with the allele fraction that was **simulated** -- the fraction of the depth the fragments spike planted actually make up. On an **additive** event (a fusion, or a DUP under `--dup-model junction`) that is the *junction* evidence: the fraction the fragments across the breakpoint make up. A junction DUP also plants interior depth copies, and those are scaled by `SIM_REQ_VAF`, not by the capped fraction, so a capped one's interior dosage is above its `SIM_VAF`: at `af=0.99` the junction gets the 0.950 recorded while every interior copy is drawn at the uncapped 0.99. The default `--dup-model full` tiles the whole tandem haplotype and has no such split
- `SIM_REQ_VAF` with the fraction that was **requested** (`af=`, or `--allele-fraction`). The two differ exactly where a mechanism moved the count off the request: the additive 0.95 cap puts `SIM_VAF` below `SIM_REQ_VAF`, the two-fragment floor puts it above. Rounding the count to a whole fragment does not: `SIM_VAF` is the request unless one of those two applied
- `SIM_RESIST` with the share of the reads over the event that spike **could not edit**: primary, mapped, non-duplicate, non-QC-fail reads at any MAPQ whose pair is not in the event's donor pool (below `--min-mapq`, not a proper pair, a mate unmapped or failing a filter). They stay in `merged.bam` as they were, so the event realised is weaker than requested by about this share. The reads counted are those over the span a DEL, DUP or INV changes, the two bases around an insertion point, a small variant's REF, and the two bases around each fusion cut. `spike validate` reports it as the advisory `resistant` row. Above 0.5 spike refuses the event unless given `--allow-resistant` (RF8), so in a truth VCF written since then a value above 0.5 means that flag was used. Under `--edit-model origin` a read in a pair `origin` may remove counts as editable too. See [Reads spike cannot edit](#reads-spike-cannot-edit)
- `SIM_DEPTH_FOLD` with the largest fold between the donor's depth in any ~1 kb bin the event's synthetic fragments are drawn from and the one depth they are all scaled by. Where the two differ, the event's depth there is off by about that fold. Under `--edit-model origin` both depths are origin depths. `spike validate` reports it as the advisory `depth_fold` row. See [One depth for the whole event](#one-depth-for-the-whole-event)
- Under `--edit-model origin` only, a header line `##spike_edit_model=origin (experimental).` followed by a note that the `SIM_RESIST` and `SIM_DEPTH_FOLD` thresholds were set for `clean` and are not yet checked for `origin`
- `SIM_GENE` with the associated gene name
- BND records for fusions (with `]`/`[` notation reflecting orientation)

The run's own `README.md` prints the same two fractions per event, as the
**Requested VAF** and **Simulated VAF** columns of its events table, so the two
files never disagree about what the reads carry.

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
          Event specification(s). Can be repeated. Formats: --event "del:chr20:30000000-30005000"           (coordinate-based) --event "del:GENE:exon4-exon8"                  (gene-based, requires --exon-bed) --event "dup:GENE:exon4-exon8"                  (gene-based duplication) --event "inv:GENE:exon4-exon8"                  (gene-based inversion) --event "fusion:GENEA:exon14:GENEB:exon2"       (fusion: adds one junction, not a balanced translocation; requires --exon-bed) --event "dup:chr20:30000000-30005000" --event "inv:chr20:30000000-30005000" --event "ins:chr20:30000000:500"                (random insertion sequence) --event "ins:chr20:30000000:ACGTACGT"           (explicit insertion sequence) --event "snp:chr20:30000000:A:T"                (SNP/small variant, POS is 1-based) --event "snp:chr20:30000000:A>T"                (alternate syntax, POS is 1-based) --event "snp:chr20:30000000:ACG:A"              (small deletion, POS is 1-based) --event "snp:chr20:30000000:A:ACGT"             (small insertion, POS is 1-based) Per-event AF (appended with ;): --event "del:GENE:exon4-exon8;af=0.15" --event "fusion:GENEA:exon14:GENEB:exon2;af=het" af=<number>: exact AF, af=het: 0.5, af=hom: 1.0

      --vcf <VCF>
          Input VCF file with variant records. Supports DEL, INS, DUP, INV, BND, and standard SNP/indel records (no SVTYPE, explicit REF/ALT alleles). Can be combined with --event. At least one of --event or --vcf required

      --vcf-info-af
          Read INFO/AF from --vcf records as the allele fraction to simulate.
          
          Off by default: in a population VCF (gnomAD, 1000G) AF is the allele frequency in the population, not the fraction of this sample's reads that should carry the allele, so using it silently produces a truth set at the wrong VAF. Without this flag only SIM_VAF and VAF are read, and a record that carries AF but has no usable SIM_VAF or VAF falls back to --allele-fraction, with a count of how many did on stderr.

      --exon-bed <EXON_BED>
          Exon BED file. Required when using gene-based --event specs (e.g. "del:GENE:exon4-exon8")

      --allele-fraction <ALLELE_FRACTION>
          Target allele fraction, in (0.0, 1.0] -- above 0 and at most 1.
          
          0 is rejected, not treated as "plant nothing": an event asked for at AF 0 would still be written to the truth VCF. NaN is rejected for the same reason.
          
          [default: 0.5]

  -o, --output <OUTPUT>
          Output directory
          
          [default: /tmp/spike]

      --seed <SEED>
          Random seed for reproducibility
          
          [default: 42]

  -t, --threads <THREADS>
          Threads for spike and for the align.sh and merge.sh scripts it writes.
          
          spike splits the work that draws no random numbers across this many threads, so a given --seed gives the same output at any count. The aligner (bwa-mem2 -t, minimap2 -t, bowtie2 -p) and samtools (sort -@, merge -@) get this many too.
          
          [default: 4]

      --region <REGION>
          Extra read extraction region (e.g. "chr19:11080000-11140000"), on top of event ± flank -- it does not replace the event window. Use this to ensure the output BAM covers the full gene/region of interest. For an event on the same chromosome, the region and the event ± flank are merged into one query when they overlap or touch, and kept as two queries when they do not, so a distant fusion partner costs one extra event-sized window rather than every read in between. A region on another chromosome than the event is ignored for that event

      --flank <FLANK>
          Flanking region (bp) to include around events. Always defines the event window; when --region is also set, the region adds another window alongside it (see --region) rather than replacing this one.
          
          Minimum 2000, and a smaller value is rejected: synthetic reads cover event +/- 2000bp, so a narrower extraction window would leave original reads the synthetic ones are meant to replace outside it.
          
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

      --into-fastq <RAW_R1> <RAW_R2>
          Spike the events into the sample's raw FASTQ pair: the whole sample, spiked, is written to <output>/<NAME>_R1.fastq.gz and <NAME>_R2.fastq.gz (NAME from --fastq-prefix, default spiked), the pair to give a pipeline that starts from FASTQ. RAW_R1/RAW_R2 are every lane concatenated, from the run whose BAM --bam is; spike runs fastq.sh on them at the end. (R1.fq.gz and R2.fq.gz hold only spike's reads around the events.)

      --fastq-prefix <NAME>
          Name of the whole-sample pair --into-fastq writes: <output>/NAME_R1.fastq.gz and NAME_R2.fastq.gz [default: spiked]. A plain name, so many samples' pairs can sit side by side; -o picks the folder

      --indel-error-rate <INDEL_ERROR_RATE>
          Indel error rate per base in synthetic reads (fraction of total error that is indel rather than substitution). Default 0.0 means substitution-only. Typical Illumina: 0.0 to 0.05
          
          [default: 0]

      --gvcf <GVCF>
          Optional gVCF/VCF with SNP calls for LOH simulation. When provided, het SNP positions are extracted from this file to determine which reads belong to the deleted haplotype (more accurate than the default pileup-based approach). Supports .vcf and .vcf.gz (requires bcftools in PATH for .vcf.gz)

      --allow-overlap
          Allow events on the same chromosome that overlap or come within 7000bp of each other.
          
          By default they are rejected to keep event effects independent and truth interpretation unambiguous. The check is on each event's replacement footprint, not its span: an event replaces reads -- removes the originals and tiles synthetic ones -- across its span grown by 3500bp on each side (2000bp of haplotype flank + 1500bp of fragment). So two spans must be at least 7000bp apart, or each event's synthetic flank writes plain reference over the other event's edit and cancels it.

      --allow-resistant
          Simulate an event even when more than half the reads over it are ones spike cannot edit.
          
          Reads below --min-mapq, not in a proper pair, or with a mate that fails a filter stay in the merged BAM as they are, so the event that reaches the reads is about VAF x (1 - that share). Above one half, spike refuses by default rather than write a truth record the reads cannot back. truth.vcf records the share as SIM_RESIST either way. Under --edit-model origin a read in a pair origin may remove counts as editable too; the one-half threshold was set for --edit-model clean and is not yet checked for origin.

      --dup-model <DUP_MODEL>
          Duplication model: "full" (default) builds a full tandem haplotype with duplicated region appearing twice, producing both junction reads and correct depth increase from a single tiling pass. "junction" uses the legacy junction-only haplotype with separate depth copies
          
          [default: full]

      --edit-model <EDIT_MODEL>
          Which original reads an event replaces. "clean" (default): only the donor pool's pairs (both mates at --min-mapq or above, a proper pair, no duplicate, secondary, supplementary or QC-fail flag) inside the event's footprint, and the duplicates of those it removes. "origin" (experimental): every primary read at the event and at its look-alikes, each removed by its chance of having come from the edited copy, read from its MAPQ and its XA tag. It needs the aligner's XA tags (bwa-mem and bwa-mem2 write them)
          
          [default: clean]

  -h, --help
          Print help (see a summary with '-h')

  -V, --version
          Print version
```


## What spike refuses

spike would rather stop than write a truth VCF that does not describe the reads
beside it, so most input problems are a non-zero exit with nothing written
rather than a warning. Every refusal below is reachable from the command line;
the message is the exact text spike prints, measured by running it.

| What | Message | Checked in |
| --- | --- | --- |
| No events at all | `no events specified (use --event or --vcf)` | `main.rs` |
| `--allele-fraction` outside `(0.0, 1.0]` — including `0`, a negative, `inf` and `nan` | `allele-fraction must be in (0.0, 1.0]` | `validate_allele_fraction` |
| `--flank` below 2000 | `--flank 500 is too small: it must be at least 2000 so every original read replaced by synthetic reads is extracted` | `validate_flank` |
| `--dup-model` other than `full`/`junction` | `invalid --dup-model 'tandem', expected 'full' or 'junction'` | `main.rs` |
| `--region` with a 0 start (it is 1-based) | `region start must be >= 1 (1-based), got 0 in 'chr20:0-1000'` | `parse_region` |
| `--region` with start after end | `region start > end (5000 > 1000) in 'chr20:5000-1000'; check your interval` | `parse_region` |
| A long-read BAM (read length above spike's max fragment length) | `input BAM's read length (1501bp, its most common) exceeds spike's max supported fragment length (1500bp); spike simulates short paired-end reads and does not support long-read (PacBio/ONT) libraries` | `validate_read_length` |
| An event beyond the end of its chromosome | `DEL event start on chr20 is at or beyond chromosome length (99000000 >= 64444167)` | `validate_interval` / `validate_point` |
| An `--event` spec whose start is past its end | `del coordinate-based spec has start > end (38422500 > 38412500); check your interval` | `main.rs` |
| A `snp:` REF that is not what the reference has there | `REF allele mismatch at chr20:38412500-38412500: specified 'A' but reference has 'G'. Check that the position is correct (1-based in event spec) and matches the reference genome.` | `validate_ref_allele` |
| Two events whose replacement footprints (span ± 3500 bp) intersect, without `--allow-overlap` | `overlapping events detected (default is to reject overlaps).`<br>`Each event replaces reads across its span grown by 3500bp on each side (2000bp of haplotype flank + 1500bp of fragment), and two such replacement footprints may not intersect: two spans on one chromosome must be at least 7000bp apart.`<br>`Use --allow-overlap to override.`<br>`  - events 1 and 2 have intersecting replacement footprints on chrT: spans 10000-11000 and 12000-13000 (footprints 6500-14500 and 8500-16500)` | `main.rs` |
| A small variant at a site the sample already carries another allele at: at least 10 reads span it, and at least a fifth of them carry another allele there ([SNPs and small indels](#snps-and-small-indels)). Every such event is listed at once. | `the sample already carries another allele at 1 small variant(s): SNV  chr1:5001 T>A (38 of 38 reads). spike edits one copy and keeps the other copy's reads, so those reads keep that allele and the requested fraction cannot be reached (a hom-alt site asked for at 0.5 comes out all ALT). Remove these events from the input` | `carried.rs` |
| An event where more than half the reads over it are ones spike cannot edit (`SIM_RESIST` above 0.5), without `--allow-resistant` ([Reads spike cannot edit](#reads-spike-cannot-edit)). Every such event is listed at once. | `The reads over 1 event would carry less than half of what truth.vcf would claim, so spike stops rather than write it (RF8):`<br>`  DEL  chr20:27100001-27110000 (10000bp): 5731 of 5749 reads over it (99.7%) are ones spike cannot edit`<br>`Those reads are below --min-mapq, not a proper pair, or have a mate that fails a filter, and they stay in the merged BAM as they are. If most of them are below --min-mapq, lowering it lets spike edit them. Otherwise remove the events from the input, or pass --allow-resistant to simulate them anyway; truth.vcf then records the share as SIM_RESIST.` | `census.rs` |
| A donor pool under 30 read pairs ([Too few donor reads](#too-few-donor-reads)) | `event DEL  chr20:38412501-38422500 (10000bp) has too few usable donor reads: 0 read pair(s) extracted from chr20:38410500-38424500, fewer than the 30 spike needs (2097 read pair(s) in those windows were dropped for unusable base qualities and are not in that count). ...` | `finish_donor_pool` |
| No donor coverage at the event's breakpoints ([No donor coverage at the breakpoint](#no-donor-coverage-at-the-breakpoint)) | `event chr20:30000000-30010000 has no donor coverage at any of its breakpoints (chr20:29999999, chr20:30010000): the pool holds 6117 read pair(s) but none of them cover that. ...` | `simulate.rs` |
| `--edit-model origin` on a BAM with no `XA` tag in its first 100,000 records ([Editing hard spots](#editing-hard-spots---edit-model-origin-experimental)) | `--edit-model origin needs the aligner's XA tags (its alternative hits), but none of the first 100000 records of sample.bam carries one. bwa-mem and bwa-mem2 write XA by default, and a later step can strip it. Re-align with one of them, or use --edit-model clean.` | `origin.rs` |
| `--edit-model origin` with no origin depth at any of the event's breakpoints | `event over chr20:14603000-14608000 has no origin depth at any of its breakpoints (chr20:14604999, chr20:14606000): no read at the spot or at its 0 look-alike region(s) could have come from there, so spike would invent the reads it plants and still write a truth VCF beside them.` | `simulate.rs` |
| A `--reference` FASTA that is gzip-compressed but not named `.gz`/`.bgz` | `misnamed.fa is gzip-compressed (starts with the gzip magic bytes 1f 8b) but is not named .gz/.bgz, so it would be read as raw uncompressed sequence; rename it to end in .gz or .bgz with a matching .gzi index, or decompress it first` | `reference.rs` |
| A gene or exon `--event` names that the `--exon-bed` has not got | `gene 'NOSUCH' not found. Available: GENEA, GENEB` | `exon.rs` |
| A `--gvcf` that cannot be read — missing, not bgzipped, no index, or no `bcftools` for a `.vcf.gz` | `could not read the --gvcf 'bad.vcf.gz'; spike stops rather than simulate without the sample's SNPs`<br>`Caused by: bcftools exited with status exit status: 255 on gVCF 'bad.vcf.gz': Failed to open bad.vcf.gz: not compressed with bgzip. The run stops here.` | `loh.rs` / `simulate.rs` |
| A read whose quality string does not match its sequence, or holds a byte outside `!`-`~` | `read <name>/1 has 150 quality byte(s) for 151 base(s); refusing to write invalid FASTQ` | `write_paired_fastq` |

Under `--edit-model origin` the resistant-reads refusal says the reads are `neither in the donor pool nor in a pair origin may remove`, and ends with the note that its threshold was set for `clean`. That text is checked by the unit test `census::tests::test_under_origin_the_messages_say_what_origin_counts_and_that_thresholds_are_unchecked`, not by a run like the rows above.

`--exon-bed` has a family of related refusals of its own — a duplicate exon
number under one gene, an exon range whose start is after its end, a malformed
BED line — each naming the gene and the line; see
[Gene/exon-based events](#geneexon-based-events).

An unreadable `--gvcf` used to be a warning: spike suppressed that region's
reads at random and exited 0, so the output could not be told from a sample
with no SNPs there. It is a refusal now (CR3 in
`docs/review/CLINICAL_SV_DESIGN_NOTES.md`). Measured on the HG002 BAM with a plain-gzip
`--gvcf`: exit 1, and the output directory is left empty.

**Still not a refusal:** without `--gvcf`, spike reads the sample's SNPs by
pileup from the BAM. If that read fails, spike logs `could not read the
sample's SNPs in <region>: ...` and goes on with no SNPs for the region, as
before. The same holds when a readable `--gvcf` has no het SNPs in the region
and the pileup it falls back to fails. That gap is known and not yet closed.
See [The sample's SNPs from a gVCF](#the-samples-snps-from-a-gvcf).

## Architecture

```
main.rs          CLI, event parsing, orchestration
types.rs         Core types: SimEvent, ReadPair, ReadPool, SimConfig
haplotype.rs     Variant haplotype construction (segment-based)
simulate.rs      Read suppression + synthetic read tiling
origin.rs        --edit-model origin: each read's chance of coming from the edited copy, look-alikes, origin depth;
                 under either model, the duplicates of a removed pair
census.rs        Reads spike cannot edit (SIM_RESIST) and the warn/refuse thresholds and messages
synth.rs         Quality-profiled synthetic read generation
extract.rs       BAM/CRAM read pair extraction
stats.rs         Fragment length distribution
loh.rs           The sample's SNPs: two copies, phasing, read assignment
exon.rs          Exon BED and event spec parsing
vcf_input.rs     VCF input parser (DEL/INS/DUP/INV/BND/SNP/indel)
truth.rs         Truth VCF output
fastq.rs         Gzipped paired FASTQ writer
reference.rs     Indexed FASTA reading + in-memory sequence store
bam_stats.rs     BAM/CRAM read length, adapter trimming, read-name shape and single-end check
read_name.rs     spike's read names: SPIKE_ in the shape of the input's names
validate.rs      `spike validate` subcommand: automated spike-in quality checks
```

## Simulation model details

All SV types share a common simulation framework: (1) build a linear variant haplotype from segments, (2) suppress original reads in the affected region at the VAF rate, (3) tile synthetic reads across the haplotype to replace the suppressed fraction. The details of what the haplotype looks like, how reads are suppressed, and what observable signals result differ by variant type.

A segment whose flank runs off the end of a contig is truncated to what the reference actually holds: the haplotype is shorter there, its reference footprint stops at the contig end, and the count of synthetic reads follows the shorter haplotype.

spike sequences each read for the input library's cycle count, the most common read length among the first 50,000 primary records, and trims reads the way that library's were. A library counts as **adapter-trimmed** when at least 1,000 of those records are full length and under 5% of them end in `A`, the TruSeq adapter's first base. On the HG001 and HG002 NovaSeq X 30x BAMs this share is 0.000; on an untrimmed HG002 BAM it is 0.300.

- **Adapter-trimmed library:** a fragment shorter than the cycle count gives two reads as long as the fragment, and the fragment model keeps every insert size up to 1500bp. When the fragment runs past both reads, each read then loses the longest 3' end that spells a prefix of `AGATCGGAAGAGC`. That is what the real reads of those BAMs show: about 61% are 151 bp and 27% are 150 bp, and 99.3-100% of the 147-150 bp ones are followed by exactly the adapter bases they lost (docs/superpowers/plans/2026-09-28-read-length.md).
- **Otherwise:** every read is the cycle count long, and no fragment is shorter than one read.

spike caps the fragment lengths it draws at 1500bp, so it does not support long-read (PacBio/ONT) libraries. If the input BAM/CRAM's read length exceeds 1500bp, spike exits with an error before doing any work rather than generating reads that don't match the library.

For the same reason spike needs a **paired-end** library: every donor read comes from an extracted pair, and the quality profile is trained on R1 and R2 separately. A single-end BAM/CRAM yields no donor material, so spike exits with an error naming the file rather than emitting a handful of synthetic pairs built on an empty quality profile. Single-end is detected from SAM flag `0x1`, which a paired library sets on every read whether or not the pair aligned properly — a paired BAM with no proper pairs in it passes this check. Reaching 50,000 primary records with no `0x1` on any of them **is** how the scan establishes that, so it stops there instead of reading the file to the end; a paired file stops there too, since the scan only needs the read length (the fragment lengths spike draws come from the extracted donor pairs, not from this scan). The message says whether the cap bound: at end of file the count it reports is every primary record in the file and the verdict is certain, while at the cap it is a window, and the message says so — a file that does hold pairs but puts 50,000 consecutive primary records without `0x1` at its head would be reported the same way. That is bounded rather than impossible: spike's input is a coordinate-sorted indexed BAM/CRAM, where a mixed library's single- and paired-end reads interleave at every locus, so no supported input is known to reach it.

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
The inserted sequence is either user-specified or randomly generated. It has no reference origin (novel sequence). A generated sequence is kept on the event, so the truth VCF spells it out: `REF` is the anchor base at `POS` and `ALT` is that base plus the inserted bases (see [Truth VCF](#truth-vcf)).

**Simulation**: Suppress-and-replace. Reads spanning left_flank→novel or novel→right_flank produce chimeric reads at the insertion point. Tiling draws each fragment start from the starts that overlap reference sequence, so no fragment is placed entirely within the novel sequence (such reads wouldn't align to the reference at all).

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

**What fusion mode is not: a balanced translocation.** It adds evidence for
**one** join and removes nothing. A balanced heterozygous translocation would
change one copy of each partner, keep the other copies, make **both**
derivative joins, and keep copy number the same. Fusion mode does none of
those four. The two BND records in the truth VCF are the two ends of the one
join, not the reciprocal one. At `af=0.5` the junction gets as many added
fragments as the local depth, on top of all the original reads, so its
evidence is stronger than a copy-neutral sample would show. Use it to test
whether a caller finds a junction; do not read its results as sensitivity or
genotyping for balanced translocations. spike says this in a warning on every
run that has a fusion event. See CR6 in `docs/review/CLINICAL_SV_DESIGN_NOTES.md`.

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

This is the default, `--edit-model clean`. The experimental `--edit-model origin` replaces this step. It decides for every primary read at the event and at its look-alikes, in the donor pool or not, and removes each by its chance of having come from the event's footprint times this same per-copy rate; see [Editing hard spots](#editing-hard-spots---edit-model-origin-experimental).

## Synthetic read tiling

The number of synthetic reads to tile is:

- **Non-additive events**: `n = round(coverage * VAF * starts / mean_fragment_length)`, where `starts` is the number of fragment start positions tiling can use: `haplotype_length - mean_fragment_length`, minus starts that would lie wholly inside inserted sequence. Starts are uniform over exactly that set, so the flanks get `VAF * coverage` synthetic depth, replacing what was suppressed.
- **Additive events** (breakpoint-only tiling): every original read is kept, so `n = round(coverage * VAF / (1 - VAF))` per breakpoint makes junction fragments a `VAF` fraction of the depth there (VAF capped at 0.95 -- a request above it is simulated at 0.95, and the truth VCF records `SIM_VAF=0.950` with the request in `SIM_REQ_VAF`).

`mean_fragment_length` is the library's own mean, taken over the same range the generator samples in -- `[read_length, 1500]` -- so the count is normalised by the distribution that is actually emitted. An observed insert size outside that range can never be generated, so it is left out of the model: on HG002 chr20 that is 38 of 4,595 donor pairs (all of them shorter than one 151 bp read), and it moves the mean from 418.6 to 421.0 and the count planted for `del:chr20:38412500-38422500` from 291 to 289. `coverage` is the donor pool's mean fragment depth in a 2 kb window around the first breakpoint, counted only on that breakpoint's own chromosome: a fusion's pool holds both partners, and reads from the far side would otherwise be added to the near side's depth. Under `--edit-model origin` it is the origin depth in the same window instead (see [Editing hard spots](#editing-hard-spots---edit-model-origin-experimental)), because inside a perfect twin the pool holds no read.

### One depth for the whole event

That one `coverage` scales every fragment the event tiles, wherever it lands. Where the donor's own depth differs from it, the event's depth there is wrong: a het DUP over a stretch at a quarter of the breakpoint's depth came out **81.09x** where a local CN2 -> CN3 is **28.13x** (`CR2` in `docs/review/REVIEW.md`). spike does not fix that yet, but it measures it (CR2 option B). Each reference interval the event's fragments are drawn from is cut into ~1 kb bins, and each bin's donor depth `D` is measured the same way as `coverage` (`C`). The event's **depth fold** is the largest `max((D+1)/(C+1), (C+1)/(D+1))` over its bins. The log prints it with the bin it came from, the run README has a **Depth fold** column, and `truth.vcf` records it as `SIM_DEPTH_FOLD`, which `spike validate` reads back as its advisory `depth_fold` row. Above **1.5** spike warns:

```
DUP  chrT:10001-28000 (18000bp): the donor's depth over chrT:17000-18000 is 25.0x, but every fragment this event tiles is scaled by the 100.0x measured at one of its breakpoints (3.88-fold). Where the donor's depth differs from that, the event's depth there is wrong by about that much; truth.vcf records the fold as SIM_DEPTH_FOLD (CR2).
```

(The depths are fragment depths, which is why the review's 75x read depth shows as 100x.) On ordinary loci the fold is not 1: measured on 40 seeded 10 kb DUPs inside the HG002 SV benchmark on chr20 (35x), it ran from **1.11** to **2.46**, median **1.27**, and **6 of the 40** warned. `D` comes from the donor pool, which holds only reads at `--min-mapq` or above, so a bin of low mappability reads thin whether or not the library is; the fold counts it anyway, and some of those six may be that. Under `--edit-model origin`, `D` and `C` are both origin depths, which count low-MAPQ reads by their chance. The 1.5 threshold was set for the pool's fold and is not yet checked for origin's.

Fragment lengths are sampled from the empirical distribution of the donor reads, which holds the observed insert sizes in `[read_length, 1500]` and nothing else: generation draws from exactly that range, so the model and the generator describe one distribution rather than two. If no donor insert size falls in it, spike warns and falls back to a 400 +/- 80 default clamped into the same range (never a distribution it cannot draw from). Each fragment is placed at a random position on the haplotype and a read pair is synthesized with quality scores from the learned Markov model. The fragment's left end is always read forward and its right end reverse (FR), and a coin flip out of the same seeded stream decides which of the two is R1: about half the pairs come out F1R2 and half F2R1, as in a real library, so read-orientation filters (Mutect2's, for one) see a balanced strand mix. Whichever mate is R1 is sampled from the R1 quality model.

Each fragment comes from one of the sample's copies and carries its alleles. Up to VAF 0.5 all fragments come from the event copy. Above 0.5 the other copy gives `max(0, 2*VAF - 1) / (2*VAF)` of them, matching what was suppressed from it.

For additive events, fragment placement is restricted to positions that cross a segment boundary (breakpoint). For non-additive events, placement is uniform across the haplotype; where the haplotype carries novel (non-reference) sequence, the start is drawn from the starts whose fragment overlaps reference, in proportion to those intervals' lengths for the fragment length drawn, so a fragment is never placed entirely within novel sequence.

## The sample's SNPs

For every event, spike reads the sample's SNPs over the haplotype's whole reference footprint (event ± 2 kb; each side of a fusion separately), at every VAF. Het SNPs are phased and assign reads to copies; het and hom-alt alleles go into the synthetic reads. Two strategies are supported:

### Pileup (default)

A pileup counts the bases (A/C/G/T) at each position (MAPQ >= threshold, excluding secondary/supplementary/duplicate/QC-fail). Positions with >= 10 reads are called:
- **het**: two alleles each at 20%-80%,
- **hom-alt**: one allele at >= 90% that differs from the reference.

A second pass records each fragment's bases at the het SNPs only, so memory grows with the number of het SNPs, not with the region's size.

### gVCF (optional, `--gvcf`)

SNPs are loaded from a pre-called VCF (e.g., DeepVariant gVCF): het (`0/1`, `1/0`) and hom-alt (`1/1`) biallelic SNPs. A het SNP whose base a deletion on the other haplotype removes (a carried deletion allele, or ALT `*`) has no copy carrying REF, so it counts as hom-alt. A BAM pass then records which allele each read carries at the het SNPs. Phased genotypes (`0|1`, `1|0`, with an optional `PS` phase set) are used for phasing (below); a phased VCF such as a GIAB/T2T benchmark gives the most realistic result. If the gVCF has no het SNPs in the region, spike falls back to pileup.

Whether the gVCF can hold SNPs for the chromosome at all is read from its own `##contig` names, not guessed from an empty result: `bcftools view -r chr20:...` on a file that calls that chromosome `20` prints nothing and exits 0, exactly like a region that genuinely has no SNPs. spike reads the header itself (plain or bgzipped, no index needed) and warns only when the names cannot match — the other convention (`20` for `chr20`) or no such chromosome at all. A region that simply has no SNPs in it is not warned about. Each message says what happens next: a fallback to pileup when the file was read, and — when it could not be read at all, e.g. an unindexed `.vcf.gz` that bcftools refuses, with bcftools' own reason quoted — that the run stops. It used to go on with that region's original reads suppressed at random; see [What spike refuses](#what-spike-refuses).

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

**A pool under 1,000 read pairs gets a warning** (the run still goes on). The
Markov levels need 30 observations *after a low-quality base* at each cycle,
and a small pool rarely has that many. Sampling then falls back to the levels
with no memory of the previous quality, so low-quality bases come out
scattered instead of in runs. On HG002 35x this was measured against held-out
reads from the same window, with a tolerance set by how much real reads from
ten other chr20 windows differ from them. Pools of 30-500 pairs missed
low-quality persistence by 0.15 against a tolerance of 0.10; at 1,000 pairs
it was 0.04. Per-cycle mean quality and the share of bases under Q20 were
already within tolerance from 60 pairs. The warning names the pool size and
the bin census, for example:

```
WARN spike::synth] Quality profile learned from 625 donor pairs; below 1000 its low-quality runs come out shorter than the sample's (measured on HG002 35x). Base-conditioned bins: 1208/1208 usable. Markov bins: base 1200/4832, cycle 445/1208 usable. Widen --flank or --region for a larger pool.
```

That is `snp:chr20:38600002:G:A --flank 2000` on the HG002 35x BAM; the
default 10 kb flank gives 3,032 pairs and no warning.

Carrying Q2 forward costs the bases *after* an `N` nothing in practice. Measured over 40 seeds at a 1 bp reference gap on chr2 (an `N`-containing read there is 98.3% real sequence), the non-`N` bases of `N`-containing reads average **Q35.621** when the chain carries Q2 and **Q35.615** when it does not — a difference of +0.006 Q against a 0.022 standard error. The reason is that real Illumina Q2 is rare enough (6 of 389,429 donor bases in that window) that the after-Q2 transition bin never reaches the 30-observation threshold, so sampling falls straight through to the non-Markov levels.

### Too few donor reads

Every simulated read is built from the donor pool extracted for its event -- the quality profile above, the fragment-length distribution and the coverage the tiling count is scaled by all come from it. A pool holding fewer than **30 read pairs** after deduplication is refused: the run exits non-zero, naming the event and the windows it searched, and writes nothing.

```
Error: event DEL  chr20:30000001-30010000 (10000bp) has too few usable donor reads: 0 read
pair(s) extracted from chr20:29990000-30020000, fewer than the 30 spike needs (0 read pair(s) in
those windows were dropped for unusable base qualities and are not in that count). ...
```

The drop count is in the message because it is one of the ways a pool empties:
a CRAM that stores its qualities as read features has every record dropped by
the quality check below, and the pool is then 0 through no fault of the region
or of `--min-mapq`. Measured on the chr20 slice with every `QUAL` set to `*`:
`0 read pair(s) extracted from chr20:38410500-38424500, fewer than the 30 spike
needs (2097 read pair(s) in those windows were dropped for unusable
base qualities ...)`. When nothing was dropped the count is `0` and the other three
causes the message lists are the ones to look at.

Without that check a starved window is silent. spike logs `Built read pool: 0 pairs`, then falls through to every substitute in turn -- the constant Q20 last resort above for every base, the default 400 +/- 80 fragment distribution (clamped into the sampled range) -- and exits **0** with a truth VCF and no reads behind it. (The 2-read tiling floor is no longer part of that fall-through: `compute_tiling_count` returns 0 outright at zero coverage rather than floored to 2, so nothing is invented -- but a truth record with nothing behind it is the same wrong answer, quieter.) An event in a zero-coverage region, an off-target panel BAM and a mistyped `--region` all reach it.

30 is the observation count the quality model itself requires before it will sample from a bin, and a pool of *n* pairs puts *n* observations in each cycle-only bin (level 4 above), so it is the smallest pool at which any level of the model is trained to its own threshold. It is a floor on "measured from this library at all", not a coverage requirement: the 30 kb window of `del:chr20:38412500-38422500` yields 4,595 pairs on the 35x HG002 BAM, so it would have to fall to roughly 0.2x before 30 pairs bound. If a real event does sit in a region that thin, widen `--flank` or `--region`, lower `--min-mapq`, or use a BAM that covers it.

### No donor coverage at the breakpoint

The 30-pair floor above counts the pool as a whole, summed over every window
the event was extracted from. That is not the same question the tiling count
asks: the number of synthetic fragments is `coverage x VAF`, and `coverage` is
measured in a 2 kb window around a breakpoint. A pool can clear 30 pairs and
still measure coverage 0 where the event is -- `--region` pointing somewhere
the event is not, or a fusion whose other partner carries the whole pool. The
tiling count then collapses to its floor of 2, and spike used to write those 2
invented pairs and a truth VCF beside them and exit **0**.

An event with no donor coverage is now refused the same way a starved pool is:
the run exits non-zero and writes nothing. Every breakpoint of the haplotype is
measured, both sides of each -- the last reference base before the cut and the
first after it -- so the verdict never depends on the order the event names its
parts. **How many of them have to be covered depends on how many places the
donor reads came from:**

- A **fusion** is extracted from two loci, one per partner, and every fragment
  spike plants spans the junction between them. A partner with no donor reads
  makes half of every planted read invention, so **both** sides of the junction
  must be covered. `fusion:GENEA:exon1:GENEB:exon2` and
  `fusion:GENEB:exon1:GENEA:exon2` over the same two loci now give the same
  answer.
- **Every other event** (DEL, DUP, INV, INS, small variants) is extracted from
  one locus and tiled across its whole haplotype, so it needs donor depth
  *somewhere* around it: it is refused only when **no** side of any of its
  breakpoints is covered. One uncovered side is ordinary input -- a thin spot,
  or the far side of an event that straddles the edge of a sliced or panel BAM.
  `del:chr20:41490000-41600000` on a BAM whose reads stop at 41,500,000 is
  simulated, and so is the mirror image whose *near* side is the uncovered one;
  the tiling count is then scaled by the first covered side's coverage.

**Keeping such an event is not silent.** The haplotype still spans every side,
so the fragments tiled across the bare part are scaled by donor depth measured
at the *other* side and land where the input BAM has no read -- a coverage
island the input does not have. spike logs a `WARN` naming the bare sides, and
the generated run README repeats it under the `Events` table, so a run whose
log has scrolled past still says so. Measured on the chr20 37.5-41.5 Mb slice
(reads start at 37,499,851), `del:chr20:37400000-37510000 --seed 1`: the
4000 bp haplotype is `[37,398,000-37,400,000) | [37,510,000-37,512,000)`, the
left half is outside the slice, and **258 of the 516 synthetic records** land
there, taking chr20:37,398,000-37,400,000 from **0x** in the input BAM to
**19.4x** in the output. The FASTQ, truth VCF and `replaced_reads.txt` are
byte-identical to the run before the warning existed.

```
WARN spike::simulate] event chr20:37400000-37510000 is kept although the donor pool
has no reads over chr20:37399999: one bare breakpoint side is the far edge of a sliced
or panel BAM, not a reason to refuse. But the haplotype spans every side, so the
fragments tiled across the bare part are scaled by the 60.7x measured elsewhere and
land where the input BAM has no read -- a coverage island the input does not have.
Extract a wider BAM, or narrow the event, if that matters downstream.
```

```
Error: event chr20:30000000-30010000 has no donor coverage at any of its breakpoints
(chr20:29999999, chr20:30010000): the pool holds 6117 read pair(s) but none of them
cover that. ...

Error: event chr20:38420200-38420200 has no donor coverage on one side of its junction
-- every read spike plants for a fusion spans the junction, so each partner needs donor
reads of its own (chr20:30005000): the pool holds 3022 read pair(s) but none of them
cover that. ...
```

When the coverage is real but thin, the floor still applies -- a haplotype
shorter than one fragment asks for 0 fragments however good the coverage is,
and planting nothing would leave a truth VCF with no reads behind it. The floor
then plants *more* support than `--allele-fraction` asked for, so spike warns
with both numbers **and** the truth VCF records what was planted: `SIM_VAF` is
the fraction the two fragments make up, `SIM_REQ_VAF` the request.

```
WARN spike::simulate] coverage 3.8x at VAF 0.030 asks for 1 tiled fragment(s); spike
emits the 2 it needs to plant the event at all, so the realized allele fraction is
above the 0.030 requested; the truth VCF records the realized fraction as SIM_VAF and
the 0.030 requested as SIM_REQ_VAF
```

That run's truth record reads `SIM_VAF=0.058;SIM_REQ_VAF=0.030`: two fragments
of the library's mean 400 bp, over 3.8x coverage across the haplotype's 3,600
usable start positions, are 5.8% of the depth there, not the 3.0% asked for.
Feeding that record back through `--vcf` asks for 0.058, gets the same two fragments and
writes `SIM_VAF=0.058` again, with no warning -- the request a thin region can
actually be spiked at.

### Missing or unusable donor base qualities

SAM's QUAL field is all-or-nothing per record: a read either has a full quality string or none at all (`*`). The two containers encode `*` differently and spike sees both — BAM keeps a per-base array with every byte `0xFF`, while noodles normalises a CRAM record's all-`0xFF` buffer to an *empty* one before spike ever sees it. A record with either shape, or carrying a raw quality above Q93 (the SAM maximum — a malformed record), is dropped during extraction: it is not included in the learned quality profile and does not contribute a read pair to the output. Dropped records are counted and logged as a warning (spike logs to stderr), e.g.:

```
WARN spike::extract] 457 record(s) considered for the chr20:38402500-38432500 donor pool had no quality scores (SAM '*'); their pairs were dropped from the pool and merge.sh removes them from the merged BAM
```

**A dropped pair is removed from the merged BAM even though nothing replaces it.** spike has no quality string to write for it, so it cannot come back through the FASTQ — but it is still listed in `replaced_reads.txt`, so `merge.sh` drops it. That costs real depth. Leaving it in costs more: the pair would sit inside every event it overlaps as reference support that no allele fraction can suppress and no synthetic read can replace, so the realised VAF would come out diluted by the unusable-quality fraction while the truth VCF still claimed the full one. Because the same pairs are dropped across the whole extraction window, not just inside the event, the loss cancels in a flank-normalised ratio; the dilution would not. Measured on an HG002 chr20 slice with 5% of donor pairs' quality stripped to `*`, over a 10 kb deletion: 147 of 2868 eligible records inside the event (5.13%) carried `*` quality and were un-suppressible; removing them takes that to 0.00%, at a cost of those same 147 records of depth.

**The count is in the generated run README, per event and in total**, not only on stderr: a run whose log has scrolled past leaves the depth dip invisible otherwise, and a localized dip is exactly what a depth-based CNV caller reads as signal. The `Events` table gains a `Dropped (unusable quality)` column, and a line under it gives the total and how many of those are removed with nothing put back in their place. Measured on the chr20 slice with `QUAL` set to `*` on every 20th record — **5.0%** of records (709 of 14,184) — over a 10 kb deletion with `--flank 2000`: the pool is 1855 pairs and **201** pairs are dropped, **9.8%** of the 2056 pairs the window would otherwise have held, because a pair is lost when *either* of its two records is stripped. The run README now says so; the FASTQ, truth VCF and `replaced_reads.txt` are byte-identical to the run before this column existed.

Without the drop, a missing quality used to decode to an invalid FASTQ quality byte (space) for kept reads and poisoned the learned quality model, so synthetic reads sampled from it came out mostly-Q0 with effectively random bases. On the CRAM path it produced a FASTQ record with a full-length SEQ line next to a zero-length QUAL line, which `samtools import` rejects outright (`truncated file`).

As a last line of defense, `write_paired_fastq` validates every pair *before* it creates either output file — so a refusal leaves no half-written `.fq.gz` in `--output` — and returns an error if a quality string's length does not match its SEQ, or if any byte falls outside the printable Phred+33 range `!`-`~` (33-126).

`write_paired_fastq` also flushes each gzip stream's returned writer explicitly after `finish()`, rather than only letting it drop. `finish()` empties flate2's own internal buffer into the `BufWriter` it hands back, but that write can land entirely inside the `BufWriter`'s own buffer without ever reaching the file — a write failure (e.g. a full disk) then surfaces only when the `BufWriter` is dropped, where `Drop`'s own flush swallows any error. Calling `write_paired_fastq` with `R1.fq.gz` symlinked to `/dev/full` and a single small read pair returned `Ok(...)` before this fix and now returns `Err("No space left on device (os error 28)")`. Both streams' `finish()`+`flush()` are computed unconditionally before either error is propagated, so an R1 failure can't short-circuit past R2's — leaving R2's `GzEncoder` to be dropped unfinished would rediscover the exact same swallowed-error bug for R2 alone; when both streams fail, the reported error now names both.

### Indel error model

When `--indel-error-rate` is set above 0, a fraction of sequencing errors are modeled as insertions or deletions (50/50 split) rather than substitutions. This keeps each read's length: insertions consume an output position without advancing the reference, and deletions skip a reference base without consuming an output position.

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
fails the run when any VAF's recall *outside the background control*
(`Recall_outside_bg`, below) is less than `F`, and `--min-gain N` is how many
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
`<outdir>/background_control/`. The version each tool reported is written to
`<outdir>/tool_versions.tsv`, one `tool<TAB>version` line per tool, before any
step runs; a tool that answers neither `--version` nor `version` is recorded as
`unknown`.

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
see `N1` in `docs/review/REVIEW.md`.

#### What the harness establishes, and what it does not

`scripts/validate_pipeline.sh` is an **integration and regression test**. What
it establishes is that spike's whole loop runs on real data, and that the
caller's result moves when the spike-in is present. It is not a measure of a
caller's sensitivity or precision, and the numbers in
`validation_summary.tsv` should not be quoted as one.

It is a deletion test with one caller. Step 2 keeps only `SVTYPE=DEL` records
that are het and 500-50,000 bp, and both `delly call` invocations run with
`-t DEL`. Nothing in the run plants or calls an INS, DUP, INV or BND, so a
green run says nothing about those types — or about any caller other than
Delly.

The verdict is a gain-over-background gate, and a gain is not an
identification. It establishes that the injection changed the caller's output;
it does not say *which* event was added. The measured example above is exactly
that gap: at `--vafs 0.001` the run gains one event over the control, and that
event is `sim_del_8`, a DEL the background already carries which re-alignment
flipped from FN to TP.

Precision against this truth set is not a clinical precision. Truvari scores
Delly's PASS calls -- the script passes `--passonly`, so non-PASS records are
dropped from both sides -- against spike's truth VCF, which holds the **added
events only**, so a call matching a real variant of the background sample is counted
FP for being absent from a truth set that never described the background.
Recall is the column that means something here — how many planted events the
caller recovered — and even that is read against the background row rather than
against zero, for the reason given above. The last three columns do that reading
for you (RF12). `N_outside_bg` is the number of truth DELs the background control
does *not* recover, `TP_outside_bg` how many of those the spiked run recovered, and
`Recall_outside_bg` the ratio. A truth DEL the background already carries leaves
both sides, because recovering it says nothing about the spike-in. Measured on the
default run (NA18488 chr20, 20 truth DELs):

| VAF | Recall | Recall_outside_bg |
| --- | --- | --- |
| background | 0.3500 (7 of 20) | n/a |
| 0.5 | 0.8500 | **0.7692** (10 of 13) |
| 0.25 | 0.7000 | **0.5385** (7 of 13) |
| 0.1 | 0.6500 | **0.4615** (6 of 13) |

The exclusion is only as good as the control's caller. The background carries more
of the truth DELs than Delly finds: 4 more (570-1527 bp, in repeats) show
reduced depth in NA18488, and HG001's long-read calls confirm a deletion at each
site, but no short-read caller tried finds them (Delly, Manta, TIDDIT,
CNVpytor). Those 4 stay in `N_outside_bg`, so `Recall_outside_bg` can still count
a background deletion as the spike-in's. On the default run it does: 2 of the 10
at VAF 0.5 and 1 of the 7 at VAF 0.25 are among those 4, and none at VAF 0.1. See
RF12 in `docs/review/REVIEW.md`. The converse is no safer: a locus at
which the caller made no call has not been shown to be variant-free. The run
measures what this caller recovered from this background, not what is in it.

Genotype accuracy is not measured at all. Matching is whatever Truvari makes of
the flags the script passes it — `--passonly -r 500 -p 0.5 -P 0.5 -s 500`
(`scripts/validate_pipeline.sh`) — and the script does not pin a Truvari
version. It records the one it ran in `<outdir>/tool_versions.tsv`, so read
those flags against that version rather than against this paragraph. What the script itself does is
fixed: it reads back `TP-base`, `FP`, `FN`, recall, precision and F1 from
Truvari's `summary.json`, and no genotype is compared.
