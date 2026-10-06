# Why real reads are soft-clipped, and what spike should model

October 2026. An offline analysis (Python, outside spike) of the soft clips in real Illumina reads. It asks which clips a read simulator has to reproduce. spike's own reads are soft-clipped far less often than real reads (0.2% against 2.4% in an earlier inspection of 20 transplanted deletions).

## Data

- **Reads:** a chr20 slice (38.5–40.2 Mb) of the GIAB HG002 NovaSeq PCR-free 35x BAM, aligned with bwa-mem2 2.2.1 to GRCh38 no-alt (no decoy sequences).
- **Counted:** primary, mapped, non-duplicate reads: 488,443. Of these, 13,442 (2.75%) carry at least one soft clip.
- **Variant list:** HG002's own large variants come from the GIAB T2T-Q100 v1.1 benchmark.

## Causes

Each clip is put in one cause by the checks in `clips/` (`classify2.py`, `foreign.py`, `foreign_sw.py`). `summary.py` then gives each read one cause; a read with clips of two kinds counts as the less ordinary one.

| Cause | How it was recognised | Clipped reads | All reads |
|---|---|---|---|
| Bad end | Clipped bases are reference sequence read badly: ≥ 50% identical where they would align, or poly-G, or 1–4 bp, or aligning within 300 bp once gaps are allowed | 44.7% | 1.23% |
| Adapter | A 3′ clip holding the TruSeq R1 or R2 adapter (or another Illumina adapter or primer 12-mer), or the read spells the adapter start where its fragment ends | 25.7% | 0.71% |
| The sample's own variant | ≥ 3 reads clip at the same place (±2 bp), or within 150 bp of a ≥ 10 bp HG002 variant | 19.7% | 0.54% |
| Chimeric fragment | The clip realigns elsewhere (MAPQ ≥ 20 or SA tag), the mate maps > 1.5 kb away or on another chromosome, or the clip is an inverted local copy | 4.7% | 0.13% |
| Not in the reference | Non-reference sequence that does not align locally and maps nowhere. Of the unexplained clips ≥ 20 bp, with relaxed bwa-mem2 settings (`-k 13 -T 20 -B 2`): 51% map only partly, 34% not at all | 3.6% | 0.10% |
| Simple repeat | One base ≥ 60% of the clip, or a 1–3 bp unit over ≥ 80% | 1.6% | 0.04% |

- **Typical clips.** First-pass bad ends: median 26 bp, mean quality 22, 84% at the 3′ end. Adapter clips: median 15 bp, mean quality 35, all at the 3′ end. "Not in the reference": median 53 bp, mean quality 23, 79% at the 3′ end.
- **First pass.** `classify.py` knew only the R1 adapter. `classify2.py` fixed that, but still left 22.6% of clipped reads as "foreign": non-reference and not mapping elsewhere. `foreign.py` and `foreign_sw.py` sort those into the rows above.

## What it means for spike

- **Bad ends.** spike's quality model is a first-order Markov chain, so a read cannot stay bad. spike makes almost no bad ends.
- **Adapter.** spike draws fragment lengths from [read length, 1500] for a library that is not adapter-trimmed, so its reads never run into adapter.
- **The sample's own variants** are a property of the sample, not of sequencing. There is nothing to simulate.
- **Chimeric fragments, sequence not in the reference and simple repeats** together come to about 0.3% of reads: library and background noise.

Proposal under discussion: reproduce only the bad ends, and keep this record of what the other clips are.

## Quality and error models for bad ends (`qmodel/`)

### Setup

- Real reads come from 25 windows of the SNV run (5,509 training pairs, 5,354 held-out pairs).
- For the soft-clip tests, reads were rebuilt from the reference at the held-out pairs' own positions, so only the quality and error model differs. They were then aligned with bwa-mem2.

### Quality strings (`e3.py`, `e7.py`), held-out real reads against models

| | Read-to-read SD of mean quality | Perfect reads | Crashed 3′ ends | R1–R2 correlation |
|---|---|---|---|---|
| Real | 2.02 | 5.7% | 1.04% | 0.55 |
| First-order chain (spike today) | 0.52 | 0.03% | 0% | 0.00 |
| Copy a real pair's two quality strings | 1.91 | 5.6% | 0.81% | 0.53 |
| fqzcomp-style context as a generator | 2.04 | 5.1% | 0.80% | 0.54 |

The fqzcomp context is htscodecs fqzcomp_qual's: last 5 qualities, position from the 3′ end, a change flag, and a per-read selector.

### Bits per quality, held out (`e5b.py`, `e8.py`)

Trained on 200k pairs. The selector's cost is included.

| Model | Bits per quality |
|---|---|
| fqzcomp context | 0.383 |
| + the base (fqzcomp5-like) | 0.379 (−1.1%) |
| Bigger contexts | 0.386 |
| Bigger contexts + 3 bases | 0.391 |

The bases carry about 1% of what can be predicted. Bigger contexts lose to sparse data.

### Errors

- At the same quality, poor reads err far more than good reads. Aligned bases only: Q25 0.07% against 0.95%.
- With clips of other causes masked from learning and scoring (`nonerr_mask.py`, `e4m.py`, `e9m.py`, `score_masked.py`):

  | | Error-type clipped reads | ≥ 5 errors in the last 20 cycles |
  |---|---|---|
  | Real | 0.97% | 0.53% |
  | fqzcomp-style + an error table by quality, read class, low-quality run and distance from the 3′ end | 1.07% (0.89–1.29) | 0.38% |
  | Copying each donor pair's error positions | 1.28% (1.08–1.52) | 0.63% |

### A mistake along the way

The first model tests counted every clipped base as a sequencing error (`clumping.py`, `score_e4.py`). They concluded that real errors clump at read ends and that copying error positions matched real clips (2.62% against 2.34%). That copy only matched by turning adapter and foreign sequence into random substitutions. See `.claude/judgment-gate-cases.md`, 2026-10-06.

The masked tests used the first-pass causes. Under the final causes, short clips and gapped local alignments also count as bad ends, so the real error-type share would be nearer 1.2% than 0.97%. The masked tests were not re-run on the final causes.

## Limits

- One sample, one 1.7 Mb region and one aligner.
- The cause rules are thresholds, listed in each script's docstring.
- The reference has no decoys, so some "not in the reference" sequence may be decoy sequence.
- Nothing here is in spike yet.

## Re-running

```bash
REF=GRCh38_no_alt.fasta SLICE=hg002_chr20_slice.bam Q100_VCF=GRCh38_HG2-T2TQ100-V1.1.vcf.gz bash clips/run.sh
REF=GRCh38_no_alt.fasta SLICE=hg002_chr20_slice.bam SNV_BAM=snv/merged.bam bash qmodel/run.sh
```

- `SNV_BAM` is the merged BAM of the 25-SNV run in `docs/presentation/scripts/runs.sh`.
- `qmodel/e4m.py` and `e9m.py` read `clips/nonerr_mask.pkl`.
- `e8.py` takes about 25 minutes.
