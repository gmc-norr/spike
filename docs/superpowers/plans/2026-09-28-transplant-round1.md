# Transplant test, round 1: do spiked deletions look like real ones?

Locked 2026-09-28, before any code, with the user's agreement. Branch `transplant`, off master `fe15a46`.

## Principle

A deletion spike plants should leave the same read evidence as the same
deletion really present in a comparable sample, sequenced the same way.

## Samples (checked 2026-09-28)

| | HG001 (NA12878) | HG002 (GM24385) |
|---|---|---|
| BAM | `Seq25-7598` 30x resample | `Seq25-7600` 30x resample (SeraCare inherited-cancer mix) |
| Run | LH00352:45:227NC2LT1, lane 2 | the same |
| Aligner | bwa-mem2 2.2.1 `mem -M -K 100000000`, MarkDuplicates | the same |
| Reference | spike's GRCh38 no_alt (195 contigs, same names and lengths) | the same |
| SV truth | Platinum Pedigree v1.2 (`NA12878_hq_v1.2.svs`) | T2T Q100 v1.1 |

Both are in `~/dev/cnv_validation/from_rv/raredisease_seracare_na12878_resample/`.

## Where

Q100's `stvar` benchmark region AND Platinum's `svs` region, autosomes only,
minus the 7 SeraCare genes (BRCA1, BRCA2, MSH2, MSH6, MLH1, PMS2, CDKN2A)
+-2 kb and +-2 kb around the cell line's MSH2-intron <-> IGL junction ends:
2,543 Mb (`scripts/transplant/count.sh`, committed with the code).

## Why aligning only spike's reads is enough (measured 2026-09-28)

Each region's own read pairs from the HG001 BAM were realigned alone with the
real command and compared read by read with the BAM:
- chr20:15-16 Mb (ordinary): 238 of 247,452 records (0.096%) moved, proper-pair
  bit 0.002%, MAPQ 0.004%.
- chr20:7.10-7.14 Mb (a repeat): 890 of 10,166 (8.8%) moved, 888 of them MAPQ 0
  before and after, all out of the region -- a re-drawn pick among equally
  good copies, one-way because the other copies' reads were not realigned.
- Pinning the insert size with `-I` changed nothing that matters.

So spike's reads are aligned alone, and the sham below measures what that
does at each event.

## Events: deletions >= 50 bp (measured: truvari bench with its default matching, `--sizemin 50 --passonly`; genotype shown is HG001's, HG002's genotype of the shared ones not yet read)

| Size | HG001 only, het | HG002 only, het | both (het in HG001) |
|---|---|---|---|
| 50-299 bp | 1,809 | 3,372 | 1,621 |
| 300-999 bp | 478 | 692 | 540 |
| 1-10 kb | 234 | 267 | 193 |
| 10 kb+ | 19 | 31 | 7 (report only) |

- **Forward:** HG001-only het deletions, spiked into HG002; real = HG001.
- **Reverse:** HG002-only het deletions, spiked into HG001; real = HG002.
- **Control (real vs real):** deletions het in both; HG001's reads vs HG002's.
- Excluded before drawing: any other truth SV (>= 50 bp) of either sample
  within 1 kb; and, for transplants, any locus where the recipient's own
  BAM already shows evidence of the deletion (more than 1 evidence read, N1).
- Drawn: up to 100 per size bin and set, seed 1. Counts after the
  exclusions: PREDICTED (not run) >= 100 except 1-10 kb.

## Spiking

- spike master `fe15a46`, `--edit-model clean`, `--seed 1`, VAF 0.5,
  without `--allow-resistant`: an event spike refuses is counted as refused
  per bin, not compared.
- Events of a set go in a few spike runs, >= 100 kb apart within a run.
- Aligned with the real command (custom `--aligner`):
  `bwa-mem2 mem -M -K 100000000 -t 16 -R '@RG\tID:sim\tPL:ILLUMINA\tSM:<sample>'`.
- No whole-BAM merge: evidence is read from the recipient BAM, skipping the
  names in `replaced_reads.txt`, plus `sim.bam`.
- **Sham** per spike run: the same replaced reads, realigned unchanged with
  the same command, read the same way. It shows what realigning alone does.

## Evidence (a separate script, not `spike validate`)

Per event, over primary, non-duplicate, non-QC-fail reads, divided by the
depth of the 1 kb on each side:
- **J (junction evidence):** read pairs showing the deletion in any way -- a
  CIGAR gap of >= 50% of its length starting within 20 bp of its start, a
  soft clip >= 10 bp within 10 bp of either end, a split alignment (SA)
  joining the two ends within 20 bp, or a discordant FR pair across it with
  insert above bwa's proper-pair bound for that BAM (878 for HG001, from
  chr1:50-60 Mb by bwa's own rule; HG002's measured the same way first).
- **E1 (depth ratio):** mean depth over the deleted bases / flank depth
  (deletions >= 300 bp).

## Pass rule (per size bin, per direction, for J and for E1)

- d = fake minus real, per transplanted event; c = HG002 minus HG001, per
  shared event.
- **Pass** if the median of d lies between the 25th and 75th percentiles of
  c, **and** the 10th-90th percentile width of d is at most 1.5x that of c.
- **Refuted** if either part fails with the broken controls failing too.

## The test must be able to fail (broken spikes)

Run on the same forward events, read the same way; each must FAIL the rule,
or the verdict is "inconclusive" (the test is blind):
- **B1:** VAF 0.25 instead of 0.5 (half the evidence).
- **B2:** each deletion planted 200 bp to the right of the truth, evidence
  read at the truth position (evidence in the wrong place).

## Order

1. **Pilot (can only stop us, not change the rule):** 10 events per set in
   the 300-999 bp bin, with B1. Stop and rethink the evidence if B1 passes
   or if c's 10th-90th width is wider than the median real J.
2. Full run, all bins.
3. `result:` commit: per bin and direction, supported / refuted /
   inconclusive, with refusals, sham and broken-control numbers.

## Not in this round

Insertions, duplications, inversions, small variants; `--edit-model
origin`; SV callers (round 2); PD-29's depth-fold threshold.
