# Transplant test, round 2: do spiked SNVs, small indels and duplications look like real ones?

Locked 2026-09-28, before any code, with the user's agreement ("1").
Branch `transplant`. Round 1b (`99dde72`)
supported deletions >= 50 bp. This round keeps its method and pass rule
(`2026-09-28-transplant-round1.md`, as changed by `...-round1b.md`) and
applies them to the variant types the user will spike at the hospital
(LDLR and similar genes): SNVs, 1-49 bp indels, and duplications.

## Gate A

1. **Principle** (unchanged): a variant spike plants should leave the same
   read evidence as the same variant really present in a comparable sample,
   sequenced the same way.
2. **What would kill it:** in a group, the median fake-minus-real difference
   falls outside the real-vs-real middle half, or is more spread, while the
   half-dose control fails. RF15 predicts one way this can happen: spike's
   20-39 bp insertions came out cleaner than HG002's own. That was measured at
   random sites against real ones at their own sites. Same-site transplants
   separate the site from the content.
3. **Already refuted here?** No. The method was supported in round 1b. RF15
   is open, and this round measures it.
4. **Simplest version:** one metric for SNVs (allele fraction). For indels, two:
   exact-form allele fraction, and evidence in any form. The second is there
   because of RF15's failing real insertions: of 60, 23 had no I in any read
   and 25 had fewer than two full-length I's. The exact count alone cannot
   tell that from a lower dose. Duplications of 50-299 bp get J only, like round 1's smallest bin.
5. **Inputs, and how each was checked:** see "Inputs" below.

## Samples, spiking, sham

As round 1: HG001 `Seq25-7598` and HG002 `Seq25-7600` 30x BAMs (same run and
lane), bwa-mem2 `mem -M -K 100000000`. spike master `fe15a46` (binary md5
`0302a168`), `--edit-model clean`, `--seed 1`, VAF 0.5, no
`--allow-resistant`; events >= 100 kb apart within a run. Aligned with the
real command, evidence read from the recipient BAM minus
`replaced_reads.txt` plus `sim.bam`, and a sham per run.

Event specs: `snp:chr:POS:REF:ALT;af=0.5` for SNVs and indels (POS 1-based,
as the truth VCF), `dup:chr:s-e;af=0.5` for duplications (s 0-based, default
`--dup-model full`). RF8 refusals are read from spike's labels,
`SNV  chr:POS REF>ALT` and `DUP  chr:s+1-e (Lbp)`, each with two spaces
after the type (`src/main.rs` `event_label`).

## Truth sets and region

| | HG001 | HG002 |
|---|---|---|
| Small variants | Platinum Pedigree v1.2.1 `latest-smallvar` | T2T Q100 v1.1 (the round 1 VCF) |
| Small-variant region | Platinum v1.2 `smallvar.bed` | GIAB `HG002_GRCh38_v5.0q_smvar.benchmark.bed` |
| Duplications | Platinum v1.2 `svs` | T2T Q100 v1.1 |

- **Small variants:** the two regions intersected, autosomes, minus the
  round-1 masks: 2,514 Mb. Both VCFs are split and left-normalised
  (`bcftools norm -m-any -c x`), and records are matched on CHROM, POS, REF
  and ALT.
- **Duplications:** round 1's SV region (2,543 Mb), with insertions matched by
  truvari bench (default matching, `--sizemin 50 --passonly`). A duplication
  is an insertion >= 50 bp whose inserted bases exactly equal the reference
  right after POS or ending at POS. Neither truth set writes any record as
  `SVTYPE=DUP`, so this is the only way to find them.

## Events (measured 2026-09-28)

Het only. Isolated: no other record of either sample within 150 bp (small),
or no other SV >= 50 bp within 1 kb (duplications, as round 1).

| Group | HG001 only | HG002 only | both |
|---|---|---|---|
| SNV | 397,275 | 413,723 | 284,638 |
| DEL 1-4 | 28,936 | 30,925 | 22,356 |
| DEL 5-19 | 3,674 | 3,836 | 2,515 |
| DEL 20-49 | 489 | 527 | 318 |
| INS 1-4 | 28,150 | 30,715 | 21,720 |
| INS 5-19 | 2,716 | 2,957 | 2,135 |
| INS 20-49 | 449 | 489 | 311 |
| DUP 50-299 | 271 | 404 | 207 (het in HG001) |

The DUP row is before the 1 kb isolation filter; after it, PREDICTED (not
run): >= 100 in each set. Larger exact duplications are too few to test:
13 / 13 / 5 at 300-999 bp and 0 / 2 / 0 at 1 kb and up.

- Forward: HG001-only into HG002, real = HG001. Reverse: HG002-only into
  HG001, real = HG002. Control c: variants in both, recipient minus donor.
- Before drawing a transplant event, the recipient's own BAM must show at most
  1 read carrying it (N1), counted as A's carriers below (J's reads for DUP).
- Drawn: up to 100 per group and set. The pilot draws with seed 1. The full
  run draws with seed 2 and cannot draw the pilot's events (as in round 1b).
- Shared events are het in both samples. A shared event's isolation ignores
  its own record in each sample.

## Evidence (a new script beside `evidence.py`)

Reads: primary, mapped, not duplicate, not QC-fail (flags 0xF04 excluded),
any MAPQ, as round 1.

- **SNV, A:** reads with ALT at POS divided by reads with any base aligned at
  POS.
- **Indel repeat region:** the deleted or inserted bases as a unit, extended
  along the reference in both directions for as long as it repeats (as
  `spike validate` does). Outside a repeat it is the event itself.
- **Indel, A (exact form):** carriers / (carriers + reference reads).
  - A carrier has an indel op of the truth's type and exact length, starting
    in the repeat region +-1 bp. For an insertion, the op must also give the
    truth's sequence there.
  - A reference read has aligned bases over the repeat region +-10 bp, with no
    indel and no clip in that span.
- **Indel, E (any form):** read pairs (by name) showing the indel in any
  way, divided by the mean depth of the 1 kb on each side of the repeat
  region. A read counts if it is a carrier; or has
  an indel of the same type, >= 50% of the length, starting within the repeat
  region +-10 bp; or has a soft clip >= 5 bp whose clip point lies within the
  repeat region +-10 bp.
- **DUP, J:** read pairs (by name) showing the duplication, divided by the
  mean depth of the 1 kb on each side of the segment [s, e). A pair counts
  for any of:
  - an insertion op >= 50% of the segment's length, within 20 bp of it;
  - a soft clip >= 10 bp within 10 bp of s or e;
  - a split alignment (SA) with one part ending within 20 bp of e and the
    other starting within 20 bp of s.

## Pass rule and verdicts (as round 1b)

- d = fake minus real per transplanted event; c = recipient minus donor per
  shared event.
- **Pass:** median d lies within c's 25th-75th percentiles, and d's 10th-90th
  width is at most 1.5x c's.
- Broken controls, on the forward events:
  - **B1** (VAF 0.25) must fail every metric in every group.
  - **B2** (planted 200 bp to the right) must fail J for DUP.
  - A metric whose must-fail control passes is inconclusive.
  - Reverse uses the forward controls.
- Per group and direction: **supported** if every metric passes; **refuted**
  if any fails; **inconclusive** otherwise.
- Refusals (RF8) and the sham are reported beside the verdicts.

## Order

1. **Pilot (can only stop us, not change the rule):** 10 events per set in
   SNV, INS 20-49 and DUP 50-299, with B1 (and B2 for DUP). Stop and rethink
   if either happens:
   - B1 passes any of those metrics;
   - spike fails on an event for any reason other than RF8.
2. Full run, all groups.
3. `result:` commit.

## Inputs, and how each was checked

- **spike takes these specs and plants them in place.** Smoke run on the
  HG002 BAM (chr20, one of each type), with 4 of 4 accepted at the requested
  positions. In `sim.bam`:
  - the SNV G>A had ALT in 13 of 25 reads;
  - the 8 bp deletion had 10 reads with an 8 bp D;
  - the 20 bp insertion had 8 reads with a 20 bp I;
  - the 100 bp duplication had no 100 bp I, but 24 reads clipped by >= 10 bp.
- **Both truth tables are complete.** The first count normalised HG002 only
  up to chr9. `bcftools norm` stopped at a REF mismatch at chr9:70701801, with
  its stderr sent to /dev/null. It now drops such records (`-c x`, 1 dropped)
  and keeps its log, and both tables hold all 22 autosomes (HG001 4,200,000
  records, HG002 6,001,476). The counts above are from the complete tables.
- **The HG002 small-variant bed fits the Q100 v1.1 VCF.** Not verified from a
  header. In chr20:10-20 Mb of the common region (9.81 Mb), the four sets hold
  similar numbers of records:

  | Set | Records |
  |---|---|
  | Platinum | 17,199 |
  | Q100 v1.1 | 16,709 |
  | GIAB v4.2.1 HG001 | 15,654 |
  | GIAB v4.2.1 HG002 | 15,519 |

  The N1 recipient check guards the draw against a variant the recipient
  carries but its VCF misses.
- **Known bias:** N1 drops sites where the recipient already shows 2 or more
  reads with the same indel, such as homopolymer stutter. The control c is not
  filtered that way.

## Not in this round

- Base quality, MAPQ and strand of the carrying reads
- MNVs and complex alleles
- Duplications of 300 bp and up (too few real ones)
- Non-tandem insertions of 50 bp and up
- `--edit-model origin`
- Callers

## Pilot (2026-09-28): STOPPED by its own check

Run with `d867cd6` (`2d68cb2`, plus a fix: the first try crashed sorting an SV
and a record with the same span, before anything was spiked), spike master
`fe15a46` (md5 `0302a168`); raw output in the session scratchpad,
`transplant2/pilot/`. 10 events per set in SNV, INS20-49 and DUP50-299.

**Pools after every filter** (`pools.tsv`): each group has >= 100 in each
set. DUP50-299 has 174 / 220 / 107, which confirms the plan's unrun
prediction. The small groups are within 1% of the plan's counts; the SV
neighbour filter removes the rest.

**The check stops the run:** B1 passes two of the pilot's metrics.

| Group | Metric | B1 median d | c 25th..75th | normal median d |
|---|---|---|---|---|
| SNV | A | -0.217 (fails, as it must) | -0.050..0.035 | 0.007 |
| INS20-49 | A | 0.000 (**passes**) | -0.042..0.133 | -0.008 |
| INS20-49 | E | -0.106 (fails, as it must) | -0.030..0.120 | -0.008 |
| DUP50-299 | J | 0.000 (**passes**) | -0.095..0.000 | 0.000 |

B2 fails DUP J on spread (width 0.417 against 1.5 x 0.244). spike refused 1
B2 duplication (RF8).

**Why, from the pilot's own rows (seen, not judged):**

- **INS20-49, A.** bwa-mem2 writes few 20-49 bp insertions as an exact I; it
  clips instead (RF11). Per event the exact carriers are 0 to 17 reads, with
  0 in 7 of 20 real events. At 10 events, halving a count that small does
  not move the median. E, which counts any form, sees B1.
- **DUP50-299, J.** J looks at the edges of one copy. When the copy lies
  inside a longer tandem array, the reads show the change at the array's
  edges, or nowhere, so J is about 0 in the real sample too. Array length
  over copy length, with real J, for the 20 transplanted events:
  - 1.00-1.05 (the copy is unique DNA): J 0.38-1.12 (7 events);
  - 1.14-1.46: 0.00-0.24 (6 events);
  - 1.47-5.46: 0.00-0.03 (7 events).
- Over the whole pools, copies in unique DNA (array < 1.10 x copy) number 61
  forward, 65 reverse and 35 shared.

Also seen, not judged: reverse INS20-49 E, fake minus real median +0.165.
Twice the real HG002 site showed 0 exact carriers where spike's reads show
6 and 10 (chr10:124378145, chr19:19863596). That is RF15's direction, on
10 events.

The full run is not started.
