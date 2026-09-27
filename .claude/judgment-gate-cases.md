# Judgment-gate case file for `spike`

This repo's own history of retractions: the specific ways this codebase has fooled its author. A
general gate catches general errors; this catches the ones that actually happen here. Each entry says
what was claimed, what the claim actually rested on, and which gate would have caught it earlier.

## 2026-09-26 — T3: "mappability dominates, 5 of 6" was 3 of 6

- **Claimed.** Of CR2's six depth-fold warnings on real chr20 duplications, five "would not have
  warned had the fold been counted at any MAPQ". A design note recommending an any-MAPQ fold was
  written on that basis.
- **What it actually rested on.** A proxy measured in the wrong units. `SIM_DEPTH_FOLD` is a ratio of
  **fragment** coverage over **proper pairs whose both mates pass `--min-mapq`**
  (`simulate::estimate_coverage_at` counts `pool.pairs`), over spike's own per-segment bins. The proxy
  `fold_any` was a ratio of `samtools depth` **read-base** coverage with **no proper-pair
  requirement**, over a uniform grid. The locked plan then declared the two equivalent — "1.5 is not a
  new number: it is `census::DEPTH_FOLD_WARN_ABOVE` … So 'mappability' means exactly *this bin would
  not have warned at any MAPQ*" — and every downstream sentence used that equivalence.
- **How it was caught.** The final whole-branch review ran the counterfactual spike's own way:
  `depth_fold` is computed from the pool before anything is planted, so `spike --min-mapq 0` *is* the
  fold at any MAPQ, in spike's units, over spike's bins. Folds became 1.62, 1.35, 1.70, 1.48, 1.20,
  1.69 — **three still warn**. The plan's own verdict rule (4 of 6) then gives *inconclusive*, whose
  outcome rule is "propose nothing".
- **The tell that was in the data all along.** Every anchor in T3's own table has
  `pool_depth / any_depth` between 1.31 and 1.40 — the insert-to-read-length ratio. If one estimator
  were a subset of the other that ratio could not exceed 1. It was printed in the result table and not
  read.
- **Which gate would have caught it.** Gate A step 5 (*what must be true of the inputs, and how was
  each verified*). The plan did record that `scaled_by` is a fragment depth — and then reasoned past
  it: "no measurement here is compared against a read depth across estimators — every comparison is a
  ratio of two windows under one estimator." True of each ratio, false of the **threshold**, which was
  borrowed from the other estimator. **Rule for this repo: when a plan borrows a threshold from
  production code, the plan must show the measured quantity is in the same units as the quantity that
  threshold was set on — not merely that the measurement is internally consistent.**
- **Cheaper check that existed and was not used.** The direct counterfactual was three spike runs with
  one flag changed. It was never considered because a bespoke estimator felt more precise.

## 2026-09-25 — T3's first census run: a silent tool failure summed to zero

- **Claimed.** Nothing, yet — the run produced `any_read_depth = 0.0` for every window.
- **What it rested on.** `--ff` is `samtools view`'s spelling of the flag filter; `samtools depth` has
  no such option. It printed usage to stderr, exited non-zero, and wrote nothing to stdout, which the
  caller summed to zero.
- **How it was caught.** The pool depths printed beside them read 55x. A depth of 0.0 next to 55x is
  impossible, and the control column is what made it visible.
- **Which gate would have caught it.** Gate D question 2 (*does my control distinguish working from
  broken?*). **Rule for this repo: every shell-out in a measurement script raises on a non-zero exit.
  A helper that returns stdout unchecked will eventually return silence and be summed.**

## 2026-09-26 -- RF8: a noise test on random spots said nothing about real SV sites

- **What the plan set out to show.** RF8's kill test K was to show that refusing above
  `SIM_RESIST` 0.5 would not stop ordinary runs. It drew 40 random 10 kb spots inside the
  HG002 benchmark on chr1, and **0 of 40** were above 0.5.
- **What it rested on.** Random placements. spike's users do not spike random spots. They
  spike known SV sites, and so does the repo's own `validate_pipeline.sh`. Those sites sit in
  repeats. On the pipeline's background, **6 of its 20** real HG002 deletions are above 0.5,
  and the new default stops its spike step outright.
- **How it was caught.** Before the docs step, every repo script that runs spike was checked
  for events the new rule would refuse. It was not a locked criterion.
- **Which gate would have caught it.** Gate A step 5 (*what must be true of the inputs*).
  **Rule for this repo: a noise test for a new default is drawn from what the rule will
  actually see. For spike, that is known SV sites, and the repo's own harnesses first. A
  random placement is a control, not the population.**

## 2026-09-26 -- a count that read 0 under mawk and 6 under gawk

- **What almost shipped.** `validate_pipeline.sh` was to log how many events have
  `SIM_RESIST` above 0.5, using `awk ... substr(...) + 0 > 0.5`.
- **What it rested on.** That `awk` reads `0.776` as a number. This machine's `LC_NUMERIC` is
  `sv_SE.UTF-8`, which has a decimal comma. `mawk` honours it, so every value became 0 and the
  count was 0, with exit 0. `gawk`, the default `awk` here, ignores it and printed 6.
- **How it was caught.** The script's own comment claimed mawk support, so the count was run
  under both. The two disagreed.
- **Which gate would have caught it.** Gate D question 2: a check must be shown to go red, or
  here, to be non-zero where it must be. **Rule for this repo: a script that parses decimals
  with `awk` runs it as `LC_ALL=C awk`. This machine's locale uses a decimal comma, and `mawk`
  turns every such number into 0 silently.**

## 2026-09-26 -- a locale bug predicted in `highest_vaf` that was not there

- **Claimed.** After the entry above, docs/review/REVIEW.md filed `highest_vaf` in `validate_pipeline.sh`
  as having the same trap, right today only because 0.5 comes first in `--vafs`. RF12's
  `--min-recall` comment said the same of its compare: "both sides become 0".
- **What it rested on.** Pattern-matching on the entry above, not a run. Both compares pass
  their decimals with `awk -v`. `mawk 1.3.4` reads `-v` values in the C locale. Only fields
  read from input, and `printf` output, go through `LC_NUMERIC`.
- **How it was caught.** The fix was to start with a failing test. Run under `mawk` with
  `LC_NUMERIC=sv_SE.UTF-8`, the function was right in all 7 orders tried. A mutant that feeds
  the same values as input picked 0.1 from `0.1 0.25 0.5`, so the harness could go red.
- **Which gate would have caught it.** Gate D question 3 (*did I measure it, or infer it?*).
  A follow-up in docs/review/REVIEW.md is a durable claim. **Rule for this repo: a locale claim says how
  the number gets into awk. Under mawk, input fields and `printf` output follow `LC_NUMERIC`;
  `-v` values do not. `LC_ALL=C awk` stays the rule, because it is right either way.**

## 2026-09-26 -- RF6: a junction probe built from the reference, on reads that carry the sample

- **What the plan set out to show.** A 31-base probe `ref[start-15, start) + ref[end, end+16)`
  would find every correct deletion's junction in its reads. The plan was locked, then the kill
  test ran.
- **What it rested on.** That the reads spike writes spell the *reference* on both sides of
  the junction. They do not. spike puts the sample's own SNPs on the event copy (het alleles
  of the chosen haplotype, and every hom-alt), so wherever one falls inside the probe, **no**
  correct read matches it exactly. That was 2 of 40 chr20 deletions. Both pass `split_reads`
  today, so the row would have *added* two false failures.
- **How it was caught.** The kill test K1 (at least 2 carriers on 40 of 40), run before any
  code. That is the gate working.
- **Which gate would have caught it earlier.** Gate A step 5: *what must be true of the
  inputs.* **Rule for this repo: any check that compares read bases to the reference must
  allow for the sample's own alleles, because spike writes them onto the simulated reads on
  purpose. This includes SNPs inside a probe, and hom-alt SNPs on both copies.**

## 2026-09-26 -- RF11: "a small hole in the check" was a rule that fails both ways

- **Claimed.** After one 40 bp insertion failed `ins_reads` because bwa-mem2 clipped all its
  reads, RF11 was described to the user as a small hole: count clips below 50 bp too.
- **What it rested on.** The false-FAIL side alone. The false-PASS side of the same rule had never
  been measured. N10 set the 50 bp clip floor from 5 empty sites at a 3 bp threshold. Nobody had
  counted how often an empty site passes at 50.
- **How it was caught.** The locked plan's null test measured today's rule beside each proposed
  floor. Today's rule passes a 50 bp insertion that is not there at 12 of 200 random and 46 of 200
  simple-repeat sites. That was confirmed with `spike validate` itself on 400 of 400 sites. Every
  lower floor passes more. Refuted before any code.
- **Which gate would have caught it earlier.** Gate B step 2 (*run the trivial baseline beside
  it*). **Rule for this repo: before changing what a check counts as evidence, measure the current
  rule's pass rate where the event is absent, at the threshold in use. A false FAIL is visible; a
  false PASS is not, and only a null test finds it.**

## 2026-09-26 -- RF13: a mismatch tolerance larger than the part of the probe that decides

- **What the plan set out to show.** A row counting spike's own reads that carry a 31-base
  junction probe, within 2 substitutions, would pass correct insertions and fail wrong ones.
  It was locked, then the kill test ran.
- **What it rested on.** That 2 substitutions are small next to the difference between a
  right and a wrong truth record. The tolerance was sized for the sample's SNPs and
  sequencing errors across the whole probe. For a 4 bp insertion, the part that tells right
  from wrong is 4 bases, and a wrong one often differs in only 2. N2 matched on 3 of 6.
- **How it was caught.** The locked K- test with wrong letters, at the shortest length in the
  set. The smoke test had been on a 40 bp insertion only.
- **Which gate would have caught it earlier.** Gate A step 5 (*what must be true of the
  inputs*), over the whole input range. **Rule for this repo: a tolerance in a matching rule
  is checked against the smallest discriminating part of the input, not the average one.
  Name the shortest case, and put it in the smoke test.**

## 2026-09-26 -- RF14: a kill test that forced spike past its own refusal

- **What the plan set out to show.** At least 41 of 42 correct deletions would have 1 or more
  of spike's reads carrying the join. It was locked, then the kill test ran.
- **What it rested on.** That every K positive was a properly planted event. To copy
  `validate_pipeline.sh`, the plan passed `--allow-resistant` on every run. That flag exists
  to override RF8's refusal of events spike cannot plant properly (`SIM_RESIST` above 0.5).
  One "clean" random site, chr20:7119236, is a low-MAPQ region at 0.93 to 0.95. spike tiled
  36 pairs over its ~4 kb haplotype, and two runs had no read across the join.
- **How it was caught.** The locked K+ (at most 1 miss). At `SIM_RESIST` 0.5 or below it was
  34 of 34.
- **Which gate would have caught it earlier.** Gate A step 5 (*what must be true of the
  inputs*). The site draw defined "clean" by the benchmark and by distance from indels, not by
  whether spike can edit the reads there. **Rule for this repo: when a kill test turns off one
  of spike's own safety defaults, the plan says what the row must do on the events that
  default exists to stop, or judges them apart. `SIM_RESIST` is only known after spike runs,
  so it cannot be filtered at draw time without running spike.**

## 2026-09-27 -- origin: "the 35x HG002 BAM has no XA tags", from one window

- **What was claimed.** The 35x HG002 BAM carries no `XA` tags, so `--edit-model origin`
  cannot run on it. It was written into the spec, the plan (`84ca4c8`), a README draft and
  the notes to the user.
- **What it rested on.** One 3 kb window at chr20:7118000-7121000, where 0 of 1070 primary
  records carry `XA` and 618 are MAPQ 0. That spot is a repeat whose MAPQ 0 reads have more
  than 5 hits, so bwa-mem leaves `XA` off by its `-h 5` rule. The file keeps the tag: its
  first 100,000 records hold 15,255 with `XA`, and chr20:10-11 Mb holds 2,841.
- **How it was caught.** A review of the plan found a MAPQ 40 read with `XA` in Task 9's
  footprint. The counts above were then measured over the whole file's start and a 1 Mb
  window.
- **What else it hid.** The spec's per-spot input check ("MAPQ 0 reads and none with `XA`:
  stop") would have stopped `origin` at exactly the repeats it exists for. The check is now
  on the file.
- **Which gate would have caught it earlier.** Gate D3 (*did I measure it, or infer it?*). A
  count from one window was stated as a fact about the file. **Rule for this repo: a claim
  about a whole file is measured on the whole file, or on more than one region of it, and
  says which. An absence at one spot is a fact about that spot.**

## 2026-09-27 -- origin: "the 35x HG002 BAM holds chr20 only", so no look-alike read was removed

- **What was claimed.** The final review of `--edit-model origin` said the 35x HG002 BAM
  holds reads on chr20 only ("idxstats: chr20 22,206,212, everything else 0"). So `origin`
  removed 0 fragments at every look-alike on real data, and the real GIAB test had to use
  whole-genome BAMs to reach that half. It went into the run's PLAN-DEFECTS (PD-27,
  escalated), the notes to the user and the plan for the GIAB test.
- **What it rested on.** An idxstats line whose command and file were not written down. The
  file does not match it: `samtools idxstats` on
  `HG002.novaseq.pcr-free.35x.bwamem2.dedup.grch38_no_alt.bam` shows reads on 195 contigs
  (chr1: 81,466,062). The 14 look-alikes of `del:chr20:7119236-7120236` hold 3518 primary
  reads.
- **What was actually going on.** `origin` still removes 0 fragments there on the whole
  genome. Of 52 look-alike fragments with a read whose `XA` hit lies in the footprint, 25 have
  a mate without one (R4 keeps them), 18 have a mate outside the regions read, and 9 pass
  (awk over `samtools view`, `XA` start positions only). A whole-genome BAM is needed for the
  look-alike half, but not enough: at this spot R4 is what keeps it small.
- **How it was caught.** The warning written for PD-27 (a look-alike region that holds no
  read) did not fire on that BAM. idxstats was then run on it.
- **Which gate would have caught it earlier.** Gate D3 (*did I measure it, or infer it?*). A
  zero was explained ("the reads are not in the file") without measuring the explanation, and
  the evidence was quoted without its command. **Rule for this repo: a claim about what a file
  holds names the file and the command beside the number, and a zero's cause is measured
  before it is written down.**
