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
