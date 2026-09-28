# Transplant test, round 3: fake-vs-real mismatch against counting noise at each site

Locked 2026-09-28, before its code and before any split-half number exists,
with the user's agreement ("1": round 3, no spike change). Branch `transplant`.
No new spike runs: this round re-reads round 2b's reads (`c00c643`, raw output
in the session scratchpad, `transplant2/full/`).

## Why a round 3

Round 2b supported 13 of 16 group-by-direction verdicts. DUP50-299 reverse was
refuted on spread, and INS20-49 was inconclusive in both directions.

After that result, a look outside the locked rule (seen, not judged) found this:
- At shared sites, two real samples differed *less* than a Poisson (J) or
  binomial (A) noise model allows: 0.63 for DUP J, 0.56 for INS20-49 A.
- Two independent samples cannot do that if the model is right. A review by
  another model (GPT 6 Astra) pointed this out.
- So the noise model is wrong for these metrics. On top of that, the yardstick
  c comes from other sites than d, sites with different depth and counts.

Round 3 measures the counting noise at each site itself, from its own reads,
and compares like with like.

## Gate A

1. **Principle.** spike's version of a variant should differ from the real one
   no more than two real samples differ at a variant they share, once each
   difference is scaled by the counting noise at its own site.
2. **What would kill it (for spike):** in a group, the fake-vs-real mismatch
   over noise is more than 2.25 times the real-vs-real one, at the lower end of
   its interval (the rule is below).
3. **Already refuted here?** Round 2b refuted DUP50-299 reverse on raw spread,
   with c taken from other sites. This round asks whether that holds once each
   site's own noise is measured. It is a new question, not a rerun.
4. **Simplest version.** Split each site's read pairs into two random halves.
   The gap between the halves measures counting noise there. There is no new
   spike run and no realigning, and the metrics are round 2's
   (`small_evidence.py` is unchanged).
5. **Inputs, and how each is checked:** see "Checks that stop the run".

## Method

**Reads.** The same reads as round 2b, as `small_evidence.py` reads them:
- real: the donor BAM;
- fake: the recipient BAM, minus the run's `replaced_reads.txt`, plus the
  run's `sim.bam`;
- shared: each sample's own BAM.

Each event's run is found the way `round2.py` made it (`transplant.batches`
over the arm's events, in `events.tsv` order). spike refused 8 event-arm pairs
(RF8, `refused.tsv`); those pairs have no fake reads and are left out.

**Split.** A read pair goes to half `blake2b(f"{k}:{name}")[0] & 1` for salt k.
Mates share a name, so they stay together.

**Per event and sample** (x is fake or real, or HG001 or HG002):
- x is the metric on all reads, as in round 2b.
- For each of K = 10 salts, x0 and x1 are the metric on each half.
- v = the mean over salts of (x0 - x1)^2 / 4.
  - A salt where either half's metric is undefined (no reads to divide by) is
    skipped.
  - An event with no usable salt, for either sample, is left out of round 3.
    The number left out is reported per group.
  - Why /4: a half has half the read pairs. For a ratio of counts, a half's
    variance is about twice the full sample's, and the difference of two
    independent halves has twice that again.

**Per set, arm, group and metric:**
- R_f = sum of d^2 over sum of (v_fake + v_real), over transplanted events,
  with d = fake - real.
- R_c = sum of c^2 over sum of (v_HG001 + v_HG002), over shared events.
- rho = R_f / R_c.
- An interval for rho comes from 2,000 bootstrap draws. Transplanted and shared
  events are resampled apart (Python `random.Random(3)`). The interval is the
  5th and 95th percentiles (`transplant.percentile`). A draw whose sum of v is
  0 is skipped.

d^2 holds bias and spread together, so one number judges both.

**Per metric:**
- **pass** if rho's 95th percentile is at most 2.25;
- **fail** if its 5th percentile is above 2.25;
- **unsure** otherwise.

2.25 is not a new number. It is 1.5^2: rounds 1b and 2 allowed d's 10th-90th
width up to 1.5 times c's, and this is the same allowance on a variance scale.

**Broken controls:** B1 (VAF 0.25), and B2 (200 bp off) for J, as round 2.
- A metric counts only when each control that must break it **fails** under
  this rule. A control that passes, or is unsure, makes the metric
  inconclusive. That is stricter than round 2, where a control only had to not
  pass.
- Reverse uses the forward controls, as before.

**Per group and direction:**
- supported if every metric passes;
- refuted if any counted metric fails;
- inconclusive otherwise.

## Checks that stop the run

1. **Same reads as round 2b.** Every event's all-reads metric must equal its
   value in round 2b's `evidence.tsv`, to its 4 written decimals (an empty cell
   for an undefined value). Any mismatch stops the run, with no verdicts.
2. **The split-half noise is right where the answer is known.** For SNVs,
   binomial noise A(1-A)/N, with N the reads counted, is right. Measured in
   the look after round 2b: the shared-site ratio was 1.01.
   - Over every real SNV measurement (forward real, reverse real, and both
     shared samples), the sum of v over the sum of A(1-A)/N must lie within
     0.8-1.25.
   - Outside that range the estimator is not trusted, and the run stops with no
     verdicts.
3. **Each event is in the run it is said to be in.** An event's position must
   lie inside one of its run's `events.bed` intervals. Any miss stops the run.

## What each outcome would mean (written before the run)

- **DUP50-299 reverse passes:** round 2b's refutation came from comparing
  different sites without their own noise. spike's per-site mismatch is within
  natural differences.
- **It fails:** spike misses something at these sites, beyond counting noise
  and beyond natural differences between two people. Next come spike changes:
  read lengths, fragment placement, and phase and indels from a phased VCF.
- **It is unsure:** 54 transplanted duplications per direction and 33 shared
  ones are too few to tell.
- **INS20-49:** whether B1 now fails decides whether the group can be judged.

Also reported beside the verdicts (seen, not judged): R_f and R_c per group and
direction, the events left out, and round 2b's verdicts next to round 3's.

## The phased-neighbour check: not run

It was chosen with this round. It was measured before locking, and it cannot be
done well with the truth sets on disk:
- Each sample's small-variant region leaves out the area around its own
  structural variants.
- Of the 142 duplication sites (54 forward, 55 reverse, 33 shared), **0** have
  their window (the copy +-150 bp) fully inside the common small-variant
  region.
- The carrier's own region covers a median of 0.60 of the window (HG001,
  forward) and 0.72 (HG002, reverse). So calls next to the duplications are not
  benchmark-trusted (`scratchpad round3/precount_bed.py`).
- For small variants the check is empty by design: no other record lies within
  150 bp of a drawn event (round 2's isolation rule).

It can come back with calls trusted around structural variants.

## Order

1. This plan (`plan:` commit).
2. Code: `scripts/transplant/round3.py` with tests written first, and a
   mutation run (`__pycache__` cleared, `python -B`, for every mutant).
3. The run, at most 16 processes. Its log is kept, its exit status is checked,
   and its last line must be there.
4. `result:` commit.

## Not in this round

- Changes to spike: read lengths, fragment placement, phase, indels.
- GPT 6 Astra's 2x2 design (needs complete local haplotypes).
- Deeper data from the same libraries.
- The seed-2 replant from the dig. It changed spike's copy choice too, so it is
  not pure read noise.
