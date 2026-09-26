> Copied unchanged from the detached run's `/home/parlar_ai/spike-codex-run/NEW-FINDINGS.md` (2026-09-25). The `STATUS.md`, `TASKS.md` and `scratch/` it names are in that folder, outside git.
>
> **Status since the copy** (the findings below are left as the run wrote them):
>
> | ID | Status |
> | --- | --- |
> | NF1 | Open. Changing `MIN_EVENTS` or the region guidance is a default change. |
> | NF2 | Fixed, `dc73bf1`: the adjacent probes pass `--allow-overlap`; the script runs to exit 0 with all 14 rows. |
> | NF3 | Fixed, `395dcc7`: the last pid-only temp name in `truth.rs` is unique per test. |
> | NF4 | Written down, `39b1aa1`: the README's Build section says one target dir per commit, and md5 the binaries. A practice, not code. |
> | NF5 | Not fixed; an advisory row added beside it, `084104f`: `split_reads_each_end` asks for two joining reads at *each* breakpoint, and failed 6 of 40 real deletions against the pooled row's 4. The pooled `split_reads` check itself is unchanged — its row and the run's exit status moved on none of the 40 — and the new row is blind when the breakpoints are 500 bp apart or less. |
> | NF6 | Fixed, `d35f31b`: a type with no expected ratio gets no verdict (`N/A`, not a pass). |
> | NF7 | Fixed in part, `2844aff`: the harness records every tool's version in `<outdir>/tool_versions.tsv`. It still pins none. |
> | RF1 | Open. Out of that run's scope; the fix is `#[cfg(debug_assertions)]` on the three tests, or a refusal checked in both profiles. |
> | RF2 | Fixed, this docs commit: both README round trips re-measured on this branch, `13/13 PASS` and `12/12 PASS`, both exit 0. The `5/6` and `3/6` "before the check existed" figures are historical and stay. |
> | RF3 | Open. The shape of `summary` was not in T1's locked plan, and changing it is a default output change for existing `--json` parsers. |
> | RF4 | Open. The `<=` comparison is what T1's plan locked, and no spike version writes a negative value. |
> | RF5 | Open. T1's plan allowed the help text one new flag line; the footer sentences were added, extending the checks table was not ratified. README.md covers both rows. |
> | RF6 | Open. The cause is in the aligner's representation, not in the check's arithmetic, and changing the default verdict is exactly the default change this run's standing choice forbids. |
> | RF7 | Open. Adding the worst bin to `truth.vcf` for every event changes what spike emits; a per-event column in the run README would be the cheap fix, and it is a design decision rather than a defect. |
> | RF8 | Open, and a decision for the human, recorded with its numbers. A coverage floor would make runs that succeed today start failing, and `--min-donor-depth` would be a new flag. |
>
> RF1–RF8 are from the run after the CR4 and CR2 census (2026-09-26) and are written up at the end
> of this file, under their own heading.

# NEW-FINDINGS — found during this run, not fixed

## NF1 — `scripts/validate_pipeline.sh` can abort on a small `--region`, independent of this run

- **Where:** `scripts/validate_pipeline.sh:619-628` (`MIN_EVENTS=5`, line 99).
- **What:** the harness aborts when fewer than `MIN_EVENTS` truth DELs survive clustering. On
  a 1 Mb `--region` there are too few to start with.
- **How I know (measured):** het DELs of 500-50000 bp in
  `data/giab_hg38/HG002/GRCh38_HG2-T2TQ100-V1.1_chr20.vcf.gz`, clustered by the script's own
  rule: whole chr20 gives 247 raw -> 156 (old gap) -> 121 (new 7000 bp gap);
  `chr20:38000000-39000000` gives 6 raw -> **3** under both gaps. 3 < `MIN_EVENTS=5`.
- **Severity:** Low. The abort is loud, names the cause after task 2's fix, and `--min-events`
  lowers it. It is not caused by task 2 — the small-region count was already 3 before.
- **Caveat:** this count skips the harness's benchmark-intersection step, so it is an upper
  bound on survivors; the real numbers are lower.
- **Not fixed** in this run: changing `MIN_EVENTS` or the region guidance is a default change.

## NF2 — `scripts/review_sv_model.py` cannot run to completion on a branch where CR1 is fixed

- **Where:** `scripts/review_sv_model.py:139` (the `del_0.5_adjacent` / `del_1_adjacent`
  probes) with `run()`'s `check=True` at line 106.
- **What:** those two probes request `del:chrT:10000-11000` + `del:chrT:12000-13000` with no
  `--allow-overlap`. From commit `76986ca` spike rejects that — which is the CR1 fix working
  as intended. `run()` raises, and the script dies before it reaches the `ins_*`, `dup_*`,
  `del_lowmap`, `dup_homdel`, `vcf_hom` and `mate_recovery` rows.
- **How I know (measured):**
  `subprocess.CalledProcessError: Command '[... '--event', 'del:chrT:10000-11000;af=0.5', '--event', 'del:chrT:12000-13000;af=0.5']' returned non-zero exit status 1`
- **Severity:** Medium for anyone re-running the review's reproduction on a fixed branch.
- **Not fixed** in this run: the script was committed verbatim as the review's durable
  reproduction, and its sha256 (`0075194c06bd40a7…`) is the run's provenance check. A patched
  copy that adds `--allow-overlap` to multi-event probes lives at
  `/home/parlar_ai/spike-codex-run/scratch/review_sv_model_allowoverlap.py`; every "after"
  measurement in this run used it, and it reproduces the base binary's rows exactly.
- **Suggested fix for the human:** add `--allow-overlap` to the two `*_adjacent` probes in the
  repo's script, with a comment saying the probe deliberately measures the pre-CR1 symptom.

## NF3 — truth.rs's test temp-file pattern can produce a false red AND a false green

- **Where:** `src/truth.rs` tests, the `spike_truth_<pid>.vcf` naming pattern.
- **What:** two tests in the same process sharing one temp filename race. Task 4's implementer
  hit this for real while developing: `438 passed` on one run and `437 passed` on the next,
  from the same code.
- **How I know (measured, reported by the implementer and confirmed by the reviewer's seven
  consecutive suite runs after the fix):** the committed tests now use
  `spike_truth_ins_<pid>_<tag>.vcf` and are stable.
- **Severity:** Low now, latent. The pre-existing test in that file still uses the
  pid-only prefix; it is distinct today, so nothing is broken.
- **Not fixed** in this run: it touches tests unrelated to any CR finding.

## NF4 — a shared `CARGO_TARGET_DIR` silently returns a stale binary across different source trees

- **Where:** any before/after measurement in this repo that builds two commits.
- **What:** three different source trees built with one shared `CARGO_TARGET_DIR` produced
  three **byte-identical** binaries, with cargo printing `Finished ... in 0.06s`. The
  measurement would have shown a zero effect for a change that has a real one.
- **How I know (measured):** task 6's implementer caught it by md5 and redid the work with
  separate target dirs; task 6's reviewer independently confirmed the tainted binary's md5
  `968579338cfd43997e5fbd30d7789f7e` equals its own fresh build of `39361fe`.
- **Severity:** High for anyone doing before/after measurement here — it produces a confident
  wrong answer, not an error.
- **Not fixed** (it is a measurement practice, not repo code). **The rule:** give every commit
  its own `CARGO_TARGET_DIR` and md5 the binaries before trusting any comparison. This run also
  hit the milder version of the same hazard once — measuring with a `target/release/spike` that
  had not been rebuilt after `cargo test` (see task 1).

## NF5 — `check_split_reads` pools its evidence across both breakpoints

- **Where:** `src/validate.rs:650-671` (`names.extend` over both breakpoint windows) and the
  `MIN_SPLIT_READS` test at `:679`; the `expected` string is built at `:677`.
- **What:** the check requires two distinct read names *in total*, so both may sit at the same
  breakpoint with nothing seen at the other. Its own `expected` string, `">=2 joining chr:pos"`,
  reads as though each end must contribute.
- **How I know:** read, and independently confirmed by task 8's reviewer against those lines.
- **Severity:** Medium for anyone reading a `split_reads` PASS as evidence of a two-ended
  junction. It is documented now (`README.md`, commit `eac6f40`) but not fixed.
- **Not fixed** in this run: changing the threshold's meaning is a validation-model change, and
  CR9's overhaul is a Phase 3 design note.

## NF6 — `coverage_ratio_result`'s non-DEL/DUP arm is unreachable

- **Where:** `src/validate.rs:622` (a 0.5 tolerance arm) versus the dispatch at
  `src/validate.rs:138-147`, which only routes DEL and DUP to the coverage check.
- **What:** the arm's tolerance can never apply in production; its only callers are `:577` and
  unit tests at `:2655-2662`.
- **How I know:** read, and independently confirmed by task 8's reviewer.
- **Severity:** Low — dead-ish code, not a wrong answer. Task 8 deliberately left its tolerance
  out of the README so as not to document an unreachable number.
- **Not fixed** in this run.

## NF7 — `scripts/validate_pipeline.sh` pins no Truvari version and records none

- **Where:** `scripts/validate_pipeline.sh:75`, `TRUVARI="$(find_tool truvari "${TRUVARI:-}")"`,
  and the `truvari bench` invocations at `:855` and `:931`.
- **What:** the harness's matching criteria are whatever the installed Truvari makes of
  `--passonly -r 500 -p 0.5 -P 0.5 -s 500`. Truvari has renamed and re-defaulted these between
  major versions. The run records no version, so a result is not reproducible from its output.
- **How I know:** read. `truvari` is **not installed** in this environment, so I could not
  execute `truvari bench --help` to confirm any flag's current meaning — which is itself the
  point of the finding, and why `README.md` now states the flags rather than glossing them.
- **Severity:** Medium for reproducibility. It is also what made the review's own R9 gloss, and
  the README's first attempt at one, unverifiable.
- **Not fixed** in this run: pinning a tool version is a harness policy change.


## The run after the CR4 and CR2 census (2026-09-26)

### RF1 — `cargo test --release` fails three tests on an untouched tree

- **Where:** `src/truth.rs`, the tests
  `test_truth_refuses_a_depth_fold_slice_that_does_not_match_the_events`,
  `test_truth_refuses_a_resistant_slice_that_does_not_match_the_events` and
  `test_truth_refuses_an_adjusted_af_slice_that_does_not_match_the_events`.
- **What:** each expects a refusal that only fires in a debug build, so a release test run
  reports `473 passed; 3 failed` on code that is not broken.
- **How I know:** measured at `985e50f` with nothing changed. `cargo test --release` →
  `test result: FAILED. 473 passed; 3 failed; 1 ignored`; `cargo test` →
  `test result: ok. 476 passed; 0 failed; 1 ignored`.
- **Severity:** Low for correctness; Medium for anyone who runs the suite in release and reads
  the failure as a real defect.
- **Not fixed:** out of that run's scope. The fix is either `#[cfg(debug_assertions)]` on the
  three tests or a refusal checked in both profiles.

### RF2 — README's `spike validate` row counts are stale measurements

- **Where:** `README.md:797`, `:798`, `:804`, `:806`.
- **What:** each quotes a `<pass>/<total> PASS` from a round trip run before the census fields
  and the advisory rows existed. The `exit 0` / `exit 1` claims beside them all still hold --
  advisory rows stay out of the exit status -- but the totals are now two rows per event low.
- **How I know:** read, after T1 measured the same shape on a real slice: a one-event truth VCF
  that used to print `Result: 5/5 PASS` now prints `7/7`.
- **Severity:** Low. A reader comparing their own output against the README sees a mismatch
  that is not a defect.
- **Not fixed** in T1's fix pass: re-measuring needs a chr20 slice and a probe run, which was
  outside that pass.

### RF3 — `--json`'s `summary.fail` can be nonzero on an exit-0 run

- **Where:** `src/validate.rs`, `print_results_json`'s `summary` object.
- **What:** `summary.total` and `summary.fail` count every row, advisory included, so a run
  whose only failing row is advisory reports `fail: 1` and still exits 0. Nothing in `--json`
  records whether `--strict` was given, so a consumer cannot tell which rows the exit status
  counted.
- **How I know:** the T1 reviewer measured `fail: 7` on a run that exited on
  `5/5 validation checks failed`.
- **Severity:** Low. `scripts/validate_pipeline.sh` reads only `passed` and `len(checks)`, so
  its verdict is unaffected. A JSON consumer keying on `summary.fail` would now disagree with
  the exit status.
- **Not fixed:** the shape of `summary` was not in T1's locked plan, and changing it is a
  default output change for existing parsers.

### RF4 — a negative `SIM_RESIST` or `SIM_DEPTH_FOLD` passes its advisory row

- **Where:** `src/validate.rs`, the two census row builders.
- **What:** the rule is `value <= threshold`, so `SIM_RESIST=-1` passes. Nothing is recomputed
  from the BAM, so a nonsense number is believed as long as it parses.
- **How I know:** read, and flagged independently by T1's implementer and its reviewer.
- **Severity:** Low. No spike version writes a negative value; only a hand-edited truth VCF can
  produce one.
- **Not fixed:** the comparison is what T1's plan locked.

### RF5 — `spike validate --help` never says what `(advisory)` means

- **Where:** `src/validate.rs`, `print_usage()`'s "Checks, by truth-event type" block.
- **What:** the block lists the checks each event type gets but not the two advisory rows a user
  now sees in the output, and nothing in `--help` defines `(advisory)` beyond the footer
  sentence saying they are left out of the exit status.
- **How I know:** read; raised by T1's reviewer.
- **Severity:** Low. README.md covers both.
- **Not fixed:** T1's plan allowed the help text one new flag line; the footer sentences were
  added in the fix pass, and extending the checks table was not ratified.

### RF6 — `split_reads` fails on one correct real deletion in ten, by default

- **Where:** `src/validate.rs`, `check_split_reads`; the check is **not** advisory, so it
  decides the exit status.
- **What:** on 40 seeded 10-kb deletions inside the HG002 T2T-Q100 SV benchmark on chr20, each
  spiked into the 35x HG002 BAM at the default AF and run through `align.sh`, `merge.sh` and
  `spike validate`, the check reported `observed 0` on **4 of the 40** -- events 15
  (`chr20:25322805-25332805`), 25 (`chr20:49921688-49931688`), 33 (`chr20:56072849-56082849`) and
  34 (`chr20:56107004-56117004`). Each of those four runs exits 1 although nothing about the
  spike-in is wrong: the other checks pass, `coverage_ratio` reads 0.46-0.54 against an expected
  0.50, and the census rows are quiet.
- **How I know:** measured in T2's C4, and reproduced on **master's binary** for event 34 to rule
  out T2 as the cause: `split_reads >=2 joining chr20:561... 0 FAIL`, `Result: 4/5 PASS`, exit 1.
- **Severity:** Medium. `spike validate`'s exit status is what `scripts/validate_pipeline.sh` and
  any CI wrapper read, and a 10% false-failure rate on correct input teaches a user to ignore it.
  The check reads only the contig and position of an `SA:Z` entry, so a junction bwa-mem2 chose
  to represent some other way produces no evidence at all.
- **Not fixed:** the cause is in the aligner's representation, not in the check's arithmetic, and
  changing the check's verdict for the default case is exactly the default change this run's
  standing choice forbids. **It is also the baseline T5 must be judged against:** a
  per-breakpoint split-read row can only be stricter than a pooled one, so T5's bar has to allow
  for these four before anything about T5's own rate can be read.

### RF7 — the depth fold's worst bin is computed for every event and printed for almost none

- **Where:** `simulate::depth_fold` (`src/simulate.rs:925`) fills `DepthFold::worst_bin` and
  `worst_depth` for every event; `census::depth_fold_warning` (`src/census.rs:117`) returns `None`
  at or below `DEPTH_FOLD_WARN_ABOVE`, and the run README's "Donor depth off the scaling depth"
  line is filtered by the same threshold (`src/main.rs:1902`).
- **What:** only `SIM_DEPTH_FOLD`, the number, reaches the truth VCF. The bin it was measured over,
  and the depth in it, exist in memory for all 40 events of a run and are discarded for the 34 that
  do not warn. A reader who wants to know whether a 1.55-fold warning is unusual for that sample
  has nothing to compare it against.
- **How I know:** it blocked T3's own control. The locked control needed all 40 events' worst bins
  and could not have them, so the statistic had to be amended to a uniform grid the measurement
  computes for itself (plan amendment `597ad72`).
- **Severity:** Low for correctness, Medium for interpreting a warning. It also makes the
  measurement T3 did impossible to repeat from a run's own output.
- **Not fixed:** adding the bin to `truth.vcf` for every event is a change to what spike emits,
  which this run's standing choice forbids. The run README is not part of that promise (it "may
  gain lines"), so a per-event bin column there would be the cheap fix, and it is a design decision
  rather than a defect.

### RF8 — spike accepts an event whose donor depth at the breakpoint is 0.2x

- **Where:** `simulate::donor_coverage_for_tiling` (`src/simulate.rs:453`). Its `covered` test is
  `!cov.is_nan() && cov > 0.0`, so **any** non-zero coverage is enough; the tiling count is then
  `cov x VAF`.
- **What:** on `chr20:27100000-27110000` (99.69% of reads below MAPQ 20, 84x of any-MAPQ depth), the
  pool's depth at the breakpoint spike scaled by is **0.2x**. Spike accepts, plants `0.2 x VAF`
  fragments, exits **0**, and writes a truth VCF claiming a 10 kb deletion — while **5731 of the
  5749 reads over the event (99.7%) stay exactly as they were**. Both censuses warn
  (`SIM_RESIST=0.997`, `SIM_DEPTH_FOLD=15.44`, that fold being a 17.2x bin scaled by a 0.2x
  anchor), so the run is not silent; but nothing refuses it, and a user who reads only the exit
  status gets a truth VCF for an event that was barely simulated.
- **How I know:** measured in T4 on the 35x HG002 BAM with master's binary and `--seed 1`. Two more
  loci behave the same way: `chr20:26700000-26710000` scaled by 0.2x at `SIM_RESIST=0.998`, and
  `chr20:28200000-28210000` scaled by 0.7x at 0.995. The three windows one rank harder are refused
  outright for *zero* breakpoint coverage, so the accept/refuse boundary on this data sits between
  0.0x and 0.2x.
- **Severity:** Medium-High for anyone spiking into difficult regions. It is the same family as
  CR4 — the pool's filter decides what can be simulated — but sharper: here the pool is not merely
  incomplete, it is two orders of magnitude thinner than the locus, and the run still succeeds.
- **Not fixed:** a coverage floor would make runs that succeed today start failing, which this
  run's standing choice forbids. The conservative alternatives are a **warning** when the scaling
  depth is far below the event's any-MAPQ depth (which `SIM_DEPTH_FOLD` already half-says), or an
  opt-in `--min-donor-depth`, which would be a new flag and so out of this run. **It is a decision
  for the human**, recorded with the numbers above.
