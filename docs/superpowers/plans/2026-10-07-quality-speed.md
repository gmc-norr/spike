# spike's quality model v2 runs fast enough, and writes exactly what it writes today

**Asked 2026-10-07.** v2 (`7f28bb0`) passed K1, K2 and K7 but **failed its own speed limit**:
1.27x against the pre-v2 binary, where the plan allowed 1.25x
(`docs/superpowers/plans/2026-10-07-quality-model-v2.md`, Result). Two claims in that result were
also wrong, and are corrected here, reported but not gated.

`$RUN` below is the run directory this work was done from (set it to wherever the evidence and the
two reference binaries were kept); nothing in the repository depends on it.

## Found

**The cost.** `draw_classes` recomputes `run_ratio(m, h, c)` for every cell of the weight table,
for every pair it makes. `READ_CLASSES = 8` (`src/quality.rs:50`) and `RUN_BINS = 4`
(`src/quality.rs:77`), so that is **64 weights and 128 `run_ratio` calls per pair**, and each call
runs `class_given_run` and `class_share`, which both sum over the whole `class_run` table
(`2 · RUN_BINS · READ_CLASSES = 64` counts). The whole table has only **64 distinct cells**, and
`class_run` never changes after the profile is learned.

**Measured before this plan** (`dup:chr20:14550000-15550000` on the HG002 35x BAM, `--seed 1
--threads 8`, 3 runs each, alternating; `$RUN/evidence/speed.txt`, `speed_probe.txt`):

| binary | wall seconds | median | against pre-v2 |
|---|---|---|---|
| pre-v2 (`spike-premaster`) | 9.613 / 9.026 / 9.373 | 9.37 | 1.00x |
| v2 (`7f28bb0`) | 12.660 / 11.859 / 11.537 | 11.86 | **1.27x** |
| pre-v2, second sitting | 9.610 / 9.102 / 9.219 | 9.22 | 1.00x |
| a throwaway probe that caches the table | 10.505 / 10.389 / 10.185 | 10.39 | **1.13x** |

The probe's output was md5-identical to v2's on two events (`dup:chr20:14550000-15550000` plus
`snp:chr20:38550000:T:C;af=0.1`, `--seed 1 --threads 8`): R1/R2 FASTQ, `truth.vcf`,
`replaced_reads.txt`, 120,099 `SPIKE_` reads.

**The probe is not the design.** It falls back to `1.0` when its table is empty, which silently
gives wrong draws for any profile built without the fill.

## Gate A

1. **Principle.** The cached table is the same arithmetic, read instead of recomputed. Nothing
   about the model changes, so **nothing about the output may change**.
2. **What would kill it.** Any byte of spike's output differing from `master` `7f28bb0` (B1).
   That is the deciding check, and it decides against the change: if identity cannot be kept, the
   code is not committed.
3. **Not a model change.** No threshold, no distribution, no random draw is touched. K1, K2, K7,
   K6 and every other v2 check keep their verdicts *by* B1, not by being re-run.
4. **Simplest thing.** 64 f64s computed once, at the one point where `class_run` becomes final.

## Design (locked)

**The table.** `QualityProfile` gets a `run_ratios` field: `run_ratio(m, h, c)` for every
`(mate, run bin, class)`, `2 · RUN_BINS · READ_CLASSES = 64` f64s, flattened as
`(m · RUN_BINS + h) · READ_CLASSES + c`.

- **Filled once, where `class_run` becomes final.** `class_run` is written in exactly one place
  (`src/quality.rs:573`, inside `learn`) and read nowhere else but `class_share` and
  `class_given_run`. `learn` is the only constructor: `from_read_pairs`, `from_donor_pairs`,
  `from_blocks` and `from_input` all route through it. So `learn` is split: a private
  `learn_counts` keeps today's body and returns the profile with the counts, and `learn` calls it
  and then fills the table. `learn_counts` is private and called from nowhere else, so there is no
  path to a profile whose table was not filled — including the early return for an empty quality
  alphabet.
- **No silent fallback.** The field is `Option<Vec<f64>>`, `None` until filled, and the read
  `expect`s it. A missing table panics in release as well as in debug. There is no neutral `1.0`
  to fall back to.
- **Read by `draw_classes`, in today's order.** The weight stays
  `joint · rr0 · rr1 · r0 · r1`, the same five factors multiplied left to right, so each f64 is
  bit-identical to today's. The run-bin index is clamped the same way `class_given_run` clamps it
  (`h.min(RUN_BINS - 1)`).
- **Allowed second step, only if S1 fails with the table alone.** Also cache the 64 weights per
  `(run bin 1, run bin 2)` per generator, with the same product order. Whether it was needed is
  recorded in `$RUN/STATUS.md`.

Nothing else changes. No new crate; `Cargo.toml` and `Cargo.lock` stay as they are.

## Checks

### B1 — byte identity. The deciding check.

`$RUN/bin/spike-base` (md5 `961e2a3c2549db04ba5119caaaa7b8a8`, built from `7f28bb0`) against the
new binary. md5 of the **decompressed** `R1.fq.gz` and `R2.fq.gz`, and of `truth.vcf`,
`replaced_reads.txt` and, where it exists, `fastq_removed_reads.txt`, on:

- **B1a** 35x BAM: `--event dup:chr20:14550000-15550000 --event "snp:chr20:38550000:T:C;af=0.1" --seed 1 --threads 8`
- **B1b** 31-value BAM: the K6 events exactly as `$RUN/evidence/k6.sh` runs them
  (`--event "snp:chr20:38550000:T:C;af=0.1" --event del:chr20:38900000-38910000 --threads 8`, **no `--seed`**)
- **B1c** 35x BAM: `--event del:chr20:38900000-38910000 --seed 7 --threads 1`

**PASS only if every file matches in all three.**

**What would trivially pass this:** two empty files have equal md5s. So each compared `R1.fq.gz`
must also be checked **non-empty and holding `SPIKE_` reads**, and the `SPIKE_` count is recorded
beside each md5. A run that failed and wrote nothing is a FAIL, not a match.

### S1 — speed. Run once; the verdict is final.

`dup:chr20:14550000-15550000` on the 35x BAM, `--seed 1 --threads 8`, wall time,
`spike-premaster` and the new binary **alternating, 5 runs each**, `uptime` recorded before every
run. **PASS when median(new) ≤ 1.25 × median(premaster).** The machine is shared (load average
11.20 / 16.39 / 11.21 on 72 cores when this plan was written), so the loads are reported with the
times. It is not rerun for a better number.

### T — tests, clippy, toolchain

- Every test passes. The count is **747 passed, 0 failed, 3 ignored** (measured on `7f28bb0`,
  one test binary) **plus the tests added here**, with nothing removed. A `test result: FAILED`
  line, or no `test result` line at all, is a FAIL.
- `cargo clippy --all-targets` is no worse than the baseline: **16** lines matching
  `^(warning|error)` — bin "spike" 12 warnings, bin "spike" test 14 (12 duplicates), 0 errors.
- `cargo +1.82 check --all-targets` is clean. `7f28bb0` already is (`Finished dev profile ... in
  8.84s`), so this is a real gate, and it forbids `is_multiple_of` (Rust 1.87).

### M — mutants. Each must turn at least one test **FAILED**.

A mutant counts as caught only when a `test result: FAILED` line names a failing test. An
`error[E…]` with no `test result` line is a compile error, not a caught mutant. The 3 ignored
tests are green by being skipped and never count.

- **M1** the table filled with `1.0` everywhere.
- **M2** the table read with the mates swapped — mate 1's row used for mate 0.
- **M3** the run-bin index clamped to 0 when the table is read.
  - The brief's M3 was "the table filled before the class/run counts are learned". With all-zero
    `class_run`, `class_given_run = (0 + RUN_PRIOR_READS · class_share) / (0 + RUN_PRIOR_READS) =
    class_share`, so every ratio is exactly `1.0`: that mutant is **bit-identical to M1**, not a
    second mutant. The clamp is locked instead, as the brief allows. The early-fill variant is
    still run, and both are reported.

If a mutant survives, that is a finding: the test that catches it is added, and both the survival
and the new test go in this plan's Result and in `$RUN/STATUS.md`.

**Tests, written first and seen red before the code exists:**

1. **The table is the formula, cell for cell.** On a learned, asymmetric profile, every one of the
   64 cells equals `run_ratio(m, h, c)` recomputed, bit for bit. The test also asserts the table
   has teeth — not `1.0` everywhere, mate 0's row differing from mate 1's, and run bin 3's row
   differing from run bin 0's — so it can go red for M1 and M3.
2. **`draw_classes` picks what the formula picks.** On the same profile, for 200 fixed seeds and
   for the asymmetric run pairs `(3, 0)` and `(0, 3)`, the pair `draw_classes` returns equals the
   pair an oracle returns, where the oracle rebuilds the weights from `run_ratio` directly with
   the same product order and consumes one `f64` from an identically seeded RNG. This is the only
   test that can see M2.

## Corrections (reported, not gated)

Two claims in v2's committed Result are wrong. `$RUN/evidence/k6_runbins.py` and `mates.py` are
re-run, their output pasted into the Result, and both scripts kept under
`docs/analysis/quality-model-v2/`:

- **K6's "a 31-value alphabet crashes too often" (24.3% against 14.4%) is withdrawn.** It compared
  spike's remade reads, which sit over the event's own template, with every read in a 20 kb window.
- **K3's "mates crash together less often" (4.74 against 11.6) is downgraded.** It rested on 5
  both-crashed pairs.

Neither is a model change, and neither is gated. v2's `2026-10-07-quality-model-v2.md` gets a
`## Correction` section appended at its end; its Result section is left exactly as committed.

## Result (2026-10-07)

**Done as planned, and the speed limit is met: 1.18x against 1.25x, with every byte of the output
unchanged.** The allowed second step (caching the 64 weights per run-bin pair) was **not needed**.

Code: `src/quality.rs` only, +172 / −2 lines (`693fe83`). `learn` is split into `learn` and a
private `learn_counts`; `learn` fills `run_ratios`, a `Option<Vec<f64>>` of
`2 · RUN_BINS · READ_CLASSES = 64` cells; `run_ratio_at` reads it and `expect`s it, so a missing
table panics in release too; `draw_classes` keeps the product order `joint · rr0 · rr1 · r0 · r1`.

### B1 — byte identity. **PASS.**

`spike-base` (md5 `961e2a3c2549db04ba5119caaaa7b8a8`, from `7f28bb0`) against the new binary
(`df547601763d864dad84be90eaf6b69d`). Five files per case; the md5s below are **both** binaries'.

| case | R1.fq.gz | R2.fq.gz | truth.vcf | replaced_reads.txt | fastq_removed_reads.txt | R1 content |
|---|---|---|---|---|---|---|
| **B1a** 35x, `dup:chr20:14550000-15550000` + `snp:chr20:38550000:T:C;af=0.1`, `--seed 1 --threads 8` | `0bcd4052fb4e9935ab92d0a3b5a42503` | `27e7b36a4bbbb648b3b757c82f5e3e36` | `3c9571b20cc99a3356e46515376ae891` | `36db9382788438eda41c33c7cad14904` | `71ff0263712f6b670219c105387a7a5c` | 758,396 lines, **120,099 `SPIKE_`** |
| **B1b** 31-value, the K6 events, **no `--seed`**, `--threads 8` | `fd005167f50e988227a139c804dde2de` | `2a563a9898bee79a225f54156c4e8ca8` | `43ecf37fa7d6a53145801104c080b75b` | `06b4805f713adc21583077cea59842b4` | `33397ca6851175fecb59f4532de765a5` | 26,428 lines, **325 `SPIKE_`** |
| **B1c** 35x, `del:chr20:38900000-38910000`, `--seed 7 --threads 1` | `df3e10b2262bdc00da5e23a7d4fcd046` | `7e329da95e8920b3ddfb6f503ed52101` | `8ef6a9249b8f710dd4722fe18cb2da6e` | `d4308d74e6d0d56d6978f410fadb4275` | `03a2b6fd09f26609c262de802f6755f3` | 15,512 lines, **303 `SPIKE_`** |

**The gate was run against input it must reject.** With `spike-premaster` in the new binary's
place: `B1a FAIL: md5 differ`, `B1b FAIL: md5 differ`, `B1c FAIL: md5 differ`, `B1 FAIL`. So the
PASS above is not two empty files matching.

### S1 — speed. **PASS, 1.18x.** Run once.

| run | premaster (s) | load before | new (s) | load before |
|---|---|---|---|---|
| 1 | 8.554 | 2.48 5.30 8.37 | 10.084 | 3.03 5.33 8.34 |
| 2 | 8.162 | 3.76 5.40 8.34 | 9.778 | 3.78 5.38 8.31 |
| 3 | 8.250 | 3.96 5.36 8.28 | 10.012 | 3.98 5.32 8.23 |
| 4 | 9.121 | 3.99 5.27 8.18 | 10.237 | 4.42 5.33 8.17 |
| 5 | 9.106 | 4.05 5.22 8.10 | 10.197 | 4.12 5.20 8.06 |

Median premaster **8.554 s**, median new **10.084 s**, ratio **1.1789 → 1.18x**. The ceiling was
1.25 × 8.554 = 10.693 s. v2 was 1.27x. The machine was quieter than when v2 was timed (load 2.5-4.4
against 11.2), so premaster is faster here too; the ratio is what the check compares, and both
binaries ran alternately under the same load.

### T — tests, clippy, toolchain. **PASS.**

- `test result: ok. 749 passed; 0 failed; 3 ignored; 0 measured; 0 filtered out` — the baseline's
  **747** plus the **2** tests added here, nothing removed.
- `cargo clippy --all-targets`: **16** lines matching `^(warning|error)`, and `diff` of the two
  lists against the baseline prints nothing: identical, warning for warning.
- `cargo +1.82 check --all-targets`: `Finished dev profile [unoptimized + debuginfo] target(s) in
  2.45s`, clean. `grep -rn is_multiple_of src/` finds none.

The two tests, both seen red before the code existed:

| step | state | what the suite said |
|---|---|---|
| RED 1 | nothing | `7 × error[E0599]: no method named run_ratio_at` and **no `test result` line** — a compile error, not yet a caught red |
| RED 1b | first `asymmetric_profile` | `only 4 distinct pairs drawn: too few to see a wrong weight` — six quality values collapse the class cuts; the helper was widened |
| RED 2 | `run_ratio_at` stubbed to 1.0 | `test result: FAILED. 748 passed; 1 failed` — `cell (mate 0, run bin 0, class 0): table 1 against the formula 0.038461538461538464` |
| RED 3 | the stub, read by `draw_classes` | `test result: FAILED. 746 passed; 3 failed` — the table test, the oracle test at `runs (3, 0), seed 0`, and the existing `synth::tests::test_reads_crash_after_a_run_of_one_base_as_the_donors_do` |
| GREEN | the real table | `test result: ok. 749 passed; 0 failed; 3 ignored` |

### M — mutants. **3 of 3 caught.** None survived.

| mutant | the edit | suite | the test that went red |
|---|---|---|---|
| **M1** table all 1.0 | `fill_run_ratios`: `.map(\|(_m, _h, _c)\| 1.0)` | `FAILED. 746 passed; 3 failed` | `test_the_run_ratio_table_is_the_formula_cell_for_cell` (`table 1 against the formula 0.038461538461538464`), `test_draw_classes_picks_what_the_run_ratio_formula_picks`, and `synth::tests::test_reads_crash_after_a_run_of_one_base_as_the_donors_do` |
| **M2** mates swapped at the read | `draw_classes`: `run_ratio_at(1, runs.0, c1)` and `run_ratio_at(0, runs.1, c2)` | `FAILED. 748 passed; 1 failed` | **only** `test_draw_classes_picks_what_the_run_ratio_formula_picks`. The oracle test is the single check that sees this, which is why it exists |
| **M3** run bin clamped to 0 | `run_ratio_at`: `(m * RUN_BINS + 0)` | `FAILED. 746 passed; 3 failed` | the table test, the oracle test, and the same `synth` run test |

Two further variants, reported not gated:

- **The brief's own M3** — the table filled before the counts are learned — prints
  `table 1 against the formula 0.038461538461538464`, M1's message character for character, and
  reddens the same three tests. With `class_run` all zero,
  `class_given_run = (0 + 8 · class_share) / (0 + 8) = class_share`, so every ratio is exactly
  `1.0`: it **is** M1. That is why the clamp was locked as M3.
- **The mate swapped inside `run_ratio_at`** instead of at the call site:
  `FAILED. 747 passed; 2 failed`, `table 0.056338028169014086 against the formula
  0.038461538461538464` — caught by both new tests.

Every mutant was restored by copying a saved clean `src/quality.rs` back; its md5 after each
restore is `108d951eaa7f8f8b7216343142b0b0d9`, equal to the saved copy. No `git checkout --` was
used on uncommitted work.

### Corrections, re-run here (reported, not gated)

Both scripts now live in `docs/analysis/quality-model-v2/`.

`python3 docs/analysis/quality-model-v2/k6_runbins.py <v2's K6 sim.bam> <the 31-value BAM> <reference>`:

```
crashed share by run bin 0 / 1 / 2 / 3 (reads in brackets)
donors in windows    R1:   4.1% (2449)    6.2% ( 403)    9.0% ( 543)   58.6% ( 263)
donors in windows    R2:  14.4% (2500)   13.0% ( 377)   18.9% ( 545)   77.7% ( 224)
spike remade reads   R1:   6.0% ( 184)    3.1% (  32)    7.1% (  28)   61.7% (  81)
spike remade reads   R2:  14.6% ( 185)   14.8% (  27)   23.1% (  39)   73.0% (  74)
donors in windows    overall 13.88% of 7304; run bin 3 holds 6.7%
spike remade reads   overall 24.31% of 650; run bin 3 holds 23.8%
expected for spike's reads at the donors' per-bin rates: 23.76%
```

`python3 docs/analysis/quality-model-v2/mates.py <prefix>…` on the K2b held-out 35x pairs and on
the K6b held-out 31-value pairs:

```
k2b/real:    pairs 62682  P(R1) 0.84%  P(R2) 1.27%  P(both) 0.083%  ratio 7.77
k2b/spike21: pairs 62682  P(R1) 0.74%  P(R2) 1.12%  P(both) 0.088%  ratio 10.55
k2b/spike22: pairs 62682  P(R1) 0.73%  P(R2) 1.11%  P(both) 0.089%  ratio 11.02
k2b/spike23: pairs 62682  P(R1) 0.74%  P(R2) 1.10%  P(both) 0.097%  ratio 11.91
k6b/real:    pairs 40232  P(R1) 6.78%  P(R2) 18.55%  P(both) 3.482%  ratio 2.77
k6b/spike:   pairs 40232  P(R1) 6.15%  P(R2) 16.43%  P(both) 2.960%  ratio 2.93
```

Because **B1b passed**, the new binary writes v2's K6 reads byte for byte: the `sim.bam` these
numbers were measured on is the new binary's output as much as v2's, so nothing above needs
re-running against this commit.

`docs/superpowers/plans/2026-10-07-quality-model-v2.md` gets a `## Correction` section at its end;
its Result section is left exactly as committed. The README's "How close" line drops the 31-value
claim.

## Whole-branch review (2026-10-07)

One review agent on the most capable model, in the foreground, against this plan and trying to break
B1. **Verdict: SOUND WITH FINDINGS.** It rebuilt the branch (`cargo build --release` gives
`df547601763d864dad84be90eaf6b69d`, equal to the binary B1 was run with), rebuilt the `7f28bb0`
baseline independently from `git archive` (`747 passed; 0 failed; 3 ignored`), re-ran every gate, ran
**12 byte-identity cases this plan never tried** — an insertion, `--edit-model origin`, a fusion with
`--exon-bed`, the 31-value BAM with `--indel-error-rate`, `--dup-model junction`, `--flank 2000`,
`--min-mapq 40` (a different profile), an indel at `af=0.02`, `--into-fastq` (7 files), `--vcf` with
`--allow-overlap`, `--threads 16`, `--region` — **all byte-identical**, and re-ran all three mutants
plus four of its own. It could not break B1, and closed the argument structurally: one write site for
`class_run`, one caller of `learn_counts`, no `Clone`/`Default`/serde, no struct literal outside
`learn_counts`, `runs` provably ≤ 3, and no fast-math in the profile, so no reassociation.

### Finding 2, Important — fixed in `e399a24`

Nothing in the suite protected this branch's purpose. Delegating `run_ratio_at` back to `run_ratio`
throws the whole speed win away and **all 749 tests stayed green**, with B1 and clippy unchanged;
only the one-off S1 timing would have noticed.

`test_draw_classes_reads_the_table_rather_than_recomputing_it` pokes
`(mate 0, run bin 3, class 5)` to `1e12` and asserts mate 0 then draws class 5 every time, and that
it did not already do so before the poke. Against the delegation it is the **only** test that
reddens: `test result: FAILED. 749 passed; 1 failed`, `poking (mate 0, run bin 3, class 5) to 1e12
left mate 0 drawing {5, 7, 4, 3}: draw_classes is not reading the table`.

Mutants with it present: **M1** 747/3, **M2** 748/2, **M3** 746/4, **M4** (the delegation) 749/1.
Tests now **750 passed, 0 failed, 3 ignored**; clippy still **16**, identical to the baseline;
`cargo +1.82 check` clean; the release binary's md5 **unchanged**, so B1 still holds.

Writing the poke index as `(0 * RUN_BINS + 3) * READ_CLASSES + c` was a clippy **error** ("this
operation will always return zero"), taking clippy from 16 to 18 with a hard error. Spelled through
a closure instead.

### Finding 1, Critical — fixed

The correction's own attribution was wrong. It named **one** T13 "about 200 bp left of the deletion
edge" as the cause of the bin-3 enrichment. Reproduced independently with `bin3_runs.py`, now in
`docs/analysis/quality-model-v2/`:

```
bin-3 spike reads: 155
     38 reads  runs in span: (('T', 12, 38910539),)
     38 reads  runs in span: (('T', 13, 38911171),)
     32 reads  runs in span: (('T', 15, 38911534),)
     30 reads  runs in span: (('T', 13, 38899797),)
     17 reads  runs in span: (('A', 17, 38549533),)
read-start clusters:
    38549303-38549529 (17)
    38899580-38899788 (30)
    38910304-38910529 (38)
    38910934-38911534 (70)
```

The named T13 (chr20:38,899,798, 202 bp left of the left edge — that part was right) carries **30 of
155**, i.e. **4.6% of the 650**, not 23.8%. 108 sit over a T12, a T13 and a T15 0.5-1.5 kb to the
*right* of the deletion, and 17 over an A17 in the SNV window. The withdrawal itself is unaffected
and stronger — five long runs under a narrow footprint against 6.7% of a 25 kb window — but naming
one run without counting the reads under it is the same class of mistake as the claim being
withdrawn, so it is written up as such in `.claude/judgment-gate-cases.md` too.

### Finding 4, Minor — fixed

The K6b held-out set is one seed with no committed recipe, where K2b spans three. The Correction now
says so and no longer leans on it as a size.

### Finding 3, Minor — declined, with the reason

The v2 plan's Result still states the withdrawn K6 number and still lists it under "Still off", with
the withdrawal 50+ lines below. The review suggested a bracketed "[withdrawn]" pointer per bullet.
**The brief for this work forbids editing that Result section — it is the committed record** — and a
bracketed insert is an edit. The `## Correction` section opens by saying the Result is left as it
stands and what is withdrawn from it, and the commit subject of this result says it too. Left for the
human to decide.
