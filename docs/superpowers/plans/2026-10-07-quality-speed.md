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
