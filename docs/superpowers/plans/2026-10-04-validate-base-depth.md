# `spike validate`'s coverage rows count sequenced bases, not read spans

**Asked 2026-10-04.** This is review finding 5 (P2, `docs/review/2026-10-04-independent-review.md`), the user's last pick of the five. We reproduced it on master `2764dbc`:
- the reviewer built a homozygous 50 bp deletion with exact alignments, at 0 aligned-base depth on every deleted base (`samtools depth -aa`);
- `coverage_ratio` observed **1.02** against an expected 0.00, and failed.

## Read in the code

- `count_depth_in_region` (`src/validate.rs`) is documented as "Mean read depth (aligned bases per base)".
- It gives each read a single span, `[start, start + ref_span)`. `ref_span` adds up M, D, N, = and X, so a read with a `D` across a deletion counts as covering every base it deletes.
- Its only users are `check_coverage_ratio`'s two rows, `coverage_ratio` and the advisory `coverage_any_mapq`. `ref_span` has no other caller.
- Which events it touches:
  - `coverage_ratio` applies only to structural DEL and DUP records. Small indels get `allele_freq`.
  - A read crosses a deletion with a `D` only when the deletion is shorter than the read, so the deletions that move are short structural ones, of about 50-150 bp.
  - The 20 deletions of `validate_pipeline.sh` (the RF12 runs in `scratchpad/rf12/pipe/`) run from 505 to 33,028 bp.

## Design (locked)

`count_depth_in_region` keeps one interval per aligned block, the M, = and X ops, instead of one per read. D and N advance the reference without adding an interval. `mean_depth` is unchanged; it sums the intervals' overlap with the region. The function then means what its doc says, and what `samtools depth` counts without `-J`. Nothing else changes: not the windows, the flanks, the expected ratios or the tolerance.

**Tests, written first and seen red:**
- a read `10M50D10M` over a region inside its `D` adds 0 to the region's depth, and 10 bases to each side's;
- a read `5S20M` adds 20 bases, at its alignment start;
- `N` adds nothing; `=` and `X` count like `M`.

**Mutation checks** (each must turn a test red; the unmutated tests are green first):
1. D counted again;
2. N counted;
3. `=`/`X` not counted;
4. the soft clip moving the reference position.

## Checks (locked)

**V1: the reviewer's probe.** `scripts/review_20261004.py` (the scratch copy from the earlier results) against the new debug binary: the `true_del` `coverage_ratio` row observes 0.00 and passes. Every other result is as before.

**V2: nothing moves where nothing should.** The three RF12 runs (`spike_vaf_{0.1,0.25,0.5}`, 20 deletions each) are validated by the new and the master binaries with the pipeline's settings.
- Every `coverage_ratio` and `coverage_any_mapq` row's observed value moves by at most 0.02, and no verdict flips.
- Every other row is identical.

**V3: what it is for, on real reads.** On the hospital BAM: 8 structural deletions, `del:chr20:P-(P+L);af=0.5`, with P = 12, 15, 20, 35, 40, 45, 50 and 55 Mb, and L = 60, 68, ..., 116 bp. They are run with `--seed 1 --threads 16 --align`.
- If spike refuses some of them (RF8), those are dropped once and the run is repeated.
- A regional merged BAM: for each event, the original reads over P ± 6,000 bp minus `replaced_reads.txt`, plus `sim.bam`, sorted and indexed.
- Both binaries validate it.

Pass if, for every event, the new `coverage_ratio` equals an independent ratio within 0.02. That ratio is `samtools depth -a -Q 20 -G UNMAP,SECONDARY,SUPPLEMENTARY,QCFAIL,DUP` (no `-J`): the mean over the event `[POS, END)` divided by the mean of the two 5,000 bp flanks.
- The master binary's values and both verdicts are reported.
- **Control:** the same samtools ratio with `-J`, which counts deletions, must differ from the new value by more than 0.1 on at least 6 of the 8. Otherwise the reads do not cross these deletions with a `D`, and V3 tests nothing.

## Result (2026-10-04): V1 and V2 pass; V3 fails its control, so it is inconclusive

Code `a7486d1`. Binaries: new release `93eb71b1`, new debug `4f8c7118`, master release `91f8fbab`. Scratch: `scratchpad/vbd/` (`v3.sh`, `v3.log`).

**Tests.** 699 pass (2 new), 0 fail, 1 ignored. The mutation runner checks the unmutated `validate::` tests are green first (137), then catches **4 of 4** mutants. Mutant 1 (D counted) is caught only by the new test: no earlier test saw the bug. clippy gives the same warnings as master.

**V1: PASS.** The reviewer's `true_del` probe now gives `coverage_ratio` expected 0.00, observed 0.00, PASS, and `coverage_any_mapq` the same. The run still exits 1, on `split_reads`, `split_reads_each_end` and `dup_rate`, which the review set aside as separate. Every other result is as in the indel-error-rate result.

**V2: PASS.** Over the three RF12 runs' 120 coverage rows (60 `coverage_ratio`, 60 `coverage_any_mapq`):
- the largest move is 0.01, with 12 rows moving by 0.01;
- no verdict flips;
- every other row is identical.

**V3: FAIL on the control, so V3 is inconclusive.** All 8 deletions ran, and none was refused.

| deletion | samtools, no `-J` | samtools, `-J` | new | master |
|---|---|---|---|---|
| 60 bp | 0.379 | 0.469 | 0.38 PASS | 0.47 PASS |
| 68 bp | 0.289 | 0.289 | 0.29 PASS | 0.29 PASS |
| 76 bp | 0.467 | 0.467 | 0.47 PASS | 0.47 PASS |
| 84 bp | 0.228 | 0.227 | 0.23 PASS | 0.23 PASS |
| 92 bp | 0.680 | 0.680 | 0.68 PASS | 0.68 PASS |
| 100 bp | 0.594 | 0.594 | 0.59 PASS | 0.59 PASS |
| 108 bp | 0.422 | 0.422 | 0.42 PASS | 0.42 PASS |
| 116 bp | 0.586 | 0.586 | 0.59 PASS | 0.59 PASS |

- The new values equal samtools without `-J` at 8 of 8.
- The control needed `-J` to differ by more than 0.1 at 6 or more. It did at **0 of 8**; the 60 bp deletion came nearest, at 0.09.

**Why, counted rather than inferred.** These are the counted reads (flags `0xF04` off, MAPQ 20 or more) overlapping each deleted span in the regional merged BAM.
- At 60 bp, 3 of spike's 5 reads there carry the deletion as a `D`.
- At 68-116 bp, none does. Of spike's reads there, the ones at 100-116 bp are soft-clipped instead (8 of 8, 18 of 20, and 23 of 24).
- No original read carries a `D` at any of them.

So the plan's statement that "the deletions that move are short structural ones, of about 50-150 bp" was a prediction, and it was wrong. bwa-mem2 writes a `D` here only up to about 60 bp. On real reads the fix moves only deletions about that short. The reviewer's 50 bp case is one of them, and so is the 60 bp case above, which moved from 0.47 to 0.38.

**What this leaves.**
- The fix is right where reads cross a deletion with a `D`: V1 on exact alignments, and the 60 bp row on real reads, both equal to samtools.
- It changes nothing that should not change (V2).
- How large its effect is on real data was not shown beyond that one row.
