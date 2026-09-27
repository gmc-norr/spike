# Threads, step 1: split the work that draws no random numbers

Locked 2026-09-27, before any code. Branch `threads`, off master `cb756ec`.

## Goal

Make big events faster on more than one core, with **byte-identical output**.
Only work that draws no random numbers is split. Tiling, the removal draw,
phasing's coin flips and the order of events stay on one thread with the one
random stream, so every `--seed` gives the same files as before (step 2 would
change that and is not part of this plan).

## Baseline, measured

`dup:chr20:14550000-17550000` (3 Mb) on the 35x HG002 BAM, `--seed 1
--allow-resistant`, master `cb756ec` plus timing log lines (binary `ecae674f`,
patch `timing.patch` in the session scratchpad). Time between log lines:

| Step | clean | origin | Random numbers? |
| --- | --- | --- | --- |
| Read extraction, pass 1 | 2.15 s | 2.16 s | no |
| origin: read the spot (892,170 records) | - | 2.83 s | no |
| origin: read 1,938 look-alike regions | - | 13.99 s | no |
| origin: `fragments()`, about six calls per site | - | ~6.5 s | no |
| Census | 1.13 s | 1.12 s | no |
| Quality profile | 3.64 s | 3.68 s | no |
| SNP pileup (count and collect passes) | 4.38 s | 4.35 s | no |
| Phasing | 2.49 s | 2.39 s | coin flips |
| Depth fold (origin: incl. building origin depth) | 0.35 s | 2.16 s | no |
| Tiling | 17.20 s | 15.32 s | **yes** |
| origin: removal chances | - | 2.76 s | no |
| origin: removal draw | - | 2.61 s | **yes** |
| Write R1 and R2 | 4.39 s | 4.32 s | no |
| **Total** | **39.88 s** | **67.82 s** | |

(`fragments()` is in the rows "origin: read the spot" to "editable built",
"origin depth at", "depth fold", "removal chances" and "origin removed"; ~6.5 s
is their sum less the other work in them, an estimate from the log.)

## What gets split

1. **A thread pool.** `rayon`, one global pool of `--threads` threads (default
   4, as now). `--threads` then sets spike's own threads as well as the
   scripts'; its help text and the README say so.
2. **Look-alike reads.** The regions are split into batches in list order;
   each batch opens its own reader. The per-region results are joined in
   region order and then go through today's loop unchanged (skips, empty
   regions, dedup by name and mate).
3. **`fragments()`.** Records are grouped by name with a stable parallel sort
   of their indices instead of a `BTreeMap`. Same fragments, same name order,
   same mate order within a name.
4. **Region reads split by where a read starts (BAM only).** A long region is
   cut into chunks; each record belongs to the chunk holding its alignment
   start (the first chunk also keeps records that start before the region).
   Joined in chunk order, that is the file's own order. Used for: the SNP count
   and collect passes, origin's spot read, extraction pass 1 and the census.
   A region shorter than one chunk is one chunk, so small events read the file
   as before. CRAM stays sequential.
5. **Quality profile.** Pairs are split into chunks, each fills its own bins,
   bins are joined and each is sorted. Sorted `u8` bins do not depend on order.
6. **FASTQ.** R1 and R2 are written by two threads, the same bytes each.
   Every pair is still checked before either file is created, and both
   streams are still finished before either error is reported.
7. **Depth fold.** Each bin's depth is computed in parallel; the worst bin is
   picked by the same loop, in bin order.
8. **Removal chances.** Each fragment's chance is computed in parallel, kept
   in fragment order.

Not split: tiling, the removal draw, phasing, the per-event loop, CRAM reads.

## Done means (locked)

1. **Byte-identical to master `cb756ec`** at `--threads` 1, 4 and 16:
   `R1.fq.gz` and `R2.fq.gz` (the gzip bytes and their content),
   `replaced_reads.txt`, and `truth.vcf` without its `##` lines, on:
   - the 4 clean and 5 origin sets of `fu-check.sh` (a 100 kb slice of the
     35x HG002 BAM, the whole 35x BAM, and the NA18488 background);
   - `del:` and `dup:chr20:14550000-17550000` on the 35x HG002 BAM, clean and
     origin;
   - one set from a CRAM of the slice, clean and origin.

   The log is identical too once timestamps are removed, apart from any line
   that names the thread count.
2. **Tests.** The suite passes. Each split piece has a test that compares its
   result at 1 thread with its result at 4 on input the split really cuts, and
   each such test is shown to fail under a deliberate mutation (chunks joined
   in the wrong order, or a record owned by two chunks).
3. **Clippy**: no new warning (12 now).
4. **Speed, reported as measured**: wall time and peak memory at `--threads`
   1, 4, 8 and 16 for the two 3 Mb events, clean and origin, and for
   `del:chr20:7119236-7120236` (1 kb), against master `cb756ec`.
   - **Kept only if** at `--threads 8` both 3 Mb events under origin are at
     least 25% faster than master, and no run at `--threads 1` is more than
     5% slower than master (or 0.2 s, whichever is larger).
   - PREDICTED (not run): at 8 threads origin about 30 s (from 68 s), clean
     about 26 s (from 38 s). Tiling (15-17 s) stays on one thread, so step 1
     cannot go below about 24 s.

If a gate fails, the branch is not offered for merge; the result says why.

## Result (2026-09-28): REFUTED by the `--threads 1` gate

Built as planned (`17fe11d`..`3db68cc`, binary md5 `d20c7c4b`; master
`cb756ec` binary `fd4ad3f5`). Gates 1-3 pass; gate 4 fails on its
`--threads 1` half, so the branch is not offered for merge.

1. **Byte-identical: PASS.** All 15 sets (the 11 small ones plus the four
   3 Mb runs) give the same signature as master at `--threads` 1, 4 and 16:
   exit code, gzip and content md5 of R1 and R2, `replaced_reads.txt`, the
   truth VCF body and the log without timestamps. The check can fail: the
   same `c1` and `o1` sets at `--seed 2` give different FASTQ md5s and logs,
   and the comparison reports them different.
2. **Tests: PASS.** 625 pass, 1 ignored (613 on master). Each split piece has a
   many-thread vs one-thread test, each shown red under a mutation.
3. **Clippy: PASS.** 12 bin / 14 test warnings, the same as master.
4. **Speed: FAIL.** One run each, one at a time, after a master warm-up.
   Master at `--threads 1` is the baseline (master's `src/` has no thread
   code, so its thread count does not change its own work).

   | Event | Model | master | t1 | t4 | t8 | t16 |
   |---|---|---|---|---|---|---|
   | del 1 kb (`chr20:7119236-7120236`) | clean | 0.4 s | 0.5 s | 0.5 s | 0.5 s | 0.5 s |
   | del 1 kb | origin | 0.6 s | 0.8 s | 0.8 s | 0.7 s | 0.7 s |
   | del 3 Mb (`chr20:14550000-17550000`) | clean | 18.2 s | 20.0 s (+9.9%) | 9.4 s | 7.8 s | 7.1 s |
   | del 3 Mb | origin | 49.0 s | 57.5 s (+17.3%) | 23.7 s | **19.0 s (-61.2%)** | 16.0 s |
   | dup 3 Mb | clean | 35.2 s | 37.4 s (+6.2%) | 25.5 s | 23.7 s | 23.0 s |
   | dup 3 Mb | origin | 64.1 s | 71.1 s (+10.9%) | 40.0 s | **36.7 s (-42.7%)** | 32.7 s |

   - At `--threads 8` both 3 Mb origin events are well over 25% faster:
     this half passes.
   - At `--threads 1` all four 3 Mb runs are more than 5% slower than
     master: this half fails. The 1 kb runs are within 0.2 s (timing
     resolution 0.1 s).
   - A second pair of runs repeated it (3 Mb del, `--threads 1`): clean
     19.2 s vs 21.8 s (+13%), origin 49.2 s vs 58.2 s (+18%). The log's
     1-second stamps put the loss in five or six steps of 1-3 s each, not in
     one place. The largest is origin's look-alike read (18 s to 21 s). The
     others are the census, the quality profile, the fragment grouping and
     the removal chances.
   - Peak memory rises at every count: 3 Mb runs +12% to +49%; the 1 kb
     origin run grows from 188 MB to 538 MB at 16 threads. Why was not
     measured (each look-alike batch opening its own reader is one guess).

   PREDICTED at t8 was origin about 30 s and clean about 26 s. Measured:
   origin 19.0 / 36.7 s, clean 7.8 / 23.7 s.

**Not done:** the README "Run time" section still describes master, as
this branch is not kept.
