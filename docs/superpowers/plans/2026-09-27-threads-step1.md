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
