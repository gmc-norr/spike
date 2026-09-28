# Threads, step 1, second attempt: master's work at one thread

Locked 2026-09-28, before any code. Branch `threads`, on top of the first
attempt (`2026-09-27-threads-step1.md`, refuted by its `--threads 1` gate in
`3831455`). Chosen by the user: "Fix the 1-thread slowdown".

## Why the first attempt failed, measured

- At `--threads 1` every 3 Mb run was 6-17% slower than master and used
  15-34% more memory.
- Three rounds on `del:chr20:14550000-17550000`, clean, `--threads 1`, one run
  at a time (session scratchpad `idc/gb/res.txt`):

  | binary | mean | runs | peak memory |
  |---|---|---|---|
  | master `cb756ec` | 19.39 s | 19.82, 18.50, 19.84 | 1230 MB |
  | pool and R1/R2 only (`17fe11d`) | 19.91 s | 20.90, 19.08, 19.76 | 1230 MB |
  | first attempt (`3db68cc`) | 20.78 s | 20.40, 20.83, 21.10 | 1643 MB |

  Most of the loss, and all of the extra memory, comes after the pool commit:
  in the pieces. Master alone spans 18.50-19.84 s (7%), so one run each cannot
  judge a 5% rule.
- At one thread the pieces still do work master did not: extra readers, the
  records of one chunk copied into a merged result, a clone of the whole spot
  (`concat`), a sort where master used a map, a list of all pass-1 records
  before the maps are filled, and one copy of every quality bin (the pool's
  `reduce` joins the only chunk into an empty identity).

## The rule

At one thread, and for any region read in one chunk, each piece does
**master's work**: the same reader, the same passes, the same data
structures, and no chunking, merging or copying on top. Output stays
byte-identical, as before. Work at more than one thread is unchanged, except
where the same copy is removed there too (marked "all counts").

`one_thread()` is `rayon::current_num_threads() == 1`.

| # | Piece | One thread / one chunk does |
|---|---|---|
| 1 | R1/R2 write | master's single pass writing both encoders in turn; both still finished before either error is raised |
| 2 | Quality profile | one set of bins filled in one pass (no `par_chunks`, no join); all counts: `reduce_with` in place of `reduce` with an empty identity |
| 3 | Look-alike reads | master's loop: the one `Source` of `gather`, region by region, each result taken before the next is read |
| 4 | Origin spot | one chunk: `Source::scan`, master's one query; all counts: the chunks are moved into one list, not cloned (`concat`) |
| 5 | `fragments()` | master's `BTreeMap` grouping |
| 6 | SNP count and collect passes | all counts: the chunks are merged into the first chunk's map, so one chunk is the result as read |
| 7 | Extraction pass 1 | each chunk fills its own maps as master did; all counts: merged into the first chunk's maps in chunk order (a later record of a name replaces an earlier one, as in one pass) |
| 8 | Census | one chunk: master's `for_each_alignment` |
| 9 | Removal chances | sequential filter and map, as master |

Not changed: the depth fold (at one thread it makes master's calls, plus a list
of its bins and their depths; no reader, copy or merge), the pool itself, and
everything not split. rayon 1.12.0's `reduce` starts each fold from
`identity()` and calls `op(identity, item)` (`src/iter/reduce.rs` lines 45-49
and 92-95), which is where piece 2's copy comes from.

## Done means (locked)

The same four gates as the first attempt; only how speed is measured changes,
since one run each cannot tell 5% from noise.

1. **Byte-identical to master `cb756ec`** at `--threads` 1, 4 and 16 on the
   same 15 sets (`idc/run.sh ... full`, `idc/sig.sh`), logs included, against
   `idc/ref.sig`. The `--seed 2` control (`c1`, `o1`) must still compare as
   different.
2. **Tests.** The suite passes. Each one-thread or one-chunk path above (1-9)
   runs in a test against an independent reference (a one-thread pool or one
   chunk, beside the many-thread result or a value from the fixture), and
   each such test is shown red under a mutation of that path.
3. **Clippy**: no new warning (12 bin / 14 test).
4. **Speed.**
   - How: per event and model, a master warm-up, then three rounds; each
     round runs master at `--threads 1`, the branch at 1, and the branch at
     8, one at a time. Means of the three are compared. The branch at 4 and
     16 runs once, reported only. Peak memory is reported as the mean.
   - Events: `del:` and `dup:chr20:14550000-17550000` and
     `del:chr20:7119236-7120236`, clean and origin, on the 35x HG002 BAM,
     `--seed 1 --allow-resistant`.
   - **Kept only if** at `--threads 8` both 3 Mb events under origin are at
     least 25% faster than master, and no event and model at `--threads 1` is
     more than 5% slower than master (or 0.2 s, whichever is larger).
   - PREDICTED (not run): `--threads 1` within 3% of master in time and 5% in
     peak memory; `--threads 8` as fast as the first attempt or faster.

If a gate fails, the branch is not offered for merge; the result says why.
