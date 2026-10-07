# spike's reads get the sample's quality and errors, learned once from a sample of the input, with the DNA's own runs

**Asked 2026-10-07.**
- The first quality model (`2026-10-07-quality-model.md`, result `bccc43a`) was refuted. One event's donor pool, about 2,800 pairs, was too thin for crashed reads and their tail errors.
- The user chose option 1 four times:
  1. learn once from a big sample of the input BAM, and speed it up;
  2. first look at what happens in the bad tails;
  3. test the parts found there cheaply;
  4. build v2 with the parts that helped, keeping the locked checks and reporting the big held-out clip test beside them.

## Found

All of this was measured before this plan, on the HG002 35x BAM unless said otherwise. Scripts and outputs are in the session scratchpad:
- `tails/`: the bad-tail analysis;
- `qm2gate/`: the throwaway model tests, with `README.txt` as the log, and the diffs `hp_all.patch`, `clip_all.patch` and `clip_all2.patch`.

**Quality differs by region** (`qm2gate/regions.py`, `mix.py`, `mix2.py`). Twenty 100 kb blocks, chr1-20 at 30% of each length:
- The crashed-read share runs 0.60-2.62% at 35x and 0.28-1.13% on the hospital BAM.
- Mean error of each block's crash share, in points:

| predicted from | 35x | hospital |
|---|---|---|
| the genome-wide share | 0.40 | 0.22 |
| the block's own read-class mix | 0.19 | 0.13 |
| the block's template one-letter runs | 0.29 | 0.13 |
| both (class × run) | 0.18 | 0.09 |

**What bad tails are** (`tails/`; the chr20 slice, sequencing orientation).
- **Quality values.** The values are 2, 11, 25 and 37, and Q2 only ever marks an `N`. Low quality (Q11) is mixed with Q25 and Q37 right to the read's end.
  - The run of low bases ending at the last base has a median of 1.
  - No crashed read stays low once low.
  - So v1's error-table bin "low qualities in a row" was almost always 0 or 1.
- **Runs of one letter in the template set off sticky, one-strand crashes.** After the run, low quality and errors rise for the rest of the read. Sixty cycles later it is still 7x for runs of 12+.

  | After a run of | Low quality | Errors |
  |---|---|---|
  | 9-11 (A/T) | 2.1-2.6x | 2.7-3.5x |
  | 12+ (A/T) | 6.9-9.4x | 14-22x |
  | C, 7-8 | 5.9x | 15x |

  - Poly-G hardly matters.
  - Reads with a run of 12+ crash 15.0% of the time, against 0.46% with no run longer than 4. About 38% of crashes come from runs of 7 or more.
  - At shared clip spots, same-strand reads go from 1.1% to 8.2% errors after the spot; the other strand stays at 0.8% and 0.8%.
- **Tails differ by read beyond their qualities.**
  - In crashed tails, Q11 bases err at 12% in mild tails and 51% in bad ones; Q37 bases err at 0.5% and 7.3%.
  - The per-read rate has SD 0.160, against 0.108 if every read erred alike.
  - Not the cause: tails are not shifted by one base (0.9% of clips), and errors barely clump within a tail (0.435 against 0.408).

**The parts, tested in a throwaway model** (`qm2gate/`).
- **Setup.** The 20 blocks, split by name hash: 120,138 pairs to train and 119,925 held out.
- **Quality strings.** Crashed share overall, and by the read's longest-run bin (none / 7-8 / 9-11 / 12+). Real: 1.14% overall; 0.53 / 0.90 / 1.64 / 14.57% by bin.

  | model | crashed | by run bin |
  |---|---|---|
  | v1 context, 120k pairs | 0.96% | 0.94 / 1.10 / 0.98 / 0.94% |
  | + run bin in the context | 0.75% | 0.65 / 0.79 / 1.06 / 2.11% |
  | + classes drawn given the run | 0.96% | 0.42 / 0.81 / 1.42 / 12.79% |
  | + lows among the last 16 | 0.98% | 0.43 / 0.90 / 1.45 / 12.81% |

  - The class is cut on the read's mean quality, which already holds the crash a run set off. Without the class draw given the run, the run bin makes things worse.
  - Three more tries did not move the overall share: 12 finer classes (0.99%), lows at every fallback level (0.96%), and v1 with 12 classes (0.99%, flat by bin).
- **Errors simulated along the real held-out qualities.** Crashed tails at Q<15 / Q15-29 / Q30+, as a share of the real rates (0.339 / 0.113 / 0.032), and the per-read spread (real SD 0.172):

  | model | crashed tails | per-read SD |
  |---|---|---|
  | v1 table (low run, class, end) | 0.80 / 0.58 / 0.44 | 0.129 |
  | lows among the last 16 instead of the low run, and the run bin | 0.87 / 0.72 / 0.62 | 0.123 |
  | + errors among the last 30 bases | 0.89 / 0.78 / 0.70 | 0.155 |

  Reads that did not crash, last 40 cycles: 1.02-1.08 throughout.
- **The run bin must be learned from the reference under each read, not its called bases.** Crashed real reads on a run-free template show a run in their called bases 22% of the time, and 7.8% reach the 12+ bin. Other reads do so 0.5% of the time.
- **Clip test** (rule written before the run, `qm2gate/README.txt`).
  - The held-out pairs with TLEN ≥ 151, made again from the reference at each mate's own 5′ position and aligned with bwa-mem2 2.2.1.
  - A clip is bad-end when it is 1-4 bp or ≥ 50% reference, with no Q100 variant within 10 bp.
  - Real reads: 1.40%. By run bin: 0.99 / 1.27 / 2.35 / 8.76%.

  | model | bad-end clips, seeds 21 / 22 / 23 | by run bin (seed 21) |
  |---|---|---|
  | v1 context | 1.13% (z −8.3) | 1.12 / 1.10 / 1.17 / 1.29% |
  | run, class, lows, error history; run from called bases | 1.26% (z −4.2) | 0.71 / 1.32 / 2.02 / 12.06% |
  | the same, run from the reference (A) | 1.25 / 1.22 / 1.28% | 0.80 / 0.86 / 2.16 / 10.47% |
  | A + slips | 1.31 / 1.27 / 1.29% (z −2.6 / −3.9 / −3.3) | 0.83 / 1.02 / 2.41 / 10.52% |

  - **Slips are small.** Learned with Q100 variants masked: runs of 12-14 slip 1.0% (A) and 1.6% (T). The 28% indel share of real reads with a 12+ run is mostly the sample's own variants. Slips add about 0.04 points.
  - **What is left:** run-free reads clip 0.8% against 0.99%, and reads with a 12+ run clip 10.5% against 8.8%.
  - For scale: master made about 0.2% bad-end clips at the hospital events (`bccc43a`, K2).
- **Startup cost** (8 threads, 240,063 pairs, v1 code):
  - reading 20 blocks takes 0.6 s in parallel (3.7 s one by one);
  - learning the qualities takes 1.4 s;
  - learning the errors takes about 19 s, on one thread with hash maps.

## Gate A

1. **Principle.** spike's reads are read as badly, and in the same way, as the sample's own reads:
   - how quality holds together along a read, differs from read to read and between mates, and goes wrong after the DNA's own runs of one letter;
   - and the errors that come with it.

   All of it is learned from the sample.
2. **What would kill it.**
   - After alignment, the share of spike's reads with a bad-end clip differs from the sample's own reads beside them (K2).
   - Or spike's read-to-read spread, crash share or perfect-read share stays apart from the sample's (K1).
3. **Refuted before.**
   - v1 learned from each event's pool, which was too thin (`bccc43a`). This learns once, from a sample about 30-50x larger (about 120,000 pairs against pools of 2,278-3,778).
   - Run memory learned from called bases was refuted in the clip test; it is learned from the reference here.
   - Slips (about 0.04 points), 12 classes, and lows at every level were measured to add nothing worth their cost, and are left out.
4. **Simplest thing.** Each part has the measurement above behind it:

   | part | measurement |
   |---|---|
   | the big sample | v1's refutation |
   | run bin and the class drawn given it | crash by run bin; clips by run bin |
   | lows among the last 16 | crash 0.96 → 1.06% alone; error ratios |
   | error history | per-read spread 0.123 → 0.155 |
   | the event's class mix | region error 0.29 → 0.18 points with the run mix |
   | the error rewrite | 19 s |

5. **Inputs, and how each was checked.**
   - Raw qualities: checked for v1, from the @PG lines of both BAMs.
   - spike already needs an indexed BAM or CRAM and the reference FASTA.
   - **That the index's leaf bins (BAI) or slices (CRAI) mark where reads are** is not checked yet. Test T8 checks it on a slice-shaped BAM, and K4 reports the blocks found on the real inputs.
   - Quality alphabets: 4 values on both NovaSeq BAMs, 31 on `HG002.GRCh38.chr20.bam` (K6).

## Design (locked)

**The startup sample** (`quality::sample_input`, once per run, after `bam_stats`).
- **Where.** The input's windows with reads, from its index:
  - BAI: leaf bins with chunks (16 kb windows);
  - CRAI: each slice's reference span.

  These windows are listed in contig order.
- **Blocks.** 20 blocks of 50 kb.
  - Block `k` starts at window `floor((k + 0.5) · W / 20)`.
  - A block that overlaps the one before it is moved to start where that one ends.
  - On a slice (one 1.7 Mb stretch of chr20) this still gives 20 blocks inside the reads.
- **Pairs.** From each block with `extract_read_pairs` (`--min-mapq`, same filters as donor pools), on the thread pool, then deduplicated by name.
  - 50 kb, not 100 kb, keeps the sample near the 120,000 pairs the clip test learned from at 35x: 100 kb blocks gave 240,063 pairs.
- **Reference.** Each block's bases are read with `reference::fetch_window`, padded by 1 kb.
- **Thin sample.** Below 50,000 pairs, the run goes on with a warning that gives the count; 50k reached 0.99% crashed in the gate, against 0.94% at 120k. Below `MIN_DONOR_PAIRS` (30) pairs, the run stops, as for a pool.

**The model** (`QualityProfile`, learned once from the sample).
- **Template.** A donor mate's template is the reference from its 5′ end (alignment start less a leading clip; or alignment end plus a trailing clip, reverse-complemented), as many bases as the read. Runs are read from it, never from the called bases.
- **Run bin.** The longest run of one base read so far, sticky:
  - A/G/T: 7-8 → 1, 9-11 → 2, 12+ → 3;
  - C: 5-6 → 2, 7+ → 3;
  - otherwise 0.
- **Lows bin.** Qualities below Q15 among the last 16: 0, 1-2, 3-5, 6-9, 10+.
- **Quality context.** v1's 6 backoff levels, with the lows bin added to levels 1-2 and the run bin to levels 1-4. Everything else is as in v1: history, position, change flag, class, `MIN_CONTEXT_OBS` 20, `N` → Q2.
- **Classes.** v1's 8 classes and quantiles. A pair's classes are drawn with weights:

  `joint(c1, c2) · Π_mates [P_m(c | h_m) / P_m(c)] · r_m(c)`

  - `h_m` is the run bin over the mate's template.
  - `P_m(c | h)` and `P_m(c)` come from the sample (with +0.5 and +4 pseudo-counts).
  - `r_m(c)` is the event's class mix (below).
- **The event's class mix.** From the event's donor pool, classes cut at the sample's cuts:

  `r_m(c) = (n_m(c) + 0.5) / (Σ_h n_m(h) · P_m(c | h) + 0.5)`

  - `n_m(c)` counts the pool's mates of class `c`; `n_m(h)` counts them by run bin, from their templates.
  - So a region's poor reads count once, and the part its runs already explain is not counted twice.
  - Without a pool (unit tests), `r = 1`.
- **Error table.**
  - Rate per (quality, class, lows bin with the base itself, end bin, run bin, errors among the last 30 bases: 0, 1, 2-3, 4+).
  - Backoff at < 200 counted bases: (quality, lows, run, errors) → (quality, class) → (quality) → 10^(-Q/10).
  - The history counts errors the read made: the donor's own when learning, the drawn ones when generating.
  - Learning counts bases as in v1: clip rule, variant mask, inserts left out.
- **Slips:** none. Errors are substitutions, with `indel_error_rate` unchanged.

**Speed.**
- **Context tables.** A fast multiply-shift hasher for the context hash maps, written in the module (no new crate).
- **Error learning.** Per block, on the thread pool. Each block's pileup and reference are flat arrays over its window.
- **No per-event learning.** Each event only counts its pool's classes and run bins.

**Wiring.**
- `main.rs` builds the profile once and passes it by reference to every event's generator. `synth_generator` adds the event's class mix.
- `generate_from_template` passes the template base to the state and records each drawn error.
- The pair generators compute both templates' run bins before drawing the classes.

**Log lines.**
- At startup: the blocks (count, span, pairs), the alphabet, the class cuts, the error census, and the time taken.
- Per event, at debug level: the class mix `r`.

**Random stream.** It changes, as in v1. Nothing here expects byte identity.

**README.** The quality and error sections are rewritten: what is learned, from where, and the startup sample. The deck in `docs/presentation` is left alone.

**Tests, written first and seen red against v1's model where the behaviour is new:**
1. **Runs.** A pool whose reads crash after a 12-base run of T and not otherwise. Generated reads crash ≥ 5x more often after such a run than on run-free templates.
2. **Run from the template.** A pool whose crashed reads, on run-free templates, call `AAAAAAAAAA` in their tails. Run-free templates get the pool's crash share within 30% (relative), and templates with a 12+ run get no more than 1.3x it.
3. **Interleaved lows.** A pool whose poor reads alternate Q11 and Q37 in their tails, and err at 40% at Q11 there, against 5% at Q11 elsewhere. Generated Q11 tail bases of poor reads err within 0.05 of 0.40.
4. **Error history.** A pool where half the crashed reads err at 60% at Q11 and half at 10%. The per-read rate's SD over generated crashed reads is at least 0.8x the pool's.
5. **The event's class mix.** A sample of mixed reads, and an event pool of only poor reads. Reads made for that event have a mean class within 1 of the pool's.
6. **v1's tests still pass:** spread, mates, errors by read kind, two-base pattern, large alphabet, N, learning errors, reversed R1, sequencing order, one vs many threads.
7. **Sampler on a slice.** A BAM with reads only in 1.2 Mb of a 60 Mb contig. `sample_input` returns 20 non-overlapping blocks, all inside that stretch.
8. **Sampler on a CRAM**, the same.

**Mutation checks.** The unmutated suite is green first, then each mutant must turn a test red:
1. The run bin left out of the context.
2. Classes drawn ignoring the run.
3. The run bin read from the called bases.
4. The lows bin replaced by v1's low run.
5. The error history left out.
6. The event's class mix ignored (`r = 1`).
7. The sampler spreading blocks over contig lengths instead of the index's windows.

## Checks (locked)

**Shared definitions.** As v1:
- A read is **crashed** when ≥ 10 of its last 20 qualities are below Q15.
- A clip is **bad-end** when it holds ≥ 5 bases matching the reference at ≥ 50% of placed positions, or is 1-4 bases, and fewer than 3 other reads clip at the same boundary (±2 bp).

**K1, quality strings.** As v1:
- **Data:** the 25-SNV run on the chr20 slice (`docs/presentation/scripts/runs.sh`). spike's reads against the real reads in the same 25 windows.
- **Pass:** spike's read-mean SD is 0.85-1.15x the real one, and the perfect-read and crash shares each have |z| < 3.
- **Control:** master's run fails at least one of the three. On `bccc43a` it failed all three.
- **Reported:** R1-R2 correlation, the share with mean below Q33, and the crash share by run bin.

**K2, soft clips after alignment.** As v1:
- **Data:** the 22 events of the read-level look, re-planted on the hospital BAM with its bwa-mem2 command (scratch `inspect/`); reads at MAPQ ≥ 20 in the flanks.
- **Pass:** spike's reads against the sample's own reads beside them, on the bad-end share: |z| < 3.
- **Control:** master's binary gives |z| ≥ 3. On `bccc43a` it gave 0.17% against 1.38%, z −8.60.
- **Reported:** the any-clip share, mismatches per 100 aligned bases, and the bad-end share by run bin.
- **Power.** With about 8,000 reads a side, a shortfall like the clip test's (8%) would not show here: PREDICTED (not run) z ≈ −0.6. The next check is reported for that reason.

**K2b, the big held-out clip test** (reported, no pass rule).
- The gate's clip test, run on this branch's own code. An ignored test, `measure_clip_share`:
  - samples the 35x BAM with `sample_input`;
  - learns from the even-hash pairs;
  - writes the odd-hash pairs (TLEN ≥ 151) as sequenced, and as spike makes them from the reference at each mate's 5′ position with `generate_from_template`.
- Aligned with bwa-mem2 2.2.1 and scored with `qm2gate/clipscore.py`, copied into `docs/analysis/quality-model-v2/`.
- Reports the bad-end share, z, and the share by run bin for seeds 21, 22 and 23.
- The gate's prototype gave 1.25-1.28% against 1.40%.

**K3, mates** (reported). P(both crash) / (P(R1) · P(R2)) at the K2 events. The sample's own value is 11.6.

**K4, the startup sample** (reported). On the slice, the 35x BAM and the hospital BAM:
- the blocks found (first and last), the pairs, and the time the sample and the learning take.

**K5, base dependence** (reported). The share below Q15 by called base at the K2 events, as v1.

**K6, a large alphabet** (reported). As v1, on `HG002.GRCh38.chr20.bam` (31 values). The run must finish.

**K7, no regressions.** The round trips of `docs/presentation` (DEL+INS, small indels, hom DEL, 25 SNVs), run as v1 did in `qtq/runs_new.sh`.
- **Allowance.** For each round trip, master is run with `--seed` 42 (the default), 2, 3, 4 and 5. Its allowance is the most `allele_freq` rows whose pass/fail differs between any two of those five runs.
- **Pass:** at `--seed` 42, every row that passes on master's seed-42 run passes on the new binary, except that `allele_freq` rows may flip up to the allowance.

**Tests and speed.**
- The whole suite passes, and the 7 mutants are caught. clippy shows no new warnings.
- **Speed:** `dup:chr20:14550000-15550000` on the 35x BAM, `--seed 1 --threads 8`, 3 runs each. The new binary's median wall time is at most 1.25x master's (8.5-9.5 s on master for v1).

## Result (2026-10-07): K1, K2, K7, tests and mutants pass; speed fails

- **Code:** `e9a6f0c`.
- **Binaries (md5):** new `961e2a3c`, master `e21d282b` (50ae8e1, as for v1).
- **Scratch** (session scratchpad):
  - `qtq2/`: `speed.sh`, `runs.sh`, `all.sh`, `k7.py`, `k1runs.py`, `k2runs.py`, `k6score.py`, and the runs;
  - `qmp2/mutate.py`.

The locked checks K1, K2 and K7 pass, as do the tests and the mutants. The speed check fails: 1.27x against a limit of 1.25x. The big held-out clip test matches the sample.

**Deviations from the design**, all in `e9a6f0c`'s message:
- **Class smoothing.** P(class | run bin) is shrunk toward P(class) by 8 reads, not +0.5/+4. With the plan's pseudo-counts, a run bin the sample barely holds drew every class alike.
- **T3 is redesigned.** Every read has one class and a dense low stretch mid-read. In the plan's version the class alone told the reads apart.
- **T3b is added**, because mutant 4 survived the planned tests.

**Tests: 747 pass** (3 ignored), on stable and on Rust 1.82. clippy's warnings are master's.
- **Mutants: 7 of 7 caught**, each by its own test:

  | Mutant | Caught by |
  |---|---|
  | 1 | T1 |
  | 2 | T1 |
  | 3 | T2 |
  | 4 | T3b; it survived T1-T5 |
  | 5 | T4 |
  | 6 | T5 |
  | 7 | T7 |

- With all seven parts back to v1's behaviour at once, 6 tests go red.

**K1: PASS.** The 25-SNV run on the chr20 slice, seed 42.

| | spike, new | spike, master | real |
|---|---|---|---|
| read-mean SD | 1.973 (ratio 0.931, pass) | 0.516 (0.248) | 2.118 |
| perfect reads | 5.52% (z +0.19, pass) | 0.02% (z −27.17) | 5.47% |
| crashed reads | 0.88% (z −2.49, pass) | 0.01% (z −11.88) | 1.16% |
| R1-R2 corr (reported) | 0.514 | 0.034 | 0.543 |
| mean below Q33 (reported) | 7.85% | 0% | 8.63% |

- The control fires: master's run fails all three.
- **Crashed share by the template's longest-run bin** (none / 7-8 / 9-11 / 12+), reported: spike 0.36 / 0.56 / 1.22 / 12.69%, real 0.52 / 0.59 / 1.33 / 15.23%.

**K2: PASS.** The 22 events, re-planted on the hospital BAM; flank reads at MAPQ ≥ 20.

| | spike, new | spike, master | sample's own |
|---|---|---|---|
| reads | 7,579 | 7,636 | 8,899 |
| bad-end clip | 0.94% (z −2.64, pass) | 0.17% (z −8.60) | 1.38% |
| any soft clip (reported) | 0.94% | 0.17% | 2.42% |
| crashed (reported) | 1.49% (z +0.04) | 0.04% | 1.48% |
| mismatches per 100 aligned bases (reported) | 0.270 | 0.243 | 0.308 |

- **Bad-end share by run bin** (reported): spike 0.26 / 0.13 / 0.50 / 8.62%, own 0.73 / 0.98 / 0.95 / 8.16%.
- On this BAM the shortfall is in reads without a long run. The reads after runs of 12+ match.
- The pass is narrow: the plan predicted z ≈ −0.6 for an 8% gap, and this gap is 32%.

**K2b, the big held-out clip test** (reported). The 35x BAM sampled by `sample_input`. The model learned from 63,040 even-hash pairs, and 62,682 odd-hash pairs were written (`docs/analysis/quality-model-v2/k2b.sh`).

| set | bad-end clip | z | any clip | crashed (z) | by run bin |
|---|---|---|---|---|---|
| real | 1.26% | | 1.95% | 1.06% | 0.86 / 1.00 / 1.98 / 9.91% |
| spike, seed 21 | 1.23% | −0.86 | 1.33% | 0.90% (−3.89) | 0.77 / 0.73 / 1.99 / 11.67% |
| spike, seed 22 | 1.25% | −0.41 | 1.34% | 0.91% (−3.79) | 0.79 / 0.75 / 2.04 / 11.84% |
| spike, seed 23 | 1.27% | +0.14 | 1.37% | 0.91% (−3.67) | 0.83 / 0.67 / 1.97 / 11.92% |

- The real share is 1.26% here against 1.40% in the gate, because 50 kb blocks are other reads than the gate's 100 kb ones.
- The crashed share stays about 15% short, as in the gate. Reads after runs of 12+ clip about 19% too often.

**K3, mates (reported).** P(both crash) / (P(R1) · P(R2)) among spike's pairs at the K2 events: 4.74 (4,780 pairs), against the sample's 11.6 and v1's 10.16. Drawing each mate's class given its own run weakens the link.

**K4, the startup sample (reported).**

| input | blocks | pairs | read | learned |
|---|---|---|---|---|
| chr20 slice | 20, chr20:38,533,461 to 40,155,624 (all inside the slice) | 142,458 | 0.3 s | 1.4 s |
| 35x HG002 | 20, chr1 to chrX | 126,228 | 0.6 s | 1.2 s |
| hospital BAM | 20, chr1 to chrX | 113,938 | 2.8 s | 1.6 s |
| `HG002.GRCh38.chr20.bam` | 20, chr20 | 143,080 | 2.4 s | 3.3 s |

**K5, base dependence (reported).** The share below Q15 by called base at the K2 events:
- spike A 0.0161, C 0.0171, G 0.0168, T 0.0167;
- own A 0.0189, C 0.0185, G 0.0156, T 0.0124.

Still flat: the base itself is not in the context, only its runs.

**K6, a large alphabet (reported).** `HG002.GRCh38.chr20.bam` (novoalign), one SNV and the 10 kb DEL.
- The run finishes, and the BAM's 31 values come out, and only those.
- spike's 650 reads against 5,127 real reads in the two windows:
  - read-mean SD 4.51 against 3.64 (ratio 1.24);
  - crashed 24.3% against 14.4% (z +6.61).
- This alphabet now over-crashes. v1 had 13.4%.

**K7: PASS.** Master was run with `--seed` 42, 2, 3, 4 and 5.
- The observed values differ between seeds, but no row's pass/fail does, so each round trip's `allele_freq` allowance is 0.
- The new binary at seed 42: DEL+INS 15/15, small indels 12/12, hom DEL 10/10 and 25 SNVs 71/78, row for row as master. No row is lost.

**Speed: FAIL.** `dup:chr20:14550000-15550000`, 35x BAM, `--seed 1 --threads 8`, alternating, 3 runs each.
- Median: new 11.86 s against master 9.37 s, 1.27x. The limit is 1.25x.
- The startup sample takes 1.8 s of that: 0.6 s to read, 1.2 s to learn.

**What this means.**
- spike's reads now carry the sample's read-to-read spread, perfect and crashed reads, and crashes after the DNA's own runs (K1).
- They clip at a bad end as often as the sample's own reads on 125,000 held-out 35x reads (K2b). On the hospital events they clip 0.94% against 1.38% (K2, inside the limit).
- **Still off:**
  - crashed reads about 15% short on 35x (K1, K2b);
  - run-free reads clip too little on the hospital BAM (K2);
  - mates crash together less often (K3);
  - a 31-value alphabet crashes too often (K6);
  - the run is 1.27x slower on a single event.
