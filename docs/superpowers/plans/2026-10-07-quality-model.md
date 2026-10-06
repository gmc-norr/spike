# spike's reads get the sample's read-to-read quality and the errors that go with it

**Asked 2026-10-07.**
- The user chose to reproduce only the soft clips that come from reads read badly ("bad ends"), and to leave adapter, chimeric and foreign clips out.
- They then chose an fqzcomp-style quality model with an error table (option 1).
- The analysis behind this is `docs/analysis/soft-clips/` on branch `soft-clip-analysis`.

## Found

**spike's quality model, `QualityProfile` (`src/synth.rs`).**
- Each base's quality is drawn from P(Q | mate, cycle, base, previous quality in 4 bins), with fallbacks.
- With a memory of one base, every read comes out average.
- Errors are `p_err = 10^(-Q/10)`, whatever the read.

**Measured on the HG002 35x chr20 slice** (real reads and spike's reads in the 25 windows of the SNV run; `docs/presentation`, Figure 12):

| | real | spike, master |
|---|---|---|
| SD of each read's mean quality | 2.08 | 0.52 |
| reads all Q37 | 5.5% | 0.02% |
| reads with mean below Q33 | 8.4% | 0% |

**Why real reads are soft-clipped** (488,443 reads in the slice, 2.75% clipped):

| Cause | Of clipped reads |
|---|---|
| bad end | 44.7% |
| adapter | 25.7% |
| the sample's own variants | 19.7% |
| chimeric fragments | 4.7% |
| not in the reference | 3.6% |
| simple repeats | 1.6% |

**Errors depend on the read, not only on Q.** Aligned bases only, good reads against poor ones: Q25 0.07% against 0.95%.

**Offline tests (Python; `docs/analysis/soft-clips/qmodel/` and scratch `e10.py`, `e11.py`).**
- Reads were rebuilt from the reference at 5,190 held-out pairs' positions, then aligned with bwa-mem2.
- With other clips masked, the real bad-end clipped share is 1.15% (0.96-1.37).

| | bad-end clipped | read-mean SD | crashed | R1-R2 corr |
|---|---|---|---|---|
| real | 1.15% | 2.02 | 1.04% | 0.55 |
| fqzcomp-style context + error table, full clip rules | 1.18% (0.99-1.40) | 2.04 | 0.80% | 0.54 |
| the same with the simple clip rule below | 1.19% (1.00-1.42) | | | |
| copying a donor pair's strings and error positions | 1.49% (1.28-1.75) | 1.91 | 0.81% | 0.53 |

**Pool size, with backoff** (`e10.py`):

| donor pairs | read-mean SD | crashed |
|---|---|---|
| 1,000 | 1.76 | 0.41% |
| 5,000 | 1.98 | 0.92% |
| 200,000 | 2.02 | 0.93% |

**The bases** carry about 1% of what can be predicted (held-out bits per quality 0.383 → 0.379). Bigger contexts were worse.

## Gate A

1. **Principle.** spike's reads are read as badly, and in the same way, as the sample's own reads: how quality holds together along a read, how it differs from read to read and between mates, and the errors that come with it. All of it is learned from the sample.
2. **What would kill it.**
   - After alignment, the share of spike's reads with a bad-end clip differs from the sample's own reads beside them (K2).
   - Or spike's read-to-read quality spread, crash share or perfect-read share stays apart from the sample's (K1).
3. **Refuted before?**
   - Copying whole quality strings (`quality-tails`, `de806dd`) fixed crashed tails but not clips. Errors stayed at 10^(-Q/10).
   - That error part is what this adds.
   - The first copy-with-error-spots comparison was confounded by adapter and foreign tails (case file, 2026-10-06). The tests above mask them.
4. **Simplest thing.**
   - Copying a donor pair's strings and error positions is simpler, but it overshoots the clip share: 1.49% (1.28-1.75) against 1.15% (0.96-1.37).
   - The context model's extra parts are the read selector, the mate link and the error table. Each was measured above to be needed.
   - The base is left out: about 1%, and no measurement demands it.
5. **Inputs, and how each was checked.**
   - **Raw qualities.** Neither BAM has been recalibrated:
     - 35x HG002: @PG is bwa-mem2 and samtools only;
     - hospital BAM: bwa-mem2, MarkDuplicates and samtools.
   - **The donor pool's MAPQ ≥ 20 filter** keeps 95.1% of bad-end reads against 99.4% of unclipped ones. In the pool, bad ends are 1.18% of reads against 1.23%.
   - **Pool sizes.** A 10 kb DEL's pool on HG002 chr20 is 4,595 pairs (README).
   - **Quality alphabets.**
     - 35x HG002: 4 values (2, 11, 25, 37).
     - hospital BAM: 4 values (2, 9, 24, 40).
     - `HG002.GRCh38.chr20.bam` (novoalign): 31 values. So the history must bin a large alphabet.
   - **ReadPair keeps no alignment today.** Mismatches against the reference need each donor read's position, strand and CIGAR.

## Design (locked)

**Alignment kept at extraction.**
- `ReadPair` gains `align: Option<[MateAlignment; 2]>`, indexed R1 then R2. `MateAlignment` holds the start (0-based), whether the read is reverse, and the CIGAR ops.
- Both the BAM and the CRAM parsers fill it.
- spike's own synthetic pairs set `None`.

**Donor errors** (`learn_errors`, with the reference, when the profile is built). For each donor mate, in sequencing order, every base gets "error", "not error" or "not counted":
- **Aligned bases** (M, =, X): error if the base differs from the reference. `N` is not counted.
- **Inserted bases** are not counted.
- **Soft-clipped ends** are placed where they would have aligned.
  - An end of 1-4 bases, or one whose bases match the reference at ≥ 50% of placed positions, is counted base by base.
  - Any other clipped end (adapter, foreign, chimeric) is not counted.
- **The sample's own variants.** Reference positions where ≥ 5 pool reads have an aligned base and ≥ 10% of them differ from the reference are not counted.

**The quality model** (`QualityProfile`, replaces the Markov bins).
- **Alphabet:** the sorted distinct quality bytes of the pool.
- **History:**
  - With ≤ 4 values, the last 5 qualities.
  - Otherwise, the last quality exactly plus the two before it in `prev_q_bin`'s four bins.
  - Cycle 0 has an empty history.
- **Position:** `min(7, cycles_left >> pshift)`, with `pshift = max(0, round(log2(rl / 8)))`. For 151 cycles, pshift is 4.
- **Change flag:** whether the quality has changed bin (`prev_q_bin`) at least twice so far in the read.
- **Read selector:** 8 classes of a read's mean quality, cut at the pool's 1, 3, 8, 20, 40, 60 and 80% quantiles (both mates pooled).
  - Each synthetic pair draws (R1 class, R2 class) from the donor pairs' joint table, so mates share state.
- **Tables:** one per mate. Counts per context and quality, with backoff at < 20 observations:
  1. full context;
  2. (2-quality history, position, selector);
  3. (last quality, position, selector);
  4. (last quality, position);
  5. (position);
  6. all of the mate.
- **Draw:** a categorical draw over the alphabet from the first level with ≥ 20 observations.
- **`N`:** an `N` the read emits reports Q2 (L18), and Q2 enters the history.

**The error table.**
- Rate per (quality, selector, low run, distance from the 3′ end):
  - low run: consecutive bases below Q15 ending at this base, in bins 0, 1, 2-3, 4-7, 8+;
  - distance: cycles to the read's end, in bins 1-10, 11-30, 31-60, 61+.
- Learned from the counted donor bases.
- Backoff when a cell has < 200 counted bases: (quality, selector), then (quality), then `10^(-Q/10)`.
- A base is wrong when a uniform draw is below the rate. The indel/substitution split (`indel_error_rate`) is unchanged.

**Removed:**
- `PREV_Q_BINS` stays (the history uses it). `MIN_MARKOV_OBS`, `MIN_BASE_OBS`, `ProfileBins`, `PROFILE_CHUNK`, `from_read_pairs_in_chunks` and `sample_quality(_inner)` go.
- `MIN_DONOR_PAIRS` stays; its comment is rewritten.

**The log line** gives:
- the pool size and alphabet;
- the selector cuts;
- the share of donor bases counted for errors, and why the rest were not;
- each mate's crash share, real and as generated over 1,000 draws.

**Thin pool.** K4 sets `MIN_PROFILE_PAIRS`.

**Random stream.** It changes: spike's reads differ from master's at any seed. Nothing here expects byte identity.

**README.** The quality and error sections are rewritten. The deck in `docs/presentation` is left as it is (it describes 50ae8e1).

**Tests, written first and seen red against master's model:**
1. **Spread.** A pool where half the pairs are all-high and half carry long low runs gives generated reads whose read-mean SD is within 15% of the pool's.
2. **Mates.** A pool where R1 and R2 are poor together gives generated pairs whose mates' mean qualities correlate above 0.5.
3. **Errors by read.** A pool whose poor reads err at 30% at Q11 and good reads at 2% gives generated errors within 0.05 of each, by class.
4. **Learning errors.**
   - A clipped end that is reference sequence with 6 mismatches in 20 bases counts 6 errors.
   - A clipped adapter end counts none.
   - A position where every pool read differs counts none.
   - An inserted base counts none.
5. **Large alphabet.** A pool with 40 quality values generates only values from the pool.
6. **N.** An `N` reports Q2, and generation goes on.
7. **Reversed R1.** R1 keeps R1's model when R1 is the reverse mate. `test_reversed_r1_keeps_the_r1_quality_model`, adapted.
8. **Sequencing order.** The reverse mate's quality decays in sequencing order. L17's test, adapted.

Tests of the removed internals go: the chunk tests, `test_markov_*`, `test_quality_sampling_distribution`, `test_base_conditioned_vs_fallback` and `test_quality_chain_carries_the_q2_an_n_reported`.

**Mutation checks.** The unmutated suite is green first, then each mutant must turn a test red:
1. Each mate draws its own selector.
2. The selector is left out of the context.
3. Errors go back to `10^(-Q/10)`.
4. Every clipped end counts as errors.
5. The history keeps one quality only.

## Checks (locked)

**Shared definitions.** A read's clip is **bad-end** when it holds ≥ 5 bases matching the reference at ≥ 50% of placed positions, or 1-4 bases, and fewer than 3 other reads clip at the same boundary (±2 bp). A read is **crashed** when ≥ 10 of its last 20 qualities are below Q15.

**K1, quality strings.**
- **Data:** the 25-SNV run (`docs/presentation/scripts/runs.sh`, HG002 35x slice), new binary. Spike's reads against the real reads in the same 25 windows.
- **Pass:**
  - spike's read-mean SD is within 0.85-1.15 times the real one;
  - the perfect-read share (all at the top quality) and the crash share each have a two-proportion |z| < 3.
- **Control:** master's run must fail at least one of the three. Master's: SD 0.52 against 2.08, all-Q37 0.02% against 5.5%.
- **Reported:** R1-R2 correlation of mean quality, and the share with mean below Q33.

**K2, soft clips after alignment.**
- **Data:**
  - the 20 events of the read-level look, re-planted on the hospital BAM with the hospital's bwa-mem2 command (as for `quality-tails`);
  - reads at MAPQ ≥ 20 in the flanks.
- **Pass:** spike's reads against the sample's own reads beside them, on the share with a bad-end clip: two-proportion |z| < 3.
- **Control:** master's binary, the same way, must give |z| ≥ 3.
- **Reported:** the same on the SNV run's windows (35x), the any-clip share, and mismatches per 100 aligned bases.

**K3, mates** (reported). P(both crash) / (P(R1) · P(R2)) at the K2 events. The sample's own is 11.6.

**K4, pool size** (reported; sets the warning).
- An ignored test, `measure_quality_pool_size`, runs on the 35x BAM window chr20:38,402,500-38,602,500. That is 200 kb: N7's 30 kb window holds about 2,400 pairs, too few for a 5,000-pair pool.
- It draws donor pools of 500, 1,000, 2,000 and 5,000 pairs (20 repeats each) and generates the held-out half's strings.
- It reports K1's three metrics against the held-out reads.
- `MIN_PROFILE_PAIRS` becomes the smallest size whose medians pass K1's rule, as do all larger sizes.

**K5, base dependence** (reported). The share below Q15 by called base, spike against the sample's own, at the K2 events.

**K6, a large alphabet** (reported).
- One SNV and one 10 kb DEL from the K7 runs are planted on `HG002.GRCh38.chr20.bam` (31 quality values).
- The run must finish, and K1's metrics are reported.

**K7, no regressions.** The round trips of `docs/presentation`, new binary: DEL+INS (15/15 on master), small indels (12/12), hom DEL (10/10) and the 25 SNVs (71/78).
- **Pass:** every row that passes on master passes, except that at most 1 `allele_freq` row may flip, since the random stream differs.

**Tests and speed.**
- The whole suite passes, and the 5 mutants are caught.
- Speed: `dup:chr20:14550000-15550000` on the 35x BAM, `--seed 1 --threads 8`, 3 runs each. The new binary's median wall time is at most 1.25x master's.
