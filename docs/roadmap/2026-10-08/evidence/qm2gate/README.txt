Gate B checks for quality-model v2 (2026-10-07), before any plan:
regions.py  - quality stats vary by region: crashed 0.60-2.62% (35x), 0.28-1.13% (hospital) over 20 chromosome blocks
mix.py      - local class mix + pooled per-class rates: mean |crash error| 0.19 vs 0.40 points (35x), 0.12 vs 0.19 (hospital)
cells200k.py- hospital, ~207k pairs: crashed-tail error cells 9/16 with >= 200 bases (rest 123-191)
Rust temp test (reverted): model trained on 20 blocks' even-hash pairs, held-out odd-hash pairs (35x):
  current context: 120k pairs crash 0.94% vs 1.14% (z -6.54); 50k 0.99% (z -5.04); 5k 0.76% (z -13.42)
  + lows-in-last-16 (5 bins) at levels 0-1: 120k crash 1.04% vs 1.14% (z -3.17); SD ratio 1.010; perfect z 1.33

Bad-tail follow-up (2026-10-07; ../tails/ for the analysis). Throwaway Rust changes, applied in order by
hp_apply.py, hp_apply2.py, hp_apply3.py, hp_apply4.py (full diff: hp_all.patch; runs: hp_run*.sh, outputs
hp_run*.out). Same 20 blocks, even-hash train (120,138 pairs) / odd-hash held out (119,925), 35x BAM.
SPIKE_QX bits: 1 lows-in-16 in quality levels 0-1; 2 longest one-letter run so far (sticky; A/G/T 7-8,
9-11, 12+, C 5-6, 7+) in quality levels 0-3; 4 lows-in-16 replaces the error table's low run; 8 run bin in
the error table (+ a quality/lows/run backoff level); 16 pair classes drawn given each mate's run bin;
32 errors in the last 30 bases in the error table; 64 lows in quality levels 2-4 too.
Quality strings, crash share (real 1.14%) and by the read's longest-run bin (real 0.53/0.90/1.64/14.57%):
  QX 0  (v1)            0.96%  by bin 0.94/1.10/0.98/0.94    SD 0.997 perfect z 2.77
  QX 1                  1.06%  1.07/1.05/0.99/0.99
  QX 2  (run, no class) 0.75%  0.65/0.79/1.06/2.11   <- class is cut on the mean, which holds the crash
  QX 16 (class only)    0.97%  0.65/1.13/1.24/7.44           perfect z 3.60
  QX 18 (run + class)   0.96%  0.42/0.81/1.42/12.79
  QX 19 (+ lows)        0.98%  0.43/0.90/1.45/12.81          perfect z 0.63
  12 classes, QX 19     0.99%  0.42/0.88/1.45/13.41          (QX 0 with 12 classes: 0.99%, flat by bin)
  QX 83 / 82 (lows at all levels) 0.96% / 0.95%
  -> the run dependence is reproduced; the overall crash share stays ~13-15% short whatever was tried.
Errors on held-out bases given the real qualities, simulated along them (ratio to observed):
  crashed tails Q<15 / Q15-29 / Q30+; per-read tail-rate SD (observed 0.172)
  QX 0   0.80 / 0.58 / 0.44; 0.129
  QX 12  0.87 / 0.72 / 0.62; 0.123
  QX 36  0.89 / 0.78 / 0.70; 0.155      (QX 44 the same; 12 classes: 0.89 / 0.74 / 0.66; 0.155)
  not crashed, last 40: 1.02-1.08 throughout; all bases 1.00.

Clip test (user "1", 2026-10-07). RULE WRITTEN BEFORE RUNNING:
  Data: the held-out (odd-hash) pairs of the same 20 blocks, 35x BAM, with |TLEN| >= 151 in the input
  (no adapter). Three FASTQ sets, the same pairs: real (the reads as sequenced), v1 (QX 0) and new
  (QX 63 = bits 1,2,4,8,16,32), fake reads made from the reference at each mate's own 5' position and
  strand, substitution errors only. All three aligned with bwa-mem2 2.2.1 mem (defaults, -t 32).
  Scored the same way (clipscore.py): primary, MAPQ >= 20, proper pair, |TLEN| >= 151, 151 bp. A clip
  is bad-end when it is 1-4 bp or matches the reference at >= 50% of its placed bases, and no Q100
  variant lies within 10 bp of its boundary. Share of reads with a bad-end clip.
  PASS: new against real, two-proportion |z| < 3. v1 is the control and is only reported.
  Also reported: any clip, crashed, mismatches per 100 aligned bases, bad-end share by the template's
  longest-run bin.
(clip test: after hp_all.patch and hp_apply5.py, READ_CLASSES set back to 8 with the v1 quantiles)
RESULT (clip_run.out; full diff clip_all.patch): FAIL by the rule.
  set    reads   bad-end clip  95% CI       z vs real  any clip  crashed (z)    mm/100
  real   237549  1.40%         1.36-1.45%              2.13%     1.11%          0.394
  v1     237410  1.13%         1.09-1.18%   -8.31      1.22%     0.96% (-4.91)  0.328
  new    237310  1.26%         1.22-1.31%   -4.24      1.39%     0.98% (-4.29)  0.282
  bad-end share by the template's longest-run bin 0/1/2/3:
    real 0.99/1.27/2.35/8.76%   v1 1.12/1.10/1.17/1.29%   new 0.71/1.32/2.02/12.06%
Why (clipdiag.out, hpsource.out, after the verdict):
  - reads with an indel in the alignment, by bin: real 1.45/3.38/4.21/27.97%, new 0.01/0.03/0.02/0.20%:
    after long runs real reads slip (or carry the sample's own homopolymer variant); spike makes
    substitutions only, which the aligner clips instead.
  - run-free reads: crashed 0.45% vs 0.55%, mismatches 0.225 vs 0.346 per 100 (real includes variants).
  - the run memory was LEARNED FROM CALLED BASES: 22% of crashed reads on a run-free template show a run
    in their called bases (7.8% reach bin 3), against 0.5% of other reads. So training taught
    "poor read -> run": too few poor reads on run-free templates, too many after long runs.
    Next would learn the run bin from the reference under each donor read (as spike generates).

Clip test, rerun with two fixes (user "1", 2026-10-07). SAME data, scoring and rule as above (real.bam reused):
  new-A  (QX 63): the run bin (quality context, class draw, error table) learned from the reference under
         each donor read, from its 5' end in sequencing order, not from its called bases.
  new-AB (QX 191 = 63 + 128): new-A plus slips: a template run of one base (length >= 3) is read one base
         short or long at the rate donor reads show, by base and length bin (3-4, 5-6, 7-8, 9-11, 12-14,
         15+): an alignment indel inside the run's reference span counts; runs within 3 bp of a Q100
         variant (q100_blocks.tsv) or not fully aligned are left out. Delete/insert share learned too.
  PASS for each: against real, two-proportion |z| < 3 on the bad-end share.
RESULT (clip_run2.out): new-A 1.25% (z -4.43) FAIL; new-AB 1.31% (z -2.62) PASS by the rule.
  But learned slip rates are small (12-14: A 1.0%, T 1.6%; 15+: A 2.3%, T 1.0%); 50,780 runs masked as
  Q100 variants: the 28% indel share of real long-run reads is mostly the sample's own variants, not slips.
  new-A vs new-AB differ by ~1.8 SE, so before trusting the pass: seed check (seed 21 above; 22, 23 below).
SEED CHECK (clip_run3.out; diff clip_all2.patch): the seed-21 pass was a lucky draw.
  bad-end clip, real 1.40%:  new-A 1.25 / 1.22 / 1.28% (z -4.4 / -5.6 / -3.6)
                             new-AB 1.31 / 1.27 / 1.29% (z -2.6 / -3.9 / -3.3)  -> FAIL on 2 of 3 seeds.
  Slips add ~0.04 points on average (AB - A: +0.06, +0.05, +0.01).
  Left over, all seeds: run-free reads 0.78-0.83% vs 0.99% (crashed 0.51 vs 0.55%, mismatches 0.23 vs
  0.35/100 with the sample's own SNVs in real); reads with a 12+ run 10.4-10.9% vs 8.76% (too many).
  Progress: master ~0.2% (K2), v1 large-sample 1.13%, now 1.25-1.29%, real 1.40%.
