# Transplant test, round 2b: the full run, after the round 2 pilot stopped

Locked 2026-09-28, before the full run and before its code, with the user's
agreement ("1"). Branch `transplant`. Everything in
`2026-09-28-transplant-round2.md` holds except the three changes below:
samples, spiking, sham, groups, evidence (A, E, J), the pass rule and the
broken controls.

## Why a round 2b

The pilot (`f1847c0`) stopped because B1 passed two metrics:
- **INS20-49 A:** too few exact 20-49 bp insertions per event to see a half
  dose in 10 events. E saw it.
- **DUP50-299 J:** J reads the edges of one copy. Inside a longer tandem
  array the reads show the change elsewhere or nowhere, so J was about 0 in
  the real sample too. Where the copy was unique DNA, real J was 0.38-1.12 in
  all 7 such events.

## Changes

1. **Duplications in unique DNA only.** A 50-299 bp duplication is drawn only
   if its copy is not part of a longer tandem array: the copy, extended along
   the reference while the reference repeats it (`repeat_region`, as for
   indels), spans less than 1.10x the copy's length.
   - Measured before leaving out the pilot's: 61 forward, 65 reverse, 35
     shared.
   - All of them are drawn, since that is fewer than 100. The DUP verdicts
     rest on about half the events the other groups get, and on about a third
     of the shared controls.
2. **No pilot stop.** The full run is not gated on the pilot again. A metric
   whose must-fail control passes is inconclusive, as `verdict` already rules
   (round 1b's rule). INS20-49 A, or any other A or J that B1 cannot move,
   therefore makes its group inconclusive, not supported.
3. **Fresh events:** seed 2, and none of the pilot's 90 events (10 per set in
   SNV, INS20-49 and DUP50-299) can be drawn.

## Verdicts

Unchanged: per group and direction, supported / refuted / inconclusive, with
refusals (RF8), the sham and the broken controls reported beside them.

## Result (2026-09-28): 13 of 16 supported, 2 inconclusive, 1 refuted

Run with `6db499b`, spike master `fe15a46` (md5 `0302a168`), 16 threads,
14:25-15:09 (44 minutes); raw output in the session scratchpad,
`transplant2/full/`. No pilot event was drawn again (0 of 90).

**Events.** 100 per group and direction. DUP50-299 took all it had: 54
forward, 55 reverse and 33 shared, from pools of 61 / 65 / 35.

**Verdicts** (`verdicts.tsv`):

| Group | Forward | Reverse |
|---|---|---|
| SNV | supported | supported |
| DEL 1-4 | supported | supported |
| DEL 5-19 | supported | supported |
| DEL 20-49 | supported | supported |
| INS 1-4 | supported | supported |
| INS 5-19 | supported | supported |
| INS 20-49 | inconclusive | inconclusive |
| DUP 50-299 | supported | **refuted** |

**B1 fails every metric it must**, except in INS20-49. B1's median d:
- SNV A: -0.256;
- indel A: -0.21 to -0.26;
- indel E: -0.14 to -0.19;
- DUP J: -0.140.

B2 fails DUP J (median d -0.300).

**Refused by spike (RF8):**
- forward: 2 (both INS20-49);
- reverse: 4 (DEL5-19, DEL20-49, INS20-49 and DUP50-299, 1 each);
- B1: 2.

**Sham:** the median A or J of the replaced reads, realigned unchanged, is
0.000 in every group and direction.

**Median A (J for DUP): real / fake.**

| Group | Forward | Reverse |
|---|---|---|
| SNV | 0.500 / 0.473 | 0.468 / 0.482 |
| DEL 1-4 | 0.520 / 0.544 | 0.531 / 0.517 |
| DEL 5-19 | 0.500 / 0.500 | 0.529 / 0.467 |
| DEL 20-49 | 0.439 / 0.444 | 0.450 / 0.429 |
| INS 1-4 | 0.500 / 0.500 | 0.516 / 0.500 |
| INS 5-19 | 0.439 / 0.459 | 0.468 / 0.444 |
| INS 20-49 | 0.136 / 0.167 | 0.208 / 0.167 |
| DUP 50-299 | 0.300 / 0.377 | 0.384 / 0.348 |

**INS20-49 is inconclusive** because B1 passes both metrics there:
- On A, B1's median d is -0.023, inside c's -0.061..0.038. On E it is
  -0.105 against c's 25th percentile of -0.109.
- Seen, not judged: B1 still moves E down in 69 events and up in 26. So the
  half dose is visible, but real-vs-real scatter at this size is too wide
  for the median rule.
- spike's own d (E) is +0.042 forward and -0.010 reverse. The pilot's hint
  of RF15's direction (+0.165 reverse, on 10 events) is not seen at 99
  events. This round cannot confirm or clear RF15.

**DUP50-299 reverse is refuted on spread, not bias.**
- The median d is -0.005, inside c's -0.062..0.048.
- d's 10th-90th width is 0.482, against at most 1.5 x 0.233 = 0.350.
- Forward passes the same test narrowly (0.327).
- Seen, not judged: in both directions, 23 of 54 d values fall outside c's
  10th-90th range, where about 11 would if d spread like c. At the
  extremes, real J is 1.18 against fake 0.48 at chr7:134190878 (real 41
  pairs over a flank depth of 35), and fake exceeds real by 0.19-0.30 at six
  sites.
- So spike's 50-299 bp duplications match real ones on average, but single
  events stray from their real counterpart more than two real samples do.
- c rests on 33 shared events, so its width is itself uncertain. The plan
  named that weakness before the run.

**Not covered:** duplications of 300 bp and up; duplications inside longer
tandem arrays (340 of the 501 isolated real 50-299 bp ones); read-level detail;
`--edit-model origin`; callers.
