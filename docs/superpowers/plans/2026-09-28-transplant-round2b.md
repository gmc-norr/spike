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
