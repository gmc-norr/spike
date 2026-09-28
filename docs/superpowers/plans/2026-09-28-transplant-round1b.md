# Transplant test, round 1b: the full run, after the pilot stopped round 1

Locked 2026-09-28, before the full run, with the user's agreement ("1").
Branch `transplant`. Everything in `2026-09-28-transplant-round1.md` holds
-- samples, region, event sets, spiking, sham, evidence J and E1, the pass
rule -- except the three changes below.

## Why a round 1b

The pilot (`6004fcf`) stopped on round 1's second pilot check: the
real-vs-real spread of J for one event (10th-90th width 0.509) was wider
than the median real J (0.409). That check compared the scatter of single
events with the signal; the pass rule judges the median over up to 100
events, and the half-dose control B1 failed clearly (J 0.25 against 0.49),
so the test does see a wrong dose. The job the check was meant to do --
catch a test too blunt to fail -- is the broken controls' job, and they stay.

## Changes

1. **No per-event spread check.** The pilot's second check is dropped.
2. **Each broken control is judged on what it breaks**, per size bin, on the
   forward events:
   - B1 (VAF 0.25) must fail the rule on J and on E1 (E1 where measured);
   - B2 (planted 200 bp to the right) must fail the rule on J. It is not
     asked to fail E1: moving a 1-10 kb deletion 200 bp barely changes the
     depth over it, which is not what B2 tests.

   A bin's verdict on a metric is "inconclusive" if a control that must fail
   that metric passes it. The same bin's reverse verdict uses the forward
   controls (none are run on reverse).
3. **Fresh events:** seed 2, and none of the pilot's 30 events (10 forward,
   10 reverse, 10 shared) can be drawn.

## Verdicts

Per size bin (50-299, 300-999, 1k-10k), per direction, per metric: pass,
fail or inconclusive, as above. A bin and direction is **supported** when
every metric measured there passes, **refuted** when any fails (not
inconclusive), and **inconclusive** otherwise. Refusals (RF8) and the sham
are reported per bin beside the verdicts.
