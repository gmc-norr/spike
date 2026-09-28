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

## Result (2026-09-28): SUPPORTED in all six bins and directions

Run with `35e38c3`, spike master `fe15a46` (binary md5 `0302a168`),
16 threads, 12:13-12:39 (25 minutes); raw output in the session scratchpad,
`transplant/full/`. No pilot event was drawn again (0 of the 30 in the new
`events.tsv`).

**Events.** 100 per bin and direction; shared 100 / 100 / 78 (only 78 1-10 kb
shared deletions pass the exclusions). To keep 100 per transplant set, 93
forward and 99 reverse candidates were dropped because the recipient's own
BAM already showed more than 1 evidence read (N1).

**Refused by spike (RF8), so not compared:**

| Bin | Forward | Reverse |
|---|---|---|
| 50-299 | 1 | 4 |
| 300-999 | 6 | 6 |
| 1k-10k | 6 | 7 |

The verdicts below are about the deletions spike accepts.

**Verdicts** (`verdicts.tsv`): every bin and direction is **supported**; every
metric passes, and every broken control fails every metric it must fail.

**The pass rule, per bin** (`judge.tsv`; d = fake minus real, c = recipient
minus donor on shared events; pass needs median d inside c's 25th-75th and
d's 10th-90th width <= 1.5x c's):

| Set | Arm | Bin | Metric | n | median d | c 25th..75th | width d / c | pass |
|---|---|---|---|---|---|---|---|---|
| fwd | spike | 50-299 | J | 99 | 0.000 | -0.030..0.059 | 0.209 / 0.254 | yes |
| fwd | B1 | 50-299 | J | 99 | -0.031 | -0.030..0.059 | 0.316 / 0.254 | no |
| fwd | B2 | 50-299 | J | 98 | -0.083 | -0.030..0.059 | 0.451 / 0.254 | no |
| fwd | spike | 300-999 | J | 94 | 0.014 | -0.082..0.104 | 0.314 / 0.404 | yes |
| fwd | B1 | 300-999 | J | 94 | -0.187 | -0.082..0.104 | 0.448 / 0.404 | no |
| fwd | B2 | 300-999 | J | 96 | -0.331 | -0.082..0.104 | 0.585 / 0.404 | no |
| fwd | spike | 300-999 | E1 | 94 | 0.042 | -0.126..0.067 | 0.338 / 0.407 | yes |
| fwd | B1 | 300-999 | E1 | 94 | 0.252 | -0.126..0.067 | 0.348 / 0.407 | no |
| fwd | spike | 1k-10k | J | 94 | 0.085 | -0.073..0.158 | 0.415 / 0.415 | yes |
| fwd | B1 | 1k-10k | J | 94 | -0.194 | -0.073..0.158 | 0.416 / 0.415 | no |
| fwd | B2 | 1k-10k | J | 94 | -0.404 | -0.073..0.158 | 0.492 / 0.415 | no |
| fwd | spike | 1k-10k | E1 | 94 | 0.019 | -0.054..0.057 | 0.184 / 0.207 | yes |
| fwd | B1 | 1k-10k | E1 | 94 | 0.245 | -0.054..0.057 | 0.184 / 0.207 | no |
| rev | spike | 50-299 | J | 96 | 0.000 | -0.059..0.030 | 0.269 / 0.254 | yes |
| rev | spike | 300-999 | J | 94 | -0.014 | -0.104..0.082 | 0.416 / 0.404 | yes |
| rev | spike | 300-999 | E1 | 94 | 0.031 | -0.067..0.126 | 0.315 / 0.407 | yes |
| rev | spike | 1k-10k | J | 93 | -0.049 | -0.158..0.073 | 0.377 / 0.415 | yes |
| rev | spike | 1k-10k | E1 | 93 | 0.020 | -0.057..0.054 | 0.209 / 0.207 | yes |

(B2 is not asked about E1; for the record it fails it in both bins, median
d 0.249 and 0.070 against c's 75th percentile 0.067 and 0.057.)

**Median evidence per event** (J; E1 for >= 300 bp):

| Bin | fwd real | fwd spike | rev real | rev spike | shared HG001 / HG002 | sham |
|---|---|---|---|---|---|---|
| 50-299 J | 0.089 | 0.120 | 0.071 | 0.121 | 0.152 / 0.174 | 0.000 |
| 300-999 J | 0.408 | 0.419 | 0.412 | 0.420 | 0.459 / 0.466 | 0.000 |
| 1k-10k J | 0.486 | 0.563 | 0.446 | 0.443 | 0.513 / 0.562 | 0.000 |
| 300-999 E1 | 0.524 | 0.560 | 0.536 | 0.534 | 0.603 / 0.577 | 0.979 / 0.975 |
| 1k-10k E1 | 0.521 | 0.524 | 0.526 | 0.540 | 0.534 / 0.529 | 1.005 / 1.012 |

The sham (spike's replaced reads realigned unchanged) leaves no junction
evidence and no depth drop: realigning alone does not make or hide a
deletion at these events.

**Weak spot: 50-299 bp.** There the half-dose control B1 fails the rule by
0.001 (median d -0.031 against c's 25th percentile -0.030). The rule's median
is diluted by ties: J is 0 in both samples for 20 of 99 forward events (22
of 96 reverse; 17 of 100 in c), mostly where the measure sees no reads at
all. Seen, not judged: B1 moves J down for 60 events and up for 15 (mean
-0.091), while spike's fake moves it up for 42 and down for 37 (mean +0.001),
so the half dose is visible, but the locked median rule only just catches
it. The 50-299 verdict stands as the rule gives it; it is weaker evidence
than the two larger bins, where B1 and B2 miss c's middle half by 0.10 to
0.33.

**Not tested here:** read-level detail (where clips land, MAPQ, insert-size
spread), a "which is fake" test, other SV types, `--edit-model origin`,
SV callers.
