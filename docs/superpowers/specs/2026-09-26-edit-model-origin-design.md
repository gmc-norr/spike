# `--edit-model origin`: edit reads by where they came from

Date: 2026-09-26. Status: design agreed with the user in four parts; revised 2026-09-27 for
five gaps a review found (R1-R5 below); not built.

## Goal

spike should change the reads the way a real variant changes them, including where the aligner
cannot tell a spot from a look-alike elsewhere. The user wants to use spike to ask whether an
event can be detected where no real example exists. If spike plants hard events as easier than
real ones, or plants nothing, that question cannot be answered.

**Principle.** A spiked event removes the reads that came from the edited copy, wherever the
aligner put them, and replaces them with reads made from the edited copy, which the aligner
then places. Every other read stays as it is.

## What happens today (measured)

- The donor pool takes a read pair only if both mates pass `--min-mapq` (default 20), form a
  proper pair and are not duplicate, secondary, supplementary or QC-fail
  (`src/extract.rs:104-109`). Only pool pairs can be removed. Everything else stays.
- **At chr20:7119236** (RF14's first kill test, 1 kb deletion, 35x HG002), 382 primary,
  non-duplicate reads lie over the event:
  - 359 are below MAPQ 20, and 242 of those are at MAPQ 0;
  - 6 are not in a proper pair;
  - 23 were editable.

  spike tiled new reads at 8.5x fragment depth, against 40.4x for the donor nearby
  (`SIM_DEPTH_FOLD` 4.37). The deletion was barely planted.
- By default spike refuses such an event (RF8, `SIM_RESIST` above 0.5).
- `--min-mapq 0` makes the MAPQ 0 reads editable. At RF8's locus it cut the uneditable reads
  from 5731 to 58 of 5749 (README). It removes them as if they certainly came from here,
  though.
- **A pool can be empty over the whole footprint and still pass `MIN_DONOR_PAIRS`.** That
  guard (30 pairs, `src/main.rs:1383`) counts the whole extraction pool. In a toy genome with
  two exact 10 kb copies, aligned with bwa-mem 0.7.19, all 2530 primary reads in a 5 kb
  footprint inside the first copy were MAPQ 0. No pool pair had a read there, while 2820 pool
  pairs had a read within 8 kb of it. Read from the code, not run: today spike would refuse
  such an event, because no breakpoint side has pool coverage (`src/simulate.rs:473-478`).
- **A read can have several placements inside one footprint.** In a toy genome with three
  copies of a 400 bp segment, 100 bp apart, plus one remote copy, bwa-mem put read
  `cluster_1` on the remote copy at MAPQ 0 with 3 `XA` entries, all three in the cluster.
- **spike's new reads are never flagged duplicate or QC-fail.** `align.sh` pipes the aligner
  into `samtools sort` (`src/main.rs:1046-1047`), and `merge.sh` only removes and merges.
- **The overlap guard compares footprint coordinates** (`src/main.rs:748-776`). Two events at
  distant look-alikes pass it.

## The physics, on paper (not measured)

Take a spot L with one perfect look-alike P, both diploid. The aligner cannot tell them
apart, so it places each read at L or P at random. Now make a het deletion of part of L.

| | depth at L | depth at P |
| --- | --- | --- |
| Real life | 3/4 | 3/4 |
| `--edit-model clean` (today) | 1 | 1 |
| `--min-mapq 0` | 1/2 | 1 |
| `--edit-model origin` | 3/4 | 3/4 |

Real life: 3 of the 4 copies are left, and their reads split evenly over L and P. The "random
placement" assumption is checked by the physics test below, not assumed in the result.

## Design

### The switch

- `--edit-model clean`: today's behavior. It is the default until the real-data test below
  supports `origin`.
- `--edit-model origin`: this design.
- It applies wherever spike removes reads today: DEL, INV, INS, small variants and the full
  DUP model. Additive events (fusion, `--dup-model junction`) remove nothing and are
  unchanged.

### Candidates at the spot

- **Which records.** Every primary record whose alignment overlaps the haplotype's reference
  footprint (event +/- `HAP_FLANK`, 2 kb). That includes any MAPQ, duplicates, QC-fail
  records, non-proper pairs, and pairs whose mate is elsewhere or unmapped.
  - Secondary and supplementary records are not judged on their own. They share their
    primary's name and go with it.
- **Placements.** A read's placements are its primary position and each `XA` entry.
- **Chance of the primary placement**, `p_here`, read from the record:
  - MAPQ of 1 or more: `1 - 10^(-MAPQ/10)`, the aligner's own estimate;
  - MAPQ 0 with an `XA` tag: `1 / (1 + number of XA entries)`;
  - MAPQ 0 without `XA`: `1/6`. bwa-mem lists alternatives in `XA` only when there are at
    most 5 (its `-h 5` default), so no `XA` means more than 5.
- **Chance of each `XA` placement:** `(1 - p_here) / k`, where `k` is the number of `XA`
  entries.
- **A read's chance of having come from the footprint** is the sum of the chances of all its
  placements inside the footprint, the primary included (R2). Examples:
  - MAPQ 0 at the spot, with 1 `XA` entry that is also inside the footprint: 1/2 + 1/2 = 1.
  - `cluster_1` above, if the cluster is the footprint: 0 + 3 x 1/4 = 3/4.
- **A pair's chance**, `p_origin`, is its surest mate's: the mate with the highest `p_here`.
  When both mates tie, the higher of their two sums. One mate pinned uniquely somewhere pins
  the fragment there.
- **Each fragment is judged once per event**, whether it was found at the spot or at a
  look-alike.

### Look-alike spots

- Gather every `XA` position named by candidates at the spot. Group positions within 1 kb of
  each other into regions (each grown by the read length).
- A region is a look-alike if at least 2 candidate reads point into it.
- Read the primary records in each look-alike region. One is a candidate if its `XA` names a
  position inside the spot's footprint.
  - Its chance is the rule above: the sum over its placements inside the footprint.
  - Example: MAPQ 0 with one `XA` entry, pointing at the spot, gives 1/2.

### Removal

- **An event's removal chance** for a candidate pair is `p_origin x copy_rate(copy, vaf)`.
  `copy_rate` is today's (`src/synth.rs:843`).
  - At the spot, `copy` is today's phase call from the sample's het SNPs.
  - At a look-alike, `copy` is unknown (`None`), so the rate is `vaf`.
- **Only pairs the new reads can replace are removable (R4).** Every mapped primary mate
  needs a placement fully inside the footprint: its primary, or an `XA` entry. An unmapped
  mate does not block. Any other candidate is kept.
  - Why: the new fragments never reach past the footprint. Removing a pair that sticks out
    would leave a depth dip where it sticks out.
  - This is today's Inside rule (`src/simulate.rs:185-189`), applied to every candidate, and
    through `XA` to reads placed at a look-alike.
- **Chances from different events add up (R5).** A fragment came from one copy, so "it came
  from event 1's edited copy" and "it came from event 2's" cannot both be true.
  - A pair's total chance is the sum over events, capped at 1.
  - All events' candidates are gathered first. Then each pair gets one draw against its
    total.
  - Example: het deletions at both copies of a perfect twin. Each event gives
    1/2 x 1/2 = 1/4, so the total is 1/2: 2 deleted copies of 4. Separate draws per event,
    as today, would give 1 - (3/4)^2 = 7/16.
- **Duplicate and QC-fail records (R3)** are removed by the same chance. A duplicate is a
  copy of one molecule, so it came from where that molecule came from.
  - A duplicate family shares one draw, compared with each member's own total. A family is
    the records with the same fragment ends: chromosome, and each mate's unclipped 5'
    position and strand.
  - They add nothing to the new reads' depth (below), because the new reads are never
    flagged.
- **Pool pairs** are re-emitted when kept, as today. Their chance is now `p_origin`, not 1.
  For MAPQ 60 the two agree to within 1e-6.
- **Non-pool candidates** are never re-emitted.
  - Removed: added to `replaced_reads.txt`. `merge.sh` removes them by name (both mates,
    every record), wherever they sit.
  - Kept: not touched at all; they stay in the BAM as they are.

### How many new reads

- Today the tiling depth is the pool's **fragment** coverage within 1 kb of a breakpoint
  (`estimate_coverage_at`, `src/simulate.rs:961`).
- `origin` builds it from where reads came from instead:
  - **Origin read depth.** For each candidate read, each of its placements inside the window
    adds that placement's chance: at the spot, and at look-alikes through their `XA`.
    Duplicate and QC-fail records add nothing (R3). This is in **read** units.
  - **Conversion to fragment units.** Multiply by `f`: the pool's summed fragment spans
    divided by its summed aligned read lengths, over the **whole pool** (R1).
    - It depends on fragment and read lengths, not on the spot.
    - The pool holds at least `MIN_DONOR_PAIRS` pairs, so `f` is never 0/0, even when no pool
      read lies in the footprint.
  - The result is in fragment coverage, the unit the tiling count already uses. This follows
    the case file's T3 rule: never compare numbers across estimators.
  - It never divides by the pool's depth at the spot, which is 0 inside a perfect twin.
  - In the twin case above, reads at L count 1/2 and so do their twins at P. The sum is L's
    true read depth.
- **The coverage refusal moves with it (R1).** Under `origin`, an event is refused for
  missing coverage only when no breakpoint side has origin depth above 0. Today's test is pool
  coverage (`src/simulate.rs:473-478`), which would refuse every event inside a perfect twin.
- The fragment-length and base-quality models are still learned from the clean pool only.

### What the user sees

- The log names each event's look-alike regions, the reads removed there, and the origin
  depth beside today's.
- `truth.vcf` keeps its fields. `SIM_RESIST` keeps its meaning (reads spike could not edit).
  Under `origin` that should be near 0, so RF8's refusal will rarely fire; that is expected.
- **Input check.** Under `origin`, if the footprint holds MAPQ 0 records and none of them
  carries `XA`, spike stops with an error that names the problem. Otherwise every MAPQ 0 read
  would silently get 1/6.
- **Known side effect, out of scope.** `spike validate`'s coverage rows expect a normal spot,
  where a het deletion halves the depth. At a look-alike the true drop is smaller (3/4 in the
  twin case), so they may flag correct events. That is RF14's kind of problem and is handled
  separately.

## Assumptions about the inputs, and how each is checked

| Assumption | How it is checked |
| --- | --- |
| The BAM keeps bwa-mem's `XA` tags | the input check above; counted on Monday's BAMs before any plan |
| MAPQ follows bwa-mem's meaning | same aligner family; the physics test measures it |
| The aligner places equally good hits at random | the physics test measures depth at L and P |
| The replacement reads' aligner (`align.sh`, bwa-mem2) acts like the one that made the donor BAM (the user's: bwa-mem) | reported on Monday; the physics test uses one aligner for both |
| Duplicates were marked by fragment ends (Picard or samtools markdup, no UMIs), so the family key above matches the marking | read from the BAM header's `@PG` lines on Monday |

## Testing

All of it follows the judgment gate: plan commit, then code, then result commit.

1. **Unit tests, written first**, each seen red:
   - `p_here` for MAPQ 60, MAPQ 3, MAPQ 0 with 1 `XA`, and MAPQ 0 without `XA`;
   - a pair's chance is its surest mate's;
   - a look-alike needs 2 pointing reads, and a read there counts only if its `XA` names
     the footprint;
   - the origin depth in the twin case, including a window where the pool has no read;
   - a read's chance sums its placements in the footprint: 3/4 for `cluster_1`'s case, and 1
     for MAPQ 0 at the spot with its one `XA` also inside (R2);
   - with no pool read in the footprint, `f` is finite and the event is not refused (R1);
   - a pair with a mate outside the footprint and no placement inside is kept, and an
     unmapped mate does not block (R4);
   - duplicate and QC-fail records add no depth, and a duplicate family gets one fate (R3);
   - two events at both copies of a twin give a removal chance of 1/2, not 7/16 (R5);
   - kept non-pool candidates are neither re-emitted nor listed for removal;
   - the missing-`XA` input check;
   - `clean` output is byte-identical to today's.
2. **The physics test** (before Monday, no real data).
   - The genome is a small made-up diploid one: unique sequence, a 10 kb segment S, more
     unique sequence, an exact copy of S, more unique sequence. The event is a het 1 kb
     deletion in the middle of the first S, so its 5 kb footprint lies inside S.
   - PREDICTED (from the toy genome above, not run): the pool holds no read in the
     footprint. That is the R1 case, so the test runs through it.
   - **Truth:** reads simulated from the genome with the deletion, aligned with bwa-mem2.
     Two seeds give the real-vs-real spread.
   - **Donor:** reads simulated from the genome without it, aligned. spike runs on it three
     ways: `clean`, `--min-mapq 0` and `origin`.
   - Measured: depth at L and P (any MAPQ), each as a ratio to the unique flanks.
   - **Pass:** `origin` falls inside the two truth seeds' spread, widened by a margin locked
     in the plan, at both L and P.
   - **The gate must reject something.** `--min-mapq 0` must fail (on paper it gives 1 at P,
     not 3/4). `clean` must fail or be refused; PREDICTED (not run): with no pool read at the
     breakpoints, today's code refuses it. If `--min-mapq 0` passes, the test cannot tell the models apart.
3. **The real test, on the user's GIAB BAMs (arriving 2026-09-28).**
   - Pick real deletions in hard spots (low MAPQ, look-alikes) that sample A carries and
     sample B lacks, from the GIAB truth sets.
   - Spike each into B three ways: `clean`, `--min-mapq 0` and `origin`.
   - Compare with A's real reads, at the spot and its look-alikes: depth at any MAPQ, MAPQ mix,
     split reads and discordant pairs.
   - Two real samples that both carry the same deletion set the real-vs-real spread.
   - **Pass:** `origin` lands inside the spread more often than both `clean` and
     `--min-mapq 0`. If `--min-mapq 0` does as well, keep the simpler one and stop.
   - Sites are fresh, and the thresholds are locked in a commit before any comparison is run.
   - Only a pass makes `origin` the default.

## Not in this design

- Changing `spike validate`'s expectations at look-alike spots.
- Regenerating whole regions from scratch (the rejected option 2): it would make the region
  all synthetic, and spike's synthetic reads are cleaner than real ones (RF15).
- Reads from look-alikes the aligner does not list (more than 5 alternatives). They get 1/6
  at most, which moves little.
- Making duplicates among the new reads (R3). In the flanks, duplicate-flagged depth falls
  by the removed share, because the new reads carry no duplicates. Unflagged depth is right,
  and callers that skip duplicates, as most do, see no change. Inside a deletion the fall is
  real: the deleted copy makes no reads at all.

## Review fixes (2026-09-27)

A review of `f5fa9ba` found five gaps. Each was checked against the code or the reviewer's
toy alignments before this revision.

| | Gap | Fix |
| --- | --- | --- |
| R1 | `f` was 0/0 when no pool read lies in the footprint; `MIN_DONOR_PAIRS` counts the whole pool | `f` over the whole pool; the coverage refusal uses origin depth |
| R2 | a read got one placement's chance when several lie in the footprint | sum over all placements in the footprint |
| R3 | duplicate and QC-fail depth became unflagged new reads | removed by origin, no depth added, one draw per duplicate family |
| R4 | non-pool pairs sticking out of the footprint were removed, leaving a dip | the Inside rule for every candidate, through `XA` at look-alikes |
| R5 | separate draws per event give 7/16 where the answer is 1/2 | chances add over events; one draw per pair |
