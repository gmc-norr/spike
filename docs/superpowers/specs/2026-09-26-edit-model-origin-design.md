# `--edit-model origin`: edit reads by where they came from

Date: 2026-09-26. Status: design agreed with the user in four parts; not built.

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
- **Chance the read came from here**, `p_here`, read from the record:
  - MAPQ of 1 or more: `1 - 10^(-MAPQ/10)`, the aligner's own estimate;
  - MAPQ 0 with an `XA` tag: `1 / (1 + number of XA entries)`;
  - MAPQ 0 without `XA`: `1/6`. bwa-mem lists alternatives in `XA` only when there are at
    most 5 (its `-h 5` default), so no `XA` means more than 5.
- **A pair's chance** is the highest `p_here` among its mates placed at the spot. One mate
  pinned uniquely here pins the fragment here.

### Look-alike spots

- Gather every `XA` position named by candidates at the spot. Group positions within 1 kb of
  each other into regions (each grown by the read length).
- A region is a look-alike if at least 2 candidate reads point into it.
- Read the primary records in each look-alike region. One is a candidate if its `XA` names a
  position inside the spot's footprint.
  - Its chance of having come from the spot is `(1 - p_there) / k`. Here `p_there` is the rule
    above applied to its own placement, and `k` is its number of `XA` entries.
  - Example: MAPQ 0 with one `XA` entry pointing at the spot gives 1/2.

### Removal

- Each candidate pair is removed with probability `p_origin x copy_rate(copy, vaf)`.
  `copy_rate` is today's (`src/synth.rs:843`).
  - At the spot, `copy` is today's phase call from the sample's het SNPs.
  - At a look-alike, `copy` is unknown (`None`), so the rate is `vaf`.
- **Pool pairs** keep today's handling: they are re-emitted when kept, and suppressed only
  when the pair lies Inside the footprint. Their chance is now `p_origin`, not 1. For MAPQ 60
  the two agree to within 1e-6.
- **Non-pool candidates** are never re-emitted.
  - Removed: added to `replaced_reads.txt`. `merge.sh` removes them by name (both mates,
    every record), wherever they sit.
  - Kept: not touched at all; they stay in the BAM as they are.
- A read removed by any event is removed. Today's union of suppressed names across events
  covers this.

### How many new reads

- Today the tiling depth is the pool's **fragment** coverage within 1 kb of a breakpoint
  (`estimate_coverage_at`, `src/simulate.rs:961`).
- `origin` builds it from where reads came from instead:
  - **Origin read depth.** The sum of `p_origin` over candidate reads covering that window: at
    the spot, and at look-alikes through their `XA`. This is in **read** units.
  - **Conversion to fragment units.** Multiply by `f`, the pool's fragment coverage divided by
    its read coverage, both measured over the whole footprint.
  - The result is in fragment coverage, the unit the tiling count already uses. This follows
    the case file's T3 rule: never compare numbers across estimators.
  - It does not divide by the pool's depth in the window, which can be 0 at a perfect twin.
    A pool empty over the whole footprint is already refused (`MIN_DONOR_PAIRS`), since the
    quality and fragment models need it.
  - In the twin case above, reads at L count 1/2 and so do their twins at P. The sum is L's
    true read depth.
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

## Testing

All of it follows the judgment gate: plan commit, then code, then result commit.

1. **Unit tests, written first**, each seen red:
   - `p_here` for MAPQ 60, MAPQ 3, MAPQ 0 with 1 `XA`, and MAPQ 0 without `XA`;
   - a pair's chance is its surest mate's;
   - a look-alike needs 2 pointing reads, and a read there counts only if its `XA` names
     the footprint;
   - the origin depth in the twin case, including a window where the pool has no read;
   - kept non-pool candidates are neither re-emitted nor listed for removal;
   - the missing-`XA` input check;
   - `clean` output is byte-identical to today's.
2. **The physics test** (before Monday, no real data).
   - The genome is a small made-up diploid one: unique sequence, a 5 kb segment S, more
     unique sequence, an exact copy of S, more unique sequence. The event is a het 1 kb
     deletion inside the first S.
   - **Truth:** reads simulated from the genome with the deletion, aligned with bwa-mem2.
     Two seeds give the real-vs-real spread.
   - **Donor:** reads simulated from the genome without it, aligned. spike runs on it three
     ways: `clean`, `--min-mapq 0` and `origin`.
   - Measured: depth at L and P (any MAPQ), each as a ratio to the unique flanks.
   - **Pass:** `origin` falls inside the two truth seeds' spread, widened by a margin locked
     in the plan, at both L and P. `clean` must fail, or the test cannot tell them apart.
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
