# Read length: sequence and trim spike's reads the way the real library's were

Locked 2026-09-28, before any code, with the user's agreement ("1": copy the
sequencer's trimming rather than draw lengths at random). Branch `read-length`,
off master `c87ab53`.

## The gap

spike gives every synthetic read one length, the input BAM's mean read length
rounded (`src/main.rs`, `bam_stats.read_length.round()`). On the HG001 and
HG002 30x BAMs that is 148 bp. It also never makes a fragment shorter than
that length: `FragmentDist::from_read_pairs` drops them, and every
`sample_in_range` call starts at the read length.

The real reads are measured below (all numbers from the HG001 and HG002 30x
BAMs of transplant rounds 1-3; scratchpad `readlen/`).

## What the real reads show

**Length shares** (HG002, primary, proper, not duplicate):

| Length | chr1:50-52 Mb (482,603 reads) | chr20:10-12 Mb (462,818 reads) |
|---|---|---|
| 151 | 0.607 | 0.605 |
| 150 | 0.269 | 0.273 |
| 149 | 0.069 | 0.067 |
| 148 | 0.022 | 0.023 |
| 147 | 0.005 | 0.006 |

**The 3' end of a read is cut where it spells the adapter's start.** For
reads in fragments of 170 bp and up (fully matched, NM 0), the reference bases
right past an L-bp read's 3' end (read orientation) are the first 151-L bases
of `AGATCGGAAGAGC`, the TruSeq adapter's start:

| | 150 bp | 149 bp | 148 bp | 147 bp |
|---|---|---|---|---|
| HG002 chr20:10-10.6 Mb | 0.997 | 0.996 | 0.993 | 0.993 |
| HG002 chr1:50-50.6 Mb | 0.997 | 0.999 | 0.997 | 0.997 |
| HG001 chr20:10-10.6 Mb | 0.998 | 0.997 | 0.998 | 1.000 |
| HG001 chr1:50-50.6 Mb | 0.997 | 0.998 | 0.996 | 1.000 |

In the same 4 sample-regions, **0.000** of the 151-bp reads end in `A`.
- The cut is the *longest* read suffix equal to an adapter prefix. A 148-bp
  read is followed by `AGA`, not just `A`.
- It happens before the FASTQ, on the sequencer's side. The pipeline's own
  fastp report for HG002
  (`raredisease_results/trimming/*_LNUMBER1.fastp.json`) gives a mean read
  length of 149 **before** fastp. fastp cut adapters from 5,690,224 of
  760,000,000 reads (0.75%).

**Short fragments give reads as long as the fragment:**
- 2.9% of pairs have a fragment under 151 bp: 7,037 of 241,277 in chr1:50-52
  Mb, and 6,656 of 231,383 in chr20:10-12 Mb.
- Of the reads in them, 95.4% and 96.1% are exactly the fragment's length.
- Of the reads under 140 bp, 99% sit in such fragments: all but 115 of
  11,387, and all but 20 of 10,809.

**Mates are cut independently.** In fragments of 151 bp and up, the joint
R1-by-R2 length table matches the product of its margins. For example, 151/151
is 90,780 observed against 91,534 expected (chr1), and 86,238 against 87,002
(chr20).

**An untrimmed BAM looks different:**
- `data/validation/hg002_novaseq_chr20.bam` holds 130,971 reads in
  chr20:10-10.5 Mb. All of them are 151 bp, and 0.300 end in `A`.
- So the signature separates a trimmed library from an untrimmed one: 0.000
  against 0.300.

## Gate A

1. **Principle.** spike's reads should be sequenced and trimmed the way the
   input library's reads were, so their lengths come out the way the real ones
   do, and from the same causes.
2. **What would kill it:** K1 below fails. That would mean spike's reads do not
   show the real length mix, the real sequence dependence, or the real share of
   short fragments.
3. **Already refuted here?** No. The history holds only L5 (`158356e`,
   `7c680f4`): reject a mean read length above 1,500 bp. That stays.
   "Fixed-length" was a design choice, not a measured one.
4. **Simplest version.** The alternative is to draw each read's length from the
   real mix. It is rejected by the measurement above, because length depends on
   the sequence (a 151-bp real read never ends in `A`). The trimming rule is no
   bigger: one suffix match.
5. **Inputs, and how each was checked:** the tables above, in two samples and
   two regions each. The whole-file head that spike samples is checked in K0.

## Design

**What spike learns from the input BAM** (`bam_stats`, the first 50,000
primary records, as today):
- `cycles`: the most common read length. It replaces the rounded mean.
  - On an untrimmed library this is the one length every read has.
  - It is 151 in the regions measured above; K0 checks the head spike reads.
- `adapter_trimmed`: at least 1,000 of the sampled reads are `cycles` long, and
  under 0.05 of those end in `A` at their 3' end, in read orientation. Measured:
  0.000 trimmed, 0.300 untrimmed.

**When `adapter_trimmed` is false:** spike behaves as today, with `cycles` as
the read length. Every read is that long, and fragments shorter than it are not
drawn.

**When `adapter_trimmed` is true:**
- **Fragments:** the fragment model keeps every positive insert size up to
  1,500 bp (`FragmentDist::from_read_pairs` with a minimum of 1), and every
  draw uses that minimum. The tiling count's mean fragment length then
  describes the fragments actually drawn (CR-FRAG).
- **Each mate** is sequenced for min(`cycles`, fragment length) bases from its
  end of the fragment. A short fragment gives two reads of the fragment's
  length, as real trimming does.
- **When the fragment is at least `cycles` long,** each mate's 3' end (FASTQ
  orientation, called bases with their errors) loses its longest suffix that
  equals a prefix of `AGATCGGAAGAGC`, with the same qualities cut.
- Everything else is as today:
  - positions, orientation and the R1/R2 coin;
  - the quality model, learned per cycle up to `cycles`; shorter pool reads
    fill their own first cycles;
  - the random stream stays in a fixed order, so output is the same for the
    same `--seed` and the same at any `--threads`.

The same `--seed` will give different reads than before this change, in every
mode.

**Where it lands:**
- `bam_stats.rs`: learn `cycles` and `adapter_trimmed`.
- `main.rs`: SimConfig, the generator, `finish_donor_pool`'s minimum,
  `validate_read_length`.
- `stats.rs`: nothing new; callers pass the minimum.
- `synth.rs`: the trim rule, and the pair generators (`generate_read_pair`,
  `generate_haplotype_read_pair`, `generate_depth_pair`).
- `simulate.rs`: `tile_haplotype_reads`'s minimum fragment.
- The README's "fixed-length" statements.

## Checks, locked

**K0 (the code):**
- Every test passes. New tests are written first and watched fail.
- A mutation run on the new code: every mutant must go red.
- `spike` logs `cycles` and `adapter_trimmed` for three BAMs:
  - the HG001 30x BAM, which must read trimmed, 151;
  - the HG002 30x BAM, which must read trimmed, 151;
  - `data/validation/hg002_novaseq_chr20.bam`, which must read **not**
    trimmed, 151. It is the case the check must reject.

**K1 (it does what it says).** spike runs round 2b's forward normal events
(HG002 recipient, 3 runs) with the new binary. Its synthetic reads (names
`ev*`, from `R1.fq.gz` and `R2.fq.gz`) are compared with the recipient BAM's
reads over the same events' windows (primary, proper, not duplicate):
- a. The share at each of 151, 150, 149 and 148 bp, and under 148 bp, lies
  within 2 percentage points of the recipient's.
- b. Under 1% of synthetic 151-bp reads end in `A` (3' end, FASTQ
  orientation).
- c. The share of synthetic pairs whose aligned fragment is under 151 bp
  (`sim.bam`, proper pairs, |TLEN|) lies within 1 point of the recipient's.
- d. Of synthetic reads under 140 bp in proper pairs, at least 90% are exactly
  |TLEN| long. Measured real: 95.4-96.1% for all reads in such fragments.

K1 passes only if a-d all hold.

**K2 (no harm).** Transplant round 2b is rerun with the new binary: the same
events (`round2.py full`, seed 2), everything else unchanged. Then
`round3.py` is run on that output.
- a. Under round 2b's rule, no group and direction goes from supported (round
  2b) to refuted.
- b. Under round 3's rule, no normal-arm metric fails.
- Seen, not judged: DUP50-299 and INS20-49 medians and widths, and R_f,
  against round 2b and round 3.

The result is **supported** if K0, K1 and K2 all pass, **refuted** if K1
fails, and **needs a look** if K2 fails with K1 passing.

## Order

1. This plan (`plan:` commit).
2. Code, test-first (`code:` commits). Separate `CARGO_TARGET_DIR`; binary
   md5 recorded.
3. K0, K1 and K2 runs, at most 16 threads, logs kept and checked for exit
   status and last line.
4. `result:` commit. Merge and push only on the user's word.

## Not in this round

- Other adapters (Nextera and others): such a library reads as not trimmed and
  keeps today's behaviour.
- Untrimmed libraries' adapter read-through in short fragments.
- fastp's own changes: pair-overlap correction, and dropping the low-quality
  2.6%.
- Illumina-style names for synthetic reads (needed for the hospital FASTQ, and
  a separate change).
