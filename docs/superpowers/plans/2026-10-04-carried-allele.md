# Refuse a small variant the sample already carries

**Asked 2026-10-04.** This is review finding 3 (P1, `docs/review/2026-10-04-independent-review.md`), the user's second pick. We reproduced it on master `2764dbc`.
- The donor was hom-alt for `chr1:5001 T>A`, and the same variant was requested at `af=0.5`.
- spike exited 0 and wrote `SIM_VAF=0.500;GT=0/1`.
- The output reads held **25 ALT and 0 REF**.

## Why it happens (read in the code)

- spike's REF check compares against the reference FASTA, never against the sample.
- The sample's alleles at the event's own bases are left out of both copies (`own_bases`, `src/simulate.rs:461`).
- For an event at AF `v`, spike replaces reads by copy (`synth::copy_rate`). The copies come from the het SNPs that phase each fragment, and a coin per phase block picks the event copy.
  - The event copy's reads go at rate `min(1, 2v)`, and the other copy's at `max(0, 2v - 1)`.
  - A fragment that no het SNP places goes at rate `v`.
  - So at `v = 0.5`, the other copy's reads stay as they were.
- So when the sample already has a non-REF allele at the event's bases, the reads left behind still show it:
  - **Hom-alt:** both copies carry it, so the result is 100% ALT.
  - **Het with the same ALT:** the result depends on the coin. It is 0.5 when the event lands on the ALT copy, and 1.0 when it lands on the REF copy.
  - **Het with another ALT:** the site ends up with two ALTs, and `GT=0/1` is wrong.
- spike's own sample model (`loh::sample_copies`) holds SNPs only. A pileup calls them at depth 10 or more: het when the top two bases are each 20-80%, hom-alt at 90% or more. So it cannot see an indel the sample has.

## Gate A

1. **Principle.** spike's truth must describe what its reads hold. When the sample already carries another allele at the bases a small variant changes, no edit spike makes gives the requested fraction, so spike refuses that event.
2. **What would kill it.** The kill test K below:
   - the rule misses sites the sample carries (K1);
   - or it refuses sites the sample does not carry, including real variants of another sample (K2).
3. **Refuted before?** No. `git log` has no commit on an existing or carried allele. The case file's RF8 entry applies: a noise test is drawn from what the rule will actually see. That is why K2(c) uses another sample's real variants.
4. **The simplest thing.** Count reads at the event's own bases, before any other work, and refuse past one threshold. The trivial baseline is spike's existing SNP rule at the event position, which covers SNVs only. K reports it beside the new rule.
5. **Inputs.**
   - The hospital BAM is the sample at the hospital: HG002 SeraCare, 30x, raredisease, bwa-mem2.
   - Its own DeepVariant calls (`GM24385_seracare_cancer_snv.vcf.gz`) are the sites it carries. FILTER is `.` on every record.
   - GIAB HG001 v4.2.1 is the other sample's variants.
   - The read filter is spike's pileup filter (`loh::count_alleles`): primary, not duplicate, not QC-fail, MAPQ at least `--min-mapq` (20).

## Design (locked)

**The site.** REF and ALT are trimmed of their common prefix, then their common suffix. That leaves the changed reference bases `[s, e)` (0-based). A pure insertion has `e = s`: its sequence goes between bases `s-1` and `s`.

**A read covers the site** when its aligned reference span (from its start through its last M, =, X, D or N) includes bases `s-1` through `e`.

**A covering read carries another allele** when any of these holds:
- an aligned base (M, =, X) at a position in `[s, e)` differs from the reference there (a read base N does not count);
- a D or N op overlaps `[s, max(e, s+1))`;
- an I op sits at a boundary `b` with `s <= b <= e` (the insertion goes between bases `b-1` and `b`).

**The rule.** The reads are the pileup's: primary, not duplicate, not QC-fail, MAPQ at least `--min-mapq`.
- With N covering reads and K carrying another allele, spike refuses the event when N >= 10 and K/N >= 0.2.
- Both numbers are the pileup's own: it needs 10 reads, and 20% is its lower bound for a het allele. They are in the same units: a share of reads, under the same filter.
- When N < 10, spike logs a warning that the site could not be checked, and goes on.

**Where.** Small variants only (`SimEvent::SmallVariant`: SNVs, MNVs and small indels given as REF/ALT). The check runs once per such event, after the overlap check and before any extraction.
- Each site logs one info line: `<chrom>:<pos> <REF>><ALT>: K of N reads carry another allele here`.
- If any event is refused, spike stops with one error that lists every refused event with its K and N. The error says the sample already carries another allele there, so no edit can give the requested fraction, and to remove those events. There is no override flag.
- Structural events (`del:`, `ins:`, `dup:`, `inv:`, fusions) are not checked. That is a limit of this fix, and the README says so.

**Tests, written first and seen red.**
- The trimmed site for an SNV, an MNV, a deletion, an insertion, and a REF/ALT with a shared suffix.
- Read classes, on hand-made records:
  - a matching read;
  - a mismatch inside the site;
  - a mismatch outside it;
  - a deletion overlapping it;
  - an insertion at each boundary, and one past it;
  - a read that stops inside the site (not covering);
  - an N base.
- On a small BAM:
  - hom-alt at the site: refused, with the count in the message;
  - het: refused;
  - REF: not refused, and the info line logged;
  - 9 covering reads: not refused, and warned.
- A run whose VCF holds one carried and one REF event stops before extraction, and its message names only the carried one.

**Mutation checks** (each must turn a named test red; the unmutated tests must be green first):
1. no trimming;
2. mismatches outside the site counted;
3. deletions ignored;
4. insertions ignored;
5. the insertion boundary as `s < b < e`;
6. coverage not required;
7. the threshold 0.2 changed to 0.5;
8. the depth floor 10 changed to 1;
9. duplicates counted;
10. MAPQ not filtered;
11. only the first refused event reported.

## Kill test K (locked; run before any Rust code)

A Python probe (`scripts/carried/probe.py`, pysam) applies the rule above to the hospital BAM.
- **The BAM:** `D24-14230_Seq25-7600_30x_sorted_md.bam`.
- **The reference:** GRCh38 no-alt.
- **Every draw:** chr20, seed 1, single-ALT records only. No two sites within 50 bp of each other.

**K1: sites the sample carries.** From its own calls, with GQ at least 20, draw 300 het SNVs, 300 hom-alt SNVs, 200 het indels and 200 hom-alt indels. Each is requested as called (same REF/ALT). Het means `0/1` or `1/0`.

**K2: sites it does not carry.** Each site has no call of the sample within 100 bp.
- (a) 500 random SNVs at random non-N positions, with a random ALT;
- (b) 250 random deletions and 250 random insertions of 1-10 bp. A deletion's REF is the anchor plus the next 1-10 bases. An insertion's sequence is random;
- (c) GIAB HG001 v4.2.1 records on chr20 inside its benchmark bed: 300 SNVs and 200 indels of up to 50 bp.

**Reported per group:** the sites checked (N >= 10), the share of those refused, and the share not checked (N < 10). For SNV groups, also the baseline's refusals: spike's pileup SNP rule at the event's position.

**Pass (all must hold, on checked sites):**
- K1: refused share at least 0.95 for het SNVs, 0.99 for hom-alt SNVs, 0.90 for het indels and 0.95 for hom-alt indels;
- K2: refused share at most 0.01 in each of (a), (b) and (c).

Every K1 miss and every K2 refusal is listed with its N, K and reads, so a cause can be named. If K fails, nothing is built. The result is reported, and the design goes back to the user.

## Checks after building (locked)

**G1: spike counts what the probe counts.** One spike run whose VCF holds every K site (they do not overlap) logs one line per site, then stops, because K1 sites are carried.
- For every site, spike's N and K equal the probe's. Pass if all decisions agree and at least 99% of sites have identical counts.
- Every difference is listed and explained.

**G2: the reviewer's probe.** `scripts/review_20261004.py` (the scratch copy from the fastq-safety result) against the new debug binary:
- `existing_event` exits non-zero, and its message says another allele is carried;
- every other result is as before.

**G3: nothing else changes.** On the hospital BAM, with `--seed 1 --threads 16`, the new and master binaries must give identical `truth.vcf`, `R1.fq.gz`, `R2.fq.gz`, `replaced_reads.txt` and `fastq_removed_reads.txt`:
- for F1's three structural events;
- for one SNV at a K2(a) site.
