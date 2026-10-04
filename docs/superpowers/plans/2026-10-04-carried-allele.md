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

## K result (2026-10-04): PASS, so build

Probe `scripts/carried/probe.py` at the commit before this one. Output: `scratchpad/carried/k/` (`probe.log`, `sites.tsv`, `sites.vcf`). The run took 66 s.

| group | sites | checked | refused | share | unchecked | baseline refused | bound |
|---|---|---|---|---|---|---|---|
| K1 het SNV | 300 | 300 | 298 | 0.993 | 0 | 298/300 | >= 0.95 PASS |
| K1 hom SNV | 300 | 299 | 299 | 1.000 | 1 | 299/299 | >= 0.99 PASS |
| K1 het indel | 200 | 200 | 195 | 0.975 | 0 | - | >= 0.90 PASS |
| K1 hom indel | 200 | 200 | 197 | 0.985 | 0 | - | >= 0.95 PASS |
| K2a random SNV | 500 | 473 | 1 | 0.002 | 27 | 0/473 | <= 0.01 PASS |
| K2b random deletion | 250 | 238 | 0 | 0.000 | 12 | - | <= 0.01 PASS |
| K2b random insertion | 250 | 229 | 0 | 0.000 | 21 | - | <= 0.01 PASS |
| K2c HG001 SNV | 300 | 297 | 0 | 0.000 | 3 | 0/298 | <= 0.01 PASS |
| K2c HG001 indel | 200 | 200 | 1 | 0.005 | 0 | - | <= 0.01 PASS |

**The misses (seen, not judged).**
- 2 het SNVs, with K/N of 16/101 and 2/28. The baseline missed the same two.
- 8 indels: 7 insertions and 1 deletion. Six have K = 0, at a called het or hom site with 22-36 reads. A K of 0 at a hom call means the reads carry the gap somewhere other than the call's own boundary. That fits bwa placing an insertion in a repeat at another position than the VCF does. Not measured.

**The two K2 refusals.**
- `chr20:31166704 T>C`: 92 of 322 reads, about 10 times the usual depth.
- HG001's `chr20:271444 TACAC>T`: 9 of 42 reads.

Both sites already show another allele in this sample's reads.

**The baseline.** spike's pileup SNP rule, compared site by site on the SNVs, gives the same answer as the new rule on all but 2 of 1,400 sites:
- the K2a site above (new: refuse, baseline: pass);
- one K2c site that the new rule leaves unchecked (fewer than 10 reads span bases `s-1` to `e`) and the baseline passes.

So it would have been enough for SNVs, but it cannot see indels.

## Result (2026-10-04): supported

Code `a4b0c41`. Binaries: new release `61f2f01b`, new debug `86c2dd65`, master release `91f8fbab`. Runs: `scratchpad/carried/g_run.sh` (log `carried/g_run.log`) and `scratchpad/carried/g2/`.

**Tests.** 696 pass (11 new), 0 fail, 1 ignored (bcftools on the PATH), and the probe's 14 pass. The mutation runner checks the unmutated tests are green first, then catches 13 of 14 mutants: all 11 planned, plus suffix trimming and an N base counted.
- The survivor removes the call from `main`. No unit test reaches `main`, which is one function. G1 and G2 run the binary and do reach it.
- The plan asked for a test that sees the info line logged at a REF site. The test log capture keeps only warnings, so that line is seen in G1 and G3 instead.
- clippy gives the same warnings as master, after two of the new code's were fixed (a large enum variant, a complex test type).

**G1: PASS.** One run of the new binary with all 2,500 K sites (`--allow-overlap`) took 14.6 s and exited 1.
- It logged one line per site, 2,500 of 2,500.
- It listed 991 refused events, the 989 K1 and 2 K2 refusals the probe found.
- It warned at 64 sites, the probe's 64 unchecked.
- The counts (N, K) are identical at **2,500 of 2,500** sites, and so are the decisions.

**G2: PASS.** The reviewer's script stops at `existing_event` because its `run` raises on a non-zero exit. A scratch copy (`carried/g2/scripts/`) records the refusal instead.
- `existing_event`: exit 1, `SNV  chr1:5001 T>A: 38 of 38 reads carry another allele here`, then `Error: the sample already carries another allele at 1 small variant(s): SNV  chr1:5001 T>A (38 of 38 reads). ...`
- Every other result is as before:
  - `hom_only` 5/0 and `hom_and_het` 0/7;
  - `asymmetric_mates` 150/150;
  - `low_fraction` 2 pairs, 0 ALT, `SIM_VAF=0.014`;
  - `indel_rate` NaN, 2.0 and -0.5 all exit 0;
  - `true_del_validate` exit 1, base depth 0;
  - clean tiling 148 reads and origin 111;
  - the two FASTQ probes refused, as after the fastq-safety fix.

**G3: PASS.** The new and master binaries gave identical `truth.vcf`, `R1.fq.gz`, `R2.fq.gz`, `replaced_reads.txt` and `fastq_removed_reads.txt`:
- for F1's three structural events, which logged no check line;
- for `snp:chr20:51218639:T:C;af=0.5`, the first K2a site with N of at least 25, which logged `0 of 33 reads carry another allele here`.
