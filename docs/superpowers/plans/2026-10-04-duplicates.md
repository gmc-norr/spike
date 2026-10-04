# Duplicates: do spiked sites look different after raredisease marks duplicates again?

**Asked 2026-10-04.** The user picked option 1: measure the duplicates gap before changing anything. Real data has duplicate reads, but spike makes none and never takes one.

## Measured before this plan

- **Real data is 4.7% duplicates.** In the hospital 30x BAM's Picard metrics (`alignment/*.MarkDuplicates.metrics.txt`), `PERCENT_DUPLICATION` is `0.047185`. Of its 17,385,340 duplicate pairs, 10,361,338 are optical. In the window the stand-in FASTQ covers (chr20:9,000,000-13,000,000), 40,129 of 965,254 primary mapped records carry the duplicate flag, or 4.16% (`samtools view -c -F 0x904`, with and without `-f 0x400`).
- **spike never takes a duplicate.** `passes_filters_bam` skips any record with the duplicate flag (`src/extract.rs:837`). So `replaced_reads.txt` never lists a duplicate, and neither does `fastq_removed_reads.txt`, which is a subset of it. `merge.sh` keeps them on purpose (the comment at `src/main.rs:2060`).
- **What raredisease runs** (from the same run's `pipeline_info` and BAM header):
  - fastp 0.23.4 with `--detect_adapter_for_pe --length_required 40 --correction --overrepresentation_analysis`, so no `--dedup`;
  - bwa-mem2 2.2.1 `mem -M -K 100000000`;
  - Picard MarkDuplicates 3.3.0, default settings plus `--MAX_SEQUENCES_FOR_DISK_READ_ENDS_MAP 50000`, `--OPTICAL_DUPLICATE_PIXEL_DISTANCE 100`;
  - DeepVariant 1.6.1 for SNVs.
- **The concern** (read from the code, not yet measured):
  - In the BAM route, the duplicates keep their flag and callers skip them.
  - In the FASTQ route, raredisease marks duplicates again from scratch. A duplicate is a copy of one fragment. If spike removed the copy Picard had kept, a leftover copy has nothing to be a copy of. It then counts as an ordinary read carrying the sample's old allele.
- **PREDICTED (not run):** about 3.5-4% of reads at a spiked site leak back this way. A hom SNV (1.0) would read about 0.96; a het SNV (0.5) about 0.49.
- **Local tools:** bwa-mem2 2.2.1 (the same version), Picard 3.3.0 (`scratchpad/names/picard-3.3.0.jar`), pysam 0.23.3. fastp here is 0.22.0, not 0.23.4. It is left out because raredisease's fastp does not deduplicate, and the stand-in's reads already went through fastp once.
- **The binary:** spike `8ae6e449`, built from the code at `master` `f15bab8`. This plan changes no spike code.
- **Real calls in the window.** In raredisease's SNV VCF over chr20:9,500,000-12,500,000, 1,593 records are 1-bp-for-1-bp with GT `1/1` and 2,727 with GT `0/1`. Every record's FILTER is `.`.

## Checks (locked before running)

**Allele counts.** At a 1-based position, count each read that is primary and mapped, has no duplicate, QC-fail, secondary or supplementary flag, has MAPQ ≥ 5, and has a base there (not a deletion) with base quality ≥ 10. The MAPQ and base-quality floors are DeepVariant 1.6.1's documented `make_examples` defaults, not run here. A read is `ref` when its base equals the reference base there.

**The events** (`draw_events.py`, seed 7):
- Draw positions uniformly from chr20:9,500,000-12,500,000. Accept one when:
  - its reference base is A, C, G or T, and there is no N within 1,000 bp;
  - raredisease's VCF has no record within 1,000 bp;
  - in the source BAM, 20 to 45 primary, mapped, non-duplicate, non-QC-fail records cover it, and at least 95% of them have MAPQ ≥ 20 and the proper-pair flag;
  - it is at least 20,000 bp from every event already accepted.
- The first 20 accepted positions are het SNVs (`SIM_VAF=0.5`), the next 20 are hom SNVs (`SIM_VAF=1.0`). ALT is the transition (A↔G, C↔T).
- Then 3 hom deletions of 300 bp (`SVTYPE=DEL`, `SIM_VAF=1.0`). A deletion is accepted when every 50th base of `[pos, pos + 300)` passes the tests above, and its span is 1,000 bp from any VCF record and 20,000 bp from every event.
- At most 20,000 candidates. If fewer are accepted, stop before running.

**The runs.**
- **Spiked:** `spike --bam <hospital 30x BAM> --reference <GRCh38 no-alt> --vcf events.vcf --seed 1 --threads 16 --align --raw-fastq <stand-in R1> <stand-in R2>`. The stand-in is F1's (`scratchpad/fastq/standin/main`, 478,815 pairs).
- **Baseline:** the stand-in itself.
- Each FASTQ pair is aligned with `bwa-mem2 mem -M -K 100000000 -t 16 -R '@RG\tID:sim\tPL:ILLUMINA\tSM:sim'`, sorted with `samtools sort`, then marked with Picard 3.3.0 using raredisease's options above.
- **Fallbacks.**
  - If spike refuses events (RF8), drop the events it names and run again, once.
  - If `fastq.sh` refuses (a removed original outside the stand-in), rebuild the stand-in from chr20:8,000,000-14,000,000 the way F1 did, once.

**D (judged): do hom SNVs leak the old allele?**
- `R_spk`: `ref` reads over all counted reads, pooled over the 20 spiked hom SNVs, in the spiked BAM.
- `R_real`: the same over real hom SNVs in the baseline BAM. These are raredisease's GT `1/1` 1-bp-for-1-bp records in chr20:9,500,000-12,500,000 that pass the same coverage test as the events, and share no position with another record.
- `L`: of the `ref` reads in `R_spk`, the share whose name carried the duplicate flag in the source BAM.
- **Matters:** `R_spk − R_real > 0.02` and `L ≥ 0.5`.
- **Does not matter:** `R_spk − R_real ≤ 0.02`. Log it and stop.
- **Something else:** `R_spk − R_real > 0.02` but `L < 0.5`. The extra old-allele reads come from something other than duplicates; look at them before anything else.
- **Why 0.02:** half the predicted leak. A leak that small is under one read per hom site at 30x.

**Controls. The run counts only if all three pass.**
- **C1, the leak is seen when it is known to be there.**
  - Mark the baseline's aligned BAM again with `--TAG_DUPLICATE_SET_MEMBERS true`.
  - With seed 7, pick 200 duplicate sets (one `DI` value) with exactly one non-duplicate pair, at least one duplicate pair, and only proper, primary pairs.
  - Drop the non-duplicate pair of the first 100 sets from the stand-in by name. Align it and mark it the same way, tagging on.
  - **Pass:** at least 95 of those 100 sets have a member pair that is no longer flagged duplicate, and at most 5 of the other 100 do.
- **C2, the counter reads the right base.**
  - **Pass:** `R_real ≤ 0.10`.
  - **Pass:** at the same real sites shifted +1 bp (dropping any that falls on a VCF record), the pooled `ref` share is ≥ 0.90.
- **C3, enough reads.** **Pass:** at least 15 of the 20 spiked hom SNVs and at least 100 real hom SNVs have ≥ 15 counted reads.

If a control fails, the result is inconclusive.

**Reported, not judged:**
- the `ref` reads at spiked hom SNVs, split four ways:
  - source duplicate;
  - spike did not take it (not `SPIKE_`, not a source duplicate, not in `replaced_reads.txt`);
  - in `replaced_reads.txt`;
  - `SPIKE_`;
- the same pooled `ref` share for the BAM route: source BAM reads not named in `replaced_reads.txt`, plus `sim.bam`;
- the pooled alt share at spiked het SNVs against real GT `0/1` sites (same tests) in the baseline;
- the share of duplicate-flagged primary records within 300 bp of each spiked SNV, spiked BAM against baseline BAM;
- for each hom deletion, the counted-filter reads (MAPQ ≥ 5, not duplicate) that lie wholly inside it, in the spiked and baseline BAMs, and how many were source duplicates;
- Picard's `PERCENT_DUPLICATION`, spiked against baseline.

**Code** (`scripts/duplicates/`, tests first, seen red):
- `draw_events.py`;
- `measure.py`: the counter, the D rule, C2, C3 and the reported rows;
- `control.py`: C1's set picking and counting;
- `run.sh`: the steps above.

All paths are arguments, so no local path goes into the repo.

**Not run here:**
- raredisease itself;
- fastp;
- DeepVariant (its read filters are copied, not run);
- `--edit-model origin`, which already gives a duplicate its original's fate (`src/origin.rs` `test_a_duplicate_shares_its_originals_fate`).

## Result (2026-10-04): supported. The gap matters.

**Setup.**
- spike `8ae6e449`, scripts at code commit `3c462dc`.
- The scripts have 24 tests, and 23 of 23 mutants are caught. The first round caught 19; four tests were added or fixed to catch the other four:
  - a site near a call or an N;
  - a deletion over a low-MAPQ stretch;
  - the margin boundary, which floats had hidden.
- The run took 4 min wall, with an 18 GB maximum resident size.
- Logs: `scratchpad/dups/{go.log,go.err,run/}`.

**The run.**
- spike took all 43 events (no fallback) and exited 0.
- `replaced_reads.txt` holds 100,938 pairs, and `fastq_removed_reads.txt` holds 13,920.
- `fastq.sh`: `removed 13920 original pairs, added 13495 of spike's pairs`.
- `truth.vcf` has 43 records.

**D: matters.**
- `R_spk` = 23/676 = **0.0340**, over the 20 spiked hom SNVs in the spiked BAM.
- `R_real` = 81/47,597 = **0.0017**, over 1,435 real hom SNVs in the baseline BAM.
- `R_spk − R_real` = **0.0323**, above the locked 0.02.
- `L` = 23/23 = **1.00**: every old-allele read at a spiked hom SNV carried the duplicate flag in the source BAM.
  - The split is 23 source duplicates, 0 not taken, 0 in `replaced_reads.txt`, 0 `SPIKE_`.
  - Per site, 0-3 old-allele reads; 14 of the 20 sites have at least one.

**Controls: all pass.**
- **C1:** 100 of 100 sets with the kept pair removed have a member that is no longer flagged duplicate, against 0 of 100 untouched sets.
- **C2:**
  - `R_real` is 0.0017.
  - At +1 bp, the `ref` share is 47,234/47,245 = 0.9998, over 1,421 sites.
- **C3:** 20 of 20 spiked hom SNVs and all 1,435 real ones have ≥ 15 reads.

**Reported, not judged:**
- **The BAM route has no leak.** At hom SNVs it reads 7/660 = 0.0106 old allele.
  - All 7 are reads spike cannot take: pairs whose mate maps to another chromosome (chr2, chr6, chr8, chr15, chrM).
  - The stand-in FASTQ lacks them, because they have no mate in its window. So the FASTQ route above shows 0 of them.
  - A real full FASTQ holds them, so there the hom share would be about 1 point higher on top of the leak. This is the `SIM_RESIST` share, not duplicates.
- **Het SNVs:** spiked alt share 314/678 = 0.4631, against real 40,113/80,802 = 0.4964 over 2,438 sites.
- **Duplicates within 300 bp of the 40 spiked SNVs:** 60/6,897 = 0.0087 in the spiked BAM, against 253/7,025 = 0.0360 in the baseline. The spiked regions lose about three quarters of their duplicates.
- **Hom deletions:** counted reads wholly inside each one.

  | Deletion | Spiked | of which source duplicates | Baseline |
  |---|---|---|---|
  | chr20:10856999-10857298 | 2 | 2 | 34 |
  | chr20:11012916-11013215 | 3 | 3 | 38 |
  | chr20:12334123-12334422 | 3 | 3 | 37 |

- **Picard `PERCENT_DUPLICATION`:** spiked 0.040578, baseline 0.041689.

**What it means.**
- After raredisease's own duplicate marking, a spiked hom SNV reads about 3.4% old allele where real hom SNVs read 0.17%.
- A spiked hom deletion keeps 2-3 reads inside where a real one keeps none.
- Every such read is a source duplicate that spike left in the FASTQ after removing the copy Picard had kept.
- The likely fix is the simplest one: spike takes duplicates too, and each one gets its kept copy's fate. `--edit-model origin` does this already. It is not built yet; that waits for the user's word.

**Not run here:**
- raredisease, fastp and DeepVariant themselves;
- a real full FASTQ, which also holds the cross-chromosome pairs above;
- `--edit-model origin`.

## Fix: a duplicate shares its kept copy's fate (plan, 2026-10-04)

The user picked option 1: fix the default (`--edit-model clean`), then re-run the check above as the after-test.

**Gate A.**
1. **Principle.** A duplicate is another read of the same molecule, so it goes wherever that molecule goes.
2. **What would kill it.** The re-run still shows the leak (K1 below).
3. **Tried before?** Nothing here refuted it. `--edit-model origin` already removes duplicates of a removed pair. README "Editing hard spots" measured 76 such duplicates among the 82 records origin adds to the list.
4. **Simplest version.** Once the run's removal list is final, read the event windows again, group the reads into duplicate families, and add the duplicates of every removed kept pair to the list. Pools, tiling, written reads and the random draws are untouched.
5. **What must be true of the inputs.**
   - The BAM is duplicate-marked. The hospital BAM is (Picard 3.3.0 `@PG`). An unmarked BAM has no flagged duplicate, so nothing is added.
   - A family key of the two mates' unclipped 5' ends with their strands (`origin`'s `five_prime` and `Fragment::family`) matches Picard's sets. Not verified; K1 and K2 measure it.
   - A removed pair's mates lie within its extraction window ± `MAX_FRAGMENT_LEN` (1,500 bp). Not verified; the log counts removed pairs whose two mates were not both seen.

**Design (locked).**
- **`origin::duplicates_of(records, removed)`.**
  - Group primary, non-QC-fail records by name. A fragment counts only when both its mates are there.
  - A family is removed when one of its fragments with no duplicate flag is in `removed`.
  - Return the names of the fragments in removed families whose mates both carry the duplicate flag.
- **`origin::removed_duplicates(bam, reference, spans, removed)`.** Read every span with origin's BAM/CRAM reader and keep only duplicate records and records named in `removed`. Return `duplicates_of` of them, plus how many names in `removed` had fewer than two mates in the spans.
- **`main`, clean mode only.**
  - `removed` = `fastq_removed_names(replaced_names, all_output_pairs)`, the originals removed and not written back.
  - The spans are each event's extraction windows (`extraction_bounds`; a fusion has two), widened by `MAX_FRAGMENT_LEN` on both sides.
  - The names returned are added to `replaced_names` before `replaced_reads.txt` and `fastq_removed_reads.txt` are written. So `merge.sh` and `fastq.sh` both drop them.
  - One log line gives the count, and the incomplete count.
- **`--edit-model origin` is not changed.**
- **Docs:** README (output rows, the paragraph on `replaced_reads.txt`, the origin paragraph that names this difference), the run README rows, and the `merge.sh` comment.

**Tests, written first and seen red:**
- a duplicate of a removed pair is returned;
- a duplicate of a kept pair is not;
- a duplicate one base off, or on the other strand, is not;
- a removed pair seen through one mate matches nothing, and is counted;
- `removed_duplicates` on a fixture BAM;
- the scan spans are the extraction windows widened by 1,500 bp (a near-0 window is clipped; a fusion gives two).

**Mutation checks (each must turn a test red):**
1. strand left out of the key;
2. any family returned, removed or not;
3. a fragment with one mate used;
4. the duplicate flag not required on the returned fragment;
5. no widening of the scan spans.

**Checks (locked before running).** Run `scripts/duplicates/run.sh` again with the new binary into a fresh directory. The same seeds give the same 43 events and the same stand-in.
- **K1: the leak is gone.**
  - **Pass:** D says "does not matter".
  - **Pass:** at most 3 source-duplicate reads remain among the old-allele reads at the 20 hom SNVs plus the reads wholly inside the 3 hom deletions. Before the fix: 23 + 8 = 31. The allowance of 3 (about 10%) is for families whose kept pair Picard and the key above disagree on.
  - **Pass:** C1, C2 and C3 pass again.
- **K2: only duplicates are added, and the right ones** (old run `scratchpad/dups/run/spike/out` against the new):
  - `R1.fq.gz` and `R2.fq.gz` decompress byte-identical. `truth.vcf` and `sim.bam`'s records (`samtools view`) are identical too.
  - `replaced_reads.txt` new = old plus added names, with at least 1 added, and `fastq_removed_reads.txt` new = old plus the same names.
  - An independent Python checker reads the source BAM over the same spans with its own key: the sorted (chrom, unclipped 5' end, strand) of both mates.
    - **Pass:** 100% of added names are duplicate pairs whose family's kept pair is in the old `fastq_removed_reads.txt`.
    - **Pass:** at least 99% of such duplicate pairs (both mates inside the spans) are added.
- **K3: origin is unchanged.** Old and new binaries, `--edit-model origin`, the same 43 events and seed. `replaced_reads.txt`, `fastq_removed_reads.txt`, `R1.fq.gz`, `R2.fq.gz` (decompressed) and `truth.vcf` must be byte-identical.

**Outcomes.**
- **Supported:** K1, K2 and K3 pass. Merge on the user's word.
- **K1 fails while K2 passes:** the leak has another source; look before anything else.
- **K2 or K3 fails:** the code is wrong; fix it with a test first.

**Reported:** how many duplicate pairs were added; the het alt share; the duplicate share near the sites; Picard's `PERCENT_DUPLICATION`; spike's wall time, old against new.

## Fix result (2026-10-04): supported

**Setup.**
- spike `5c4cd9a2`, from code commit `a93cc7f`.
- 669 tests. 9 of 9 planned-style mutants are caught, with the unmutated suite checked green first: 5 planned and 4 more (QC-fail filter, dedup across spans, the flag on both mates, a fusion's second side).
- Logs: `scratchpad/dups/{go-fix.log,go-fix.err,run-fix/,k2.log,k3.log}`.

**The run.**
- The same 43 events (`events.vcf` byte-identical to the first run). spike exited 0.
- The log line: `Duplicates of the removed originals: 638 pair(s), removed with them (...); 0 removed pair(s) had a mate outside the windows read`.
- `fastq.sh`: `removed 14558 original pairs, added 13495 of spike's pairs` (before: 13,920 removed).
- spike's own log spans 45 s, against 43 s before.

**K1: passes. The leak is gone.**
- D says **does not matter**. `R_spk` = 0/653 = 0.0000, against `R_real` 0.0017.
- **0** source-duplicate reads are left: 0 old-allele reads at the 20 hom SNVs, and 0 reads of any kind wholly inside the 3 hom deletions. Before the fix there were 23 + 8 = 31.
- C1 (100/100 against 0/100), C2 and C3 pass again.

**K2: passes. Only duplicates are added, and the right ones.**
- `R1.fq.gz`, `R2.fq.gz` and `truth.vcf` are identical, and so are `sim.bam`'s 201,212 records.
- `replaced_reads.txt` is 100,938 → 101,576 (638 added), and `fastq_removed_reads.txt` is 13,920 → 14,558, with the same 638 names.
- The independent Python key (sorted mates' unclipped 5' end and strand, read from the source BAM):
  - all 638 added names are duplicates of a kept pair in the old `fastq_removed_reads.txt`, so 100%;
  - 638 of the 638 such duplicate pairs inside the spans are added, so 1.0000.
- **Controls: both fail, as they must.**
  - The checker on old against old (nothing added) fails the 99% rule: 0 of 638.
  - The new lists plus 20 source duplicates spike did not list fail the 100% rule: 638 of 658.

**K3: passes. Origin is unchanged.** With `--edit-model origin`, the old and new binaries give byte-identical results:
- `replaced_reads.txt` (101,541 lines);
- `fastq_removed_reads.txt` (14,418);
- `truth.vcf`;
- `R1.fq.gz` and `R2.fq.gz` (100,605 records each, decompressed).

The new binary logs no clean-mode duplicate line in origin mode.

**Reported, not judged:**
- **Het SNVs:** spiked alt share 314/666 = 0.4715, against real 0.4964. Before the fix: 0.4631.
- **Duplicates within 300 bp of the spiked SNVs:** 55/6,702 = 0.0082 in the spiked BAM, against 0.0360 in the baseline.
  - spike still makes no duplicates of its own reads, so after marking, the spiked regions hold about a quarter of the duplicates real ones do.
  - Callers skip duplicates, so this shows only in duplicate counts and in a viewer that shows duplicates. It was not part of this fix.
- **Picard `PERCENT_DUPLICATION`:** spiked 0.040446, baseline 0.041689.
- **The BAM route at hom SNVs:** 7/660 old allele, as before. These are the cross-chromosome pairs spike cannot take (the `SIM_RESIST` share), not duplicates.

**Not run here:**
- raredisease, fastp and DeepVariant themselves;
- a real full FASTQ;
- a CRAM input (the reader is origin's, which its own tests cover on CRAM);
- an input BAM that was never duplicate-marked.
