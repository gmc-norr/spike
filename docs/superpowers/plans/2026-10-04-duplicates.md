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
