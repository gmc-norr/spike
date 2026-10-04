# origin: add reads back at the chance used to remove them

**Asked 2026-10-04.** This is review finding 2 (P1, `docs/review/2026-10-04-independent-review.md`). The user picked it ("finish the review") after fixes 1, 7, 3, 8, 5 and 4. `--edit-model origin` is experimental; the hospital demo uses `clean`.

**Reproduced on master `2ab9fc7`** (debug `38ade71d`). This used the reviewer's probe, from the patched scratch copy that survives the refusals added since. Its `discordant_confidence.bam` has 408 pairs. Every R1 is MAPQ 60. Every R2 is MAPQ 0 with one `XA` hit, at chr1:9001, outside the footprint. The event is `snp:chr1:5001:T:A;af=1`, run with `--min-mapq 0 --flank 2000 --seed 42 --threads 1`.
- `clean`: removed 149 pairs, added 148 (`Tiling 148 ... cov=12.0`).
- `origin`: removed 149 fragments, added **111**. The log says `origin depth at chr1:4999: 9.0x (the donor pool's there: 12.0x)`. 9.0 / 12.0 = 0.75.

The review also mentions a run "with the default MAPQ and flank settings, retaining high-quality donors outside the footprint". Its BAM is not in the script. On this BAM at default settings, both models refuse: 0 pairs pass `--min-mapq 20`. That run is not repeated here.

## Read in the code

- **Removal** (`src/origin.rs:231-248`) uses a fragment's `p_origin`, the chance that it came from the footprint. That is its surest mate's: the highest primary chance, with ties going to the higher in-footprint sum. The design (`docs/superpowers/specs/2026-09-26-edit-model-origin-design.md:101`) says why: "One mate pinned uniquely somewhere pins the fragment there." `decide` (`:566-589`) then draws once per duplicate family, against its highest member's total.
- **Origin depth** (`:392-410`, summed in `read_coverage_at`, `:437-468`) counts every unflagged read of a removable fragment. Each of its placements adds **that read's own** chance there (spec `:160-161`).
- **The tiling count** is computed from a coverage: the origin depth at the first covered breakpoint, times `f` to turn it into fragment depth (`src/simulate.rs:587-592`). `compute_tiling_count` takes that coverage, the VAF and the mean fragment length (`:970`). `SIM_DEPTH_FOLD` uses the same depth (`src/simulate.rs:286-291`).

So a pair with a MAPQ 60 mate and a MAPQ 0 mate with one hit outside is removed at about 1. It is counted at about 1 + 1/2 over its two reads, though, where `clean` counts 2. That is the 0.75 above. The spec defines the two rules in separate sections, and nothing checks that they agree.

## Gate A

1. **Principle.** A fragment came from one place, and both its reads came from there. So spike counts what it adds back at a place by the same chance it uses to take reads away from there.
2. **What would kill it.** On the reviewer's probe, origin still adds back clearly fewer pairs than it removes. Or the physics test (2026-09-27 plan, Task 12) leaves its bands.
3. **Refuted before?** No.
   - The surest-mate rule came in `a745304`. `git log -L` shows no later change to its lines.
   - The per-read depth came in `c0eeac1`. Its later commits, `41a5573` and `eda06a0`, made it faster with bit-equal output; they did not change the rule.
   - No case-file entry is about either.
4. **Simplest thing.** Change the depth side only. The removal rule, `decide`, and therefore which fragments a seed removes, stay as they are.
5. **Inputs, and how each was checked.**
   - Both mates of a removable fragment are read. This is R4 (`removable`, `:254-261`); a mate that was not read must be unmapped.
   - The probe's mates differ as described: counted from its SAM (408 R1 at MAPQ 60, 408 R2 at MAPQ 0, all 408 R2 with `XA`).
   - The hospital BAM carries `XA`. Counted on two regions: chr20:10-11 Mb has 247,023 primary records, 2,192 of them with `XA`; chr1:50-51 Mb has 249,853, with 2,757.
   - In the physics test every read over the footprint is MAPQ 0 with `XA` (its harness line on master: 2,066 of 2,066). So both mates of every pair there are equally unsure.

## Design (locked)

In `OriginSite::depth`:
- Each removable fragment gets a weight `W`: the highest `p_origin` among the removable fragments of its duplicate family. This is the chance `decide` draws the family's fate against, before the copy rate.
- For each unflagged read of that fragment, let `c` be the read's own in-footprint sum (`chance_within`). Each of its placements **inside** the footprint adds `chance x W / c` instead of `chance`. Its placements outside the footprint are unchanged. If `c` is 0, nothing is scaled.

So each unflagged read adds `W` inside the footprint, the fragment's chance of having come from there. When `c == W`, as for the surest mate, or for two equally sure mates, `W / c` is exactly 1.0. Those depths stay bit-equal.

Removal, `decide`, the copy rate, look-alike regions, R4 and the census rules do not change. The README's origin bullet on depth says what is counted.

**Tests, written first and seen red** (on the `origin.rs` test helpers). The reviewer's four cases:
1. **A sure mate and an unsure mate.** R1 MAPQ 60 at the spot; R2 MAPQ 0 at the spot, with one hit outside the footprint. Read depth over R2's spot is 1 - 1e-6, not 1/2.
2. **Both mates unsure.** Each read is MAPQ 0 with one hit outside, so each counts 1/2, as on master. This one is a guard: it is green on master by design.
3. **The unsure mate's hit on another contig.** As in 1, but R2's hit is on chr2. R2 counts 1 - 1e-6 at the spot, and its chr2 placement keeps 1/2.
4. **A duplicate family whose flagged member is surer.** The unflagged member's mates are both MAPQ 0 with one hit outside (1/2). The flagged member's R1 is MAPQ 60. The unflagged reads count 1 - 1e-6 each, the family's chance.

**Mutation checks.** Each must turn a test red. The runner first checks that the unmutated `origin::` tests are green.
1. No scaling (master's rule).
2. `W` is the fragment's own `p_origin`, not its family's highest.
3. `W` is the family's lowest.
4. Placements outside the footprint are scaled too.
5. `c` is summed over all of the read's placements, not only those inside the footprint.

## Checks (locked)

**K: the reviewer's probe.** The same patched script, against the new debug binary:
- `origin` removes 149 fragments (removal is unchanged);
- `origin` adds between 145 and 151 pairs, which is `clean`'s 148 ± 3. Master adds 111, which fails this.
- PREDICTED, not run: 148, since the origin depth becomes 12.0x.

**P: the physics test still passes.** `scripts/origin_physics.py` with the new release binary must give `VERDICT: SUPPORTED` under its locked rule. On master today it gives L 0.739 and P 0.763, against bands L [0.717, 0.820] and P [0.762, 0.865].

**B: unchanged where nothing should change.**
- **B1: `clean`.** The clean path never builds an origin depth. On the probe, `clean`'s `R1.fq.gz`, `R2.fq.gz`, `truth.vcf` and `fastq_removed_reads.txt` must be identical to master's. On the 100 kb slice below, its `replaced_reads.txt` must be identical too.
- **B2: the physics test's `origin` run.** Its `R1.fq.gz`, `R2.fq.gz`, `truth.vcf` and `replaced_reads.txt` must be identical to master's, because both mates of every pair there are equally unsure (Gate A, step 5). If they are not, the result names the fragments whose weight changed before anything else is concluded.

**D: the README's origin numbers, re-measured (reported).** Both were run on master today, and the README's figures reproduce except one:
- **The 100 kb slice.** `samtools view -b` of the 35x HG002 BAM over chr20:14,500,000-14,600,000, with `del:chr20:14530000-14531000` and default flags. `origin` lists 2,748 names and `clean` 2,736; 25 are `origin`'s alone and 13 are `clean`'s alone. These match the README.
- **RF14's site.** `del:chr20:7119236-7120236` on the 35x HG002 BAM, default flags, `--edit-model origin`. It logs `origin depth at chr20:7119235: 33.7x (the donor pool's there: 8.5x)`, `SIM_RESIST=0.016` and `SIM_DEPTH_FOLD=1.30`, all matching the README. But it removes **201** fragments at the spot and 0 at the look-alikes, where the README says 219.

The README takes the new binary's numbers, says which ones moved, and corrects the 219.

**R: how much it moves origin on real data (reported, not judged).**
- **Sites.** 150 random SNV sites on chr20 of the hospital BAM. They are drawn with seed 1, uniform on [100,000, length − 100,000), with REF from the reference (skip anything not ACGT). ALT is the next of A, C, G, T after REF.
- **Runs.** One run per site, with `--edit-model origin --allow-resistant --seed 1 --threads 16`, under master and the new release binary.
- **Report:**
  - the sites each binary refuses;
  - the new / master ratio of the logged origin depth and of the tiled count (min, median, max);
  - how many sites change at all;
  - the five largest changes, with their MAPQ 0 share in the footprint.

**T: speed.** `dup:chr20:10000000-11000000` on the 35x HG002 BAM, with `--edit-model origin --allow-resistant --seed 1 --threads 16`. Run master, new, master, new. The new binary's best time must be at most 1.25 x master's best. If it is slower, the speed is fixed before the result is written.

## Result (2026-10-04): supported

Code `6ef9c19`. Binaries: new release `8fc328b8`, new debug `8f913663`; master `2ab9fc7` release `920d2b7c`, debug `38ade71d`. Scratch: `scratchpad/ocp/` (logs, `mutate.py`, `lookprobe.py`).

**Tests.** 705 pass (4 new), 0 fail, 1 ignored. clippy gives 16 lines, the same as master.
- **Mutation.** The runner first checks that the 53 unmutated `origin::` tests are green. It then catches **5 of 5** mutants:
  - 1 (no scaling) reddens the sure/unsure, other-contig and duplicate-family tests;
  - 2 (own chance) and 3 (family lowest) redden the duplicate-family test;
  - 4 (outside scaled) reddens the other-contig test;
  - 5 (own chance over all placements) reddens 7 tests.
- **A changed test.** `simulate::tests::test_origin_depth_sums_what_summing_every_placement_sums` failed on the new code. Its 1,500 random reads used their random start as their 5' end, so reads starting at the same base formed duplicate families, and those now count their family's chance. The test is about the depth index, so each read now gets its own 5' end. Its brute-force side is unchanged.

**K: PASS.** `origin` removed 149 and added **148**, inside [145, 151]. Its log reads `origin depth at chr1:4999: 12.0x (the donor pool's there: 12.0x)`. Master added 111. Every other probe result is identical to master's, apart from timestamps and temp paths.

**P: PASS.** `VERDICT: SUPPORTED`, with L 0.739 and P 0.763, the same as master.

**B1: PASS.** Each probe run makes its own temp folder, so the two `truth.vcf` files differ in their `##reference` line only. B1 was therefore rerun on one set of input files, with both debug binaries. `clean`'s `R1.fq.gz`, `R2.fq.gz`, `truth.vcf`, `fastq_removed_reads.txt` and `replaced_reads.txt` are identical. On the slice, `clean`'s `replaced_reads.txt` is identical too.

**B2: PASS.** The physics runs also differ in `##reference` only. So `origin` was rerun on master's `donor.bam` and `ref.fa`, with the harness's arguments (`del:chrT:24500-25500 --seed 7 --allow-resistant`). `R1.fq.gz`, `R2.fq.gz`, `truth.vcf`, `replaced_reads.txt` and `fastq_removed_reads.txt` are identical. Master's rerun equals the harness's own run.

**D.**
- **The slice** is unchanged: 2,748 against 2,736 names, 25 and 13 alone. `origin`'s `replaced_reads.txt` is identical to master's.
- **RF14's site:** origin depth 33.7x, 142 tiled reads and 201 fragments removed, all as on master. `SIM_DEPTH_FOLD` moves from 1.30 to **1.29**.
- **The README** now says 1.29 and 201. The 219 it said before was already 201 on master.

**R (reported).** 150 sites.
- **Refused:** 5 sites, the same 5 under both binaries. 4 are too-few-reads refusals near the centromere (chr20:27.2-28.4 Mb). One, chr20:1598012, has no origin depth at any breakpoint.
- **Changed:** 9 of the other 145. Tiled reads, new / master: min 0.9746, median 1.0000, max 1.0089. Origin depth: 0.9744 to 1.0111.
- **Largest changes:**

| site | MAPQ 0 share | tiled reads | origin depth | pool depth |
|---|---|---|---|---|
| chr20:30264835 | 0.005 | 197 → 192 | 43.0 → 41.9 | 41.3 |
| chr20:29486139 | 0.083 | 294 → 296 | 63.3 → 63.8 | 52.8 |
| chr20:61987027 | 0.020 | 224 → 226 | 46.2 → 46.7 | 45.1 |
| chr20:26259627 | 0.012 | 224 → 225 | 48.7 → 48.7 | 46.5 |
| chr20:47827272 | 0.105 | 225 → 226 | 48.9 → 49.1 | 42.4 |

**Why the depth also goes down (measured, `lookprobe.py`, chr20:30264835).** At the spot, 32 reads are scaled up, by 0.62 in total, and 16 down, by 11.20. All 16 reads scaled down (listed by `lookprobe_all.py`) have a mate that is MAPQ 43-60 and placed at a look-alike, so its share in the footprint is at most `(1 - p_here)`, below 1e-4. Removal already treats such a pair as almost surely from elsewhere. Master's depth still counted the pair's other read (MAPQ 0-30) at the spot, at its own 0.17-1.00.

**T: PASS.** The 1 Mb DUP took 25.6 s (cold cache) and 11.79 s on master, and 11.76 s and 12.10 s on the new binary. Best against best is 0.998. Its outputs are identical to master's.

**What this means.** `origin` now adds back what it removes. On the reviewer's probe, the quarter it fell short is gone. On real data the change is small: 9 of 145 random hospital sites, each by at most 2.5%, and in both directions. `clean` does not change.
