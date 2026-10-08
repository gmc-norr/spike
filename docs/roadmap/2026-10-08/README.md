<!-- Written by the improvement-roadmap workflow (wf_db17a130-4b0, 2026-10-07/08): 6 map readers, 7 idea lenses, a merge, 2 skeptics per idea, a completeness critic and this synthesis. The workflow script is workflow/spike-improvement-roadmap.js; every agent's raw result, keyed by its label (map:*, ideas:*, merge, reality:*, value:*, critic), is in workflow/results.json. -->

# spike improvement roadmap, 2026-10-07

**Status.** Every idea now has both skeptic verdicts: M1–M49 from the lenses and C1–C7 from the critic, 56 ideas in all.
- Killed: 4 (M24, M32, C1, C6).
- Surviving: 52. Almost all of them were weakened, with a smaller scope, corrected evidence or a different falsifier.
- No idea reached value 5. The highest mean value is 3.5.

**Code state.**
- Local master is f18c7f5. Its code equals ba9f4e7; the commits after it are docs only.
- Public `spike/master` is 4a5672d. It still has the v2 decoy-CRAM abort and the old LDLR exon 1. The fixes are only in local ba9f4e7. See the open questions.

**Conventions**
- Evidence labels: **meas.** = measured (the file or script is named), **doc.** = documented (the doc is named), **inf.** = inferred.
- Value and cost are the mean of the two skeptic scores, on a 1–5 scale. Cost letters: ≤2.0 = S, 2.5–3.0 = M, 3.5–4.0 = L, ≥4.5 = XL. Where a skeptic re-scoped an idea, the reduced cost is also given.
- `scratch/` = the session scratchpad the evidence was made in (`/tmp/claude-1066/-home-parlar-ai-dev-spike/23a35c29-9600-4da2-a2c9-618425a1e0d3/scratchpad/`), which does not last. Its scripts, logs and small results (each at most 100 kB) were copied on 2026-10-08 to [`evidence/`](evidence/) beside this file, at the same paths: `scratch/X` and a bare `X/...` are `evidence/X`. Large files (BAMs, FASTQs, reference copies, spike run folders, anything over 100 kB) were not kept; [`evidence/SKIPPED.tsv`](evidence/SKIPPED.tsv) lists each with its size and reason. Baselines that lived only in those files have to be rebuilt with the scripts.
- Samples:
  - **35x** = HG002 NovaSeq 6000 PCR-free 35x BAM. The K2b held-out sets are in `scratch/qtq2/k2b`.
  - **hospital HG002** = D24-14230_Seq25-7600. It is GM24385/HG002 in a SeraCare inherited-cancer mix, even though the folder name says `na12878`.
  - **hospital HG001** = D25-7403_Seq25-7598, which is NA12878.
  - **DV** = DeepVariant. The hospital runs 1.6.1; only 1.9.0 is cached locally.
- Where a skeptic corrected a claim, only the corrected figure appears below.

---

## 1. Top 10

| # | Idea | Why (evidence kind) | Falsifying measurement (short) | Cost | Value |
|---|---|---|---|---|---|
| 1 | **M41** Background-integrity yardstick | Spiking breaks the sample's own indels in 214/752 events. Hom indels fall from median AF 0.934 to 0.479; a sham run changes 0/431 records (meas., `scratch/valid_lens/cr3_bg.py`, `m41_skeptic/sham_bg.py`). | Draw ≥50 random events (fixed seed). Run DV on the spiked view and on the sham. Kill if the after-minus-sham GT-change rate is <0.05 per event. | S | 3.5 |
| 2 | **M28** ClinVar-native VCF ingest | ClinVar records ≥50 bp have no SVTYPE. Fed by VCF, LDLR DUP 251140 plants 313 synthetic pairs with SIM_ALT_FRAGS 0; the coordinate spec plants 1,540 (meas., `scratch/cv/`). | Routing census of 3,996 P/LP records against locked counts. Byte identity with the del:/dup:/ins: specs. A shifted record and a GRCh37 record must be refused. | S | 3.5 |
| 3 | **M48 (+M46 simple form)** Per-variant demo receipt and SPIKE_ survival audit | No SNV, indel or exon DEL/DUP spike has ever been scored at caller level. DV on the unspiked hospital HG002 run matches GIAB het SNVs 1,845/1,850 (meas.). | Mutant shams: ALT swap ≤1% called, GT flip ≤2% matched, multi-allelic ≥95% after normalising. Per-event survival compared with SIM_ALT_FRAGS. Footprint diff against the unspiked run. | M | 3.5 |
| 4 | **M42** DeepVariant caller-equivalence (first step) | Real het call rates sit below the ceiling for INS20-49 (78/100) and DUP50-299 (24/54). At sites DV misses in the real reads, spiked reads carry more evidence (meas., `m42value/proxy.py`). | Paired McNemar per group and direction, with the shared-site discordance as margin. The old fe15a46 views must shift GQ (power control). | S first step, M full | 3.5 |
| 5 | **M10** In-read SNV/MNV editing | 92.3% of remade SNV pairs never touch the site. 36/100 SNV footprints damage a background indel (meas.). | DV on 100 footprints. Master must exceed the sham by ≥5 background GT changes, otherwise kill. The prototype must stay within sham + 1. | M | 3.5 |
| 6 | **M12** Male-X guard (reduced) | A DMD DEL at the default af on male HG002 exits 0 with SIM_VAF 0.500, GT 0/1 and no warning. The pipeline calls every male non-PAR X record 1/1 (meas.). | DV in haploid mode on 20 af=0.5 and 20 af=1 SNVs, compared with real hemizygous calls. Mutants: a PAR-blind version and an idxstats-based version must fail. | S reduced, L full | 3.5 |
| 7 | **M20** SMN1 deletion via `--edit-model origin`, scored by SMNCopyNumberCaller | Today's binary already removes 0.243 of SMN1+SMN2 reads (0.25 expected), and paralog-informative (PSV) reads only at SMN1. The caller has never been run on the output (meas., `scratch/m20/`). | Het: Total_CN_raw ≈2.88 ±0.2, SMN1 1, SMN2 2. Four-bin per-locus placement test. Must reject clean `--min-mapq 0 --allow-resistant`. | S | 3.5 |
| 8 | **M36** Provenance, simplest form | Pasting the README command runs seed 42 into /tmp/spike. A refused run leaves the previous run's truth and FASTQ in place (meas., `robust/stale`). | The pasted command gives byte-identical truth. A stale `--into-fastq` folder is cleared. `validate_pipeline.sh` step 3 still runs. | S | 3.0 |
| 9 | **M43** Manta equivalence, paired design | Real Alu-Alu het DELs ≥300 bp: Manta calls 34/77, some caller calls 62/77. The spiked ClinVar 4845403 is missed by Manta but found by TIDDIT and CNVnator (meas.). | Paired Δ across both directions with a bootstrap CI against a locked δ of 0.10–0.14. B1 must fail the PR+SR ratio. | M (reduced ≈2) | 3.5 |
| 10 | **M22** Streamed merge.sh | merge.sh took 1 h 42 min on B1. On a 387 MB slice the stream took 4.04–4.25 s against 49.8 s, with the same records (meas.). | Same record multiset, order and total on the 2.7 GB chr20 BAM and on B1. ≤20 min at 16 threads. A mutant without `^` must fail. | S | 3.0 |

**Why this order.** Items 1–5 serve the hospital LDLR/raredisease goal directly:
- background damage around every planted variant;
- ClinVar as the natural input;
- the first caller-level evidence;
- the cheapest fix for the most common class.

Items 6–7 are cheap guards or tests for pipeline steps that clinicians see: X-linked calls and SMN. Items 8–10 make results traceable and make the whole-BAM SV route usable, now that the slice route has been refuted for Manta and CNVnator.

---

## 2. Quick wins (cost S, value ≥3)

| Idea | Cost | Value | First step |
|---|---|---|---|
| M41 background integrity | S (1.5) | 3.5 | Add the read-AF column with its sham null to every transplant run (3 s per 750 events, meas.). Then the 50-event DV gate. |
| M28 ClinVar ingest | S (2.0) | 3.5 | Map CLNVC to SVTYPE for records ≥50 bp, keeping the REF and DUP-copy checks. Refuse Microsatellite, Indel, complex and N-containing records, and count them. |
| M20 SMN measurement | S (2.0) | 3.5 | No code: run today's origin binary for het and hom SMN1 DELs at ≥3 seeds, then SMNCopyNumberCaller and the four-bin placement count. |
| M22 streamed merge.sh | S (2.0) | 3.0 | Run the one-pass `-u -U -` form against today's script on the chr20 BAM with a real sim.bam. |
| M36 provenance | S (2.0) | 3.0 | build.rs commit stamp; shell-quoted command; record the effective seed and defaults; delete spike's own file names at startup; write truth via temp file and rename. |

**Just below the bar** (cost S, value 2–2.5). These are cheap enough to batch:
- C4: reword the RF8 `--min-mapq` hint.
- M31: refuse on a contig length (LN) mismatch.
- M26: delete the reference cache that is never hit.
- M29 leftovers: relabel the demo VCF, cap the "gene not found" list, add a MANE recipe.
- M34: drop pairs with a read holding >5 N.
- M2: emit Q2 as N.
- M16: log the chosen copy per event.
- M14: warn on random-base `INS:ME`.
- M11 step 1: a README sentence saying the sample's indels are dropped from remade reads.

---

## 3. Big bets

No surviving idea scored value 5, so no idea qualifies as a big bet by the rule. These are the largest bets with value ≥3 and cost L, each with the Gate B experiment that should come before any build.

- **M21** Coverage-preserving remap, scoped to DUPs (value 3.5, cost L as proposed, ≈2 when DUP-only).
  - **Gate B:** no build. Run master `--dup-model junction` and `--dup-model full` on the 1 Mb chr20 DUP and on the LDLR B3 DUP. Regress the added reads on donor depth per 1 kb bin, excluding bins within 1.5 kb of a breakpoint.
  - **Predictions:** uniform tiling gives slope ≈ −v; an image map gives slope ≈ +v.
  - **Rules:** master must fail the slope test (locked precondition). The remap passes if the 95% CI lies inside [0.4, 0.6] at v=0.5.
- **M11** Carry the sample's indels onto remade copies, CR3 option A (value 3.0, L).
  - **Gate B (V3):** hap.py on the spiked footprints against Q100 plus spike truth, beside the sham.
  - **Rule:** build only if background FP+FN beyond the sham is ≥0.1 per event (95% CI), measured after M10/M21 have been tried.
- **M46** Whole-sample acceptance for the hospital route (value 3.0, L as proposed, 2–3 in the simple form).
  - **Gate B:** a 2-minute check at the hospital: run fastp 0.23.4 with `--thread 6` twice on the head of the raw FASTQ and compare read order. If it reorders, the noise floor is an unspiked rerun, not a sham FASTQ.
- **M12** full ploidy and sex model (value 3.5, L).
  - **Gate B (M2 of its skeptic):** DV 1.6.1 with `--haploid_contigs chrX,chrY` on 20 af=0.5 and 20 af=1 SNVs in non-PAR chrX of hospital HG002, against real Q100 GT=1 sites.
  - **Rule:** change defaults or truth bytes only if af=0.5 gives RefCall, filtered or out-of-distribution GQ while af=1 matches the real calls.

---

## 4. Ideas by category

### 4.1 Read realism

#### M2 Substitution spectrum, Q2→N, per-mate rates (value 2.5, cost S, kill 0)
- **Problem.** The error table has no mate dimension (quality.rs:253-264). Wrong bases are uniform among the other three (synth.rs:719-728). The called base is never stored (quality.rs:1014). N is left out of learning (quality.rs:1006), so a drawn Q2 falls back to the nominal 0.631 error rate and emits a called base.
- **Evidence** (meas., `scratch/phys/errproc_35x.out`, `errproc_hosp.out`):
  - Q<15 top-alt share, R1/R2: 35x real 0.747/0.692; hospital 0.837/0.819; spike21 0.337/0.338.
  - Hospital at Q24: 0.864/0.833. At Q40: 0.625/0.592.
  - Q2: 35x real 495 bases, all N; hospital 782, all N; spike21 441, none N.
  - R1/R2 mismatch ratio: 35x real 0.433 (Q25) and 0.542 (Q37); hospital 0.680 (Q24) and 0.483 (Q40); master 0.873/0.795.
- **Corrections.**
  - Spike does write N at Q2 for reference N, contig padding and deletion padding (synth.rs:189-203, 245-248). Only called bases never become N.
  - Mate pooling explains the R1 excess only at Q25. At Q37 both mates over-err: pooled real is 0.0119%, against spike R1 0.0186% and R2 0.0234%. In the last 10 cycles of class ≥35 reads, spike is 3.77e-4 against real 9.01e-5 (meas., hiq_35x.out). A mate split moves this excess between the mates but does not remove it.
  - For callers, the relevant hospital bands are Q24 and Q40. Q9 is below DV's minimum base quality (tool docs, not checked here).
  - The spectrum alone would take ≥3-read sites from 2 to about 8 and ≥2-read sites from 329 to about 600, against real 181 and 1,247 (inf.).
  - The fastp N-filter impact does not hold.
  - The pattern is per instrument or run, not per library.
- **Proposal.**
  1. Q2→N: about 5 lines.
  2. A 4×4 emission matrix per Q bin, pooled over mates, learned from the called base. Use one matrix per quality symbol, or at least a Q<30/Q≥30 split; a split at Q15 merges different spectra.
  3. A per-(mate, symbol) multiplier on the cell rate, instead of a mate bit that halves the 10,240 cells.
  4. Widen the variant mask only around CIGAR I/D seen in the pool.
- **Simplest.** Q2→N alone, gated on the N check plus K7.
- **Falsifier.**
  - Before building the matrix: on the odd-numbered blocks, redraw spike21's alt counts from a matrix learned on the even blocks. Build only if the redraw closes ≥25% of the ≥2-read gap (329 → 1,247).
  - After building: on held-out windows, not the learning blocks, the Q<15 top-alt share per true base is within ±0.05 of real. The R1/R2 ratio is within ±15% of the targets above. A transposed matrix, swapped mates and master must all fail.
  - N: zero non-N Q2 bases, and an N rate within 2× of 2.6e-5.
  - Guards: K1, K2 (passed at z −2.64) and K7.
- **Gate B.** Learnability on even vs odd blocks of real reads.
- **Cost.** About one day with K7 reruns, mutants and guards (cost 2).
- **Risks.** Thinner cells could push K2 past its limit. The output bytes change. A spectrum at Q≥30 may contain unmasked variants.
- **Skeptic notes.** The gaps are real and cheap to close. The effect on calls is small. Drop the clip-boundary mask: it contradicts bad-tails-20261007, where shared clip spots were mostly sequence-triggered bad ends (27% vs 3% shuffled).

#### M1 Site-specific strand errors (value 2.0, cost L, kill 0)
- **Problem.** The error draw has no genomic position (synth.rs:208; quality.rs:894-921). block_errors masks positions with ≥5 reads and ≥10% differing (quality.rs:963). Real recurrent same-strand errors disappear inside remade footprints.
- **Evidence** (meas., 35x K2b, `phys/recur.py` → `recur_35x.out`, reproduced):
  - Sites with ≥3 same-strand same-alt reads: real 181, spike 2 and 2.
  - Candidates (≥2 alt reads, AF ≥0.12): real 453 (500/Mb, 366 one-strand); spike 85 and 74.
  - Split-half (meas., `phys/split_both.out`): P(B≥2 | A≥3) = 0.372 (n=191) against a base rate of 0.000228. The other strand repeats 3.7%.
  - Runs: 48/181 sites have a run ≥7 nearby; 111/181 have no run ≥5.
- **Corrections.**
  - At BQ≥20 the 35x numbers are real 33 candidates and 13 sites, against spike 1 and 0 (meas., `sk/recur_other.out`).
  - Hospital at BQ≥10: 36 sites with ≥3 reads, 137 with ≥2, and 49 candidates (51.8/Mb). 85–93% of hospital recurrence rides on Q9 bases (meas., `skeptic_sites/`).
  - The "~5× more candidates" figure was measured on 35x, not on hospital data.
  - Redrawing spike's error bases from the real substitution matrix alone takes candidates from 85 to 154–158 (meas., `matrix_inject.out`).
  - Matrix-independent metric: same-strand any-alt ≥3 is real 231 vs spike 16.
  - Hospital DV over the 20 blocks: 1,351 true SNVs and exactly one isolated false positive (chr19:19835502, 3 reverse-strand Q40 alt reads). That is about 1 site-error call per Mb, so a footprint loses about 0.004 calls per event (inf.).
  - The dupcheck figure "690/757" has no output file and is unverified.
- **Proposal.** A per-site, per-strand excess error map from the pool, shrunk toward the iid rate, emitted on reference-origin bases.
- **Simplest.** Pool sites with ≥2 same-strand same-alt reads and 0 on the other strand. Rate = (k − expected)/n, capped at 0.5. Flank reads only.
- **Falsifier.**
  - Value test first: DV 1.6.1 on real vs master-remade reads over ≥20–50 Mb of hospital bench. Kill for the raredisease goal if master's non-truth call rate is not significantly below real.
  - Otherwise, a hospital K2b-style set at BQ≥10 and BQ≥20: within 0.7–1.3× of real. Site overlap must reach ≥0.6× the split-half oracle, not a fixed 40%, which sits above the oracle.
  - A matrix-only arm is the baseline. Shuffled positions and swapped strands must fail. Emitted errors must keep their low-Q coupling.
- **Gate B.** Offline injection into spike21 with the existing `phys/` scripts, at real pool depth.
- **Cost.** L. A faithful version must also bias the quality draw at the site.
- **Risks.** Thin pools (about 15–17 reads per strand); copying unlisted variants.
- **Skeptic notes.** The gap is real and new. For germline DV at the hospital the value is near zero. It matters for low-VAF or somatic use and for the IGV look. Do M2's matrix first.

#### M5 Last-cycle quality bin and physical cycle axis (value 2.0, cost S, kill 0)
- **Problem.** Position bins are 16 cycles wide (quality.rs:547, 718) and end_bin lumps the last 10 (162-169). Learning takes cycles_left from qual.len() (737, 1023), while generation runs on physical cycles and then trims (synth.rs:147). So on trimmed libraries learning and generation already disagree.
- **Evidence** (meas., errproc, `skeptic_lastcycle/fqend.py`, `skeptic_cyc/`):
  - 35x R1 cycle 151 mismatch 0.538% vs 0.324% at 146–150 (1.66×); R2 1.57×. Spike is flat (0.392 vs 0.414).
  - Cycle-151 low share, R1: real 5.76% vs 4.17%; spike 3.95% vs 4.01%. R2: 7.46% vs 5.66%.
  - Hospital (full-length reads): R1 2.27% at cycle 151 vs 1.58% at 150; R2 2.72% vs 1.83%. The last base of 150-bp reads has no dip (1.53%/1.59%).
  - Hospital read lengths: 61.0% 151 bp, 36% 147–150 bp (adapter-prefix trimmed), 2.6% cut to their fragment.
- **Corrections.**
  - The 4.2× excess of aligned high-Q errors in the last 10 cycles is a clumping gap, not a per-cycle gap. Over the last 4 cycles, placed errors are 160 vs real 169. But real has 24 reads with ≥2 errors there, against 0.09 expected; real clips 251/334 of its last-10 errors, spike 68/421. Bins would leave spike at about 3.4× real (inf.).
  - The 3' short-clip gap disappears once clips within 10 bp of a Q100 record are masked: 0.131% vs 0.122%.
  - The "272 high-Q clip errors" count covers both ends and is mostly the sample's variants.
  - 5' qualities already match (cycle 1–5 mean Q 36.21 vs 36.20).
  - A qual.len()-based last-cycle bin would be about 39% diluted on hospital data.
- **Proposal.** A last-cycle position in the quality key, plus physical cycles_left = read_length − c in count_read, block_errors and generation.
- **Simplest.** As above, about 20 lines.
- **Falsifier.**
  - The excess low share at physical cycle 151 over 150, per mate, must be within ±25% of real on both BAMs, over 3 seeds. Targets: 35x R1 1.59 and R2 1.78 points; hospital R1 0.69 and R2 0.89.
  - The last base of 150-bp hospital reads stays within 0.15 points of physical cycle 150.
  - Master gives about 0. A qual.len() version gives about 58% of the hospital excess. Both must fail.
  - The K2b harness must keep trimmed reads (synth.rs:2752, 2759 drop them today).
  - Report the placed error rate per cycle, but do not gate on the aligned 4.2×.
- **Gate B.** The share of cycle-1 bases whose full error cell reaches 200 bases. The backoff drops the end dimension (quality.rs:900-904).
- **Cost.** S (cost 2 with the K1/K2/K2b/K7 reruns).
- **Risks.** Thin terminal cells.
- **Skeptic notes.** No caller outcome depends on it; the value is detectability. Clumping deserves its own idea: reads with ≥2 high-Q errors in the last 4 cycles, real 3.3e-4 vs spike 1.4e-4.

#### M3 Homopolymer slips per library (value 2.0, cost M, kill 0)
- **Problem.** The default `--indel-error-rate` is 0 (main.rs:178-182). When enabled, the indel branch inserts a random base anywhere downstream (synth.rs:210-226), not ±1 at the run.
- **Evidence** (meas.; slips per pass, with the bench BED or the chr17 block removed):
  - 35x: 8–9 bp runs 0.23% (12/5,113); 10–11 0.72% (10/1,394); 12–14 0.53–0.72%; 15+ 0.95% (6/632).
  - Hospital: 0.32% (28/8,821); 1.84–1.85%; 2.44–2.51%; 5.26–5.35%.
  - spike21: 0 at every run ≥8.
  - Library ratio on the same runs: 2.7×, 4.6× and 5.6×. "Library" here mixes instrument, trimming and batch.
- **Corrections.**
  - The unrestricted 35x 15+ figure (2.2%) is inflated by a missed sample variant at chr17:36270887 (A×18, 6/12 passes) in a block with no bench-BED interval.
  - Slips belong to the fragment: overlapping pairs give 0 one-mate slips against about 12.6 expected (hospital) and about 3.2 (35x). A per-read model has the wrong shape. Re-measure on a whole hospital chromosome first.
  - A fix-B prototype already exists (`qm2gate/hp_apply6.py`). It was judged only on clip share (+0.04 points).
  - Impact at coding runs is small: 6–9 bp runs slip 0.02–0.3% per pass. PanelApp green exons hold 29 runs ≥10 bp (2.9/Mb), against 372/Mb in the sample blocks.
  - The "rate ratio collapses at `--indel-error-rate` 0.05" claim is inferred.
- **Proposal.** A per-fragment ±1 template edit at runs, drawn once per fragment and read by both mates. Mask a run when one signed change is ≥30% of its passes. Fit log p = a + bL.
- **Simplest.** Per-fragment ±1 at runs ≥10 with a two-parameter p(L).
- **Falsifier.**
  - Step 0, value, on real reads: remove stutter fragments from Q100 1-bp indels in runs ≥10. Value exists only if DV genotype, FILTER or GQ bin changes at ≥5% of sites beyond random-removal controls. Add a low-AF arm.
  - Build: hospital rates for 10–14 and 15+ inside the real Poisson 95% CI. Ins share within ±0.15. One-mate pairs ≤20%. A per-read mutant, master and `--indel-error-rate 0.05` must fail.
- **Gate B.** Hash-split halves, ≥30 events per long bin, halves agreeing within 1.5×.
- **Cost.** M. Generation must fetch one template for the whole fragment.
- **Risks.** CR3 (M11) erases the sample's own indels at these runs anyway, so fix that first.
- **Skeptic notes.** Frame this as low-AF homopolymer realism, not germline.

#### M6 Dead-tail read state (value 2.0, cost L, kill 0)
- **Problem.** Spike's aligned high-Q mismatch rate is 1.76× real.
- **Evidence** (meas.):
  - Rate: 1.188e-4 vs 2.097e-4.
  - Reads with ≥1 aligned high-Q mismatch: 1.42% vs 2.21%.
  - K2b crash share 0.90–0.91% vs 1.06% (doc.).
- **Corrections** (meas., `deadtail_check/where.py`, `where_hiq.py`, `crash_by_var.py`).
  - Spike counts more high-Q errors than real (4,839 vs 3,811, 1.27×). It does not "learn the right number".
  - Only 16% of real high-Q errors sit in crashed tails.
  - About 76% of the excess is spike erring again within 30 cycles of an earlier error in reads that did not crash (5.3× and 3.6× real). The backward error memory is too sticky.
  - The crash shortfall sits at reads over Q100 indels: 5.73% vs 3.61%/3.37% (z −5.8/−6.6). Reads with no Q100 record match (0.83% vs 0.78/0.82%). That points to CR3, not a missing state.
  - The read class already is a per-read state; only the onset cycle is missing.
  - The hospital excess has never been shown.
- **Proposal.** A latent dead/alive state with an onset cycle.
- **Simplest.** A two-state flag with an empirical onset per class.
- **Falsifier.** Cheap ablations first:
  1. Keep class and end bin in the error-history fallback (quality.rs:919 drops both).
  2. Mask clusters near Q100 variants.
  3. Step 0: draw errors along real quality strings.
  - Build the state only if these leave the aligned not-crashed high-Q errors outside ±20% of real (1,450 on variant-free templates) and the hospital shows the excess.
  - A Q≥30×0.6 rescale must fail the clip-share criterion.
- **Gate B.** Step 0 decomposition (1–2 h, existing qm2gate tooling).
- **Cost.** L+.
- **Risks.** Every previous latent-state attempt failed as a quality model.
- **Skeptic notes.** Kill condition: drop this if the cheaper fixes close the excess. Value is low; the extra errors are random, not recurrent.

#### M4 Molecule errors shared by both mates (value 2.0, cost M, kill 0)
- **Evidence** (meas., `phys/mateconc2.py`, `molerr/strat.py`, `skeptic_mol/`):
  - Both mates Q≥30, same wrong base: 35x 44 sites in 913,655 bases (4.8e-5); hospital 141 in 1,454,067 (9.7e-5); spike21 0.
  - Shared share of erroneous high-Q positions: 49% (35x), 77% (hospital), 0% (spike).
  - 43/44 high-Q 35x sites are singletons.
- **Corrections.**
  - It is not "decisive for low VAF". At AF 0.05 about 95% of footprint reads stay real (inf. from the README low-AF table), so they keep their real molecule errors.
  - At AF 0.5–1 the effect is about 0.5–1 two-read positions per footprint at AF ~0.07, below DV's 0.12 (inf.).
  - The low-Q half of the 35x concordance (50/104 sites) fits the biased spectrum, which is M2.
  - The spectrum order is noise past T>C.
  - Spike's single-mate high-Q errors in overlaps are already 229+241 vs real 23+23, so the proposed "subtract from the table" step is ill-posed.
  - Qualities must be conditioned on the error.
  - REVIEW "M7" is seed determinism, not duplicates; spike makes no synthetic duplicates.
- **Proposal.** Mutate the fragment once, before cutting both mates.
- **Simplest.** One rate from the both-Q≥30 stratum plus a 12-cell table.
- **Falsifier.** Value test first: inject shared errors post hoc into the spike21 FASTQs, realign, run DV/HC plus Mutect2 or LoFreq. Build only if injection closes at least half of a significant false-positive gap. If built, score per quality stratum, each with its own denominator.
- **Gate B.** As above (about 2 h).
- **Cost.** M.
- **Risks.** Fix the high-Q 3'-end excess first.
- **Skeptic notes.** Low priority.

#### M9 Adapter read-through for untrimmed libraries (value 2.0, cost M, kill 0)
- **Evidence.**
  - 35x chr20:10–11 Mb: 0.831% of proper R1 reads have |TLEN|<151; 0.473% have <140 (meas.).
  - The adapter consensus is exactly TruSeq for 33 bp, minimum per-position agreement 0.895/0.932 (meas., `adapter_check/cons.py`).
  - Adapter clips are 0.71% of reads on another region (doc., soft-clips README; an earlier pass gave 0.48%).
- **Corrections.**
  - The user chose to leave adapter clips out (quality-model plan :4), and the read-length plan deferred them. Asking first is required.
  - The hospital library is sequencer-trimmed, and raredisease fastp trims real adapters anyway.
  - The fragment-length floor is already wired (types.rs:211, main.rs:1930).
  - The real costs are re-baselining every 35x byte-identity case, and a misfire on software-trimmed libraries: bam_stats would call them untrimmed and a hard-coded TruSeq would add adapter.
  - The "|TLEN|<151 ±0.2" gate passes by construction.
- **Proposal.** Learn each mate's tail only from short pairs whose reads are longer than TLEN. With no read-through seen, fall back to trimmed-mode truncation.
- **Simplest.** Park it with a README limitation note.
- **Falsifier.** Step 0 counterfactual on real 35x chr20: remove or cut the short pairs, then rerun Manta and DV/HC. Park if no Manta call ≥50 bp changes and F1 moves <0.001.
- **Gate B.** As above.
- **Cost.** M.
- **Risks.** 35x RNG baselines move.

#### M7 Synthetic duplicate families (value 1.5, cost M, kill 1)
- **Evidence.**
  - Hospital duplicate rate 4.71%, of which 59.6% optical (meas., MarkDuplicates metrics).
  - Duplicate share within 300 bp of 40 chr20 stand-in SNVs: 0.0082 vs 0.0360 (meas., `dups/run-fix`).
  - BAM route, ±3 kb of the 22 K2 deletions: 0.0276 vs 0.0370 (meas.).
- **Corrections.**
  - raredisease has no per-region duplicate QC.
  - Whole-library PERCENT_DUPLICATION moved only 0.041689 → 0.040446 even on the dense stand-in.
  - FastQC never sees SPIKE_ reads; they are appended at the end of the file.
  - 93% of duplicate sets are pairs.
  - The pool holds no duplicates (extract.rs:876/934).
  - Nothing in germline raredisease uses family-shared errors.
  - BAM-route copies would be unflagged, because merge.sh never marks duplicates (origin's R3 defect again).
  - The optical fraction varies by window: 0.703 in the chr20 window; 0.43–0.59 in other hospital runs.
  - The naming question was already settled by read-names K1b (50c75b6).
- **Proposal.** FASTQ-route only, with one copy probability per pair.
- **Simplest.** Same; decisions taken from a name hash.
- **Falsifier.** Step 0 value test: awk-copy about 4.7% of SPIKE_ pairs, then run DV, Manta, CNVnator and TIDDIT. The idea stays closed if nothing changes. The same copies with 0x400 stripped serve as positive control.
- **Gate B.** Not needed (naming settled).
- **Cost.** M.
- **Risks.** RNG stream; byte-identity gates.
- **Skeptic notes.** Deferred on purpose twice. One kill vote on value.

#### M8 Fragment-end composition motif (value 1.5, cost M, kill 1)
- **Evidence** (meas., `phys/startcomp2.py`, re-run).
  - Hospital K2 cycle 1: sample reads C .337, G .291; SPIKE_ C .226, G .222, flat.
  - The bias runs to cycle 10. It holds genome-wide in both mates and in both hospital libraries. 35x is mild.
- **Corrections.**
  - The chemistry is unknown. The adapter is TruSeq-type; only 1/6,137 full-length reads end in A.
  - There is a third uniform start draw at simulate.rs:937.
  - A 3-bp PWM fails its own cycle 4–6 test.
  - The "per-read tell" contradicts the owner's distinguishable-SPIKE_ rule. The classifier AUC is only 0.591–0.624.
  - FastQC cannot see the SPIKE_ block, and the module already fails on the real sample.
  - "GC bias 0.88–1.08" is unverified.
  - Gate B passes by construction.
- **Proposal / simplest.** Only behind an opt-in flag.
- **Falsifier.** Value gate first: a ReadPosRankSum/GQ comparison against real look-alikes, or V1 FastQC/fastp status. Kill if every status matches.
- **Gate B.** Already answered.
- **Cost.** M.
- **Risks.** Changes the RNG stream.
- **Skeptic notes.** Record as "measured, real, left alone". M10/M21 would bring real starts for most reads for free.

#### M10 In-read base editing for SNVs/MNVs (P1b) (value 3.5, cost M, kill 0)
- **Problem.** Clean mode remakes every event-copy pair over ±2 kb (simulate.rs:250-263, 817). Real qualities, duplicates and the sample's indels are lost there.
- **Evidence.**
  - 92.3% of 16,089 remade pairs over 80 SNV events never reach the site; per event 3.6–11.9% do (meas., transplant2 run0).
  - Under SNV events: 35 hom background indels fall from median 0.955 to 0.462; 38 het from 0.464 to 0.271.
  - 36/100 SNV footprints damage ≥1 background indel. 178 of 214 damaged footprints are under indel or DUP events, which P1b does not fix (meas., `p1b_value/bg_by_group.py`).
  - 9/18 LDLR exon footprints hold a hom Q100 non-SNP (meas.).
- **Corrections.**
  - fastp corrects a mismatch only when one base is ≥Q30 and the other ≤Q14 (basecorrector.cpp). Hospital bases are about 94–96% Q40, and both mates cover the site in only 5.6% of covering pairs (meas.). So the proposed "one-mate edit loses ALT after fastp" control would not fail.
  - The "0.17% ref reads at real hom SNVs" figure is not a measured gap: P(0 of 653) = 0.33.
  - KS is invalid on 4 tied quality values.
  - Duncavage et al. is J Mol Diagn 2023;25(1):3-16, and that category also describes spike's current method.
  - merge.sh does no duplicate re-marking.
  - pysam is available in system python3 (0.23.3), but not under `python3 -I`.
  - SNV site evidence is already supported, so the gain is in the flanks.
- **Proposal.** Rewrite the base in covering event-copy fragments at copy_rate, in both mates and in the duplicate family. Keep qualities and errors. Rename to SPIKE_.
- **Simplest.** For both routes: add an edited copy of each covering pair, built from the BAM primary record put back in sequencing orientation, with the original quality and tile:x:y, and list the original as removed. Caveat: BAM sequence is post-fastp.
- **Falsifier (Gate B′)**, on round 2b's 100 forward SNV events; the master arm is already on disk.
  - DV 1.6.1 on footprint slices for unspiked, sham, master and prototype. Slice calls must first equal whole-BAM calls.
  - Master must exceed the sham by ≥5 background GT-changed records, otherwise KILL. The prototype must stay within sham + 1.
  - Mate concordance at the site against real carriers (Fisher); a one-mate mutant must fail.
  - Edited duplicate families: 0 REF representatives after Picard.
  - Round 3 SNV A metric, both directions.
  - Do not use "no flank metric moves" as a kill rule.
- **Gate B.** As above (hours).
- **Cost.** M or more: MNVs, clips, overlapping events, LOH, origin mode.
- **Risks.** Conflicts with the convention "spike never edits a real read in place" (README:1427); this needs the owner. RNG changes.
- **Skeptic notes.** Weigh it against the cheaper alternative of re-synthesising only site-covering pairs. Build it before M11, M15 and M16 groups.

### 4.2 Biological realism

#### M11 Carry the sample's indels onto remade copies (CR3 option A) (value 3.0, cost L, kill 0)
- **Problem.** SampleCopies is SNP-only (loh.rs:34-43). The parser drops non-SNPs (554-557). Truth carries no phase (truth.rs:232, 468-474).
- **Evidence.**
  - R11: 17 hom indels fall from median 0.95 to 0.48; 0/211 synthetic reads carry them (meas., `cr3probe/out80.tsv`). 16/17 are 1/1 in hospital DV.
  - Het: 41 distinct positions, 0 synthetic carriers.
  - LDLR: 18/18 exon footprints hold a non-SNP, 9/18 a hom one (meas., MANE BED).
  - bcftools 1.9 mini Gate B: 0/9 recipient 1/1 sites stay 1/1 (meas., `cr3probe/gateb`). bcftools mis-genotypes homopolymer indels, though.
- **Corrections.**
  - `--gvcf` already reads phased GT|PS (47f0f64). Extend it rather than adding a flag.
  - T7's result was "supported: a warning would fire on 38/40", not refuted.
  - Q100 and Platinum carry no PS, so phase is contig-wide.
  - Clinical visibility at LDLR is about nil. raredisease drops gnomAD AF >0.70 before ranking: all 8 LDLR Q100 hom indels are DV 1/1 but absent from the ranked VCFs. 9 of the 17 R11 indels reach research-ranked with RankScore −4 to −10 (meas.).
  - R11 skipped indels <300 bp from the event, so the effect on the planted call is unmeasured. 3/18 LDLR exons have a non-SNP within 150 bp.
  - Value gate (a) as written cannot fail.
- **Proposal and order.**
  1. Now: a README note that the sample's indels are dropped.
  2. M10 first.
  3. If V3 still shows residual error: carry hom-alt indels onto every remade copy, with no phasing.
  4. Then phased het indels and phased truth, designed together with CR7.
- **Simplest.** SNPs plus indels ≤50 bp, small-variant events only, inside the truth BED, byte-identical without the flag.
- **Falsifier.**
  - V1: hospital outcome through the gnomAD filter and the clinical list.
  - V2: planted-call outcome, variants with a background indel within 150 bp against matched ones without.
  - V3: hap.py footprint FP+FN per event beyond the sham. Build only if ≥0.1 per event remains after M10/M21.
  - Mechanism gate: exact allele counting after `bcftools norm`; a paired test against the recipient's real carry rate; a cis/trans phase check; a fixed event copy.
- **Gate B.** V3 (DV 1.9.0 is cached).
- **Cost.** L. Background segments must stay out of junctions(), SIM_ALT_FRAGS, own_bases and the planted probes.
- **Risks.** Mask the 7 SeraCare genes ±2 kb on D24-14230. Truth errors become spike errors.
- **Skeptic notes.** Value lies with benchmarkers and realism purity more than with the demo VCF.

#### M12 Sex and ploidy awareness on chrX/chrY (value 3.5, cost L full / S reduced, kill 0)
- **Evidence** (meas.).
  - chrX/autosome density: hospital 0.499, 35x 0.501, strobealign 0.506, HG001 0.977. PAR1 depth ratio 0.910.
  - Pipeline: male non-PAR X 105,196/105,196 calls 1/1; PAR1 3,794 and PAR2 279 also 1/1.
  - 0.45% of male non-PAR SNVs have ref fraction 0.2–0.8.
  - ClinVar P/LP: 7.20% on non-PAR X.
  - del:chrX:31500000-31502000 at af=het: exit 0, SIM_VAF 0.500, GT 0/1. At af=hom: SIM_VAF 1.000, GT 1/1.
- **Corrections.**
  - A chr20-only slice gives a ratio of 0.007 and would be called male, so use windowed index counts with an "unknown" guard.
  - Defaulting to GT '1' contradicts the pipeline. Keep '1/1' plus SIM_PLOIDY=1; make '1' opt-in.
  - Haploid PAR is a pipeline setting (`--haploid_contigs chrX,chrY`, `par_bed: null`), not biology. Spike's PAR stays diploid.
  - Skipping phasing is moot at af=1.
  - The full bundle is L.
  - The R12 hemizygous AF numbers are not re-verified.
- **Proposal (reduced).**
  - A warning on non-PAR X/Y when the background looks male, from 5 windows of 1 Mb at MAPQ≥20, with a built-in GRCh38 PAR table and byte-identical output.
  - A README line telling users to use af=hom there.
  - Optional `--sex male`: omitted af defaults to 1.0, plus SIM_PLOIDY.
  - Drop `--strict-ploidy`, `--ploidy-bed`, sex inference as a decision-maker, and CR7(a) riding along.
- **Simplest.** The warning plus the README line.
- **Falsifier.**
  - M1 (minutes): warns on HG002 35x and hospital HG002 at the default af. Silent on af=hom, on the HG001 CRAM and on PAR1. "Unknown" on the chr20 slice. PAR-blind and idxstats mutants fail.
  - M2: the DV haploid-mode test described in §3.
  - M3: a DMD DEL through Manta and gCNV against real Q100 chrX DELs ≥50 bp.
- **Gate B.** M1 and M2.
- **Cost.** S reduced, L full.
- **Risks.** The af default change on male X needs a release note.
- **Skeptic notes.** Tell the hospital that `par_bed` is null.

#### M13 Homology-aware breakpoints (value 2.5, cost M, kill 0)
- **Evidence** (meas., `m13skeptic/probe.py`).
  - Demo LDLR VCF: 15/22 records have both ends on multiples of 100. Demo DELs have median homology 0 bp; 12/19 are blunt.
  - 19 sequence-resolved ClinVar LDLR DELs: median 2 bp; 11/19 ≥2 bp; 3/19 ≥20 bp.
  - Rounding 15 precise Alu-Alu DELs drops their 3–40 bp homology to 0–2 bp.
- **Corrections.**
  - "Alu pair within 500 bp of both ends" holds by chance for 81% of random events.
  - Real LDLR DEL ends are not Alu-enriched: 6/19 have both ends in an Alu, against 26–29% random.
  - The low-MAPQ claim is wrong. Real Alu-Alu carriers' evidence reads have MAPQ median 60. What separates them is few clipped reads: median 2 vs 15 for no-repeat events.
  - 4845403's homology is disputed: 0 bp under the VCF-anchor convention vs 33 bp after a 1-bp shift. M49's read check found about 26 bp of near-identical AluSx at the junction.
  - No rounded LDLR spike exists, because B1–B3 used ClinVar coordinates.
  - The SIMPLEST rule snaps 22/22 demo records to 10–37 bp homology, moves ends up to 827 bp, recovers ≤1/19 real junctions and changes the exons in 2/22.
- **Proposal.** Opt-in snapping only, with an empirical junction-class mixture. Use an rmsk-free longest exact match ≥ about 20 bp (chance is about 11–12 bp, inf.). Record CIPOS/CIEND and SIM_SNAP.
- **Simplest.** No code. Replace or relabel the demo VCF (M29 part 2), document that exon specs cut at exon edges, and use precise ClinVar coordinates.
- **Falsifier.**
  - Gate 2 first: Manta and TIDDIT on B1 precise vs B1 rounded, and on two small demo exon DELs. Build only if the junction class flips a call.
  - Gate 0: a snap must restore the true junction ±300 bp in ≥10/15 rounded real events, and must not invent Alu-Alu junctions in >3/6 no-repeat events.
- **Gate B.** As above.
- **Cost.** M (S for DEL-only).
- **Risks.** Snapping imposes NAHR rather than discovering it.

#### M14 Mobile-element insertion content (value 2.0, cost M, kill 0)
- **Evidence.**
  - SeraCare AluY: 26 clipped reads. 15 of 26 have mates off-contig, 14/15 of those at MAPQ 0 (meas., `mei_evidence.tsv`; 15/26, not 15/15).
  - Hospital Manta calls INS within ±300 bp for 300/991 HG002 Q100 Alu-head insertions of 250–400 bp (30.3%) (meas., `m14/classify.py`).
  - ME calling is skipped at the hospital (meas.).
- **Corrections.**
  - The SeraCare AluY is an engineered construct, n=1 with unknown AF, and no caller calls it.
  - ClinVar's VCF has 0 symbolic ALTs and 0 INS:ME records, so a clinical MEI is hand-written and can carry explicit sequence today.
  - Explicit-sequence ingest already works.
  - The ME subtype is dropped silently (README:349 says "identically").
  - Most real insertion content (483/765) is local copy, not ME.
- **Proposal.** A warning on random-base INS:ME and length-only `ins:`. Fix README:349. Document the recipe, plus an optional helper writing AluY + poly-A + TSD.
- **Simplest.** The warning and the doc fix.
- **Falsifier.** Gate 0: about 100 HG002-only Alu-head insertions planted into HG001, explicit sequence vs random bases, with whole-genome Manta (not `--callRegions`). The explicit arm must fall within the binomial CI of 30.3%. If the random arm falls within the explicit arm's CI, build only the warning.
- **Gate B.** Gate 0 (about half a day of Manta runs).
- **Cost.** M for the library; 1–2 for the warning.
- **Risks.** L1/SVA need fragments wholly inside the insert.

### 4.3 New capabilities

#### M20 SMN1/SMN2 copy-number events via origin (value 3.5, cost S reframed, kill 0)
- **Problem.** Reframed: the capability already exists under origin. What is missing is the caller-level measurement.
- **Evidence** (meas., `scratch/m20/`, existing binary, hospital HG002, del:chr5:70924941-70953015 af=het, origin mode):
  - Exit 0, SIM_RESIST 0.008.
  - Removed: SMN1 1,868/6,733, SMN2 1,449/6,911, combined 0.243.
  - MAPQ≥20: SMN1 456/970 removed, SMN2 0/1,143.
  - c.840: 12/32 SMN1 C-reads removed, 0/30 SMN2 T-reads.
  - Clean mode refuses (85.6% uneditable).
  - 98.6% of MAPQ-0 SMN1 reads carry a single XA hit at SMN2, which is the twin regime.
- **Corrections.**
  - "Refused by RF8" holds for clean mode only.
  - Clean `--allow-resistant` probably gives a false SMN1=1/SMN2=3 call (inf.).
  - Origin has run on real BAMs, but its locked real-truth validation is still open (REVIEW.md:6385-6394).
  - The hospital HG002 is the SeraCare GM24385.
  - No SMA carrier BAM is on disk. Total_CN_raw at CN4 ranges 3.683–3.947 over 8 samples.
  - Judge per-site values against a CN1 spread, not a CN2 one.
- **Proposal.** No code first. Sequence-informed origin is likely redundant. Gene conversion is a separate item (needs M15).
- **Simplest.** Het and hom SMN1 DELs with today's origin, no `--allow-resistant`, at ≥3 seeds, plus an unspiked control through the same path.
- **Falsifier.**
  - Het: Total_CN_raw ≈2.88–2.92 (±0.15–0.2), SMN1 1, SMN2 2.
  - Hom: ≈1.92–1.95, SMN1 0, isSMA true.
  - Placement after/before ratios: MAPQ0 at SMN1 and SMN2 0.75 ±0.04; MAPQ≥20 at SMN1 0.50 ±0.06, at SMN2 1.00 ±0.06.
  - Clean `--min-mapq 0 --allow-resistant` must fail the placement test.
  - Test origin's 50:50 MAPQ-0 premise on an unequal-CN sample (Seq25-16020/16025 have SMN2=3, at the hospital) or on 1000G CN1 slices.
- **Gate B.** Install SMNCopyNumberCaller; run spike and the caller.
- **Cost.** About 2.
- **Risks.**
  - SIM_DEPTH_FOLD 4.52: the upstream flank holds about 8× raw depth while tiling is scaled to 31.4×, so a coverage bump is likely (inf.).
  - 28 small off-target look-alikes.
  - Origin is experimental.
- **Skeptic notes.** SMN results reach clinicians through the hospital's Scout fork, but the owner has never named SMN, so its demo value is inferred.

#### M48 Per-variant demo receipt (value 3.5, cost M, kill 0)
- **Problem.** No tool scores a spiked run's calls variant by variant.
- **Evidence** (meas., `scratch/m48/`, GIAB v4.2.1 against hospital DV, chr19:10–13 Mb):
  - Het SNV 1,845/1,850; het indel 402/406; hom SNV 793/794; hom indel 199/205; multi-allelic 0/29 by exact key.
- **Corrections.**
  - Spike scripts do read raredisease DV calls (`scripts/duplicates/`, `gvcf/reach.py`), but only to choose sites.
  - The clinical ranked VCF holds 573,084 records, so "present" says little. Rank is driven mostly by annotation.
  - The DV VCF is already split and normalised, FILTER is '.', and it contains './1' GTs.
  - The owner's 2026-10-05 decision puts caller loops, including scoring, in cnv_validation. A scorer in spike needs a new decision.
  - Existing pieces: `score.py` (slice-safety 8b0c4eb), `caller_attribution.py`, truvari in `transplant.py`/`validate_pipeline.sh`.
  - Gate B against score.py is circular.
  - Write `truth.vcf.gz` in addition to `truth.vcf`; about 30 scripts read the plain file.
- **Proposal.**
  - Small variants: `bcftools norm -m-any` plus exact keys.
  - SVs: truvari.
  - `bcftools query` for RankScore and svdb_origin.
  - A per-event footprint diff against the unspiked run.
- **Simplest.** One script; build it once a spiked hospital run is scheduled.
- **Falsifier.**
  - Mutant shams: ALT swap ≤1% called; GT flip ≤2% matched; ≥95% match on normalised multi-allelic and right-shifted indels; random positions ≤1% called.
  - SVs: truvari on Q100 DEL/DUP ≥50 bp on chr19 against the svdb-merged VCF split by origin. A sham-planted SV must come out missed.
  - Realism claims go to M42/M43. A call-rate band has power only for gross failures: at n=100, 60% power at a 2% miss rate (inf.).
- **Gate B.** Correctness shams (minutes).
- **Cost.** About half a day.
- **Risks.** The SM copied into truth touches M38.
- **Skeptic notes.** Fold into M46's scoring step.

#### M18 Repeat-expansion (STR) events (value 3.0, cost M, kill 0)
- **Evidence and corrections** (meas.):
  - HG001 SAMD12 "76 units" is EH noise. Platinum TR truth there is 121/111 bp against a 106 bp reference. RAPGEF2 calls are noise in both samples.
  - EH's in-repeat reads at HG002 RFC1 are all mate-anchored (8 reads, mates MAPQ 60 on chr4). Today's tiling already makes that class, so "ADIR ~0 today" is wrong for alleles up to about one fragment.
  - The excluded share rises with insert length: 4.6% of reads with an in-repeat base for a 383 bp tract, 21% at 520 bp, 74% at 1.94 kb (inf. from the measured TLEN distribution).
  - Only FMR1 and C9ORF72 have off-target regions in the v5.0.0 catalog.
  - EH undercalls RFC1: truth ≈116 units, EH 51 (CI 43–85). Compare EH(spiked) with EH(real), never with the planted count.
  - HTT and ATXN pathogenic sizes are shorter than a fragment and can be planted today as explicit `ins:`.
  - 48/66 EH records differ between the samples.
- **Proposal.** SIM_RU/SIM_REPCN plus a documented `ins:` recipe. Admit novel-only starts only if Gate 0 fails at fragment length or longer.
- **Simplest.** As above (S).
- **Falsifier (Gate 0, no code).**
  - Plant RFC1 and the other long pure STR insertions from Q100 (47 ≥150 bp) into the hospital HG001 resample BAM. Align against the whole genome, never a slice reference. Run EH v5.0.0 with a custom catalog.
  - Pass per locus: ADIR and ADFL within the Poisson 95% band of real; REPCI overlapping real.
  - A random-sequence control must fail.
  - Use EH replicate agreement as the noise floor.
  - Use Platinum's 239 HG001 TR alleles ≥1 kb for the longer-than-fragment half.
  - Check that stranger flags explicit inserts above pathogenic thresholds.
- **Gate B.** Gate 0.
- **Cost.** S for Gate 0 and the recipe; M–L for the full spec.
- **Risks.** GC-extreme CGG/GGGGCC dropout is not modelled. FMR1 on male X needs M12. CR3 replaces the event copy's own STR allele.

#### M19 Mitochondrial handling (value 2.5, cost M, kill 0)
- **Evidence** (meas.):
  - Pileup phasing fires on chrM at the rCRS N placeholder 3107, so the haplotype writes C there.
  - 5,234–5,338 synthetic reads carry C at 3107 and none carries the deletion; 12,563–12,795 kept real reads carry the deletion and none a base. This is 136 bp from m.3243.
  - 19% (35/187) of P/LP MT sites lie within 2 kb of 3107.
  - The homoplasmic 310 T>TC (AF 0.999) is also dropped (CR3).
  - copy_rate(None, v) = v already gives heteroplasmy semantics: planted 0.3, realised 0.291.
- **Corrections.**
  - The proposed minimal mode would emit N at Q2 instead of C, still unreal.
  - A real 20–80% heteroplasmy exists: HG001 16023 G>A at 0.694–0.697.
  - The 7 kb rule allows 3 point events per chrM.
  - The example sites are HG001-only.
  - 6210 T>C is a doubtful yardstick (recurs in another lineage).
  - Truth GT 0/1 already matches Mutect2.
  - `skip_mt_annotation` is true at the hospital.
- **Proposal.** Never call a het at a non-ACGT reference base. Carry homoplasmic indels (both fixes help genome-wide). Use unknown copy and no phasing on chrM. Warn near the chrM ends.
- **Simplest.** As above (cost about 2).
- **Falsifier.**
  - A no-code mixture built by read count, at the evidence level per site: AF, strand, base quality, read position.
  - A read-form check at 3107/302/310. Today's binary must fail it.
  - Plant within 300 bp of 16023 to measure distortion from the copy coin.
- **Gate B.** As above.
- **Cost.** About 2 (rescoped).
- **Risks.** Defer the circular haplotype and the eKLIPse presets.

#### M16 Haplotype-addressed events, hap=1|2 (value 2.5, cost S, kill 0)
- **Evidence** (meas., `hapcheck/kmers2.py`, existing binary, HG001 CRAM with Platinum as `--gvcf`):
  - One SNV at seeds 1–4 lands on hap1, hap2, hap1, hap1. truth.vcf is byte-identical apart from the date.
  - Two LDLR SNVs at seeds 1–6: 5 cis, 1 trans.
  - 98–100% of SPIKE k-mers agree with one Platinum haplotype.
  - 101/153 LDLR exon pairs are ≥7 kb apart.
- **Corrections.**
  - The hospital pipeline cannot see phase: genmod runs without `--phased`, 0 GTs are phased, and every case is a singleton.
  - Trio patterns are already possible at the genotype level.
  - 6 Seq25-160* folders, not 7.
  - Honouring input GT must be opt-in. Refuse cross-PS cis/trans requests.
  - The memory note that `--gvcf` "takes positions only (not phase)" is out of date.
- **Proposal / simplest.** At cost 1: log the chosen copy and write SPIKE_HAP plus a phased GT when the block holds a phased SNP, with no RNG change. Full hap= is S–M.
- **Falsifier.**
  - Gate 0: cis vs trans seeds through DV and genmod (singleton, `--whole_gene`). Pass only if an output differs (predicted identical).
  - Gate 1: per event, ≥90–95% of informative ALT fragments carry the requested haplotype, with ≥5 informative fragments. md5-identical without hap=. Refusals for an unphased gVCF and for cross-PS groups. The master coin build fails with probability 1 − 2⁻²⁰ over 20 events.
- **Gate B.** As above.
- **Cost.** 1–2.
- **Risks.** Allele1/allele2 orientation bugs on 1|0.

#### M15 Grouped events on one haplotype pair (value 2.5, cost L, kill 0)
- **Evidence** (meas.):
  - 52/153 LDLR exon pairs are <7 kb apart.
  - 49.4% of LDLR P/LP SNV pairs are refused (295,948/598,965).
- **Corrections.**
  - At the hospital phase is invisible: 263/338 records near LDLR are already AR_comp, so a single spiked het already pairs.
  - Homozygous FH is a single AF 1.0 event.
  - Under `--allow-overlap` with matching coins: depth 1.25–1.5× and ALT 0.33–0.40. With differing coins, trans comes out nearly right (inf.).
  - A cis close pair can likely be planted today as one MNV (inf., untested).
  - The SIMPLEST trans group erases a real het DEL that DV calls (chr19:11104275 CA>C, QUAL 41).
- **Proposal.** SNV groups on top of M10's in-read editing (about M), with one per-fragment copy assignment.
- **Simplest.** Document the MNV workaround. Use M16 hap= to force trans for pairs ≥7 kb.
- **Falsifier.**
  - Step 0: plant pairs ≥7 kb apart and 300 bp apart (matching and differing coins). Measure AF, depth bins, co-carry and background AF.
  - Use count-based bounds, not "2× a near-zero rate".
  - Whether 1.25–1.5× depth over about 7.5 kb triggers a CNV call is the strongest unmeasured value question, so make it Gate B.
- **Gate B.** As above.
- **Cost.** L (M on M10).
- **Risks.** CR3 is magnified across two remade copies.

#### M17 Trio-consistent spiking (value 1.5, cost XL, kill 1)
- **Corrections.**
  - genmod reads unphased genotypes, and every hospital run is a singleton.
  - Per-member runs already give genotype-level inherited, de novo and compound designs.
  - Father and child BAMs are index-only here; the mother is from another batch.
  - 45% of covering fragments are informative at 60 HG001 het SNVs (meas.).
  - Random inserted bases differ between runs on different BAMs (inf.), so inherited `ins:` needs explicit sequence.
- **Proposal / simplest.** Drop the XL subcommand. Write a README recipe and use M16 plus `whatshap --ped` if phase is ever wanted.
- **Falsifier.** An A/B test: haplotype-consistent vs wrong-copy plants through genmod with the hospital's options. Predicted identical, which would kill it.
- **Gate B.** As above.
- **Cost.** XL as proposed; recipe about 0.
- **Risks.** Family FASTQs are 7-lane (M39).

### 4.4 Algorithms and performance

#### M22 Streamed merge.sh (value 3.0, cost S, kill 0)
- **Evidence** (meas.):
  - B1 merge.sh took 1 h 42 min on 741,918,308 records.
  - `-U` output is single-threaded.
  - On a slice, `view -u -N ^list | merge --write-index` gave the same record md5 in 4.04–4.25 s against 49.8 s. A subset reproduced this (md5 066116b5).
- **Corrections.**
  - The cause is the compression of `-U`. The one-pass form `view -@T -u -N list -o removed.bam -U - ORIGINAL | merge …` took 0.82 s with the same md5 and needs no `^`. Its guard fires after the merge.
  - Keep BAM output for CRAM input.
  - Disk is not binding (10 TB free).
  - `-N ^` works in samtools 1.21–1.23; README "≥1.13" needs updating.
  - RG collision renames differ only with the stand-in.
  - Value is for BAM-route loops, not the FASTQ demo route: about 1 h 45 min becomes about 10 min per event (inf.).
- **Proposal.** Stream the selection into merge with `--write-index`; drop the re-sort; take the total from idxstats; write under a temp name and mv merged.bam and its .bai.
- **Simplest.** As above, plus rewriting the stub-samtools tests (main.rs:3626-3800).
- **Falsifier.**
  - On the chr20 2.7 GB BAM with a real sim.bam, then on B1: same record multiset, non-decreasing (tid,pos), total = orig − removed + sim, BAI magic.
  - The wrong-BAM guard still aborts. A mutant without `^` fails.
  - ≤20 min, cold cache.
  - Report the in-order md5 separately.
- **Gate B.** Done on the slice.
- **Cost.** S.
- **Risks.** Tie order on non-samtools-sorted inputs. Check that sim.bam's @HD says SO:coordinate.

#### M21 Coverage-preserving remap, DUP-scoped (value 3.5, cost L as proposed / ≈2 scoped, kill 0)
- **Evidence** (meas.):
  - 1 Mb chr20 DUP: 63,390 pairs suppressed, 120,046 tiled, SIM_DEPTH_FOLD 4.69, SIM_ALT_FRAGS 17, 19.2 s, 704 MB.
  - LDLR B3 DUP: synthetic depth flat at about 296 reads/kb whatever the donor depth; the thin bin gets ratio 1.73 against 1.5.
- **Corrections.**
  - The 4.69-fold bin (chr20:14862000-14863000) is a real GC 0.126 dip: 103 reads against 225–317 in neighbouring bins, 82/103 at MAPQ≥20. The proposed mappability filter would drop it.
  - Depth shape is diluted, not discarded.
  - The BAM route re-emits 66,621 kept pairs.
  - 81.09× is from the Codex base binary; master's probe gives fold 3.88. The sham cannot run through the CLI (AF must be >0).
  - SIMPLEST's ownership rule is wrong at the DUP ends. A tandem DUP is purely additive.
  - The performance claim is overstated.
  - DEL interiors are already local. SNVs belong to M10.
- **Proposal.** DUP only: keep every original; add one jittered depth copy per event-copy interior pool fragment (the existing junction mode); add junction fragments at about v·D_local.
- **Simplest.** As above.
- **Falsifier.** Slope of added reads on donor depth, as in §3. Report the any-MAPQ undershoot in poorly mappable bins as a known CR4-like limit. Boundary bins are a separate locked check.
- **Gate B.** As in §3.
- **Cost.** About 2 scoped.
- **Risks.** Copies come only from the MAPQ≥20 pool. LDLR value is modest: B3 was found by Manta through its junction, and CNVnator missed it anyway.

#### M23 Parallel whole-sample FASTQ filter (value 2.0, cost M, kill 0)
- **Evidence.**
  - 2 M pairs: 11.95 s at 8 threads (logged).
  - Per-mate filters plus a cmp over name FIFOs: 8.10 s vs 19.47 s on a loaded machine, and 2.75 s vs 4.5–5.5 s on 500k pairs, md5-identical, catching a swapped pair at record 100001 (meas., `m23check`, `m23v`).
- **Corrections.**
  - The baseline is about 38–42 min, not 33 (inf.).
  - At the default `--threads 4`, compression is the bottleneck, so stage 1 gains about 1.09× at most.
  - Stage 2 needs BGZF block copy; without it, about 1.6× over stage 1 (fails its own rule). Its cost is L.
  - The chr20 FASTQ is a re-extract, and the perf input is not BGZF.
  - The scratch awk streamed only kept names; every record's name must be streamed.
  - A raredisease run takes hours, so the saving is modest.
- **Proposal.** Stage 1 only.
- **Simplest.** Same.
- **Falsifier.**
  - Value gate at the hospital first: time fastq.sh, check the raw R1 for BGZF, and build only if fastq.sh is ≥10% of spike-to-VCF time and THREADS ≥8.
  - Compressed md5 identical at THREADS 4 and 16.
  - The 13 existing fastq tests pass.
  - Mutants (no cmp; kept names only) fail.
  - Fault injection leaves no partial output.
  - Interleaved timing medians.
- **Gate B.** As above.
- **Cost.** S for stage 1.
- **Risks.** mawk hosts are slower. Drop slot interleaving: appended reads were measured harmless, 1,296/1,296.

#### M25 Per-record cost in the duplicate scan and LOH pileup (value 2.0, cost S, kill 0)
- **Evidence** (meas., `m25v/run1.log` with millisecond stamps, `m25prof/` gdb sampling):
  - Pass 2 0.34 s; duplicate scan 2.38 s (7× pass 2); LOH pileup 1.94 s; tiling 4.33 s; startup sample about 9.7 s.
  - Inflate is a minority of each scan (2/22 samples in the duplicate scan).
- **Corrections.**
  - The original log had whole-second stamps.
  - The labels mix steps.
  - The duplicate scan does a full RecordBuf conversion per record.
  - The "simplest" filter saves at most about 45% of record building (54.5% of records pass), so it misses its own ≥50% gate.
  - perf is blocked here.
  - The one-pass plan is deprioritised, not refuted.
- **Proposal.** A lean duplicate scan (flags and borrowed name on the lazy record, no XA). A dense Vec pileup. Optionally the noodles libdeflate feature as a separate byte-identical A/B.
- **Simplest.** As above.
- **Falsifier.**
  - Each step falls to ≤1.5× the pass-2 interval.
  - User time drops ≥2.5 s on a 3 Mb DEL at 1 thread, A/B/A/B.
  - md5 identity on the 15 sets, a CRAM set and a hom event.
  - Mutants (skip duplicates; ignore qc_fail) fail.
- **Gate B.** Done.
- **Cost.** S.
- **Risks.** Small gains; ranks after M22/M23. The decision "no pass-2 strips, no one-pass" is free, because neither was on the owner's backlog.

#### M26 Reference-cache peak memory (value 2.0, cost S, kill 0)
- **Evidence** (meas., `scratch/m26/`):
  - Load peaks are about 3× one chromosome, and 2.07× for 23 chromosomes (6,122,108 kB at 2,890 MB loaded). RSS falls to about 2.9–3.0 GB after the load.
  - The cache never gets a hit in production.
- **Corrections.**
  - Flagged twice before (REVIEW.md:2834; independent review :139).
  - Streaming outputs would break md5 identity, cross-event dedup and origin; drop it.
  - The kept-pair clone is about 9% of one big event (inf.).
  - Peaks for big events and chrM come from the pool, not the reference.
- **Proposal / simplest.** Delete the cache. Optionally uppercase in place.
- **Falsifier.** MAXRSS/loaded ≤1.25 (vs 2.02); no load spike above the plateau plus the largest chromosome; md5-identical on the 23-SNV panel and on 3 Mb DEL/DUP at 1 and 4 threads; a mutant that keeps the cache fails.
- **Gate B.** Done.
- **Cost.** S.
- **Risks.** None found.
- **Skeptic notes.** Value is low; an FH panel peaks at about 1.1–1.4 GB (inf.).

### 4.5 Robustness and correctness

#### M28 ClinVar-native VCF ingest (value 3.5, cost S, kill 0)
- **Evidence** (meas., `scratch/cv/`, `skeptic_ldlr39.txt`, `skeptic_gw.tsv`):
  - ClinVar's '19' contig is refused.
  - 251409 is labelled 'SNV', `sim_var_1`, SIM_GENE unknown.
  - 251140 via VCF plants 313 pairs with SIM_ALT_FRAGS 0; the coordinate spec plants 1,540 pairs with fold 1.51.
  - Adding SVTYPE alone makes R1 byte-identical.
  - 0 of 6,329 ClinVar records ≥50 bp carry SVTYPE.
  - LDLR: 31 DEL, 3 DUP, 3 INS, 2 Indel.
  - Genome-wide P/LP ≥50 bp: about 2,600 Deletions (2,602–2,603 depending on the query), all pure; 283/285 Duplications are exact adjacent copies.
- **Corrections.**
  - 2 of 3 LDLR insertions contain N runs and would become Q2 N runs; genome-wide 353/448 Insertions and 107/305 Microsatellites carry N.
  - 3767245 (8,118 bp, no anchor) must be refused. 430754 (a 63 bp delinsA) must stay on the small path.
  - DEL reads are nearly the same either way, but validate treats a record with no SVTYPE as a SNP, so del_planted never runs.
  - Microsatellites must be refused until M18 (168 of 198 plantable ones are DMPK).
  - Routing drops the REF check (main.rs:1634-1643), so with a contig alias a GRCh37 record would plant at the wrong place. Keep a REF and DUP-copy check.
  - RF8's 0/35 covered the DELs plus 2 DUPs.
  - Re-run every gate on a master build.
- **Proposal.**
  - Map CLNVC to DEL/DUP/INS for records ≥50 bp without SVTYPE.
  - Refuse and count the rest, naming the reason.
  - Add SIM_SRC_ID, SIM_GENE (GENEINFO) and SIM_HGVS.
  - Warn on coordinate DEL/DUP/INV specs <50 bp: the cnv_validation BED holds the 4 LDLR DUPs as 1 bp intervals.
  - The contig alias is optional, comes with the REF check, and refuses when both spellings exist.
  - Drop `--vcf-ids` and HGVS specs.
- **Simplest.** CLNVC mapping, SIM_SRC_ID and the refusals.
- **Falsifier.**
  - Routing census of 3,996 records against the locked counts (2,602 DEL; 283 DUP; 95 INS; the rest refused).
  - md5 identity for 251409, 4845403, 251140, 976329, 403639 and 1074682 against their specs.
  - A 500-record sample <50 bp keeps R1/R2 md5-identical.
  - Negative controls are refused: shifted 251140, a GRCh37 record, a mutated DUP copy, a FASTA with both contig spellings.
  - Replace the 66 h whole-BAM gate with validate on the existing B2/B3 merged BAMs.
- **Gate B.** The classifier over the 39 LDLR records (2 already done).
- **Cost.** S.
- **Risks.** CLNVC=Duplication for 403639 is left-aligned.

#### M36 Provenance and stale outputs (value 3.0, cost S, kill 0)
- **Evidence** (meas.):
  - B1 README reads 'spike version: 0.1.0' with an unquoted command.
  - A paste test with an echo stand-in passed only `--event …`, then bash printed "`--seed`: command not found" after the run with seed 42 into /tmp/spike.
  - R47: a refused second run leaves run 1's truth.vcf (a5538b6f2fa9) and R1 (84851d3db490) in place.
- **Corrections.**
  - Exit code 2 is clap's usage code; use 3.
  - RUN_OUTPUT_FILES already includes sim.bam. It lacks merged.bam(.bai), merge leftovers and the whole-sample pair.
  - Any failure after `create_dir_all` (main.rs:489) leaves stale output, not only refusals.
  - Refusing a dirty `-o` breaks validate_pipeline.sh:757-759 and repeated default runs; delete spike's own names instead.
  - `##spike_*` lines in truth would break the cross-binary md5 gates (rf8_compare strips only fileDate), so put provenance in the README.
  - Every harness already uses `rm -rf` and checks exit status; the realistic victim is a hospital `--into-fastq` loop.
- **Proposal / simplest.**
  - build.rs with `git describe --dirty` and rerun-if-changed on .git/HEAD and the branch ref.
  - README: commit, the command quoted through sh_word, the parsed Args.
  - Delete spike's own file names at startup; write truth via temp file and rename.
  - Exit 3 for refusals.
  - Defer run.json, run id, @CO, `--overwrite` and the input hash.
- **Falsifier.**
  - Paste test from another directory, after checking the chosen site is not refused as carried.
  - Stale `--into-fastq` case, plus a SIGKILL case: no truth.vcf may remain.
  - validate_pipeline step 3, slice_loop on 2 events, and two default runs all still pass.
  - The version changes after an empty commit.
  - Cross-binary md5 with a documented strip rule.
  - Exit codes: refusal ≠1/2; missing flag = 2.
- **Gate B.** Done for the stale case.
- **Cost.** S simplest, M full.
- **Risks.** Tarball builds have no git.

#### M30 `--gvcf` sample identity (value 2.5, cost S, kill 0)
- **Evidence** (meas., `robust/gv2_*`, `wrong_gvcf.py` re-run):
  - HG001's VCF on HG002 35x: exit 0, 0 WARN. At chr20:45000494, 15/15 SPIKE reads carry REF where HG002 is 1/1. Depth there goes 29 → 47.
  - 64.5% of windows get ≥1 wrong het. 195/678 hets are shared. 13.5–14% of windows lose a hom-alt.
- **Corrections.**
  - The hard-coded column was noted at REVIEW.md:2835.
  - Trio multi-sample VCFs are hypothetical: all 20 local raredisease VCFs are single-sample with names matching their BAM.
  - The LDLR test footprint came out md5-identical with the wrong VCF.
  - GIAB names never match lab SMs, so a name rule can only warn.
- **Proposal.**
  - Pick the column whose name equals SM. Refuse a multi-sample VCF without a match. Warn on a single sample with another name.
  - A run-wide pooled concordance over het and hom-alt sites (startup blocks plus footprints, ≥50 sites). Refuse below a safe floor (<0.5) with an override.
  - Drop the per-footprint refusal.
- **Simplest.** The column and name part.
- **Falsifier.**
  - Correct pairs ≥0.95; the wrong pair ≤0.5 (about 0.33 predicted).
  - Both the chr20 test run and the LDLR run with the HG001 VCF are refused, run-wide.
  - Correct runs stay byte-identical.
  - A merged HG001+HG002 VCF picks the HG002 column.
  - A sites-only VCF is refused. The hospital D24/D25 swap is caught.
- **Gate B.** bcftools mpileup at about 100 footprints' sites.
- **Cost.** S (S–M for the pooled gate).

#### M31 Refuse a BAM aligned to another build (value 2.5, cost S, kill 0)
- **Evidence** (meas.):
  - CHM13 BAM with the GRCh38 FASTA: exit 0, 0 WARN, 4,483 hom-alt SNPs in 6 kb, 69.3% of the quality sample masked (0.22% on HG002).
  - samtools merge also accepts the LN mismatch silently.
- **Corrections.**
  - CRAM is already protected by noodles' slice MD5.
  - An M5 refusal would falsely refuse UCSC hg38 BAMs: 18/24 primary contigs have a different M5 at equal LN; on chr19 the difference is only the centromere N.
  - The alt-aware warning already exists (INFO, since the M32 fix).
  - The harness Delly check is still needed.
  - "Event contigs only" misses chrM.
  - Masked share is WARN at most.
- **Proposal / simplest.** An LN check over all contig names the input and the .fai share; refuse naming the contig and both lengths.
- **Falsifier.**
  - CHM13 is refused before extraction.
  - HG002 35x, hospital HG002, NA18488 and the HG001 CRAM pass md5-identical (0/195 LN differences).
  - Mutants fail: check skipped; full set equality; a chrM-only fixture.
  - An M5-only change on an NA18488 header is not refused.
- **Gate B.** Done.
- **Cost.** S (about 15 lines).

#### M33 Warn on duplicate-unmarked input (value 2.5, cost S, kill 0)
- **Evidence** (meas.):
  - The strobealign BAM has no markdup @PG and 0/149,768 flagged reads, against 12.10% on 35x. Re-marking finds 13.09%.
  - Two more unmarked BAMs (novoalign and bwamem2 chr20) have only 2.8% 5'-end families.
- **Corrections.**
  - The 0.535 bias is an upper bound. With per-pair removal, about 0.516–0.532 at v=0.5 (inf. from measured families).
  - NA18488 is 9.16%, not 8.47%.
  - The hospital BAM is marked.
  - The quality sample drops duplicates before any count, so a detector must count while skipping.
- **Proposal / simplest.** One warning reporting the pool's family share and the "up to" bias. No refusal, no flag.
- **Falsifier.**
  - Paired design: a 35x slice against its twin with flags stripped. Mechanism ratio ≈1.12–1.14 ±0.03.
  - 300–500 SNVs: harm is real if Δ ≥+0.012 with a CI excluding 0, dead if the upper bound is <+0.008.
  - Depth and DUP readouts.
  - The detector fires on the unmarked BAMs and stays silent on marked and remove-duplicates BAMs.
- **Gate B.** As above (minutes).
- **Cost.** S.
- **Skeptic notes.** Low value for the demo.

#### M27 Per-event random streams (value 2.5, cost M, kill 0)
- **Evidence** (meas.):
  - R39: adding event A 13 kb away changes all 276 of B's SPIKE_ records, SIM_ALT_FRAGS 22 → 19 and suppressed 292 → 278, with an identical pool.
  - R49: only 11/273 sequences are shared.
  - Tiling of a 3 Mb DUP takes 17.20 s at 1 thread (doc., threads-step1).
- **Corrections.**
  - Refused transplant attempts write nothing, and the B1/B2 arms re-draw under any key.
  - K7 is already handled by v2's seed-spread allowance; master flipped 0/115 over 5 seeds.
  - "Panel equals singles" holds only in clean mode with non-interacting pools.
  - It needs the owner's word, because the byte-identity rule was set for speed-ups (deferred 2026-09-27).
- **Proposal / simplest.**
  - Seed each event's StdRng from splitmix64(seed, FNV-1a(canonical spec)), with three sub-streams: suppression, placement, synthesis.
  - Origin gets its own stream.
  - Names and IDs unchanged. Drop stage 2 and stable names.
- **Falsifier.**
  - Value check first: an 8-event LDLR panel at 3 seeds; count threshold-crossing flips.
  - Invariance across {B}, {A,B}, {B,A}, {A,B,C}, with A inside the 10 kb flank, over every RNG consumer. Master and an index-seeded version fail.
  - K7-style equivalence at seeds 42, 2, 3, 4, 5.
  - A synthesis-only change leaves the suppressed sets identical.
- **Gate B.** Code reading (largely done).
- **Cost.** About 2 simplest.
- **Risks.** Every output changes once, and baselines must be re-cut.

#### C4 Reword the RF8 `--min-mapq` hint (value 2.0, cost S, kill 0)
- **Evidence.**
  - Already measured: on an exact twin, clean mode at `--min-mapq 0` gives L 0.271 / P 1.296 against truth 0.767–0.770 / 0.812–0.815 (doc., REVIEW.md physics test). README:516 documents it.
  - The refusal text (census.rs:232-238) and README:507 still promote lowering `--min-mapq`.
- **Corrections.**
  - Only PRODH crosses RF8's 0.5 in the 400-exon sample. CYP21A2 (0.106; 0.068 on 35x) is never refused.
  - PRODH's look-alike is chr22_KI270734v1_random.
  - An XA-majority gate fails at the README's own chr20:27.1 Mb locus: 5,322 MAPQ-0 reads, 1 with XA.
  - The kill rule is ill-posed for diverged paralogs.
- **Proposal / simplest.** Always attach the caveat to the hint, and name `--edit-model origin` as experimental (checked only on an exact twin). README bullet and table row, test update (census.rs:483), and optionally a truth header line when `--min-mapq` <20.
- **Falsifier.** A unit test on the message. No Gate B run.
- **Cost.** XS.

#### M29 Exon-spec leftovers after the MANE fix (value 2.5, cost S, kill 0)
- **Evidence.**
  - Already done locally: exon 1 = 11089462-11089615 and exon 18 end 11133820 (1103835, merged in ba9f4e7, with a MANE test). Public 4a5672d still has the old exon 1.
  - Open items (meas.):
    - the known-deletions VCF has round-number breakpoints named like published alleles; only one record matches ClinVar (251409, shifted 193 bp);
    - "gene not found" writes a 13 MB log (5.1 MB error line, 37,491 WARN lines);
    - a one-line awk MANE→BED recipe reproduces the bundled BED and gives 19,062 genes with correct minus-strand numbering.
- **Corrections.**
  - 4 ClinVar records sat in the old exon 1, not 5.
  - The Scout BED is coding-only and merged across transcripts. Do not support it.
  - The README gVCF example misses MANE exon 1.
- **Proposal / simplest.** Push. Relabel or remove the demo VCF. Cap the near-match list at about 10. Add the README recipe and a note that exon specs give idealised edge breakpoints. `--exon-gff` only if needed (M).
- **Falsifier.** Recipe checks: genome-wide strand-order numbering, the CDS start in exon 1, a DMD exon45-50 run.
- **Cost.** S.

#### M34 N-aware tiling (value 1.5, cost S, kill 0)
- **Evidence** (meas.):
  - At chr20:61000-61400 with `--allow-resistant`: 123/644 SPIKE reads are all N, 149 at least half N. RF8 refuses by default (0.814).
  - Near N runs: 6/65,142 PanelApp green exons (SHANK2, C1R, GRK1); LDLR 0.
- **Corrections.**
  - The cited sampler serves only INS haplotypes, so the simplest fix misses the DEL case.
  - Contig ends are already clamped.
  - Depth over real sequence is not inflated.
  - Excluding whole fragments is unrealistic: 51 real pairs span a 100-N gap.
  - The N-rate gate is unsound.
  - At C1R, RF8 would let the event through.
  - fastp's default n_base_limit (5) drops these reads on the hospital route.
- **Proposal / simplest.** After the tiling loop, drop pairs with a read holding >5 N. Warn on ≥10 bp of reference N. Leave N bins out of SIM_DEPTH_FOLD.
- **Falsifier.** Reads with >5 N go 149 → 0; depth bins within Poisson noise of master; fixtures with isolated N and 9-bp runs unchanged; identity sets byte-identical.
- **Cost.** About 1.
- **Risks.** 94 IUPAC positions in the reference need their own check.

#### M35 EventStat consolidation, with the audit replaced (value 2.0, cost S, kill 0)
- **Corrections.**
  - EventStat already exists (main.rs:1823-1854).
  - The four name invariants hold by construction and are unit-tested.
  - RF1 was a debug-vs-release assert defect, not an index misalignment ("a real run cannot trip them").
  - The proposed mutants cannot test an in-memory audit.
  - run.json does not exist.
- **Proposal.** Add alt_frags to EventStat when the event loop is next touched, under a byte-identity gate. Instead of the audit:
  - (a) count SPIKE_ pairs before and after dedup_by_name; a mutant with duplicated name prefixes must refuse;
  - (b) have fastq.sh check ADDED1==ADDED2 and R1/R2 name order before the raw pass (<10 s).
- **Falsifier.** As stated in each item.
- **Cost.** S.
- **Skeptic notes.** A separate inferred risk: re-spiking an already-spiked sample may collide SPIKE_ names. Test it.

### 4.6 Usability and clinical deployment

#### M37 Fail fast on the wrong raw FASTQ (value 2.5, cost S, kill 0)
- **Corrections** (meas.):
  - The header check misses the likeliest mistake: hospital D24-14230 and D25-7403 share RUN:FLOWCELL:LANE (LH00352:45:227NC2LT1:2).
  - The SIMPLEST version falsely refuses multi-flowcell samples (the HiSeq HG002 BAM has 5 RUN:FLOWCELL values).
  - R1/R2 mismatches already fail at record 1.
  - 42 min is predicted, not measured.
- **Proposal.** A name-overlap probe of the first ~1M raw records against the 50k names bam_stats samples. Measured: 51 hits for the right pair, 0 for a wrong one, 1.7 s; hospital expectation about 67 (inf.). And/or a fastq.sh deadline after M = ceil(ln(1e6)·N_est/k) records (about 15–80 s at the hospital, inf.).
- **Simplest.** The probe; warn rather than refuse on coordinate-ordered or subset FASTQs.
- **Falsifier.** D24 BAM with a D25 FASTQ is refused in <2 min. The HiSeq multi-flowcell pair passes. The right pair stays byte-identical. A mutant without the deadline fails.
- **Cost.** S.
- **Skeptic notes, spin-off.** Hospital FASTQs hold exactly 760,000,000 and 800,000,000 reads ("30x_resample"). A downsampled BAM given its full-depth FASTQ would silently dilute VAF (inf.). Add an END check of the raw pair count against the BAM's pairs (hospital ratio 1.03).

#### C5 Hospital-installable spike (value 2.5, cost M, kill 0)
- **Evidence** (meas.):
  - The dev binary fails to start on glibc 2.31 and 2.35 containers ("GLIBC_2.39 not found", exit 1). That container run overrides a symbol-table reading that put the floor at 2.34.
  - The build is not pure Rust: bzip2-sys and lzma-sys are C, and the binary links liblzma.so.5.
- **Corrections.**
  - fastq.sh uses no samtools or bwa-mem2. Only `mv -T` fails late. BAM-route tools fail at their first command.
  - The hospital @PG tools ran in apptainer.
  - Data tests are already `#[ignore]`; CI needs bcftools.
  - A different libm could change RNG-driven paths, so a portable build must be md5-identical to the dev build.
- **Proposal.** Gate B first. Then one portable binary (musl with static lzma, cargo-zigbuild with old glibc, or an apptainer image), with SHA256. Make fastq.sh fail early on `mv -T`. CI (`cargo test --locked` with bcftools, MSRV 1.82) is separate hygiene.
- **Simplest.** The host check plus a README "installing at the hospital" paragraph.
- **Falsifier.**
  - The portable binary runs `--into-fastq` end to end on the host.
  - md5-identical to the dev build on the identity sets; within about 10% wall time.
  - The dev binary must fail on rockylinux:8.
  - A stub `mv` without `-T` makes fastq.sh exit before it reads the raw file.
- **Gate B.** 10 minutes on the host: copy the binary and run `--version`; `ldd --version`; `command -v` for the tools; network to crates.io and github.
- **Cost.** M for the full package, XS for the check.

#### M38 Share-safe outputs, reduced (value 2.0, cost S, kill 0)
- **Corrections.**
  - The README already says R1/R2 hold kept originals (2,435 kept vs 228 synthetic in B1).
  - The align.sh and sim.bam hits come from SM on purpose (README:489).
  - The symlink test resolves through canonicalize.
  - truth `##reference` and bwa's @PG carry absolute paths.
  - Script permissions protect nothing: the -o directory itself follows the umask.
  - D24-14230 is the SeraCare GM24385/HG002 sample, not NA12878; its IDs are already public and the owner accepted that (doc., github-push.md).
  - With a patient background, about 91% of R1/R2 are real reads, and the read names carry flowcell IDs. A "share-safe" label would mislead.
- **Proposal / simplest.**
  - Warn when -o is the /tmp/spike default or already holds spike outputs.
  - An opt-in that leaves paths out of the README, truth and script defaults (the scripts then require $1).
  - A README warning that bundles from a patient sample must not leave the lab whatever their labels.
  - Drop `--portable` SM relabelling.
- **Falsifier.** A path-token and SM-token grep over all outputs on a real file under a FAKEID directory; scripts still run when given arguments.
- **Cost.** S.

#### M39 Multi-lane FASTQ route (value 2.0, cost M, kill 0)
- **Corrections** (meas.):
  - The 6 CEPH samples are 7-lane, but the mother NA12878 is single-lane. All were run as singletons with empty parent IDs.
  - There is no LB in the read groups, MarkDuplicates sees "Unknown Library", tagging is DontTag, and DV ignores read groups. So a lane layout can change little beyond bwa batch noise.
  - Under the concatenation workaround, the "SPIKE_ name claims a missing lane" problem does not arise.
- **Proposal / simplest.** A recipe: run fastq.sh on the concatenated pair, split the output by the RUN:FLOWCELL:LANE fields of each name, and write a 7-row samplesheet (S).
- **Falsifier.**
  - Whole-genome comparison: per-lane vs concatenated, against a null of two concatenation orders. Build the full feature only if call differences exceed the null.
  - A runtime check: build only if splitting saves ≥3 h end to end. Per-lane BWA ran 58–65 min in parallel, against 5 h 45–6 h 06 for a single pair (meas.), but whole-run spans are confounded.
- **Cost.** M full, S recipe.

### 4.7 Validation methodology

#### M41 Background-integrity yardstick (value 3.5, cost S, kill 0)
- **Problem.** No yardstick checks what a spike does to the sample's own variants in the remade event ±2 kb.
- **Evidence** (meas., `valid_lens/cr3_bg.py` re-run, `m41_skeptic/`):
  - 297/754 forward footprints hold ≥1 HG002 non-SNP.
  - Hom indels: median 0.934 → 0.479, 187/198 below 0.75 afterwards (28 before).
  - Het indels: median 0.448 → 0.240; 109 dropped by >0.2.
  - 214/752 events damage ≥1 background indel.
  - Sham: 0/198 hom and 0/233 het changed.
  - The damage does not depend on size (1 bp 0.932 → 0.500; 5–19 bp 1.000 → 0.474) and is milder at 1.5–2 kb (0.937 → 0.588). There is no synthetic share beyond 2 kb.
  - LDLR: 27 HG002 non-SNP records within ±2 kb of the 18 exons (6 hom).
- **Corrections.**
  - Footprints falling from ≥0.75 to <0.75: 121 (159 records), not 137.
  - Hospital DV 1.6.1 already calls these 1/1: 159/159 in that set; overall 197/200 hom 1/1 and 230/234 het 0/1. So the caller-level gate will likely confirm the harm.
  - The SNP control capped at about 400 SNPs.
  - T7 was "supported", and a size threshold for its warning has no basis.
- **Proposal.** A standing read-AF column with the sham null in every realism round. A DV caller column in the hospital acceptance run (M46). It becomes the locked falsifier for M10, M11 and M21.
- **Simplest.** The read-AF column (3 s per 750 events).
- **Falsifier.**
  - Pass rule for any CR3 fix: the count of hom indels below 0.75 lies within a binomial 95% margin of the sham.
  - Caller gate: ≥50 events drawn at random with a fixed seed. DV on the after view and on the sham. Kill if the after-minus-sham GT-change rate is <0.05 per event, or if the sham shows the same rate.
  - Mutant: emptying SampleCopies must flip ≥90% of hom SNPs (today 0/198).
- **Gate B.** As above (about 1 h with the cached DV 1.9.0).
- **Cost.** S.
- **Risks.** The crude carrier detector. The round 2b binary is old, so re-measure on master.
- **Skeptic notes.** One skeptic's before-call lookup used 41 min of CPU, over budget.

#### M42 DeepVariant caller-equivalence round (value 3.5, cost M / S first step, kill 0)
- **Evidence** (meas., `m42chk/rates.py`, reproduced):
  - Forward het call rates: SNV 100, DEL1-4 99, DEL5-19 100, DEL20-49 95, INS1-4 98, INS5-19 97, INS20-49 78, DUP 24/54.
  - Shared real-vs-real discordance: 0–6 per 100.
  - Median |dGQ| 3–6.
- **Corrections.**
  - DUP misses are mostly genuine no-calls (23/30 forward have nothing nearby), so normalisation recovers little.
  - GLnexus keeps DV's GT and GQ (revise_genotypes false).
  - Shared sites are easier (89/100 called vs 78/100 forward), so the yardstick is biased against spike.
  - At real-uncalled forward INS20-49 sites the donor shows 0/21 sites with ≥2 carriers against 9/21 for fake; at called sites they match. Stratify by donor support.
  - Fisher at n=100 has no power (78% vs 87%, p=0.136). Use paired McNemar with the full pools.
  - GQ touches rank only through ModelScore bins <10 and 10–20.
  - B2 for small variants must be new.
- **Proposal.** The same DV binary on donor, recipient, spiked, sham, old fe15a46 and shared arms. Per-group McNemar with a three-way rule. Ceiling groups judged on dGQ/dAD-VAF and GQ-bin crossings. Pre-set donor-evidence strata. Record background GT changes (M41) for free.
- **Simplest.** DV 1.9.0 on round 2b's existing views for INS20-49 and SNV, seen-not-judged (2–3 h).
- **Falsifier.**
  - The sham reproduces the recipient's calls exactly (locked stop).
  - The old fe15a46 views (read-mean quality SD 0.248× real) must shift GQ, otherwise "no power".
  - A new B2 must lose the call.
- **Gate B.** DV 1.9.0 against the raredisease 1.6.1 calls on donor sites.
- **Cost.** S first step, M full.
- **Risks.** Version mismatch. Expect little power for SNV and DEL/INS 1–19.

#### M43 Manta equivalence at real DELs, stratified by Alu (value 3.5, cost M / ≈2 reduced, kill 0)
- **Evidence** (meas., `m43_reality/rates.py`, `valid_lens/`):
  - Real Manta rates per size bin: forward 49/65/76, reverse 48/73/71; any caller 50/78/86 and 49/80/86.
  - Alu-Alu DELs ≥300 bp: Manta 34/77 (counting SVDB Intersection records), some caller 62/77. No-Alu 187/257; one-Alu 64/66.
  - By direction: HG001 forward 13/37, HG002 reverse 19/40. At shared Alu-Alu sites, HG002 calls 6 that HG001 misses and 0 the other way, which is a library effect.
  - Spiked 4845403: missed by Manta FULL, found by TIDDIT and CNVnator FULL.
- **Corrections.**
  - Events in one run are ≥100 kb apart, not 7 kb.
  - Runs 1–2 hold 1–10 events each.
  - The old sim.bam has 148 bp reads and predates several fixes.
  - The callRegions route was probed but not locked (10 extra INS calls on chr19).
  - The Fisher rule cannot detect the idea's own example (p=1.0).
  - Caller loops belong in cnv_validation.
- **Proposal.** Re-plant forward run0, reverse run0, B1 run0 and a sham on current master (about 30 min). Build merged BAMs and run Manta `--callRegions`.
- **Simplest.** As above.
- **Falsifier.**
  - Paired Δ = ½[(spiked in HG002 − real HG001) + (spiked in HG001 − real HG002)], per bin, with Alu-Alu pooled, bootstrap CI against a locked δ of 0.10–0.14.
  - Continuous alt PR+SR per depth.
  - B1 must fail the continuous metric (expected log ratio ≈ −0.69). The sham adds 0 DEL calls.
  - Skip TIDDIT until its noise is measured.
- **Gate B.** As above.
- **Cost.** About 2 (3–4 BAM builds; faster with M22).
- **Skeptic notes.** The LDLR question is mostly answered: spiked B1 behaves like the typical real Alu-Alu DEL.

#### M46 Whole-sample acceptance for the hospital route (value 3.0, cost L / 2–3 simple, kill 0)
- **Corrections** (meas., `m46_order/`, `m46value/`):
  - The TIDDIT "~2,400" was FULL runs on four different BAMs. Manta (8,864 calls) and CNVnator (2,546) were identical away from the event on the BAM route.
  - Order sensitivity comes from bwa's tie-break hash on input position (bwamem.c:553/1224/1230), not from insert size. Moving 1% of pairs changed 622/17,198 records, all MAPQ 0. Rerunning the same input changed 0.
  - fastp 0.22.0 with `--thread 6` reorders about 71% of reads per run. 0.23.4 is untested.
  - The published baseline lacks per-caller SV VCFs.
  - One raredisease run takes about 14 h.
  - fastp on SPIKE_ reads is benign (C3).
- **Proposal (simple).**
  - Per sample: one spiked run with:
    - a survival audit against SIM_ALT_FRAGS (fail any event below the real 5th percentile);
    - at-event calls;
    - an outside-footprint BAM diff excluding the mate loci of removed pairs (provisional ≤0.01% changed, inf.);
    - a broken-build control on the local F1 stand-in.
  - Once per pipeline version: one unspiked rerun with per-caller VCFs as the noise floor.
- **Simplest.** As above.
- **Falsifier.** As above. An empty-removal-list build must fail the at-event and diff checks, with an old-allele share above 0.0017 at hom SNVs.
- **Gate B.** The fastp 0.23.4 order test.
- **Cost.** 2–3 simple.
- **Risks.** Watch SMNCopyNumberCaller, CNVnator and EH for MAPQ-0 re-draws.

#### M40 Read-realism panel and "can you tell" classifier (value 2.5, cost S, kill 0)
- **Corrections.**
  - The hospital positive control (0.958 vs 0.58) compares spike with HG001 reads over 13.7× larger windows, which is not composition-controlled. Its top features differ from K2b's.
  - Pre-v2's 0.958 rests on sdQ and meanQ, the gap K1 already caught. So far the classifier has found no new gap.
  - K2b AUC 0.571–0.576 over seeds 21/22/23; nulls 0.495–0.505. Removing the last-cycle features moves AUC only 0.571 → 0.564, so "fix the top feature lowers AUC" cannot fail reliably.
  - Pre-v2 binaries cannot build K2b sets.
  - The error table has no mate key to swap.
  - Bands from a real-vs-real split run at half depth.
  - v2 already fails most rows.
  - The harness drops trimmed reads.
- **Proposal.** Commit the phys scripts and disc.py with one runner reporting distances with seed SDs.
- **Simplest.** As above (2–3 h; runs in about 1 min).
- **Falsifier.**
  - Gate on change: the target row's distance shrinks by >3 seed SDs, and no other row grows by >3 SDs.
  - Mutants v2 currently passes must go red: both mates using mate 1's tables; templates forced to run bin 0.
  - A caller-facing row: DV-rule candidates per Mb (real 500.1, spike 80.7–92.7).
  - The classifier serves only as a regression guard with a fixed random_state.
- **Cost.** S.

#### M44 Re-run the transplant yardsticks on master (value 2.5, cost M, kill 0)
- **Corrections.**
  - Rounds 2b/3 were already re-run with the trimming binary 858f13f (38 m 41 s): 7 supported cells, not 5. Only round 1b is still from the 148 bp binary.
  - The transplant metrics did not react to the read-length change (0/6 cells), so the falsifier "separates the 148 bp binary" is vacuous.
  - The "future plans cite this" claim refers to roadmap ideas, not committed plans.
  - round2.py raises on non-RF8 refusals.
  - Baselines live in the volatile scratchpad: move them first.
  - SEED is hard-coded (transplant.py:49).
- **Proposal / simplest.** A one-off no-harm re-run on master at seeds 1 and 2, plus 858f13f at seed 2.
- **Falsifier.**
  - Round 3's rule; any supported → refuted cell is a regression.
  - The paired rank rule: "new wider" in ≥4 of 6 junction/indel cells, p ≈ 0.009.
  - Then make "transplant no-harm" a standard K in plans that change read synthesis. No per-merge suite.
  - Cheaper pre-merge guard: an input-matrix smoke test (decoy CRAM, chr-less names, CSI-only index, 31-value BAM). It would have caught M32.
- **Cost.** 40–80 min of compute.

#### M45 Dose sensitivity of the transplant verdicts (value 2.5, cost M, kill 0)
- **Corrections.**
  - The B1 lists are stale.
  - round3.py cannot be used unchanged.
  - An offset ladder makes sense only for DUP J.
  - A non-monotone curve is expected where spike sits above real.
  - The design is already paired.
  - The 35x cross-library comparison must stay separate.
  - Dropping E after the fact would flip verdicts and hide the RF15 signal.
  - The rho bias threshold does not depend on n.
- **Evidence** (inf., interpolation from existing normal and B1 arms, `m45check/`, `skeptic_mde/`):
  - "Supported" bounds the dose error only to about 20–35%.
  - Paired signed bias, seen not judged, 28 cells: forward SNV A −0.025 [−0.049, −0.002]; reverse DEL5-19 A −0.045 [−0.072, −0.016]; reverse DEL5-19 E −0.052 [−0.080, −0.025]; forward INS20-49 E +0.059 [+0.029, +0.090].
- **Proposal / simplest.** Step 1, no spike runs: an in-silico dose curve plus a locked, Holm-corrected bias rule (±0.05 A or ±20% dose) beside rho. Step 2, optional: VAF 0.40 and 0.60 arms on one binary with locked predictions.
- **Falsifier.** Observed mean d inside its interpolated interval, and monotone across arms.
- **Cost.** S for step 1.

#### M47 Composite benchmark truth (value 2.5, cost M, kill 0)
- **Corrections.**
  - truth.vcf already has contig, ALT, INFO and FORMAT headers; only bgzip/tabix is missing.
  - events.bed is ±10 kb, not the footprint. The footprint is span ±2000.
  - Removing footprints from the confident region drops the planted windows; add ±50 bp back.
  - Keep background SNPs (2,129 of 2,576 run0 footprint records).
  - Use Q100 v1.1 with the v5.0q bed.
  - Mask SeraCare.
  - bedtools, hap.py and rtg are absent.
  - The concept is a documented review recommendation (CLINICAL_SV_REVIEW.md:131-133).
  - Demo-scale artefacts are about 0.4 per event (inf.).
- **Proposal / simplest.** truth.vcf.gz plus .tbi with SM (shared with M48). Spike writes planted.bed and footprints.bed. A README bcftools recipe; truvari for SVs; no normalisation inside spike.
- **Falsifier.**
  - Composer correctness: GIAB used as the query gives FN = planted-in-confident, FP = 0, footprint non-SNPs UNK.
  - A mutant without the exclusion re-admits 159 hom indels.
  - Keep masking only if naive spurious FP+FN per event is >0.1 after M41's DV slices.
  - The "rest unchanged" claim goes to M46.
- **Cost.** S simplest.

#### M49 Per-read truth file (value 2.0, cost S, kill 0)
- **Evidence** (meas., `skeptic_m49/gate_b.py` on ld1/B1, hospital HG002):
  - Flank reads: 438/438 within 10 bp of truth.
  - Junction reads: 7/18 within 10 bp. 11 are off by exactly 2,135 bp through about 26 bp of near-identical AluSx. Those are equivalent representations, mostly a 1-bp deletion plus a short clip and no SA.
  - Manta had about 10 MAPQ-60 discordant pairs but no candidate, so the miss is in the caller (inf.).
- **Corrections.**
  - PairSpans holds haplotype coordinates only.
  - locate.py is chr20-only and DEL-only, and independent by design.
  - DEL haplotypes have no reverse-complement segments, so the RC mutant cannot fail.
  - A 10-bp scorer is unsound at homologous junctions.
  - Attribution is already possible from read names.
- **Proposal / simplest.** Optional TSV of haplotype spans, segment-to-reference mapping, source copy and ALT flag. Or a generalised standalone locate.py.
- **Falsifier.** Value gate: build only if some documented miss changes attribution with per-read truth. Validate the writer on INV, DUP, INS and fusions. Any scorer must be homology-aware.
- **Cost.** S (TSV).

#### C2 dup_planted / inv_planted rows (value 3.0, cost M, kill 0)
- **Evidence** (meas., replica of split_reads, approximate):
  - Real HG002 tandem DUPs ≥300 bp on 35x pass split_reads 9/38.
  - NA12878 INVs pass 31/39 on the CRAM (8/36 fail at ≥300 bp in a second peek).
  - Spike's own events pass: LDLR B3 DUP 6 joining reads; chr20 DUP 12, INV 32/25, BND 3.
- **Corrections.**
  - Real INVs come from Platinum CEPH (HG001), not HG002 (Q100 has 0); use END=POS+SVLEN.
  - round2's DUP candidates stop at 299 bp.
  - Windows overlap below 1 kb (a 106-bp INV passed 38/38), so a null side is needed.
  - INV is not a demo target, and no BND truth exists.
  - The advisory line is validate.rs:622.
- **Proposal / simplest.** Measure the demo DUPs first. Then dup_planted with RF14's matcher. inv_planted only after the real-INV comparison. Drop BND.
- **Falsifier.**
  - (A) 4 LDLR ClinVar DUPs at af 0.5 and 0.1 with 3 seeds. Build no dup_planted if ≤1/12 fail at 0.5.
  - (B) HG001 INVs with SUPP ≥2, spiked into HG002 at the same site.
  - (C) Kill arms: END shifted 1 kb, the unspiked BAM, a DUP relabelled INV.
- **Cost.** M.

#### C3 Hospital read-model gap (value 2.5, cost S, kill 0)
- **Evidence** (meas., `c3fastp/`, re-scored existing K2 BAMs):
  - Bin-0 bad-end clips: spike 14/5,430 vs own 46/6,313 (z −3.57), and −3.38 after Q100 masking.
  - About 80% of the gap sits in Alu-overlapping reads (0.42% vs 1.45%, z −3.58; non-Alu z −1.15).
  - About 60% of the bin-0 excess is in 150-bp reads.
  - The fastp double pass is closed: SPIKE_ reads drop 0.48% as low quality and 0.08% as too short, with 0.015% of bases corrected; originals 0.
- **Corrections.**
  - The harness drops trimmed reads (only 35.9% of pairs are both 151 bp).
  - The kill rule cannot fire, because 35x already fails bins 1 and 3.
  - Ranking M1–M9 by the clip gap is ill-posed for 5 of them.
  - Crash share differs within Alu strata.
  - Counts are small, and 16/22 events were chosen as Alu-flanked.
- **Proposal.** Adapt the harness for trimmed libraries. Then a value gate: count clusters of ≥2 bad-end clips within ±5 bp in real Alu flanks. Under 1 per 10 kb means no build is justified.
- **Falsifier.** Kill a "hospital-specific gap" if the Alu-stratum ratio's CI overlaps 35x's range (0.90–0.97).
- **Cost.** S–M.

#### C7 SeraCare background note (value 1.5, cost S, kill 1)
- **Evidence.** Mostly done:
  - The manifest is on disk (cnv_validation/seracare_inherited_cancer_truthset_v1.0, 23 rows) (doc.).
  - Transplant masks already cover the 7 genes ±2 kb plus IGL.
  - The falsifier was run: 28 green-exon non-Q100 het calls, 23 of them in the 7 SeraCare genes (meas., `scratch/c7/`).
- **Corrections.**
  - The "230 pathogenic" figure is mostly baseline (NA12878 has 187).
  - A GQ≥30 filter would miss manifest variants; use a region mask.
- **Proposal.** One sentence in the M41/M46/M47/M48 specs and in the demo notes: use count.sh's mask or the truthset, or use Seq25-7598 as the demo background. Do not hard-code the vendor manifest in spike.
- **Cost.** About 0.

---

## 5. Killed or already done

- **M24** Reuse the quality profile across runs. Killed: Gate B fails on existing logs, about 0.5–2 min saved per transplant campaign against a 10 min bar. Slice loops never hit a file-keyed cache, and the design repeats a stale-cache failure (case file 2026-09-28). Side issue: slice-route runs learn quality only from the slice.
- **M32** Decoy CRAM abort in the quality sample. Done: fixed in 1103835 … ba9f4e7 (C1 exit 0, 10/10 md5-identical). The only action left is to push it.
- **C1** Large-DUP realism from CEPH carriers. Killed: its own Gate B yields n=1 against the 20 needed. Truth DUPs are mostly satellite or VNTR with flank identity median 0.956, while spike plants exact copies. Reframe it as a comparison with confirmed carriers at the hospital.
- **C6** R_c < 1 noise estimator. Killed: Monte Carlo bias 1.037–1.068 is under the 1.1 kill line, and analytic noise changes 0/6 unsure cells. The INS20-49 A anomaly remains as a separate question.
- **fastp second pass over SPIKE_ reads.** Closed by measurement (C3).
- **LDLR exon 1 in intron 1** (M29 part 1). Fixed locally, not pushed.

**Dropped sub-proposals** (do not re-propose):
- M23: stage 2 and slot interleaving.
- M26: streaming outputs.
- M35: the four-invariant audit.
- M38: `--portable` SM relabelling and script-umask modes.
- M17: the XL trio subcommand.
- M27: stage 2 and stable names.
- M12: GT '1' default, `--strict-ploidy`, `--ploidy-bed`.
- M2: clip-boundary mask and mate bit.
- M1: fixed 40% overlap threshold.
- M7: BAM-route copies.
- M8: any default-on build.
- M9: hard-coded TruSeq adapter.
- M13: Alu-first SIMPLEST snap.
- M14: bundled ME library before Gate 0.
- M15: both-copies remake before M10/M11.
- M19: emitting N at 3107.
- M22: CRAM-in/CRAM-out.
- M25: pass-2 strips and one-pass reads (deprioritised).
- M31: M5 refusal and masked-share refusal.
- M34: excluding whole fragments over N.
- M36: refusing a dirty `-o`, and exit code 2.
- M37: a RUN:FLOWCELL-only check.
- M40: using the classifier ranking to choose work.
- M44: a per-merge suite.
- M45: the full ladder, and dropping E.
- M47: normalisation inside spike.
- M49: a 10-bp scorer.
- C4: an XA-majority gate.
- C5: `--check-tools`.
- M10: the fastp one-mate control.
- M41: the outcome-selected 137 footprints and the AND kill rule.

**Refuted claims** (do not repeat):
- The 4845403 Manta miss is due to low-MAPQ junction reads.
- The 35x 15+ slip rate of 2.2% is all slips.
- Real LDLR DEL ends are Alu-enriched.

---

## 6. Open questions for the owner

1. **Push.** Should local ba9f4e7 (plus docs 3211d55 and f18c7f5) go to public `spike`? Public 4a5672d aborts on decoy-aware CRAMs and ships the old LDLR exon 1.
2. **Where scoring lives.** Keep caller loops and scoring in cnv_validation (the 2026-10-05 decision), or allow a generic scorer in public spike (M48, M43 runner)?
3. **Byte identity.** May outputs change once for per-event RNG streams (M27) and for read-model fixes (M2 Q2→N, M5)?
4. **In-place edits.** May the rule "spike never edits a real read in place" (README:1427) be relaxed for SPIKE_-renamed in-read SNV edits (M10)?
5. **Adapter clips (M9).** In scope, or does the earlier "bad ends only" choice stand?
6. **Demo background.** Seq25-7600 (HG002 SeraCare mix, needs masks) or Seq25-7598 (plain NA12878)?
7. **Demo content.** Are X-linked genes (DMD etc.) in the demo (M12)? Is SMN in scope (M20)? Are trio or family demos wanted (M15–M17, M39)? Are low-VAF or somatic uses in scope (M1, M3, M4)?
8. **Hospital-side facts and compute.**
   - Host glibc, network, apptainer, and whether it can build from source (C5 Gate B).
   - fastp 0.23.4 order test and one extra raredisease rerun (M46).
   - Does the lab hold MLPA- or array-confirmed DUP carriers (C1 reframing)?
9. **Pipeline config.** Tell the hospital that raredisease runs with `par_bed: null` and calls male PAR haploid (M12)?
10. **T7 / CR3 option B.** The pending decision on a footprint warning. M41 shows the damage does not depend on size, so a size threshold has no basis.
11. **Origin mode.** Point users to `--edit-model origin` in the RF8 message as "experimental" (C4)? Use exit code 3 for refusals (M36)?