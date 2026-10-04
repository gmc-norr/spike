# Keep the gVCF's hom-alt calls when it has no het call in a region

**Asked 2026-10-04.** This is review finding 4 (P2, `docs/review/2026-10-04-independent-review.md`), picked by the user after fixes 1, 7, 3, 8 and 5. We reproduced it on master `2764dbc`.
- The probe has a hom-alt SNP at 8 reads' depth in an event's replacement flank. Given a gVCF holding only that hom-alt call, the synthetic reads held **5 REF and 0 ALT** there.
- With one unrelated het record added to the same gVCF, they held **0 REF and 7 ALT**.

## Read in the code

`loh::sample_copies` (`src/loh.rs:165`) gets the sample's SNPs over the event's footprint. For an SNV that is the event plus 2 kb on each side.
- With `--gvcf`, it reads the gVCF's single-base SNV calls there: het, and hom-alt. Then:
  - **When the gVCF has at least one het call there,** it uses the gVCF alone.
  - **When it has none,** it logs `no het SNPs from gVCF, trying pileup fallback` and throws the whole gVCF result away, hom-alt calls included. The pileup then calls SNPs on its own: depth at least 10, het when the top two bases are each 20-80%, hom-alt at 90% or more.
- So a gVCF hom-alt the pileup does not call (fewer than 10 reads at MAPQ 20, or under 90%) is lost. The synthetic reads then carry REF where the sample is hom-alt.
- The fallback has been in the code since the first commit, with no recorded reason. The README's review notes give it one job: a gVCF that names the chromosome differently (`20` against `chr20`) gives nothing for the region, and the pileup stands in.

## Design (locked)

When the gVCF has no het call in a region, the pileup still runs as now, and the gVCF's hom-alt calls are then added on top of the pileup's result. Where the two disagree, the gVCF wins:
- every gVCF hom-alt call goes in as hom-alt with the gVCF's allele, replacing a pileup hom-alt with another allele at that position;
- a pileup het at a position the gVCF calls hom-alt is removed;
- everything else the pileup found stays as it is.

The gVCF-with-hets path and the no-`--gvcf` path do not change. A gVCF with no call at all in the region (the renamed-chromosome case) adds nothing, so that case does not change either. The info line says how many hom-alt SNPs came from the gVCF.

This is the smallest change that keeps the gVCF's calls. It leaves alone the pileup's hets in such regions, which master already uses to phase.

**Tests, written first and seen red:**
- the merge, on hand-made sets:
  - a gVCF hom-alt the pileup lacks is added;
  - a pileup het at a gVCF hom-alt position is removed;
  - a pileup hom-alt with another allele there takes the gVCF's;
  - pileup hets and hom-alts elsewhere stay;
- `sample_copies` on a small BAM with 8 reads carrying T at one position, and a plain VCF calling that position `1/1` with no het calls: both copies carry T there. On master they carry nothing there.

**Mutation checks** (each must turn a test red; the unmutated tests are green first):
1. the merge not called;
2. the pileup's het kept at a gVCF hom-alt position;
3. the pileup's allele kept over the gVCF's;
4. the pileup's other hets dropped;
5. the pileup's other hom-alts dropped.

## Checks (locked)

**R (reported, not judged): how often this bites at the hospital.** Use the hospital BAM and its own DeepVariant calls (`GM24385_seracare_cancer_snv.vcf.gz`, the same raredisease run, given as `--gvcf`). Draw 500 random 4,001 bp windows on chr20 (seed 1). Report:
- the share of windows with at least one hom-alt SNV call and no het SNV call;
- in those windows, the share of hom-alt calls that spike's pileup rule misses.

The rule is applied in Python: primary, not duplicate, not QC-fail, MAPQ at least 20, at least 10 reads, and the top base at 90% or more and not REF. These misses are the sites where master writes REF onto synthetic reads.

**H1: the reviewer's probe.** `scripts/review_20261004.py` (the scratch copy from the fastq-safety and later results) against the new debug binary:
- `hom_only` gives 0 REF and at least 5 ALT;
- `hom_and_het` stays at 0 REF and 7 ALT;
- every other result is as before.

**H2: unchanged where nothing is lost.** On the hospital BAM, `--seed 1 --threads 16`, with one SNV event, the new and master binaries give identical `truth.vcf`, `R1.fq.gz`, `R2.fq.gz`, `replaced_reads.txt` and `fastq_removed_reads.txt`:
- (a) without `--gvcf`, at the first of R's windows;
- (b) with `--gvcf`, at the first of R's windows that has a het call;
- (c) with `--gvcf`, at the first of R's windows that has hom-alt calls and no het call, where the pileup misses none of them.

**H3: changed where it should.** With `--gvcf`, at the first of R's windows where the pileup misses a gVCF hom-alt, at that site:
- every synthetic read (`SPIKE_`) from the new binary that covers it carries the ALT;
- master's do not.

If R finds no such window, H3 uses a window found by widening the draw to the whole of chr20, scanned in order. That is said in the result.

## Addendum, before H2 and H3 ran: R's output and how the sites are chosen

**R (run, `scripts/gvcf/reach.py`, output `scratchpad/gvcf/r.log` and `r.tsv`).**
- Of 500 windows: 336 have a het call, 88 (0.176) have hom-alt calls and no het call, and 76 have no SNV call.
- In the 88, the pileup misses 17 of 377 hom-alt calls (0.045), in 11 windows.
- 15 of the 17 lie in chr20:26.7-30.0 Mb, with 0-7 reads at MAPQ 20. The two others are chr20:4767379 (26 reads, not 90% ALT) and chr20:53717174 (7 reads).

**Site rule for H2(c) and H3, fixed now.** Before H3 runs, it is already clear that the first missed sites hold 0-2 reads, where spike may refuse to run at all.
- Each window's event is an SNV at the window's centre (start + 2,000, 0-based), so its footprint is the window. If the centre is one of the sample's calls, the event moves 1 bp right until it is not.
- REF is the reference base. ALT is the next of A, C, G, T after it.
- H3 takes R's missed sites in order. It uses the first at which the master binary runs without refusing, and lists the ones it skipped and why.
- H2(c) takes the first window with hom-alt calls, no het call and no miss, under the same refusal rule.

## Result (2026-10-04): supported

Code `7d6e605`. Binaries: new release `920d2b7c`, new debug `b553389e`, master `658baa4` release `93eb71b1`. Scratch: `scratchpad/gvcf/` (`h_run.py`, `h_run.log`, `h/`).

**Tests.** 701 pass (2 new), 0 fail, 1 ignored. The mutation runner checks the unmutated `loh::`/`gvcf` tests are green first (27), then catches **5 of 5** mutants. clippy gives the same warnings as master.

**R** is in the addendum: 88 of 500 windows are affected, and the pileup misses 17 of their 377 hom-alt calls.

**H1: PASS.** In the reviewer's script (scratch copy), `hom_only` gives **0 REF and 5 ALT** (master: 5 REF and 0 ALT), and `hom_and_het` stays at 0 REF and 7 ALT. Every other result is as in the validate-base-depth result.

**H2: PASS.** Each case gave the same 5 files as master, byte for byte:
- (a) without `--gvcf`, at R's first window, with `snp:chr20:10019032:A:C`;
- (b) with `--gvcf`, at the first window with a het call. That is the same window, and the same event;
- (c) with `--gvcf`, at `snp:chr20:5237028:C:G`, the first window with hom-alt calls, no het call and no miss, where master ran.

**H3: PASS.** R's missed sites, taken in order:
- The first 8 (chr20:26.7-28.2 Mb) were skipped, because master refused their events. Four times this was RF8. Three times the donor reads were too few: 3, 5 and 6 pairs, against the 30 needed. The plan's rule says to skip those. So at those sites, spike does not run an event at all.
- At the 9th, `chr20:4767379 A>G` (window 181, event `snp:chr20:4768371:A:C`), the pileup sees 26 reads, not 90% G, so it misses the call. There, every one of spike's 13 reads that cover the site carries **G** with the new binary. With master, all 13 carry **A**.

**What this means at the hospital.** With `--gvcf`, a hom-alt call that the pileup misses no longer turns into REF on spike's reads. In R's sample that happened at 17 of 377 hom-alt calls in affected windows. Most of them lie where spike refuses to run anyway, but not all: chr20:4767379 is one where it does run. Without `--gvcf`, nothing changes.
