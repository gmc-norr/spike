> Copied unchanged from the detached run's `/home/parlar_ai/spike-codex-run/DESIGN-NOTES.md` (2026-09-25). The `STATUS.md`, `TASKS.md` and `scratch/` it names are in that folder, outside git.
>
> **Done since the copy** (the notes below are left as the run wrote them): the first three
> items of its suggested order. **CR7(c)**, `8dca8d0`: `af=het` is exactly 0.5. **CR3's
> fail-closed half**, `f08228d`: an unreadable `--gvcf` stops the run, with no flag; a failed
> pileup still only warns, and the footprint scan is not done. **CR6 relabel**, `a1ba176`: a
> warning, the help line and the README say what fusion mode is; no rename. Details and
> measurements are in `REVIEW.md` under each CR.

# DESIGN-NOTES — the Codex findings that change spike's model or its defaults

These six are **not** fixed in this run. Each changes what spike simulates, or what it does by
default, and TASKS.md reserves those decisions for the human. Every number below was measured in
Phase 1 on the base binary `4efa0f4`, with `scripts/review_sv_model.py`, and the raw output is
in `STATUS.md`.

A note on the recommendations: they are mine, made with the measurements in hand and nothing
else. Where I think the honest answer is "restrict what spike claims" rather than "build the
model", I have said so, because a smaller correct tool is worth more here than a larger one
whose limits are undocumented.

---

## CR2 — donor-aware depth

### The problem

Spike measures fragment depth **once**, in a 2 kb window at the first covered breakpoint
(`src/simulate.rs:480`), and applies that single number across the whole variant haplotype
(`src/simulate.rs:548`, `:725`).

Measured: a heterozygous tandem duplication of `chrT:10000-28000` over a donor whose interior
section runs at **18.75x** makes that section **81.09x** — a **4.32x** rise, where a locally
proportional CN2->CN3 predicts **28.13x**. The 75x section of the same duplication becomes
**110.92x** against an expected 112.5x, so the error is confined to the section whose depth
differs from the anchor's. On a uniformly covered donor both sections behave as expected
(110.85x and 109.81x against 112.5x).

The rule spike implements inside a duplication is
`D_out(x) = (1 - v)·D_donor(x) + 2v·D_anchor`, not `(1 + v)·D_donor(x)`.

### The governing principle

Dosage must be a local property of the donor, because that is what a depth-based caller reads.

### The options

**Option A — project a learned spatial fragment-start intensity onto the altered haplotype.**
Learn a binned start-rate along the donor, then generate each synthetic fragment from the
intensity at the reference position its start maps back to, rather than from one scalar.

**Option B — refuse the claim instead of fixing the model.** Keep the single estimate, but
measure the donor's depth profile across the event's footprint and **reject or loudly label**
any event whose footprint depth varies beyond a threshold against the anchor window.

**My recommendation: A, with B shipped first as a gate.** B is days of work and immediately
stops spike from silently producing the 4.32x error; A is the real fix and can land behind it.
Doing B first also gives A its acceptance test for free: the events B rejects are exactly the
ones A must get right.

The review's own caution applies to A and should be written into its design: do **not** treat
low aligned depth as a molecular sampling bias. Some of it is mappability, and it should
reappear on its own when simulated molecules are aligned. A naive intensity learned from
aligned depth would bake mappability into the molecules and then double-count it.

### The measurement that decides it

**Works:** on the `dup_variable` probe, the 18.75x interior section reads within sampling noise
of **28.13x** (and the 75x section stays near 112.5x). Then the same, binned, on the real HG002
chr20 BAM: output/donor depth ratio per 1 kb bin across the event footprint should be flat at
`1 + v` in well-mappable bins, rather than rising where the donor is thin.

**Does not work:** the binned ratio is flat on the synthetic probe but not on HG002 — that
would mean the intensity is learning mappability rather than library sampling, and option A has
become a way to hide the problem rather than model it. A second falsifier: a **control** run
with no event must come back with a binned ratio of 1.0 everywhere; if the intensity model
perturbs an unedited genome, it is wrong.

### Rough size

B: small — one pass over the donor pool, a threshold, an error, a flag to override, and tests.
A: large. It touches the count formula, placement, and the depth-copy path, and it needs a new
learned object with its own serialisation and tests. Expect it to be the biggest single item on
this list.

### What it breaks

B changes accepted input: events over heterogeneous coverage start being rejected. That will
bite panel and exome users hardest, and anyone spiking near a segmental duplication.
A changes the emitted reads for every event on a non-uniform donor, so every existing
regression baseline moves. Both need a release note.

---

## CR3 — the sample's own indels and SVs, and failing closed

### The problem

`SampleCopies` stores **one base per reference position** (`src/loh.rs:34`,
`HashMap<u64, u8>`), so an indel cannot be represented at all. The gVCF path drops non-SNP
alleles outright (`src/loh.rs:513`: `if ref_allele.len() != 1 || alt_allele.len() != 1 { return; }`).
Synthetic haplotypes are built from the reference and then patched with SNP alleles
(`src/simulate.rs:156`), so the sample's own indels and rearrangements are replaced with
reference sequence.

Measured: a donor carrying a **homozygous** 2 bp deletion, with a heterozygous tandem
duplication spanning it, comes out at **AF 0.360** (32 deletion-supporting reads against 57
reference-supporting). All three resulting copies should still carry it; the two synthetic
copies resurrect the reference bases.

Separately, when the sample-variant input cannot be read at all, `src/simulate.rs:93` logs a
warning and proceeds with **empty** copies — so a failed load is indistinguishable in the
output from a sample with no variants.

### The governing principle

A synthetic copy must be built from the sample's sequence, not from the reference; and where it
cannot be, spike must say so rather than quietly substitute the reference.

### The options

**Option A — a real sample-haplotype model.** Replace `HashMap<u64, u8>` with an ordered
variant list per copy (SNVs, indels, and existing SVs), built from phased calls or an assembly,
and construct each haplotype by applying it to the reference before the event is composed.

**Option B — fail closed, and label.** Keep the SNP-only model, but (i) make an unreadable
sample-variant input a hard error under a strict mode, and (ii) scan the event footprint for
non-SNP variation in the supplied calls and reject or explicitly label any footprint containing
it.

**My recommendation: B now, A when there is a real need for it.** B is cheap and removes the
silent wrong answer, which is the actual harm: today a user gets a truth VCF that looks fine.
A is a large piece of work whose value depends on having phased, sequence-resolved sample calls
to feed it — if the lab does not have those, A buys nothing. Do A when a concrete use case
supplies the input.

The fail-closed half of B is the part I would not defer at all. It is small and the current
behaviour is indefensible on its own terms.

### The measurement that decides it

**Works (A):** the `dup_homdel` probe returns to **AF 1.0** (today 0.360), and a heterozygous
background indel under the same duplication comes out at the fraction its copy count implies
rather than diluted. On real data: hom-alt SNPs *and* indels inside an event footprint keep
AF 1.0 in the merged BAM.

**Works (B):** an unreadable `--gvcf` exits non-zero with a message naming the file, and an
event whose footprint contains a non-SNP call in the supplied VCF is rejected by default.

**Does not work:** `dup_homdel` reaches AF 1.0 but phasing breaks — the background indel and
the event end up on copies that contradict the input phase. That is the failure mode to watch,
and the test for it is a phased het SNP plus a phased background indel in one footprint: their
haplotype assignment in the output must match the input's.

### Rough size

B: small for the fail-closed half, medium for the footprint scan (it needs the variant reader to
report non-SNPs it currently discards).
A: large, and it overlaps CR7's ploidy work — do not design them separately.

### What it breaks

B changes accepted input twice over: runs that today proceed on an unreadable gVCF will stop,
and events over known indels will be rejected. Both are currently silent wrong answers, so the
breakage is the point — but it will fail existing pipelines, and it needs a flag and a release
note.
A changes the emitted sequence wherever the sample differs from the reference, so every
baseline moves.

---

## CR4 — training eligibility versus replacement eligibility

### The problem

One filter decides two different things. Only proper pairs passing the MAPQ and flag filters
enter the donor pool (`src/extract.rs:104`, `:777`), and that same set is the complete
population of molecules the event can edit. Everything else — low-MAPQ pairs, discordant pairs,
unmapped mates, flagged records — survives the event untouched in the merged BAM.

Measured: with half the donor pairs at MAPQ 0, a deletion requested at **AF=1** leaves
**37.5x** inside the deletion, exactly the MAPQ-0 half of a 75x input. That is not a homozygous
deletion. `spike validate` reports `coverage_ratio expected 0.00, observed 0.00, pass` — because
its own default `--min-mapq` of 20 hides precisely the reads that survived.

### The governing principle

Which molecules are good enough to *learn a library model from* is a different question from
which molecules the *edit applies to*, and using one answer for both makes the altered genome
depend on the aligner's confidence.

### The options

**Option A — separate the two populations.** Keep the filtered pool for training, but track
**all** molecules overlapping the footprint — their mates and their supplementary and secondary
records — and decide explicitly what the altered genome does to each. Emit an exclusion census
alongside the truth, quantifying event-resistant depth.

**Option B — census and gate only.** Do not change which molecules are edited. Measure the
resistant fraction inside every event footprint, write it into the truth VCF and the run README,
and fail the run when it exceeds a threshold.

**My recommendation: B, then A only for the classes with a defensible answer.** A sounds right
and is partly a trap: for a multimapper there is no correct answer without modelling assignment
uncertainty, and blindly deleting one would be worse than leaving it. B makes the problem
visible and quantified today, which is what a validation study actually needs. Then extend A to
the cases that are unambiguous — a proper pair failing only on MAPQ, and a mate outside the
query — and declare the rest unsupported at that locus.

The review says the same thing in its own words: *"Avoid blindly deleting unrelated
multimappers: assignment uncertainty should be modeled or the locus declared unsupported."*

### The measurement that decides it

**Works (B):** on the `del_lowmap` probe, the run reports a resistant fraction of **0.50** and,
by default, fails rather than writing an `AF=1` truth record beside 37.5x of surviving reads.

**Works (A):** the same probe's deletion interior reads **0.00x** with `--min-mapq 0` on the
merged BAM, not just with the default filter.

**Does not work:** the interior reaches 0.00x but reads elsewhere in the genome have gone
missing — deleting a multimapper's other copies. The falsifier is a whole-BAM read-count and
per-chromosome depth comparison against an unedited control: only the footprint may change.

### Rough size

B: medium. The census is a second pass over the unfiltered records in the footprint, plus a new
truth field and a gate.
A: large, and it interacts with CR2 (which molecules exist) and CR3 (what sequence they carry).

### What it breaks

B's gate rejects events at difficult loci that run today. `spike validate` should also learn to
report the resistant fraction rather than inheriting the default MAPQ filter that hides it —
otherwise the tool that is supposed to catch this remains blind to it.

---

## CR6 — balanced translocation

### The problem

Fusion mode is additive. `is_additive` is true for `SimEvent::Fusion` (`src/simulate.rs:144`),
which short-circuits all suppression, and the count is
`n = coverage · v/(1-v) · breakpoints.len()` (`src/simulate.rs:584`). At `v = 0.5` that is `C`
junction fragments added on top of `C` retained originals. `from_fusion`
(`src/haplotype.rs:345`) builds **one** join with flanks — not a derivative chromosome pair.

A balanced heterozygous reciprocal translocation should alter one copy at each participating
chromosome, retain the other copies, produce **both** derivative adjacencies, and conserve copy
number. Spike does none of those four. Two BND mate records describe the two ends of *one*
adjacency; they are not the reciprocal derivative.

Local breakpoint evidence is therefore much stronger than in a copy-neutral germline sample.

### The governing principle

A balanced rearrangement conserves copy number, so a model that only adds evidence cannot
represent one.

### The options

**Option A — derivative-chromosome paths.** Let an event name both derivatives explicitly, with
dosage and phase, and suppress/replace across both partners so copy number is conserved.

**Option B — relabel, and stop claiming otherwise.** Keep the additive mode, rename it in the
CLI and the docs to what it is — additive junction evidence — and state plainly that it does not
establish sensitivity or genotyping for balanced translocations.

**My recommendation: B now, A only if translocations are in the validation scope.** B costs a
day and removes a claim that is currently false. A is a genuine model extension and should not
be started until CR1's grouped-event composition exists, because a reciprocal translocation is
by definition two linked events that must be generated once, together — exactly the problem CR1
identified and this run only gated rather than solved.

The review adds one more thing worth acting on cheaply: **the legacy junction DUP mode shares
the additive dosage problem.** Use the full tandem model for simple DUP work.

### The measurement that decides it

**Works:** on a synthetic two-chromosome probe with equal coverage, a balanced heterozygous
translocation leaves both partners at **1.0x** their donor depth across the whole span (copy
number conserved), while both derivative adjacencies carry junction evidence at the fraction the
genotype implies — roughly half the local depth, not `C` on top of `C`.

**Does not work:** the adjacencies appear at the right strength but a depth ratio somewhere on
either partner has moved off 1.0 — the rearrangement is no longer balanced, and the model has
turned a translocation into an unbalanced one.

### Rough size

B: small — a rename, docs, and a warning.
A: large, and **blocked on CR1's grouped-event composition**. Sequence it after that.

### What it breaks

B renames a mode, which the run's rules treat as a CLI change, so it needs an alias and a
deprecation period.
A changes what a fusion event emits entirely; every fusion baseline moves.

---

## CR7, the rest — input GT, ploidy and copy number, and what `af=het` means

Three separate decisions that happen to share a finding. The two bounded parts of CR7 — the
insertion sequence and the AF caps — are **fixed** in this run (tasks 4 and 5).

### The problem

**a) Input GT is not the simulated genotype.** Measured: an input DEL with `GT=1/1` and no
simulation-specific AF produces `SIM_VAF=0.500; GT=0/1`. `genotype_from_vaf`
(`src/truth.rs:321`) is `if vaf >= 0.9 {"1/1"} else {"0/1"}` — the genotype is a function of AF
alone, and the input GT is never read. FILTER and homozygous-reference records are not
selection gates either.

**b) No ploidy or absolute copy number.** Both the simulation and the emitted GT assume two
starting copies. A male non-PAR X or Y event cannot have its haploid truth represented; baseline
CNVs and CN>4 gains have no model. `AF=1` obtains complete deletion of *eligible* reads (see
CR4) but still writes diploid `1/1`.

**c) `af=het` changes the underlying event fraction, not the observation.**
`Beta(40,40)` (`src/main.rs:430`) is sampled once into `resolved_af`, which drives both copy
suppression and generation. A value above 0.5 makes the second copy contribute. That is a dosage
mixture, not a fixed diploid heterozygote observed with sampling noise.

### The governing principle

The truth file must describe the genome that was built, and the genome must be described by a
genotype and a copy number rather than inferred from a read fraction.

### The options

**Option A — a genotype-and-ploidy model.** An event carries a genotype and a baseline copy
number; ploidy comes from a sex/karyotype declaration or a contig-ploidy map; the emitted GT is
what was built, and the read fraction follows from it rather than the reverse. `af=het` becomes
`GT=0/1` with sampling noise at the observation, not a moved event fraction.

**Option B — two explicit modes.** Keep the AF-driven path as an explicitly named
*specification* mode, and add a *genotype-reproduction* mode that reads input GT, honours FILTER
and hom-ref, and refuses events on contigs whose ploidy it has not been told. The review asks
for exactly this distinction.

**My recommendation: B, which is A arriving in a shippable order.** The two are not really
rivals — B's genotype-reproduction mode *is* A, restricted to where it can be supported, and
naming the existing behaviour honestly is most of the value. The AF path is legitimate for
titration work and should not be removed.

Do (c) first regardless of which is chosen: separating `af=het`'s Beta from the event fraction
is small, and until it is done every "heterozygote" spike has a randomised dosage. The review's
practical advice stands in the meantime — *use an explicit 0.5*.

### The measurement that decides it

**Works:** the `vcf_hom` probe round-trips — an input `GT=1/1` produces a truth record with
`GT=1/1` and a merged BAM with essentially no reference support at the locus. A male non-PAR X
deletion writes a haploid GT and removes all of the (one) copy. `af=het` over many seeds gives a
constant event fraction of exactly 0.5 with the *observed* VAF scattering around it, rather than
the event fraction itself scattering.

**Does not work:** GT round-trips but the dosage does not — `GT=1/1` beside residual reference
support that is not explained by CR4's resistant fraction. That would mean the genotype became a
label rather than a description, which is the current defect wearing a new face.

### Rough size

(c): small, and independent — do it now.
(a): medium — a mode flag, an input-GT path, FILTER and hom-ref gating, and truth changes.
(b): large — ploidy declaration, per-contig ploidy, CN in the truth, and it overlaps CR3's
sample-haplotype model. Design those two together.

### What it breaks

(c) changes emitted reads for every `af=het` run — deliberately, and every `af=het` baseline
moves.
(a) and (b) change the truth VCF's GT column for VCF-driven runs, which is what anyone comparing
against it will notice first. All three need release notes; (a) should default to today's
behaviour for a release, with the new mode opt-in, then flip.

---

## CR9, the rest — the validation overhaul

The documentation half of CR9 is **done** in this run (task 8): `README.md` now states, per
check, what it establishes and what it does not, and says that `scripts/validate_pipeline.sh` is
an integration and regression test rather than a measure of caller sensitivity or precision. The
overhaul below is what remains.

### The problem

Several checks are evidence-presence checks wearing the name of event validation.

Measured: CR2's **4.32x** interior amplification **passes** the DUP event-average
`coverage_ratio` at an observed 1.32 against an expected 1.50. CR4's **37.5x** of resistant
low-MAPQ depth **passes** at an observed 0.00, because validate's own `--min-mapq` default of 20
hides it. `check_split_reads` counts SA entries within 500 bp of the partner and reads only an
entry's contig and position (`src/validate.rs:2237-2243`) — never strand, CIGAR, sequence or
allele fraction — and it pools its two required reads **across both breakpoints**, so both may
sit at one end (NF5). `check_ins_reads` counts `I` operations and, above `min_len` 50, soft
clips near POS, and never reads which bases were inserted. Every global check is a fixed library
heuristic with no comparison to the measured donor, so a probe's nonzero validator exit comes
from `insert_size 400+/-0` and `dup_rate no dup flags` rather than from anything being wrong
with the event.

### The governing principle

Four different things are being conflated — the intended genome, the molecules emitted, the
alignment evidence, and the caller's output — and a check should say which one it tests.

### The options

**Option A — a full four-layer validator.** Genome truth (sequence and CN at the locus),
molecular truth (what spike emitted, recorded at generation time), alignment evidence (what the
merged BAM shows), and caller output, each with its own checks, matched by event identity, and
each check declaring which layer it speaks for.

**Option B — strengthen the existing checks in place.** Binned local depth instead of one
event-average; both junctions and their orientations; sequence identity for INS (now possible —
task 4 put the inserted sequence in the truth); require evidence at each breakpoint rather than
pooled; report the resistant fraction from CR4; and compare the globals against the measured
donor rather than fixed thresholds.

**My recommendation: B, and then A only if a formal validation study needs the provenance.**
Every item in B is individually small, immediately useful, and directly closes a measured blind
spot — the 1.32 pass and the 0.00 pass are the two that matter most. A is a re-architecture
whose value is mostly in provenance and reporting, not in catching more defects, and it should
follow the model work (CR2, CR3, CR4) rather than lead it. Validating a model that is still
wrong is effort spent on the wrong layer.

One item of A is worth doing early and cheaply: **record the molecular truth at generation
time** — how many fragments were planted, on which copy, with what sequence. Spike knows all of
it and currently throws it away, and several of B's checks become trivial once it is written
down. Task 5 already took the first step by recording the fraction actually simulated.

### The measurement that decides it

**Works:** re-run the review's four probes. The `dup_variable` probe must now **FAIL** its depth
check (today it passes at 1.32), and `del_lowmap` must **FAIL** (today it passes at 0.00). A
correct control run must still pass everything. That pair — a true positive on the two known
blind spots and no false positive on a good run — is the whole test.

**Does not work:** the two probes fail but so does the control, or so do real correct spike-ins
on HG002 chr20. A validator that cries wolf is worse than one with known blind spots, because
the blind spots are documented and the noise is not. Measure the false-failure rate on at least
twenty correct real-data events before changing any default.

### Rough size

B: medium overall, but it is six independent small items and can land one at a time. The INS
sequence check is now nearly free.
A: large. It needs a provenance record, an identity-matching layer, and a report format.
The molecular-truth record alone: small-to-medium, and it unblocks much of B.

### What it breaks

B makes `spike validate` **fail on runs that pass today** — that is its purpose, and it will
look like a regression to anyone with a green pipeline. Stage it behind a flag for one release,
report both verdicts, then flip the default.
Comparing the globals against the measured donor also changes what `validate` needs as input: it
would have to be given the donor BAM, which today it is not.

---

## Sequencing, if the human wants one list

1. **CR7(c)** — separate `af=het`'s Beta from the event fraction. Small, independent, and every
   heterozygote spike is randomised until it is done.
2. **CR3 fail-closed** — an unreadable sample-variant input becomes an error. Small; the current
   behaviour is a silent wrong answer.
3. **CR6 relabel** — say what fusion mode is. Small; it removes a false claim.
4. **CR2 option B** and **CR4 option B** — gate and census. Both make a measured, currently
   silent error visible.
5. **CR9 option B** — strengthen the checks, starting with binned depth and the resistant
   fraction, so that 4 has something to report into.
6. **CR1's grouped-event composition** — not on this list, because it is a Phase 2 item this run
   only *gated*. It blocks CR6 option A and is the prerequisite for anything multi-event.
7. **CR2 option A**, **CR3 option A**, **CR7(a)/(b)**, **CR9 option A** — the model work, in
   whatever order the validation scope demands. Design CR3's sample haplotypes and CR7's ploidy
   together; they are one object.
