export const meta = {
  name: 'spike-improvement-roadmap',
  description: 'Map spike, propose improvements through 7 lenses, adversarially verify each, and rank a roadmap',
  phases: [
    { title: 'Map', detail: '6 readers map subsystems, history and known gaps' },
    { title: 'Ideate', detail: '7 lenses propose improvements with evidence and a falsifying measurement' },
    { title: 'Merge', detail: 'dedupe and normalise the ideas' },
    { title: 'Verify', detail: '2 skeptics per idea: reality/history and value/measurement' },
    { title: 'Critic', detail: 'completeness critic finds what is missing; new ideas verified' },
    { title: 'Synthesize', detail: 'ranked roadmap' },
  ],
}

const A = args || {}
const REPO = A.repo
const MEM = A.memory
const DATE = A.date

const GROUND = `
You are analysing the Rust tool **spike** at ${REPO} (read-only). spike plants synthetic variants
(SNV, small indels, DEL, DUP, INV, INS, fusions/translocations, exon-level events, LOH, VCF input) into a real
sample's sequencing reads: it extracts donor read pairs around each event from a real BAM/CRAM, edits the
haplotype, re-synthesises reads with a learned quality/error model, and writes FASTQ + truth VCF + scripts to align
and merge back (or a whole-sample spiked FASTQ). Users benchmark variant callers with it; the owner's end goal is
spiking clinically relevant variants (e.g. LDLR) into hospital demo data run through the nf-core/raredisease pipeline.

Project culture you must respect:
- "Model what really happens": real data is the realism yardstick; claims rest on measurements against real reads.
- Judgment gate: every proposed mechanism needs a falsifying measurement; many past mechanisms were REFUTED.
  The case file ${REPO}/.claude/judgment-gate-cases.md and the owner's notes in ${MEM}/*.md record what was
  measured, refuted, merged and still open. Do not re-propose a refuted mechanism without new evidence.
- Label every number as measured (cite file/command), documented (cite doc), or inferred.

Hard rules: READ-ONLY. Do not edit, create, commit, checkout or build anything inside ${REPO} or any git worktree
(including /home/parlar_ai/quality-speed-run, where an unattended run is working). No cargo build/test.
Cheap read-only commands (grep, git log/show, samtools view/idxstats on data under ${REPO}/data) are fine;
keep anything heavy under ~1 minute of CPU. Use absolute paths.`

const MAP_SCHEMA = {
  type: 'object',
  properties: {
    subsystem: { type: 'string' },
    summary: { type: 'string', description: 'how it works, 150-300 words' },
    algorithms: { type: 'array', items: { type: 'object', properties: {
      name: { type: 'string' }, how: { type: 'string' }, where: { type: 'string', description: 'file:line' } },
      required: ['name', 'how', 'where'] } },
    limitations: { type: 'array', items: { type: 'object', properties: {
      what: { type: 'string' }, evidence: { type: 'string' }, evidence_kind: { type: 'string', enum: ['measured', 'documented', 'inferred'] },
      where: { type: 'string' } }, required: ['what', 'evidence', 'evidence_kind'] } },
    open_issues: { type: 'array', items: { type: 'string' } },
    refuted_or_done: { type: 'array', items: { type: 'string' }, description: 'mechanisms already tried and refuted, or already built, relevant to this area' },
  },
  required: ['subsystem', 'summary', 'algorithms', 'limitations', 'open_issues', 'refuted_or_done'],
}

const READERS = [
  { key: 'events', prompt: `Map the EVENT MODEL and haplotype editing: src/simulate.rs, src/haplotype.rs, src/origin.rs (edit model "origin" vs "clean"), src/exon.rs, src/loh.rs, src/carried.rs, src/vcf_input.rs, src/types.rs. How each variant type is planted (which reads are remade, how VAF/zygosity/ploidy are realised, breakpoints, inserted sequence, look-alike copies, refusals like SIM_RESIST/RF8). What is unrealistic or missing vs real biology (phasing with the sample's own variants, complex/compound events, STRs, MEIs, microhomology, CNV allelic balance, sex chromosomes, mosaicism).` },
  { key: 'donors', prompt: `Map DONOR READ HANDLING: src/extract.rs, src/reference.rs, src/bam_stats.rs, src/read_name.rs, src/census.rs and the relevant parts of src/main.rs. How pairs are selected (MAPQ, proper pairs, duplicates, supplementary), windows, coverage estimation (SIM_DEPTH_FOLD), library characterisation (read length, adapter trimming), CRAM, contig naming. Where selection bias or edge cases could distort the output (low mappability, duplicates, non-proper pairs, chimeras, tiny/huge events, high depth).` },
  { key: 'synth', prompt: `Map READ SYNTHESIS and the QUALITY/ERROR model: src/synth.rs, src/quality.rs, src/stats.rs, docs/analysis/quality-model-v2/, docs/analysis/soft-clips/, README section "Quality profile", plan docs docs/superpowers/plans/*quality*. Fragment lengths, read tiling, base errors, indel errors, adapters, read classes, runs, mates, N, sequencing order. List what real reads show that spike's reads do not (and which of those were measured), e.g. indel errors after homopolymers, GC/coverage bias, PCR duplicates, optical duplicates, chimeric fragments, index hopping, strand bias, cycle effects, platform differences (NovaSeq binned qualities, other instruments).` },
  { key: 'outputs', prompt: `Map OUTPUTS and PIPELINE INTEGRATION: src/fastq.rs, src/truth.rs, src/main.rs (CLI, output scripts align.sh/merge.sh/fastq.sh/README written into the output dir, --into-fastq, JSON output), README sections "Output files", "CLI reference", "What spike refuses", "Pipeline validation examples". How users run it end-to-end (BAM merge vs whole-sample FASTQ; raredisease). Ergonomics, failure modes, performance hot spots (main.rs is 6k lines), determinism (--seed, --threads), safety for clinical data (names, leakage).` },
  { key: 'validate', prompt: `Map VALIDATION and measurement tooling: src/validate.rs (8.6k lines), src/census.rs, src/truth.rs checks, scripts/ (realism_probe.py, real_events*, origin_physics.py, transplant/validation scripts), docs/review/*. What "spike validate" checks, which checks are advisory vs gating, which realism yardsticks exist (transplant HG001<->HG002, real-variant comparisons), and what is still unmeasured (e.g. caller-level equivalence: does a caller call spiked variants like real ones?).` },
  { key: 'history', prompt: `Build the HISTORY and KNOWN-GAPS inventory. Read ${REPO}/.claude/judgment-gate-cases.md fully, every note in ${MEM}/ (MEMORY.md index first), docs/review/REVIEW.md, docs/review/2026-10-04-independent-review.md, docs/review/CLINICAL_SV_*.md, and \`git -C ${REPO} log --oneline | head -300\` subjects (plus 'result:' commit messages, which hold verdicts). Output: (a) every open issue / unfixed finding still listed as open (RF15, PD-27, PD-29, CR3, review finding 6, untested SNVs/indels/dups at the hospital, slices refuted for Manta/CNVnator, etc.), (b) every REFUTED mechanism with one line on why, (c) the owner's stated goals and preferences, (d) what is merged/pushed. Put (a) in limitations/open_issues and (b) in refuted_or_done.` },
]

const IDEA_ITEM = {
  type: 'object',
  properties: {
    title: { type: 'string' },
    category: { type: 'string', enum: ['realism-reads', 'realism-biology', 'algorithm-performance', 'validation-methodology', 'usability-clinical', 'robustness-correctness', 'new-capability'] },
    problem: { type: 'string', description: 'the gap in spike today, concretely' },
    evidence: { type: 'string', description: 'why we believe the gap exists; cite file:line, doc, memory note or command output' },
    evidence_kind: { type: 'string', enum: ['measured', 'documented', 'inferred'] },
    proposal: { type: 'string' },
    simplest_version: { type: 'string', description: 'the dumbest version that could work (Gate A q4)' },
    falsifying_measurement: { type: 'string', description: 'the specific measurement and outcome that would kill the idea (Gate A q2), against real data' },
    cheapest_gate_b: { type: 'string', description: 'an hours-not-days experiment to run before building' },
    impact: { type: 'string', description: 'who benefits and how (caller benchmarking, hospital demo, speed...)' },
    cost: { type: 'string', enum: ['S', 'M', 'L', 'XL'] },
    risks: { type: 'string' },
    history_check: { type: 'string', description: 'was this or something like it tried/refuted/built here? cite' },
    files: { type: 'array', items: { type: 'string' } },
  },
  required: ['title', 'category', 'problem', 'evidence', 'evidence_kind', 'proposal', 'simplest_version', 'falsifying_measurement', 'cheapest_gate_b', 'impact', 'cost', 'risks', 'history_check'],
}
const IDEAS_SCHEMA = { type: 'object', properties: { ideas: { type: 'array', items: IDEA_ITEM } }, required: ['ideas'] }

const LENSES = [
  { key: 'reads', prompt: `LENS: realism of the READS as sequencing physics. Compare spike's synthetic reads with what real Illumina (and other platform) data does: substitution/indel error processes, homopolymer slips, quality binning, cycle/strand/tile effects, GC and coverage bias, fragment-length distribution, adapters, PCR and optical duplicates, chimeric fragments, index hopping, soft-clip causes, mapping artifacts after realignment. Prefer gaps with measured evidence in the notes; propose how to learn each from the input sample itself.` },
  { key: 'biology', prompt: `LENS: realism of the VARIANTS and genome biology. Phasing with the sample's own nearby variants, compound/complex events, STR expansions, MEIs (Alu/L1/SVA with TSDs and polyA), breakpoint microhomology/templated insertions, CNV allelic balance and B-allele frequency, segmental duplications and look-alike copies, low-mappability regions, sex chromosomes/ploidy, mosaic/somatic VAF and tumour purity, mitochondrial heteroplasmy, de novo trios, realistic clinical variant sets (ClinVar).` },
  { key: 'algorithms', prompt: `LENS: ALGORITHMS and PERFORMANCE. Time and memory complexity of the main paths (extraction, tiling, synthesis, quality model, validation), threading and determinism, IO (BAM/CRAM decode, FASTQ gzip), whole-sample FASTQ route scalability to 30-40x WGS, startup costs (quality sample), caching, data structures. Better algorithms (e.g. coverage-preserving tiling, exact read-placement models, streaming). Measure where cheap (e.g. read existing timing notes; do not run heavy jobs).` },
  { key: 'validation', prompt: `LENS: VALIDATION METHODOLOGY. How would we know spike's output is indistinguishable from real data at the level that matters (caller outputs)? Transplant designs, real-twin comparisons, caller-equivalence tests (call rates, genotype quality, SV evidence types), statistical power, regression harnesses/CI that guard realism metrics, discriminator tests (can a classifier tell spiked reads from real?), benchmark-ready outputs (GA4GH truth sets, stratifications).` },
  { key: 'clinical', prompt: `LENS: USABILITY and the CLINICAL WORKFLOW. The owner's goal: spike LDLR and other clinically relevant variants into hospital demo data, full FASTQ, through nf-core/raredisease. What blocks or slows that today: untested variant classes at the hospital, CLI ergonomics, event specs (HGVS input? ClinVar IDs? gene/exon names), reports, provenance, privacy (no patient data leakage, read names), docs, failure messages, reproducibility.` },
  { key: 'robustness', prompt: `LENS: ROBUSTNESS and CORRECTNESS. Edge cases and silent-failure risks: contig naming (chr vs no chr), alt/decoy/HLA contigs, chrM, sex chromosomes, CRAM reference mismatch, events near contig ends/gaps/Ns, overlapping events, very small/large events, high-depth regions, empty pools, multi-sample/multi-read-group BAMs, long reads input, odd read lengths. Code health: very large files (main.rs 6k, validate.rs 8.6k lines), test gaps, property/fuzz tests, error handling. Cite concrete file:line.` },
  { key: 'landscape', prompt: `LENS: the LANDSCAPE and NEW CAPABILITIES. Compare spike with BAMSurgeon, VarSim, NEAT, ART, InSilicoSeq, SURVIVOR, Sim-it, VISOR, Mason, simuG, pIRS and others (load WebSearch with ToolSearch "select:WebSearch" if useful). What do they do that spike does not, and what does spike do uniquely? Propose new capabilities that fit spike's design (e.g. long-read support, somatic/tumour mode, trio de novo, RNA-seq/fusions, targeted panels/exome capture bias, UMI data, methylation) only where the owner's goals make them valuable.` },
]

const VERDICT_SCHEMA = {
  type: 'object',
  properties: {
    verdict: { type: 'string', enum: ['keep', 'weaken', 'kill'] },
    gap_is_real: { type: 'boolean' },
    already_done_or_refuted: { type: 'string', description: 'cite if this exists or was refuted; empty if not' },
    corrected_claims: { type: 'string', description: 'any claim in the idea that is wrong, with the correct fact and its source' },
    value: { type: 'integer', minimum: 1, maximum: 5 },
    cost: { type: 'integer', minimum: 1, maximum: 5, description: '1 cheap ... 5 very expensive' },
    measurement_is_sound: { type: 'boolean' },
    better_measurement: { type: 'string' },
    notes: { type: 'string' },
  },
  required: ['verdict', 'gap_is_real', 'already_done_or_refuted', 'corrected_claims', 'value', 'cost', 'measurement_is_sound', 'notes'],
}

const VERIFY_LENSES = [
  { key: 'reality', prompt: `You are a skeptic checking REALITY and HISTORY. Try to refute the idea: (1) does the gap really exist in spike's current code at ${REPO} (master; read the cited files and lines yourself)? (2) Is it already implemented, partly implemented, or already REFUTED here (grep ${REPO}/.claude/judgment-gate-cases.md, ${MEM}/*.md, git log, docs/superpowers/plans)? (3) Are its evidence and numbers correct? Default to verdict=kill if the gap is not real or the mechanism was refuted without new evidence; weaken if partly true.` },
  { key: 'value', prompt: `You are a skeptic checking VALUE and MEASUREMENT. Try to refute the idea on worth: (1) would it change an outcome that matters (a caller's calls on spiked vs real variants, the hospital LDLR/raredisease demo, runtime that blocks users)? (2) Is the falsifying measurement well-posed, against real data, with a pass/fail rule that a broken version would fail? Propose a better one if not. (3) Is the cost estimate honest; is there a simpler version? Default to weaken if the value is speculative.` },
]

function ideaText(i) {
  return `TITLE: ${i.title}\nCATEGORY: ${i.category}\nPROBLEM: ${i.problem}\nEVIDENCE (${i.evidence_kind}): ${i.evidence}\nPROPOSAL: ${i.proposal}\nSIMPLEST: ${i.simplest_version}\nFALSIFYING MEASUREMENT: ${i.falsifying_measurement}\nGATE B: ${i.cheapest_gate_b}\nIMPACT: ${i.impact}\nCOST: ${i.cost}\nRISKS: ${i.risks}\nHISTORY: ${i.history_check}\nFILES: ${(i.files || []).join(', ')}`
}

// ---------------- Map ----------------
phase('Map')
const maps = (await parallel(READERS.map(r => () =>
  agent(`${GROUND}\n\nTASK: ${r.prompt}\nBe concrete and cite file:line. Return the structured map.`,
    { label: `map:${r.key}`, phase: 'Map', schema: MAP_SCHEMA })))).filter(Boolean)
log(`mapped ${maps.length}/${READERS.length} areas`)
const MAPTEXT = maps.map(m => `### ${m.subsystem}\n${m.summary}\nALGORITHMS:\n${m.algorithms.map(a => `- ${a.name}: ${a.how} (${a.where})`).join('\n')}\nLIMITATIONS:\n${m.limitations.map(l => `- ${l.what} [${l.evidence_kind}: ${l.evidence}] ${l.where || ''}`).join('\n')}\nOPEN: ${m.open_issues.join(' | ')}\nREFUTED/DONE: ${m.refuted_or_done.join(' | ')}`).join('\n\n')

// ---------------- Ideate ----------------
phase('Ideate')
const lensIdeas = (await parallel(LENSES.map(l => () =>
  agent(`${GROUND}\n\nHere is a map of spike written by readers of its code, docs and history:\n\n${MAPTEXT}\n\n${l.prompt}\n\nThink deeply. Propose 6-10 improvements through this lens, most valuable first. Check the code yourself before asserting a gap. For each, fill every field; the falsifying measurement must be against real data and must fail for a broken version. Skip anything the map lists as refuted unless you have new evidence (say what).`,
    { label: `ideas:${l.key}`, phase: 'Ideate', schema: IDEAS_SCHEMA }).then(r => r ? r.ideas.map(i => ({ ...i, lens: l.key })) : [])))).filter(Boolean).flat()
log(`${lensIdeas.length} raw ideas from ${LENSES.length} lenses`)

// ---------------- Merge ----------------
phase('Merge')
const MERGED_SCHEMA = {
  type: 'object',
  properties: {
    ideas: { type: 'array', items: { ...IDEA_ITEM, properties: { ...IDEA_ITEM.properties, id: { type: 'string' }, lenses: { type: 'array', items: { type: 'string' } } }, required: [...IDEA_ITEM.required, 'id', 'lenses'] } },
    merged_notes: { type: 'string', description: 'which raw ideas were merged into which, and any dropped as out of scope, with reason' },
  },
  required: ['ideas', 'merged_notes'],
}
const raw = lensIdeas.map((i, k) => `[R${k + 1} lens=${i.lens}]\n${ideaText(i)}`).join('\n\n')
const merged = await agent(`${GROUND}\n\nBelow are ${lensIdeas.length} raw improvement ideas from 7 lenses. Merge duplicates and near-duplicates into one idea each (keep the strongest evidence and the sharpest falsifying measurement; list all source lenses). Do not drop distinct ideas; only drop an idea if it is plainly out of scope for spike, and say so in merged_notes. Give each merged idea an id M1, M2, ... Keep every field.\n\n${raw}`,
  { label: 'merge', phase: 'Merge', schema: MERGED_SCHEMA })
const ideas = merged ? merged.ideas : lensIdeas.map((i, k) => ({ ...i, id: `R${k + 1}`, lenses: [i.lens] }))
log(`${ideas.length} merged ideas`)

// ---------------- Verify ----------------
const verifyIdea = (idea, phaseName) => parallel(VERIFY_LENSES.map(v => () =>
  agent(`${GROUND}\n\n${v.prompt}\n\nThe idea:\n${ideaText(idea)}\n\nReturn your verdict.`,
    { label: `${v.key}:${idea.id}`, phase: phaseName, schema: VERDICT_SCHEMA })))
  .then(vs => ({ idea, verdicts: vs.filter(Boolean).map((x, k) => ({ ...x, lens: VERIFY_LENSES[k] ? VERIFY_LENSES[k].key : '?' })) }))

phase('Verify')
const judged = (await parallel(ideas.map(i => () => verifyIdea(i, 'Verify')))).filter(Boolean)

const score = j => {
  const vs = j.verdicts
  const kills = vs.filter(v => v.verdict === 'kill').length
  const val = vs.length ? vs.reduce((s, v) => s + v.value, 0) / vs.length : 0
  const cost = vs.length ? vs.reduce((s, v) => s + v.cost, 0) / vs.length : 0
  return { kills, val, cost }
}
const survivors = judged.filter(j => score(j).kills < 2 && j.verdicts.length > 0)
const killed = judged.filter(j => !survivors.includes(j))
log(`${survivors.length} survive, ${killed.length} killed`)

// ---------------- Critic ----------------
phase('Critic')
const summaryLine = j => { const s = score(j); return `${j.idea.id} [${j.idea.category}] ${j.idea.title} (value ${s.val.toFixed(1)}, cost ${s.cost.toFixed(1)}, kills ${s.kills})` }
const critic = await agent(`${GROUND}\n\nA panel produced and verified these improvement ideas for spike:\nSURVIVING:\n${survivors.map(summaryLine).join('\n')}\nKILLED:\n${killed.map(summaryLine).join('\n')}\n\nMap of spike:\n${MAPTEXT}\n\nYou are the COMPLETENESS CRITIC. What important opportunities are missing entirely - a lens not taken, a subsystem not examined, an open issue in the history not addressed, a cheap high-value win overlooked, a risk to the owner's hospital goal nobody named? Check the code and notes yourself. Propose up to 8 missing ideas with every field filled (ids C1, C2, ...). Do not repeat the ideas listed.`,
  { label: 'critic', phase: 'Critic', schema: IDEAS_SCHEMA })
const extra = critic ? critic.ideas.map((i, k) => ({ ...i, id: `C${k + 1}`, lenses: ['critic'] })) : []
const judgedExtra = (await parallel(extra.map(i => () => verifyIdea(i, 'Critic')))).filter(Boolean)
const survivorsExtra = judgedExtra.filter(j => score(j).kills < 2 && j.verdicts.length > 0)
const killedExtra = judgedExtra.filter(j => !survivorsExtra.includes(j))
log(`critic added ${extra.length}; ${survivorsExtra.length} survive`)

// ---------------- Synthesize ----------------
phase('Synthesize')
const allSurv = [...survivors, ...survivorsExtra]
const allKilled = [...killed, ...killedExtra]
const dossier = j => {
  const s = score(j)
  return `${ideaText(j.idea)}\nID: ${j.idea.id}  LENSES: ${(j.idea.lenses || []).join(',')}\nMEAN VALUE ${s.val.toFixed(1)} / MEAN COST ${s.cost.toFixed(1)} / KILL VOTES ${s.kills}\n` +
    j.verdicts.map(v => `  - ${v.lens}: ${v.verdict}; real=${v.gap_is_real}; done/refuted: ${v.already_done_or_refuted || '-'}; corrections: ${v.corrected_claims || '-'}; measurement sound=${v.measurement_is_sound}; better: ${v.better_measurement || '-'}; notes: ${v.notes}`).join('\n')
}
const report = await agent(`${GROUND}\n\nWrite the final ROADMAP for improving spike, dated ${DATE}, as GitHub-flavoured markdown. Inputs are verified ideas with skeptic verdicts. Apply every skeptic correction (never repeat a claim a skeptic corrected). Structure:
1. "Top 10" table: rank, idea, why (one line, with its evidence kind), falsifying measurement (short), cost, value. Rank by value/cost and by how directly it serves real-data realism and the hospital LDLR/raredisease goal.
2. "Quick wins" (cost S, value >= 3).
3. "Big bets" (value 5, cost L/XL), each with its Gate B experiment.
4. Sections by category with every surviving idea: problem, evidence (kind + source), proposal, simplest version, falsifying measurement with pass/fail, Gate B, cost, risks, skeptic notes.
5. "Killed or already done": one line each with the reason (so nobody re-proposes them).
6. "Open questions for the owner": decisions only the owner can make.
Be precise; keep each idea tight. Do not invent numbers: every number must come from the inputs, labelled as in them.

SURVIVING IDEAS:\n\n${allSurv.map(dossier).join('\n\n')}\n\nKILLED IDEAS:\n\n${allKilled.map(dossier).join('\n\n')}`,
  { label: 'roadmap', phase: 'Synthesize' })

return {
  report,
  counts: { raw: lensIdeas.length, merged: ideas.length, survived: survivors.length, killed: killed.length, critic: extra.length, critic_survived: survivorsExtra.length },
  merged_notes: merged ? merged.merged_notes : '',
  survivors: allSurv.map(j => ({ id: j.idea.id, title: j.idea.title, category: j.idea.category, ...score(j) })),
  killed: allKilled.map(j => ({ id: j.idea.id, title: j.idea.title, reasons: j.verdicts.map(v => `${v.lens}:${v.verdict}:${v.already_done_or_refuted || v.notes}`.slice(0, 300)) })),
}
