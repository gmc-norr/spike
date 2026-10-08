"""Mutation run for the carried-allele plan; scores only when the unmutated tests are green."""
import os, re, subprocess, sys
WT = '/tmp/claude-1066/-home-parlar-ai-dev-spike/23a35c29-9600-4da2-a2c9-618425a1e0d3/scratchpad/gh'
TARGET = '/tmp/claude-1066/-home-parlar-ai-dev-spike/23a35c29-9600-4da2-a2c9-618425a1e0d3/scratchpad/target-gh'
C, M = 'src/carried.rs', 'src/main.rs'
L = 'src/loh.rs'
MUTANTS = [
    ('1 merge not called', L, 'call_snps(&counts, region_start, ref_seq).with_hom_alt_from(gvcf_hom_alt)', 'call_snps(&counts, region_start, ref_seq)'),
    ('2 pileup het kept at a gVCF hom-alt', L, '        self.het.retain(|s| !hom_alt.contains_key(&s.pos));\n', ''),
    ('3 pileup allele kept over the gVCF', L, '        self.hom_alt.extend(hom_alt);\n', '        for (p, a) in hom_alt { self.hom_alt.entry(p).or_insert(a); }\n'),
    ('4 pileup other hets dropped', L, '        self.het.retain(|s| !hom_alt.contains_key(&s.pos));\n', '        self.het.clear();\n'),
    ('5 pileup other hom-alts dropped', L, '        self.hom_alt.extend(hom_alt);\n', '        self.hom_alt = hom_alt;\n'),
]
MUTANTS_OLD = [
    ('1 no trimming', C, 'let prefix = ref_allele.iter().zip(alt_allele).take_while(|(r, a)| r == a).count();', 'let prefix = 0;'),
    ('1b no suffix trimming', C, 'let suffix = r.iter().rev().zip(a.iter().rev()).take_while(|(x, y)| x == y).count();', 'let suffix = 0;'),
    ('2 a mismatch past the site counted', C, 'for p in ref_pos.max(start)..(ref_pos + len as u64).min(end) {', 'for p in ref_pos.max(start)..(ref_pos + len as u64).min(end + 1) {'),
    ('3 deletions ignored', C, 'if ref_pos < end.max(start + 1) && ref_pos + len as u64 > start {', 'if false {'),
    ('4 insertions ignored', C, 'if start <= ref_pos && ref_pos <= end {', 'if false {'),
    ('5 insertion boundary s < b < e', C, 'if start <= ref_pos && ref_pos <= end {', 'if start < ref_pos && ref_pos < end {'),
    ('6 spanning not required', C, 'match (align_start < start && ref_pos > end, other) {', 'match (true, other) {'),
    ('7 share 0.2 -> 0.5', C, 'pub const MIN_SHARE: f64 = 0.2;', 'pub const MIN_SHARE: f64 = 0.5;'),
    ('8 floor 10 -> 1', C, 'pub const MIN_READS: u32 = 10;', 'pub const MIN_READS: u32 = 1;'),
    ('9 duplicates counted', C, '        || flags.is_duplicate()\n', ''),
    ('10 MAPQ not filtered', C, '&& mapq.map_or(0, u8::from) >= min_mapq', '&& mapq.map_or(0, u8::from) >= 0'),
    ('11 only the first refused event reported', M, 'refused.join("; ")', 'refused[0].clone()'),
    ('12 an N base counted', C, 'if matches!(base, Some(b) if b != b\'N\' && Some(b) != want) {', 'if matches!(base, Some(b) if Some(b) != want) {'),
    ('13 the check not called', M, '    refuse_carried_alleles(&events, &config.bam_path, &config.ref_path, config.min_mapq, &shared_ref)?;\n', ''),
]
def run():
    env = dict(os.environ, CARGO_TARGET_DIR=TARGET)
    p = subprocess.run(['cargo', 'test', '--locked', '--offline', '--bin', 'spike', '--', 'loh::', 'gvcf'],
                       cwd=WT, env=env, capture_output=True, text=True)
    out = p.stdout + p.stderr
    if 'error[' in out or 'could not compile' in out:
        return None, out[-1500:]
    return sorted(set(re.findall(r'^test (\S+) \.\.\. FAILED', out, re.M))), (re.findall(r'test result: .*', out) or [''])[-1]
orig = {f: open(os.path.join(WT, f)).read() for f in (C, M, L)}
try:
    red, summary = run()
    print(f'unmutated: red={red} {summary}', flush=True)
    if red is None or red:
        sys.exit('unmutated tests are not green; not scoring')
    caught = 0
    for label, f, old, new in MUTANTS:
        n = orig[f].count(old)
        if n != 1:
            print(f'{label}: anchor found {n} times; NOT RUN', flush=True); continue
        open(os.path.join(WT, f), 'w').write(orig[f].replace(old, new))
        red, summary = run()
        open(os.path.join(WT, f), 'w').write(orig[f])
        if red is None:
            print(f'{label}: COMPILE ERROR {summary[-300:]}', flush=True); continue
        caught += bool(red)
        print(f'{label}: {"CAUGHT" if red else "SURVIVED"} by {red}', flush=True)
    print(f'caught {caught} of {len(MUTANTS)}')
finally:
    for f, text in orig.items():
        open(os.path.join(WT, f), 'w').write(text)
