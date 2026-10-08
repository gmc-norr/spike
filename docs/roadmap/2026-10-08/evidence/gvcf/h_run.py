"""Locked checks H2 and H3 of docs/superpowers/plans/2026-10-04-gvcf-hom-alt.md (site rule: the addendum)."""
import csv, filecmp, os, shutil, subprocess, sys
from collections import Counter
import pysam
S = '/tmp/claude-1066/-home-parlar-ai-dev-spike/23a35c29-9600-4da2-a2c9-618425a1e0d3/scratchpad'
D = f'{S}/gvcf/h'
NEW, OLD = f'{S}/target-gh/release/spike', f'{S}/target-master/release/spike'
BAM = os.path.expanduser('~/dev/cnv_validation/from_rv/raredisease_seracare_na12878_resample/D24-14230_Seq25-7600_30x_resample/raredisease_results/alignment/D24-14230_Seq25-7600_30x_sorted_md.bam')
REF = os.path.expanduser('~/dev/spike/data/giab_hg38/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta')
CALLS = '/home/parlar_ai/dev/cnv_validation/from_rv/raredisease_seracare_na12878_resample/D24-14230_Seq25-7600_30x_resample/raredisease_results/call_snv/genome/GM24385_seracare_cancer_snv.vcf.gz'
ALIGNER = "bwa-mem2 mem -M -K 100000000 -t 16 -R '@RG\\tID:sim\\tPL:ILLUMINA\\tSM:sim'"
FILES = ['truth.vcf', 'R1.fq.gz', 'R2.fq.gz', 'replaced_reads.txt', 'fastq_removed_reads.txt']
seq = pysam.FastaFile(REF).fetch('chr20').upper()
called = {rec.start for rec in pysam.VariantFile(CALLS).fetch('chr20')}

rows = list(csv.DictReader(open(f'{S}/gvcf/r.tsv'), delimiter='\t'))
windows = {}
for r in rows:
    w = windows.setdefault(int(r['window']), {'start': int(r['start']), 'het': int(r['het_calls']), 'hom': int(r['hom_calls']), 'missed': []})
    if r['master'] == 'MISSED':
        w['missed'].append((int(r['hom_pos']), r['alt']))
order = sorted(windows)

def event(w):
    p = windows[w]['start'] + 2000
    while p in called:
        p += 1
    ref = seq[p]
    alt = 'ACGT'['ACGT'.index(ref) + 1 - 4] if ref in 'ACGT' else None
    return f'snp:chr20:{p + 1}:{ref}:{alt};af=0.5'

def run(tag, who, ev, gvcf, align=False):
    out = f'{D}/{tag}/{who}'
    shutil.rmtree(out, ignore_errors=True); os.makedirs(out)
    cmd = [NEW if who == 'new' else OLD, '--bam', BAM, '--reference', REF, '--seed', '1', '--threads', '16',
           '--aligner', ALIGNER, '--event', ev, '-o', 'out'] + (['--gvcf', CALLS] if gvcf else []) + (['--align'] if align else [])
    with open(f'{out}/spike.log', 'w') as log:
        rc = subprocess.run(cmd, cwd=out, stdout=log, stderr=subprocess.STDOUT).returncode
    return rc, out

def compare(tag, ev, gvcf):
    rcs = {who: run(tag, who, ev, gvcf)[0] for who in ('new', 'old')}
    same = [f for f in FILES if rcs['new'] == rcs['old'] == 0 and filecmp.cmp(f'{D}/{tag}/new/out/{f}', f'{D}/{tag}/old/out/{f}', shallow=False)]
    ok = rcs['new'] == rcs['old'] == 0 and len(same) == len(FILES)
    print(f'{tag} {ev} gvcf={gvcf}: exit new {rcs["new"]} master {rcs["old"]}; identical {len(same)} of {len(FILES)}: {"PASS" if ok else "FAIL"}', flush=True)
    return ok

def last_error(out):
    lines = [l for l in open(f'{out}/spike.log') if 'Error' in l or 'refus' in l.lower()]
    return lines[-1].strip()[:160] if lines else '(no error line)'

results = []
# H2(a): the first window, no --gvcf.
results.append(compare('h2a', event(order[0]), False))
# H2(b): the first window with a het call, --gvcf.
results.append(compare('h2b', event(next(w for w in order if windows[w]['het'])), True))
# H2(c): the first hom-alt-only window with no miss where master runs.
for w in order:
    x = windows[w]
    if x['het'] or not x['hom'] or x['missed']:
        continue
    rc, out = run('h2c_try', 'old', event(w), True)
    if rc:
        print(f'h2c skip window {w}: master refused: {last_error(out)}')
        continue
    results.append(compare('h2c', event(w), True))
    break
# H3: R's missed sites in order; the first where master runs (with --align).
for w in order:
    for pos1, alt in windows[w]['missed']:
        ev = event(w)
        rc_old, out_old = run('h3', 'old', ev, True, align=True)
        if rc_old:
            print(f'h3 skip chr20:{pos1} (window {w}): master refused: {last_error(out_old)}')
            continue
        rc_new, out_new = run('h3', 'new', ev, True, align=True)
        counts = {}
        for who, out in (('new', out_new), ('old', out_old)):
            c = Counter()
            bam = pysam.AlignmentFile(f'{out}/out/sim.bam')
            for col in bam.pileup('chr20', pos1 - 1, pos1, truncate=True, stepper='nofilter', min_base_quality=0, max_depth=100000):
                for pr in col.pileups:
                    if pr.alignment.query_name.startswith('SPIKE_') and not pr.is_del and not pr.is_refskip:
                        b = pr.alignment.query_sequence[pr.query_position]
                        c['alt' if b == alt else 'ref' if b == seq[pos1 - 1] else 'other'] += 1
            counts[who] = dict(c)
        ok = rc_new == 0 and counts['new'].get('alt', 0) > 0 and set(counts['new']) == {'alt'} and counts['old'].get('alt', 0) == 0
        print(f'h3 chr20:{pos1} {seq[pos1 - 1]}>{alt} (window {w}, event {ev}): SPIKE_ reads new {counts["new"]} master {counts["old"]}: {"PASS" if ok else "FAIL"}', flush=True)
        results.append(ok)
        break
    else:
        continue
    break
print('H2/H3: ' + ('PASS' if all(results) and len(results) == 4 else 'FAIL'))
