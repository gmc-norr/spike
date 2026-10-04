#!/usr/bin/env python3
"""Reproduce the 2026-10-04 independent review on disposable synthetic data.

Requires Python 3.9+, samtools on PATH, and target/debug/spike (cargo build).
Uses no real sequencing files. The input/output-alias probe deliberately tests
file destruction using only a generated throwaway FASTQ inside a fresh /tmp dir.
This is a diagnostic snapshot, not a caller benchmark or an assertion that the
current bugs should be preserved. Results and inputs stay in the printed folder.
"""
import gzip, json, pathlib, random, subprocess, shutil, tempfile

ROOT = pathlib.Path(tempfile.mkdtemp(prefix='spike-review-20261004-'))
ROOT.mkdir(exist_ok=True)
SPIKE = str(pathlib.Path(__file__).resolve().parents[1] / 'target/debug/spike')
SAMTOOLS = shutil.which('samtools')
if SAMTOOLS is None:
    raise SystemExit('samtools is required on PATH')
print(f'Review artifacts: {ROOT}', flush=True)
random.seed(142)
REF = ''.join(random.choices('ACGT', k=12000))
REFPATH = ROOT / 'ref.fa'
REFPATH.write_text('>chr1\n' + REF + '\n')
subprocess.run([SAMTOOLS, 'faidx', str(REFPATH)], check=True)

def rc(s): return s.translate(str.maketrans('ACGT', 'TGCA'))[::-1]
def alt(b): return next(x for x in 'ACGT' if x != b)
def fq_write(path, records):
    with gzip.open(path, 'wt') as f:
        for name, seq in records:
            f.write(f'@{name}\n{seq}\n+\n' + '~'*len(seq) + '\n')
def fq_read(path):
    with gzip.open(path, 'rt') as f:
        lines=f.read().splitlines()
    return [(lines[i][1:].split()[0].removesuffix('/1').removesuffix('/2'), lines[i+1]) for i in range(0,len(lines),4)]

def fixture(name, step=25, variants=None, fragment=300, read_lengths=(100,100)):
    variants=variants or {}
    seq=list(REF)
    for pos,base in variants.items(): seq[pos]=base
    seq=''.join(seq)
    sam=ROOT/f'{name}.sam'
    raw=[[],[]]
    records=[]
    for i,start in enumerate(range(500,11000-fragment,step)):
        n=f'pair_{i:06d}'
        l1,l2=read_lengths
        for mate,(pos,flag,length,mpos,tlen) in enumerate([
            (start,99,l1,start+fragment-l2,fragment),
            (start+fragment-l2,147,l2,start,-fragment)]):
            bases=seq[pos:pos+length]
            records.append((pos, f'{n}\t{flag}\tchr1\t{pos+1}\t60\t{length}M\t=\t{mpos+1}\t{tlen}\t{bases}\t'+ '~'*length + '\tRG:Z:rg\n'))
            raw[mate].append((n+f'/{mate+1}', rc(bases) if mate else bases))
    sam.write_text('@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:12000\n@RG\tID:rg\tSM:sample\n'+''.join(row for _,row in sorted(records)))
    bam=ROOT/f'{name}.bam'
    subprocess.run([SAMTOOLS, 'view','-b','-o',str(bam),str(sam)],check=True)
    subprocess.run([SAMTOOLS,'index',str(bam)],check=True)
    for mate in [1,2]: fq_write(ROOT/f'{name}_R{mate}.fq.gz',raw[mate-1])
    return bam

def vcf(name, rows):
    p=ROOT/f'{name}.vcf'
    p.write_text('##fileformat=VCFv4.3\n##contig=<ID=chr1,length=12000>\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n'+''.join(f'chr1\t{pos+1}\t.\t{REF[pos]}\t{a}\t.\tPASS\t.\tGT\t{gt}\n' for pos,a,gt in rows))
    return p

def run(name,bam,event,extra=()):
    out=ROOT/name
    cmd=[SPIKE,'--bam',str(bam),'--reference',str(REFPATH),'--event',event,'--flank','2000','--threads','1','--seed','42','-o',str(out),*map(str,extra)]
    p=subprocess.run(cmd,text=True,capture_output=True)
    (ROOT/f'{name}.log').write_text(p.stdout+p.stderr)
    if p.returncode: raise RuntimeError(f'{name}: {p.stderr}')
    return out

def kmer_counts(out,pos,variant):
    refprobe=REF[pos-12:pos+13]
    altprobe=REF[pos-12:pos]+variant+REF[pos+1:pos+13]
    counts={'ref':0,'alt':0}
    for m in [1,2]:
        for n,seq in fq_read(out/f'R{m}.fq.gz'):
            if not n.startswith('SPIKE_'): continue
            counts['ref']+=refprobe in seq or rc(refprobe) in seq
            counts['alt']+=altprobe in seq or rc(altprobe) in seq
    return counts

results={}
bg=4000
bam=fixture('low',variants={bg:alt(REF[bg])})
gvcf=vcf('hom_only',[(bg,alt(REF[bg]),'1/1')])
control=vcf('hom_and_het',[(bg,alt(REF[bg]),'1/1'),(6000,alt(REF[6000]),'0/1')])
event=f'snp:chr1:5001:{REF[5000]}:{alt(REF[5000])};af=1'
for name,g in [('hom_only',gvcf),('hom_and_het',control)]:
    out=run(name,bam,event,['--gvcf',g])
    results[name]=kmer_counts(out,bg,alt(REF[bg]))

# Request an existing homozygous ALT at 50%: replacement must remove original ALT support too.
existing=fixture('existing',step=5,variants={5000:alt(REF[5000])})
out=run('existing_event',existing,f'snp:chr1:5001:{REF[5000]}:{alt(REF[5000])};af=0.5')
refprobe=REF[4988:5013]
altprobe=REF[4988:5000]+alt(REF[5000])+REF[5001:5013]
counts={'ref':0,'alt':0}
for m in [1,2]:
    for n,seq in fq_read(out/f'R{m}.fq.gz'):
        counts['ref']+=refprobe in seq or rc(refprobe) in seq
        counts['alt']+=altprobe in seq or rc(altprobe) in seq
results['existing_event']=counts

# Asymmetric paired-end cycles are conflated into one shared cycle count.
asym=fixture('asym',step=5,read_lengths=(100,150))
out=run('asym_out',asym,event)
results['asymmetric_mates']={str(m):sorted({len(s) for n,s in fq_read(out/f'R{m}.fq.gz') if n.startswith('SPIKE_')}) for m in [1,2]}

# Generated FASTQ script should not report success on mismatched mate order.
script=ROOT/'hom_only'/'fastq.sh'
for mate,records in [(1,[('a/1','ACGT'),('b/1','TGCA')]),(2,[('b/2','TGCA'),('a/2','ACGT')])]:
    fq_write(ROOT/f'wrong_R{mate}.fq.gz',records)
testdir=ROOT/'route'
testdir.mkdir(exist_ok=True)
shutil.copy(script,testdir/'fastq.sh')
(testdir/'fastq_removed_reads.txt').write_text('')
for mate in [1,2]: fq_write(testdir/f'R{mate}.fq.gz',[(f'SPIKE_ev0001_hap_000000/{mate}','ACGT')])
outs=[testdir/f'out_R{mate}.fq.gz' for mate in [1,2]]
p=subprocess.run(['bash',str(testdir/'fastq.sh'),str(ROOT/'wrong_R1.fq.gz'),str(ROOT/'wrong_R2.fq.gz'),*map(str,outs)],text=True,capture_output=True)
results['mismatched_mates']={'exit':p.returncode,'names':[[n for n,s in fq_read(path)] for path in outs]}

# Output aliasing an input truncates the input; probe only on expendable fixtures.
alias=ROOT/'alias_R1.fq.gz'
shutil.copy(ROOT/'wrong_R1.fq.gz',alias)
p=subprocess.run(['bash',str(testdir/'fastq.sh'),str(alias),str(ROOT/'wrong_R2.fq.gz'),str(alias),str(testdir/'alias_R2.fq.gz')],text=True,capture_output=True)
results['input_output_alias']={'exit':p.returncode,'raw_input_still_exists':alias.exists(),'stderr':p.stderr.strip()}

(ROOT/'results.json').write_text(json.dumps(results,indent=2)+'\n')
print(json.dumps(results,indent=2))


import gzip, json, pathlib, subprocess, shutil
R=ROOT

ST=SAMTOOLS
ref=''.join((R/'ref.fa').read_text().splitlines()[1:])
def alt(b): return next(x for x in 'ACGT' if x!=b)
def fq(path):
    with gzip.open(path,'rt') as f: ls=f.read().splitlines()
    return [(ls[i][1:].split()[0],ls[i+1]) for i in range(0,len(ls),4)]
def rc(s): return s.translate(str.maketrans('ACGT','TGCA'))[::-1]
def run(name,event,extra=()):
    out=R/name
    p=subprocess.run([SPIKE,'--bam',str(R/'low.bam'),'--reference',str(R/'ref.fa'),'--event',event,'--flank','2000','--threads','1','--seed','42','-o',str(out),*extra],text=True,capture_output=True)
    (R/f'{name}.log').write_text(p.stdout+p.stderr)
    return p,out
results={}
# The floor claims to plant an event but counts pairs across the whole footprint.
p,out=run('low_fraction',f'snp:chr1:5001:{ref[5000]}:{alt(ref[5000])};af=0.001')
probe=ref[4988:5000]+alt(ref[5000])+ref[5001:5013]
synth=[(n,s) for m in (1,2) for n,s in fq(out/f'R{m}.fq.gz') if n.startswith('SPIKE_')]
results['low_fraction']={'exit':p.returncode,'synth_pairs':len(synth)//2,'alt_probe_reads':sum(probe in s or rc(probe) in s for n,s in synth),'truth':[l for l in (out/'truth.vcf').read_text().splitlines() if not l.startswith('#')]}

# Values that are not probabilities are accepted for indel_error_rate.
for value in ['NaN','2.0','-0.5']:
    p,out=run('indel_rate_'+value,f'snp:chr1:5001:{ref[5000]}:{alt(ref[5000])}',[f'--indel-error-rate={value}'])
    results['indel_rate_'+value]={'exit':p.returncode,'wrote_truth':(out/'truth.vcf').exists()}

# A correct 50-base deletion represented with CIGAR D, with uniform haplotype sampling.
start,end=5000,5050
hap=ref[:start]+ref[end:]
records=[]
def map_pos(p):return p if p<start else p+end-start
for i,hstart in enumerate(range(500,10700,5)):
    n=f'SPIKE_ev0001_hap_{i:06d}'
    span=map_pos(hstart+299)+1-map_pos(hstart)
    for pos,flag,mate,tlen in [(hstart,99,hstart+200,span),(hstart+200,147,hstart,-span)]:
        if pos<start<pos+100:
            a=start-pos
            cigar=f'{a}M{end-start}D{100-a}M'
        else:cigar='100M'
        row=f'{n}\t{flag}\tchr1\t{map_pos(pos)+1}\t60\t{cigar}\t=\t{map_pos(mate)+1}\t{tlen}\t{hap[pos:pos+100]}\t'+ '~'*100+'\n'
        records.append((map_pos(pos),row))
sam=R/'true_del.sam'
sam.write_text('@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:12000\n'+''.join(v for _,v in sorted(records)))
bam=R/'true_del.bam'
subprocess.run([ST,'view','-b','-o',str(bam),str(sam)],check=True)
subprocess.run([ST,'index',str(bam)],check=True)
truth=R/'true_del.vcf'
truth.write_text('##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'+f'chr1\t5000\tsim_del_1\t{ref[4999]}\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=5050;SIM_VAF=1\n')
p=subprocess.run([SPIKE,'validate','--bam',str(bam),'--reference',str(R/'ref.fa'),'--truth',str(truth),'--json'],text=True,capture_output=True)
(R/'true_del_validate.json').write_text(p.stdout)
(R/'true_del_validate.log').write_text(p.stderr)
results['true_del_validate']={'exit':p.returncode,'stdout':p.stdout}
p=subprocess.run([ST,'depth','-aa','-r','chr1:5001-5050',str(bam)],text=True,capture_output=True,check=True)
results['true_del_base_depth']=sorted({int(l.split()[2]) for l in p.stdout.splitlines()})
(R/'more_results.json').write_text(json.dumps(results,indent=2)+'\n')
print(json.dumps(results,indent=2))


import pathlib, subprocess, json, gzip
R=ROOT; ST=SAMTOOLS

rows=[]
for line in (R/'low.sam').read_text().splitlines():
    if line.startswith('@'): rows.append(line); continue
    f=line.split('\t')
    if int(f[1])&128:
        f[4]='0'
        f.append('XA:Z:chr1,+9001,100M,0;')
    rows.append('\t'.join(f))
sam=R/'discordant_confidence.sam'; sam.write_text('\n'.join(rows)+'\n')
bam=R/'discordant_confidence.bam'
subprocess.run([ST,'view','-b','-o',str(bam),str(sam)],check=True)
subprocess.run([ST,'index',str(bam)],check=True)
results={}
for model in ['clean','origin']:
    out=R/f'confidence_{model}'
    p=subprocess.run([SPIKE,'--bam',str(bam),'--reference',str(R/'ref.fa'),'--event','snp:chr1:5001:T:A;af=1','--flank','2000','--threads','1','--seed','42','--min-mapq','0','--edit-model',model,'-o',str(out)],text=True,capture_output=True)
    (R/f'confidence_{model}.log').write_text(p.stdout+p.stderr)
    relevant=[l for l in p.stderr.splitlines() if any(w in l for w in ['origin depth at','Tiling ','origin removed','simulate_event:'])]
    results[model]={'exit':p.returncode,'log':relevant}
    if p.returncode==0:
        results[model]['removed_pairs']=len((out/'fastq_removed_reads.txt').read_text().splitlines())
        with gzip.open(out/'R1.fq.gz','rt') as f: ls=f.read().splitlines()
        results[model]['synthetic_pairs']=sum(l.startswith('@SPIKE_') for l in ls[::4])
    else: results[model]['stderr']=p.stderr
(R/'origin_results.json').write_text(json.dumps(results,indent=2)+'\n')
print(json.dumps(results,indent=2))
