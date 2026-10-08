import statistics as st
rows=[l.rstrip('\n').split('\t') for l in open('../transplant2/full/evidence.tsv')]
h=rows[0]; rows=[dict(zip(h,r)) for r in rows[1:]]
H={'forward':{tuple(l.split('\t')[:4]) for l in open('H1.calls.tsv')},'reverse':{tuple(l.split('\t')[:4]) for l in open('H2.calls.tsv')}}
for s in ('forward','reverse'):
  for g in ('INS20-49','DEL20-49','DUP50-299','INS1-4'):
    for which in ('real','fake_normal','fake_B1'):
      called=[];missed=[]
      for r in rows:
        if r['set']!=s or r['group']!=g or r['which']!=which: continue
        k=(r['chrom'],r['pos'],r['ref'],r['alt'])
        m='J' if g.startswith('DUP') else 'E'
        try: v=float(r[m])
        except: continue
        (called if k in H[s] else missed).append(v)
      if called or missed:
        print(s,g,which,'metric',m,'donor-called n',len(called),'med',round(st.median(called),3) if called else None,'| donor-missed n',len(missed),'med',round(st.median(missed),3) if missed else None, 'missed<0.1:',sum(1 for x in missed if x<0.1))
