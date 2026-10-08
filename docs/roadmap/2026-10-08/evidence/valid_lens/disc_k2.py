# Quality-only discriminator on the hospital K2 events: SPIKE_ reads vs real reads over the same footprints.
# pre-v2 (inspect/, bdd568a) vs v2 (qtq2/k2v2/). Seen-not-judged look for an idea list.
import os, glob, numpy as np, pysam
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.model_selection import cross_val_predict, GroupKFold
from sklearn.metrics import roc_auc_score
S='/tmp/claude-1066/-home-parlar-ai-dev-spike/23a35c29-9600-4da2-a2c9-618425a1e0d3/scratchpad'
def feats(q, seq, rl):
    q=np.asarray(q,dtype=np.int16); n=len(q); low=q<15
    lr=cur=0
    for v in low:
        cur=cur+1 if v else 0; lr=max(lr,cur)
    return [q.mean(), q.std(), low.sum(), (q==2).sum(), seq.count('N'), q[-1], q[0], q[:10].mean(), q[-10:].mean(),
            np.count_nonzero(np.diff(q)), lr, int(np.argmax(low)) if low.any() else n, low[-10:].sum(), n]
def collect(d):
    X=[];y=[];g=[]
    for ev in sorted(glob.glob(d+'/aa*')):
        bed=open(ev+'/run/events.bed').readline().split('\t')
        c,s,e=bed[0],int(bed[1])-2000,int(bed[2])+2000
        for path,lab in ((ev+'/spiked.bam',1),(ev+'/real.bam',0)):
            for a in pysam.AlignmentFile(path).fetch(c,max(0,s),e):
                if a.is_secondary or a.is_supplementary or a.is_duplicate or a.is_unmapped or a.reference_end is None: continue
                isspk=a.query_name.startswith('SPIKE_')
                if lab==1 and not isspk: continue
                if lab==0 and isspk: continue
                if a.reference_start<s or a.reference_end>e: continue
                q=list(a.query_qualities); sq=a.query_sequence
                if a.is_reverse: q=q[::-1]
                X.append(feats(q,sq,len(q))); y.append(lab); g.append(os.path.basename(ev))
    return np.array(X,float),np.array(y),np.array(g)
for name,d in (('pre-v2 bdd568a',S+'/inspect'),('v2',S+'/qtq2/k2v2')):
    X,y,g=collect(d)
    p=cross_val_predict(HistGradientBoostingClassifier(max_iter=150),X,y,groups=g,cv=GroupKFold(5),method='predict_proba')[:,1]
    print(name,'spike reads',int(y.sum()),'real reads',int((1-y).sum()),'events',len(set(g)),'AUC %.3f'%roc_auc_score(y,p))
