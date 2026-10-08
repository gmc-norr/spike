# Quick "can you tell" gate: quality-string-only features, real vs spike (K2b sets), name-free.
import sys, numpy as np
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import cross_val_predict, GroupKFold
from sklearn.metrics import roc_auc_score
D='/home/parlar_ai/quality-speed-run/evidence/k2b/'
SYM=None
def feats(qs, seq):
    q=np.frombuffer(qs,dtype=np.uint8)-33
    n=len(q)
    low=q<15
    tr=np.count_nonzero(np.diff(q))
    # longest low run
    lr=0;cur=0
    for v in low:
        cur=cur+1 if v else 0
        lr=max(lr,cur)
    first_low = int(np.argmax(low)) if low.any() else n
    lastlow10 = low[-10:].sum()
    return [q.mean(), q.std(), low.sum(), (q==2).sum(), seq.count(b'N'), q[-1], q[0], q[:10].mean(), q[-10:].mean(), tr, lr, first_low, lastlow10,
            ((q==2)&(np.frombuffer(seq,dtype=np.uint8)!=78)).sum()]
def load(path, maxn):
    X=[]
    with open(path,'rb') as f:
        while True:
            h=f.readline()
            if not h or len(X)>=maxn: break
            s=f.readline().rstrip(); f.readline(); q=f.readline().rstrip()
            X.append(feats(q,s))
    return np.array(X,float)
N=int(sys.argv[1]) if len(sys.argv)>1 else 20000
names=['meanQ','sdQ','nlow','nQ2','nN','lastQ','firstQ','first10','last10','transitions','lowrun','firstlow','lastlow10','Q2notN']
for spk in ['spike21']:
    for mate in ['R1','R2']:
        r=load(D+f'real_{mate}.fq',N); s=load(D+f'{spk}_{mate}.fq',N)
        X=np.vstack([r,s]); y=np.r_[np.zeros(len(r)),np.ones(len(s))]
        g=np.r_[np.arange(len(r)),np.arange(len(s))]  # pair index as group: same template in both sets
        cv=GroupKFold(5)
        p=cross_val_predict(HistGradientBoostingClassifier(max_iter=200),X,y,groups=g,cv=cv,method='predict_proba')[:,1]
        print(spk,mate,'HGB AUC all features %.4f'%roc_auc_score(y,p))
        keep=[i for i,n in enumerate(names) if n not in ('nQ2','Q2notN','nN')]
        p2=cross_val_predict(HistGradientBoostingClassifier(max_iter=200),X[:,keep],y,groups=g,cv=cv,method='predict_proba')[:,1]
        print(spk,mate,'HGB AUC without Q2/N features %.4f'%roc_auc_score(y,p2))
        for i,nm in enumerate(names):
            a=roc_auc_score(y,X[:,i]); print('   %-12s single-feature AUC %.4f  real mean %.3f spike mean %.3f'%(nm,max(a,1-a),r[:,i].mean(),s[:,i].mean()))
