import sys
for f in sys.argv[1:]:
    n=0; s=[0]*12; low=[0]*12; s5=[0]*11; low5=[0]*11
    with open(f) as fh:
        for i,line in enumerate(fh):
            if i%4!=3: continue
            q=line.rstrip('\n'); L=len(q); n+=1
            for k in range(1,11):
                v=ord(q[L-k])-33; s[k]+=v; low[k]+= v<15
            for c in range(1,11):
                v=ord(q[c-1])-33; s5[c]+=v; low5[c]+= v<15
    print(f.split('/')[-1], n, '3p mean Q cl1..10:', ' '.join('%.2f'%(s[k]/n) for k in range(1,11)))
    print('   3p low%% cl1..10:', ' '.join('%.2f'%(100*low[k]/n) for k in range(1,11)))
    print('   5p mean Q c1..10:', ' '.join('%.2f'%(s5[c]/n) for c in range(1,11)), ' low%:', ' '.join('%.2f'%(100*low5[c]/n) for c in range(1,11)))
