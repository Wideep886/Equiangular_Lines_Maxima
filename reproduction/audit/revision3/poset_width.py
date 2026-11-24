"""Root-shell bounds from Stanley's Sperner theorem for parabolic quotients.

Usage: python3 poset_width.py PATH/TO/THE/EXTRACTED/joint23 [host names joined with _].
No chain certificate is read. Input module ade.py regenerates the Niemeier
glue code and the norm-four and norm-six dominant weights.
"""
from functools import cache
from collections import Counter
from itertools import product
from math import lcm,prod
from pathlib import Path
import sys
root=Path(sys.argv[1])
sys.path.insert(0,str(root/'thirteen'))
from ade import factors,glue,shapes

HOSTS=['A6_A6_A6_A6','A7_A7_D5_D5','A8_A8_A8','A9_A9_D6','A11_D7_E6','A12_A12','A15_D9','A17_E7','A24','D4_D4_D4_D4_D4_D4','D6_D6_D6_D6','D8_D8_D8','D10_E7_E7','D12_D12','D16_E8','D24','E6_E6_E6_E6','E8_E8_E8']

@cache
def degrees(f,I):
    left=set(I);out=[]
    while left:
        sub={left.pop()};todo=list(sub)
        while todo:
            a=todo.pop()
            nbs={i for i in left if f.C[a][i]==-1};sub.update(nbs);left.difference_update(nbs);todo.extend(nbs)
        n=len(sub);edges={a:{i for i in sub if f.C[a][i]==-1} for a in sub}
        branching=[a for a in sub if len(edges[a])==3]
        if not branching:out.extend(range(2,n+2));continue
        a,=branching;arms=[]
        for b in edges[a]:
            i=a;d=1
            while len(edges[b])>1:
                c,=edges[b]-{i};i,b=b,c;d+=1
            arms.append(d)
        if sorted(arms)[:2]==[1,1]:out.extend([n]+list(range(2,2*n,2)))
        else:
            assert sorted(arms)==[1,2,n-4]
            out.extend({6:(2,5,6,8,9,12),7:(2,6,8,10,12,14,18),8:(2,8,12,14,18,20,24,30)}[n])
    return tuple(sorted(out))

@cache
def quotient(num,den):
    # Poincare polynomial of a parabolic quotient; exact division
    def poly(degs):
        p=[1]
        for n in degs:
            q=[0]*(len(p)+n-1)
            for j,x in enumerate(p):
                for k in range(n):q[j+k]+=x
            p=q
        return p
    a,b=poly(num),poly(den);n=len(a)-len(b)+1;q=[]
    for i in range(n):
        u=a[i];q.append(u)
        for j,v in enumerate(b):a[i+j]-=u*v
    assert not any(a) and min(q)>=1 and q==q[::-1]
    assert max(q)==q[len(q)//2]
    return q[len(q)//2]

@cache
def full_orbit(f,x):
    import numpy as np
    return np.array(f.orbit(x),dtype=np.int64)

@cache
def local(f,w,x,mark,L):
    import numpy as np
    # records for parabolic-dominant members of W x: pairing, degree factors
    orbit=full_orbit(f,x)
    I=tuple(i for i,k in enumerate(w) if not k and all(f.C[j][i]==0 for j in mark))
    if mark: orbit=orbit[np.all(orbit[:,list(mark)]==0,axis=1)]
    if I:orbit=orbit[np.all(orbit[:,I]>=0,axis=1)]
    b=np.asarray(f.Z,dtype=np.int64)@np.array(w)*(L//f.d)
    q=Counter()
    for row,m in zip(orbit,(orbit@b).tolist()):
        K=tuple(i for i in I if not row[i])
        q[(m,degrees(f,I),degrees(f,K))]+=1
    return tuple((m,a,b,k) for (m,a,b),k in q.items())

def bound(fs,w,sh,L,mark):
    tot=0;num_orbits=0
    for x in sh:
        pieces=[local(f,h,g,mt,L) for f,h,g,mt in zip(fs,w,x,mark)]
        for row in product(*pieces):
            if sum(t[0] for t in row)!=3*L: continue
            ds=tuple(sorted(a for t in row for a in t[1]))
            ks=tuple(sorted(a for t in row for a in t[2]))
            k=prod(t[3] for t in row)
            num_orbits+=k;tot+=k*quotient(ds,ks)
    return tot

def go(name):
    fs=factors(name.split('_'));code=glue(fs);xs=shapes(fs,code,4);hs=shapes(fs,code,6);L=lcm(*(f.d for f in fs))
    null=tuple(() for _ in fs);marks_tested=0;unmarked=0;bds=[];weak=[];pairs=0
    for k,w in enumerate(hs):
        b=bound(fs,w,xs,L,null);unmarked=max(unmarked,b)
        if b>=116:
            weak.append(k); caps=[]
            for u,(f,v) in enumerate(zip(fs,w)):
                for i in range(f.n):
                    for j in range(i):
                        if v[i] or v[j] or f.C[i][j]!=-1:continue
                        marks_tested+=1
                        mt=list(null);mt[u]=(j,i)
                        caps.append(bound(fs,w,xs,L,tuple(mt)))
            pairs+=sum(t>=116 for t in caps)
            b=min(b,max(caps,default=0))
        bds.append(b)
    print(name,len(hs),'unmarked max',unmarked,'unmarked unresolved',len(weak),'marks tested',marks_tested,'unresolved marks',pairs,'effective max',max(bds),'unresolved vectors',sum(b>=116 for b in bds),flush=True)

if __name__=='__main__':
    for name in sys.argv[2:] or HOSTS:go(name)
