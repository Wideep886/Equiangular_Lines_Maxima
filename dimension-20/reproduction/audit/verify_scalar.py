#!/usr/bin/env python3
"""Verify the orbit-energy inequality for the 20 general Niemeier models.
NO full candidate-vector shells, graphs, colorings, global marking vectors,
 or supplied global Weyl-orbit representatives. Only scalar norm profiles,
 local parabolic-dominant tables, and one local histogram maximization.
"""
import sys
from pathlib import Path
from orbit_models import *
from energy_bound import combine,cap
if not __debug__:raise RuntimeError('Run with assertions enabled')
BASE=Path(__file__).resolve().parent

@lru_cache(None)
def histograms(ty,N,a,s):
    f=factor(ty);H=[]
    for M,j in f.byclass[a]:
        if M!=N:continue
        H.append({b:Counter((n,n2,tuple(sorted(es)),m,u==v) for n,n2,u,v,es,m in R)
                  for b,R in f.options(j,s).items()})
    assert H;return H

def rowtable(H):
    return {b:tuple((n,n2,es,m,selfpart,k) for (n,n2,es,m,selfpart),k in R.items()) for b,R in H.items()}

def union_max(H):
    E=defaultdict(Counter)
    for h in H:
        for b,R in h.items():E[b]|=R
    return E

def score(H,C,s):
    # Half the signed-orbit sum, double the contribution of a partner-fixed
    # signed orbit. All summands are nonnegative. Linear in any one histogram.
    tabs=list(map(rowtable,H));twice=0
    for c in C:
        for P in combine([t.get(b,()) for t,b in zip(tabs,c)],4*s):
            es=sum((t[2] for t in P),());m=prod(t[3] for t in P);selfpart=all(t[4] for t in P)
            h=5-sum(e for d,e,m in es);b=cap(es,h,m,selfpart)
            twice+=(2 if selfpart else 1)*b*prod(t[5] for t in P)
    return twice

def check(obj):
    L,C,s=validated_general(obj);patterns=0;envelope_pass=0;local_max=0;upper=0;exact_upper=0;ambiguous=0
    for c in C:
        T=[sorted({(N*s//f.d,N) for N,j in f.byclass[a]}) for f,a in zip(L,c)]
        for N in budgets(T,6*s):
            patterns+=1
            H=[histograms(f.ty,n,a,s) for f,n,a in zip(L,N,c)]
            J=[j for j,Lj in enumerate(H) if len(Lj)>1]
            assert len(J)<=1 # follows also from the minimum glue norm >=4
            ambiguous+=bool(J)
            E=[union_max(h) for h in H]
            b=score(E,C,s)//2
            # Exact maximum over a single LOCAL histogram; for an
            # unambiguous profile E is already the actual local data.
            if J:
                j=J[0];A=[]
                for h in H[j]:
                    EE=E.copy();EE[j]=h;k=score(EE,C,s);assert k%2==0;A.append(k//2)
                q=max(A);assert q<=b
            else:
                assert score(E,C,s)%2==0;q=b
            exact_upper=max(exact_upper,q)
            if b<=90:envelope_pass+=1;upper=max(upper,b)
            else:
                local_max+=1;assert J and q<=90,(obj['idx'],c,N,b,q)
                upper=max(upper,q)
    assert envelope_pass+local_max==patterns and exact_upper<=90 and upper<=90
    print('SCALAR',obj['idx'],'patterns',patterns,'ambiguous',ambiguous,'envelope_pass',envelope_pass,'local_max',local_max,'bound',exact_upper,'envelope_or_refined',upper,flush=True)
    return patterns,exact_upper,envelope_pass,local_max
if __name__=='__main__':
    D=json.loads((BASE/'glue_models.json').read_text());assert [o['idx'] for o in D]==list(range(20))
    # Assert coverage of all 20 distinct rooted systems as specified in TYPES.
    assert len({tuple(tuple(t) for t in o['types']) for o in D})==20
    ids=list(map(int,sys.argv[1:])) if len(sys.argv)>1 else list(range(20))
    assert len(set(ids))==len(ids) and set(ids)<=set(range(20))
    R=[check(D[j]) for j in ids]
    print('SCALAR VERIFIED',len(ids),'types;',sum(r[0] for r in R),'scalar patterns;',sum(r[2] for r in R),'envelope passes;',sum(r[3] for r in R),'single-local refinements; maximum',max(r[1] for r in R),flush=True)
