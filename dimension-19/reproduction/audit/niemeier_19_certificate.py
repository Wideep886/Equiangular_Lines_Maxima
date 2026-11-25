#!/usr/bin/env python3
"""Exact check of the marked A3^8 and A2^12 Niemeier sections.

Requires Python 3.10+, NumPy, NetworkX. No floating-point rank or LP solver.
Run: python niemeier_19_certificate.py

The seven A2 cases flagged ANALYTIC are exactly those with a marked
component of norm 8/3 and five other components of norm 2/3. They are
excluded by the shortened ternary Golay code / K6 saturation proof,
not by this program alone. All other representatives have bound <=72
for a hypothetical selected set of 73 axes.

The conditional support deletion uses the published N(18)<=59.
The analytic K6 argument uses N(17)=48.
"""
from collections import Counter, defaultdict
from functools import lru_cache
from itertools import product, combinations
from math import isqrt
import numpy as np
import networkx as nx

H4 = ((3,3,2,3,1,0,0,0),(3,0,3,2,3,1,0,0),
      (3,0,0,3,2,3,1,0),(3,0,0,0,3,2,3,1))
P3 = ((0,1,1,1,1,1),(1,0,1,1,2,2),(1,1,0,2,1,2),
      (1,1,2,0,2,1),(1,2,1,2,0,1),(1,2,2,1,1,0))
# Each (p,s) acts by (g(c))_i=s_i*c_{p_i}. All are verified below.
G4 = (
 ((2,3,1,7,0,6,5,4),(1,-1,-1,1,-1,1,1,-1)),
 ((7,2,3,4,1,6,0,5),(-1,1,1,-1,1,-1,-1,1)),
)
G3 = (
 (tuple(range(12)),(-1,)*12),
 ((0,2,4,1,5,3,6,8,10,7,11,9),(1,)*12),
 ((0,2,1,4,3,5,6,8,7,10,9,11),(1,)*12),
 ((0,1,2,3,11,10,7,6,9,8,5,4),(1,1,1,1,-1,-1,1,1,1,1,-1,-1)),
)

def dot(a,b):
    return sum(x*y for x,y in zip(a,b))

def norm(a):
    return dot(a,a)

def glue(m):
    if m==4:
        C={c for c in product(range(4),repeat=8)
           if all(dot(h,c)%4==0 for h in H4)}
    else:
        C={u+tuple(sum(u[i]*P3[i][j] for i in range(6))%3
                   for j in range(6)) for u in product(range(3),repeat=6)}
    n=24//(m-1)
    assert len(C)==m**(n//2)
    # The minimum norm of discriminant class c in A_{m-1} is c(m-c)/m.
    # Even norms for every coset imply integrality by polarization.
    sums=Counter(sum(c*(m-c) for c in v) for v in C)
    assert sums[0]==1
    assert all(t%(2*m)==0 for t in sums)
    assert all(t==0 or t>=4*m for t in sums)
    # Index |C| in the extension gives determinant 1. No new roots.
    for p,s in (G4 if m==4 else G3):
        assert sorted(p)==list(range(n)) and set(s)<=set((-1,1))
        assert all(tuple(s[i]*c[p[i]]%m for i in range(n)) in C for c in C)
    return sorted(C)

@lru_cache(None)
def local(m,cap,dominant=False):
    """All scaled local dual vectors v: sum v_i=0; v_i congruent mod m."""
    out=[]; R=isqrt(cap)
    for a in product(range(-R,R+1),repeat=m-1):
        v=a+(-sum(a),); N=norm(v)
        if N>cap or any((x-v[0])%m for x in v):continue
        if dominant and v!=tuple(sorted(v,reverse=True)):continue
        out.append(((-v[0])%m,N,v))
    return out

def profiles(m,C):
    """All Weyl-normalized (w,marked component,equal coordinate value)."""
    n=len(C[0]); cap=6*m*m; D=defaultdict(list); out=set()
    for c,N,v in local(m,cap,True): D[c].append((N,v))
    for c in C:
        minimum=[min(N for N,v in D[a]) for a in c]
        rem=[sum(minimum[i:]) for i in range(n+1)]
        if rem[0]>cap:continue
        def visit(i,N,w):
            if i==n:
                if N==cap:
                    for j,v in enumerate(w):
                        for a,count in Counter(v).items():
                            if count>=2:out.add((w,j,a))
                return
            for t,v in D[c[i]]:
                if N+t+rem[i+1]<=cap:visit(i+1,N+t,w+(v,))
        visit(0,0,())
    return out

def transform(record,g):
    w,j,a=record; p,s=g; k=p.index(j)
    return (tuple(tuple(sorted((s[i]*x for x in w[p[i]]),reverse=True))
                  for i in range(len(w))),k,a*s[k])

def representatives(m,records):
    remaining=set(records); ans=[]; G=G4 if m==4 else G3
    while remaining:
        p=min(remaining); remaining.remove(p); queue=[p]
        for x in queue:
            for g in G:
                y=transform(x,g)
                assert y in records
                if y in remaining:remaining.remove(y);queue.append(y)
        ans.append((p,len(queue)))
    assert sum(k for p,k in ans)==len(records)
    return ans

@lru_cache(None)
def local_options(m,w,mark):
    D=defaultdict(list)
    for c,N,x in local(m,4*m*m):
        if mark and x[mark[0]]!=x[mark[1]]:continue
        D[c].append((N,dot(x,w),x))
    return D

def shell(m,C,w,mark):
    """Scaled x, norm 4*m^2, x.w=3*m^2, marked entries equal."""
    n=len(w); b,ij=mark; cap=4*m*m; pair=3*m*m
    tabs=[local_options(m,v,ij if i==b else ()) for i,v in enumerate(w)]
    wr=[sum(norm(v) for v in w[i:]) for i in range(n+1)]
    ans=[]
    for c in C:
        blocks=[tabs[i][a] for i,a in enumerate(c)]
        if any(not B for B in blocks):continue
        minimum=[min(N for N,d,v in B) for B in blocks]
        rem=[sum(minimum[i:]) for i in range(n+1)]
        if rem[0]>cap:continue
        def visit(i,N,d,x):
            if i==n:
                if N==cap and d==pair:ans.append(x)
                return
            for t,h,v in blocks[i]:
                if N+t+rem[i+1]>cap:continue
                if (pair-d-h)**2>(cap-N-t)*wr[i+1]:continue
                visit(i+1,N+t,d+h,x+v)
        visit(0,0,0,())
    assert len(ans)==len(set(ans))
    W=sum(w,()); S=set(ans)
    assert all(tuple(a-b for a,b in zip(W,x)) in S for x in S)
    return np.array(ans,dtype=np.int64).reshape(-1,m*n)

def compatibility(m,w,mark,X):
    """One representative per antipodal pair and two incompatibility tests."""
    b,ij=mark; Y=2*X-np.array(w,dtype=np.int64).flatten()
    idx=np.argmax(Y!=0,axis=1)
    Y=Y[Y[np.arange(len(Y)),idx]>0]
    assert len(Y)*2==len(X)
    assert np.all(np.sum(Y*Y,axis=1)==10*m*m)
    free=[]
    for i,wi in enumerate(w):
        for a in sorted(set(wi)):
            block=[i*m+j for j in range(m) if wi[j]==a and
                   not (i==b and j in ij)]
            if block:free.append(block)
    fixed=[i for i in range(Y.shape[1]) if i not in sum(free,[])]
    A0=(np.abs(Y@Y.T)==2*m*m); A=A0.copy()
    def cone(D):
        ok=np.ones(len(D),dtype=bool);pos=ok.copy();neg=ok.copy()
        if fixed:ok&=np.all(D[:,fixed]==0,axis=1)
        for g in free:
            T=np.cumsum(D[:,g],axis=1)
            ok&=(T[:,-1]==0);pos&=np.all(T>=0,axis=1);neg&=np.all(T<=0,axis=1)
        return ok&(pos|neg)
    for i in range(len(Y)):
        A[i]&=~(cone(Y[i]+Y)|cone(Y[i]-Y))
    assert np.array_equal(A,A.T)
    groups=[np.array(sorted(g),dtype=int) for g in
            nx.connected_components(nx.from_numpy_array(~A0))]
    return Y,A,groups

def omega(A):
    """Exact maximum clique, using a greedy coloring bound in each recursion."""
    n=len(A)
    if n==0:return 0
    masks=[sum(1<<j for j in np.flatnonzero(row).tolist()) for row in A]
    best=0
    def expand(P,s):
        nonlocal best
        if not P:best=max(best,s);return
        order=[];colors=[];R=P;k=0
        # Partition current vertices into independent color classes.
        while R:
            k+=1;Q=R
            while Q:
                bit=Q&-Q;v=bit.bit_length()-1
                order.append(v);colors.append(k);R^=bit;Q&=~(bit|masks[v])
        # Reverse color order; remaining prefix has at most c colors.
        for v,c in zip(reversed(order),reversed(colors)):
            if s+c<=best:return
            if (P>>v)&1:
                expand(P&masks[v],s+1);P^=1<<v
    expand((1<<n)-1,0)
    return best

def test(m,w,mark,X):
    Y,A,groups=compatibility(m,w,mark,X)
    mats=[A[np.ix_(g,g)] for g in groups]
    @lru_cache(None)
    def local_bound(j,B):
        ids=[i for i in range(len(groups[j])) if B>>i&1]
        return omega(mats[j][np.ix_(ids,ids)]) if ids else 0
    def bits(a):return sum(1<<i for i in np.flatnonzero(a).tolist())
    def bound(S):return sum(local_bound(j,bits(S[g])) for j,g in enumerate(groups))
    U=np.ones(len(Y),dtype=bool); base=bound(U); deleted=[]
    cell_types=Counter((len(g),local_bound(j,(1<<len(g))-1)) for j,g in enumerate(groups))
    if base>72:
        for b in range(len(w)):
            V=Y[:,m*b:m*(b+1)]; E=np.eye(m,dtype=np.int64)
            fs=list(E)+[E[i]-E[j] for i,j in combinations(range(m),2)]+[w[b]]
            if m==3:fs.extend(np.cross([1,1,1],v) for v in np.unique(V,axis=0))
            for f in fs:
                T=(V@f!=0)&U
                if np.any(T):
                    q=bound(T)
                    if q<=13:
                        U[T]=False;deleted.append((b,tuple(int(a) for a in f),q))
            if bound(U)<=72:break
    return base,bound(U),cell_types,deleted

def run(m):
    C=glue(m); records=profiles(m,C); reps=representatives(m,records)
    expected={4:(13192,29),3:(43980,88)}
    assert (len(records),len(reps))==expected[m]
    print(f'A{m-1}^{24//(m-1)}: {len(records)} normalized records; {len(reps)} subgroup orbits')
    checked=analytic=0; max_bound=0
    for k,((w,j,a),size) in enumerate(reps):
        ij=tuple(i for i,x in enumerate(w[j]) if x==a)[:2];mark=(j,ij)
        X=shell(m,C,w,mark); norms=sorted(norm(v) for v in w)
        # The rank-two support propagation proof handles precisely this shape.
        if m==3 and norm(w[j])==24:
            assert norms==[0]*6+[6]*5+[24] and len(X)==572
            analytic+=1; print(k,size,len(X),'ANALYTIC K6',flush=True);continue
        if len(X)<=144:
            base=final=len(X)//2;cells={};deleted=[]
        else:base,final,cells,deleted=test(m,w,mark,X)
        assert final<=72,(m,k,base,final,w,mark)
        checked+=1;max_bound=max(max_bound,final)
        print(k,size,len(X),base,final,flush=True)
        if m==4 and len(X) in (472,600) and (norm(w[j])==48 or norms==[0,0,12,12,12,12,16,32] and norm(w[j])==0):
            print('  cell types (size, clique number):',sorted(cells.items()))
        if deleted:print('  support deletions (component,functional,upper bound):',deleted)
    assert analytic==(0 if m==4 else 7)
    print(f'CHECKED {checked}; ANALYTIC {analytic}; MAX COMPUTED BOUND {max_bound}\n')

if __name__=='__main__':
    run(4)
    run(3)
