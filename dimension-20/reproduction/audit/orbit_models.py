"""Integer/rational local weight tables for the orbit-energy proof.
No full local or global norm-four shell is enumerated for the general models.
Only dominant weights and parabolic-dominant (w,x) pairs are generated.
"""
from fractions import Fraction as F
from itertools import product, permutations
from functools import lru_cache
from collections import Counter, defaultdict
from math import prod, lcm, factorial
from pathlib import Path
import json
# The elementary Cartan routines and the small-code models are copied from v4.
from cartan_utils import TYPES, cartan, inverse_and_det, dot, matvec, weyl_order, budgets, small_local, H4, P3, G3, G4

@lru_cache(None)
def dominant(K,cap):
    """All nonnegative integer u with u^T K u <= cap, K positive entrywise."""
    n=len(K);assert all(a>0 for row in K for a in row);out=[]
    def go(u,N):
        i=len(u)
        if i==n:
            u=tuple(u);assert dot(u,matvec(K,u))==N;out.append((N,u));return
        cross=2*sum(K[i][j]*a for j,a in enumerate(u));b=0
        while N+cross*b+K[i][i]*b*b<=cap:
            go(u+[b],N+cross*b+K[i][i]*b*b);b+=1
    go([],0);return tuple(out)

def connected(C,ids):
    todo=set(ids);out=[]
    while todo:
        i=min(todo);todo.remove(i);q=[i]
        for j in q:
            t={k for k in todo if C[j][k]};todo-=t;q.extend(sorted(t))
        out.append(tuple(sorted(q)))
    return tuple(out)

def submatrix(C,I):return tuple(tuple(C[i][j] for j in I) for i in I)

@lru_cache(None)
def parabolic(C,Z):
    ans=[]
    for I in connected(C,Z):
        A=submatrix(C,I);d,K=inverse_and_det(A);order=weyl_order(A)
        records=[]
        # On a global valid centered vector the projected energy is <=5.
        # Its pairings on these roots are twice those of x, divided by sqrt(2).
        for N,u in dominant(K,5*d//2):
            assert 2*N<=5*d
            J=tuple(i for i,a in enumerate(u) if not a)
            den=weyl_order(submatrix(A,J));assert order%den==0
            e=F(2*N,d);m=order//den
            records.append((u,e,m))
        ans.append((I,tuple(records)))
    return tuple(ans)

class Factor:
    def __init__(self,ty):
        self.ty=ty;self.n=ty[1];self.C=cartan(ty);d,K=inverse_and_det(self.C);self.d,self.K=d,K
        D=dominant(K,7*d);R=sorted({tuple(a%d for a in matvec(K,u)) for N,u in D});assert len(R)==d
        self.idx={r:i for i,r in enumerate(R)};self.dom=[(N,self.idx[tuple(a%d for a in matvec(K,u))],u) for N,u in D]
        self.min=[min(N for N,c,u in self.dom if c==j) for j in range(d)]
        assert all(sum(N==self.min[j] and c==j for N,c,u in self.dom)==1 for j in range(d))
        assert sum(N==2*d and c==0 for N,c,u in self.dom)==1
        self.add=[[self.idx[tuple((a+b)%d for a,b in zip(x,y))] for y in R] for x in R]
        self.byclass=defaultdict(list);self.hist=defaultdict(Counter);order=weyl_order(self.C)
        for j,(N,c,u) in enumerate(self.dom):
            self.byclass[c].append((N,j))
            if N<=6*d:
                den=weyl_order(submatrix(self.C,tuple(i for i,a in enumerate(u) if not a)));assert order%den==0
                self.hist[c][N]+=order//den
        self.neg=[];lookup={u:j for j,(N,c,u) in enumerate(self.dom)}
        for N,c,u in self.dom:
            v=self.move(tuple(-a for a in u),tuple(range(self.n)));self.neg.append(lookup[v])
    def move(self,v,Z):
        # Reflection until dominant with respect to the parabolic subgroup.
        # Terminates because a vector in a finite Weyl orbit strictly increases
        # its pairing with a regular dominant element at every step.
        while any(v[i]<0 for i in Z):
            i=next(i for i in Z if v[i]<0);v=tuple(b-v[i]*c for b,c in zip(v,self.C[i]))
        return v
    @lru_cache(None)
    def options(self,j,s):
        """Parabolic-dominant x, directly from its pairings with simple roots."""
        N,_,w=self.dom[j];W=matvec(self.K,w);d,K=self.d,self.K
        Z=tuple(i for i,a in enumerate(w) if not a);J=tuple(i for i,a in enumerate(w) if a)
        P=parabolic(self.C,Z);out=defaultdict(list)
        for U in product(*(R for I,R in P)):
            if sum(e for u,e,m in U)>5:continue
            es=tuple((len(I),e,m) for (I,R),(u,e,m) in zip(P,U) if e)
            size=prod(m for u,e,m in U)
            a=[0]*self.n
            for (I,R),(u,e,m) in zip(P,U):
                for i,b in zip(I,u):a[i]=b
            for v in product(range(-2,3),repeat=len(J)):
                for i,b in zip(J,v):a[i]=b
                u=tuple(a);Ku=matvec(K,u);n=dot(u,Ku);n2=N+n-2*dot(u,W)
                if not (0<=n<=4*d and 0<=n2<=4*d):continue
                y=tuple(2*b-c for b,c in zip(u,w));q=F(dot(y,matvec(K,y)),2*d)
                if q>5:continue
                assert q>=sum(e for _,e,_ in es)
                c=self.idx[tuple(b%d for b in Ku)];partner=self.move(tuple(b-a for a,b in zip(u,w)),Z)
                out[c].append((n*s//d,n2*s//d,u,partner,es,size))
        for c,L in out.items():assert len({t[2] for t in L})==len(L)
        return out

@lru_cache(None)
def factor(ty):return Factor(ty)

def validated_general(obj):
    idx=obj['idx'];ty=tuple(tuple(t) for t in obj['types']);assert ty==tuple(TYPES[idx]);L=tuple(factor(t) for t in ty);s=lcm(*(x.d for x in L))
    C=tuple(tuple(c) for c in obj['glue']);CS=set(C)
    assert len(C)==len(CS) and len(C)**2==prod(x.d for x in L) and (0,)*len(L) in CS
    for c in C:
        assert len(c)==len(L) and all(type(a)==int and 0<=a<x.d for x,a in zip(L,c))
        n=sum(x.min[a]*(s//x.d) for x,a in zip(L,c));assert n%(2*s)==0 and (n>=4*s or not any(c))
    assert all(tuple(x.add[a][b] for x,a,b in zip(L,c,d)) in CS for c in C for d in C)
    # Shell cardinalities via local Weyl orbit sizes, not shell enumeration.
    hist=Counter()
    for c in C:
        H={0:1}
        for x,a in zip(L,c):
            new=Counter()
            for n,b in H.items():
                for m,k in x.hist[a].items():
                    v=n+m*(s//x.d)
                    if v<=6*s:new[v]+=b*k
            H=new
        hist.update(H)
    n2,n4,n6=(hist[j*s] for j in (2,4,6))
    def roots(t):
        a,n=t;return n*(n+1) if a=='A' else 2*n*(n-1) if a=='D' else {6:72,7:126,8:240}[n]
    assert n2==sum(map(roots,ty)) and n4==196560-24*n2 and n6==16773120+252*n2
    return L,C,s
