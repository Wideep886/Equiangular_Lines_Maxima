#!/usr/bin/env python3
"""Extra independent small checks for the expanded audit. Requires SymPy and NumPy.
Run with --source pointing at the folder containing the original verifier modules.
None of these unit tests replaces the exhaustive proof or the cited classifications.
"""
from itertools import combinations, product
from collections import Counter
from fractions import Fraction as F
from math import comb, factorial
from pathlib import Path
import sys, argparse, hashlib
import sympy as sp
import numpy as np


def golay():
    # Binary cyclic (23,12) code, extended with parity. Coefficients indexed from 0.
    g=sum(1<<i for i in (0,1,5,6,7,9,11))
    code=[]
    for u in range(1<<12):
        v=0
        for i in range(12):
            if u>>i&1:v^=g<<i
        v|=(v.bit_count()%2)<<23
        code.append(v)
    assert len(set(code))==4096
    hist=Counter(x.bit_count() for x in code)
    assert hist=={0:1,8:759,12:2576,16:759,24:1},hist
    octads=[b for b in code if b.bit_count()==8]
    count=Counter(tuple(k for k in range(24) if (mask>>k)&1)
        for b in octads for mask in [] )  # empty; actual count below
    for b in octads:
        for I in combinations([i for i in range(24) if b>>i&1],5):count[I]+=1
    assert len(count)==comb(24,5) and set(count.values())=={1}
    assert all((a&b).bit_count() in (0,2,4) for a,b in combinations(octads,2))
    print('Binary Golay: weight distribution',dict(sorted(hist.items())), '; all 42504 five-subsets checked')
    # Asche's construction, as in Example 5.19 of Greaves et al., JCTA 138 (2016).
    B1,B2=next((a,b) for a,b in combinations(octads,2) if not (a&1 or b&1) and (a&b)==6)
    A=[b for b in octads if b&1 and not b&6 and (b&B1).bit_count()==(b&B2).bit_count()==2]
    def vec(b):return tuple(4*((b>>i)&1)-4*int(i==0)-1 for i in range(24))
    X=[vec(b) for b in A]
    dot=lambda a,b:sum(x*y for x,y in zip(a,b))
    normals=[tuple(4*int(i==0)+1 for i in range(24)),tuple(int(i==0)-int(i==1) for i in range(24)),tuple(int(i==0)-int(i==2) for i in range(24)),vec(B1),vec(B2)]
    assert len(X)==72 and len(set(X))==72
    assert all(dot(x,x)==80 and all(dot(x,n)==0 for n in normals) for x in X)
    assert all(abs(dot(a,b))==16 for a,b in combinations(X,2))
    assert sp.Matrix(normals).rank()==5 and sp.Matrix(X).rank()==19
    print('Asche lower bound: 72 integer rows, squared norm 80, pairings +/-16, exact rank 19')
    print('Fixed octads B1, B2 (one-based):',*[tuple(i+1 for i in range(24) if b>>i&1) for b in (B1,B2)])
    return X


def cartan(t):
    a,n=t;C=2*sp.eye(n)
    if a=='A':edges=[(i,i+1) for i in range(n-1)]
    elif a=='D':edges=[(i,i+1) for i in range(n-2)]+[(n-3,n-1)]
    else:edges=[(i,i+1) for i in range(n-2)]+[(2,n-1)]
    for i,j in edges:C[i,j]=C[j,i]=-1
    return C


def root_set(C):
    n=C.rows;cc=[list(map(int,C.row(i))) for i in range(n)]
    Q=[tuple(int(i==j) for i in range(n)) for j in range(n)];S=set(Q)
    for a in Q:
        for j in range(n):
            b=list(a);b[j]-=sum(a[i]*cc[i][j] for i in range(n));b=tuple(b)
            if b not in S:S.add(b);Q.append(b)
    assert all((sp.Matrix([a])*C*sp.Matrix(a))[0]==2 for a in S)
    return S


def weyl_order(C):
    n=C.rows;unseen=set(range(n));ans=1
    while unseen:
        i=min(unseen);unseen.remove(i);Q=[i]
        for v in Q:
            for j in list(unseen):
                if C[v,j]==-1:unseen.remove(j);Q.append(j)
        d=len(Q);degs={j:sum(C[j,k]==-1 for k in Q) for j in Q}
        if max(degs.values())<=2:ans*=factorial(d+1)
        else:
            b=next(j for j in Q if degs[j]==3);arms=[]
            for j in Q:
                if C[b,j]!=-1:continue
                u,v=b,j;k=1
                while degs[v]>1:
                    u,v=v,next(w for w in Q if w!=u and C[v,w]==-1);k+=1
                arms.append(k)
            arms.sort()
            if arms[:2]==[1,1]:ans*=2**(d-1)*factorial(d)
            else:assert arms==[1,2,d-4];ans*={6:51840,7:2903040,8:696729600}[d]
    return ans


def minuscule_data(slots):
    types=[('A',i) for i in range(1,20)]+[('D',i) for i in range(4,11)]+[('E',6),('E',7),('E',8)]
    for t in types:
        a,n=t; C=cartan(t); R=root_set(C); M=C.inv();h=len(R)//n
        assert len(R)==n*h
        pos=[x for x in R if all(v>=0 for v in x)]
        high=tuple(max(r[j] for r in pos) for j in range(n));assert high in R and all(v>0 for v in high)
        # For a dominant integral weight u, pairing with the highest root must be <=1.
        inds=[i for i,x in enumerate(high) if x==1]
        us=[sp.zeros(n,1)]+[sp.eye(n)[:,i] for i in inds]
        def cls(u):return tuple(F(x)%1 for x in M*u)
        def w(u):return F((2*u.T*M*u)[0])
        W={cls(u):w(u) for u in us};assert len(W)==len(us)==int(C.det())
        rec=[];groups={}
        for i,u in zip(inds,us[1:]):
            s=cls(u);v=w(u);delta=v-W[cls(2*u)]/4
            z=[j for j in range(n) if j!=i]
            orbit=weyl_order(C)//weyl_order(C.extract(z,z))
            rec.append((v,delta,orbit))
            # A type large subset and D spin group share one slot; D4 slots distinct.
            if a=='A':group='s' if delta==1 else 'l'
            elif a=='D' and n==4:group=str(i)
            elif a=='D':group='v' if i==0 else 'l'
            else:group='l'
            groups[s]=group
        if t==('E',8):assert not rec;continue
        expected={(x[0],F(x[1]),x[3]) for x in slots(t)}
        assert set(rec)==expected,(t,rec,expected)
        for u,v in product(us[1:],repeat=2):
            width=w(u)+w(v)-(W[cls(u-v)]+W[cls(u+v)])/2
            assert width>=2,(t,u,v,width)
            if groups[cls(u)]==groups[cls(v)]!='s':assert width>=4,(t,u,v,width)
        print(f'{a}{n}: roots={len(R)}, minuscule classes={len(us)-1}, data and width bounds PASS')
    print('Minuscule data: all 29 irreducible types checked, including E8')


def chains():
    # Obtain symmetric chains of a product of two chains recursively.
    def rect(a,b):
        if a==0 or b==0:return []
        first=[(0,j) for j in range(b)]+[(i,b-1) for i in range(1,a)]
        return [first]+[[(i+1,j) for i,j in c] for c in rect(a-1,b-1)]
    cs=[]
    for c in rect(2,3):
        cs.extend([[c[i]+(j,) for i,j in d] for d in rect(len(c),3)])
    assert len(cs)==5 and set(x for c in cs for x in c)==set(product(range(2),range(3),range(3)))
    assert sum(map(len,cs))==18 and all(all(all(a<=b for a,b in zip(x,y)) for x,y in zip(c,c[1:])) for c in cs)
    print('Explicit 5 chains of [2]x[3]x[3]:',cs)


def clique_unit(omega):
    count=0
    for n in range(1,6):
        ed=list(combinations(range(n),2))
        for b in range(1<<len(ed)):
            A=np.zeros((n,n),dtype=bool)
            for k,(i,j) in enumerate(ed):
                if b>>k&1:A[i,j]=A[j,i]=True
            brute=max(len(S) for k in range(n+1) for S in combinations(range(n),k) if all(A[i,j] for i,j in combinations(S,2)))
            assert omega(A)==brute
            count+=1
    print('Clique recursion unit test: all',count,'labeled simple graphs of orders 1..5 PASS')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,default=Path(__file__).resolve().parent.parent/'equiangular19_revised');p.add_argument('--output',type=Path,default=Path(__file__).resolve().parent/'Asche72_integer_vectors.txt');args=p.parse_args()
    sys.path.insert(0,str(args.source));from root19_mixed_deficit import slots
    from niemeier_19_certificate import omega
    X=golay();args.output.write_text('\n'.join(' '.join(map(str,x)) for x in X)+'\n')
    minuscule_data(slots);chains();clique_unit(omega)
    print('ALL ADDITIONAL CHECKS PASSED')
