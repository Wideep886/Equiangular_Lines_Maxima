#!/usr/bin/env python3
"""Check the short incidence table in the evaluation-code lemma (no evaluation-word search).
All loops range over integer intersection counts, at most five at a time
for the only undetermined triples of counts. See the manuscript formulas.
"""
from functools import lru_cache
from fractions import Fraction as F
from math import comb
from itertools import product

assert __debug__, 'Run without -O'
L = (253,77,21,5,1)
assert L == tuple(comb(23-j,4-j)//comb(7-j,4-j) for j in range(5))
assert 2**12*sum(comb(23,j) for j in range(4))==2**23

@lru_cache(None)
def A(k,u,v):
    if k==0:
        assert v==0
        return {L[u]}
    M1=(k-v)*L[u+1]+v*L[u]
    if k==7:
        return {2*L[u]-M1+4*(u==v)}
    assert k in (8,11)
    a=0 if k==8 else 1
    if u<=2:
        M2=2*comb(k-v,2)*L[u+2]+(2*v+1)*(k-v)*L[u+1]+v*v*L[u]
        val=F(L[u])-F(-M2+(2*a+4)*M1-a*(a+4)*L[u],2)
        assert val.denominator==1
        return {int(val)}
    assert u==3
    I=(a,a+2,a+4)
    out=set()
    for n0,n1 in product(range(6),repeat=2):
        n2=5-n0-n1
        if n2<0: continue
        ns=(n0,n1,n2)
        if any(n and (i<v or 7-i<u-v) for n,i in zip(ns,I)): continue
        if sum(n*i for n,i in zip(ns,I))!=M1: continue
        out.add(5-2*n1)
    assert out
    return out

expected = {
 (0,0):((0,),), (0,1):((176,),), (0,2):((112,),), (0,3):((128,),),
 (7,0):((148,),), (7,1):((112,),(164,)), (7,2):((148,),(128,),(100,)), (7,3):((128,),(132,),(144,),(148,)),
 (8,0):((176,),), (8,1):((128,),(100,)), (8,2):((128,),(148,),(144,)), (8,3):((112,128,144),(132,),(128,),(148,)),
 (11,0):((176,),), (11,1):((132,),(112,)), (11,2):((128,),(148,),(144,)), (11,3):((148,),(112,128,144),(132,),(144,))}
possible=set(); minimum_types=set()
for (k,t),expected_row in expected.items():
    row=[]
    for r in range(min(k,t)+1):
        sums={0}
        for u in range(t+1):
            for v in range(max(0,u-(t-r)), min(u,r)+1):
                for _ in range(comb(r,v)*comb(t-r,u-v)):
                    sums={x+(-2)**u*y for x in sums for y in A(k,u,v)}
        assert all((253-(-1)**t*x)%2==0 for x in sums)
        vals=tuple(sorted(k+(253-(-1)**t*x)//2 for x in sums))
        row.append(vals)
        possible.update(vals); possible.update(276-w for w in vals)
        if 100 in vals: minimum_types.add((k,t,r))
        if 176 in vals: minimum_types.add((23-k,t,t-r))
    assert tuple(row)==expected_row, (k,t,row)
    print('INCIDENCE', k,t,tuple(row))
assert possible=={0,100,112,128,132,144,148,164,176,276}
assert minimum_types=={(7,2,2),(8,1,1),(12,0,0),(15,0,0),(23,1,1)}
for k,t,r in minimum_types:
    assert t==r and (k-2*t)%3==0
    assert F(k,4)-F((k-2*t)**2,108)==F(5,3)
# For k=8, u=v, A=total-2*(number of intersections of size two).
assert [(L[u]-next(iter(A(8,u,u))))//2 for u in range(3)]==[168,42,6]
assert 168-2*42+6==90
print('LEECH INCIDENCE TABLE AND MINIMUM LIFTS VERIFIED; LOWER BOUND 90')
