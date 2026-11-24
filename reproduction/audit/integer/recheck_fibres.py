"""Read codes.txt and a labelled cuts list, checking every cut without automorphisms."""
from collections import Counter
from functools import reduce
from itertools import combinations, product
from operator import and_
from pathlib import Path
import sys

def components(blocks):
    ids=list(range(len(blocks)))
    def root(i):
        while ids[i]!=i: i=ids[i]
        return i
    for i,j in combinations(range(len(blocks)),2):
        overlap=(blocks[i]&blocks[j]).bit_count()
        assert overlap in (1,3)
        if overlap==3:
            ids[root(i)]=root(j)
    return [[b for j,b in enumerate(blocks) if root(j)==k] for k in set(map(root,range(len(ids))))]

def check_switch(B):
    # A component with four supports is a point times two disjoint binary choices of pairs.
    masks=Counter(sum(((b>>j)&1)<<i for i,b in enumerate(B)) for j in range(24))
    del masks[0]
    assert masks.pop(15)==1 and len(masks)==4
    assert all(k.bit_count()==2 and masks[k^15]==2 and n==2 for k,n in masks.items())

def check_stars(B):
    multiplicity=[sum((b>>i)&1 for b in B) for i in range(24)]
    assert sorted(k for k in multiplicity if k)==[4]*6+[6]*6
    core=sum(1<<i for i,k in enumerate(multiplicity) if k==6)
    stars={b&core for b in B}
    leaves={b&~core for b in B}
    assert len(stars)==4 and len(leaves)==3
    assert all(p.bit_count()==2 for p in leaves) and all(a&b==0 for a,b in combinations(leaves,2))
    assert all(s.bit_count()==3 for s in stars)
    edges={a&b for a,b in combinations(stars,2)}
    assert len(edges)==6 and all(e.bit_count()==1 for e in edges)
    assert set(B)=={a|b for a,b in product(stars,leaves)}

def signed_maximum(B):
    # Complete Bron--Kerbosch maximal-clique enumeration with pivoting.
    v=[]
    for b in B:
        indices=[j for j in range(24) if ((b>>j)&1)]
        v.extend((b,sum(s<<j for j,s in zip(indices[1:],signs))) for signs in product(range(2),repeat=4))
    nb=[set() for x in v]
    for j,(b,m) in enumerate(v):
        for i,(a,n) in enumerate(v[:j]):
            s=a&b
            if abs(s.bit_count()-2*(s&(m^n)).bit_count())==1:nb[j].add(i);nb[i].add(j)
    def expand(P,X):
        if not P:return -len(v) if X else 0
        u=max(P|X,key=lambda k:len(P&nb[k]))
        val=0
        for i in list(P-nb[u]):
            val=max(val,1+expand(P&nb[i],X&nb[i]))
            P.remove(i);X.add(i)
        return val
    return expand(set(range(len(v))),set())

def main():
    codes={}
    for row in Path(sys.argv[1]).read_text().splitlines():
        key,n,*words=row.split()
        codes[key]=[int(w) for w in words]
        assert len(codes[key])==int(n)
    # Tetrad components with their vertex and tetrad counts certify inequivalence.
    fingerprints=set()
    for words in codes.values():
        tetrads=[w for w in words if w.bit_count()==4]
        unseen=set(range(24));pieces=[]
        while unseen:
            start=min(unseen);component=1<<start;prev=0
            while component!=prev:
                prev=component
                for t in tetrads:
                    if t&component:component|=t
            unseen.difference_update(j for j in range(24) if (component>>j)&1)
            pieces.append((component.bit_count(),sum(bool(t&component) for t in tetrads)))
        fingerprints.add(tuple(sorted(pieces)))
    assert len(fingerprints)==len(codes)==9
    numbers=Counter();representatives={}
    for line in Path(sys.argv[2]).read_text().splitlines():
        code,one,zero=line.split();one,zero=int(one),int(zero)
        block=[w^one for w in codes[code] if w&(one|zero)==one]
        assert all(b.bit_count()==5 for b in block)
        pieces=components(block)
        shape=tuple(sorted(map(len,pieces)))
        if shape==(8,16):
            assert all(reduce(and_,p) for p in pieces)
        elif shape==(2,2,2,4,4,4):
            for p in pieces:
                if len(p)==4:check_switch(p)
                else:assert (p[0]&p[1]).bit_count()==3
                representatives[len(p)]=p
        else:
            assert shape==(12,12,12)
            for p in pieces:check_stars(p)
        numbers[shape]+=1
    assert numbers=={(8,16):2880,(2,2,2,4,4,4):11520,(12,12,12):32}
    print(sorted(numbers.items()))
    caps=[signed_maximum(p) for p in [[31],representatives[2],representatives[4]]]
    assert caps==[5,8,10]
    print(caps)
if __name__=='__main__':main()
