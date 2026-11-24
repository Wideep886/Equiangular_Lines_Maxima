"""Certify every marked polarization representative against a freshly enumerated shape universe.
Place beside exceptional_shells.py and pass EL18_58_search/joint23 as sole argument.
"""
from collections import defaultdict
from itertools import product, permutations
from functools import cache
from pathlib import Path
import sys
from exceptional_shells import cases, blocks, codes

@cache
def normalize(n,mark,block,v):
    if n==6 and block==4:
        # The D4 Weyl group permutes coordinates and changes an even number
        # of signs. With an A2 mark on the first triple, the first three
        # signs can be flipped together with the fourth sign.
        w=([abs(v[0])]*3+[abs(v[3])]) if mark==block else sorted(map(abs,v))
        if all(v) and sum(t<0 for t in v)%2:
            w[-1]*=-1
        return tuple(w)
    if mark==block:
        return v[:3]+tuple(sorted(v[3:]))
    return tuple(sorted(v))


def universe(n,mark):
    c=codes(n);d=len(c[0]);S=set()
    def coordinate_patterns(i,label):
        length=4 if n==6 and i==4 else n
        return set((e,normalize(n,mark,i,v)) for e,v in blocks(n,length,label,mark==i))
    parts={(i,k):coordinate_patterns(i,k) for i in range(d) for k in set(x[i] for x in c)}
    target=6*n*n
    for word in c:
        p=[parts[(i,k)] for i,k in enumerate(word)]
        tail=[0]*(d+1)
        if any(not s for s in p):continue
        for i in range(d-1,-1,-1):tail[i]=tail[i+1]+min(e for e,v in p[i])
        def rec(i,e,row):
            if i==d:
                if e==target:S.add(row)
            elif e+tail[i]<=target:
                for a,v in p[i]:
                    if e+a+tail[i+1]<=target:rec(i+1,e+a,row+(v,))
        rec(0,0,())
    return S


def make_transforms(n,wlist,marked):
    code=set(wlist); d=len(wlist[0])
    if n==3:
        gens=[((2,5,3,8,4,6,9,7,0,1,10),(1,-1,1,1,1,1,1,-1,1,-1,1)),
              ((2,5,3,0,6,8,7,1,4,10,9),(1,-1,-1,-1,1,-1,1,1,1,1,1))]
        maps=[]
        for p,s in gens:
            inv=tuple(p.index(i) for i in range(d))
            maps.append((inv,tuple(s[j] for j in inv),0))
        maps.append((tuple(range(d)),(-1,)*d,0))
    else:
        maps=[];a=d-(n==6)
        # Columns of a real isometry: arbitrary signs on A factors and
        # for D4 an optional change of the last coordinate sign.
        init=([w for w in wlist if w[:4].count(0)==3 and w[:4].count(1)==1]
              if n==4 else sorted(code))
        if n==4:assert len(init)==4
        for perm in permutations(range(a)):
            p=perm+((a,) if n==6 else ())
            if marked is not None and p[marked]!=marked:continue
            if n==4:
                # Signs are determined on the odd positions of a generator
                # by its image. This filter reduces the direct loop without
                # assuming a list of automorphisms.
                v=tuple(init[0][j] for j in p)
                feasible=[u for u in wlist if all(min(x,4-x)==min(y,4-y) for x,y in zip(u,v))]
                opts=set()
                for u in feasible:
                    sets=[(1,-1) if y%2==0 else ((1,) if x==y else (-1,)) for x,y in zip(u,v)]
                    opts.update(product(*sets))
            else:opts=product((1,-1),repeat=a)
            for s in opts:
                signs=s+((1,) if n==6 else ())
                for f in range(2 if n==6 else 1):
                    def ok(x):
                        y=tuple(signs[i]*x[p[i]]%n for i in range(a))
                        if n==6:y+=(x[-1]^(f if x[-1]>=2 else 0),)
                        return y in code
                    if all(ok(x) for x in init):maps.append((p,signs,f))
    # Verify each accepted operation on every word, independently of the
    # early generating-row filter above.
    for p,s,f in maps:
        for x in code:
            a=d-(n==6)
            y=tuple(s[i]*x[p[i]]%n for i in range(a))
            if n==6:y+=(x[-1]^(f if x[-1]>=2 else 0),)
            assert y in code
    return maps


def full(root):
    for n in (6,5,4,3):
        rows=list(cases(root,n)); code=codes(n)
        for m in sorted({r[0] for r in rows},key=str):
            maps=make_transforms(n,code,m)
            S=universe(n,m);pending=set(S)
            def image(row,iso):
                p,s,f=iso;result=[]
                for i,j in enumerate(p):
                    v=row[j]
                    if s[i]<0:v=tuple(-z for z in v)
                    if n==6 and i==4 and f:v=v[:-1]+(-v[-1],)
                    result.append(normalize(n,m,i,v))
                return tuple(result)
            reps=0
            for mark,k,w,num in rows:
                if mark!=m:continue
                w=tuple(normalize(n,m,i,v) for i,v in enumerate(w))
                assert w in pending
                if n==3:
                    orbit={w};stack=[w]
                    while stack:
                        t=stack.pop()
                        for g in maps:
                            t1=image(t,g)
                            if t1 not in orbit:stack.append(t1);orbit.add(t1)
                else:orbit={image(w,g) for g in maps}
                assert w in orbit and orbit<=pending
                pending-=orbit;reps+=1
            assert not pending
            print(n,'mark',m,'number of maps',len(maps),'component shapes',len(S),'disjoint marked orbits',reps)
if __name__=='__main__':full(Path(sys.argv[1]))
