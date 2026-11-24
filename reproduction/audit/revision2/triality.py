"""Enlarge the subgroup on A5^4 D4 by all six D4 triality actions, lifting each on the fixed mark.
Run beside exceptional_{shells,orbits}.py with joint23 as argument.
"""
from itertools import product,permutations
from pathlib import Path
from exceptional_shells import codes,cases
from exceptional_orbits import normalize
import sys


def action(a,b):
    # Matrices are doubled, including every scaled D4 component vector.
    y=tuple(sum(r*t for r,t in zip(row,b)) for row in a)
    assert not any(t%2 for t in y)
    return tuple(t//2 for t in y)

def compose(a,b):
    return tuple(zip(*(action(a,x) for x in zip(*b))))

def label(v):
    assert all((z-v[0])%6==0 for z in v)
    return ((sum(v)//6)%2)+(2 if v[0]%6 else 0)

def generators():
    I=tuple(tuple(2*int(i==j) for j in range(4)) for i in range(4))
    F=I[:3]+(tuple(-z for z in I[3]),)
    H=((1,1,1,1),(1,1,-1,-1),(1,-1,1,-1),(1,-1,-1,1))
    roots=[v for v in product((-1,0,1),repeat=4) if sum(t*t for t in v)==2]
    W=[]
    for p in permutations(range(4)):
        for sg in product((-1,1), repeat=4):
            if sg.count(-1)%2==0:
                W.append(tuple(tuple(2*int(j==k)*s for j in range(4)) for k,s in zip(p,sg)))
    def sub(a):
        return tuple(label(action(a,v)) for v in ((0,0,0,0),(6,0,0,0),(3,3,3,3),(3,3,3,-3)))
    reps={sub(I):I};q=[I]
    while q:
        a=q.pop()
        for u in (H,F):
            b=compose(a,u)
            if sub(b) not in reps:q.append(b);reps[sub(b)]=b
    assert len(reps)==6
    T={v for v in roots if v[-1]==0 and sum(v[:3])==0}
    lift={}
    for act,rep in reps.items():
        for w in W:
            a=compose(w,rep)
            if {action(a,r) for r in T}==T:
                assert {action(a,r) for r in roots}==set(roots)
                assert sub(a)==act
                lift[act]=a;break
    assert len(lift)==6
    return lift


def full(root):
    n=6;words=set(codes(n));D=generators();rows=list(cases(root,6))
    for m in (0,4):
        maps=[]
        for p in permutations(range(4)):
            if m==0 and p[0]!=0:continue
            for s in product((-1,1),repeat=4):
                for d,mat in D.items():
                    if all(tuple(a*x[j]%6 for a,j in zip(s,p))+(d[x[-1]],) in words for x in words):
                        maps.append((p,s,mat))
        def key(x):
            return min(tuple(normalize(6,m,i,tuple(z*s[i] for z in x[j])) for i,j in enumerate(p))
                  +(normalize(6,m,4,action(mat,x[4])),) for p,s,mat in maps)
        orbits={}
        for j,k,x,c in rows:
            if j==m:orbits.setdefault(key(x),[]).append(int(k))
        print('marked factor',m,'actions',len(maps),'reduced orbit count',len(orbits),'example identified indices',next((group for group in orbits.values() if len(group)>1),None),'terminal class' ,next(g for g in orbits.values() if 32 in g))
if __name__=='__main__':full(Path(sys.argv[1]))
