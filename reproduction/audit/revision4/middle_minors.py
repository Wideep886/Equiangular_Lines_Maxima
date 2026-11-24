"""Verify the middle-layer determinants and two extra dimensions for A3 #6 and #10.
Invocation: python3 middle_minors.py path/to/eleven/a3
"""
import sys
from itertools import combinations
from collections import Counter
from pathlib import Path
sys.path.insert(0,str(Path(sys.argv[1]).resolve()))
from chamber import data, shells
from saturated import maximum_sets, basis
cases=[
    (((-1,-1,-1,3),(-3,1,1,1),(-3,1,1,1),(-2,-2,2,2),
      (-2,-2,2,2),(-2,-2,2,2),(0,0,0,0),(-1,-1,-1,3)),
     [(2,5,6,7,10,17),(1,8,9,11,14,15)], [3,4,12,13,16,18], 6, []),
    (((0,0,0,0),(-5,-1,3,3),(-3,1,1,1),(-3,1,1,1),
      (-1,-1,-1,3),(-2,-2,2,2),(0,0,0,0),(0,0,0,0)),
     [(9,11,12,14),(9,12,15,16),(11,14,15,16)], [5,6,7],1,[(-1,3,-1,-1)])
]
def determinant(a):
    a=[row[:] for row in a];n=len(a);sign=1;previous=1
    for k in range(n-1):
        i=next(i for i in range(k,n) if a[i][k])
        if i!=k:a[k],a[i]=a[i],a[k];sign=-sign
        p=a[k][k]
        for i in range(k+1,n):
            for j in range(k+1,n):
                t=a[i][j]*p-a[i][k]*a[k][j]
                assert t%previous==0
                a[i][j]=t//previous
        previous=p
    return sign*a[-1][-1]
for case,(h,families,extras,factor,line) in zip((6,10),cases):
    w=sum(h,())
    pts,ch,_,_,capacity=data(case,h,shells(4))
    assert capacity==59
    y=[tuple(2*a-b for a,b in zip(x,w)) for x in pts]
    sets=[]
    for ids,graph in ch:
        k,items=maximum_sets(graph)
        sets.append([tuple(ids[i] for i in item) for item in items])
    def proj(r):
        return y[r][4*factor:4*factor+4]
    assert not set.intersection(*(set(blocks) for blocks in families))
    for blocks in families:
        rows=[]
        for bi in blocks:
            matches,=sets[bi]
            rows.extend(sorted(matches))
        assert all((not any(proj(i))) if not line else proj(i)==line[0] for i in rows)
        chosen=[]
        for i in rows:
            if len(basis([y[t] for t in chosen+[i]]))>len(chosen):chosen.append(i)
        R=basis([y[i] for i in chosen])
        cols=[next(i for i,a in enumerate(r) if a) for r in R]
        assert len(chosen)==17
        det=determinant([[y[r][c] for c in cols] for r in chosen])
        print('mark',case,'blocks',*blocks,'rows',*chosen,'cols',*cols,'det / 2^v',end=' ')
        v=0
        while det%2==0:v+=1;det//=2
        print(det,'v',v)
    assert len(extras)>=2
    total=0
    for block in extras:
        for r,s in sets[block]:
            assert len(basis([proj(r),proj(s)]+line))==2+len(line)
            total+=1
    print('mark',case,'extra full pairs verified',total)
