"""Regenerate the seven out-neighbour sets from the archived A3^8 coordinates.
Invocation: python3 mark5_table.py path/to/eleven/a3
"""
from fractions import Fraction
from pathlib import Path
import sys
sys.path.insert(0,str(Path(sys.argv[1]).resolve()))
from octacode import cells, dot, shells, vertices
h=((-1,-1,-1,3),(-3,1,1,1),(-3,1,1,1),(-3,1,1,1),
   (-1,-1,-1,3),(-3,1,1,1),(-3,1,1,1),(-3,1,1,1))
arcs={}
shapes=list(cells(h,shells(4)))
assert len(shapes)==42
for c in shapes:
    E=tuple(i for i,g in enumerate(c) if len(g)>1)
    assert len(E)==2 and 0 not in E and all(len(c[i])==3 for i in E)
    for i in E:
        ds={dot(tuple(a-b for a,b in zip(x,y)),tuple(a-b for a,b in zip(x,y)))
            for x in c[i] for y in c[i] if x!=y}
        assert ds=={32} # coordinates scaled by four
    vs=vertices(c)
    n=len(vs)
    total=[2*sum(col)-n*b for col,b in zip(zip(*vs),sum(h,()))]
    coef=[]
    for i,b in enumerate(h):
        m=total[4*i:4*i+4]; t=Fraction(dot(m,b), n*dot(b,b))
        assert tuple(m)==tuple(t*n*x for x in b)
        coef.append(t)
    t0=coef[0]
    assert abs(t0)==1
    if t0==1:continue
    assert all(abs(t)==(Fraction(1,3) if j in E else 1) for j,t in enumerate(coef))
    i,j=sorted(E,key=lambda j:-coef[j])
    assert coef[i]==-coef[j]==Fraction(1,3) and E not in arcs
    arcs[E]=(i,j)
assert len(arcs)==21 
for i in range(1,8):
    nei=sorted(b for a,b in arcs.values() if a==i)
    assert len(nei)==3
    print(i, *nei)
