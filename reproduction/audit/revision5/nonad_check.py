import sys
from itertools import product, permutations
from collections import defaultdict,Counter
from nonad_cells import codewords,dot,shell
words=list(codewords())
mini={s:sorted(set(permutations((s,s,-2*s)))) for s in (-1,1)}
w=next(w for w in words if sum(map(abs,w))==9)
h=[mini[s][0] if s else (0,0,0) for s in w]
o=[i for i,s in enumerate(w) if s==0]+[i for i,s in enumerate(w) if s!=0]
blacks=defaultdict(list)
signed=set()
selected=set()
for v in shell(words,mini):
    if sum(dot(vv,ww) for vv,ww in zip(v,h))!=27: continue
    x=tuple(tuple(2*a-b for a,b in zip(v[i],h[i])) for i in o)
    assert x not in signed
    signed.add(x)
    cen=tuple(dot(xx,h[i])//6 for xx,i in zip(x[2:],o[2:]))
    idx=tuple(i for i,e in enumerate(cen) if not e)
    if len(idx)==1:
        j=next(j for j,s in enumerate(x[:2]) if any(s))
        if min(x[j])!=-4: continue
    else:
        assert len(idx)==3
        j=2
        if next(z for z in cen if z)<0:continue
    selected.add(x);blacks[(j,idx)].append((x,cen))
neg=lambda x:tuple(tuple(-c for c in y) for y in x)
assert len(signed)==408 and len(selected)==204 and signed==selected|{neg(x) for x in selected}
tri={}; steps={}
for (j,idx),data in blacks.items():
    points,cs=list(zip(*data));N=len(points); cen=cs[0];assert all(c==cen for c in cs)
    for i in range(11):
        loc={x[i] for x in points}
        if i<2:
            assert [sum(x[i][m] for x in points) for m in range(3)]==[0]*3
            if i==j: assert tri.setdefault(i,loc)==loc
            continue
        u=i-2
        assert all(sum(x[i][m] for x in points)==N*cen[u]*h[o[i]][m] for m in range(3))
        if u in idx:
            a,b=sorted(loc);diff=tuple(y-z for y,z in zip(a,b))
            assert steps.setdefault(u,diff)==diff
            assert dot(diff,h[o[i]])==0
assert len(tri)==2 and len(steps)==9
print("PASS 408 distinct signed points, 204 axes paired by reflection, all barycenters exactly in H, nine global binary steps and two consistent ternary triangles")
