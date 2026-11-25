"""Independent verification on actual shell permutations, without composing abstract normal maps."""
from pathlib import Path
from collections import Counter
import sys
if len(sys.argv)!=3:raise SystemExit('usage: python3 normal_anchor_check.py path/to/audit/revision2 path/to/EL18_58_search/joint23')
audit=Path(sys.argv[1]).resolve();root=Path(sys.argv[2]).resolve()
sys.path.insert(0,str(audit))
from exceptional_shells import codes,cases
from exceptional_orbits import make_transforms
from odd_sections import signed_shell, make
import numpy as np
def closure(gens):
    e=tuple(range(1,12));todo=[e];seen={e}
    while todo:
        perm=todo.pop()
        for step in gens:
            composition=tuple(step[abs(v)-1]*(1 if v>0 else -1) for v in perm)
            if composition not in seen:seen.add(composition);todo.append(composition)
    return seen

def local(g,w):
    slots=[]
    for target,old in enumerate(g):
        src=abs(old)-1;sign=1 if old>0 else -1
        values=[x*sign for x in w[src]]
        if sorted(values)!=sorted(w[target]):return None
        buckets={v:iter(sorted([u for u in range(3) if values[u]==v],reverse=sign<0)) for v in values}
        slots.extend(src*3+next(buckets[v]) for v in w[target])
    return slots

def tables(model):
    rows=model.splitlines();_,n=map(int,rows[0].split())
    cells=[]
    for j in range(n):
        _,*ws=map(int,rows[2*j+1].split())
        adj=list(map(int,rows[2*j+2].split()))
        cells.append((ws,adj))
    return cells

C=codes(3);generators=make_transforms(3,C,None)
grp=closure([tuple((j+1)*s for j,s in zip(p,sign)) for p,sign,_ in generators])
def binary(v):
    return sum((x%2)<<i for i,x in enumerate([x for i,x in enumerate(v) if i%3!=2]))
for m,k,w,n in cases(root,3):
    if k not in(2,4):continue
    flat=sum(w,()); allx=signed_shell(3,C,None,w)
    xs=sorted(v for v in allx if v<tuple(x-y for x,y in zip(flat,v)))
    stab=[g for g in grp if local(g,w)]
    # Each gluing symmetry fixes the actual w after choice of local reflection and maps
    # each simple root (successive equal slots) to a positive simple root.
    roots=set()
    for i,block in enumerate(w):
        for val in block:
            occ=[3*i+j for j,v in enumerate(block) if v==val]
            roots.update((a,b) for a,b in zip(occ,occ[1:]))
    def apply(g,bs,v):
        return tuple((1 if g[i//3]>0 else -1)*v[old] for i,old in enumerate(bs))
    z=binary(flat);p=z.bit_length()-1
    def q(x):
        y=binary(x)
        if (y>>p)&1:y^=z
        return ((y>>(p+1))<<p) | (y%(1<<p))
    Z=[q(v) for v in xs];ids={z:i for i,z in enumerate(Z)}
    data=tables(make(3,m,w,C,n))
    assert len(ids)==len(xs)
    cells={frozenset(words) for words,adj in data}
    edges={frozenset((words[i],words[j])) for words,adj in data for i in range(len(words)) for j in range(i) if adj[i]>>j&1}
    axis_perms=[]
    for g in stab:
        b=local(g,w)
        assert apply(g,b,flat)==flat
        changed=set()
        for a,c in roots:
            mid=b.index(a);end=b.index(c);assert mid//3==end//3
            changed.add((mid,end) if g[mid//3]>0 else (end,mid))
        assert roots==changed
        Y=[apply(g,b,v) for v in xs]
        assert all(y in allx for y in Y)
        # Every actual axis and signed-pair relation has corresponding F2 image.
        P={a:q(y) for a,y in zip(Z,Y)}
        axis_perms.append(tuple(ids[P[z]] for z in Z))
        assert set(P.values())==set(Z)
        assert {frozenset(P[z] for z in c) for c in cells}==cells
        assert {frozenset(P[z] for z in pair) for pair in edges}==edges
    print(k,'stabilizer',len(stab),'unique axes',len(xs),'simple root count',len(roots),'cells',len(cells),'allowed internal edges',len(edges),flush=True)

    # Fourier inversion of the point indicator gives the size of every hyperplane section.
    spectrum=np.zeros(2**21,dtype=np.int16)
    spectrum[Z]=1
    k2=1
    while k2<len(spectrum):
        row=spectrum.reshape(-1,2*k2);x=row[:,:k2].copy();y=row[:,k2:].copy()
        row[:,:k2]=x+y;row[:,k2:]=x-y;k2*=2
    counts=(len(Z)+spectrum)//2
    candidates=np.flatnonzero((7*counts>=3*len(Z)+232)&(np.arange(len(counts))>0))
    def clique(mask,adj,cache):
        if mask in cache:return cache[mask]
        b=mask&-mask
        cache[mask]=max(clique(mask^b,adj,cache),1+clique(mask&adj[b.bit_length()-1],adj,cache))
        return cache[mask]
    cdata=[(words,adj,{0:0}) for words,adj in data]
    def bnd(normal):
        return sum(clique(sum(1<<i for i,wrd in enumerate(ws) if (normal&wrd).bit_count()%2==0),adj,cache) for ws,adj,cache in cdata)
    strong=[int(norm) for norm in candidates if bnd(int(norm))>=58]
    masks={sum(1<<i for i,z in enumerate(Z) if (x&z).bit_count()%2==0) for x in strong}
    assert len(masks)==len(strong)
    sizes=[]
    while masks:
        section=min(masks)
        block={sum(1<<pr[i] for i in range(len(Z)) if section>>i&1) for pr in axis_perms}
        assert block<=masks
        sizes.append(len(block));masks-=block
    print(k,'Fourier count',len(candidates),'strong',len(strong),'actual section orbits',len(sizes),dict(sorted(Counter(sizes).items())),flush=True)
