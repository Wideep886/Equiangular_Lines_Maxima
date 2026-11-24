"""Modular independent check of the 70 anchors of marked A5^4 D4.

Input the candidate cells from a5/bound.py, but reconstruct both the
positive cone and the cell graphs. Only modular rank and pair constraints
are used after the anchor, so a rank increase under lifting is harmless.
"""
from functools import cache
from itertools import combinations
from pathlib import Path
import sys

if len(sys.argv) != 3:
    raise SystemExit('python3 -B a5_anchors_mod.py PATH/TO/joint23 PRIME')
sys.path.insert(0, str(Path(sys.argv[1]) / 'a5'))
from bound import ws, orbs

prime = int(sys.argv[2])
assert prime in (101, 103)
w = ws(4)[32]
h = sum(w, ())
cells = orbs(w, 4)
v = [x for cell in cells for x in cell]
y = [tuple(2*a-b for a,b in zip(x,h)) for x in v]
blocks = [tuple(j for j in range(6*i, 6*i+6) if h[j]==t)
          for i in range(4) for t in sorted(set(h[6*i:6*i+6]))]
assert all(len(b)==3 for b in blocks)

def in_cone(z):
    if any(z[-4:]) or any(sum(z[j] for j in b) for b in blocks):
        return False
    prefixes = [z[b[0]] for b in blocks]+[z[b[0]]+z[b[1]] for b in blocks]
    return (any(prefixes) and (min(prefixes)>=0 or max(prefixes)<=0))

def compatible(i,j):
    if i==j or sum(a*b for a,b in zip(v[i],v[j])) not in (36,72):
        return False
    return all(not in_cone(tuple(a+s*b for a,b in zip(y[i],y[j]))) for s in (-1,1))

active = [i for i,z in enumerate(y) if not in_cone(z)]
m = len(v)
adj = [sum(1<<j for j in active if compatible(i,j)) for i in range(m)]
groups=[]
pos=0
for cell in cells:
    groups.append(((1<<len(cell))-1)<<pos)
    pos+=len(cell)

@cache
def capacity(mask):
    if not mask:
        return 0
    bit=mask & -mask
    i=bit.bit_length()-1
    rest=mask^bit
    return max(capacity(rest),1+capacity(rest&adj[i]))

def full(mask,k):
    if not k:
        yield 0
        return
    while bin(mask).count("1")>=k:
        bit=mask&-mask
        mask ^= bit
        for extra in full(mask&adj[bit.bit_length()-1],k-1):
            yield bit|extra

def red(basis,row):
    row=list(row)
    for j,b in basis:
        q=row[j]%prime
        if q:
            row=[(z-q*x)%prime for z,x in zip(row,b)]
    return tuple(z%prime for z in row)

def add(basis,row):
    u=red(basis,row)
    if any(u):
        j=next(i for i,z in enumerate(u) if z)
        f=pow(u[j],-1,prime)
        return basis+[(j,tuple(z*f%prime for z in u))]
    return basis

def projective(v):
    s=next((t for t in v if t),0)
    if not s:
        return v
    f=pow(s,-1,prime)
    return tuple(t*f%prime for t in v)

active_mask=sum(1<<i for i in active)
large=[]
for g in groups:
    cap=capacity(g&active_mask)
    if bin(g).count("1")==27:
        forced=list(full(g&active_mask,7))
        assert cap==7 and len(forced)==1
        large.append(forced[0])
    else:
        assert bin(g).count("1")==2 and cap<=2
assert len(large)==8 and len(groups)==11
from collections import Counter
ranks, bounds, totals = Counter(), Counter(), Counter()
for chosen in combinations(large,4):
    anchor=sum(chosen)
    basis=[]
    elig=active_mask
    for j in range(m):
        if anchor>>j&1:
            elig &= adj[j]
            basis=add(basis,y[j])
    assert all(((anchor^(1<<j))&~adj[j])==0 for j in range(m) if anchor>>j&1)
    assert 16<=len(basis)<=18
    ranks[len(basis)]+=1
    residue={}
    for j in range(m):
        if elig>>j&1:
            u=projective(red(basis,y[j]))
            residue[u]=residue.get(u,0)|(1<<j)
    z=(0,)*len(h)
    constant=residue.pop(z,0)
    directions=list(residue)
    sections={constant}
    if len(basis)<18:
        sections.update(constant|mask for mask in residue.values())
    if len(basis)==16:
        for u in directions:
            first=add([],u)
            quotient={}
            for t in directions:
                if t==u:
                    continue
                q=projective(red(first,t))
                quotient[q]=quotient.get(q,0)|residue[t]
            sections.update(constant|residue[u]|mask for mask in quotient.values())
    upper=max(28+sum(capacity(s&g) for g in groups) for s in sections)
    totals[len(basis)] += len(sections)
    bounds[upper]+=1
assert ranks == {16:40,17:6,18:24}
assert totals == {16:4360,17:78,18:24}
assert bounds == {56:70}
print('prime',prime,'anchors by rank',sorted(ranks.items()),'sections',sorted(totals.items()),
      'upper bounds WITHOUT any root test',sorted(bounds.items()))
