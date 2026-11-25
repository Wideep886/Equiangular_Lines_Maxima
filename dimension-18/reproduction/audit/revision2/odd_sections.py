"""Build all 5 A2^12 and the 2 terminal A4^6 binary inputs without column saturation.
For each mark use the scaled x rows modulo 2, quotienting by the nonzero w.
Requires exceptional_shells.py in the same folder. Usage: joint23 output_dir.
"""
from collections import defaultdict
from functools import cache
from itertools import combinations
from pathlib import Path
import sys
from exceptional_shells import blocks, codes, cases

def signed_shell(n,words,mark,w):
    target,inner=4*n*n,3*n*n
    out=set()
    for word in words:
        hist=[]
        for i,(label,h) in enumerate(zip(word,w)):
            groups=defaultdict(list)
            for e,v in blocks(n,n,label,mark==i):
                if e<=target:groups[e,sum(x*y for x,y in zip(v,h))].append(v)
            hist.append(groups)
        if any(not s for s in hist):continue
        @cache
        def tails(j,e,t):
            if j==len(word):return ((),) if (e,t)==(0,0) else ()
            rows=[]
            for (energy,ip),items in hist[j].items():
                if energy>e:continue
                for tail in tails(j+1,e-energy,t-ip):
                    rows.extend(v+tail for v in items)
            return tuple(rows)
        out.update(tails(0,target,inner));tails.cache_clear()
    return out


def make(n,mark,w,C,expected):
    flat=sum(w,());xset=signed_shell(n,C,mark,w)
    def partner(x):return tuple(a-b for a,b in zip(flat,x))
    assert len(xset)==expected and all(partner(x) in xset and partner(x)!=x for x in xset)
    xs=sorted(x for x in xset if x<partner(x))
    ys=[tuple(2*a-b for a,b in zip(x,flat)) for x in xs]
    blocks_=tuple(tuple(i*n+j for j in range(n) if z==h[j] and (i!=mark or (j<3)==flag))
            for i,h in enumerate(w) for z in sorted(set(h)) for flag in ((False,True) if i==mark else (False,)))
    blocks_=tuple(v for v in blocks_ if v)
    # For the quadratic pairing test, intersections of stabilizer orbits
    # need not be considered across cells since we only seek an upper bound.
    def cone(v):
        # Coordinates against dual fundamental weights of each equal-value A block.
        signs=set()
        for row in blocks_:
            t=0
            for j in row:
                t+=v[j]
                if t:signs.add(t>0)
            if t or len(signs)>1:return False
        return len(signs)==1
    ids=list(range(len(xs)))
    def disjoint(i):
        while ids[i]!=i:i=ids[i]
        return i
    compatible=[0]*len(xs)
    for i,j in combinations(range(len(xs)),2):
        a,b=ys[i],ys[j]
        if abs(sum(x*y for x,y in zip(a,b)))!=2*n*n:
            x,y=disjoint(i),disjoint(j);ids[x]=y;continue
        if not cone(tuple(x-y for x,y in zip(a,b))) and not cone(tuple(x+y for x,y in zip(a,b))):
            compatible[i]|=1<<j;compatible[j]|=1<<i
    cells=defaultdict(list)
    for i,x in enumerate(xs):
        def sort_key(v):
            return tuple(tuple(sorted(v[j] for j in B)) for B in blocks_)
        k=disjoint(i) if n==3 else min(sort_key(x),sort_key(partner(x)))
        cells[k].append(i)
    assert sum(map(len,cells.values()))==len(xs) and max(map(len,cells.values()))<=64
    # Each unmarked A factor has the last coordinate determined by its zero sum.
    # For A4 marked triple (a,a,a,b,-3a-b), use a and b.
    kept=[i*n+j for i in range(len(w)) for j in ((0,3) if i==mark else range(n-1))]
    assert len(kept)==22
    def bits(v):return sum((v[k]&1)<<i for i,k in enumerate(kept))
    h=bits(flat);assert h
    def quotient(v):
        z=bits(v);j=h.bit_length()-1
        if (z>>j)&1:z^=h
        return ((z>>(j+1))<<j)+(z% (1<<j))
    assert all(quotient(x)==quotient(partner(x)) for x in xs)
    vwords=list(map(quotient,xs));piv={}
    for z in vwords:
        while z:
            j=z.bit_length()-1
            if j not in piv:piv[j]=z;break
            z^=piv[j]
    assert len(piv)==21
    lines=['21 '+str(len(cells))]
    for c in cells.values():
        lines.append(' '.join(map(str,[len(c)]+[vwords[i] for i in c])))
        lines.append(' '.join(str(sum(((compatible[i]>>j)&1)<<u for u,j in enumerate(c))) for i in c))
    return '\n'.join(lines)+'\n'

def run(root,folder):
    folder.mkdir(parents=True,exist_ok=True)
    for n in (3,5):
        for mark,k,w,num in cases(root,n):
            if n==5 and k not in (26,27):continue
            data=make(n,mark,w,codes(n),num)
            (folder/f'a{n-1}-row{k}.in').write_text(data)
            print('A'+str(n-1),'index',k,'signed',num,'axes',num//2,'binary rank',data.splitlines()[0].split()[0],
                  'cells',data.splitlines()[0].split()[1],flush=True)
if __name__=='__main__':run(Path(sys.argv[1]),Path(sys.argv[2]))
