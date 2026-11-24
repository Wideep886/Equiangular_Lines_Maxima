"""Second calculation for the four exceptional hosts, using distances in reflection graphs.
Usage: python poset_extra_orbits.py joint23 revision2 [table.txt]
Compare its per-mark table with poset_extra.py. Both depend on the same local
lattice model, glue and list of polarizations. This code determines the local
orbit polynomial from distances to a dominant point in a simple-reflection
graph; it makes no calls to the sign-imbalance calculation or poset_width.py.
"""
from collections import Counter, deque
from functools import cache
from itertools import combinations
from pathlib import Path
import sys
sys.path.insert(0,sys.argv[2])
from exceptional_shells import blocks,codes,cases

def simples(n,w,mark):
    d=len(w);D=(n==6 and d==4);good=set()
    for i,j in combinations(range(d),2):
        for s in ((1,-1) if D else (-1,)):
            r=tuple(int(t==i)+s*int(t==j) for t in range(d))
            if mark and len(set(r[:3]))>1:continue
            if not sum(a*b for a,b in zip(w,r)):good.add(r)
    return [r for r in sorted(good) if not any(tuple(a-b for a,b in zip(r,v)) in good for v in good)]

def reflect(v,r):
    t=sum(a*b for a,b in zip(v,r))
    return tuple(a-t*b for a,b in zip(v,r))

@cache
def local(n,w,mark,tag):
    roots=simples(n,w,mark)
    cand={v for e,v in blocks(n,len(w),tag,mark) if e<=4*n*n}
    dom=[v for v in cand if all(sum(a*b for a,b in zip(v,r))>=0 for r in roots)]
    tab=Counter()
    for v in dom:
        assert v in cand
        cand.remove(v); queue=deque([(v,0)]); lens=Counter()
        while queue:
            x,t=queue.popleft();lens[t]+=1
            for r in roots:
                b=reflect(x,r)
                if b in cand:
                    cand.remove(b);queue.append((b,t+1))
        poly=tuple(lens[k] for k in range(1+max(lens)))
        assert poly==poly[::-1] and max(poly)==poly[len(poly)//2]
        tab[sum(x*x for x in v),sum(a*b for a,b in zip(v,w)),poly]+=1
    assert not cand
    return [(k[0],k[1],k[2],t) for k,t in sorted(tab.items())]

@cache
def mult(f,g):
    q=[0]*(len(f)+len(g)-1)
    for i,a in enumerate(f):
        for j,b in enumerate(g):q[i+j]+=a*b
    return tuple(q)

def bound(n,w,mark,words):
    H=4*n*n;dot=3*n*n
    total=cap=0
    for word in words:
        tabs=[local(n,h,i==mark,g) for i,(h,g) in enumerate(zip(w,word))]
        floors=[min((r[0] for r in part),default=H+1) for part in tabs]
        if sum(floors)>H:continue
        def walk(i,e,ip,poly,m):
            nonlocal total,cap
            if i==len(tabs):
                if (e,ip)==(H,dot):
                    total+=m*sum(poly); cap+=m*poly[len(poly)//2]
                return
            ceiling=H-sum(floors[i+1:])
            for f,d,q,k in tabs[i]:
                if e+f<=ceiling:walk(i+1,e+f,ip+d,mult(poly,q),m*k)
        walk(0,0,0,(1,),1)
    return total,cap

def main():
    records=[]
    for n in (6,5,4,3):
        words=codes(n); alive=[];table=[]
        for mark,k,w,total in cases(Path(sys.argv[1]),n):
            mass,b=bound(n,w,mark,words);assert mass==total
            records.append((n,-1 if mark is None else mark,int(k),total,b))
            table.append(b)
            if b>=116:alive.append((mark,int(k),b))
        print(n,len(table),max(table),len(alive),alive,flush=True)
    if len(sys.argv)>3:
        Path(sys.argv[3]).write_text(''.join(' '.join(map(str,t))+'\n' for t in sorted(records)))
if __name__=='__main__':main()
