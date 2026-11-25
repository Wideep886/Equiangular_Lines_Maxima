"""Check a rank-three configuration of nine nonzero step vectors on AG(2,3).
Python 3.9+. No coordinate-shell data are needed for this finite table.
"""
from collections import Counter, defaultdict
from itertools import combinations, product
from math import gcd

POINTS=list(product(range(3), repeat=2))
LINES=[L for L in combinations(range(9),3)
       if all(sum(POINTS[i][j] for i in L)%3==0 for j in (0,1))]
PENCILS=[tuple([L]+[M for M in LINES if set(L).isdisjoint(M)]) for L in LINES]
PENCILS=[p for p in PENCILS if p[0]==min(p)]
SIGNS=list(product((-1,1),repeat=9))
SMALL_SIGNS=list(product((-1,1),repeat=3))

def partitions(n, prefix=()):
    if not n:
        yield prefix
    else:
        for c in range(max(prefix, default=-1)+2):
            yield from partitions(n-1, prefix+(c,))

def canonical_names(p):
    names={}
    return tuple(names.setdefault(x,len(names)) for x in p)

def affine_group():
    group=[]
    for a,b,c,d in product(range(3),repeat=4):
        if (a*d-b*c)%3:
            for s,t in POINTS:
                group.append([3*((a*x+b*y+s)%3)+(c*x+d*y+t)%3 for x,y in POINTS])
    assert len(group)==432
    return group

def relevant_partitions():
    group=affine_group()
    seen=set()
    for p in partitions(9):
        if p in seen: continue
        seen.update(canonical_names([p[i] for i in g]) for g in group)
        ns=sorted(Counter(p).values(), reverse=True)
        if len(ns)<4: continue
        mono=sum(len({p[i] for i in line})==1 for line in LINES)
        if 42+2*sum(ns[:3])+mono>=58:
            assert len(ns)<=5
            yield p
    assert len(seen)==21147

def normalize(v):
    if not any(v): return v
    g=gcd(*v)
    if next(x for x in v if x)<0: g=-g
    return tuple(x//g for x in v)

def exterior(v,w):
    return normalize(tuple(v[i]*w[j]-v[j]*w[i] for i,j in combinations(range(len(v)),2)))

def table(p):
    k=max(p)+1
    def row(T,signs):
        return tuple(sum(s for a,s in zip(T,signs) if a==j) for j in range(k))
    # A full pencil gives a signed relation using all nine occurrences.
    fullrow=[normalize(row(p,s)) for s in SIGNS]
    relations=set(fullrow)
    triangles={tr:{normalize(row(tr,s)) for s in SMALL_SIGNS} for tr in combinations(range(k),3)}
    for rr in triangles.values():relations.update(rr)
    cell_vectors={L:[row([p[i] for i in L],s) for s in SMALL_SIGNS] for L in LINES}
    difference={L:[[normalize(tuple(a-b for a,b in zip(x,y))) for y in vals] for x in vals] for L,vals in cell_vectors.items()}
    for rows in difference.values():
        for rr in rows:relations.update(rr)
    bad={row([j],[1]) for j in range(k)}
    for i,j in combinations(range(k),2):
        for s in (-1,1):bad.add(row((i,j),(1,s)))
    # Precompute which cube sign index is induced by each nine-position sign.
    local_choice={L:[SMALL_SIGNS.index(tuple(s[i] for i in L)) for s in SIGNS] for L in LINES}
    weights=[p.count(j) for j in range(k)]
    # Flip the direction of each representative so that this first relation is nonnegative.
    for v in sorted(w for w in set(fullrow) if min(w)>=0 and w not in bad):
        wed={w:exterior(v,w) for w in relations}
        zero=exterior(v,v)
        forbidden={exterior(v,w) for w in bad}
        # In rank three, at most k-3 independent relations can occur.
        projections=[None]
        if k==5:projections+=sorted(set(wed.values())-forbidden-{zero})
        for extra in projections:
            belongs={w for w,g in wed.items() if g==zero or g==extra}
            pair_count=max([sum(sorted(weights,reverse=True)[:2])]+
                            [sum(weights[j] for j in tr) for tr,rr in triangles.items() if rr & belongs])
            good_signs=[i for i,w in enumerate(fullrow) if w in belongs]
            fibers={L:[sum(d in belongs for d in D) for D in rr] for L,rr in difference.items()}
            cap=[]
            for pencil in PENCILS:
                two=sum(sorted([max(fibers[L]) for L in pencil],reverse=True)[:2])
                three=max(sum(fibers[L][local_choice[L][i]] for L in pencil) for i in good_signs)
                cap.append(max(two,three))
            yield 18+2*pair_count+sum(cap),v,extra,pair_count,cap

def main():
    orbits=defaultdict(list);high=[]
    for p in relevant_partitions():
        results=list(table(p))
        shape=tuple(sorted(Counter(p).values(),reverse=True))
        orbits[shape].append(max([0]+[r[0] for r in results]))
        for cap,*rest in results:
            if cap>=58:high.append((cap,p,*rest))
    assert sum(map(len,orbits.values()))==31 and len(high)==1
    data=high[0]
    assert data[:5]==(58,(0,0,0,0,0,0,1,2,3),(0,1,1,1),None,7)
    assert sorted(data[5])==[6,6,6,8]
    for shape,scores in sorted(orbits.items(),reverse=True):
        print(''.join(map(str,shape)),len(scores),max(scores))
    print('31 partition orbits; only capacity 58: six equal steps and a line of three steps with sum zero.')
if __name__=='__main__':main()
