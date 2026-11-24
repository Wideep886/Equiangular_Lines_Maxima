"""Enumerate rank-eleven single-quintet support templates with or without a distinguished perfect colour.

A colour is a matching of K6 identified by its edge mask. The iteration chooses
the four other matchings as an unordered tuple after fixing the first perfect
matching. Canonical keys are the actual transformed five-matchings, not just
partition-refinement invariants. Run with Python 3.11+; no data files are used.
"""
from collections import Counter
from itertools import combinations, permutations, product


def mask(vertices):
    return sum(1 << i for i in vertices)


def binary_rank(matrix):
    pivot = {}
    for row in matrix:
        while row:
            i = row.bit_length() - 1
            if i not in pivot:
                pivot[i] = row
                break
            row ^= pivot[i]
    return len(pivot)


def field_rank(matrix, prime):
    work = [[v % prime for v in row] for row in matrix]
    i = 0
    for j in range(len(matrix[0])):
        k = next((k for k in range(i, len(work)) if work[k][j]), None)
        if k is None: continue
        work[k], work[i] = work[i], work[k]
        inv = pow(work[i][j], -1, prime)
        for k in range(i+1,len(work)):
            fac = work[k][j]*inv % prime
            if fac:
                work[k] = [(b-fac*a)%prime for a,b in zip(work[i],work[k])]
        i+=1
    return i


def check_signs(triples):
    """Switch every spanning-forest edge positive, then try all remaining signs."""
    n=len(triples)
    fixed,free=[],[]
    components=list(range(n))
    for i,j in combinations(range(n),2):
        if triples[i]&triples[j]: continue
        x,y=components[i],components[j]
        if x==y: free.append((i,j))
        else:
            fixed.append((i,j))
            components=[x if u==y else u for u in components]
    G=[[2*(i==j) for j in range(n)] for i in range(n)]
    for i,j in fixed:G[i][j]=G[j][i]=1
    stats=Counter()
    for labels in product((-1,1),repeat=len(free)):
        for (i,j),v in zip(free,labels):G[i][j]=G[j][i]=v
        stats['signings']+=1
        if field_rank(G,101)>7:continue
        # A modular rank of seven proves at least seven real dimensions.
        assert field_rank(G,101)==7
        stats['possible']+=1
        total_weight=0
        for i,j in product(range(5),range(-1,6)):
            mask=1<<i
            if j>=0: mask |= 1<<(5+j)
            if any((mask&t).bit_count()>1 for t in triples):continue
            indices=[k for k,t in enumerate(triples) if not mask&t]
            for signs in product((-1,1), repeat=max(0,len(indices)-1)):
                v=[0]*n
                for k,eps in zip(indices,(1,)+signs):v[k]=eps
                M=[row+[x] for row,x in zip(G,v)]
                M.append(v+[5-mask.bit_count()])
                if field_rank(M,101)<=7:
                    total_weight+=mask.bit_count()
        stats['max_weight']=max(stats['max_weight'],total_weight)
    return stats


def main():
    edges = tuple(combinations(range(6), 2))
    edge_id = {e:i for i,e in enumerate(edges)}
    matchings = {}
    for s in range(4):
        for js in combinations(range(15),s):
            used = set(v for j in js for v in edges[j])
            if len(used) == 2*s:
                matchings[mask(js)] = js
    assert len(matchings)==76
    full = mask(edge_id[e] for e in ((0,1),(2,3),(4,5)))
    acts = []
    for p in permutations(range(6)):
        edge_perm = [edge_id[tuple(sorted((p[i],p[j])))] for i,j in edges]
        acts.append({m:mask(edge_perm[j] for j in js) for m,js in matchings.items()})
    stab = [g for g in acts if g[full]==full]
    assert len(stab)==48
    rooted = set()
    unanchored = set()
    scanned = 0
    letters=[m for m in matchings if not (m & full)]
    for other in combinations(letters,4):
        if (full|other[0]|other[1]|other[2]|other[3]).bit_count()!=3+sum(len(matchings[m]) for m in other):
            continue
        ell=3+sum(len(matchings[m]) for m in other)
        if ell<11:continue
        scanned+=1
        rooted_key=min(tuple(sorted(g[m] for m in other)) for g in stab)
        if rooted_key in rooted: continue
        rooted.add(rooted_key)
        all_colours=(full,)+other
        unanchored.add(min(tuple(sorted(g[m] for m in all_colours)) for g in acts))
    size = lambda rows:sum(len(matchings[m]) for m in rows)
    rooted_counts = Counter(3+size(rows) for rows in rooted)
    counts = Counter(size(rows) for rows in unanchored)
    empty_counts = Counter(size(rows) for rows in unanchored if 0 in rows)
    ranks=Counter()
    total=Counter()
    for rows in sorted(unanchored):
        triples=[(1<<c) | ((1<<5)<<a) | ((1<<5)<<b) for c,m in enumerate(rows) for j in matchings[m] for a,b in [edges[j]]]
        gramm2 = [mask(j for j,t in enumerate(triples) if s&t==0) for s in triples]
        ranks[len(triples),binary_rank(gramm2)] += 1
        if binary_rank(gramm2)<=7:
            checks=check_signs(triples)
            for k in ('signings','possible'):total[k]+=checks[k]
            total['surviving_orbits']+=int(bool(checks['possible']))
            total['max_weight']=max(total['max_weight'],checks['max_weight'])
    summary=lambda ss: [ss[k] for k in range(11,16)]
    print('labelled decorations after fixing one perfect matching',scanned)
    print('anchored orbits for 11..15 occupied edges',summary(rooted_counts))
    print('unanchored orbits',summary(counts),'with empty colour',summary(empty_counts))
    print('rank mod 2, by number of triples', sorted(ranks.items()))
    print('mod 101 across the ten cases', sorted(total.items()))
    assert scanned == 9341
    assert summary(rooted_counts) == [177,53,13,1,1]
    assert summary(counts) == [132,31,7,1,1]
    assert summary(empty_counts) == [3,1,0,0,0]
    assert sum(v for (n,r),v in ranks.items() if r>7) == 162
    assert total == {'max_weight':4,'possible':5,'signings':1600,'surviving_orbits':4}

if __name__=='__main__':
    main()
