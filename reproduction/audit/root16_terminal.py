"""Recheck rank-16 terminal capacities from saved word domains by explicit local coordinates.
Usage: python3 root16_terminal.py /path/to/EL18_58_search
The nine root types, 21 domains, weights, and attaining label lists are read from
root_profiles/orbit16.json. The claimed weights and weighted maxima are tested.
The completeness of the supplied word domains is a separate earlier computation.
"""
from fractions import Fraction
from functools import lru_cache
from itertools import combinations, product
from math import lcm
from pathlib import Path
import json
import sys

@lru_cache(None)
def component(t, k):
    # Coordinates below are divided by the denominator, then multiplied by sqrt(2).
    n = int(t[1:])
    if t.startswith('A'):
        N = n + 1
        return N, [tuple(N * int(j in J) - k for j in range(N))
                   for J in combinations(range(N), k)]
    assert t == 'D4'
    if k == 0:
        return 2, [(0, 0, 0, 0)]
    if k == 1:
        return 2, [tuple(2 * sign * int(i == j) for j in range(4))
                   for i in range(4) for sign in (-1, 1)]
    return 2, [v for v in product((-1, 1), repeat=4) if v.count(-1) % 2 == k - 2]

@lru_cache(None)
def interval(t, k, j, denominator):
    a, U = component(t, k)
    b, V = component(t, j)
    # Fixing a single representative in U is justified by Weyl transitivity.
    scores = [2 * denominator * sum(x*y for x,y in zip(U[0], v)) for v in V]
    assert all(s % (a*b) == 0 for s in scores)
    return min(scores)//(a*b), max(scores)//(a*b)

def cutoff(neighbors, weights, target, domain=None):
    # Greedy independent colour classes yield a weighted clique upper bound.
    if domain is None:
        domain = (1 << len(weights)) - 1
    calls = 0
    def search(P, need):
        nonlocal calls
        calls += 1
        if need <= 0:
            return True
        if not P:
            return False
        available = P
        order = []
        bound = 0
        while available:
            choices = available
            start = len(order)
            maximum = 0
            while choices:
                bit = choices & -choices
                i = bit.bit_length()-1
                available ^= bit
                choices &= ~bit & ~neighbors[i]
                maximum = max(maximum, weights[i])
                order.append((i,0))
            bound += maximum
            for q in range(start,len(order)):
                order[q] = order[q][0],bound
        for i,b in reversed(order):
            if b < need:
                return False
            if search(P & neighbors[i], need - weights[i]):
                return True
            P ^= 1 << i
        return False
    return search(domain,target), calls

@lru_cache(None)
def within(shape, weight):
    ds, parts = zip(*(component(t,k) for t,k in shape)) if shape else ((),())
    m = lcm(1,*ds)
    points = [tuple(coord*(m//d) for d,v in zip(ds,word) for coord in v)
              for word in product(*parts)]
    adj = [0]*len(points)
    for i,v in enumerate(points):
        for j in range(i):
            d = sum((x-y)**2 for x,y in zip(v,points[j]))
            if d in (4*m*m,6*m*m):
                adj[i] |= 1 << j
                adj[j] |= 1 << i
    # Vertex transitivity permits fixing point 0.
    w = [1]*len(points)
    yes,calls1 = cutoff(adj,w,weight-1,adj[0])
    no,calls2 = cutoff(adj,w,weight,adj[0])
    assert yes and not no
    return calls1 + calls2

def verify(label, record):
    name, _, p = label.split(':')
    p = int(p)
    types = name.split('.')
    bound, weights, words, witness = record
    # A common scaling of the component dot products, determined independently.
    m = lcm(4,*(int(t[1:])+1 for t in types if t.startswith('A')))
    D = m*m
    scale = None
    for row,number in zip(words,weights):
        shape = tuple(sorted((t,min(k,int(t[1:])+1-k) if t.startswith('A') else k)
                             for t,k in zip(types,row) if k))
        within(shape,number)
        h = 5*D - sum(interval(t,k,k,D)[1] for t,k in zip(types,row))
        if h:
            y,b = row[-2:]
            s = Fraction(D*(y*y+p*b*b),h)
            assert s>0 and (scale is None or scale==s)
            scale=s
        else:
            assert row[-2:]==[0,0]
    assert scale is not None
    N,Q=scale.numerator,scale.denominator
    adj = [0]*len(words)
    for i,v in enumerate(words):
        for j,u in enumerate(words[:i]):
            inn = D*Q*(u[-2]*v[-2] + p*u[-1]*v[-1])
            low=high=inn
            for t,k,l in zip(types,v,u):
                vmin,vmax = interval(t,k,l,D)
                low += N*vmin
                high += N*vmax
            if low in (-D*N,D*N) and high in (-D*N,D*N):
                adj[i] |= 1 << j
                adj[j] |= 1 << i
    assert sum(weights[i] for i in witness)==bound
    assert all((adj[i] >> j)&1 for i,j in combinations(witness,2))
    no,nodes=cutoff(adj,weights,bound+1)
    assert not no, label
    return len(words),nodes

def main():
    root=Path(sys.argv[1])
    tab=json.loads((root/'root_profiles'/'orbit16.json').read_text())
    maxima={}
    total=nodes=0
    for key,val in sorted(tab.items()):
        count,steps=verify(key,val)
        name=key.split(':')[0]
        maxima[name]=max(maxima.get(name,0),val[0])
        total+=count
        nodes+=steps
    print('domains',len(tab),'word occurrences',total,'distinct local weight problems',within.cache_info().currsize)
    print('weighted cutoff nodes',nodes,'maxima',sorted(maxima.values()))
if __name__=='__main__':
    main()
