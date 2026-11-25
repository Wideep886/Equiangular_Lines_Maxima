"""Rebuild the full norm-four shells and check the exceptional chain certificates.

Usage: python3 verify_chains.py EL18_58_search/joint23 audit/revision2 exceptional_chains.json.gz
Requires NumPy. The prescribed mark list and glue words come from revision2.
New code independently rebuilds every point and checks each directed edge by
coordinate subtraction. No reflection order, middle-layer bound, or matching
algorithm is used to certify the upper bounds.
"""
import argparse
import gzip
import json
from collections import defaultdict
from functools import cache
from itertools import product
from pathlib import Path
import sys

import numpy as np


@cache
def local_shell(n, size, label, marked, d_factor):
    limit = 4*n*n
    if d_factor:
        residue = (n//2) * int(label >= 2)
    else:
        residue = (label if n in (3,4) else -label) % n
    choices = [a for a in range(-2*n,2*n+1) if a % n == residue]
    groups = defaultdict(list)
    for xs in product(choices, repeat=size):
        if marked and len(set(xs[:3])) != 1:
            continue
        if d_factor:
            if sum(xs) % (2*n) != n*(label%2):
                continue
        elif sum(xs) != 0:
            continue
        norm = sum(a*a for a in xs)
        if norm <= limit:
            groups[norm].append(xs)
    return groups


def full_shell(n, mark, codes):
    """List all norm-four points for the indicated mark before imposing w."""
    limit = 4*n*n
    dimension = lambda i: 4 if n == 6 and i == 4 else n
    points = []

    def append_parts(groups, index, budget, parts, floors):
        if index == len(groups):
            if budget == 0:
                points.extend(sum(row, ()) for row in product(*parts))
            return
        for energy, local in groups[index].items():
            if energy + floors[index+1] <= budget:
                parts.append(local)
                append_parts(groups,index+1,budget-energy,parts,floors)
                parts.pop()

    for word in codes:
        groups = [local_shell(n, dimension(i), label, i==mark, n==6 and i==4)
                  for i,label in enumerate(word)]
        if any(not g for g in groups):
            continue
        floors = [0]*(len(groups)+1)
        for i in reversed(range(len(groups))):
            floors[i] = floors[i+1] + min(groups[i])
        if floors[0] <= limit:
            append_parts(groups,0,limit,[],floors)
    assert len(set(points)) == len(points)
    table = np.asarray(sorted(points), dtype=np.int32)
    assert np.all((table*table).sum(axis=1)==limit)
    expected = {3:98208,4:98154,5:98100,6:98046}
    assert len(table)==expected[n]
    return table


def check_record(n, mark, w, points, links):
    """Each difference is a positive multiple of one allowed coordinate root."""
    flat_w = sum(w,())
    assert len(flat_w)==len(points[0])
    component = [(j, k) for j,block in enumerate(w) for k in range(len(block))]
    left, right = set(), set()
    for i,j in links:
        assert 0<=i<len(points) and 0<=j<len(points) and i!=j
        assert i not in left and j not in right
        left.add(i);right.add(j)
        diff = [(k, b-a) for k,(a,b) in enumerate(zip(points[i],points[j])) if a!=b]
        assert len(diff)==2
        (u,mu),(v,nu)=diff
        block, first = component[u]
        block2, second = component[v]
        assert block==block2 and first<second and mu>0 and abs(nu)==mu
        assert nu==-mu or (n==6 and block==4)
        assert flat_w[u]*mu+flat_w[v]*nu==0
        assert block!=mark or first>=3
    return len(points)-len(links)


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('joint',type=Path)
    ap.add_argument('revision2',type=Path)
    ap.add_argument('certificate',type=Path)
    args=ap.parse_args()
    if not __debug__:
        raise SystemExit('Enable Python assertions.')
    sys.path.insert(0,str(args.revision2.resolve()))
    from exceptional_shells import cases, codes
    with gzip.open(args.certificate,'rt') as stream:
        records=json.load(stream)
    keys = [(r['scale'],r['mark'],r['index']) for r in records]
    lookup=dict(zip(keys,records))
    assert len(keys)==len(lookup)==216
    # Only these records are sent to the separate geometric arguments.
    unresolved = {
        3: {(-1, i) for i in (1,2,3,4)},
        4: {(0, i) for i in (0,1,5,6,10,14,15,16)},
        5: {(0,26),(0,27)},
        6: {(4,32)},
    }
    used=set()
    for n in (6,5,4,3):
        common=defaultdict(list)
        for mark,i,w,size in cases(args.joint,n):
            common[mark].append((int(i),tuple(map(tuple,w)),size))
        complete=[]
        for mark,specs in common.items():
            all_points=full_shell(n,mark,codes(n))
            for i,w,size in specs:
                height = np.asarray(sum(w,()),dtype=np.int32)
                points = all_points[(all_points@height)==3*n*n].tolist()
                key = (n,-1 if mark is None else mark,i)
                cert=lookup[key]
                assert len(points)==size==cert['points']
                cap=check_record(n,mark,w,points,cert['links'])
                assert (cap>=116) == ((key[1],key[2]) in unresolved[n]), key
                complete.append((size,len(cert['links']),cap))
                used.add(key)
        print(f'scale {n}: {len(complete)} records, {sum(s for s,_,_ in complete)} signed points, '
              f'{sum(l for _,l,_ in complete)} directed edges, {sum(c<116 for _,_,c in complete)} bounds below 116.')
    assert used==set(lookup)
    print('All 216 chain covers verified from integer differences on the full coordinate shells.')
if __name__=='__main__':
    main()
