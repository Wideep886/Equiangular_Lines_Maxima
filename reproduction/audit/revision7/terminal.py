"""Build a binary input directly from the coordinates of each terminal shell."""
from collections import defaultdict
from itertools import combinations
from pathlib import Path
import argparse
import gzip
import json
import sys

import numpy as np


def coordinate_groups(n, mark, w):
    """Simple-root strings from equal components of w; other coordinates stay fixed."""
    out = []
    offset = 0
    for i,block in enumerate(w):
        if n == 6 and i == 4:
            # This is always the marked D4; its marked A2 uses the first 3 coordinates.
            assert i == mark
            out.extend((offset+j,) for j in range(4))
        else:
            keys = defaultdict(list)
            for j,h in enumerate(block):
                key = ('single',j) if i == mark and j < 3 else ('tie',h)
                keys[key].append(offset+j)
            out.extend(tuple(items) for items in keys.values())
        offset += len(block)
    return out


def chart(v, groups):
    level, coeff = [], []
    for group in groups:
        total = 0
        for j in group[:-1]:
            total += v[j]
            coeff.append(total)
        level.append(total+v[group[-1]])
    return tuple(level), tuple(coeff)


def comparable(a,b):
    if a[0] != b[0] or a[1] == b[1]:
        return False
    return all(x<=y for x,y in zip(a[1],b[1])) or all(x>=y for x,y in zip(a[1],b[1]))


def build(n,mark,w,shell,size,destination,emit):
    flat = sum(w,())
    h = np.array(flat,dtype=np.int32)
    assert h@h == 6*n*n
    points = shell[shell@h==3*n*n].astype(int)
    assert len(points)==size and np.all(np.sum(points**2,axis=1)==4*n*n)
    pointset = set(map(tuple,points.tolist()))
    assert all(tuple(a-b for a,b in zip(flat,x)) in pointset for x in pointset)
    groups=coordinate_groups(n,mark,w)
    centers=[]
    charts=[]
    for x in sorted(pointset):
        opposite=tuple(a-b for a,b in zip(flat,x))
        if x>=opposite:
            continue
        y=tuple(2*a-b for a,b in zip(x,flat))
        by=chart(y,groups)
        my=chart(tuple(-a for a in y),groups)
        centers.append(y)
        charts.append((by,my))
    N=len(centers)
    assert 2*N == len(pointset) == size
    ymat=np.array(centers,dtype=np.int32)
    gram=ymat@ymat.T
    assert all(gram[i,i]==10*n*n for i in range(N))
    parents=list(range(N))
    def root(i):
        while parents[i]!=i:
            parents[i]=parents[parents[i]]
            i=parents[i]
        return i
    compatible=[0]*N
    for i,j in combinations(range(N),2):
        if abs(gram[i,j])!=2*n*n:
            parents[root(i)]=root(j)
        elif not any(comparable(charts[i][0],c) for c in charts[j]):
            compatible[i]|=1<<j
            compatible[j]|=1<<i
    groups_by_conflict=defaultdict(list)
    for i in range(N):
        groups_by_conflict[root(i)].append(i)
    cells=list(groups_by_conflict.values())
    print(f'scale {n}, axes {N}, cell sizes {sorted(set(map(len,cells)))}',flush=True)
    assert max(map(len,cells))<=36
    with destination.open('w') as stream:
        d,steps=emit(centers,cells,compatible,stream)
    print(f'  dimension {d}, saturation steps {steps}',flush=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('joint',type=Path)
    parser.add_argument('audit',type=Path,help='directory containing revision2 and revision6')
    parser.add_argument('destination',type=Path)
    parser.add_argument('--scale',type=int,default=0)
    a=parser.parse_args()
    if not __debug__:
        raise RuntimeError('Enable Python assertions.')
    sys.path[:0]=[str(a.audit.resolve()/'revision2'),str(a.audit.resolve()/'revision6'),str(a.joint.resolve()/'normal_sieve')]
    from exceptional_shells import cases,codes
    from verify_chains import full_shell
    from pack import emit
    a.destination.mkdir(parents=True,exist_ok=True)
    # Derive the complete frontier from the certified list partitions.
    certificate = a.audit / 'revision6' / 'exceptional_chains.json.gz'
    with gzip.open(certificate, 'rt', encoding='ascii') as src:
        records=json.load(src)
    assert len(records)==216
    frontier=defaultdict(list)
    for rec in records:
        if rec['points']-len(rec['links'])>=116:
            key=(rec['scale'], None if rec['mark']==-1 else rec['mark'])
            frontier[key].append(rec['index'])
    assert sum(map(len,frontier.values()))==15
    for (n,group),indices in sorted(frontier.items(),reverse=True):
        if a.scale and a.scale!=n:
            continue
        shell=full_shell(n,group,codes(n))
        data={int(i):(w,size) for m,i,w,size in cases(a.joint,n) if m==group}
        for i in sorted(indices):
            w,size=data[i]
            target=a.destination/f'scale{n}-record{i:02d}.in'
            print(target.name,flush=True)
            build(n,group,tuple(map(tuple,w)),shell,size,target,emit)


if __name__=='__main__':
    main()
