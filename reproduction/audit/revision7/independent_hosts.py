"""Verify binary labels, graph edges and cells against a separate shell generator.

Usage: python3 -B independent_hosts.py joint23 audit generated_input_directory
The three arguments are directories. Requires NumPy for section_isometries.chart.
"""
import argparse
from itertools import combinations
from pathlib import Path
import sys


def integral_labels(points):
    """Find the integer row lattice using unimodular Euclidean row operations."""
    rows = [list(p) for p in points]
    pivots = []
    for j in range(len(rows[0])):
        i = len(pivots)
        first = next((k for k in range(i,len(rows)) if rows[k][j]),None)
        if first is None:
            continue
        rows[i],rows[first] = rows[first],rows[i]
        for k in range(i+1,len(rows)):
            while rows[k][j]:
                q = rows[i][j] // rows[k][j]
                rows[i] = [v-q*w for v,w in zip(rows[i],rows[k])]
                rows[k],rows[i] = rows[i],rows[k]
        pivots.append(j)
    d = len(pivots)
    assert all(not any(r) for r in rows[d:])
    words = []
    for point in points:
        r = list(point)
        word = 0
        for i,p in enumerate(pivots):
            q,rem = divmod(r[p],rows[i][p]); assert rem == 0
            r = [x-q*v for x,v in zip(r,rows[i])]
            word |= (q%2)<<i
        assert not any(r)
        words.append(word)
    return d,words


def rank(words):
    pivots = {}
    for w in words:
        while w:
            i = w.bit_length()-1
            if i not in pivots:
                pivots[i]=w
                break
            w ^= pivots[i]
    return len(pivots)


def allowed_roots(scale,mark,blocks):
    roots=[]
    offset=0
    for i,block in enumerate(blocks):
        if not (scale==6 and i==4):
            for val in set(block):
                indices=[j+offset for j,v in enumerate(block) if v==val and (i!=mark or j>=3)]
                roots.extend(zip(indices,indices[1:]))
        offset += len(block)
    return roots


def in_signed_cone(v,roots):
    # Eliminate one ordered root at a time, solving the equation v=sum c_i r_i.
    r = list(v)
    sg = 0
    for i,j in roots:
        c = r[i]
        if c*sg<0:
            return False
        if c:
            sg = 1 if c>0 else -1
        r[j] += c; r[i] = 0
    return bool(sg) and not any(r)


def check_record(joint,table,n,mark,h,make_shell,directory):
    _,ys,_=make_shell(joint,n,mark,h)
    points=[tuple(map(int,p)) for p in ys if tuple(p)<tuple(-p)]
    roots=allowed_roots(n,mark,h)
    dim,words=integral_labels(points)
    conflict = [set() for _ in points]
    for i,j in combinations(range(len(points)),2):
        if abs(sum(a*b for a,b in zip(points[i],points[j]))) != 2*n*n:
            conflict[i].add(j);conflict[j].add(i)
    cells=[]
    remaining=set(range(len(points)))
    while remaining:
        stack=[min(remaining)]
        cell=set(stack)
        while stack:
            a=stack.pop()
            new=conflict[a]-cell
            cell.update(new);stack.extend(new)
        cells.append(sorted(cell));remaining-=cell
    source=(directory/f'scale{n}-record{table:02}.in').read_text().splitlines()
    D,G=map(int,source[0].split())
    assert D==dim and len(cells)==G and len(source)==1+2*G
    stored,derived=[],[]
    for k,indices in enumerate(cells):
        size,*labels=map(int,source[1+2*k].split())
        assert size==len(indices)==len(labels) and all(0<=w<1<<dim for w in labels)
        stored.extend(labels);derived.extend(words[i] for i in indices)
        expected=[0]*size
        for i,j in combinations(range(size),2):
            u,v=indices[i],indices[j]
            if u in conflict[v]:continue
            a,b=points[u],points[v]
            if any(in_signed_cone((x+sg*y for x,y in zip(a,b)),roots) for sg in (1,-1)):continue
            expected[i]|=1<<j;expected[j]|=1<<i
        assert expected==list(map(int,source[2+2*k].split()))
    assert rank(stored)==rank(derived)==rank([(x<<dim)|y for x,y in zip(stored,derived)])==dim
    print(f'scale {n} index {table}: axes {len(points)}, integral rank {dim}, cells {len(cells)}; all edges and binary kernels agree',flush=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    for k in ('joint','audit','inputs'):
        parser.add_argument(k,type=Path)
    a=parser.parse_args()
    if not __debug__:raise RuntimeError('Enable Python assertions.')
    sys.path.insert(0,str(a.audit.resolve()))
    from section_isometries import RECORDS, chart, cases
    done=set()
    for n,indices in RECORDS.items():
        for mark,k,h,_ in cases(a.joint,n):
            k=int(k)
            if k not in indices or (n==6 and mark!=4):continue
            assert (n,k) not in done
            check_record(a.joint,k,n,mark,h,chart,a.inputs)
            done.add((n,k))
    assert len(done)==15

if __name__=='__main__':
    main()
