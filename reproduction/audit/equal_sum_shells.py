"""Optional coordinate check for the two eight-block cases. Argument: extracted joint23."""
from ast import literal_eval
from collections import Counter, defaultdict
from fractions import Fraction
from itertools import combinations, product
from pathlib import Path
import sys


def rank(rows):
    basis = {}
    for row in rows:
        row = list(map(Fraction, row))
        for i, p in sorted(basis.items()):
            q = row[i]
            if q:
                row = [u-q*v for u, v in zip(row, p)]
        if any(row):
            pivot = next(i for i, u in enumerate(row) if u)
            basis[pivot] = [u/row[pivot] for u in row]
    return len(basis)


def verify(name, h, signed, scale, triples, key):
    points = set(signed)
    assert len(points) == 444 and sum(a*a for a in h) == 6*scale**2
    assert all(sum(u*u for u in x) == 4*scale**2 and
               sum(u*v for u, v in zip(h,x)) == 3*scale**2 for x in points)
    small, groups = [], defaultdict(list)
    for x in sorted(points):
        y = tuple(v-u for v,u in zip(h,x))
        assert y != x and y in points
        k = key(x)
        if k is None:
            assert key(y) is None
            if x < y:
                small.append(tuple(2*u-v for u,v in zip(x,h)))
        elif k < key(y):
            groups[k].append(x)
    assert len(small) == 6 and sorted(map(len, groups.values())) == [27]*8
    planes, centers, edges = defaultdict(list), [], []
    for group in groups.values():
        centers.append(tuple(2*sum(col)-27*a for col,a in zip(zip(*group),h)))
        support = []
        for i, idx in enumerate(triples):
            values = sorted({tuple(x[j] for j in idx) for x in group})
            if len(values) == 1:
                continue
            assert len(values) == 3
            support.append(i)
            for a,b in combinations(values,2):
                diff = [0]*len(h)
                for j,u,v in zip(idx,a,b):
                    diff[j] = u-v
                assert sum(u*u for u in diff) == 2*scale**2
                planes[i].append(diff)
        assert len(support) == 3
        all_idx = {j for i in support for j in triples[i]}
        assert all(len({x[j] for x in group}) == 1 for j in range(len(h)) if j not in all_idx)
        edges.append(set(support))
    assert len(planes) == 8 and all(rank(rows) == 2 for rows in planes.values())
    assert all(len(e & f) <= 1 for e,f in combinations(edges,2))
    assert sorted(Counter(i for edge in edges for i in edge).values()) == [3]*8
    directions = [row for rows in planes.values() for row in rows]
    assert rank(centers) == 3 and rank(directions) == 16 and rank(directions+centers) == 19
    assert all(rank(directions+centers+[x]) == 20 for x in small)
    print(name, '444 signed, 6 + 8*27 axes; degrees 3^8, pair overlap <= 1; ranks 16 + 3, each small raises rank to 20')


def main(joint):
    sys.path.insert(0, str(joint/'a2'))
    from golay import vectors, polar, shortened, dot
    c = shortened()
    tag, h = list(polar(c))[1]
    assert tag == 'three_roots' and h == (-3,0,3)*3+(0,)*24
    pts = (x for x in vectors(c) if dot(x,h)==27)
    triples = [range(3*i,3*i+3) for i in range(3,11)]
    verify('A2^12 (three_roots)', h, pts, 3, triples,
           lambda x: tuple(x[::3][i]%3 for i in range(11)) if x[0]%3 else None)
    sys.path.insert(0, str(joint/'a5'))
    from a5 import groups
    row, = [s for s in (joint/'a5/all.out').read_text().splitlines() if s.startswith('4 32 ')]
    mark, which, count, _, data = row.split(maxsplit=4)
    w = literal_eval(data)
    assert (mark,which,count) == ('4','32','222') and w == ((-3,)*3+(3,)*3,)*4+((0,)*4,)
    pts = (sum(z,()) for g in groups(w,4) for z in product(*g))
    triples = [range(3*i,3*i+3) for i in range(8)]
    verify('A5^4D4 (4,32)', sum(w,()), pts, 6, triples,
           lambda x: tuple(-x[i]%6 for i in range(0,24,6)) if not any(x[24:]) else None)


if __name__ == '__main__':
    main(Path(sys.argv[1]).resolve())
