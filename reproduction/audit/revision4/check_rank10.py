"""Hand-table replay for no-five-rank10.tex. Python 3 standard library, no input or seed."""
from fractions import Fraction
from itertools import combinations, combinations_with_replacement, product, permutations
from math import comb


def packing(s):
    if s == 5:
        return 2
    return s * ((s - 1) // 2) // 3


def profiles():
    for u in range(9):
        for a, b, c in combinations_with_replacement(range(11 - u), 3):
            if u + a + b + c == 10 and b:
                yield (u, a, b, c)


def caps(n):
    u, a, b, c = n
    h = min(a*b, b*c, c*a)
    if n == (1, 3, 3, 3):
        h = min(h, 7)
    return [packing(u)] + [min(comb(j, 2), u * (j // 2)) for j in n[1:]] + [h]


def pair_capacities(n, counts):
    u, a, b, c = n
    l0, la, lb, lc, h = counts
    return (2*(sum(comb(x, 2) for x in n) - 3*l0 - la - lb - lc),
            2*(u*a + b*c - 2*la - h), 2*(u*b + c*a - 2*lb - h),
            2*(u*c + a*b - 2*lc - h))


def colour_table():
    rows = {}
    for n in profiles():
        bounds = caps(n)
        if sum(bounds) < 10:
            continue
        for m in range(10, sum(bounds) + 1):
            words = [p for p in product(*(range(x+1) for x in bounds)) if sum(p) == m]
            maxima = tuple(map(max, zip(*(pair_capacities(n, p) for p in words))))
            rows[n, m] = maxima
    expected = {
        ((1, 2, 3, 4), 10): (12, 12, 6, 0), ((1, 3, 3, 3), 10): (12, 6, 6, 6),
        ((2, 2, 2, 4), 10): (6, 12, 12, 0), ((2, 2, 3, 3), 10): (8, 14, 8, 8),
        ((2, 2, 3, 3), 11): (6, 10, 4, 4), ((3, 0, 1, 6), 10): (12, 12, 6, 0),
        ((3, 0, 3, 4), 10): (0, 24, 6, 0), ((3, 1, 2, 4), 10): (0, 18, 12, 0),
        ((3, 1, 3, 3), 10): (0, 18, 6, 6), ((3, 2, 2, 3), 10): (0, 12, 12, 6)}
    assert rows == expected
    return rows


def verify_even():
    third = Fraction(22**2, 18)
    for (n, m), ds in colour_table().items():
        for k, maxpairs in enumerate(ds):
            for t in range(40, 4*m+1):
                for d in range(min(maxpairs, 58-t)+1):
                    for o in range(59-t-d):
                        e = 58-t-d-o
                        if (e and k) or d+2*t < 95:
                            continue
                        # For t=39 the universal first energy bound receives an extra 18.
                        lower = 368+third+4*e+8*o-8*t+10*(d % 2)
                        if lower > 84:
                            continue
                        assert ((m,t,d,e,o) in [(10,40,18,0,0),(11,44,10,0,4)] and k==1)
                        l0, l1, l2, l3, h = caps(n)
                        assert l0 + l1 + l2 + l3 + h == m
                        nu = list(map(lambda v: 4*v, [3*l0+l1+l2+l3, 2*l1+h, 2*l2+h, 2*l3+h]))
                        nu[0] += o
                        # k=1 means 01 and 23 pairs. Account separately for both capacities.
                        cap01 = 2*(n[0]*n[1] - 2*l1)
                        cap23 = 2*(n[2]*n[3] - h)
                        trials = []
                        for d01 in range(cap01+1):
                            d23 = d-d01
                            if 0 <= d23 <= cap23:
                                q = [s + v for s, v in zip(nu, [d01,d01,d23,d23])]
                                cost = 0
                                for size, total in zip(n,q):
                                    if size == 0:
                                        assert total == 0
                                        continue
                                    lo, rem = divmod(total, size)
                                    cost += (size-rem)*(lo-14)**2 + rem*(lo-13)**2
                                trials.append((q,cost))
                        expected = {(3,0,3,4):([48,0,42,66],38), (3,1,2,4):([46,14,28,68],42),
                                    (3,1,3,3):([42,18,48,48],40), (2,2,3,3):([28,36,46,46],44)}
                        assert trials == [expected[n]]
                        Q = o + 2*d + 3*t
                        E = 600+4*o-4*d-12*t + trials[0][1] + Fraction((162-Q)**2, 8)
                        assert E > 84
                        print('even', n, 'm', m, '(e,o,d,t)', (e,o,d,t), 'sums/variance', trials[0], 'E >=', E)


def check_latin():
    # Enumerate labelled order-three squares and all eight-cell deletions.
    symbols = tuple(range(3))
    squares = []
    for rows in product(permutations(symbols), repeat=3):
        if any(len({r[j] for r in rows}) != 3 for j in symbols):
            continue
        cells = {(i,j,rows[i][j]) for i,j in product(symbols, repeat=2)}
        squares.append(cells)
        def disjoint(x,y):
            return all(a != b for a,b in zip(x,y))
        for v in cells:
            neighbours = [w for w in cells if disjoint(v,w)]
            assert len(neighbours) == 2 and disjoint(*neighbours)
    assert len(squares) == 12
    partial = {frozenset(square - {v}) for square in squares for v in square}
    assert len(partial) == 108
    for pos in combinations(product(symbols, repeat=2),8):
        # The two completely filled rows permit only 3! times 2 possibilities.
        complete = sorted(set(symbols)-{next(i for i in symbols if sum(x==i for x,y in pos)==2)})
        missingrow = next(i for i in symbols if i not in complete)
        for row1,row2 in product(permutations(symbols), repeat=2):
            if any(x==y for x,y in zip(row1,row2)):continue
            rows = dict(zip(complete,(row1,row2)))
            cells = {(i,j,rows[i][j] if i!=missingrow else 3-row1[j]-row2[j]) for i,j in pos}
            assert frozenset(cells) in partial
    print('Latin squares/deletions', len(squares),len(partial))


def check_petersen():
    bounds = tuple(packing(10-a)+min(comb(a,2),(10-a)*(a//2)) for a in range(1,11))
    assert bounds==(12,9,10,10,12,13,10,8,4,0)
    points = set(range(6))
    pairs = tuple(combinations(sorted(points),2))
    mats = [frozenset(m) for m in combinations(pairs,3) if {x for pair in m for x in pair}==points]
    assert len(mats)==15
    factors=0
    for fac in combinations(mats,5):
        if len(set.union(*(set(m) for m in fac))) !=15:continue
        factors+=1
        unchosen = set(mats)-set(fac)
        omitted = {m:frozenset(i for i,p in enumerate(fac) if p.isdisjoint(m)) for m in unchosen}
        assert len(set(omitted.values())) == 10 and all(len(x)==2 for x in omitted.values())
        for x,y in combinations(unchosen,2):
            assert (not x.isdisjoint(y)) == omitted[x].isdisjoint(omitted[y])
        for deleted in fac:
            edges = {edge:frozenset(m for m in unchosen if edge in m) for edge in set(pairs)-deleted}
            assert all(len(v)==2 for v in edges.values())
            covered = set(list(unchosen)[:1])
            for _ in unchosen:
                for ends in edges.values():
                    if covered & ends:covered.update(ends)
            assert len(covered)==10
            for x,y in combinations(edges,2):
                samecol=any(x in p and y in p for p in fac)
                assert bool(edges[x] & edges[y]) == (set(x).isdisjoint(y) and not samecol)
    assert factors==6
    print('odd branch check weights', bounds, 'labelled one-factorizations', factors)


if __name__=='__main__':
    verify_even()
    check_latin()
    check_petersen()
