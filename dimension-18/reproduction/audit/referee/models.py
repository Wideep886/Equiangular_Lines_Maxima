"""Glue and marks before polarizations. Usage: models.py joint23 audit/revision2.
Python 3.9+ standard library. The two paths connect the independent construction to
the shortened-code and marked-polarization models consumed by the referee replay.
"""
from collections import Counter
from itertools import permutations, product
from math import prod
from pathlib import Path
import sys


def paley():
    chi = (0, 1, -1, -1, 1)
    c = [[0 if i == j else 1 if i*j == 0 else chi[(i-j) % 5]
          for j in range(6)] for i in range(6)]
    assert all(sum(c[i][k]*c[k][j] for k in range(6)) == 5*(i == j)
               for i in range(6) for j in range(6))
    return c


def generators(n):
    if n == 6:
        t = [(2, 2, 2, 0, 0), (0, 2, 4, 2, 0)]
        return t + [tuple(3*int(i in (0,j)) for i in range(4)) + (j,)
                    for j in (1,2,3)]
    if n == 5:
        return [tuple(pow(x,i,5) for x in range(5)) + (int(i == 2),) for i in range(3)]
    b = paley() if n == 3 else ((2,1,1,1),(-1,2,1,-1),(-1,-1,2,1),(-1,1,-1,2))
    return [tuple(int(i == j) for j in range(len(b))) + tuple(x % n for x in row)
            for i,row in enumerate(b)]


def subgroup(n, rows):
    a = len(rows[0]) - (n == 6)
    zero = (0,) * len(rows[0])
    seen, todo = {zero}, [zero]
    while todo:
        v = todo.pop()
        for r in rows:
            w = tuple((v[i]+r[i]) % n for i in range(a))
            if n == 6: w += (v[-1] ^ r[-1],)
            if w not in seen:
                seen.add(w)
                todo.append(w)
    return seen


def actions(n):
    # Each listed signed permutation lifts by permuting and negating whole A factors.
    # For D4 only a last-coordinate sign flip is used; it swaps the spinor labels.
    if n == 6:
        return [((1,0,2,3),(1,1,1,-1),1), ((2,3,0,1),(1,1,-1,-1),0)]
    if n == 5:
        return [((1,2,3,4,0,5),(1,)*6,0), ((5,1,3,2,4,0),(1,1,-1,-1,1,1),0)]
    if n == 4:
        return [((2,0,1,3,6,4,5,7),(1,-1,-1,1,1,-1,-1,1),0),
                ((4,0,1,2,7,3,6,5),(1,-1,-1,1,1,-1,1,-1),0)]
    pairs = [((0,2,3,4,5,1),(1,)*6), ((1,0,2,4,3,5),(1,1,1,-1,-1,1))]
    both = [(p+tuple(i+6 for i in p), s+s, 0) for p,s in pairs]
    swap = tuple(range(6,12)) + tuple(range(6))
    return both + [(swap,(1,)*6+(-1,)*6,0)]


def component_orbit(n, words):
    acts = actions(n)
    a = len(acts[0][0])
    for p,s,flip in acts:
        assert sorted(p) == list(range(a)) and all(abs(t) == 1 for t in s)
        def image(x):
            y = tuple(s[i]*x[p[i]] % n for i in range(a))
            return y + ((x[-1] ^ (flip if x[-1] >= 2 else 0)),) if n == 6 else y
        assert {image(x) for x in words} == words
    orbit = {0}
    while True:
        new = orbit | {p[i] for p,s,flip in acts for i in orbit}
        if new == orbit: break
        orbit = new
    assert orbit == set(range(a))
    return [len(orbit),1] if n == 6 else [len(orbit)]


def local_mark_count(d, kind):
    # Weyl actions on the unordered sets of the three root axes of an A2.
    # A roots are e_i-e_j. D roots are +/-e_i+/-e_j with i different from j.
    roots = {tuple((i == k) + s*(j == k) for k in range(d))
             for i,j in permutations(range(d),2) for s in ((-1,) if kind == 'A' else (-1,1))}
    roots |= {tuple(-v for v in r) for r in roots}
    def axis(r): return min(r,tuple(-v for v in r))
    marks = {frozenset((axis(r), axis(t), axis(tuple(i+j for i,j in zip(r,t)))))
             for r in roots for t in roots if sum(i*j for i,j in zip(r,t)) == -1}
    source = next(iter(marks))
    allowed_signs = [(1,)*d] if kind == 'A' else [s for s in product((1,-1),repeat=d) if prod(s) == 1]
    orbit = {frozenset(axis(tuple(s[i]*r[p[i]] for i in range(d))) for r in source)
             for p in permutations(range(d)) for s in allowed_signs}
    assert orbit == marks
    return len(orbit)


def main(joint, audit):
    sys.path.insert(0,str(audit))
    from exceptional_shells import codes
    sys.path.insert(0,str(joint/'a2'))
    from golay import generator
    for n in (6,5,4,3):
        gen = generators(n)
        words = subgroup(n,gen)
        recorded = set(codes(n))
        a = len(gen[0]) - (n == 6)
        rank = (n-1)*a + 4*(n == 6)
        discr = n**a * (4 if n == 6 else 1)
        assert rank == 24 and len(words)**2 == discr
        # Scale by n. The minimum in A_(n-1) class k is k(n-k)/n,
        # in every D4 nonzero class it is 1. Norms in each class
        # differ by even integers since the class lies in the dual.
        minima = [sum(k*(n-k) for k in w[:a]) + n*int(n == 6 and w[-1] != 0)
                  for w in words]
        assert all(m % (2*n) == 0 for m in minima)
        assert minima.count(0) == 1 and min(m for m in minima if m) == 4*n
        if n != 3:
            assert recorded == words
        else:
            # Orthogonality to the entire marked A2 block makes that coordinate zero.
            assert subgroup(3,[tuple(z%3 for z in g) for g in generator()]) == words
            short = {w[:-1] for w in words if w[-1] == 0}
            assert len(short) == 3**5 and short == recorded
            punct = {w[:-1] for w in words}
            assert len(punct) == 3**6
            assert all(sum(i*j for i,j in zip(w,g))%3 == 0 for w in short for g in gen)
            # As dimensions are 5 and 6 in length 11, punct = short-perp.
            h = Counter(sum(k != 0 for k in w) for w in words)
            sh = Counter(sum(k != 0 for k in w) for w in short)
            assert h == {0:1,6:264,9:440,12:24} and sh == {0:1,6:132,9:110}
            print('Golay weights full',sorted(h.items()),'short',sorted(sh.items()),'punctured',len(punct))
        roots = a*n*(n-1) + 24*(n == 6)
        marks = [local_mark_count(n,'A')] + ([local_mark_count(4,'D')] if n == 6 else [])
        print('scale',n,'glue',len(words),'det(root)',discr,'coset min 4, roots',roots,
              'component orbits',component_orbit(n,words),'marks per component',marks)
    print('PASS all four rank-24 even unimodular models, full root systems and A2 mark orbits')


if __name__ == '__main__':
    main(*(Path(arg).resolve() for arg in sys.argv[1:]))
