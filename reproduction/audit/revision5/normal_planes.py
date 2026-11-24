"""Count symmetry classes of necessary binary normal planes in the two A2 shells.
Arguments: directory of revision2 scripts; extracted joint23 directory.
The candidate words and the chamber graphs are regenerated from the marked shell.
"""
from collections import Counter
from functools import cache
from pathlib import Path
import sys

import numpy as np


def closure(identity, generators, compose):
    group = [identity]
    seen = {identity}
    for current in group:
        for generator in generators:
            product = compose(current, generator)
            if product not in seen:
                group.append(product)
                seen.add(product)
    return group


def local_slots(signed_permutation, w):
    """Choose the local permutations that fix w and the chosen root chamber."""
    slots = []
    for block, entry in enumerate(signed_permutation):
        source = abs(entry) - 1
        sign = 1 if entry > 0 else -1
        old = tuple(sign * x for x in w[source])
        target = w[block]
        if sorted(old) != sorted(target):
            return None
        for k, value in enumerate(target):
            choices = [j for j, x in enumerate(old) if x == value]
            if sign < 0:
                choices.reverse()
            slots.append(3 * source + choices[target[:k].count(value)])
    return slots


def apply_linear(word, column_images):
    result = 0
    while word:
        bit = word & -word
        result ^= column_images[bit.bit_length()-1]
        word ^= bit
    return result


def actions(w, code, transforms):
    def compose(g, h):
        return tuple(g[abs(t)-1] * (1 if t > 0 else -1) for t in h)

    identity = tuple(range(1, 12))
    generators = [tuple((i+1)*s for i,s in zip(p,sign)) for p,sign,_ in transforms]
    group = closure(identity, generators, compose)
    stabilizer = sorted(g for g in group if local_slots(g,w) is not None)
    known = {identity}
    generators = []
    for g in stabilizer:
        if g not in known:
            generators.append(g)
            known = set(closure(identity, generators, compose))
    assert known == set(stabilizer)

    def bits(point):
        return sum((a & 1) << j for j,a in enumerate(point))

    def shorten(point):
        return bits([point[3*i+j] for i in range(11) for j in (0,1)])

    h = shorten(sum(w, ()))
    pivot = h.bit_length() - 1

    def quotient(point):
        if (point >> pivot) & 1:
            point ^= h
        return ((point >> (pivot+1)) << pivot) | (point % (1 << pivot))

    primal = []
    for g in generators:
        assert {tuple(((1 if t>0 else -1)*row[abs(t)-1]) % 3 for t in g)
                for row in code} == set(code)
        slots = local_slots(g,w)
        columns = []
        for j in range(21):
            word = 1 << (j + (j >= pivot))
            point = []
            for i in range(11):
                x, y = (word >> (2*i)) & 1, (word >> (2*i+1)) & 1
                point.extend((x, y, x ^ y))
            columns.append(quotient(shorten(tuple(point[s] for s in slots))))
        primal.append(columns)

    identity = tuple(1 << j for j in range(21))
    dual_generators = [tuple(sum(((e & col).bit_count() % 2) << j
                                 for j,col in enumerate(cols)) for e in identity)
                       for cols in primal]
    dual = closure(identity, dual_generators,
                   lambda g,h: tuple(apply_linear(x,h) for x in g))
    assert len(stabilizer) % len(dual) == 0
    return primal, dual, len(stabilizer)


def read_cells(model):
    lines = model.splitlines()
    d, n = map(int, lines[0].split())
    assert d == 21 and len(lines) == 1+2*n
    cells = []
    for a,b in zip(lines[1::2],lines[2::2]):
        size,*words = map(int,a.split())
        adj = list(map(int,b.split()))
        assert size == len(words) == len(adj) and size <= 16
        cells.append((words,adj))
    return cells


def check_actions(cells, maps):
    """Every action must permute all vertices, whole cells, and every within-cell edge."""
    index = {z:(i,j) for i,(words,_) in enumerate(cells) for j,z in enumerate(words)}
    assert len(index) == sum(len(ws) for ws,_ in cells)
    for cols in maps:
        images = {z:apply_linear(z,cols) for z in index}
        assert set(images.values()) == set(index)
        for words,adj in cells:
            dest = [index[images[z]] for z in words]
            i = dest[0][0]
            assert {a for a,b in dest} == {i}
            assert {b for a,b in dest} == set(range(len(cells[i][0])))
            for j,(_,v) in enumerate(dest):
                for k,(_,u) in enumerate(dest):
                    assert (adj[j] >> k & 1) == (cells[i][1][v] >> u & 1)


def cell_table(words, adjacency):
    @cache
    def clique(mask):
        if not mask:
            return 0
        v = mask & -mask
        return max(clique(mask^v), 1 + clique(mask & adjacency[v.bit_length()-1]))

    size = len(words)
    lookup = np.fromiter((clique(i) for i in range(1<<size)), dtype=np.uint8)
    columns = [sum((z >> j & 1) << i for i,z in enumerate(words)) for j in range(21)]
    return lookup, columns, (1<<size)-1


def hyperplanes(tables):
    count = np.zeros(1<<21,dtype=np.int16)
    bound = np.zeros_like(count)
    pop = np.fromiter((m.bit_count() for m in range(1<<16)),dtype=np.uint8)
    # All entries are <= 210. Every local mask has at most 16 bits.
    for lookup,columns,full in tables:
        mask = np.empty(1<<21,dtype=np.uint16)
        mask[0] = full
        for j,col in enumerate(columns):
            k = 1 << j
            mask[k:2*k] = mask[:k] ^ col
        count += pop[mask]
        bound += lookup[mask]
    return count,bound


def pair_bounds(u, others, tables):
    capacities = np.zeros(len(others),dtype=np.int16)
    for lookup,columns,full in tables:
        a = full
        b = np.full(len(others),full,dtype=np.uint16)
        for j,col in enumerate(columns):
            if u >> j & 1:
                a ^= col
            b ^= (others >> j & 1).astype(np.uint16) * col
        capacities += lookup[b & a]
    return capacities.tolist()


def verify(case, w, code, transforms, model):
    cells = read_cells(model)
    primal,group,order = actions(w,code,transforms)
    check_actions(cells,primal)
    tables = [cell_table(*cell) for cell in cells]
    n,bound = hyperplanes(tables)
    good = np.flatnonzero((bound >= 58) & (np.arange(len(n)) != 0))
    high = set(map(int,good[7*n[good] >= 3*n[0]+232]))
    first = []
    while high:
        u = min(high)
        first.append(u)
        orbit = {apply_linear(u,g) for g in group}
        assert orbit <= high
        high.difference_update(orbit)

    seen_planes = set()
    capacities = Counter()
    for u in first:
        v = good[(good < (good^u)) & (bound[good^u] >= 58) & (good != u)]
        section_size = (n[u]+n[v]+n[u^v]-n[0])//2
        v = v[(3*section_size >= n[u]+116) & (n[v] <= n[u]) & (n[u^v] <= n[u])]
        representatives = []
        for x in map(int,v):
            plane = tuple(sorted((u,x,u^x)))
            if plane in seen_planes:
                continue
            representatives.append(x)
            for g in group:
                a,b = apply_linear(u,g),apply_linear(x,g)
                seen_planes.add(tuple(sorted((a,b,a^b))))
        if representatives:
            bounds = pair_bounds(u,np.array(representatives,dtype=np.int64),tables)
            capacities.update(bounds)
            if case == 2:
                for x, capacity in zip(representatives, bounds):
                    size = (int(n[u]) + int(n[x]) + int(n[u^x]) - int(n[0])) // 2
                    print('plane', tuple(sorted((u,x,u^x))), 'eligible axes', size,
                          'clique bound', capacity)
    total_high = sum(7*n[int(u)] >= 3*n[0]+232 for u in good)
    print('index',case,'group',order,'axes',int(n[0]),'cell bound',int(bound[0]))
    print('eligible nonzero normals',len(good),'first normals',int(total_high),
          'first representatives',len(first))
    print('distinct plane orbits',sum(capacities.values()),'max cell bound',max(capacities),
          'capacities',dict(sorted(capacities.items())))
    assert max(capacities) < 58


if __name__ == '__main__':
    audit, joint = map(Path,sys.argv[1:])
    sys.path.insert(0,str(audit.resolve()))
    from exceptional_shells import cases, codes
    from exceptional_orbits import make_transforms
    from odd_sections import make

    code = codes(3)
    transforms = make_transforms(3,code,None)
    for mark,case,w,n in cases(joint,3):
        if case in (2,4):
            verify(case,w,code,transforms,make(3,mark,w,code,n))
