"""Verify the A3^8 records 1 and 15 using coordinate spans over Q.

Usage: python3 -B full_cells.py PATH/TO/EL18_58_search
Only the archived octacode coordinates, shell construction and cell graphs
are shared with the normal-space verification. The rank test uses integers.
"""
import argparse
from ast import literal_eval
from collections import Counter
from functools import lru_cache, reduce
from itertools import combinations
from math import comb, gcd, lcm
from pathlib import Path
import sys


def primitive(row):
    divisor = reduce(gcd, row, 0)
    if not divisor:
        return ()
    if next(x for x in row if x) < 0:
        divisor = -divisor
    return tuple(x // divisor for x in row)


def clear_column(row, pivot_row, column):
    a = row[column]
    if not a:
        return row
    b = pivot_row[column]
    divisor = gcd(a, b)
    return primitive(tuple((b // divisor)*x - (a // divisor)*y
                           for x, y in zip(row, pivot_row)))


def row_span(rows, initial=()):
    """Unique reduced echelon form over Q, each row primitive integral."""
    basis = initial
    for row in rows:
        row = primitive(row)
        for j, pivot_row in basis:
            if not row:
                break
            row = clear_column(row, pivot_row, j)
        if row:
            k = next(j for j, x in enumerate(row) if x)
            reduced = [(j, clear_column(old, row, k)) for j, old in basis]
            basis = tuple(sorted(reduced + [(k, row)]))
    return basis


def perpendicular(basis, dimension):
    pivots = dict(basis)
    multiple = lcm(*(row[j] for j, row in basis))
    vectors = []
    for k in range(dimension):
        if k not in pivots:
            v = [0] * dimension
            v[k] = multiple
            for j, row in basis:
                v[j] = -(multiple // row[j]) * row[k]
            assert all(sum(x*y for x, y in zip(row, v)) == 0 for _, row in basis)
            vectors.append(v)
    result = row_span(vectors)
    assert len(basis) + len(result) == dimension
    return result


def shell_record(mark, h):
    from chamber import (cells, vertices, shells, graph, components, indgraph,
                         chamber_chart, cmp)
    from saturated import maximum_sets
    w = sum(h, ())
    partner = lambda x: tuple(t-a for t, a in zip(w, x))
    chart = chamber_chart(h)
    points = [x for cell in cells(h, shells(4)) for x in vertices(cell)
              if x < partner(x) and not cmp(chart(x), chart(partner(x)))]
    assert len(set(points)) == len(points)
    rows = [(chart(x), chart(partner(x))) for x in points]
    adjacency = graph(points)
    universe = (1 << len(points)) - 1
    clusters = components([universe ^ (row | 1 << j)
                           for j, row in enumerate(adjacency)])
    centered = [tuple(2*a-b for a, b in zip(x, w)) for x in points]
    columns = [j for j, _ in row_span(centered)]
    vectors = [tuple(x[j] for j in columns) for x in centered]
    dimension = len(columns)
    blocks = []
    graphs = []
    for ids in clusters:
        g = indgraph(adjacency, ids)
        for j, k in combinations(range(len(ids)), 2):
            if any(cmp(rows[ids[j]][0], t) for t in rows[ids[k]]):
                g[j] &= ~(1 << k)
                g[k] &= ~(1 << j)
        graphs.append(g)
        width, maxima = maximum_sets(g)
        spans = {row_span(vectors[ids[i]] for i in clique) for clique in maxima}
        annihilators = row_span(v for s in spans
                                for _, v in perpendicular(s, dimension))
        common = perpendicular(annihilators, dimension)
        blocks.append((width, tuple(sorted(spans)), common, ids))
    assert (len(points), len(blocks), sum(b[0] for b in blocks)) == {
        1: (280, 10, 62), 15: (214, 20, 62)}[mark]
    return points, clusters, graphs, vectors, blocks, dimension


def block_actions(h, points, clusters, graphs, actions):
    from octacode import C, reflected
    w = sum(h, ())
    reversed_h = reflected(h)
    index = {x: i for i, x in enumerate(points)}
    cell_index = {frozenset(ids): i for i, ids in enumerate(clusters)}
    actions_on_cells = set()
    for perm, signs in actions:
        if tuple((reversed_h if signs >> i & 1 else h)[j]
                 for i, j in enumerate(perm)) != h:
            continue
        assert {tuple((-word[j] if signs >> i & 1 else word[j]) % 4
                      for i, j in enumerate(perm)) for word in C} == C
        for opposite_chamber in (False, True):
            image = []
            for point in points:
                result = []
                for i, j in enumerate(perm):
                    part = point[4*j:4*j+4]
                    if signs >> i & 1:
                        part = [-a for a in (part if j == 0 else reversed(part))]
                    if i and opposite_chamber:
                        part = [a for v in sorted(set(h[i])) for a in reversed(
                            [x for x, t in zip(part, h[i]) if t == v])]
                    result.extend(part)
                axis = min(tuple(result), tuple(b-a for b, a in zip(w, result)))
                image.append(index[axis])
            assert len(set(image)) == len(points)
            cp = tuple(cell_index[frozenset(image[x] for x in ids)] for ids in clusters)
            for j, ids in enumerate(clusters):
                other = clusters[cp[j]]
                mapped = [other.index(image[x]) for x in ids]
                assert all((graphs[j][a] >> b) & 1 ==
                           (graphs[cp[j]][mapped[a]] >> mapped[b]) & 1
                           for a, b in combinations(range(len(ids)), 2))
            actions_on_cells.add(cp)
    # Check that minimization really computes the orbits of a group.
    assert all(tuple(x[i] for i in y) in actions_on_cells
               for x in actions_on_cells for y in actions_on_cells)
    return actions_on_cells


def verify(mark, h, actions):
    pts, clusters, graphs, vectors, blocks, dimension = shell_record(mark, h)
    group = block_actions(h, pts, clusters, graphs, actions)
    loss = sum(b[0] for b in blocks) - 58
    n = len(blocks)
    reps = Counter(min(tuple(sorted(g[i] for i in missing)) for g in group)
                   for missing in combinations(range(n), loss))
    counts = Counter()

    @lru_cache(None)
    def direction_bound(span):
        assert len(span) in (17, 18)
        normals = [r for _, r in perpendicular(span, dimension)]
        keys = [primitive(tuple(sum(x*y for x, y in zip(v, r)) for r in normals))
                for v in vectors]
        in_span = 0
        additional = Counter()
        for width, _, _, ids in blocks:
            sizes = Counter(keys[i] for i in ids)
            k = min(width, sizes.pop((), 0))
            in_span += k
            if len(span) == 17:
                for direction, size in sizes.items():
                    additional[direction] += min(size, width-k)
        return in_span + max(additional.values(), default=0)

    @lru_cache(None)
    def check(span, full):
        if len(span) >= 17:
            bound = direction_bound(span)
            if bound < 58:
                counts['direction', len(span)] += 1
                return bound
        if not full:
            raise ValueError('No remaining forced block at this span.')
        choice = None
        for i in full:
            children = {s for a in blocks[i][1]
                        if len(s := row_span((v for _, v in a), span)) <= 18}
            if span in children:
                return check(span, tuple(j for j in full if j != i))
            if not children:
                counts['rank', len(span)] += 1
                return 0
            if choice is None or len(children) < len(choice[1]):
                choice = i, children
        i, children = choice
        return max(check(s, tuple(j for j in full if j != i)) for s in children)

    initial = 0
    values = []
    for missing in sorted(reps):
        full = tuple(i for i in range(n) if i not in missing)
        span = row_span(v for i in full for _, v in blocks[i][2])
        if len(span) in (17, 18) and direction_bound(span) < 58:
            initial += 1
        values.append(0 if len(span) > 18 else check(span, full))
    assert max(values) < 58 and sum(reps.values()) == comb(n, loss)
    print('mark', mark, 'axes', len(pts), 'span dimension', dimension,
          'blocks', n, 'total width', 58+loss, 'omission sets', sum(reps.values()),
          'block group order', len(group), 'orbits', len(reps),
          'immediate direction cuts', initial, 'largest conditional bound', max(values))
    print('terminal states (kind, dimension, count):',
          sorted((kind, r, v) for (kind, r), v in counts.items()))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('checkpoint', type=Path)
    args = parser.parse_args()
    if not __debug__:
        raise RuntimeError('Run with assertions enabled')
    source = args.checkpoint.resolve()/'joint23'/'eleven'/'a3'
    sys.path.insert(0, str(source))
    from octacode import automorphisms
    actions = automorphisms()
    lines = (source/'octacode.out').read_text().splitlines()
    for mark in (1, 15):
        record = next(line for line in lines if line.startswith(f'{mark} h '))
        h = literal_eval(record.split(' h ')[1].split(' multiplicity')[0])
        verify(mark, h, actions)


if __name__ == '__main__':
    main()
