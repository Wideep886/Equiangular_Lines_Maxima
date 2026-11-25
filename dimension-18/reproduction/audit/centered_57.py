"""Reconstruct the 45+12 line witness from Section 7. Uses only the standard library."""
from fractions import Fraction
from itertools import combinations, permutations


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def centered(indices):
    return tuple(10 * (j in indices) - len(indices) for j in range(10))


def schur_inertia(gram):
    matrix = [[Fraction(v) for v in row] for row in gram]
    rank = 0
    while matrix:
        pivot_index = next((i for i, row in enumerate(matrix) if row[i]), None)
        if pivot_index is None:
            assert all(x == 0 for row in matrix for x in row)
            return rank, len(matrix)
        matrix[0], matrix[pivot_index] = matrix[pivot_index], matrix[0]
        for row in matrix:
            row[0], row[pivot_index] = row[pivot_index], row[0]
        pivot = matrix[0][0]
        assert pivot > 0
        matrix = [[row[j] - row[0] * matrix[0][j] / pivot
                   for j in range(1, len(matrix))] for row in matrix[1:]]
        rank += 1
    return rank, 0


def tetrahedral_quintets():
    """Return the Lin--Munemasa--Taniguchi--Yoshino blocks in the earlier labeling."""
    vertices = (5, 6, 7, 9)
    # Each triple of successors is ordered lexicographically; the first vertex is fixed.
    cycle_labels = (8, 1, 4, 0, 2, 3)
    cycle_arcs = {}
    for label, tail in zip(cycle_labels, permutations(vertices[1:])):
        cycle = (vertices[0],) + tail
        cycle_arcs[label] = set(zip(cycle, cycle[1:] + cycle[:1]))
    all_points = set(range(10))
    assert set(vertices).isdisjoint(cycle_labels)
    assert set(vertices) | set(cycle_labels) == all_points
    representatives = []
    for arc in permutations(vertices, 2):
        i, j = arc
        block = (set(vertices) - {i}) | {c for c, arcs in cycle_arcs.items() if (i, j) in arcs}
        # Replacing a block by its complement reverses the centered vector.
        representatives.append(tuple(sorted(block if 0 in block else all_points - block)))
    return sorted(representatives)


def main():
    quintets = list(map(set, tetrahedral_quintets()))
    previous = "01239 01278 01346 01458 02345 02368 02489 02567 03478 03579 04679 05689".split()
    assert ["".join(map(str, sorted(q))) for q in quintets] == previous
    assert all(len(a) == 5 for a in quintets)
    assert all(len(a & b) in (2, 3) for a, b in combinations(quintets, 2))
    rows = [centered(a) + (0,) * 10 for a in quintets]
    matchings = []
    for i in range(9):
        matching = [{i, 9}] + [{(i - j) % 9, (i + j) % 9} for j in range(1, 5)]
        assert len(set.union(*matching)) == 10
        matchings.extend(frozenset(edge) for edge in matching)
        rows.extend(centered({i}) + centered(edge) for edge in matching)
    assert len(set(matchings)) == 45
    gram = []
    for i, a in enumerate(rows):
        scaled = [dot(a, b) for b in rows]
        assert scaled[i] == 250
        assert all(abs(x) == 50 for j, x in enumerate(scaled) if j != i)
        gram.append([x // 50 for x in scaled])
    assert len(gram) == 57 and schur_inertia(gram) == (18, 39)
    print("The directed K4 construction reproduces all twelve earlier quintets in order, up to orientation.")
    print("45 edges, 12 quintets, 66 pair intersections and 1596 off-diagonal pairs verified.")
    print("Scaled coordinates YY^T=50G. Rational Schur inertia (positive,negative,zero)=(18,0,39).")


if __name__ == "__main__":
    main()
