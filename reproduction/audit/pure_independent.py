"""Check the last fifteen pure rank-twelve obstructions by integer determinants.

Run with the root of the extracted EL18_58_search directory as the argument.
No sign normalization, floating point, or third-party software is used.
"""
from itertools import combinations, product
from pathlib import Path
import sys


def determinant(matrix):
    a = [row[:] for row in matrix]
    n = len(a)
    sign, previous = 1, 1
    for k in range(n - 1):
        pivot_row = next((i for i in range(k, n) if a[i][k]), None)
        if pivot_row is None:
            return 0
        if k != pivot_row:
            a[k], a[pivot_row] = a[pivot_row], a[k]
            sign *= -1
        pivot = a[k][k]
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                value = a[i][j] * pivot - a[i][k] * a[k][j]
                quotient, remainder = divmod(value, previous)
                assert remainder == 0
                a[i][j] = quotient
        previous = pivot
    return sign * a[-1][-1]


def main(root):
    directory = root / 'joint23' / 'pure12-independent'
    templates = [list(map(int, row.split())) for row in
                 (directory / 'supp.txt').read_text().splitlines()]
    logs = (directory / 'sign.out').read_text().splitlines()[:-1]
    assert len(templates) == len(logs) == 15
    total = 0
    for row_number, (template, line) in enumerate(zip(templates, logs), start=1):
        index, six, last, required, _ = line.split()
        assert int(index) == row_number
        length, *members = template
        assert len(members) == length
        masks = list(map(int, six.split(','))) + [int(last)]
        for mask in masks:
            members.remove(mask)
        assert len(masks) == 7 and 5 - bin(masks[-1]).count('1') == int(required)
        positions = [(i, j) for i, j in combinations(range(7), 2)
                     if not masks[i] & masks[j]]
        determinants = set()
        for signs in product((-1, 1), repeat=len(positions)):
            gram = [[0] * 7 for mask in masks]
            for i, mask in enumerate(masks):
                gram[i][i] = 5 - bin(mask).count('1')
            for (i, j), sign in zip(positions, signs):
                gram[i][j] = gram[j][i] = sign
            determinants.add(determinant(gram))
        assert 0 not in determinants
        total += 2 ** len(positions)
        print(row_number, 2 ** len(positions), sorted(determinants))
    print('templates', len(logs), 'signed seven-row determinants', total, 'zero', 0)


if __name__ == '__main__':
    main(Path(sys.argv[1]))
