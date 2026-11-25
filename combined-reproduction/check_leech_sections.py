#!/usr/bin/env python3
"""Check the sizes and rational ranks of the nested Leech sections."""

from fractions import Fraction
from itertools import combinations
from pathlib import Path
from runpy import run_path
import sys

from replay import REPOSITORY, verify_manifest


REPRODUCTION = REPOSITORY / "dimension-20" / "reproduction"
GOLAY_CHECK = Path("audit/golay_design_check.py")
EXTERIOR = range(1, 24)


def dot(left, right):
    return sum(x * y for x, y in zip(left, right, strict=True))


def rational_rank(rows):
    basis = {}
    for vector in rows:
        row = tuple(map(Fraction, vector))
        for column in range(len(row)):
            leading = row[column]
            if leading == 0:
                continue
            if column not in basis:
                basis[column] = tuple(entry / leading for entry in row)
                break
            row = tuple(x - leading * y for x, y in zip(row, basis[column], strict=True))
    return len(basis)


def scaled_normal(evaluations):
    """Return 36 times the normal with the given point evaluations."""
    total = sum(evaluations)
    return (total, *(18 * value - total for value in evaluations))


def main():
    if sys.version_info < (3, 10):
        sys.exit("Use Python 3.10 or newer.")
    if not __debug__:
        sys.exit("Enable Python assertions.")

    verify_manifest(REPRODUCTION, required=(GOLAY_CHECK,))
    # The existing Golay checker stores the verified extended octad supports in O.
    octads = [frozenset(support) for support in run_path(str(REPRODUCTION / GOLAY_CHECK))["O"]]
    heptads = [support - {0} for support in octads if 0 in support]
    octad = next(support for support in octads if 0 not in support)
    a, b = sorted(octad)[:2]

    # Use integer representatives 4 p_j and 4 q_H of the axes in the manuscript.
    points = [(3, *(8 * (i == j) - 1 for i in EXTERIOR)) for j in EXTERIOR]
    heptad_axes = [tuple(4 * (i in support) - 1 for i in range(24)) for support in heptads]
    axes = points + heptad_axes
    host_normal = (5, *(1 for _ in EXTERIOR))
    assert all(dot(axis, host_normal) == 0 and dot(axis, axis) == 80 for axis in axes)
    assert all(abs(dot(left, right)) == 16 for left, right in combinations(axes, 2))

    normals = [
        scaled_normal([1 - 2 * (i == a) for i in EXTERIOR]),
        scaled_normal([-int(i not in octad) for i in EXTERIOR]),
        scaled_normal([1 - 2 * (i == b) for i in EXTERIOR]),
    ]
    assert rational_rank(normals) == 3
    for normal in normals:
        assert dot(normal, host_normal) == 0
        assert all(dot(normal, axis) % 144 == 0 for axis in axes)

    ranks = [rational_rank(axes)]
    sizes = [len(axes)]
    section = axes
    for normal in normals:
        section = [axis for axis in section if dot(axis, normal) == 0]
        ranks.append(rational_rank(section))
        sizes.append(len(section))

    assert ranks == [23, 22, 21, 20], ranks
    assert sizes == [276, 176, 126, 90], sizes
    print("Leech section dimensions:", ", ".join(map(str, ranks)))
    print("Leech section sizes:", ", ".join(map(str, sizes)))


if __name__ == "__main__":
    main()
