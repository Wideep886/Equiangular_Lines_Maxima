#!/usr/bin/env python3
"""Check the manuscript's fixed-angle examples using exact arithmetic.

The interval endpoints and extension inequalities are checked for the six
tabulated cases. The sufficient threshold comparisons and the two small
cosine-1/7 Gram examples are also checked. The general proofs are in the text.
"""

from fractions import Fraction
import sys

from check_leech_sections import rational_rank


SEEDS = {7: (28, 10), 23: (276, 56)}
INTERVALS = (
    (19, 7, 2, 0, 49, 51),
    (111, 23, 2, 0, 529, 539),
    (321, 7, 34, 9, 945, 952),
    (5601, 7, 562, 9, 15729, 16703),
    (5601, 23, 102, 55, 28129, 28155),
    (6001, 23, 108, 23, 29785, 30167),
)


def below_coherence_cutoff(line_count, increment, excess):
    """Compare line_count with L_k * increment without evaluating a radical."""
    difference = (excess + 1) * (excess + 2) * increment - 2 * line_count
    radical_squared = line_count**2 * excess**2 * (excess + 3)
    return difference > 0 and difference**2 > radical_squared


def cutoff_floor(increment, excess):
    # L_k has denominator greater than two, and L_k * increment is irrational.
    lower = 0
    upper = (excess + 1) * (excess + 2) * increment // 2 + 1
    while upper - lower > 1:
        middle = (lower + upper) // 2
        if below_coherence_cutoff(middle, increment, excess):
            lower = middle
        else:
            upper = middle
    return lower


def check_interval(q, excess):
    size, coefficient = SEEDS[excess]
    increment = q + 1
    block_size = 2 * ((increment + 2 * coefficient - 1) // (2 * coefficient))
    twice_degree = coefficient * block_size - increment
    assert twice_degree % 2 == 0
    degree = twice_degree // 2
    assert 0 <= degree < block_size
    base_size = size * block_size
    lower = base_size - excess
    upper = cutoff_floor(increment, excess) - excess - 1
    assert lower <= upper
    assert below_coherence_cutoff(upper + excess + 1, increment, excess)
    assert not below_coherence_cutoff(upper + excess + 2, increment, excess)

    # Evaluate the circulant adjacency matrix on the alternating character.
    offsets = {
        sign * step % block_size
        for step in range(1, degree // 2 + 1)
        for sign in (-1, 1)
    }
    if degree % 2:
        offsets.add(block_size // 2)
    assert len(offsets) == degree and 0 not in offsets
    assert increment - 2 * degree > 0
    eigenvalue = increment + 2 * sum(1 - 2 * (step % 2) for step in offsets)
    assert eigenvalue >= increment - 2 * degree
    # At most upper-lower lines are added. The constant-direction Schur
    # eigenvalue at that endpoint is checked after multiplication by eigenvalue.
    assert (q - 1) * eigenvalue > (base_size - eigenvalue) * (upper - lower)
    return block_size, degree, lower, upper


def check_threshold(excess, first_q, margin):
    size, coefficient = SEEDS[excess]
    ratio = Fraction(size, coefficient) + margin
    assert below_coherence_cutoff(ratio.numerator, ratio.denominator, excess)
    assert margin * (first_q + 1) - 2 * size > 1
    block_size, _, _, _ = check_interval(first_q, excess)
    assert block_size >= coefficient
    print(f"Sufficient threshold for k={excess}: q >= {first_q}")


def check_cosine_one_seventh(block_sizes, dimension):
    blocks = [index for index, size in enumerate(block_sizes) for _ in range(size)]
    size = len(blocks)
    # For 8I - 2 diag(J_B), block constants have eigenvalue 8 - 2|B|,
    # and the contrasts have eigenvalue 8. The final J term is positive semidefinite.
    assert all(8 - 2 * width >= 0 for width in block_sizes)
    gram = [
        [8 * (row == column) - 2 * (blocks[row] == blocks[column]) + 1 for column in range(size)]
        for row in range(size)
    ]
    assert all(gram[i][i] == 7 for i in range(size))
    assert all(abs(gram[i][j]) == 1 for i in range(size) for j in range(i))
    assert rational_rank(gram) == dimension
    assert 48 * dimension // (49 - dimension) == size  # Relative bound.
    print(f"Cosine 1/7 in dimension {dimension}: {size} lines")


def main():
    if sys.version_info < (3, 10):
        sys.exit("Use Python 3.10 or newer.")
    if not __debug__:
        sys.exit("Enable Python assertions.")

    for q, excess, block_size, degree, lower, upper in INTERVALS:
        assert check_interval(q, excess) == (block_size, degree, lower, upper)
        print(f"q={q}, k={excess}: dimensions {lower} through {upper}")
    check_threshold(7, 319, Fraction(9, 50))
    check_threshold(23, 5599, Fraction(1, 10))
    check_cosine_one_seventh((4, 4, 4), 10)
    check_cosine_one_seventh((4, 4, 4, 1), 11)


if __name__ == "__main__":
    main()
