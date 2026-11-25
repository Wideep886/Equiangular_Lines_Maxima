# Verifying the increasing lists

For each marking, exceptional_chains.json.gz contains the integer scale, marked factor, marking index, total number of candidates, and pairs of zero-based indices in lexicographic order on the marked shell. A pair i,j is a directed link. Different records may describe overlapping families of line configurations.

The verifier enumerates every norm-four point in the marked-factor model, filters by its inner product with the specified w, and checks each link against the original integer coordinates. The difference is a positive multiple of a root perpendicular to the mark and w. Each endpoint has at most one outgoing and one incoming link. Positivity rules out cycles, so subtracting the number of links from the number of candidates counts all lists, including isolated points.

Run from the package root, with the original checkpoint extracted:

    python3 -B audit/revision6/replay.py EL18_58_search
    python3 -B audit/revision6/replay.py EL18_58_search --regenerate

The first command validates the supplied lists and requires NumPy. The second additionally recreates the links by matching and compares the lists. Both methods use the marking list and local coordinate and glue definitions from revision 2. The verifier enumerates the full marked shells before imposing the condition on w. The counts are compared with revision 3, which calculates middle layers independently of the links.

| Lattice | Records | Candidate occurrences | Links | Records below 116 |
| --- | ---: | ---: | ---: | ---: |
| A₅⁴D₄ | 166 | 72,534 | 58,656 | 165 |
| A₄⁶ | 28 | 11,616 | 9,012 | 26 |
| A₃⁸ | 17 | 7,290 | 5,398 | 9 |
| A₂¹² | 5 | 1,824 | 1,242 | 1 |

The separate historical certificates for the eighteen other lattices allow a nonzero sum of positive roots for a link. Their verifiers are in joint23/{thirteen,four-chains,marked-chains} in the original checkpoint.
