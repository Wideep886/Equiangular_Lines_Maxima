# Exact checks for the revised proof

Run `python3 -B audit/revision5/replay.py EL18_58_search` from the top level after extracting the checkpoint. Python 3.11 or later and NumPy suffice for this entry point. All runs are deterministic. Each test compares its output with the supplied transcript and halts on any discrepancy. The other four entry points in the manuscript guide cover the antecedent reductions and the legacy proof.

## Coverage

| Calculation | Inputs | Meaning of the output |
| --- | --- | --- |
| Integral full-rank codes | The nine length-24 Type II codeword lists from the checkpoint, and the displayed coordinate symmetries | Four cases have direct bounds. For the other five codes 562 symmetry representatives cover 13,459,600 ordered choices of two disjoint triples. The old enumeration over all nine codes used 24,227,280 labelled choices. |
| Four new A₃⁸ cases | Octacode and the four marked vectors | Upper capacities are 63, 62, 61 and 62 for indices 0, 1, 14 and 15. The appendix gives the dimension contradiction under the assumption of 58. The separate full-cell calculation tests ranks over the rationals for indices 1 and 15. |
| Nine-point support | Ternary Golay coordinate generator | Five displayed generator rows yield the affine quadratic coordinates used to derive the cells without a shell search. The 408 signed points pair to 204 axes. Two coordinate reconstructions check the shapes and centers. Direct incidence counting handles all color multiplicities. Two optional partition checks reproduce the unique equality case on the 31 affine partition orbits, one using exterior coordinates and the other by rational elimination. |
| Binary planes | Shortened ternary code, actual axes and incompatibility graphs | For index 2 the three admissible plane orbits have bounds 51, 55 and 55. Index 4 has an independent binary check. The first-normal actions are rebuilt by a separate program acting directly on axes. The outputs record the symmetry group used. |

The three bounds for index 2 apply after assuming a 58-axis selection and selecting a normal plane satisfying the two averaging inequalities. They make no assertion about other 18-dimensional sections. Distinct markings or auxiliary planes need not represent distinct configurations.

## A possible further geometric reduction

`hexad16.py` verifies the center ranks and coordinate actions in the isometric A₃⁸ index 16 model of A₂¹² index 2. The proof of the main theorem uses the three certified plane bounds.
