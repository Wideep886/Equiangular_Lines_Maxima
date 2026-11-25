# Replay of the fourth revision

Use Python 3.9 or later with NumPy installed for the full-shell isometries. All scripts use integer or rational arithmetic. Run the commands from the directory containing the main manuscript and its checkpoint archive.

```sh
unzip -q EL18_58_Search_Checkpoint.zip -d checkpoint
python3 -B audit/revision4/replay.py checkpoint/EL18_58_search
```

The driver checks each program's output against the accompanying expected output. The default transcript is in `replay.out`. To include the earlier independent weighted clique calculation on the same 21 root-rank-sixteen domains, use

```sh
python3 -B audit/revision4/replay.py --legacy-rank16 checkpoint/EL18_58_search
```

| File | What it checks |
| --- | --- |
| `check_rank10.py` | Rank-ten color, multiplicity and variance tables used in the proof. |
| `rebuild_support16.py` | From the incoming rank-sixteen frontier, rebuilds the 21 domains and compares all 2,307 labelled word occurrences, in order, with the archive. |
| `support16.py` | Counts component shapes and the residual positions used by the rank-sixteen support argument. |
| `octacode_arcs.py`, `mark5_table.py` | Obtain the same tournament for row 5 of the eight-copy A3 host from the 64 shortened codewords and from the norm-four shell. |
| `middle_minors.py` | Rebuilds the chamber blocks of A3 rows 6 and 10. Checks the five integer minors and the projection condition for every qualifying pair. |
| `verify_a4_26.py` | Reconstructs the 17 blocks in A4 row 26 and checks the ranks of the specified block centers. |
| `verify_a2_3.py` | Reconstructs the 16 blocks in A2 row 3 and the colored-edge table for its ten large blocks. |
| `../section_isometries.py` | Checks the 15 signed shells, comparing points and perpendicular roots in each of the four indicated isometries. |

The proof of each counting or rank bound is in the manuscript. The older terminal calculations remain available in the checkpoint.
