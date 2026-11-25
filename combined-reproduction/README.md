# Reproduce the finite checks for the combined manuscript

[`replay.py`](replay.py) runs the existing verification drivers for dimensions 18, 19 and 20. It checks the SHA-256 manifests and uses temporary working directories for extracted certificates and generated outputs. A mismatch in the inputs or expected outputs causes a nonzero exit status.

## Run all three dimensions

Use Python 3.11 or newer, Node.js (`node` on `PATH`) and a C++17 compiler (`c++` on `PATH`). Run these commands from the repository root.

```sh
python3 -m venv .venv
. .venv/bin/activate
python -m pip install -r dimension-19/reproduction/requirements.txt
unset PYTHONOPTIMIZE
python -B combined-reproduction/replay.py
```

The dependency list installs NumPy, NetworkX and SymPy. After installation, the checks use the files in this repository and run offline. Allow at least 200 MB of temporary disk space for the dimension-18 checkpoint and its working files.

## Run one dimension

Pass `--dimension 18`, `--dimension 19` or `--dimension 20`. Dimensions 19 and 20 also support Python 3.10. Dimension 18 needs only NumPy from the Python dependency list, plus Node.js and the compiler. Dimension 19 uses all three Python packages and dimension 20 uses only the Python standard library.

| Dimension | Combined manuscript | Checks | Details |
| --- | --- | --- | --- |
| 18 | Section 2.8 | The exact 57-line witness, all three exclusion stages and the second terminal reconstruction | [Reproduction guide](../dimension-18/reproduction/README.md) |
| 19 | Section 3.7 | All three stages, including the regenerated 72-vector witness and the root and lattice certificates | [Reproduction guide](../dimension-19/reproduction/README.md) |
| 20 | Section 4.7 | Four finite checks and comparison with the five archived transcripts | [Reproduction guide](../dimension-20/reproduction/README.md) |

### Dimension nineteen certificate dependencies

For Theorem 3.1, the saturation construction in Section 3.2 associates a root type and a nonnegative vector of label counts with any hypothetical 73-line system having a dual unit vector. Sections 3.3 and 3.4 explain why those counts would satisfy the rows reconstructed by [`root19_mixed_deficit.py`](../dimension-19/reproduction/audit/root19_mixed_deficit.py).

| Verifier step | Mathematical source in the combined manuscript |
| --- | --- |
| `cap`, `patterns` and `allowed` | The minuscule weights, residual squares, multiplicities and fiber bounds in Section 3.3 determine the candidate variables. The orbit extrema and block-trace lower bound in Section 3.4.1 give the zero-residual filter. |
| `setup` counting and packing rows | The total count is 73. The orbit extrema, support intersections and singleton residual bounds in Section 3.3 give the slot and support-packing bounds. Section 3.4.1 adds the subspace bounds and the averaged projection rows. |
| `setup` trace and branch rows | The frame identity and the block-trace bounds in Sections 3.1 and 3.4.1 constrain the component and complementary traces. The one-dimensional residual argument gives the square-class row for a branched certificate. Section 3.4.2 gives the additional rows for the specified types. |
| `deficit` and `tangent` | The edge and fiber counts of Section 3.3 bound the deficit of each label. The root-system deficit lemma in Section 3.2 gives the total deficit bound and the tangent rows in Section 3.4.1. |
| `verify_dual` | Section 3.4.3 uses the stored nonnegative rational multipliers to combine the reconstructed rows into an infeasibility certificate. |

For Theorem 3.2, the [small-model verifier](../dimension-19/reproduction/audit/niemeier_19_certificate.py) enumerates 88 marking orbits in `A2^12` and computes bounds for 81. It identifies the seven remaining orbits by the marked shape of the norm-six vector and reports them as `ANALYTIC`. The ternary Golay argument in Section 3.6.2 excludes these seven using common height and rank-two projection support.

### Nested Leech sections

The separate [Leech section checker](check_leech_sections.py) verifies the construction in Section 4.6. From the repository root, run with Python 3.10 or newer and assertions enabled.

```sh
python3 -B combined-reproduction/check_leech_sections.py
```

The checker constructs the 276 axes, verifies their pairings, and uses the three normals in the manuscript to recover nested sections of dimensions 23, 22, 21, and 20 with 276, 176, 126, and 90 axes, respectively. It checks the ranks over the rationals. The generalized Hamming-weight argument in Sections 4.5 and 4.6 proves that these counts are maximal within the 276-axis configuration.
