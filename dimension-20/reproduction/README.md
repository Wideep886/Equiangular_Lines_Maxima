# Reproduce the finite calculations in dimension twenty

This supplement accompanies the dimension-twenty manuscript. The [replay driver](audit/referee/replay.py) runs the four supplied verifiers in a temporary directory and compares their outputs with the five [reference transcripts](audit/referee/expected/). The [exact glue-model data](audit/glue_models.json) are included with the code.

## Requirements and commands

Use Python 3.10 or newer with assertions enabled. The checks use the Python standard library and run offline. From the repository root, run

```sh
cd dimension-20/reproduction
shasum -a 256 -c MANIFEST.sha256
unset PYTHONOPTIMIZE
python3 -B audit/referee/replay.py
```

The replay reports `PASS dimension-twenty replay (five transcripts matched)` after it checks the four detailed transcripts and the master summary. The four stage completion messages can appear in any order. The replay leaves the supplied inputs and reference transcripts in place.

## Evidence in the manuscript

| Paper | Computation and supporting files |
| --- | --- |
| Section 3, the twenty general rooted types | The [scalar verifier](audit/verify_scalar.py) checks 1,885 scalar profiles, comprising 1,464 envelope checks and 421 single-local refinements. The maximum refined bound is 87. The [Cartan utilities](audit/cartan_utils.py), [local orbit routines](audit/orbit_models.py), and [orbit capacities](audit/energy_bound.py) use integer or exact rational arithmetic. |
| Section 3, the two small rooted types | The [small-model verifier](audit/verify_small.py) covers 3,840 and 1,256 profiles for `A2^12` and `A3^8`, grouped into 22 and 10 scalar orbits with 24 and 13 local alternatives, respectively. The bounds for the contradiction argument are 90 and 88. |
| Section 4, the binary Golay model | The [Golay verifier](audit/golay_design_check.py) checks the punctured and extended binary Golay codes, the 759-octad Steiner system `S(5,8,24)`, and ten partition bounds. |
| Sections 5–6, Leech-lattice incidence and the construction | The [incidence verifier](audit/leech_incidence_check.py) checks 16 heptad-incidence rows, the five minimum-lift types, and the lower-bound count $168-2\cdot 42+6=90$. |

The [proof dependencies](AUDIT.md) identify the mathematical arguments and the classical lattice results used in addition to these finite checks. The [SHA-256 manifest](MANIFEST.sha256) covers the verifier sources, the exact input, the expected transcripts, and this guide.
