# Reproduce the finite calculations

This supplement matches the method and commands in the manuscript's verification paragraph. The main verifier is [`audit/referee/replay.py`](audit/referee/replay.py). The [checkpoint archive](EL18_58_Search_Checkpoint.zip) and the accompanying `audit/` files provide its inputs, integer witnesses, and expected results.

## Requirements and commands

Use Python 3.11 or newer with NumPy, Node.js, and a C++17 compiler available as `c++`. On macOS, the Xcode Command Line Tools provide the compiler. Run the commands below from the repository root.

```sh
cd dimension-18/reproduction
shasum -a 256 -c MANIFEST.sha256
python3 -m venv .venv
source .venv/bin/activate
python3 -m pip install numpy
unset PYTHONOPTIMIZE
python3 -m zipfile -e EL18_58_Search_Checkpoint.zip .
python3 -B audit/referee/replay.py EL18_58_search --second
```

Use `--only A`, `--only B`, or `--only C` to run one stage. The `--second` flag also rebuilds the terminal graphs and binary kernels and checks the section counts by a Walsh transform. The two checks use the same lattice models. Compare the output with the [full expected transcript](audit/referee/replay.out).

## Evidence in the manuscript

| Paper | Computation and supporting files |
| --- | --- |
| Section 1.1 | [`audit/centered_57.py`](audit/centered_57.py) checks the 57-line construction in exact rational arithmetic. Run `python3 -B audit/centered_57.py`. |
| Sections 2–4, stage A | The checkpoint supplies `joint23/`, `root_profiles/`, and `integer_frame/`. Checks in [`audit/revision2/`](audit/revision2/), [`audit/revision4/`](audit/revision4/), and [`audit/revision5/`](audit/revision5/) cover the root bounds, rank transfers, 21 rank-16 domains, and 562 integral-code orbits. |
| Section 6, stage B | The checkpoint rebuilds the candidate shells and verifies chain bounds for 18 ordinary lattice types. [`audit/revision6/exceptional_chains.json.gz`](audit/revision6/exceptional_chains.json.gz) contains the explicit partitions for 216 exceptional markings. |
| Section 7, stage C | [`audit/revision7/`](audit/revision7/) reconstructs the inputs for the 15 remaining marked records and checks all 472,704 ordered-first plane tests. The run reports zero passing tests. |

The [SHA-256 manifest](MANIFEST.sha256) records the digest of the certificate archive. Stage A also checks the checkpoint's internal `joint23/SHA256.json` manifest.

## Rebuild the thirteen unmarked chain certificates

The checkpoint's `joint23/thirteen/` directory covers 376 dominant choices of the squared-norm-six vector across thirteen lattice types. The largest bound in these types is 90 signed points. Starting in `reproduction/` after extraction and environment setup, verify these certificates with `python3 -B EL18_58_search/joint23/thirteen/verify.py`. To regenerate them, copy `ade.py`, `charts.py`, `sizes.py`, and `verify.py` from that directory to an empty directory, run `python3 -B charts.py` with NumPy installed, then run `python3 -B verify.py`. The generator writes its certificates alongside the scripts.
