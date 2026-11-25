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

| Dimension | Checks | Details |
| --- | --- | --- |
| 18 | The exact 57-line witness, all three exclusion stages and the second terminal reconstruction | [Reproduction guide](../dimension-18/reproduction/README.md) |
| 19 | All three stages, including the regenerated 72-vector witness and the root and lattice certificates | [Reproduction guide](../dimension-19/reproduction/README.md) |
| 20 | Four finite checks and comparison with the five archived transcripts | [Reproduction guide](../dimension-20/reproduction/README.md) |
