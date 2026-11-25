# Reproduce the finite calculations in dimension nineteen

This supplement accompanies the computer-assisted proof of $N(19)=72$. The [replay driver](audit/referee/replay.py) checks the certificate archive's internal digests, reruns the checks in a temporary directory, and compares their output with the [reference transcripts](audit/referee/expected/). The [certificate archive](EL19_Certificates.zip) contains the rational certificates and the 72-vector lower-bound witness.

## Requirements and commands

Use Python 3.10 or newer with NumPy, NetworkX, and SymPy. Run the following commands from the repository root.

```sh
cd dimension-19/reproduction
python3 -m venv .venv
source .venv/bin/activate
python3 -m pip install -r requirements.txt
unset PYTHONOPTIMIZE
shasum -a 256 -c MANIFEST.sha256 && \
python3 -B audit/referee/replay.py
```

Use `--only A`, `--only B`, or `--only C` to run one stage. The replay requires assertions to be enabled. Stage A uses the Python standard library, including exact rational arithmetic for certificate verification. Stages B and C use the packages listed in [requirements.txt](requirements.txt). The rational verifier uses the stored rational multipliers.

## Evidence in the manuscript

| Paper | Computation and supporting files |
| --- | --- |
| Sections 2–4, stage A | [Coverage](audit/check_coverage.py) checks the exact disjoint decomposition of 6,709 allowed root types into 6,565 mixed types, 140 complementary types, and four exceptional types. The [rational verifier](audit/root19_mixed_deficit.py) checks 6,571 mixed certificates, including the residual branches, and 140 complementary certificates. The [orthogonal-root verifier](audit/root_rank19_pure_A1_certificate.py) checks the three pure `A1` cases. |
| Section 6, stage B | The [small-type verifier](audit/niemeier_19_certificate.py) checks the `A3^8` and `A2^12` models, with 81 computed records and seven analytic records in the ternary model. The [independent driver](audit/verify_niemeier20.py) checks twenty further lattice types, all 540 remaining marked records in the second pass, and their norm-2, norm-4, and norm-6 shell counts. The number in these script names refers to the number of lattice types in this computation. |
| Section 1 and supplementary checks, stage C | The [elementary audit](audit/audit_elementary.py) regenerates Asche's 72 vectors from the binary Golay code and checks their inner products and exact rank. It also checks minuscule data, a product-chain calculation, and the clique routine on all labeled simple graphs of orders one through five. |

The complementary transcript shows `FULL MIXED COVERAGE False` because it checks the 140 complementary types. The mixed transcript reports full coverage of the 6,565 mixed types. The separate coverage check verifies that these two groups and the four exceptional types form a disjoint, complete list.

The upper-bound proof also uses the geometric reductions and analytic cases in the manuscript. The `A2^7` argument and the seven records marked `ANALYTIC` are handled by proofs in the text.

The `shasum` command above checks the outer [SHA-256 manifest](MANIFEST.sha256), which covers the archive, this guide, the programs, the requirements file, and the expected outputs.
