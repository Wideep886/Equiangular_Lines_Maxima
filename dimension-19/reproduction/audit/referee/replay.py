#!/usr/bin/env python3
"""Replay the finite checks accompanying the dimension-nineteen manuscript."""

import argparse
import difflib
import hashlib
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
from zipfile import ZipFile


HERE = Path(__file__).resolve().parent
AUDIT = HERE.parent
ARCHIVE = AUDIT.parent / "EL19_Certificates.zip"
EXPECTED = HERE / "expected"
MIXED = "root19_all_mixed_certificates.jsonl"
HOMOGENEOUS = "root19_homogeneous_DE_recheck.jsonl"
ASCHE = "Asche72_integer_vectors.txt"
CERTIFICATES = {MIXED, HOMOGENEOUS, ASCHE}


def unpack_certificates(directory):
    with ZipFile(ARCHIVE) as archive:
        members = archive.namelist()
        if len(members) != len(CERTIFICATES) + 1 or set(members) != CERTIFICATES | {"MANIFEST.sha256"}:
            raise ValueError("Unexpected entries in the certificate archive.")
        manifest = {}
        for line in archive.read("MANIFEST.sha256").decode("ascii").splitlines():
            digest, name = line.split(maxsplit=1)
            if name in manifest:
                raise ValueError(f"Duplicate certificate digest for {name}.")
            manifest[name] = digest
        if set(manifest) != CERTIFICATES:
            raise ValueError("The certificate manifest does not cover the expected inputs.")
        for name, digest in manifest.items():
            content = archive.read(name)
            if hashlib.sha256(content).hexdigest() != digest:
                raise ValueError(f"Certificate digest mismatch for {name}.")
            (directory / name).write_bytes(content)


def verify(directory, script, expected, *options):
    print(f"Running {expected}", flush=True)
    result = subprocess.run(
        [sys.executable, "-B", script, *options],
        cwd=directory,
        stdout=subprocess.PIPE,
        text=True,
    )
    if result.returncode:
        print(result.stdout, end="", file=sys.stderr)
        result.check_returncode()
    reference = (EXPECTED / f"{expected}.out").read_text().splitlines(keepends=True)
    observed = result.stdout.splitlines(keepends=True)
    if observed != reference:
        difference = difflib.unified_diff(reference, observed, f"{expected}.out", "actual", n=2)
        sys.stderr.writelines(difference)
        raise RuntimeError(f"Output mismatch in {script}.")
    print(f"PASS {expected}", flush=True)


def stage_a(directory):
    verify(directory, "check_coverage.py", "coverage")
    verify(directory, "root19_mixed_deficit.py", "mixed", "--verify", MIXED, "--complete")
    verify(directory, "root19_mixed_deficit.py", "homogeneous", "--verify", HOMOGENEOUS)
    verify(directory, "root_rank19_pure_A1_certificate.py", "pure-a1")


def stage_b(directory):
    verify(directory, "niemeier_19_certificate.py", "small-niemeier")
    verify(directory, "verify_niemeier20.py", "niemeier-twenty")


def stage_c(directory):
    witness = directory / "recomputed-Asche72.txt"
    verify(directory, "audit_elementary.py", "elementary", "--source", str(directory), "--output", str(witness))
    if witness.read_bytes() != (directory / ASCHE).read_bytes():
        raise RuntimeError("The regenerated Asche witness does not match the 72 archived rows.")
    print("PASS 72-vector witness", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--only", choices=("A", "B", "C"), help="Run a single stage.")
    args = parser.parse_args()
    if not __debug__ or os.environ.get("PYTHONOPTIMIZE"):
        parser.error("Enable Python assertions and clear PYTHONOPTIMIZE.")
    if sys.version_info < (3, 10):
        parser.error("Use Python 3.10 or newer.")
    with tempfile.TemporaryDirectory(prefix="el19-replay-") as scratch:
        work = Path(scratch)
        unpack_certificates(work)
        for source in AUDIT.glob("*.py"):
            shutil.copy2(source, work / source.name)
        stages = {"A": stage_a, "B": stage_b, "C": stage_c}
        for name, run in stages.items():
            if args.only is None or args.only == name:
                print(f"Stage {name}", flush=True)
                run(work)
    print("PASS dimension-nineteen replay", flush=True)


if __name__ == "__main__":
    main()
