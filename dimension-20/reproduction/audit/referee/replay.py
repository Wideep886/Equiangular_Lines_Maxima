#!/usr/bin/env python3
"""Replay the finite checks accompanying the dimension-twenty manuscript."""

import difflib
import hashlib
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile


HERE = Path(__file__).resolve().parent
AUDIT = HERE.parent
SUPPLEMENT = AUDIT.parent
EXPECTED = HERE / "expected"
TRANSCRIPTS = (
    "scalar_verification.log",
    "small_verification.log",
    "golay_verification.log",
    "leech_incidence_verification.log",
)


def verify_manifest(required):
    checked = set()
    for line in (SUPPLEMENT / "MANIFEST.sha256").read_text(encoding="ascii").splitlines():
        digest, name = line.split("  ", 1)
        relative = Path(name)
        if relative.is_absolute() or ".." in relative.parts or relative in checked:
            raise ValueError(f"Invalid or duplicate manifest path: {name}")
        actual = hashlib.sha256((SUPPLEMENT / relative).read_bytes()).hexdigest()
        if actual != digest:
            raise ValueError(f"SHA-256 mismatch: {name}")
        checked.add(relative)

    missing = {path.relative_to(SUPPLEMENT) for path in required} - checked
    if missing:
        raise ValueError("Files missing from the manifest: " + ", ".join(map(str, sorted(missing))))


def compare(name, actual):
    expected = (EXPECTED / name).read_text(encoding="utf-8")
    if actual == expected:
        return
    difference = difflib.unified_diff(
        expected.splitlines(keepends=True),
        actual.splitlines(keepends=True),
        fromfile=f"expected/{name}",
        tofile=f"replayed/{name}",
    )
    sys.stderr.writelines(difference)
    raise RuntimeError(f"Output mismatch in {name}.")


def stable_master(output):
    # Independent stages run concurrently, so completion order can vary.
    lines = output.splitlines(keepends=True)
    return "".join(sorted(lines[: len(TRANSCRIPTS)]) + lines[len(TRANSCRIPTS) :])


def main():
    if sys.version_info < (3, 10):
        sys.exit("Use Python 3.10 or newer.")
    if not __debug__ or os.environ.get("PYTHONOPTIMIZE"):
        sys.exit("Enable Python assertions and clear PYTHONOPTIMIZE.")

    inputs = (*sorted(AUDIT.glob("*.py")), AUDIT / "glue_models.json")
    transcripts = (EXPECTED / name for name in (*TRANSCRIPTS, "master.log"))
    try:
        verify_manifest((*inputs, *transcripts, HERE / "replay.py"))
    except (OSError, ValueError) as error:
        sys.exit(f"Cannot verify the supplement: {error}")

    print("Running the four dimension-twenty finite checks", flush=True)
    with tempfile.TemporaryDirectory(prefix="el20-replay-") as scratch:
        work = Path(scratch)
        for source in inputs:
            shutil.copy2(source, work)

        result = subprocess.run(
            [sys.executable, "-B", "verify_all.py"],
            cwd=work,
            stdout=subprocess.PIPE,
            text=True,
            encoding="utf-8",
        )
        print(result.stdout, end="", flush=True)
        result.check_returncode()

        for name in TRANSCRIPTS:
            compare(name, (work / name).read_text(encoding="utf-8"))

        expected = (EXPECTED / "master.log").read_text(encoding="utf-8")
        if stable_master(result.stdout) != stable_master(expected):
            raise RuntimeError("Output mismatch in master.log (ignoring stage completion order).")

    print("PASS dimension-twenty replay (five transcripts matched)", flush=True)


if __name__ == "__main__":
    main()
