#!/usr/bin/env python3
"""Replay the finite checks supporting the combined manuscript."""

import argparse
import hashlib
import os
from pathlib import Path
import shutil
import subprocess
import sys
from tempfile import TemporaryDirectory
from zipfile import BadZipFile, ZipFile


REPOSITORY = Path(__file__).resolve().parents[1]
DIMENSIONS = (18, 19, 20)
REFEREE_REPLAY = Path("audit/referee/replay.py")
CHECKPOINT_ARCHIVE = Path("EL18_58_Search_Checkpoint.zip")
WITNESS = Path("audit/centered_57.py")


def verify_manifest(supplement, required=(), snapshot=None):
    manifest = supplement / "MANIFEST.sha256"
    checked = set()
    for line in manifest.read_text(encoding="ascii").splitlines():
        expected, name = line.split("  ", 1)
        relative = Path(name)
        source = supplement / relative
        if (
            relative.is_absolute()
            or ".." in relative.parts
            or relative in checked
            or not source.resolve().is_relative_to(supplement)
        ):
            raise ValueError(f"Invalid or duplicate manifest path {name!r}.")

        if snapshot is not None:
            target = snapshot / relative
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(source, target)
            source = target

        digest = hashlib.sha256()
        with source.open("rb") as file:
            for block in iter(lambda: file.read(1024 * 1024), b""):
                digest.update(block)
        if digest.hexdigest() != expected:
            raise ValueError(f"SHA-256 mismatch for {name}.")
        checked.add(relative)

    missing = set(required) - checked
    if missing:
        names = ", ".join(path.as_posix() for path in sorted(missing))
        raise ValueError(f"Files missing from the manifest include {names}.")


def replay_dimension_18(supplement):
    with TemporaryDirectory(prefix="el18-combined-replay-") as scratch:
        work = Path(scratch)
        witness_log = WITNESS.with_suffix(".out")
        referee_log = REFEREE_REPLAY.with_suffix(".out")
        verify_manifest(
            supplement,
            required=(CHECKPOINT_ARCHIVE, WITNESS, witness_log, REFEREE_REPLAY, referee_log),
            snapshot=work,
        )
        checkpoint = work / "checkpoint"
        with ZipFile(work / CHECKPOINT_ARCHIVE) as archive:
            archive.extractall(checkpoint)

        for script, transcript, arguments in (
            (WITNESS, witness_log, ()),
            (REFEREE_REPLAY, referee_log, (checkpoint / "EL18_58_search", "--second")),
        ):
            result = subprocess.run(
                [sys.executable, "-B", script, *arguments],
                cwd=work,
                stdout=subprocess.PIPE,
                text=True,
                check=False,
            )
            print(result.stdout, end="", flush=True)
            result.check_returncode()
            if result.stdout != (work / transcript).read_text(encoding="utf-8"):
                raise ValueError(f"Output mismatch for {transcript.as_posix()}.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--dimension", type=int, choices=DIMENSIONS, help="Replay only this dimension."
    )
    args = parser.parse_args()
    dimensions = DIMENSIONS if args.dimension is None else (args.dimension,)
    minimum_python = (3, 11) if 18 in dimensions else (3, 10)
    if sys.version_info < minimum_python:
        parser.error(f"Use Python {minimum_python[0]}.{minimum_python[1]} or newer.")
    if not __debug__ or os.environ.get("PYTHONOPTIMIZE"):
        parser.error("Enable Python assertions and clear PYTHONOPTIMIZE.")

    for dimension in dimensions:
        supplement = REPOSITORY / f"dimension-{dimension}" / "reproduction"
        print(f"Replaying dimension {dimension}", flush=True)
        try:
            if dimension == 18:
                replay_dimension_18(supplement)
            else:
                verify_manifest(supplement, required=(REFEREE_REPLAY,))
                subprocess.run([sys.executable, "-B", REFEREE_REPLAY], cwd=supplement, check=True)
        except (BadZipFile, OSError, ValueError, subprocess.CalledProcessError) as error:
            sys.exit(f"Dimension {dimension} replay failed. {error}")
        print(f"PASS combined reproduction for dimension {dimension}", flush=True)


if __name__ == "__main__":
    main()
