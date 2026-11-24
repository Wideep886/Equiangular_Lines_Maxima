"""Recompute and compare the signed shell bounds in the structural revision."""
import argparse
import os
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory

HERE = Path(__file__).resolve().parent

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('checkpoint', type=Path)
    args = parser.parse_args()
    joint = args.checkpoint.resolve() / 'joint23'
    if not __debug__ or os.environ.get("PYTHONOPTIMIZE"):
        raise SystemExit("Python assertions must be enabled")
    if not (joint / 'replay.py').is_file():
        raise SystemExit("Supply the extracted EL18_58_search directory")

    def python(script, *args, expected=None):
        cmd = [sys.executable, '-B', HERE / script, joint, *args]
        result = subprocess.run(cmd, text=True, stdout=subprocess.PIPE, check=True)
        if expected and result.stdout.splitlines() != (HERE / expected).read_text().splitlines():
            raise SystemExit(f"Summary mismatch in {script}")

    def compare(file, expected):
        if file.read_bytes() != (HERE / expected).read_bytes():
            raise SystemExit(f"Detailed record mismatch in {expected}")

    with TemporaryDirectory() as work:
        tmp = Path(work)
        table = tmp / 'table.txt'
        for method in ('direct', 'orbit'):
            python('poset_records.py', method, table)
            compare(table, 'ordinary-4010.txt')
            print(f'PASS {method}: 1364 polarizations + 2646 marks', flush=True)
        bfs = tmp / 'bruhat'
        subprocess.run([os.environ.get('CXX', 'c++'), '-O3', '-std=c++17',
                        HERE / 'binary_bruhat.cpp', '-o', bfs], check=True)
        python('binary_bruhat.py', bfs, '--compare', HERE / 'poset_width.py',
               expected='binary_bruhat_checked.out')
        print('PASS reflection-graph distances on all 4010 records', flush=True)
        for name in ('poset_extra', 'poset_extra_orbits'):
            python(name + '.py', HERE.parent / 'revision2', table, expected='poset_extra.out')
            compare(table, 'exceptional-216.txt')
            print(f'PASS {name}: 216 marked records, 15 require terminal tests', flush=True)

if __name__ == '__main__':
    main()
