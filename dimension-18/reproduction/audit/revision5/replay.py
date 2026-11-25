"""Replay the November 24 reductions using the extracted EL18_58_search."""
import argparse
import os
from pathlib import Path
import subprocess
import sys

HERE = Path(__file__).resolve().parent

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('checkpoint', type=Path)
    root = parser.parse_args().checkpoint.resolve()
    if not __debug__ or os.environ.get('PYTHONOPTIMIZE'):
        raise SystemExit('Run with Python assertions enabled.')
    codes = root / 'integer_frame' / 'codes.txt'
    joint = root / 'joint23'
    a3 = joint / 'eleven' / 'a3'
    tables = HERE.parent / 'revision2'
    if not codes.is_file() or not (a3 / 'octacode.py').is_file():
        raise SystemExit('The argument must be the extracted EL18_58_search directory.')
    checks = [
        ('integer_cuts', 'full-rank integral codes and five residual orbit lists', [codes]),
        ('mark0', 'A3 index 0 cells', []),
        ('mark1', 'A3 index 1 cells', []),
        ('mark14_table', 'A3 index 14 centers and minor', [a3]),
        ('mark15_table', 'A3 index 15 centers and minor', [a3]),
        ('full_cells', 'independent rational rank test at A3 indices 1 and 15', [root]),
        ('nonad_glue_check', 'five Paley rows in affine quadratic coordinates', []),
        ('nonad_basis', 'second affine-code check on all 243 words and pencil signs', []),
        ('nonad_cells', 'A2 index 4 affine cells, steps and centers from the ternary code', []),
        ('nonad_check', 'all affine centers and steps in the original coordinates', []),
        ('nonad_partitions', 'A2 index 4 affine partition and rational relation tests', []),
        ('nonad_quotient', 'independent rational quotient test of the affine table', []),
        ('normal_planes', 'binary plane orbits and section capacities', [tables, joint]),
        ('hexad16', 'a geometric full-cell anchor for the last section', []),
        ('normal_anchor_check', 'first normal orbits from the coordinate actions', [tables, joint]),
    ]
    for script, label, arguments in checks:
        command = [sys.executable, '-B', str(HERE / (script + '.py')), *map(str, arguments)]
        data = subprocess.run(command, check=True, stdout=subprocess.PIPE, text=True).stdout
        if data.splitlines() != (HERE / (script + '.out')).read_text(encoding='utf-8').splitlines():
            raise SystemExit(f'Output mismatch for {script}.')
        print('PASS', label, flush=True)

if __name__ == '__main__':
    main()
