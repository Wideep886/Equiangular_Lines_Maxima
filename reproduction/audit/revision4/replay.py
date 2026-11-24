"""Replay the fourth revision against the extracted EL18_58_search checkpoint."""
import argparse
import os
from pathlib import Path
import subprocess
import sys

HERE = Path(__file__).resolve().parent

def check(name, script, expected, *args):
    result = subprocess.run([sys.executable, '-B', str(script), *map(str, args)],
                            check=True, stdout=subprocess.PIPE, text=True)
    if result.stdout.splitlines() != expected.read_text(encoding='utf-8').splitlines():
        raise SystemExit(f'Output mismatch: {name} ({expected.name})')
    print('PASS', name, flush=True)

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('checkpoint', type=Path, help='extracted EL18_58_search directory')
    parser.add_argument('--legacy-rank16', action='store_true',
                        help='also compare with the earlier weighted clique calculation')
    args = parser.parse_args()
    if not __debug__ or os.environ.get('PYTHONOPTIMIZE'):
        raise SystemExit('Run with Python assertions enabled.')
    root = args.checkpoint.resolve()
    joint = root / 'joint23'
    profiles = root / 'root_profiles'
    a3 = joint / 'eleven' / 'a3'
    for file in (a3 / 'octacode.py', joint / 'a5' / 'all.out',
                 profiles / 'frontier_current.json', profiles / 'orbit16.json'):
        if not file.is_file():
            raise SystemExit(f'Missing checkpoint file: {file}')
    tests = [
        ('rank 10 hand tables', 'check_rank10.py', 'check_rank10.out', ()),
        ('rank 16 domains from the incoming frontier', 'rebuild_support16.py', 'rebuild_support16.out', (profiles,)),
        ('rank 16 shapes and residual positions', 'support16.py', 'support16.out', (profiles / 'orbit16.json',)),
        ('A3 row 5 shortened-code tournament', 'octacode_arcs.py', 'mark5_table.out', ()),
        ('A3 row 5 tournament from the full shell', 'mark5_table.py', 'mark5_table.out', (a3,)),
        ('A3 rows 6 and 10 five integer minors', 'middle_minors.py', 'middle_minors.out', (a3,)),
        ('A4 row 26 product blocks', 'verify_a4_26.py', 'verify_a4_26.out', ()),
        ('A2 row 3 colored edges', 'verify_a2_3.py', 'verify_a2_3.out', ()),
        ('15 signed shells and 4 root-preserving isometries', '../section_isometries.py',
         '../section_isometries.out', (joint,)),
    ]
    if args.legacy_rank16:
        tests.append(('earlier rank 16 weighted clique bounds', '../root16_terminal.py',
                      '../root16_terminal.out', (root,)))
    for label, script, output, arguments in tests:
        check(label, HERE / script, HERE / output, *arguments)

if __name__ == '__main__':
    main()
