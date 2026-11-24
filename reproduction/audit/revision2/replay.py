"""Run the second review against an extracted version-23 checkpoint."""
import argparse
import os
from pathlib import Path
import subprocess
import sys
import tempfile

folder = Path(__file__).resolve().parent


def check(name, command, cwd=None):
    result = subprocess.run(command, capture_output=True, text=True, cwd=cwd)
    if result.returncode:
        raise SystemExit(f'{name} failed ({result.returncode})\n{result.stdout}\n{result.stderr}')
    expected = (folder / (name + '.out')).read_text().splitlines()
    if result.stdout.splitlines() != expected:
        raise SystemExit(f'{name}: transcript mismatch\n{result.stdout}\n{result.stderr}')
    print(f'PASS {name}', flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('checkpoint', type=Path)
    top = parser.parse_args().checkpoint.resolve()
    joint = top / 'joint23'
    if not __debug__ or os.environ.get('PYTHONOPTIMIZE'):
        raise SystemExit('Run with Python assertions enabled.')
    if not (joint / 'replay.py').is_file():
        raise SystemExit('Expected the extracted EL18_58_search directory.')
    def python_test(name, *args):
        check(name, [sys.executable, '-B', folder / (name + '.py'), *args])
    with tempfile.TemporaryDirectory() as tmp:
        work = Path(tmp)
        flags = [os.environ.get('CXX','c++'), '-std=c++17', '-O3']
        area = work / 'area'
        deletion = work / 'deletion'
        subprocess.run([*flags, folder / 'high_area.cpp', '-o', area], check=True)
        subprocess.run([*flags, joint / 'normal_sieve/audit/independent.cpp', '-o', deletion], check=True)
        python_test('discriminant_boundary')
        python_test('pure_colour_orbits')
        python_test('midrank_bridge', top)
        python_test('high_domains', top)
        python_test('high_area', top, area)
        python_test('check_frontiers', top)
        fresh_lists = work / 'fresh-lists'
        reference_lists = work / 'original-lists'
        fresh_orbits = work / 'fresh-orbits'
        reference_orbits = work / 'original-orbits'
        python_test('labelled_supports', top, fresh_lists)
        check('labelled_baseline', [sys.executable, '-B', joint / 'dual-audit/verify_only.py',
              folder / 'labelled_baseline.py', reference_lists, reference_orbits], cwd=top / 'root_profiles')
        def list_bytes(directory):
            return {p.name: p.read_bytes() for p in directory.glob('*.txt')}
        computed, original = list_bytes(fresh_lists), list_bytes(reference_lists)
        if len(computed) != 5 or computed != original:
            raise SystemExit('Labelled support lists differ.')
        print('PASS labelled_support_sets (5 files compared byte-for-byte)', flush=True)
        python_test('labelled_orbits', top, fresh_lists, fresh_orbits)
        computed, original = list_bytes(fresh_orbits), list_bytes(reference_orbits)
        if len(computed) != 5 or computed != original:
            raise SystemExit('Labelled orbit maps differ.')
        print('PASS labelled_orbit_maps (5 files compared byte-for-byte)', flush=True)
        for name in ['exceptional_shells','exceptional_orbits','triality']:
            python_test(name, joint)
        inputs = work / 'inputs'
        python_test('odd_sections', joint, inputs)
        rows = sorted(inputs.glob('*.in'))
        if len(rows) != 7:
            raise SystemExit('Expected seven regenerated binary section inputs.')
        check('odd_deletion', [deletion, *rows])


if __name__ == '__main__':
    main()
