"""Verify the four exceptional chain tables against the extracted checkpoint."""
import argparse
import gzip
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

HERE = Path(__file__).resolve().parent

def check_script(name, arguments):
    result = subprocess.run([sys.executable, '-B', str(HERE / (name + '.py')), *map(str, arguments)],
                            check=True, capture_output=True, text=True)
    if result.stdout.splitlines() != (HERE / (name + '.out')).read_text().splitlines():
        raise ValueError(f'{name} output differs from the archived transcript')

def records(path):
    with gzip.open(path, 'rt', encoding='ascii') as stream:
        return json.load(stream)

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('checkpoint', type=Path, help='extracted EL18_58_search directory')
    parser.add_argument('--regenerate', action='store_true', help='also reconstruct the links by matching')
    args = parser.parse_args()
    if not __debug__ or os.environ.get('PYTHONOPTIMIZE'):
        raise SystemExit('Use Python with assertions enabled.')
    certificate = HERE / 'exceptional_chains.json.gz'
    paths = [args.checkpoint.resolve() / 'joint23', HERE.parent / 'revision2']
    check_script('verify_chains', [*paths, certificate])
    print('PASS full coordinate shells and all 74308 links in 216 records', flush=True)
    given = records(certificate)
    counts = [(r['scale'], r['mark'], r['index'], r['points'], r['points']-len(r['links'])) for r in given]
    expected = (HERE.parent / 'revision3' / 'exceptional-216.txt').read_text().splitlines()
    assert sorted(counts) == sorted(tuple(map(int, line.split())) for line in expected)
    print('PASS all 216 list counts agree entry by entry with the independent middle-layer table', flush=True)
    if args.regenerate:
        with tempfile.TemporaryDirectory() as temp:
            rebuilt = Path(temp) / 'chains.json.gz'
            check_script('generate_chains', [*paths, rebuilt])
            assert records(rebuilt) == given
        print('PASS regenerated all 216 partitions from the marked models')

if __name__ == '__main__':
    main()
