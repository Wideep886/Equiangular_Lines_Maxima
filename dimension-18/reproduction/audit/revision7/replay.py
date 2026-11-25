"""Replay the single normal-plane test on the entire terminal frontier."""
import argparse
import os
from pathlib import Path
import re
import subprocess
import sys
from tempfile import TemporaryDirectory

HERE = Path(__file__).resolve().parent


def checked(command, **kwargs):
    return subprocess.run(command, check=True, capture_output=True, **kwargs).stdout


def transcript(name, data):
    expected = (HERE / name).read_text().splitlines()
    if data.splitlines() != expected:
        raise ValueError(f'Different result in {name}')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('checkpoint', type=Path, help='extracted EL18_58_search')
    parser.add_argument('--audit', type=Path, default=HERE.parent,
                        help='parent of revision2 and revision6')
    parser.add_argument('--second', action='store_true',
                        help='also run a Walsh implementation on the same coordinates')
    args = parser.parse_args()
    if not __debug__ or os.environ.get('PYTHONOPTIMIZE'):
        raise ValueError('Keep Python assertions enabled.')
    joint = args.checkpoint.resolve() / 'joint23'
    sources = joint / 'normal_sieve'
    cpp = [os.environ.get('CXX','c++'), '-std=c++17', '-O3']
    with TemporaryDirectory() as tmp:
        root = Path(tmp)
        inputs = root / 'inputs'
        data = checked([sys.executable, '-B', HERE/'terminal.py', joint,
                        args.audit.resolve(), inputs], text=True)
        transcript('coordinates.out', data)
        rows = sorted(inputs.glob('*.in'))
        assert len(rows)==15
        print('PASS all 15 coordinate inputs regenerated from the 216 list covers', flush=True)
        checker = root / 'checker'
        checked([*cpp, sources/'audit'/'independent.cpp', '-o', checker])
        counts = checked([checker, *rows], text=True)
        transcript('normal.out', counts)
        summaries = [re.search(r' pairs=(\d+) twoPass=(\d+) triples=(\d+)$', s)
                     for s in counts.splitlines()]
        assert len(summaries)==15 and all(summaries)
        assert all(int(s[2])==int(s[3])==0 for s in summaries)
        print('PASS all 15 direct normal tests; '
              f'{sum(int(s[1]) for s in summaries)} ordered-first plane tests; 0 pass',
              flush=True)
        if args.second:
            check_rows = checked([sys.executable, '-B', HERE/'independent_hosts.py', joint,
                                  args.audit.resolve(), inputs], text=True)
            transcript('independent_hosts.out', check_rows)
            print('PASS all cells, edges and binary kernels rebuilt using an integer row basis', flush=True)
            second = root / 'second'
            checked([*cpp, sources/'sieve.cpp', '-o', second])
            lines = []
            for row in rows:
                output=checked([second], input=row.read_bytes()).decode().strip()
                assert 'no eligible kernel' in output and 'pair passes 0' in output
                lines.append(f'{row.name} {output}')
            transcript('walsh.out','\n'.join(lines))
            print('PASS all 15 Walsh counts on the same inputs', flush=True)


if __name__=='__main__':
    main()
