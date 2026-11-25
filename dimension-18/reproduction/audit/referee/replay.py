"""Replay the three finite stages of the article. Run on an extracted checkpoint."""
import argparse
import os
from pathlib import Path
import subprocess
import sys

HERE = Path(__file__).resolve().parent
AUDIT = HERE.parent
PY = [sys.executable, '-B']

def check(cmd, folder=None, transcript=None):
    run = subprocess.run([str(c) for c in cmd], cwd=folder, check=True, stdout=subprocess.PIPE, text=True)
    if transcript is not None and run.stdout.splitlines()!=transcript.read_text().splitlines():
        raise RuntimeError(f'Output mismatch for {cmd[2]}')
    return run.stdout

def roots(top,args):
    joint=top/'joint23'
    command=PY+[joint/'replay.py','--group','roots','--manifest']
    if args.node:
        command+=['--node',args.node]
    check(command, top)
    for name in ('midrank_bridge','check_frontiers'):
        check(PY+[AUDIT/'revision2'/(name+'.py'),top],transcript=AUDIT/'revision2'/(name+'.out'))
    print('PASS A: root certificates, rank transfers and terminal frontiers',flush=True)
    profiles=top/'root_profiles'
    supplements=[
        ('revision4','check_rank10',[]),
        ('revision4','rebuild_support16',[profiles]),
        ('revision4','support16',[profiles/'orbit16.json']),
        ('revision5','integer_cuts',[top/'integer_frame/codes.txt']),
    ]
    for revision,name,inputs in supplements:
        folder=AUDIT/revision
        check(PY+[folder/(name+'.py'),*inputs],transcript=folder/(name+'.out'))
    print('PASS A: orthogonal rank 10 tables, 21 rank 16 domains, 562 integral-code orbits',flush=True)

def heights(top,args):
    joint=top/'joint23'
    folder=joint/'thirteen'
    check(PY+[folder/'verify.py'],folder,folder/'verify.out')
    print('PASS B: direct chain bounds for 13 ordinary lattice types',flush=True)
    for name in ('four-chains','marked-chains'):
        folder=joint/name
        check(PY+[folder/'verify.py'],folder,folder/'verify.out')
    print('PASS B: list bounds for the remaining five ordinary lattice types',flush=True)
    check(PY+[HERE/'models.py',joint,AUDIT/'revision2'],transcript=HERE/'models.out')
    check(PY+[AUDIT/'revision2/exceptional_orbits.py',joint],transcript=AUDIT/'revision2/exceptional_orbits.out')
    check(PY+[AUDIT/'revision6/verify_chains.py',joint,AUDIT/'revision2',AUDIT/'revision6/exceptional_chains.json.gz'],
          transcript=AUDIT/'revision6/verify_chains.out')
    print('PASS B: models, marked-factor coverage, 216 marked shells, 15 terminal records',flush=True)

def terminals(top,args):
    command=PY+[AUDIT/'revision7/replay.py',top]
    if args.second:
        command.append('--second')
    report=check(command)
    lines=(AUDIT/'revision7/replay.out').read_text().splitlines()
    assert report.splitlines()==(lines if args.second else lines[:2])
    print(report,end='',flush=True)
    print('PASS C: every required plane fails on all 15 original marked records',flush=True)

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('checkpoint',type=Path)
    parser.add_argument('--only',choices=['A','B','C'],help='replay one proof stage')
    parser.add_argument('--node',help='path to Node.js, if it is not on PATH')
    parser.add_argument('--second',action='store_true',help='also check the terminal graph inputs and Walsh counts independently')
    args=parser.parse_args()
    if not __debug__ or os.environ.get('PYTHONOPTIMIZE'):
        parser.error('Enable Python assertions and clear PYTHONOPTIMIZE.')
    top=args.checkpoint.resolve()
    if not (top/'joint23/replay.py').is_file():
        parser.error('Extract the checkpoint and pass its EL18_58_search directory.')
    for stage,run in [('A',roots),('B',heights),('C',terminals)]:
        if args.only in (None,stage):
            run(top,args)
if __name__=='__main__':
    main()
