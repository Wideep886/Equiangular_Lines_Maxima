#!/usr/bin/env python3
"""Reproduce all finite arithmetic inputs in the manuscript."""
from pathlib import Path
import subprocess, sys, concurrent.futures
if not __debug__: raise RuntimeError('Run with assertions enabled')
BASE=Path(__file__).resolve().parent
CHECKS=(('scalar','verify_scalar.py'),('small','verify_small.py'),
        ('golay','golay_design_check.py'),('leech_incidence','leech_incidence_check.py'))
def run(item):
    name,script=item
    p=subprocess.run([sys.executable,'-u',script],cwd=BASE,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
    (BASE/f'{name}_verification.log').write_text(p.stdout+f'EXIT_STATUS {p.returncode}\n')
    if p.returncode: raise RuntimeError((name,p.returncode,p.stdout[-2000:]))
    print('PASS',name,flush=True)
    return name,p.stdout
if __name__=='__main__':
    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
        out=dict(pool.map(run,CHECKS))
    assert 'SCALAR VERIFIED 20 types; 1885 scalar patterns; 1464 envelope passes; 421 single-local refinements; maximum 87' in out['scalar']
    assert 'SMALL VERIFIED 3 scalar patterns 3840 scalar orbits 22 local alternatives 24 bound 90' in out['small']
    assert 'SMALL VERIFIED 4 scalar patterns 1256 scalar orbits 10 local alternatives 13 bound 88' in out['small']
    assert 'S(5,8,24) verified' in out['golay'] and 'Normal-column partition bounds verified:' in out['golay']
    assert 'LEECH INCIDENCE TABLE AND MINIMUM LIFTS VERIFIED; LOWER BOUND 90' in out['leech_incidence']
    print('ALL REQUIRED FINITE VERIFICATIONS PASSED.',flush=True)
    print('The analytic arguments and classical theorems in the manuscript are also essential.',flush=True)
