"""Audit flag widths without degrees, polynomial division or archived chain lists.

Usage: python3 binary_bruhat.py PATH/TO/joint23 PATH/TO/BFS_EXEC [host names]
The two programs use exact integers. Norm-shell input is regenerated using
the archived root-coordinate Weyl traversal, a separate routine from
the Dynkin-coordinate enumeration of poset_width.py.
"""
import importlib.util
from math import lcm
from pathlib import Path
from tempfile import TemporaryFile
import subprocess,sys
J=Path(sys.argv[1]).resolve()
sys.path.insert(0,str(J/'thirteen'))
from ade import factors,glue,shapes
from verify import all_for_h

types='A6_A6_A6_A6 A7_A7_D5_D5 A8_A8_A8 A9_A9_D6 A11_D7_E6 A12_A12 A15_D9 A17_E7 A24 D4_D4_D4_D4_D4_D4 D6_D6_D6_D6 D8_D8_D8 D10_E7_E7 D12_D12 D16_E8 D24 E6_E6_E6_E6 E8_E8_E8'.split()
compare_module=None
if len(sys.argv)>4 and sys.argv[3]=='--compare':
    module=importlib.util.spec_from_file_location('flag_width',Path(sys.argv[4]).resolve())
    compare_module=importlib.util.module_from_spec(module)
    module.loader.exec_module(compare_module)

def audit(name):
    fs=factors(name.split('_')); G=glue(fs); xs=shapes(fs,G,4); hs=shapes(fs,G,6)
    L=lcm(*(f.d for f in fs));offset=0;edges=[]
    for f in fs:
        edges.extend((offset+i,offset+j) for i in range(f.n) for j in range(i) if f.C[i][j]<0)
        offset+=f.n
    with TemporaryFile(mode='w+t') as stream:
        for h in hs:
            rows=all_for_h(fs,xs,h,L)
            print(len(rows),*sum(h,()),len(edges),*(i for e in edges for i in e),file=stream)
            for row in rows: print(*row,file=stream)
        stream.flush();stream.seek(0)
        proc=subprocess.run([sys.argv[2]],stdin=stream,capture_output=True,text=True,check=True)
    answers=[line.split() for line in proc.stdout.splitlines()]
    assert len(answers)==len(hs)
    if compare_module is not None:
        unmarked=tuple(() for _ in fs)
        for h,fields in zip(hs,answers):
            b=compare_module.bound(fs,h,xs,L,unmarked)
            assert b==int(fields[0])
            marked={tuple(map(int,x.split(','))) for x in fields[1:]}
            expected=set()
            if b>=116:
                shift=0
                for u,(factor,v) in enumerate(zip(fs,h)):
                    for i in range(factor.n):
                        for j in range(i):
                            if v[i] or v[j] or not factor.C[i][j]: continue
                            chosen=list(unmarked);chosen[u]=(i,j)
                            expected.add((i+shift,j+shift,compare_module.bound(fs,h,xs,L,tuple(chosen))))
                    shift+=factor.n
            assert marked==expected,(name,h,marked,expected)
    caps=[min(int(a[0]),max((int(x.split(',')[2]) for x in a[1:]),default=0)) if int(a[0])>=116 else int(a[0]) for a in answers]
    print(name,len(hs),sum(len(a)-1 for a in answers),max(caps),flush=True)
    if len(sys.argv)>3 and sys.argv[3]=='--details':
        for i,a in enumerate(answers):print(i,*a,flush=True)
if __name__=='__main__':
    start=5 if compare_module is not None else 4 if len(sys.argv)>3 and sys.argv[3]=='--details' else 3
    for name in sys.argv[start:] or types:audit(name)
