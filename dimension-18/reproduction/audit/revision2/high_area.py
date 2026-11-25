"""Export the Cartan-derived residue-group arithmetic and compare all saved area sectors.

Usage: python3 high_area.py EL18_58_search/path COMPILED_HIGH_AREA_EXECUTABLE
"""
from pathlib import Path
from fractions import Fraction
from functools import lru_cache
from math import lcm
import json, subprocess, sys, tempfile
from high_domains import classes, configurations, valid

coefs=((0,2),(0,3),(1,1),(1,-1),(1,2),(1,-2),(2,1),(2,-1),(3,1),(3,-1))

@lru_cache(None)
def component(tag,S):
    weights=classes(tag);items=sorted(weights)
    representatives={}
    for i,a in enumerate(items):representatives.setdefault(weights[a][0],i)
    lines=[]
    for a in items:
        rows=[]
        for b in items:
            energies= [int(S*weights[b][0]), weights[b][1]]
            for m,n in coefs:
                index=tuple((m*x+n*y)%1 for x,y in zip(a,b))
                energies.append(int(S*weights[index][0]))
            rows.append(' '.join(map(str,energies)))
        lines.append(' '.join(sorted(set(rows))))
    return representatives, [len(x.split())//12 for x in lines], lines

def main():
    p=Path(sys.argv[1])/'root_profiles/plane16.json'
    decoded=json.loads(p.read_text()); ncases=0
    with tempfile.TemporaryFile(mode='w+t') as stream:
        print(len(decoded),file=stream)
        for name,record in sorted(decoded.items()):
            ts=name.split('.'); S=lcm(*(v[0].denominator for tag in ts for v in classes(tag).values()))
            entries={tag:component(tag,S) for tag in set(ts)}
            print(len(ts),S,file=stream)
            for tag in ts:
                _,counts,rows=entries[tag];print(len(rows),file=stream)
                for k,line in zip(counts,rows):print(k,line,file=stream)
            allprofiles=sorted((prof,x) for prof,x in configurations(ts) if x[0]==5 or x[1]<=3)
            cases=[(i,prof,x) for i,(prof,x) in enumerate(allprofiles) if
                   (h:=5-x[0])>0 and valid(x[2]+4*h,0) and valid(x[3]+9*h,1)]
            assert len(cases)==len(record);print(len(cases),file=stream)
            for i,prof,x in cases:
                a=int(S*(5-x[0]));indices=[entries[tag][0][z] for tag,z in zip(ts,prof)]
                expected=sorted(map(int,record[str(i)]))
                print(a,*indices,len(expected),*expected,file=stream);ncases+=1
        stream.flush();stream.seek(0)
        subprocess.run([sys.argv[2]],stdin=stream,check=True)
    print(ncases,'anchor records compared')

if __name__=='__main__':main()
