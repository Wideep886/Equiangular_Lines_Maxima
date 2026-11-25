"""Reconstruct all root-factor isometry orbits of the five labelled support domains.

Run labelled_supports.py ROOT OUT before this file.  Standard library Python 3.11+.
"""
from itertools import permutations, product
from pathlib import Path
from collections import Counter
from hashlib import sha256
import sys,json

def label(x):
    return sum(t*3**j for j,t in enumerate(x))

def words(cols):
    t=len(cols[0])
    return {tuple(sum(x*y for x,y in zip(v,m))%3 for v in cols)
            for m in product(range(3),repeat=t)}

def group_actions(a,b,cols,keys):
    # Enumerate the images of each dual basis vector in the full ternary dual code.
    # For each coordinate permutation derive *every* compatible signing. This
    # produces all signed permutations preserving this code, including the kernel.
    allowed=words(cols)
    images=list(product(sorted(allowed),repeat=len(cols[0])))
    left=list(permutations(range(a)))
    dims=(2,)*a+(3,)*b
    decoded={k:tuple(k//3**i%3 for i in range(a+b)) for k in keys}
    powers=tuple(3**i for i in range(a+b))
    full=set()
    for p in permutations(range(b)):
        source=tuple(cols[j] for j in p)
        for rows in images:
            signs=[]
            for j,col in enumerate(source):
                options=[s for s in (1,-1) if all((s*v-r[j])%3==0 for v,r in zip(col,rows))]
                if not options:break
                signs.append(options)
            else:
                for eps in product(*signs):
                    for q in left:
                        ix=q+tuple(a+j for j in p)
                        sign=(1,)*a+eps
                        perm=[]
                        for k in keys:
                            y=tuple((s*decoded[k][j])%d for s,j,d in zip(sign,ix,dims))
                            ny=tuple((-z)%d for z,d in zip(y,dims))
                            perm.append(label(min(y,ny)))
                        assert set(perm)==set(keys)
                        full.add(tuple(perm))
    return [dict(zip(keys,g)) for g in full]

def classify(raw,a,b,cols,output=None):
    data=set(tuple(map(int,r.split())) for r in raw.splitlines())
    block=[]
    if not data:
        if output is not None:output.write_text('')
        return [],0
    keys=sorted(set().union(*map(set,data)))
    maps=group_actions(a,b,cols,keys)
    rep=[]
    while data:
        seed=min(data)
        orbit={tuple(sorted(g[i] for i in seed)) for g in maps}
        assert orbit<=data and seed==min(orbit)
        if output is not None:
            block.extend((s,seed) for s in orbit)
        data.difference_update(orbit)
        rep.append((seed,len(orbit)))
    if output is not None:
        output.write_text(''.join(' '.join(map(str,k))+'|'+' '.join(map(str,v))+'\n' for k,v in sorted(block)))
    return rep,len(maps)

def main():
    root=Path(sys.argv[1])/'root_profiles'
    p=Path(sys.argv[2])
    outdir=Path(sys.argv[3]) if len(sys.argv)>3 else None
    if outdir is not None:outdir.mkdir(parents=True,exist_ok=True)
    codes=json.loads((root/'full_dual.json').read_text())
    old=json.loads((root/'high_normals.json').read_text())
    for r in (13,14):
        for typ,records in json.loads((root/f'fixedaxes_{r}.json').read_text()).items():
            types=typ.split('.');a,b=types.count('A1'),types.count('A2')
            for key,rec in records.items():
                if not rec['live']:continue
                cols=codes[f'3:{b}:{18-r}'][int(key)]['cols']
                filename=f'{a}_{b}_{key}.txt'
                reps,ng=classify((p/filename).read_text(),a,b,cols,outdir/filename if outdir is not None else None)
                canon=''.join(' '.join(map(str,x))+' #'+str(w)+'\n' for x,w in reps)
                prev=old[typ+':'+key]['stats']
                assert len(reps)==prev.get('orbits',0) and sum(w for x,w in reps)==prev.get('labelled',0)
                print(a,b,key,len(reps),ng,dict(sorted(Counter(w for x,w in reps).items())),sha256(canon.encode()).hexdigest(),flush=True)
if __name__=='__main__':main()
