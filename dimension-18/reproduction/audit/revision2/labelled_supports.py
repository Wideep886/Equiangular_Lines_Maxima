"""Regenerate all five labelled intermediate-rank high-support domains.
Requires the last three fixed dual-code data from the checkpoint and Python 3.11+.
Does not import the checkpoint's enumeration, fibre-weight or symmetry code.
"""
from collections import Counter
from hashlib import sha256
from itertools import product
from pathlib import Path
import json
import sys

def label(word):
    return sum(x*3**i for i,x in enumerate(word))

def digest(rows,filename=None):
    payload=''.join(' '.join(map(str,r))+'\n' for r in sorted(rows))
    if filename is not None:
        filename.write_text(payload, encoding='ascii')
    return sha256(payload.encode()).hexdigest()

def enumerate_raw(a,b,cols,alive,filename=None):
    # Supports consist of a binary coordinate in each A1 and a ternary in each A2.
    dims=(2,)*a+(3,)*b; n=len(dims); c=18-a-2*b
    shapes=sorted({tuple(sorted(x[:a])+sorted(x[a:])) for x in product(*([0,3],)*a,*([0,4],)*b)
                   if sum(x)<=15 and (sum(x)==15 or sum(bool(y) for y in x)<=3)})
    # A zero residual would give a nonzero dual word of weight three.
    assert all(sum(sum(u*v for u,v in zip(col,m)) % 3 != 0 for col in cols)>=4
               for m in product(range(3),repeat=len(cols[0])) if any(m))
    fibres=[]
    single=[0]*n; empty=0
    for x in product(*(range(d) for d in dims)):
        e=[(3 if j<a else 4)*bool(t) for j,t in enumerate(x)]
        h=15-sum(e)
        if h<=0 or sum(bool(y) for y in e)>3:continue
        if shapes.index(tuple(sorted(e[:a])+sorted(e[a:]))) not in alive:continue
        if any(sum(u*v[j] for u,v in zip(x[a:],cols))%3 for j in range(len(cols[0]))):continue
        if x>tuple((-u)%d for u,d in zip(x,dims)):continue
        I=tuple(j for j,t in enumerate(x) if t)
        if not I: empty=24*c//(25-c)
        elif len(I)==1: single[I[0]]=c+(I[0]>=a)
        else:
            # Every two sign/vertex choices must differ in >=2 coordinates.
            # Projection onto all but any one coordinate is injective.
            alph=sorted(dims[i] for i in I); cap=alph[0]*(alph[1] if len(alph)==3 else 1)
            fibres.append((sum(1<<i for i in I),I,cap,label(x)))
    # Enumerate cliques by increasing index. There is no colouring or group reduction.
    adj=[sum(1<<j for j,(t,*_) in enumerate(fibres) if j>i and (s&t).bit_count()<2)
         for i,(s,*_) in enumerate(fibres)]
    total=0; hits=[0]*n; strength=[0]*n
    cutsets=[]; chosen=[]
    def walk(P,w):
        nonlocal total
        total+=1
        slack=[max(0,cap-u) for u,cap in zip(hits,single)]
        if w+empty+sum(slack)>=58 and all(t+slack[i]>=(8 if i<a else 16) for i,t in enumerate(strength)):
            cutsets.append(tuple(sorted(chosen)))
        while P:
            u=P&-P;j=u.bit_length()-1;P^=u
            _,I,k,key=fibres[j]
            chosen.append(key)
            for t in I:hits[t]+=1; strength[t]+=k
            walk(P&adj[j],w+k)
            for t in I:hits[t]-=1; strength[t]-=k
            chosen.pop()
    walk((1<<len(fibres))-1,0)
    return len(fibres),total,len(cutsets),dict(sorted(Counter(map(len,cutsets)).items())),digest(cutsets,filename)

def main():
    root=Path(sys.argv[1])/'root_profiles'
    outdir=Path(sys.argv[2]) if len(sys.argv)>2 else None
    if outdir is not None: outdir.mkdir(parents=True, exist_ok=True)
    full=json.loads((root/'full_dual.json').read_text())
    # Fingerprints of sorted labelled tuples from the historical skeleton enumerator.
    reference=[(11424, 'bea793138deda0948d6cde3ef5fae8c7c58173c04ca4a46566a0c00fc381c5a9'), (30, 'bb0e46a8eebef8e8910dfdf269d5fd86371a1c34d2021ba5d18dde0f011460ad'), (2, '47911bf57b270914a73b24b027df4a25bae027942891239e21b207967a647ccd'), (0, 'e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855'), (128, 'd2f8ced955bf502ea8822089b6bf946cc68396662a4900cd5236d18cb52cab77')]
    result=[]
    for r in (13,14):
        for typ,records in json.loads((root/f'fixedaxes_{r}.json').read_text()).items():
            types=typ.split('.');a,b=types.count('A1'),types.count('A2')
            for key,rec in records.items():
                if not rec['live']:continue
                cols=full[f'3:{b}:{18-r}'][int(key)]['cols']
                filename=outdir/f'{a}_{b}_{key}.txt' if outdir is not None else None
                summary=enumerate_raw(a,b,cols,rec['live'],filename)
                result.append((summary[2],summary[-1]));print(a,b,key,*summary,flush=True)
    assert result==reference

if __name__=='__main__':main()
