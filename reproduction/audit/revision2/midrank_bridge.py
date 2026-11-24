"""Audit the complete rank 6..15 first-anchor census from archived integer cuts.
Usage: python3 -B midrank_bridge.py EL18_58_search  (NumPy and Python 3.11+).
No module from the checkpoint is imported. The output filenames supply certificates,
not the set of ADE types, admissible label tuples, congruence rows or comparison tests.
"""
from collections import defaultdict
from fractions import Fraction
from functools import lru_cache
from itertools import product
from math import lcm
from pathlib import Path
import json
import sys

import numpy as np


@lru_cache(None)
def partitions(total, floor=1):
    if total == 0:
        return ((),)
    return tuple((n,) + p for n in range(floor, total+1)
                 for p in partitions(total-n, n))


def root_types(rank):
    types = set()
    for partition in partitions(rank):
        options = [['A'] + (['D'] if r >= 4 else [])
                   + (['E'] if r in (6, 7, 8) else []) for r in partition]
        for letters in product(*options):
            types.add(tuple(sorted(letter+str(r) for r, letter in zip(partition, letters))))
    return sorted(types)


def lengths_and_depth(name):
    n = int(name[1:])
    if name[0] == 'A':
        return [(Fraction(2*j*(n+1-j), n+1), min(j, n+1-j)) for j in range(n+1)]
    if name[0] == 'D':
        vec = (2,2)
        spin = (Fraction(n,2),n//2)
        return [(0,0),spin,vec,spin] if n%2 else [(0,0),vec,spin,spin]
    return {'E6': [(0,0),(Fraction(8,3),2),(Fraction(8,3),2)],
            'E7': [(0,0),(3,3)], 'E8': [(0,0)]}[name]


def first_anchor_check(components, cert):
    tables = [lengths_and_depth(f) for f in components]
    scale = lcm(*(w.denominator for tab in tables for w,_ in tab))
    # Keep depth <= 3 only for a positive residual; pure-root-span rows can have arbitrary depth.
    allwords = [([],0,0)]
    for tab in tables:
        allwords = [(prefix+[j],a+b,c+d) for prefix,a,c in allwords
                    for j,(weight,d) in enumerate(tab) if a+(b:=int(scale*weight))<=5*scale]
    words = np.array([word for word,a,c in allwords if a==5*scale or c<=3],dtype='int64')
    lengths = [np.array([int(scale*w) for w,_ in tab],dtype='int64') for tab in tables]
    energies = np.column_stack([w[words[:,j]] for j,w in enumerate(lengths)])
    heights = 5*scale-energies.sum(axis=1)
    groups = [[i for i,t in enumerate(components) if t==T] for T in dict.fromkeys(components)]
    canonical = np.column_stack([np.sort(energies[:,g],axis=1) for g in groups])
    keys,rep,tag = np.unique(canonical,axis=0,return_index=True,return_inverse=True)

    def weight_sum(u,m,n):
        result = np.zeros(len(words),dtype='int64')
        for j,(typ,w) in enumerate(zip(components,lengths)):
            if typ[0]=='D' and int(typ[1:])%2==0:
                idx = ((m%2)*u[j]) ^ ((n%2)*words[:,j])
            else:
                idx = (m*u[j]+n*words[:,j]) % len(w)
            result += w[idx]
        return result

    def allowed(v,par):
        return ((v==0)|(v>=5*scale)) & (v%(4*scale)==(par%2)*scale)
    cyclic = allowed(weight_sum(words[0],0,2)+4*heights,0) & allowed(weight_sum(words[0],0,3)+9*heights,1)
    active = cyclic.copy()
    # Track all positive-height keys which pass the one-class tests before the first cut.
    live = set(map(int,tag[cyclic&(heights>0)]))
    assert all(np.all((tag!=k)|cyclic) for k in live)
    all_checks = 0
    for key,b in cert['cuts']:
        assert key in live and isinstance(b,list) and all(isinstance(x,int) for x in b)
        j = rep[key]; u = words[j]; h = int(heights[j])
        blocks = defaultdict(list)
        for i,t in enumerate(components):
            blocks[t,int(energies[j,i])].append(i)
        dim = [len(g)*int(t[1:]) for (t,e),g in blocks.items()] + [1,17-sum(int(t[1:]) for t in components)]
        assert sum(dim)==18 and len(b)==len(dim)
        xwidth = 4*h*scale
        score = np.zeros(len(words),dtype='int64')
        for s,g in zip(b,blocks.values()):
            score += s*4*h*energies[:,g].sum(axis=1)
        assert 5*xwidth*max(map(abs,b)) < 2**56
        # congruence for the sum fixes the residue of y = 2 S <z_a,z_b>.
        plus = weight_sum(u,1,1)+h+heights
        residue = np.remainder(-plus,4*scale)
        probes = [(m,n,weight_sum(u,m,n)+m*m*h+n*n*heights)
                  for m,n in ((1,1),(-1,1),(2,1),(-2,1),(1,2),(-1,2),(3,1),(-3,1))]
        mn = None
        for k in range(-3,3):  # 0<=residue<4S and |y|<=10S.
            y = residue+4*scale*k
            d = 4*h*heights - y*y
            ok = (d>=0)&active
            for m,n,v in probes:
                ok &= allowed(v+m*n*y,m+n)
            if not np.any(ok):
                continue
            val = int((score[ok]+b[-2]*y[ok]**2+b[-1]*d[ok]).min())
            mn = val if mn is None else min(mn,val)
        assert mn is not None
        tr=sum(t*z for t,z in zip(dim,b))
        var=18*sum(t*z*z for t,z in zip(dim,b))-tr*tr
        gap=18*58*mn-290*tr*xwidth
        assert gap>0 and gap*gap>1508*var*xwidth*xwidth
        active &= (tag!=key)
        live.remove(key)
        all_checks += 1
    assert set(cert['live'])==live
    return bool(live),all_checks


def fsum(dim):
    n=dim
    def bound(r): return 24*r//(25-r)
    return bound(n+1)+bound(18-n)+4+8//(n-1)


def main():
    root=Path(sys.argv[1])/'root_profiles'
    archived=json.loads((root/'frontier9.json').read_text())['types']
    counts = [0]*4
    for rank in range(6,16):
        records=json.loads((root/f'lower_{rank}.json').read_text())
        names={'.'.join(t):t for t in root_types(rank)}
        assert set(names)==set(records)
        remaining=set();checks=0
        for name,typ in names.items():
            live,nb=first_anchor_check(typ,records[name])
            checks+=nb
            if live: remaining.add(name)
        used=set(archived[str(rank)])
        missing=remaining-used
        pure='.'.join(['A1']*rank)
        exceptional={pure}
        if rank in range(6,11):exceptional.add(f'A{rank}')
        if rank in (6,7):exceptional.add(f'D{rank}')
        if rank==6:exceptional.update(['A1.A1.A1.A1.A2','A1.A1.A2.A2'])
        assert missing==exceptional and used==remaining-exceptional
        print(rank,len(names),len(remaining),len(used),checks,flush=True)
        counts=[a+b for a,b in zip(counts,[len(names),len(remaining),len(used),checks])]
    print('TOTAL input types, live, sent downstream, certified cuts:',counts)
    assert counts==[1794,1097,1078,8001]
    for expr in ([('A1',4),('A2',1)],[('A1',2),('A2',2)]):
        radius=sum(int(t[1:])*(int(t[1:])+1)*n//2 for t,n in expr)
        score=sum(int(t[1:])*n*max(0,60-8*(radius-2*(int(t[1:])+1)+3)) for t,n in expr)
        assert score>84;print(expr,score)
    for d in (6,7,8,9,10):
        assert fsum(d)<58
        print('A',d,'upper',fsum(d))


if __name__=='__main__':
    main()
