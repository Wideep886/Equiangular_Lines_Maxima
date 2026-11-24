"""Independently check signed candidate counts by multiplying local theta polynomials.
Usage: python3 exceptional_shells.py PATH/TO/EL18_58_search/joint23
Only the recorded polarization coordinates and claimed counts are read from the
certificate. Local coordinate boxes, congruence conditions, glue and convolutions
use no code from the certificate. All energies and coordinates are integers.
"""
from collections import Counter, defaultdict
from functools import cache
from itertools import product
from math import isqrt
from pathlib import Path
from ast import literal_eval
import json, re, sys


def codes(n):
    if n==6:
        seq=[]
        for b in product(range(2),repeat=4):
            if sum(b)%2:continue
            for i,j in product(range(3),repeat=2):
                v=tuple((3*x+2*y)%6 for x,y in zip(b,(i,i+j,i-j,j)))
                seq.append(v+(b[1] ^ (2*b[2]) ^ (3*b[3]),))
        return seq
    if n==5:
        return [tuple((a+b*x+c*x*x)%5 for x in range(5))+(c,)
                for a,b,c in product(range(5),repeat=3)]
    if n==4:
        g=((2,1,1,1),(-1,2,1,-1),(-1,-1,2,1),(-1,1,-1,2))
        return [v+tuple(sum(a*b for a,b in zip(col,v))%4 for col in zip(*g))
                for v in product(range(4),repeat=4)]
    chi=(0,1,-1,-1,1)
    c=[[0 if i==j else 1 if not i or not j else chi[(i-j)%5] for j in range(6)] for i in range(6)]
    full=[v+tuple(sum(a*b for a,b in zip(col,v))%3 for col in zip(*c))
          for v in product(range(3),repeat=6)]
    return [v[:-1] for v in full if v[-1]==0]


@cache
def blocks(n,dim,tag,mark):
    """Enumerate the scaled vectors in one coset through norm six."""
    D=(n==6 and dim==4)
    if D:
        residue=3 if tag>=2 else 0
    else:
        residue=(tag if n in (3,4) else -tag)%n
    lim=6*n*n; maximum=isqrt(lim); vals=tuple(v for v in range(-maximum,maximum+1) if (v-residue)%n==0)
    out=[]
    for part in product(vals,repeat=(dim if D else dim-1)):
        v=part if D else part+(-sum(part),)
        if not D and (v[-1]-residue)%n:continue
        if mark and not(v[0]==v[1]==v[2]): continue
        if D:
            h=sum(v)//6
            if h%2 != tag%2:continue
        E=sum(x*x for x in v)
        if E<=lim:out.append((E,v))
    return tuple(out)


@cache
def local(n,dim,tag,mark,w):
    if w is None:
        return tuple(Counter((E,0) for E,v in blocks(n,dim,tag,mark)).items())
    bound=4*n*n
    return tuple(Counter((E,sum(a*b for a,b in zip(v,w))) for E,v in blocks(n,dim,tag,mark) if E<=bound).items())


def theta(n,words,mark,ws):
    factors=[4 if n==6 and i==4 else n for i in range(len(words[0]))]
    result=Counter()
    lim=(6 if ws is None else 4)*n*n
    tar=3*n*n
    for g in words:
        factors_hist=[local(n,d,t,i==mark,None if ws is None else ws[i]) for i,(t,d) in enumerate(zip(g,factors))]
        if any(not v for v in factors_hist):continue
        minima=[min(k[0] for k,m in piece) for piece in factors_hist]
        if sum(minima)>lim:continue
        previous={(0,0):1}
        left=sum(minima)
        for hist,lowest in zip(factors_hist,minima):
            left-=lowest
            current=defaultdict(int)
            for (n0,t0),a in previous.items():
                for (n1,t1),b in hist:
                    if n0+n1+left<=lim:
                        current[n0+n1,t0+t1]+=a*b
            previous=current
        for e in ([2,4,6] if ws is None else [4]):
            result[e]+=previous.get((e*n*n,0 if ws is None else tar),0)
    return result


def cases(root,n):
    if n==6:
        for line in (root/'a5/all.out').read_text().splitlines():
            m=re.fullmatch(r'([04]) (\d+) (\d+) (\d+) (\(.*\))',line)
            if m:
                j,k,naxes,_,s=m.groups()
                yield int(j),k,literal_eval(s),2*int(naxes)
    elif n==5:
        for k,row in enumerate(json.loads((root/'a4/orbits.json').read_text())):
            yield 0,k,tuple(map(tuple,row['w'])),row['count']
    elif n==4:
        for line in (root/'eleven/a3/octacode.out').read_text().splitlines():
            m=re.fullmatch(r'(\d+) h (\(.*\)) multiplicity \d+ cells \d+ points (\d+) bound \d+',line)
            if m:yield 0,m[1],literal_eval(m[2]),int(m[3])
    else:
        words=codes(3)
        zero=(0,0,0)
        def singleton(i,v):
            return tuple(v if j==i else zero for j in range(11))
        yield None,0,singleton(0,(6,-3,-3)),120
        yield None,1,((3,-3,0),)*3+(zero,)*8,444
        def build(word,tag):
            w=[zero if u==0 else (u,u,-2*u) if u==1 else (-1,-1,2) for u in word]
            if tag==3:
                i=next(i for i,t in enumerate(word) if t)
                w[i]=tuple(-2*x for x in w[i])
            if tag==2:
                i=word.index(0)
                w[i]=(3,-3,0)
            return tuple(w)
        w6=next(w for w in words if sum(bool(a) for a in w)==6)
        w9=next(w for w in words if sum(bool(a) for a in w)==9)
        for k,g,number in ((2,w6,420),(3,w6,432),(4,w9,408)):
            yield None,k,build(g,k),number


def main(root):
    for n,expected in [(6,166),(5,28),(4,17),(3,5)]:
        g=codes(n);rows=list(cases(root,n));assert len(rows)==expected
        for j in sorted({v[0] for v in rows},key=str):
            raw=theta(n,g,j,None)
            print('scale',n,'marked factor',j,'full shell sizes [2,4,6]',[raw[t] for t in (2,4,6)])
        c=Counter()
        for j,k,w,num in rows:
            assert sum(sum(a*a for a in v) for v in w)==6*n*n
            assert j is None or len(set(w[j][:3]))==1
            count=theta(n,g,j,w)[4]
            assert count==num and count%2==0,(n,j,k,count,num)
            c[num]+=1
        print('marks',len(rows),'oriented minima/maxima',min(c),max(c),'sum over representatives',sum(k*v for k,v in c.items()),flush=True)
if __name__=='__main__':
    main(Path(sys.argv[1]))
