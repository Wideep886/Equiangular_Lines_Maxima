"""Reconstruct initial rank-16/17/18 root domains directly from Cartan matrices.
Input: path to the extracted EL18_58_search directory. Python 3.9+; stdlib only.
The generator never imports the original root-system or certificate code.
"""
import json
import sys
from collections import Counter
from fractions import Fraction as Rat
from functools import lru_cache
from itertools import combinations_with_replacement
from math import isqrt
from pathlib import Path

def unit(i,n):
    return tuple(int(j==i) for j in range(n))

def inverse(mat):
    n=len(mat)
    a=[list(map(Rat,row+list(unit(i,n)))) for i,row in enumerate(mat)]
    det=Rat(1)
    for j in range(n):
        p=a[j][j]
        assert p>0
        det*=p
        a[j]=[x/p for x in a[j]]
        for i in range(n):
            if i!=j and a[i][j]:
                t=a[i][j]
                a[i]=[x-t*y for x,y in zip(a[i],a[j])]
    return [row[n:] for row in a],int(det)

def diagram(typ):
    n=int(typ[1:])
    if typ[0]=='A': edges=[(i,i+1) for i in range(n-1)]
    else:
        arms=[1,1,n-3] if typ[0]=='D' else [1,2,n-4]
        edges=[]; idx=1
        for k in arms:
            parent=0
            for i in range(k):
                edges.append((parent,idx));parent=idx;idx+=1
        assert idx==n
    mat=[[2*int(i==j) for j in range(n)] for i in range(n)]
    for i,j in edges:mat[i][j]=mat[j][i]=-1
    return mat

@lru_cache(None)
def classes(typ):
    mat=diagram(typ);n=len(mat); inv,det=inverse(mat)
    roots=set(unit(i,n) for i in range(n));todo=list(roots)
    while todo:
        x=todo.pop()
        for j,row in enumerate(mat):
            y=list(x);y[j]-=sum(a*b for a,b in zip(row,x));y=tuple(y)
            if y not in roots: roots.add(y);todo.append(y)
    theta=max(roots,key=sum)
    # The minuscule classes are the trivial class and the mark-one nodes.
    zero=(Rat(0),)*n
    cls={zero:(Rat(0),0)}
    for j,mark in enumerate(theta):
        if mark!=1:continue
        dyn=list(unit(j,n))
        while any(x>0 for x in dyn):
            k=next(i for i,x in enumerate(dyn) if x>0)
            t=dyn[k];dyn=[x-t*m for x,m in zip(dyn,mat[k])]
        span=inv[j][j]-sum(v*y for v,y in zip(inv[j],dyn))
        assert span.denominator==1
        key=tuple(q%1 for q in inv[j])
        cls[key]=(2*inv[j][j],int(span))
    assert len(cls)==det
    return cls

@lru_cache(None)
def table(typ):
    cls=classes(typ)
    result={}
    for label,(energy,width) in cls.items():
        doubl=tuple((2*a)%1 for a in label);tripl=tuple((3*a)%1 for a in label)
        v=(width,cls[doubl][0],cls[tripl][0])
        if energy in result: assert result[energy]==v
        result[energy]=v
    return tuple(sorted((e,*v) for e,v in result.items()))

def dynkin_types(rank):
    atoms=sorted([f'{t}{r}' for r in range(1,rank+1) for t in 'ADE'
                  if t=='A' or (t=='D' and r>=4) or (t=='E' and r in (6,7,8))])
    def choose(left, start):
        if not left:
            yield ()
        for j in range(start,len(atoms)):
            t=atoms[j];r=int(t[1:])
            if r<=left:
                for tail in choose(left-r,j): yield (t,)+tail
    return choose(rank,0)

@lru_cache(None)
def block(t,m):
    out=[]
    for choose in combinations_with_replacement(table(t),m):
        e=sum(v[0] for v in choose)
        if e<=5:
            out.append((tuple(x[0] for x in choose),tuple(map(sum,zip(*choose)))))
    return out

def configurations(types):
    items=list(Counter(types).items())
    def visit(g,signature,values):
        if g==len(items):
            yield signature,values
            return
        for profile,u in block(*items[g]):
            n=tuple(x+y for x,y in zip(values,u))
            if n[0]<=5:
                yield from visit(g+1,signature+profile,n)
    return visit(0,(),(0,0,0,0))

@lru_cache(None)
def squarefree(h):
    v=h.numerator*h.denominator
    for k in range(isqrt(v),0,-1):
        if v%(k*k)==0:return v//(k*k)

def valid(z,p):
    return z%4==p and (z==0 or z>=5)

def main():
    p=Path(sys.argv[1])/'root_profiles'
    for rank,filename in [(16,'plane16.json'),(17,'codim1.json'),(18,'anchor_all.json')]:
        data=json.loads((p/filename).read_text())
        names=list(dynkin_types(rank))
        data={'.'.join(sorted(k.split('.'))):v for k,v in data.items()}
        assert len(names)==len(data) and set(data)=={'.'.join(t) for t in names}
        total=0; cuts=0; cells=0
        for typ in names:
            name='.'.join(typ); rec=data[name]
            configs=list(configurations(typ))
            if rank==18:
                cnt=sum(x[0]==5 for _,x in configs)
                assert rec['orbits']==cnt
                total+=cnt
            elif rank==17:
                active=[h for _,(u,w,d,t) in configs if (h:=5-u)>0 and valid(d+4*h,0)]
                qs={squarefree(Rat(z)) for z in active}
                assert qs==set(map(int,rec))
                total+=len(qs)
                cuts+=sum(row is not None for row in rec.values())
            else:
                # The ranks in plane16.json are indices in the order of the full
                # energy profiles, including zero-residual and pre-cyclic profiles.
                c=sorted((prof,x) for prof,x in configs if x[0]==5 or x[1]<=3)
                live={j for j,(_,(u,w,d,t)) in enumerate(c) if (h:=5-u)>0
                      and valid(d+4*h,0) and valid(t+9*h,1)}
                assert set(map(int,rec))==live,(name,live,set(map(int,rec)))
                total+=len(live);cells+=sum(len(v) for v in rec.values())
        print(rank,'types',len(data),'initial records',total,'rank17 trace cuts',cuts,'stored rank16 area sectors',cells)

if __name__=="__main__":
    main()
