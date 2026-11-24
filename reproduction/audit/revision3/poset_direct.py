"""Enumerate J-dominant local weights in an ellipsoid without expanding Weyl orbits.
Use Python 3.11+; first argument is the original joint23 directory. No saved
polarizations or chains are used. Requires poset_width.py in this directory.
"""
from collections import Counter, defaultdict
from functools import cache
from math import lcm, isqrt
import sys
from poset_width import factors, glue, shapes, degrees, quotient, HOSTS

@cache
def prefix(f):
    n=f.n;d=1;adj=[];out=[]
    for i in range(n):
        c=f.C[i][:i]
        v=[sum(a*b for a,b in zip(row,c)) for row in adj]
        e=2*d-sum(a*b for a,b in zip(c,v))
        cols=[(j,y) for j,y in enumerate(v) if y]
        out.append((d,e,cols))
        adj=[[((e*adj[j][k]+v[j]*v[k])//d) for k in range(i)]+[-v[j]] for j in range(i)]+[[-y for y in v]+[d]]
        d=e
    assert all(adj[i][j]*f.d==d*f.Z[i][j] for i in range(n) for j in range(n))
    return out

@cache
def local_weights(f,jset,mark):
    # A state retains the fixed Dynkin labels and the minimum square along
    # the fibre, exactly the squared norm of their projection on this prefix.
    n=f.n;out=[];p=prefix(f);u=[]
    def walk(i,h):
        if i==n:
            z=tuple(sum(t*s for t,s in zip(row,u)) for row in f.Z)
            c=f.cls.index(tuple(v%f.d for v in z))
            energy=sum(a*b for a,b in zip(z,u))
            assert energy*p[-1][1]==h*f.d
            k=tuple(j for j in jset if not u[j])
            out.append((c,energy,z,degrees(f,k)))
            return
        d,e,col=p[i]
        t=sum(u[j]*v for j,v in col)
        rem=(4*d-h)*e
        a=isqrt(rem)
        lo=max(-2,-((a-t)//d))
        hi=min(2,(t+a)//d)
        if i in jset:lo=max(0,lo)
        if i in mark:lo,hi=max(0,lo),min(0,hi)
        for x in range(lo,hi+1):
            numerator=h*e+(d*x-t)**2
            assert numerator%d==0
            u.append(x);walk(i+1,numerator//d);u.pop()
    walk(0,0)
    return tuple(out)

@cache
def groups(f,w,mark,L):
    js=tuple(i for i,c in enumerate(w) if c==0 and all(f.C[j][i]==0 for j in mark))
    result=[Counter() for _ in f.cls]
    for c,norm,z,den in local_weights(f,js,mark):
        dot=sum(k*b for k,b in zip(w,z))
        result[c][(norm*(L//f.d),dot*(L//f.d),den)]+=1
    return degrees(f,js),[tuple((e,ip,qs,m) for (e,ip,qs),m in v.items()) for v in result]

@cache
def add(a,b):return tuple(sorted(a+b))

def bound(fs,code,w,mark,L):
    gs=[groups(f,h,m,L) for f,h,m in zip(fs,w,mark)]
    num=tuple(sorted(n for d,_ in gs for n in d))
    total=0
    for c in code:
        options=[g[1][k] for g,k in zip(gs,c)]
        minima=[min((v[0] for v in group),default=5*L) for group in options]
        if sum(minima)>4*L:continue
        tails=[sum(minima[i+1:]) for i in range(len(options))]
        def walk(i,e,t,ds,m):
            if i==len(options):return m*quotient(num,ds) if (e,t)==(4*L,3*L) else 0
            return sum(walk(i+1,e+f,t+s,add(ds,d),m*n) for f,s,d,n in options[i] if e+f+tails[i]<=4*L)
        total+=walk(0,0,0,(),1)
    return total

def go(name):
    fs=factors(name.split('_'));code=glue(fs);hs=shapes(fs,code,6);L=lcm(*(f.d for f in fs));n=len(fs)
    null=((),)*n;weak=marks=miss=0;bds=[];ub=0; masses=[]
    for w in hs:
        b=bound(fs,code,w,null,L);ub=max(b,ub)
        if b>=116:
            weak+=1;caps=[]
            for k,(f,h) in enumerate(zip(fs,w)):
                for i in range(f.n):
                    for j in range(i):
                        if h[i] or h[j] or f.C[i][j]!=-1:continue
                        mt=list(null);mt[k]=(j,i)
                        caps.append(bound(fs,code,w,tuple(mt),L))
            marks+=len(caps);miss+=sum(t>=116 for t in caps)
            b=min(b,max(caps,default=0))
        bds.append(b)
    print(name,len(hs),'unmarked max',ub,'unmarked unresolved',weak,'marks tested',marks,'unresolved marks',miss,'effective max',max(bds),'unresolved vectors',sum(b>=116 for b in bds),flush=True)

if __name__=='__main__':
    for name in sys.argv[2:] or HOSTS:go(name)
