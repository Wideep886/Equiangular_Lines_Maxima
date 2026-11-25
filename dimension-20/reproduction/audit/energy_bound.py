"""The common orbit-energy bound and optional local fiber bounds."""
from fractions import Fraction as F
from math import prod

def cap(es,h,size,selfpart,*,fibers=False):
    # es contains (root-space dimension, squared projected norm, orbit size)
    assert h>=0 and h+sum(e for d,e,t in es)==5
    assert prod(t for d,e,t in es)==size
    assert not selfpart or (h==0 and size%2==0)
    r=sum(d for d,e,t in es)+(h>0);assert 1<=r<=23
    f=2 if selfpart else 1
    E=h*h+sum(e*e/d for d,e,t in es);assert E>1
    b=min(size//f,24//(E-1))
    if h>1:b=min(b,4//(h-1))
    if fibers:
        for d,e,t in es:
            if e<4 and (h>0 or len(es)>1):b=min(b,4//(4-e)*(size//t)//f)
        sig=sorted(es)
        if selfpart and sig==[(1,F(1),2)]*5:b=min(b,5)
        if selfpart and sig==[(1,F(1),2),(3,F(2),6),(3,F(2),6)]:b=min(b,8)
    return b

def combine(tabs,cap):
    """Exact two norm budgets for local signature records; no vector products."""
    if any(not t for t in tabs):return
    n=len(tabs);tail1=[0]*(n+1);tail2=tail1.copy()
    for j in range(n-1,-1,-1):tail1[j]=tail1[j+1]+min(t[0] for t in tabs[j]);tail2[j]=tail2[j+1]+min(t[1] for t in tabs[j])
    def visit(i,a,b,P):
        if i==n:
            if a==b==cap:yield tuple(P)
            return
        for t in tabs[i]:
            c,d=t[:2]
            if a+c+tail1[i+1]<=cap and b+d+tail2[i+1]<=cap:yield from visit(i+1,a+c,b+d,P+[t])
    yield from visit(0,0,0,[])
