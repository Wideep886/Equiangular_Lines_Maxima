"""Independent coordinate quotient and fibre counts for the 31 affine-plane partitions."""
from collections import Counter
from itertools import product, combinations
from math import gcd

P=list(product(range(3),repeat=2))
# first four pencils consist of fixed x and fixed y-mx, m=0,1,2
pencils=[]
for coef in [(1,0),(0,1),(2,1),(1,1)]:
    pencils.append([tuple(i for i,(x,y) in enumerate(P) if (coef[0]*x+coef[1]*y)%3==c) for c in range(3)])
lines=sum(pencils,[])
def partitions(pre,unused):
    if not unused:yield tuple(pre); return
    i=unused[0];rest=unused[1:]
    for size in range(len(rest)+1):
        for tail in combinations(rest,size):
            block=(i,)+tail
            yield from partitions(pre+[block],[v for v in rest if v not in tail])

def image_partition(p,g):
    return tuple(sorted(tuple(sorted(g[i] for i in B)) for B in p))

def trim(w):
    if not any(w):return w
    n=gcd(*w)
    if next(v for v in w if v)<0:n=-n
    return tuple(x//n for x in w)

def quotient(w,v):
    pivot=next(i for i,x in enumerate(v) if x)
    return tuple(v[pivot]*x-w[pivot]*y for i,(x,y) in enumerate(zip(w,v)) if i!=pivot)

def plus(x,y):return tuple(a+b for a,b in zip(x,y))
def minus(x,y):return tuple(a-b for a,b in zip(x,y))
def neg(x):return tuple(-a for a in x)
def allsum(columns):
    start={(0,)*len(columns[0]):1}
    for x in columns:
        new=Counter()
        for s,c in start.items():
            new[plus(s,x)]+=c;new[minus(s,x)]+=c
        start=new
    return start
def inspect(p,cnt):
    k=len(p); labels=[next(j for j,B in enumerate(p) if i in B) for i in range(9)]
    unit=[tuple(int(i==j) for i in range(k)) for j in range(k)]
    F=set(allsum([unit[i] for i in labels]))
    local=[list(allsum([unit[labels[i]] for i in line])) for line in lines]
    rel=set(F)
    for C in local:rel.update(minus(x,y) for x,y in combinations(C,2))
    for B in combinations(unit,3):rel.update(allsum(B))
    best=0;high=set()
    for v in {trim(f) for f in F if min(f)>=0}:
        cols1=[quotient(x,v) for x in unit]
        if any(not any(x) for x in cols1):continue
        # Each relevant vector can supply the second independent equation.
        extra={None}
        if k==5:extra.update(trim(q) for a in rel for q in [quotient(a,v)] if any(q))
        for w in extra:
            cols=[quotient(x,w) for x in cols1] if w is not None else cols1
            # Classes mean equality up to sign, not proportionality.
            if any(not any(x) for x in cols) or len({min(x,neg(x)) for x in cols})<k:continue
            cube=[allsum([cols[labels[i]] for i in line]) for line in lines]
            rect=max(sum(sorted(cnt,reverse=True)[:2]),max([0]+[sum(cnt[i] for i in inds) for inds in combinations(range(k),3) if (0,)*len(cols[0]) in allsum([cols[i] for i in inds])]))
            value=18+2*rect
            for j in range(4):
                # Occupied centres from one parallel class sum to zero.
                X,Y,Z=cube[3*j:3*j+3]
                cap=sum(sorted(max(C.values()) for C in [X,Y,Z])[-2:])
                for x,m in X.items():
                    for y,n in Y.items():
                        s=neg(plus(x,y))
                        if s in Z:cap=max(cap,m+n+Z[s])
                value+=cap
            best=max(best,value)
            if value>=58:
                assert w is None and k==4 and v==tuple(0 if len(B)==6 else 1 for B in p) and rect==7
                high.add((value,tuple(sorted((len(B) for B in p),reverse=True)),
                          any(all(len(p[labels[i]])==1 for i in line) for line in lines)))
    return best,high
def main():
    affine=[[3*((a*x+b*y+i)%3)+(c*x+d*y+j)%3 for x,y in P]
            for a,b,c,d,i,j in product(range(3),repeat=6) if (a*d-b*c)%3]
    seen=set();caps=Counter();limits={};qual=set();highn=0
    for B in partitions([],list(range(9))):
        key=image_partition(B,list(range(9)))
        if key in seen:continue
        seen.update(image_partition(B,g) for g in affine)
        if len(B)<4:continue
        mult=[len(g) for g in B];hist=tuple(sorted(mult,reverse=True))
        mono=sum(any(set(L)<=set(g) for g in B) for L in lines)
        if 42+2*sum(hist[:3])+mono<58:continue
        assert len(B)<=5
        score,passed=inspect(B,mult);limits[hist]=max(score,limits.get(hist,0));caps[hist]+=1
        if passed:qual.update(passed);highn+=1
    assert len(seen)==21147
    assert sum(caps.values())==31 and highn==1 and qual=={(58,(6,1,1,1),True)}
    for h in sorted(limits,reverse=True):
        print(''.join(map(str,h)),caps[h],limits[h])
    print('Independent quotient: one surviving partition, six equal steps and three collinear singleton positions.')
if __name__=='__main__':main()
