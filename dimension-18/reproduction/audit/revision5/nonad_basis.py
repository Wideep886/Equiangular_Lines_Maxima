from itertools import product,permutations
from collections import Counter
order=(4,9,0,1,2,3,5,6,7,8,10)
sgn=(1,1,-1,-1,-1,-1,-1,-1,1,1,-1)
X=(1,1,0,2,0,0,1,2,2)
Y=(1,0,1,0,2,0,2,1,2)
F=list(zip(X,Y))
assert set(F)==set(product(range(3),repeat=2))
q0=lambda x,y:(x*x+y*y)%3
q1=lambda x,y:(2*x*x+x*y+y*y)%3
urows=[(1,-1,-1,-1,1,-1),(-1,1,-1,1,0,-1),(-1,-1,-1,-1,0,-1),(-1,-1,0,1,0,0),(-1,0,-1,0,0,1)]
def paley(u):
    vals=[sum(u[1:])]+[u[0]+sum(u[i]*(0 if i==j else 1 if (i-j)%5 in (1,4) else -1) for i in range(1,6)) for j in range(1,6)]
    c=tuple(u)+tuple(vals)
    return c[-1]%3, tuple(s*c[j]%3 for s,j in zip(sgn,order))
assert all(paley(u)==(0,tuple([int(r==0),int(r==1)]+[[q0(x,y),q1(x,y),1,x,y][r] for x,y in F])) for r,u in enumerate(urows))
actual={val for u in product(range(3),repeat=6) for check,val in [paley(u)] if check==0}
expected={(a,b)+tuple((a*q0(x,y)+b*q1(x,y)+c+d*x+e*y)%3 for x,y in F) for a,b,c,d,e in product(range(3),repeat=5)}
assert len(actual)==len(expected)==243 and actual==expected
# Ordered cube-height sign conventions. These are center entries equal to ±1 or zero (a midpoint).
mini={s:sorted(set(permutations((s,s,-2*s)))) for s in (1,-1)}
h=(-2,1,1)
axis={};
for c in expected:
    if c.count(0)!=5:continue
    iword=c[2:];z=sum(t!=0 for t in c[:2]);p=iword.count(1);n=iword.count(2)
    # Native units: host vectors times 3. Target dot 27.
    opts={s:tuple(sorted({sum(x*y for x,y in zip(r,h)) for r in val})) for s,val in [(0,[(0,0,0)]),(1,mini[1]),(2,mini[-1])]}
    possible=any(sum(v)==27 for v in product(*(opts[s] for s in iword)))
    assert possible==(z==0 or (z==1 and p==4 and n==1))
    if not possible:continue
    assert sum(max(opts[s]) for s in iword)==27 and (p,n) in [(3,3),(4,1)]
    if z==1 and 1 in c[:2]:
        j=c[:2].index(1);i=iword.index(2);xx,yy=F[i];q=(q0,q1)[j]
        assert iword==tuple((q(x-xx,y-yy)-1)%3 for x,y in F)
        axis[j,i]=tuple(1 if ss==1 else -1 if ss==0 else 0 for ss in iword)
assert len(axis)==18
dirs=((1,0),(0,1),(1,1),(1,-1))
for j,q in enumerate([q0,q1]):
    eps=[1 if q(dx,dy)==2 else -1 for dx,dy in dirs]
    for i,hc in enumerate(F):
        center=[sum(eps[k] for k,(dx,dy) in enumerate(dirs) if (dx*(y-hc[1])-dy*(x-hc[0]))%3==0) for x,y in F]
        assert tuple(center)==axis[j,i]
    print("rectangle block signs",j,eps)
print("PASS all five Paley rows, all 243 shortened words and Cartesian coordinates, z=0,1,2, unique centers and pencil signs")
