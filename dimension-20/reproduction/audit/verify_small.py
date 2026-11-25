#!/usr/bin/env python3
"""Small glue models: scalar patterns, local product bounds and projections.
No global candidate vectors, no compatibility graphs or search certificates.
"""
from cartan_utils import *
from energy_bound import cap,combine
from fractions import Fraction as F
from math import gcd
import sys
if not __debug__:raise RuntimeError('Run with assertions enabled')
@lru_cache(None)
def smoptions(m,w):
    from math import factorial
    Z=defaultdict(list)
    for i,a in enumerate(w):Z[a].append(i)
    O=prod(factorial(len(I)) for I in Z.values());out=defaultdict(list)
    for a,S in small_local(m,4*m*m).items():
        for n,u in S:
            if any([u[i] for i in I]!=sorted((u[i] for i in I),reverse=True) for I in Z.values()):continue
            v=tuple(b-a for a,b in zip(u,w));n2=dot(v,v)
            if n2>4*m*m:continue
            y=tuple(2*b-c for b,c in zip(u,w));es=[]
            for I in Z.values():
                p=[y[i] for i in I];e=(F(dot(p,p))-F(sum(p)**2,len(I)))/(2*m*m)
                if e:es.append((len(I)-1,e,factorial(len(I))//prod(factorial(k) for k in Counter(u[i] for i in I).values())))
            size=O//prod(factorial(k) for I in Z.values() for k in Counter(u[i] for i in I).values())
            v=list(v)
            for I in Z.values():
                for i,b in zip(I,sorted((v[i] for i in I),reverse=True)):v[i]=b
            out[a].append((n,n2,u,tuple(v),tuple(es),size))
    return out

@lru_cache(None)
def local_ys(w,u):
    blocks=defaultdict(list)
    for i,x in enumerate(w):blocks[x].append(i)
    opts=[sorted(set(permutations(u[i] for i in I))) for I in blocks.values()]
    Y=set()
    for A in product(*opts):
        v=[0]*len(w)
        for I,a in zip(blocks.values(),A):
            for i,b in zip(I,a):v[i]=b
        Y.add(tuple(2*a-b for a,b in zip(v,w)))
    return tuple(sorted(Y))

def ray(y):
    a=gcd(*y);assert a>0
    y=tuple(x//a for x in y)
    if next(x for x in y if x)<0:y=tuple(-x for x in y)
    return y

class Cell:
    def __init__(self,w,P,selfpart):
        self.selfpart=selfpart
        self.es=sum((t[4] for t in P),())
        self.h=5-sum(e for d,e,t in self.es)
        self.size=prod(t[5] for t in P)
        self.b=cap(self.es,self.h,self.size,selfpart,fibers=True)
        self.local=tuple(local_ys(z,t[2]) for z,t in zip(w,P))
        assert all(len(S)==t[5] for S,t in zip(self.local,P))
        self.es_by_factor=tuple(t[4] for t in P)
        self.proj=tuple(bool(any(S[0])) for S in self.local)
        assert all(len({dot(y,y) for y in S})==1 for S in self.local)
        # Every cell projection is either identically zero or always nonzero.
        self.linear=tuple(dot(S[0],z) for S,z in zip(self.local,w))
        assert all(len({dot(y,z) for y in S})==1 for S,z in zip(self.local,w))
        self.rays=[]
        for j,S in enumerate(self.local):
            B=Counter(ray(y) for y in S if any(y))
            if not B:self.rays.append({});continue
            # Fix a nonzero projected vector in this factor; the other factors
            # still range over their complete product orbits.
            es=sum((t[4] for i,t in enumerate(P) if i!=j),())
            h=self.h+sum(e for d,e,t in P[j][4])
            size=self.size//len(S)
            b=cap(es,h,size,False,fibers=True)
            if selfpart:
                assert set(B.values())=={2}
                # Orient each axis to have its projection in the positive ray.
                self.rays.append({r:min(self.b,b) for r in B})
            else:
                self.rays.append({r:min(self.b,b*k) for r,k in B.items()})

def cells(m,C,w):
    tabs=[smoptions(m,z) for z in w];R={}
    for c in C:
        for P in combine([t.get(a,[]) for t,a in zip(tabs,c)],4*m*m):
            u=tuple(t[2] for t in P);v=tuple(t[3] for t in P)
            assert u not in R;R[u]=(v,P)
    ans=[]
    for u,(v,P) in R.items():
        assert v in R and R[v][0]==u
        if u>v:continue
        ans.append(Cell(w,P,u==v))
    return ans

def raycap(C,j):
    D=Counter()
    for c in C:D.update(c.rays[j])
    return max(D.values(),default=0)

def solve(C,n):
    start=sum(c.b for c in C);deletions=[]
    while sum(c.b for c in C)>90:
        # Invariant functional y -> (y_j,w_j).
        for j in range(n):
            T=[c for c in C if c.linear[j]];b=sum(c.b for c in T)
            if T and b<=14:
                deletions.append(('linear',j,b));C=[c for c in C if not c.linear[j]];break
        else:
            for j in range(n):
                T=[c for c in C if c.proj[j]];b=sum(c.b for c in T)
                if not T:continue
                R=raycap(T,j)
                if b<=14 or (R<=14 and b<=30):
                    deletions.append(('factor',j,b,R));C=[c for c in C if not c.proj[j]];break
            else:break
    b=sum(c.b for c in C)
    if b<=90:return start,b,deletions,None
    J=[j for j in range(n) if raycap(C,j)<=14 and sum(c.b for c in C if not c.proj[j])<=90]
    t=max((sum(c.proj[j] for j in J) for c in C),default=0)
    if 31*len(J)>91*t:return start,90,deletions,(len(J),t)
    return start,b,deletions,('UNSOLVED',len(J),t)

def scalar_patterns(m):
    assert m in (3,4);n=24//(m-1)
    if m==3:C=sorted({a+tuple(sum(a[i]*P3[i][j] for i in range(6))%3 for j in range(6)) for a in product(range(3),repeat=6)})
    else:C=[a for a in product(range(4),repeat=8) if all(dot(h,a)%4==0 for h in H4)]
    CS=set(C);assert len(C)==m**(n//2)
    for c in C:
        N=sum(a*(m-a) for a in c);assert N%(2*m)==0 and (N>=4*m or not any(c))
    G=G3 if m==3 else G4
    for p,s in G:assert sorted(p)==list(range(n)) and set(s)<={-1,1} and all(tuple(s[i]*c[p[i]]%m for i in range(n)) in CS for c in C)
    D=small_local(m,6*m*m,True);H=defaultdict(list)
    for a,S in D.items():
        for N,v in S:H[a,N].append(v)
    W=set();total=0
    for c in C:
        for N in budgets([sorted({(N,N) for N,v in D[a]}) for a in c],6*m*m):
            assert (c,N) not in W;W.add((c,N))
            total+=prod(sum(factorial(m)//prod(factorial(k) for k in Counter(v).values()) for v in H[a,b]) for a,b in zip(c,N))
    assert total==16773120+6048*m
    todo=W.copy();ans=[]
    while todo:
        z=min(todo);todo.remove(z);q=[z]
        for c,N in q:
            for p,s in G:
                z2=(tuple(s[i]*c[p[i]]%m for i in range(n)),tuple(N[p[i]] for i in range(n)))
                assert z2 in W
                if z2 in todo:todo.remove(z2);q.append(z2)
        ans.append((z,len(q)))
    assert sum(k for z,k in ans)==len(W)
    return C,H,len(W),ans

def check_small(m):
    C,H,num,patterns=scalar_patterns(m);upper=0;tests=0
    groups=defaultdict(list)
    for i,((c,N),orbit_size) in enumerate(patterns):
        local=[H[a,b] for a,b in zip(c,N)];assert sum(len(A)>1 for A in local)<=1
        # Only local choices of the one possible ambiguous factor are tested.
        for w in product(*local):
            cells0=cells(m,C,w);row=solve(cells0,len(w));tests+=1
            assert row[1]<=90 and not (row[3] and row[3][0]=='UNSOLVED'),(m,i,c,N,row)
            upper=max(upper,row[1]);groups[tuple(sorted(x for x in N if x))].append(row[1])
            print('SMALL',m,'scalar_orbit',i,'signed_pair_cells',len(cells0),'marking_count',orbit_size,'result',row,flush=True)
    for n,B in sorted(groups.items()):print('SHAPE',m,n,'bound',max(B),flush=True)
    print('SMALL VERIFIED',m,'scalar patterns',num,'scalar orbits',len(patterns),'local alternatives',tests,'bound',upper,flush=True)
    return upper

if __name__=='__main__':
    for m in (map(int,sys.argv[1:]) if len(sys.argv)>1 else (3,4)):check_small(m)
