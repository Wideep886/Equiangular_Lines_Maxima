"""Elementary exact Cartan and local coordinate utilities, from the v4 verifier."""
from fractions import Fraction
from itertools import product, permutations, combinations
from functools import lru_cache
from collections import Counter, defaultdict
from math import prod, lcm, isqrt, factorial

TYPES=[[('D',24)],[('D',16),('E',8)],[('E',8)]*3,[('A',24)],[('D',12)]*2,[('A',17),('E',7)],[('D',10),('E',7),('E',7)],[('A',15),('D',9)],[('D',8)]*3,[('A',12)]*2,[('A',11),('D',7),('E',6)],[('E',6)]*4,[('A',9)]*2+[('D',6)],[('D',6)]*4,[('A',8)]*3,[('A',7)]*2+[('D',5)]*2,[('A',6)]*4,[('D',4)]*6,[('A',5)]*4+[('D',4)],[('A',4)]*6]

def dot(a,b):return sum(x*y for x,y in zip(a,b))
def matvec(A,v):return tuple(dot(a,v) for a in A)
def cartan(ty):
    t,n=ty;C=[[2*(i==j) for j in range(n)] for i in range(n)]
    edges=[(i,i+1) for i in range(n-1 if t=='A' else n-2)]
    if t=='D':edges.append((n-3,n-1))
    elif t=='E':edges.append((2,n-1))
    for i,j in edges:C[i][j]=C[j][i]=-1
    return tuple(map(tuple,C))

def inverse_and_det(C):
    n=len(C);A=[[Fraction(x) for x in row]+[Fraction(i==j) for j in range(n)] for i,row in enumerate(C)];d=Fraction(1)
    for j in range(n):
        p=next(i for i in range(j,n) if A[i][j])
        if p!=j:A[p],A[j]=A[j],A[p];d=-d
        x=A[j][j];d*=x;A[j]=[a/x for a in A[j]]
        for i in range(n):
            if i!=j and A[i][j]:
                x=A[i][j];A[i]=[a-x*b for a,b in zip(A[i],A[j])]
    assert d.denominator==1 and d>0
    d=int(d);K=tuple(tuple(int(x*d) for x in row[n:]) for row in A)
    assert all(all(x*d==int(x*d) for x in row[n:]) for row in A)
    assert all(dot(C[i],tuple(K[k][j] for k in range(n)))==d*(i==j) for i in range(n) for j in range(n))
    return d,K

def weyl_order(C):
    n=len(C);todo=set(range(n));out=1
    while todo:
        p=min(todo);todo.remove(p);q=[p]
        for i in q:
            s={j for j in todo if C[i][j]==-1};todo-=s;q.extend(s)
        m=len(q);deg={i:sum(C[i][j]==-1 for j in q) for i in q};b=[i for i in q if deg[i]>2]
        if not b:out*=factorial(m+1);continue
        assert len(b)==1 and deg[b[0]]==3
        b=b[0];arms=[]
        for i in q:
            if C[b][i]!=-1:continue
            prev=b;k=1
            while deg[i]>1:
                j=next(j for j in q if j!=prev and C[i][j]==-1);prev,i=i,j;k+=1
            arms.append(k)
        arms.sort()
        if arms[:2]==[1,1]:out*=2**(m-1)*factorial(m)
        else:
            assert arms==[1,2,m-4] and m in (6,7,8);out*={6:51840,7:2903040,8:696729600}[m]
    return out

def budgets(tabs,cap1,cap2=None):
    """Assemble exact budgets from nonnegative local weights; no CS pruning."""
    if any(not T for T in tabs):return
    if cap2 is None:
        # tables: (norm, index) tuples
        tail=[0]*(len(tabs)+1)
        for j in range(len(tabs)-1,-1,-1):tail[j]=tail[j+1]+min(n for n,i in tabs[j])
        def visit(j,a,P):
            if j==len(tabs):
                if a==cap1:yield tuple(P)
                return
            for n,i in tabs[j]:
                if a+n+tail[j+1]<=cap1:yield from visit(j+1,a+n,P+[i])
        yield from visit(0,0,[])
    else:
        tail1=[0]*(len(tabs)+1);tail2=tail1.copy()
        for j in range(len(tabs)-1,-1,-1):
            tail1[j]=tail1[j+1]+min(a for a,b in tabs[j]);tail2[j]=tail2[j+1]+min(b for a,b in tabs[j])
        def visit(j,a,b,P):
            if j==len(tabs):
                if a==cap1 and b==cap2:
                    for t in product(*P):yield sum(t,())
                return
            for (s,t),A in tabs[j].items():
                if a+s+tail1[j+1]<=cap1 and b+t+tail2[j+1]<=cap2:yield from visit(j+1,a+s,b+t,P+[A])
        yield from visit(0,0,0,[])

H4=((3,3,2,3,1,0,0,0),(3,0,3,2,3,1,0,0),(3,0,0,3,2,3,1,0),(3,0,0,0,3,2,3,1))
P3=((0,1,1,1,1,1),(1,0,1,1,2,2),(1,1,0,2,1,2),(1,1,2,0,2,1),(1,2,1,2,0,1),(1,2,2,1,1,0))
G4=(((2,3,1,7,0,6,5,4),(1,-1,-1,1,-1,1,1,-1)),((7,2,3,4,1,6,0,5),(-1,1,1,-1,1,-1,-1,1)))
G3=((tuple(range(12)),(-1,)*12),((0,2,4,1,5,3,6,8,10,7,11,9),(1,)*12),((0,2,1,4,3,5,6,8,7,10,9,11),(1,)*12),((0,1,2,3,11,10,7,6,9,8,5,4),(1,1,1,1,-1,-1,1,1,1,1,-1,-1)))
@lru_cache(None)
def small_local(m,cap,dominant=False):
    R=isqrt(cap);out=defaultdict(list)
    for c in range(m):
        opts=[a for a in range(-R,R+1) if (a+c)%m==0]
        for a in product(opts,repeat=m-1):
            b=-sum(a)
            if b not in opts:continue
            v=a+(b,);n=dot(v,v)
            if n<=cap and (not dominant or v==tuple(sorted(v,reverse=True))):out[c].append((n,v))
    return out

