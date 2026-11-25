#!/usr/bin/env python3
"""Necessary-condition relaxation for 73 equiangular lines and a saturated ADE root system.
All inequalities and all exclusion certificates checked with Fraction arithmetic.
SciPy is used only to FIND multipliers. `--verify` only uses the standard library.
"""
from fractions import Fraction as F
from itertools import combinations, combinations_with_replacement, product
from collections import Counter
from math import comb, isqrt, prod
from functools import lru_cache
import json, argparse

N=73
TYPES=[('A',i) for i in range(1,20)]+[('D',i) for i in range(4,11)]+[('E',6),('E',7)]
def B(d): return min(24*d//(25-d),{14:28,16:40,17:48,18:59}.get(d,10000))
def gen(r=0,i=0,A=()):
    if 11<=r<=19 and len(set(A))>1 and any(t=='A' for t,d in A): yield A
    for j in range(i,len(TYPES)):
        if r+TYPES[j][1]<=19: yield from gen(r+TYPES[j][1],j,A+(TYPES[j],))

def cox(t):
    a,d=t
    return d+1 if a=='A' else 2*d-2 if a=='D' else 12 if d==6 else 18

@lru_cache(None)
def slots(t):
    a,d=t
    # w, delta, slot group, orbit size (or a safe upper bound), fixed-other fiber bound
    if a=='A':return tuple((F(2*k*(d+1-k),d+1),k,'s' if k==1 else 'l',comb(d+1,k),comb(d+1,k-1)//k) for k in range(1,(d+1)//2+1))
    if a=='D':
        if d==4:return tuple((F(2),2,k,8,2) for k in ('v','l','l2'))
        return ((F(2),2,'v',2*d,2),(F(d,2),d//2,'l',2**(d-1),2**d//(d+1)))
    return ((F(8,3),2,'l',27,27),) if d==6 else ((F(3),3,'l',56,56),)

def info(p):return tuple(slots(t)[i] for t,i in p)

@lru_cache(None)
def cap(p):
    I=info(p); w=sum(x[0] for x in I); d=sum(t[1] for t,i in p); h=5-w
    if h<0 or (h>0 and sum(x[1] for x in I)>3):return 0
    a=min(N,B(d+int(h>0)))
    if h>1:a=min(a,int(4/(h-1)))
    for x in I:a=min(a,x[4]*prod(y[3] for y in I)//x[3])
    if all(t==('A',1) for t,i in p) and len(p)==5:a=min(a,5)
    return a

def max_pairs(m,L):
    a,b=divmod(m,L);return a*comb(L,2)+comb(b,2)

@lru_cache(None)
def deficit(p,m):
    I=info(p); w=[x[0] for x in I]; T=sum(w)*m
    h=5-sum(w)
    # A pair may contribute at most 2 in total by the tight-frame identity.
    # For a pair with projection difference a sum of two roots in different
    # components, weights are 2/h_j and 2/h_l. Within an A component,
    # a difference of square 8 has exactly four root lines of pairing 2.
    cross=[F(2,cox(t)) for t,i in p]
    within=[(F(8,cox(t)) if t[0]=='A' else F(2)) if x[1]>=2 else F(0) for (t,i),x in zip(p,I)]
    paircap=max([F(0)]+within+[cross[i]+cross[j] for i,j in combinations(range(len(p)),2)])
    bj=[max(within[i],cross[i] if len(p)>1 else 0) for i in range(len(p))]
    lower=max(F(0), T-paircap*comb(m,2),sum(max(F(0),m*w[i]-bj[i]*comb(m,2)) for i in range(len(p))))
    # For small A weights and a consistently oriented label, only pairs at
    # Hamming distance 2 contribute; exception: self-inverse, residual zero.
    if all(x[1]==1 and t[0]=='A' for (t,i),x in zip(p,I)) and (h>0 or any(t!=('A',1) for t,i in p)):
        q=[t[1]+1 for t,i in p]; ub=F(0); per=[F(0)]*len(p)
        for i,j in combinations(range(len(p)),2):
            v=max_pairs(m,min(q[i],q[j])); ub+=(cross[i]+cross[j])*v
            per[i]+=cross[i]*v;per[j]+=cross[j]*v
        lower=max(lower,T-ub,sum(max(F(0),m*w[i]-per[i]) for i in range(len(p))))
    return lower

@lru_cache(None)
def patterns(A):
    cnt=Counter(A); out=[]
    opts=[(t,i) for t in cnt for i,x in enumerate(slots(t)) if x[0]<=5]
    def go(p,i,w):
        if p and cap(p):out.append(p)
        for j in range(i,len(opts)):
            o=opts[j];t,k=o;ww=slots(t)[k][0]
            if w+ww<=5 and sum(x[0]==t for x in p)<cnt[t]:go(p+(o,),j,w+ww)
    go((),0,F(0));return out

def square_class(q):
    q=F(q); assert q>0
    n=q.numerator*q.denominator; out=1; p=2
    while p*p<=n:
        parity=0
        while n%p==0:n//=p;parity^=1
        if parity:out*=p
        p+=1
    return out*n

def setup(A,sc=None):
    cnt=Counter(A); ts=list(cnt);r=sum(t[1] for t in A);c=19-r
    # Variable x_{p,m}: number of distinct labels of this type carrying m lines.
    # Empty support counted as one-line variables (no deficit).
    PP=patterns(A)
    # If a zero-residual label meets a second label in an A component,
    # their orbit extrema must be -1 and +1. This requires one singleton
    # and one half-size class. The latter may occur in only one label.
    def too_low(d,T):
        x=F(365*d,19)-T
        return x>0 and x*x>F(70*19-16,19)*F(d*(19-d),19)
    def allowed(p,m):
        if sum(x[0] for x in info(p))!=5:return True
        for t,i in p:
            typ,d=t; w,delta,_,_,_=slots(t)[i]
            if typ!='A' or 2*delta==d+1:continue
            extra=0
            if delta==1 and (d+1)%2==0:
                # Use a separate pattern, with support intersection exactly
                # the specified factor. Ignoring other restrictions is safe.
                cc=Counter(x for x,j in p); cc[t]-=1
                for q in PP:
                    if (t,(d+1)//2-1) not in q:continue
                    qc=Counter(x for x,j in q);qc[t]-=1
                    if any(cc[z]+qc[z]>cnt[z]-int(z==t) for z in cnt):continue
                    extra=max(extra,F(d+1,2)*cap(q))
            if too_low(d,m*w+extra):return False
        return True
    V=[((),1)]+[(p,m) for p in PP for m in range(1,cap(p)+1) if allowed(p,m)]
    n=len(V); C=[Counter(t for t,i in p) for p,m in V]
    WS=[sum(x[0] for x in info(p)) for p,m in V]
    E=[[sum(slots(x)[i][0] for x,i in p if x==t)*m for p,m in V] for t in ts]
    dims=[cnt[t]*t[1] for t in ts]
    D=[deficit(p,m) if p else F(0) for p,m in V]
    rows=[];rhs=[]
    def add(row,b): rows.append(tuple(row));rhs.append(F(b))
    count=[m for p,m in V]; add(count,73);add([-x for x in count],-73)
    if sc is not None:
        assert c==1
        add([int(w<5 and square_class(5-w)!=sc) for w in WS],0)
    add([int(not p) for p,m in V], B(c))
    for a,b in combinations_with_replacement(ts,2):
        add([comb(cc[a],2) if a==b else cc[a]*cc[b] for cc in C],comb(cnt[a],2) if a==b else cnt[a]*cnt[b])
    for t in ts:
        a,d=t
        for key in sorted(set(x[2] for x in slots(t))):
            if key=='s':
                add([m*int(len(p)==1 and p[0][0]==t and p[0][1]==0) for p,m in V],cnt[t]*(c if d==1 else (4+8//(d-1))))
            else:
                add([sum(x==t and slots(x)[i][2]==key for x,i in p) for p,m in V],cnt[t])
    # At each component of a given type, all other components in its
    # nonsingleton supports are distinct. Apply integer packing floors.
    for t in ts:
        for u in ts:
            avail=cnt[u]-int(u==t)
            for k in range(2,6):
                if avail<k:continue
                add([cc[t]*int(cc[u]-int(u==t)>=k) for cc in C],cnt[t]*(avail//k))
    # Any three supports containing three specified factors of a five-set
    # would use >=6 distinct factors. This improves the degree bound on 5.
    for t in ts:
        if cnt[t]==5:add([int(cc[t]>=3) for cc in C],2)
    # A2 triple labels (all support factors of A2 type) have unit residuals.
    # If two triples are disjoint, all others meet both and on seven
    # factors there are at most three such triples. Otherwise the residuals
    # form equiangular lines with cosine 1/3 in dimension c.
    if cnt[('A',2)]==7 and c<5:
        add([int(len(p)==3 and all(t==('A',2) for t,i in p)) for p,m in V],max(5,8*c//(9-c)))
    if cnt[('A',1)]==1 and cnt[('A',2)]==7 and r==15:
        # Six pure triples leave a triangle of unused A2 pairs.
        # The A2 pairs in mixed A1+A2+A2 labels form a matching.
        add([2*int(len(p)==3 and cc[('A',2)]==3)+int(len(p)==3 and cc[('A',1)]==1 and cc[('A',2)]==2) for cc,(p,m) in zip(C,V)],13)
    if A==(("A",1),)*3+(("D",6),("D",10)):
        # For A1 + D6 vector + D10 vector, orient the A1 coordinate
        # positively. Fix a D6 coordinate and its sign; at most two
        # representatives share it (same D10 coordinate, two signs).
        add([int(p==((('A',1),0),(('D',6),0),(('D',10),0)) and m>24) for p,m in V],0)
    if A==(("A",5),)*3+(("D",4),):
        # For A5 half-class + D4 nonzero class, group by the D4
        # coordinate after orienting its sign positively. Each group has
        # at most four lines, or at most two if a complementary pair occurs.
        # Distinct groups cannot share an antipodal A5 triple direction.
        add([int(len(p)==2 and (('A',5),2) in p and any(t==('D',4) for t,i in p) and m>12) for p,m in V],0)
    comp=[(5-w)*m for w,(p,m) in zip(WS,V)]
    if c:E.append(comp);dims.append(c)
    else:add(comp,0)
    # Zero-residual lines lie in rank r.
    add([m*int(w==5) for w,(p,m) in zip(WS,V)],B(r))
    # Symmetrized subspace bounds: omit an arbitrary subset of irreducible factors.
    for a in product(*(range(cnt[t]+1) for t in ts)):
        d=sum(x*t[1] for x,t in zip(a,ts))
        if not d:continue
        coeff=[]
        for cc,(p,m) in zip(C,V):
            pr=F(m)
            for x,t in zip(a,ts):
                if x:pr*=F(comb(cnt[t]-cc[t],x),comb(cnt[t],x))
            coeff.append(pr)
        add(coeff,B(19-d))
        if c:add([v*int(w==5) for v,w in zip(coeff,WS)], B(r-d))
    add([4*x for x in D],70)
    # Tangent inequalities at scalar average along every block split.
    for sz in range(1,len(E)):
        for J in combinations(range(len(E)),sz):
            d=sum(dims[i] for i in J);row=[sum(E[i][j] for i in J) for j in range(n)]
            q=F(70*19-16,19)*F(d*(19-d),19); u=F(isqrt(q.numerator*10000//q.denominator)+1,100);mid=F(365*d,19)
            add(row,mid+u);add([-x for x in row],-mid+u)
    return V,rows,rhs,E,dims,D

def tangent(E,dims,D,g):
    row=tuple(4*x+sum(2*gi*Ei[j] for gi,Ei in zip(g,E)) for j,x in enumerate(D))
    b=70+sum(d*(38*x+x*x) for d,x in zip(dims,g))
    return row,b

def verify_dual(rows,rhs,cert):
    # Every variable is nonnegative, and their sum is <=73 (count equality).
    n=len(rows[0]); y=[(i,F(v)) for i,v in cert]; assert all(v>=0 for i,v in y)
    assert all(0<=i<len(rows) for i,v in y)
    val=sum(rhs[i]*v for i,v in y)
    lo=min(F(0),*(sum(rows[i][j]*v for i,v in y) for j in range(n)))
    return val<73*lo

def search(A,sc=None):
    import numpy as np
    from scipy.optimize import linprog
    V,rows,rhs,E,dims,D=setup(A,sc);n=len(V); gs=[]
    for k in range(60):
        mat=np.array(rows,dtype=float);b=np.array(rhs,dtype=float)
        # Phase I; common nonnegative relaxation variable s.
        res=linprog(np.r_[np.zeros(n),1.], A_ub=np.c_[mat,-np.ones(len(rows))],b_ub=b,bounds=(0,None),method='highs')
        if not res.success:return None,'solver'
        if res.fun>1e-8:
            for S in (10**6,10**9,10**12):
                y=[(i,str(F(round(-v*S),S))) for i,v in enumerate(res.ineqlin.marginals) if round(-v*S)>0]
                if verify_dual(rows,rhs,y):return {'g':[[str(x) for x in g] for g in gs], 'y':y},'excluded'
            return None,'rounding'
        t=np.array(E,float)@res.x[:-1];delta=t-19*np.array(dims);energy=sum(delta**2/np.array(dims))+4*(np.array(D,float)@res.x[:-1])
        if energy<=70+1e-7:
            return None,'relaxation'
        g=[F(round(x*100000),100000) for x in delta/np.array(dims)];row,b=tangent(E,dims,D,g)
        if row in rows:return None,'tangent repeat'
        rows.append(row);rhs.append(b);gs.append(g)
    return None,'iteration limit'

def key(A):return '+'.join(t+str(d) for t,d in A)
def parse(s):return tuple((x[0],int(x[1:])) for x in s.split('+'))
def verify(record):
    A=parse(record['type']);V,rows,rhs,E,dims,D=setup(A,record.get('square_class'))
    for g in record['g']:
        row,b=tangent(E,dims,D,[F(x) for x in g]);rows.append(row);rhs.append(b)
    assert verify_dual(rows,rhs,record['y']),record['type']

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify');ap.add_argument('--complete',action='store_true');ap.add_argument('--small',action='store_true');ap.add_argument('--output',default='mixed_deficit_certificates.jsonl');ap.add_argument('--max-types',type=int);ap.add_argument('--large-only',action='store_true');args=ap.parse_args()
    if args.verify:
        ordinary=set();branches={};stats=Counter();ncert=0
        for line in open(args.verify):
            rec=json.loads(line);verify(rec);ncert+=1
            t=rec['type'];sc=rec.get('square_class')
            if sc is None:
                assert t not in ordinary and t not in branches
                ordinary.add(t)
            else:
                assert t not in ordinary and sc not in branches.setdefault(t,set())
                branches[t].add(sc)
        for t,scs in branches.items():
            A=parse(t);assert sum(d for _,d in A)==18
            V,*_=setup(A)
            expected={square_class(5-sum(x[0] for x in info(p))) for p,m in V if sum(x[0] for x in info(p))<5}
            assert scs==expected,(t,scs,expected)
        seen=ordinary|set(branches)
        for t in seen:stats[sum(d for t,d in parse(t))]+=1
        expected={key(A) for A in gen()}
        print('VERIFIED certificates',ncert,'distinct root types',len(seen),'by rank',dict(sorted(stats.items())))
        print('FULL MIXED COVERAGE',expected<=seen,'missing',len(expected-seen))
        if args.complete:assert seen==expected,(len(seen),len(expected),expected-seen,seen-expected)
        raise SystemExit
    from time import monotonic
    import sys
    A0=list(gen());A0.sort(key=lambda A:(sum(d for t,d in A),A))
    if args.small:A0=[A for A in A0 if all(t=='A' and d<=3 for t,d in A)]
    if args.large_only:A0=[A for A in A0 if any(t!='A' or d>=4 for t,d in A)]
    stats=Counter();start=monotonic()
    with open(args.output,'w') as f:
        for j,A in enumerate(A0):
            if args.max_types and j>=args.max_types:break
            rec,why=search(A);stats[sum(d for t,d in A),why]+=1
            if rec:
                rec['type']=key(A);f.write(json.dumps(rec)+'\n');f.flush()
            else:print('SURV',key(A),why,flush=True)
            if (j+1)%25==0:print('PROGRESS',j+1,'excluded',sum(v for (r,w),v in stats.items() if w=='excluded'),'secs',round(monotonic()-start,1),file=sys.stderr,flush=True)
    print('STATS',sorted(stats.items()),flush=True)
