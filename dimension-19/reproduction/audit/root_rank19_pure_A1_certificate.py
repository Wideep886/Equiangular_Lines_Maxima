#!/usr/bin/env python3
"""Exact finite checks for pure A1^r root systems, 73 lines in dimension 19.
Standard library only. All inequalities are rational and ranks over F_3.
The reductions and formulas are explained in the accompanying conversation.
"""
from fractions import Fraction as F
from itertools import combinations, product

N = 73

def contribution(k, q):
    # Minimum possible contribution from q quintets, with k lines in total.
    return min(sum(-4*a*(a-1) for a in aa)
               for aa in product(range(1,6), repeat=q) if sum(aa)==k)

def q1_min(r):
    s,c=r-5,19-r
    best=None
    for a in range(1,6):
        for t in range(min(4*min(s*(s-1)//2,5*(s//2)),N-a)+1):
            for d in range(N-a-t+1):
                if (d+1)//2+2*((t+3)//4)>5*s:
                    continue
                U=d+2*t
                L=(72*r+292-4*a*a-8*U+F((4*a-12)**2,5)
                   +F((U-17*s)**2,s)+F((19*r-69-4*a-U)**2,c))
                best=L if best is None else min(best,L)
    return best

def q2_min(r):
    s,c=r-9,19-r
    best=None
    for k in range(2,11):
        qcon=contribution(k,2)
        for p3 in range(4*(s//2)+1):
            for p2 in range(2*s+1):
                if (p2+1)//2+2*((p3+3)//4)>s:
                    continue
                for ab3 in range(min(64,16*s,N-k-p2-p3)+1):
                    for ab2 in range(min(32,N-k-p2-p3-ab3)+1):
                        if (ab2+1)//2+(ab3+3)//4>16:
                            continue
                        p1=N-k-p2-p3-ab2-ab3
                        T0=k+p1+p2+p3
                        T1=4*k+2*(ab2+ab3)
                        T2=p2+2*p3+ab3
                        Q=T0+T1+T2
                        L=(72*r+qcon+4*p1-4*(p2+ab2)-12*(p3+ab3)
                           +(T0-17)**2+F((T1-136)**2,8)+F((T2-17*s)**2,s)
                           +F((4+19*r-Q)**2,c))
                        best=L if best is None else min(best,L)
    return best

def q3_distinct_min(r):
    s,c=r-12,19-r
    best=None
    for k in range(3,16):
        qcon=contribution(k,3)
        for d in range(19):
            for t in range(4*min(9,3*s)+1):
                u=N-k-d-t
                if not 0<=u<=36 or (d+1)//2+(t+3)//4>9:
                    continue
                T0=2*k+d+t
                T1=3*k+d+t+3*u
                T2=t
                Q=T0+T1+T2
                L=(72*r+qcon-4*d-12*(t+u)+F((T0-51)**2,3)
                   +F((T1-153)**2,9)+(F((T2-17*s)**2,s) if s else 0)
                   +F((4+19*r-Q)**2,c))
                best=L if best is None else min(best,L)
    return best

def q3_concurrent_min():
    best=None
    for k in range(3,16):
        qcon=contribution(k,3)
        for t in range(min(64,N-k)+1):
            p=N-k-t
            T0=k+p
            T1=4*k+3*t
            Q=T0+T1
            L=(72*13+qcon+4*p-12*t+(T0-17)**2+F((T1-204)**2,12)
               +F((4+19*13-Q)**2,6))
            best=L if best is None else min(best,L)
    return best

# Rank 11 without quintets: at least 13 colored edges of K_6,
# each of five colors a matching of size <= 3. The nonzero residual
# inner products occur exactly between disjoint edges of distinct colors.
E=list(combinations(range(6),2))
EM=[sum(1<<i for i in e) for e in E]
MATCH={k:[m for m in combinations(range(15),k)
          if sum(EM[i].bit_count() for i in m)==
          (sum(EM[i] for i in m)).bit_count()]
       for k in (1,2,3)}
# The popcount equality above holds iff the edges are vertex-disjoint;
# to make this independent of that arithmetic observation, verify explicitly.
for k in MATCH:
    assert MATCH[k]==[m for m in combinations(range(15),k)
                      if len({v for e in m for v in E[e]})==2*k]

def matchings(sz):
    def rec(A,used):
        if len(A)==5:
            yield tuple(A)
            return
        i=len(A);k=sz[i]
        for m in MATCH[k]:
            # Normalize one full matching by S_6; sort other equal sizes.
            if i==0 and m!=MATCH[3][0]:
                continue
            if i>1 and k==sz[i-1] and m<A[-1]:
                continue
            mask=sum(1<<j for j in m)
            if not mask&used:
                yield from rec(A+[m],used|mask)
    yield from rec([],0)

def rank3(mat):
    a=[list(row) for row in mat]
    n=len(a);k=0
    for j in range(n):
        p=next((i for i in range(k,n) if a[i][j]%3),None)
        if p is None: continue
        a[k],a[p]=a[p],a[k]
        a[k]=[(x*a[k][j])%3 for x in a[k]]
        for i in range(k+1,n):
            v=a[i][j]%3
            if v:
                a[i]=[(x-v*y)%3 for x,y in zip(a[i],a[k])]
        k+=1
    return k

def check_coloring(A):
    labels=[(i,e) for i,m in enumerate(A) for e in m]
    n=len(labels)
    edges=[]
    adj=[set() for _ in range(n)]
    for a,b in combinations(range(n),2):
        i,e=labels[a];j,f=labels[b]
        if i!=j and not EM[e]&EM[f]:
            edges.append((a,b));adj[a].add(b);adj[b].add(a)
    # Switch signs so all edges in a spanning forest are positive.
    seen=set();forest=set()
    for a in range(n):
        if a in seen: continue
        seen.add(a);stack=[a]
        while stack:
            u=stack.pop()
            for v in sorted(adj[u]):
                if v not in seen:
                    forest.add(tuple(sorted((u,v))))
                    seen.add(v);stack.append(v)
    free=[e for e in edges if e not in forest]
    index={e:i for i,e in enumerate(free)}
    eq=[]
    for tri in combinations(range(n),3):
        pairs=list(combinations(tri,2))
        if all(b in adj[a] for a,b in pairs) and any(
             len(adj[v].intersection(tri))==1 for v in range(n) if v not in tri):
            # Such a triangle must have positive sign product by PSD.
            eq.append(sum(1<<index[e] for e in pairs if e in index))
    minrank=n;count=0
    for mask in range(1<<len(free)):
        if any((mask&e).bit_count()%2 for e in eq):
            continue
        count+=1
        mat=[[0]*n for _ in range(n)]
        for i in range(n): mat[i][i]=2
        for a,b in edges:
            e=(a,b)
            mat[a][b]=mat[b][a]=2 if e in index and mask>>index[e]&1 else 1
        rk=rank3(mat)
        assert rk>=9,(A,mask,rk)
        minrank=min(minrank,rk)
    assert count>0
    return count,minrank

def main():
    for r in (11,12,13):
        for q,fun in [(1,q1_min),(2,q2_min)]:
            v=fun(r)
            print(f'r={r} q={q}: {v}')
            assert v is not None and v>70
    for r in (12,13):
        v=q3_distinct_min(r)
        print(f'r={r} q=3 distinct: {v if v is not None else "empty"}')
        assert v is None or v>70
    v=q3_concurrent_min()
    print(f'r=13 q=3 concurrent: {v}')
    assert v>70
    for sz in [(3,3,3,3,3),(3,3,3,3,2),(3,3,3,3,1),(3,3,3,2,2)]:
        records=signs=0;minrank=15
        for A in matchings(sz):
            records+=1
            count,rk=check_coloring(A)
            signs+=count;minrank=min(minrank,rk)
        assert records>0 and minrank>=9
        print(f'matching sizes {sz}: {records} records, {signs} signs, minimum rank mod 3 = {minrank}')
    print('PASS')

if __name__=='__main__':
    main()
