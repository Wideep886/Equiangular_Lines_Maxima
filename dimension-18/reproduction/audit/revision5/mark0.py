"""Recompute the nine product blocks and three-dimensional center ranks of A3^8, mark 0."""
from itertools import product
from fractions import Fraction
from functools import lru_cache
from math import prod

B = ((2,1,1,1),(-1,2,1,-1),(-1,-1,2,1),(-1,1,-1,2))
h = ((-2,-2,-2,6),) + ((-3,1,1,1),)*4 + ((0,0,0,0),)*3

@lru_cache(None)
def local(word, position, polarization):
    forms = {}
    for a,b,c in product(range(-8+word,9,4),repeat=3):
        t = (a,b,c,-a-b-c)
        if sum(x*x for x in t)>64 or (position==0 and len(set(t[:3]))!=1):
            continue
        y = tuple(2*x-w for x,w in zip(t,polarization))
        if position==0:
            key=y
        else:
            key=tuple(tuple(sorted(v for v,w in zip(y,polarization) if w==t))
                      for t in sorted(set(polarization)))
        forms.setdefault(key,[]).append(t)
    return tuple((values,sum(x*x for x in values[0]),sum(x*w for x,w in zip(values[0],polarization)))
                 for values in forms.values())

def rank(rows):
    basis={}
    for row in rows:
        row=list(map(Fraction,row))
        for i,prev in sorted(basis.items()):
            row=[x-row[i]*y for x,y in zip(row,prev)]
        if any(row):
            i=next(i for i,x in enumerate(row) if x)
            basis[i]=[x/row[i] for x in row]
    return len(basis)

def cells(polarization=h):
    for tag in product(range(4),repeat=4):
        word=tag+tuple(sum(a*b for a,b in zip(tag,col))%4 for col in zip(*B))
        if sum(k*(4-k) for k in word)>16:
            continue
        pieces=[local(k,i,polarization[i]) for i,k in enumerate(word)]
        def complete(i,remain,pair,choice):
            if i==8:
                if remain==pair==0:yield choice
                return
            for values,q,p in pieces[i]:
                if q<=remain:
                    yield from complete(i+1,remain-q,pair-p,choice+[values])
        yield from complete(0,64,48,[])

def main():
    rows=[]
    for C in cells():
        y=[tuple(2*z-w for z,w in zip(g[0],h[i])) for i,g in enumerate(C)]
        assert y[0]==(0,0,0,0)
        center=tuple(-y[i][0] for i in range(1,5))
        if center[0]<0:continue
        assert center[0] and sum(center)==0
        variable={i:len(g) for i,g in enumerate(C) if len(g)>1}
        for i,g in enumerate(C):
            values=[tuple(2*z-w for z,w in zip(x,h[i])) for x in g]
            if i>=5:
                mean=[sum(x) for x in zip(*values)]
                assert mean==[0]*4
                if len(g)>1:
                    v=values[0];counts=sorted(v.count(a) for a in set(v))
                    assert counts==([1,3] if len(g)==4 else [2,2])
            if 1<=i<=4 and len(g)>1:
                assert len(g)==3 and len(set(values[0][1:]))==2
        rows.append((tuple(variable.items()),center))
    assert len(rows)==9 and sum(prod(n for _,n in L) for L,c in rows)==288
    groups=[({i for i,n in L if n==3},next(i for i,n in L if n==4),C)
            for L,C in rows if len(L)==3]
    small=[(next(i for i,n in L if n==6),next(i for i,n in L if n==4)) for L,C in rows if len(L)==2]
    print('large blocks (edge; colour; center coordinates)')
    for edge,k,c in sorted(groups,key=lambda v:sorted(v[0])):
        print(''.join(map(str,sorted(edge))),k,c)
    print('small blocks (octahedron, tetrahedron)', sorted(small))
    assert sorted({k for _,k,_ in groups})==[5,6,7]
    for k in (5,6,7):
        assert rank([c for _,j,c in groups if j!=k])==3
    costs={}
    edge,col,_=groups[0]
    for mask in range(1,1<<7):
        S={j+1 for j in range(7) if mask&(1<<j)}
        if S&(edge|{col}):continue
        a=len(S&{1,2,3,4});b=len(S)-a
        cost=0
        for U,k,_ in groups:
            count=len(S&(U|{k}))
            cost+=0 if not count else 2 if count==1 and k not in S else 1
        cost+=sum(len(S&{a,b})==1 for a,b in small)
        costs[a,b]=min(cost,costs.get((a,b),100))
    assert costs=={(0,1):4,(0,2):6,(1,0):6,(1,1):8,(1,2):8,(2,0):9,(2,1):9,(2,2):7}
    print('center ranks after dropping one colour: 3 3 3')
    print('minimum deficits by active vertices, colours: - 4 6 / 6 8 8 / 9 9 7')

if __name__=='__main__':
    main()
