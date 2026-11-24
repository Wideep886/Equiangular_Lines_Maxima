"""Direct coordinate check for the A4^6 norm-six marking number 26. Python 3.8+."""
from itertools import combinations, permutations, product
from math import gcd
from functools import reduce
from fractions import Fraction as Q

w = ((0,0,0,0,0),(-5,0,0,0,5),(-4,1,1,1,1),(-2,-2,-2,3,3),
     (-4,1,1,1,1),(-2,-2,-2,3,3))  # Coordinates in units of 1/5
# The first three coordinates in component zero must be equal.
letters = (('a',(0,3,4)), ('b',(1,1,2,3)), ('c',(2,1,2,3,4)),
           ('d',(3,0,1,2)), ('e',(3,3,4)), ('f',(4,1,2,3,4)),
           ('g',(5,0,1,2)), ('h',(5,3,4)))
shapes = ('c','de','f','gh','F g','cf','C d','ceg','dfh','bdg','bch','bef',
          'ad','ag','ab','aeh','')
shapes = tuple(s.replace(' ','') for s in shapes)
capacity = (1,2,1,2,4,4,4,6,6,7,6,6,2,2,2,3,1)
def residue(x): return (x+2)%5-2
def dot(a,b): return sum(x*y for x,y in zip(a,b))
def local(v): return sorted(set(permutations(v)))
def flatten(x):return sum(x,())
def invert(x):return tuple(tuple(u-v for u,v in zip(wi,xi)) for wi,xi in zip(w,x))
def scaled(row,k):
    return tuple(5*t-k for t in row)
def classes(a,b,c):
    return tuple(residue(a+b*i+c*i*i) for i in range(5))+(c,)
def key(x):
    rows=[list(row) for row in x]
    for _,(i,*inds) in letters:
        vals=sorted(rows[i][j] for j in inds)
        for j,v in zip(inds,vals): rows[i][j]=v
    return tuple(map(tuple,rows))
def canonical(x):return min(key(x),key(invert(x)))
def rank(mat):
    A=[list(map(Q,row)) for row in mat]
    n=0
    for j in range(len(A[0])):
        k=next((k for k in range(n,len(A)) if A[k][j]),None)
        if k is None:continue
        A[k],A[n]=A[n],A[k]
        for k in range(n+1,len(A)):
            t=A[k][j]/A[n][j]
            A[k]=[u-t*v for u,v in zip(A[k],A[n])]
        n+=1
    return n
def make_groups():
    adj={0:[((0,)*5,0)]}
    for k in (-2,-1,1,2):
        z=tuple((1 if k>0 else -1) if j<abs(k) else 0 for j in range(5))
        adj[k]=[(v,dot(v,v)) for v in local(scaled(z,k))]
    for r in [(1,0,0,0,-1),(1,1,0,-1,-1)]:
        adj[0].extend((v,dot(v,v)) for v in local(scaled(r,0)))
    groups={}
    for a,b,c in product(range(-2,3),repeat=3):
        word=classes(a,b,c)
        if sum(abs(k)*(5-abs(k)) for k in word)>20:continue
        blocks=[]
        for i,(k,wi) in enumerate(zip(word,w)):
            blocks.append([(v,en,dot(v,wi)) for v,en in adj[k] if i or v[0]==v[1]==v[2]])
        # Exact suffix energy and scalar-product possibilities.
        suffix=[{(0,0)}]
        for block in reversed(blocks):
            suffix.append({(s+en,t+d) for s,t in suffix[-1] for _,en,d in block if s+en<=100})
        def rec(i,u,en,d):
            if (en,d) not in suffix[6-i]:return
            if i==6:
                x=tuple(u); g=canonical(x)
                groups.setdefault(g,set()).add(x)
                return
            for v,e,p in blocks[i]:rec(i+1,u+[v],en-e,d-p)
        rec(0,[],100,75)
    return groups
def main():
    groups=make_groups()
    assert len(groups)==17
    # Representatives. The same nonzero word gives the first four cells,
    # distinguished by which one of its four local pairings is submaximal.
    words=((0,-2,2),)*4+((0,2,-1),(0,-1,2),(0,1,0),(0,-2,-2),(0,0,1),
       (0,-1,-2),(0,1,1),(0,2,0),(-1,0,0),(-1,2,-2),(-1,1,2),(-1,1,1),(-2,-2,0))
    chosen=[]; centers=[]; sizes=[]
    for i,(abc,shape) in enumerate(zip(words,shapes)):
        word=classes(*abc)
        opts=[]
        for k,G in groups.items():
            pts=[x for x in G if tuple(residue(-row[0]) for row in x)==word]
            if i<4:pts=[x for x in pts if dot(w[2+i],x[2+i])<dot(w[2+i],w[2+i])]
            if pts:opts.append((k,pts))
        assert len(opts)==1,(i,len(opts))
        k,G=opts[0];chosen.append(k)
        assert len(groups[k])==2*len(G)
        sizes.append(len(G))
        expected={'a':2,'b':3,'c':4,'d':3,'e':2,'f':4,'g':3,'h':2,'C':6,'F':6}
        count=1
        for name,(n,*J) in letters:
            vals={tuple(x[n][j] for j in J) for x in G}
            letter=name.upper() if name.upper() in shape else name
            assert len(vals)==(expected[letter] if letter in shape else 1),(i,name)
            if len(vals)>1:
                assert all(len(set(v))==2 for v in vals)
                assert all(sorted(v)==sorted(next(iter(vals))) for v in vals)
                count*=len(vals)
        assert count==len(G)
        c=[sum(2*flatten(x)[j]-flatten(w)[j] for x in G) for j in range(30)]
        centers.append(c)
    assert len(set(chosen))==17 and sum(capacity)==59
    assert rank([centers[i] for i in (1,3,4,5,12)])==5
    assert rank([centers[i] for i in (1,3,4,5,7,12)])==6
    print('The 17 cell sizes:',sizes,'total',sum(sizes))
    print('Cap total',sum(capacity),'six occupied centers independent')
if __name__=='__main__':main()
