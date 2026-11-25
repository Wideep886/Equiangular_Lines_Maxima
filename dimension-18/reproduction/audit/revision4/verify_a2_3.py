"""Independent local-coordinate table for A2^12, mark 3. No subspace or clique search."""
from collections import Counter,defaultdict
from itertools import permutations,product,combinations

def dot(a,b):return sum(x*y for x,y in zip(a,b))
def signed(x):return (x+1)%3-1
def rows(x):return tuple(sorted(set(permutations(x))))
zero=(0,0,0)
root=rows((3,0,-3))
mini={1:rows((1,1,-2)),-1:rows((-1,-1,2))}
h=(zero,)*5+((4,-2,-2),(-2,1,1),(-2,1,1),(-1,-1,2),(-1,-1,2),(-2,1,1))
def code():
    C=[[0 if i==j==0 else 1 if i==0 or j==0 else (0,1,-1,-1,1)[(i-j)%5]
        for j in range(6)] for i in range(6)]
    for u in product((-1,0,1),repeat=6):
        v=tuple(signed(dot(u,col)) for col in zip(*C))
        if v[-1]==0:yield u+v[:-1]
def shell():
    for i,j in combinations(range(11),2):
        for a,b in product(root, repeat=2):
            x=[zero]*11;x[i],x[j]=a,b;yield x
    for c in code():
        if sum(map(abs,c))==6:
            yield from product(*(mini[s] if s else [zero] for s in c))
def neg(p):return tuple(tuple(-z for z in row) for row in p)
def key(p):
    return tuple(tuple(tuple(sorted(x[j] for j in range(3) if wi[j]==q))
            for q in sorted(set(wi))) for x,wi in zip(p,h))
def main():
    table=defaultdict(set)
    for x in shell():
        if sum(dot(xi,wi) for xi,wi in zip(x,h))!=27:continue
        p=tuple(tuple(2*z-v for z,v in zip(r,s)) for r,s in zip(x,h))
        table[min(key(p),key(neg(p)))].add(p)
    big=[]; small=[]; folded=0
    for k,S in table.items():
        W=[p for p in S if key(p)==k]
        data=[]
        for j in range(11):
            opts={p[j] for p in W}
            if len(opts)>1:
                assert len(opts) in (2,3)
                assert all(sorted(v)==sorted(next(iter(opts))) for v in opts)
                diffs={tuple(a-b for a,b in zip(u,v)) for u,v in combinations(opts,2)}
                assert all(sorted(v)==[-6,0,6] for v in diffs)
                data.append((j,len(opts)))
        n=1
        for _,s in data:n*=s
        assert len(W)==n
        centers=[{dot(p[j],h[j]) for p in W} for j in range(5,11)]
        assert all(len(s)==1 for s in centers)
        c=[next(iter(s))//6 for s in centers]
        if n==18:
            assert len(S)==36
            (i,k1),(j,k2),(y,k3)=data
            assert j<5 and y>5 and [k1,k2,k3]==[3,3,2]
            assert c[0]==c[y-5]==0 and sorted(c[1:])==[-1,-1,0,1,1]
            t=c[1:]; s=next(v for v in t if v)
            big.append((i,j,y-5,tuple(v*s for v in t)))
        elif n==4:
            assert len(S)==8 and data[0]==(5,2) and data[1][1]==2 and c[0]!=0
            small.append(data[1][0]-5)
        else:
            assert n==32 and len(S)==32 and data==[(i,2) for i in range(6,11)]
            assert all(all(sum(p[j][k] for p in W)==0 for k in range(3)) for j in range(11))
            folded+=1
    assert folded==1 and sorted(small)==[1,2,3,4,5]
    assert sorted((i,j) for i,j,k,s in big)==list(combinations(range(5),2))
    for y in range(1,6):
        a,b=[(i,j,t) for i,j,k,t in big if y==k]
        assert len(set(a[:2]+b[:2]))==4 and a[2]!=b[2]
    print('edge, color, centered heights on the last five blocks')
    for i,j,k,t in sorted(big):print(str(i)+str(j),k,''.join('0+-'[z] for z in t))
    print('candidate axes',sum(map(len,table.values()))//2,'groups',len(table),'bound',10*5+5*2+4)
if __name__=='__main__':main()
