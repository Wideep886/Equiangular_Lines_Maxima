"""Reconstruct the affine plane and the centers of every cell from the ternary code."""
from collections import defaultdict
from itertools import combinations, permutations, product

def dot(x, y):
    return sum(a*b for a, b in zip(x, y))

def codewords():
    def signed(x): return (x+1)%3-1
    C = [[0 if i==j==0 else 1 if i==0 or j==0 else (0,1,-1,-1,1)[(i-j)%5]
          for j in range(6)] for i in range(6)]
    for u in product((-1,0,1), repeat=6):
        v = tuple(signed(dot(u, col)) for col in zip(*C))
        if v[-1] == 0:
            yield u + v[:-1]

def shell(words, mini):
    zero = (0,0,0)
    roots = set(permutations((-3,0,3)))
    for i,j in combinations(range(11),2):
        for a,b in product(roots,repeat=2):
            v=[zero]*11
            v[i],v[j]=a,b
            yield v
    for c in words:
        if sum(map(abs,c))==6:
            yield from product(*(mini[s] if s else [zero] for s in c))

def check():
    words = list(codewords())
    assert len(words)==243
    mini = {s:sorted(set(permutations((s,s,-2*s)))) for s in (-1,1)}
    c = next(w for w in words if sum(map(abs,w))==9)
    order = [i for i,s in enumerate(c) if not s] + [i for i,s in enumerate(c) if s]
    h = tuple((0,0,0) if s==0 else mini[s][0] for s in c)
    cells = defaultdict(list)
    signed=set()
    oriented=set()
    for v in shell(words,mini):
        if sum(dot(p,q) for p,q in zip(v,h))!=27: continue
        raw = [tuple(2*a-b for a,b in zip(v[i],h[i])) for i in order]
        point=tuple(raw)
        assert point not in signed
        signed.add(point)
        heights = tuple(dot(r,h[i]) for r,i in zip(raw[2:],order[2:]))
        assert all(y%6==0 for y in heights)
        center = tuple(y//6 for y in heights)
        active = tuple(i for i,y in enumerate(center) if y==0)
        if len(active)==1:
            j, = [t for t,r in enumerate(raw[:2]) if any(r)]
            if sorted(raw[j])!=[-4,2,2]: continue
            cell = (j,active[0])
        else:
            assert len(active)==3 and not any(any(r) for r in raw[:2])
            if next(y for y in center if y)<0:continue
            cell = (2,active)
        cells[cell].append((raw,center))
        oriented.add(point)
    assert len(cells)==30
    def opposite(point):return tuple(tuple(-z for z in r) for r in point)
    assert len(signed)==408 and len(oriented)==204
    assert signed==oriented|{opposite(x) for x in oriented}
    centers={}
    triangles={}
    steps={}
    for cell,rows in cells.items():
        cc={c for raw,c in rows}
        assert len(cc)==1
        centers[cell]=cc.pop()
        fields=[{raw[j] for raw,c in rows} for j in range(11)]
        data=[(j,len(f)) for j,f in enumerate(fields) if len(f)>1]
        for i,f in enumerate(fields):
            sums=tuple(sum(raw[i][m] for raw,_ in rows) for m in range(3))
            if i<2:
                assert sums==(0,0,0)
                if len(f)>1:assert triangles.setdefault(i,f)==f
            else:
                assert sums==tuple(len(rows)*centers[cell][i-2]*v for v in h[order[i]])
                if len(f)==2:
                    x,y=sorted(f)
                    diff=tuple(a-b for a,b in zip(x,y))
                    assert steps.setdefault(i,diff)==diff and dot(diff,h[order[i]])==0
        j, key=cell
        assert data==([(j,3),(2+key,2)] if j!=2 else [(2+i,2) for i in key])
        assert len(rows)==(6 if j!=2 else 8)
        for i,_ in data:
            for a,b in combinations(fields[i],2):
                assert sorted(x-y for x,y in zip(a,b))==[-6,0,6]
    assert len(steps)==9 and len(triangles)==2
    lines=[set(s) for (j,s) in centers if j==2]
    pencils=[]
    for line in lines:
        family=[line]+[l for l in lines if line.isdisjoint(l)]
        assert len(family)==3 and len(set.union(*family))==9
        if not any(line in p for p in pencils):pencils.append(family)
    assert len(pencils)==4
    assert sorted(tuple(sorted(ab)) for line in lines for ab in combinations(line,2))==list(combinations(range(9),2))
    for family in pencils:
        for i,line in enumerate(family):
            c=centers[2,tuple(sorted(line))]
            h1=tuple(int(j in family[(i+1)%3])-int(j in family[(i+2)%3]) for j in range(9))
            assert c==h1 or c==tuple(-y for y in h1)
    signs=[]
    for j in (0,1):
        eps=[]
        for family in pencils:
            values={centers[j,i][k] for L in family for i in L for k in L if i!=k}
            assert len(values)==1
            eps.append(values.pop())
        assert sorted(eps)==[-1,-1,1,1]
        signs.append(eps)
        for i in range(9):
            L=[next(l for l in family if i in l) for family in pencils]
            assert centers[j,i]==tuple(sum(s for s,l in zip(eps,L) if k in l) for k in range(9))
    assert len(set(zip(*signs)))==4
    target={frozenset(l) for l in combinations(range(9),3) if all(sum((i//3,i%3)[k] for i in l)%3==0 for k in (0,1))}
    iso=next(f for g in permutations(range(1,9)) for f in [(0,)+g] if {frozenset(f[i] for i in l) for l in lines}==target)
    print('204 axes; 18 rectangles and 12 affine-line cubes; all cell products, steps and centers verified.')
if __name__=='__main__':check()
