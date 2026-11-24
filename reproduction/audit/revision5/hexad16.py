"""A3^8 mark 16: all 24 centered cells, center ranks, and chamber-free anchor equivalence."""
from itertools import combinations, permutations, product
from fractions import Fraction
from mark0 import B, cells, rank

h=((0,)*4,)+((-3,1,1,1),)*4+((-2,-2,2,2),)*3

def tables():
    rows={}
    for entry in cells(h):
        means=tuple(tuple(Fraction(sum(2*x[k]-h[i][k] for x in values),len(values))
                          for k in range(4)) for i,values in enumerate(entry))
        c=tuple(z[0] for z in means)
        if next(a for a in c if a)>0:continue
        key=''
        for i,values in enumerate(entry):
            if len(values)==1:continue
            if i<=4:
                assert 1<=i<=4 and len(values)==3 and len({x[0] for x in values})==1
                key+=str(i)
            else:
                assert len(values)==2
                char='abc'[i-5]
                key+=char if values[0][0]!=values[1][0] else char.upper()
        assert key not in rows
        rows[key]=c
    return rows

def main():
    rows=tables()
    special='ABCabc'
    assert set(rows)==set(special)|{'%d%d'%(a,b) for a,b in combinations(range(1,5),2)}|set(
        '1AC 1Bc 1ab 2AB 2aC 2bc 3Ab 3BC 3ac 4Ac 4aB 4bC'.split())
    ns=[c for tag,c in rows.items() if len(tag)>1]
    assert rank(ns)==6 and all(rank(ns+[rows[t]])==7 for t in special)
    print('6 singleton blocks, 6 triangle-pair blocks, 12 mixed blocks')
    for i in range(1,5):print('triangle',i,sorted(k[1:] for k in rows if k[0]==str(i) and len(k)==3))
    print('nonsingleton centers span 6; any singleton raises the rank to 7')

    # Permutations of the four equal ternary factors and three binary factors,
    # with optional sign on each binary factor. The compensating local
    # coordinate permutations for a sign swap the two equal-coordinate pairs.
    generators=tuple(tuple(int(i==j) for j in range(4))+B[i] for i in range(4))
    code={v+tuple(sum(a*b for a,b in zip(v,col))%4 for col in zip(*B))
          for v in product(range(4),repeat=4)}
    def morph(t,p,s):return tuple(s[j]*t[p[j]]%4 for j in range(8))
    images=set();stabilizers=0
    for pp in permutations(range(1,5)):
        for qq in permutations(range(5,8)):
            p=(0,)+pp+qq
            for signs in product((-1,1),repeat=4):
                s=(signs[0],1,1,1,1)+signs[1:]
                if not all(morph(t,p,s) in code for t in generators):continue
                stabilizers+=1
                # Image of cell 1ab: new triangle i with p[i]=1; plus the two
                # output binary pairs belonging originally to factors 5,6.
                key=str(p.index(1))
                for j in range(5,8):
                    if p[j] in (5,6):
                        a='abc'[j-5];key+=a if s[j]==1 else a.upper()
                assert key in rows
                images.add(key)
    assert len(images)==12
    print('stabilizer actions:',stabilizers,'images of 1ab:',len(images))

if __name__=='__main__':
    main()
