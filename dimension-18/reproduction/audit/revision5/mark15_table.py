"""Recover the constant and moving parts of row 15 from the octacode shell."""
from collections import Counter
from fractions import Fraction
from itertools import permutations
from math import prod
from pathlib import Path
import sys
sys.path.insert(0,str(Path(sys.argv[1]).resolve()))
from octacode import cells, shells, vertices
h=((0,0,0,0),(-4,0,0,4),(-2,-2,2,2),(-2,-2,2,2),
   (0,0,0,0),(0,0,0,0),(-2,-2,2,2),(-2,-2,2,2))
w=sum(h,())
# Capitals A,B,C,D and their lowercase mates label different binary blocks.
blocks={'X':(1,(1,2)),'S':(4,(0,1,2,3)),'T':(5,(0,1,2,3))}
for label,i in zip('ABCD',(2,3,6,7)):
    blocks[label]=i,(0,1)
    blocks[label.lower()]=i,(2,3)
centers={}
for factors in cells(h,shells(4)):
    pts=[tuple(2*x-a for x,a in zip(row,w)) for row in vertices(factors)]
    n=len(pts)
    center=tuple(Fraction(sum(col),n) for col in zip(*pts))
    if next((x for x in center if x),0)>0:continue
    labels=[]
    for letter,(i,idx) in blocks.items():
        local={tuple(p[4*i+j] for j in idx) for p in pts}
        if len(local)<=1:continue
        t=next(iter(local))
        assert len(set(t))==2 and max(t)-min(t)==8 and set(permutations(t))==local
        if len(local)==6:
            assert letter in ('S','T')
            labels.append('P'+letter)
        else:labels.append(letter)
    names=tuple(sorted(labels))
    assert n==prod(6 if k.startswith('P') else len(blocks[k][1]) for k in names)
    assert names not in centers
    centers[names]=center
def name(s):return tuple(sorted(s))
def names(ss):return {name(s) for s in ss.split()}
short='Ad aD BC bc Aa Bb Cc Dd'.split()
assert set(centers)==names('Ad aD BC bc Aa Bb Cc Dd ACT AbS BdT CDS aBS acT cdS bDT XABcD XabCd')|{('PS','X'),('PT','X')}
for j in ('XABcD','XabCd'):
    assert not any(centers[name(j)])
for labels in [('Aa','Ad','aD','Dd'),('Bb','BC','bc','Cc')]:
    print('four-cycle',*labels)
for letter in ('S','T'):
    print('tetrahedron',letter,sorted(''.join(k for k in lab if k!=letter) for lab in centers if letter in lab))
small=[centers[name(s)] for s in short]
def det(a):
    return sum((-1)**i*c*det([s[:i]+s[i+1:] for s in a[1:]]) for i,c in enumerate(a[0])) if a else 1
def relation(v):
    assert not any(sum(a*b for a,b in zip(col,v)) for col in zip(*small))
relation((-1,1,0,0,0,-1,1,0))
relation((0,0,1,-1,-1,0,0,1))
# Search only for a displayed nonzero six-by-six minor.
from itertools import combinations
for cols in combinations((0,4,5,8,12,24,28),6):
    m=[[v[j]/2 for j in cols] for v in small[:6]]
    delta=det(m)
    if delta:break
assert delta
print('first six centers divided by two, columns',*cols)
for row in m:
    assert all(a.denominator==1 for a in row)
    print(*map(int,row))
print('determinant',delta)
