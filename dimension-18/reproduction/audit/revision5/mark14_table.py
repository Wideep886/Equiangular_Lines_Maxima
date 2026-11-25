"""Regenerate the quotient proof table for row 14 of the eight-copy A3 shell."""
import sys
from fractions import Fraction
from itertools import combinations, permutations, product
from math import prod
from pathlib import Path

sys.path.insert(0,str(Path(sys.argv[1]).resolve()))
from octacode import cells, shells, vertices

h=((0,0,0,0),(-4,0,0,4),(-3,1,1,1),(-3,1,1,1),
   (-2,-2,2,2),(0,0,0,0),(-3,1,1,1),(-1,-1,-1,3))
w=sum(h,())
# Each letter specifies its coordinate block in 4w, numbered from zero.
positions={'A':(2,(1,2,3)),'B':(3,(1,2,3)),'C':(6,(1,2,3)),
           'D':(7,(0,1,2)), 'X':(1,(1,2)), 'Y':(4,(0,1)),
           'Z':(4,(2,3)), 'S':(5,(0,1,2,3))}
centers={}
for factors in cells(h,shells(4)):
    points=[tuple(2*x-a for x,a in zip(row,w)) for row in vertices(factors)]
    mean=tuple(Fraction(sum(col),len(points)) for col in zip(*points))
    if next(x for x in mean if x)>0:continue  # choose one orientation per orbit
    labels=[]
    for label,(i,indices) in positions.items():
        options={tuple(point[4*i+j] for j in indices) for point in points}
        if len(options)==1:continue
        token=next(iter(options))
        assert set(permutations(token))==options
        assert len(set(token))==2 and max(token)-min(token)==8
        two_positions=sorted([token.count(value) for value in set(token)])==[2,2]
        labels.append('P' if two_positions else label)
        assert not two_positions or label=='S'
    assert len(points)==prod(6 if label=='P' else len(positions[label][1]) for label in labels)
    tag=''.join(sorted(labels))
    x0=mean[0];x1=mean[7]-mean[4]
    group=0 if x0 else 1 if x1 else 2
    assert not(x0 and x1)
    assert (group,tag) not in centers
    centers[group,tag]=mean
for group in range(3):
    print('center group',group,sorted(tag for j,tag in centers if j==group))
assert set(centers)=={(0,t) for t in 'A B C D XZ'.split()}|{(1,t) for t in 'A B C D YZ'.split()}|{(2,t) for t in 'AS BS CS DS ABX CDX ACZ BDZ ADY BCY PZ SXY'.split()}
# Scaled values of the first coordinates in factors 1,2,3,4,6.
cols=(4,8,12,16,24)
scales=(Fraction(1,2),1,1,Fraction(1,2),1)
table=[tuple(centers[2,tag][i]*k for i,k in zip(cols,scales)) for tag in ('BDZ','ADY','ACZ','BCY','CDX')]
def det(a):
    if not a: return 1
    return sum((-1)**i*c*det([row[:i]+row[i+1:] for row in a[1:]]) for i,c in enumerate(a[0]))
assert abs(det(table))==128
print('minor:')
for row in table:
    assert all(v.denominator==1 for v in row)
    print(' '.join(str(int(v)) for v in row))
print('determinant',det(table))
