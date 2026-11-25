"""Check five Paley generators and the 12+18 cells of the nonad without scanning norm shells."""
from itertools import product
from collections import Counter
pts=list(product(range(3),repeat=2))
C=[[0 if i==j==0 else 1 if i==0 or j==0 else (0,1,-1,-1,1)[(i-j)%5] for j in range(6)] for i in range(6)]
order=(4,9,6,2,5,1,0,7,3,8,10)
sign=(1,1,-1,-1,-1,-1,-1,1,-1,1,-1)
linear=[(-1,-1,-1,-1,0,-1),(-1,-1,0,1,0,0),(-1,0,-1,0,0,1)]
quad=[(1,-1,-1,-1,1,-1),(-1,1,-1,1,0,-1)]
values=[[1,x,y,x*x+y*y,2*x*x+x*y+y*y] for x,y in pts]
for i,u in enumerate(linear+quad):
    w=u+tuple(sum(x*y for x,y in zip(u,c)) for c in zip(*C))
    assert w[-1]%3==0
    ww=tuple(z*w[j]%3 for j,z in zip(order,sign))
    assert ww==(int(i==3),int(i==4))+tuple(v[i]%3 for v in values)
    if i>=3:assert Counter(ww[2:])=={0:1,1:4,2:4}
cube=set();rect=set()
for a,b,c,s,t in product(range(3),repeat=5):
    f=tuple(sum(q*v for q,v in zip((a,b,c,s,t),row))%3 for row in values)
    z=int(bool(s))+int(bool(t));r=f.count(1)
    if z+9-f.count(0)!=6 or (6-z+r)<9:continue
    assert 6-z+r==9
    if z==0:cube.add(tuple(i for i,w in enumerate(f) if w==2))
    elif s+t==1:rect.add((int(bool(t)),f.index(2)))
assert len(cube)==12 and len(rect)==18 and len(cube)*8+len(rect)*6==204
print('Five Paley basis rows and their affine quadratic evaluations agree; 12 cube lines, 18 rectangles, 204 axes.')
