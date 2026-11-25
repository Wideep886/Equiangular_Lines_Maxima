"""The 64 words in the shortened octacode give the tournament for A3^8 #5."""
from itertools import product
from collections import Counter
P = ((2,1,1,1),(-1,2,1,-1),(-1,-1,2,1),(-1,1,-1,2))
h = (3,1,1,1,3,1,1,1)
arcs = []
for b,c,d in product(range(4), repeat=3):
    u=(0,b,c,d)
    w=u+tuple(sum(x*y for x,y in zip(col,u))%4 for col in zip(*P))
    matching = [i for i in range(8) if w[i]==h[i]]
    opposite = [i for i in range(8) if w[i]==-h[i]%4]
    even = [i for i,a in enumerate(w) if a==2]
    if len(matching)==3 and len(opposite)==len(even)==1:
        arcs.append((even[0],opposite[0]))
assert len({tuple(sorted(e)) for e in arcs})==len(arcs)==21
for i in range(1,8):
    adj=sorted(v for u,v in arcs if u==i)
    assert len(adj)==3
    print(i,*adj)
