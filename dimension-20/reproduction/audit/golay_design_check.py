#!/usr/bin/env python3
"""Exact finite inputs for the analytic A1^24 proof."""
from itertools import combinations
from collections import Counter
g=sum(1<<i for i in (0,1,5,6,7,9,11))
C=[0]
for i in range(12): C += [a^(g<<i) for a in C]
assert len(set(C))==4096
assert Counter(a.bit_count() for a in C)=={0:1,7:253,8:506,11:1288,12:1288,15:506,16:253,23:1}
C=[(a<<1)|(a.bit_count()%2) for a in C]
assert len(set(C))==4096
assert Counter(a.bit_count() for a in C)=={0:1,8:759,12:2576,16:759,24:1}
O=[tuple(i for i in range(24) if a>>i&1) for a in C if a.bit_count()==8]
F=Counter(T for S in O for T in combinations(S,5))
assert len(F)==42504 and set(F.values())=={1}
# The code is linear (parity extension of a binary span); all words have weight divisible by 4.
# This implies the corresponding Construction-A lattice is even and unimodular.
print('Golay code weight enumerator and S(5,8,24) verified')

# Small integer-partition bound in the octad-plus-coordinate proof.
from math import comb

def b(s):
    return min(comb(s,2)//3, s*((s-1)//2)//3, 35-7*s+comb(s,2))

def partitions(n, low=3):
    if n == 0:
        yield ()
    for s in range(low,n+1):
        for p in partitions(n-s,s):
            yield (s,)+p

maxima = [max(sum(b(s) for s in p)+2*b(z)
              for p in partitions(15-z) if len(p)>=2)
          for z in range(10)]
assert maxima == [18,14,11,12,10,14,13,18,16,18]
assert max(maxima) <= 18
assert max(b(s) for s in range(7)) == 4 and b(7)==7
print('Normal-column partition bounds verified:', maxima)
