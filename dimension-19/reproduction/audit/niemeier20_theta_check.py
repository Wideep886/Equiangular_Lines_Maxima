#!/usr/bin/env python3
"""Independent shell-count check using Weyl group orders and the Niemeier theta series."""
import math, collections as co, sys
from niemeier20_independent import Local,TYPES,glue

def order(C):
 n=len(C); unseen=set(range(n)); ans=1
 while unseen:
  i=min(unseen);unseen.remove(i);q=[i]
  for j in q:
   for k in list(unseen):
    if C[j,k]==-1:unseen.remove(k);q.append(k)
  m=len(q);deg={j:sum(C[j,k]==-1 for k in q) for j in q}
  if max(deg.values(),default=0)<=2:ans*=math.factorial(m+1)
  else:
   b=[j for j in q if deg[j]==3];assert len(b)==1; b=b[0];arms=[]
   for j in q:
    if C[b,j]!=-1:continue
    p=b;k=j;v=1
    while deg[k]>1:
     nxt=[x for x in q if x!=p and C[k,x]==-1];assert len(nxt)==1;p,k=k,nxt[0];v+=1
    arms.append(v)
   arms.sort()
   if arms[:2]==[1,1]:ans*=2**(m-1)*math.factorial(m)
   else:assert arms==[1,2,m-4] and m in (6,7,8);ans*={6:51840,7:2903040,8:696729600}[m]
 return ans

def local_hist(x,scale):
 W=order(x.C);out=co.defaultdict(co.Counter); check=co.Counter()
 for N,c,u in x.dom:
  if N>6*x.den:continue
  z=[i for i,v in enumerate(u) if not v];s=order(x.C[z,:][:,z]);assert W%s==0;n=W//s
  out[c][N*(scale//x.den)]+=n
  if N<=4*x.den:check[N,c]+=n
 computed=co.Counter()
 for N,c,U,B in x.x:computed[N,c]+=len(U)
 assert computed==check
 return out

def run(ids):
 cache={}
 for idx in ids:
  L=[]
  for t in TYPES[idx]:
   if t not in cache:cache[t]=Local(t)
   L.append(cache[t])
  scale=math.lcm(*(x.den for x in L)); C=glue(L,scale);H=[local_hist(x,scale) for x in L];counts=co.Counter()
  for a in C:
   v={0:1}
   for h,c in zip(H,a):
    w=co.Counter()
    for i,m in v.items():
     for j,n in h[c].items():
      if i+j<=6*scale:w[i+j]+=m*n
    v=w
   counts.update(v)
  n2,n4,n6=(counts[d*scale] for d in [2,4,6]);assert n2==sum(len(x.root) for x in L)
  assert n4==179280-24*(n2-720);assert n6==16954560+252*(n2-720)
  assert all(k%(2*scale)==0 for k in counts)
  print(idx,TYPES[idx],n2,n4,n6,flush=True)
if __name__=='__main__':run(list(map(int,sys.argv[1:])) or range(20))
