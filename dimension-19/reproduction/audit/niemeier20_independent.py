#!/usr/bin/env python3
"""Exact independent fiber bounds for marked Niemeier sections."""
import sys, math, itertools as it, collections as co
import numpy as np
import sympy as sp
from niemeier_19_certificate import omega

def cartan(t):
 a,n=t; C=2*np.eye(n,dtype=np.int64)
 edges=[(i,i+1) for i in range(n-1)] if a=='A' else ([(i,i+1) for i in range(n-2)]+[(n-3,n-1)] if a=='D' else [(i,i+1) for i in range(n-2)]+[(2,n-1)])
 for i,j in edges:C[i,j]=C[j,i]=-1
 return C

class Local:
 def __init__(self,t):
  self.t=t;self.n=t[1];self.C=cartan(t);S=sp.Matrix(self.C);self.den=int(S.det());self.K=np.array(S.inv()*self.den).astype(np.int64)
  assert self.den>0 and np.array_equal(self.C@self.K,self.den*np.eye(self.n,dtype=np.int64));assert np.all(self.K>0);assert self.den<=25 and np.max(self.K)<200
  self.dom=[]
  def visit(i,u,N):
   if i==self.n:
    u=tuple(u);b=self.K@u;res=tuple(int(a%self.den) for a in b);self.dom.append((int(N),res,u));return
   cross=2*sum(int(self.K[i,j])*u[j] for j in range(i));a=int(self.K[i,i]);k=0
   while N+cross*k+a*k*k<=7*self.den:visit(i+1,u+[k],N+cross*k+a*k*k);k+=1
  visit(0,[],0)
  R={}
  for N,r,u in self.dom:R[r]=min(N,R.get(r,N))
  self.res=sorted(R);assert len(self.res)==self.den
  self.min=[R[r] for r in self.res];idx={r:i for i,r in enumerate(self.res)}
  self.add=[[idx[tuple((x+y)%self.den for x,y in zip(a,b))] for b in self.res] for a in self.res]
  self.dom=[(N,idx[r],u) for N,r,u in self.dom]
  self.x=[];self.root=[];neighbors=[np.flatnonzero(self.C[j]==-1).tolist() for j in range(self.n)]
  for N,c,u in self.dom:
   if N>4*self.den:continue
   seen={u};queue=[u]
   for v in queue:
    for j,s in enumerate(v):
     if s:
      z=list(v);z[j]-=2*s
      for k in neighbors[j]:z[k]+=s
      z=tuple(z)
      if z not in seen:seen.add(z);queue.append(z)
   U=np.array(queue,dtype=np.int64);B=U@self.K;assert np.all(np.sum(B*U,axis=1)==N)
   self.x.append((N,c,U,B))
   if N==2*self.den and c==0:self.root.extend((B//self.den).tolist())
  self.root=np.array(self.root,dtype=np.int64);self.positive=self.root[np.all(self.root>=0,axis=1)];assert len(self.root)==2*len(self.positive);assert np.max(np.abs(self.root))<=8
  self.byclass=co.defaultdict(list)
  for i,(N,c,u) in enumerate(self.dom):self.byclass[c].append((N,i))
  self.cache={}
 def markings(self,i):
  u=self.dom[i][2];zero=set(j for j,x in enumerate(u) if not x);out=[]
  while zero:
   j=min(zero);out.append(j);zero.remove(j);q=[j]
   for k in q:
    for l in list(zero):
     if self.C[k,l]==-1:zero.remove(l);q.append(l)
  return out
 def options(self,i,mark=-1):
  key=(i,mark)
  if key in self.cache:return self.cache[key]
  W=np.array(self.dom[i][2],dtype=np.int64);wB=W@self.K
  roots=self.positive[self.positive@W==0]
  if mark>=0:roots=roots[(roots@self.C)[:,mark]==0]
  rootset={tuple(x) for x in roots};simples=[]
  for a in sorted(rootset,key=lambda x:sum(x)):
   if not any(tuple(x-y for x,y in zip(a,b)) in rootset for b in rootset if sum(b)<sum(a)):simples.append(a)
  R=np.array(simples,dtype=np.int64).reshape(-1,self.n)
  if len(R):
   M=sp.Matrix(R@self.C@R.T);det=int(M.det());assert 0<det<=2**24;inv=np.array(M.inv()*det).astype(np.int64);assert np.max(np.abs(inv))<=2**23;T=R.T@inv
  else:T=None
  out=co.defaultdict(list)
  for N,c,U,B in self.x:
   if mark>=0:
    mask=U[:,mark]==0;U=U[mask];B=B[mask]
   if not len(U):continue
   dot=U@wB;wnorm=self.dom[i][0];mask=(3*self.den-dot)**2<=(4*self.den-N)*(6*self.den-wnorm)
   U=U[mask];B=B[mask];dot=dot[mask]
   if not len(U):continue
   V=U.copy()
   for _ in range(1000):
    changed=False
    for a in R:
     vals=V@a;mask=vals<0
     if mask.any():V[mask]-=vals[mask,None]*(a@self.C)[None,:];changed=True
    if not changed:break
   else:raise AssertionError('reflection descent failed')
   groups=co.defaultdict(list)
   for k,v in enumerate(V):groups[tuple(v)].append(k)
   for v,ix in groups.items():
    A=U[ix];D=B[ix];assert len(set(dot[ix]))==1
    b=self.fiber_bound(N,A,D,T);out[c].append((N,int(dot[ix[0]]),len(ix),b))
  self.cache[key]=out;return out
 def fiber_bound(self,N,A,D,T):
  n=len(A)
  if n<=1:return n
  if n>500:return n
  P=A@D.T;diff=2*N-2*P;G=(diff==4*self.den)|(diff==6*self.den)|(diff==10*self.den)
  if T is not None:
   coords=A@T
   for j in range(n):G[j]&=~(np.all(coords[j]>=coords,axis=1)|np.all(coords[j]<=coords,axis=1))
  assert np.array_equal(G,G.T)
  if n<=120:return omega(G)
  remaining=set(range(n));count=0
  while remaining:
   count+=1;avail=set(remaining)
   while avail:
    j=min(avail);remaining.remove(j);avail.remove(j);avail.difference_update(np.flatnonzero(G[j]).tolist())
  return count

TYPES=[[('D',24)],[('D',16),('E',8)],[('E',8)]*3,[('A',24)],[('D',12)]*2,[('A',17),('E',7)],[('D',10),('E',7),('E',7)],[('A',15),('D',9)],[('D',8)]*3,[('A',12)]*2,[('A',11),('D',7),('E',6)],[('E',6)]*4,[('A',9)]*2+[('D',6)],[('D',6)]*4,[('A',8)]*3,[('A',7)]*2+[('D',5)]*2,[('A',6)]*4,[('D',4)]*6,[('A',5)]*4+[('D',4)],[('A',4)]*6]

def glue(L,scale):
 dims=[x.den for x in L];size=math.isqrt(math.prod(dims));assert size**2==math.prod(dims)
 def add(a,b):return tuple(x.add[i][j] for x,i,j in zip(L,a,b))
 zero=(0,)*len(L);good=set()
 for a in it.product(*(range(d) for d in dims)):
  N=sum(x.min[i]*(scale//x.den) for x,i in zip(L,a))
  if a==zero or (N>=4*scale and N%(2*scale)==0):good.add(a)
 gs=sorted(good)
 def visit(H):
  if len(H)==size:return H
  for g in gs:
   if g in H:continue
   x=g;cosets=set(H)
   while x not in H:
    C={add(x,h) for h in H}
    if not C<=good:break
    cosets|=C;x=add(x,g)
   else:
    if len(cosets)>len(H) and size%len(cosets)==0:
     z=visit(cosets)
     if z:return z
  return None
 H=visit({zero});assert H is not None;assert len(H)==size and all(add(a,b) in H for a in H for b in H)
 return sorted(H)

def profiles(L,C,scale):
 out=[]
 for a in C:
  opt=[[(N*(scale//x.den),i) for N,i in x.byclass[c]] for x,c in zip(L,a)]
  rem=[sum(min(N for N,i in X) for X in opt[j:]) for j in range(len(L)+1)]
  def go(j,N,P):
   if j==len(L):
    if N==6*scale:out.append(P)
    return
   for n,i in opt[j]:
    if N+n+rem[j+1]<=6*scale:go(j+1,N+n,P+(i,))
  go(0,0,())
 return out

def section(L,C,scale,w,mark):
 opt=[x.options(i,mark[1] if j==mark[0] else -1) for j,(x,i) in enumerate(zip(L,w))];ans=0;points=0
 for a in C:
  O=[d.get(c,[]) for d,c in zip(opt,a)]
  if not all(O):continue
  fac=[scale//x.den for x in L];rem=[sum(min(o[0] for o in O[k])*fac[k] for k in range(j,len(L))) for j in range(len(L)+1)]
  def go(j,N,p,S,ratio):
   nonlocal ans,points
   if j==len(L):
    if N==4*scale and p==3*scale:points+=S;ans+=S*ratio[0]//ratio[1]
    return
   f=fac[j]
   for n,d,s,b in O[j]:
    if N+n*f+rem[j+1]>4*scale:continue
    r=(b,s) if b*ratio[1]<ratio[0]*s else ratio;go(j+1,N+n*f,p+d*f,S*s,r)
  go(0,0,0,1,(1,1))
 return points,ans

def main(ids):
 cache={}
 for idx in ids:
  ty=TYPES[idx];print('START',idx,ty,flush=True);L=[]
  for t in ty:
   if t not in cache:cache[t]=Local(t)
   L.append(cache[t])
  scale=math.lcm(*(x.den for x in L));C=glue(L,scale);W=profiles(L,C,scale);print('GLUE',len(C),'PROFILES',len(W),flush=True)
  G=[];CS=set(C)
  for i,j in it.combinations(range(len(L)),2):
   if ty[i]==ty[j]:
    p=list(range(len(L)));p[i],p[j]=p[j],p[i]
    if all(tuple(a[k] for k in p) in CS for a in C):G.append(p)
  S=set(W);all_profiles=set(W);assert len(W)==len(S);reps=[]
  while S:
   w=min(S);S.remove(w);q=[w];reps.append(w)
   for v in q:
    for p in G:
     z=tuple(v[k] for k in p)
     assert z in all_profiles
     if z in S:S.remove(z);q.append(z)
  mx=0;pass_max=0;count=0;bad=[];print('PROFILE REPRESENTATIVES',len(reps),flush=True)
  for w in reps:
   pts,b=section(L,C,scale,w,(-1,-1))
   if b<=144:mx=max(mx,b);pass_max=max(pass_max,b);continue
   for j,x in enumerate(L):
    for m in x.markings(w[j]):
     pts,b=section(L,C,scale,w,(j,m));count+=1;mx=max(mx,b)
     if b>144:bad.append((w,j,m,pts,b))
     else:pass_max=max(pass_max,b)
  print('RESULT',idx,'max bound',mx,'explicit markings',count,'remaining',len(bad),'passed bound',pass_max,flush=True)
  print('REMAINING RECORDS',repr(bad),flush=True)
if __name__=='__main__':main(list(map(int,sys.argv[1:])) or list(range(20)))
