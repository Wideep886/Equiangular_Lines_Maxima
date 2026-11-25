#!/usr/bin/env python3
"""Second pass for unresolved marked profiles from niemeier20_independent.py."""
import ast,sys,itertools as it,collections as co,math
import numpy as np
import sympy as sp
from niemeier20_independent import Local,TYPES,glue,omega

def root_basis(x,i,mark):
 W=np.array(x.dom[i][2]);roots=x.positive[x.positive@W==0]
 if mark>=0:roots=roots[(roots@x.C)[:,mark]==0]
 rs={tuple(a) for a in roots};R=[]
 for a in sorted(rs,key=sum):
  if not any(tuple(u-v for u,v in zip(a,b)) in rs for b in rs if sum(b)<sum(a)):R.append(a)
 return np.array(R,dtype=np.int64).reshape(-1,x.n)

def buckets(x,i,mark,scale):
 W=np.array(x.dom[i][2]);WB=W@x.K; wnorm=x.dom[i][0]; out=co.defaultdict(lambda:co.defaultdict(list)); f=scale//x.den
 for N,c,U,B in x.x:
  dot=U@WB;mask=(3*x.den-dot)**2<=(4*x.den-N)*(6*x.den-wnorm)
  if mark>=0:mask&=U[:,mark]==0
  for d in np.unique(dot[mask]):
   out[c][(int(N)*f,int(d)*f)].append(U[mask&(dot==d)])
 return {c:{nd:np.concatenate(A) for nd,A in o.items()} for c,o in out.items()}

def shell(L,C,scale,w,mark,limit=10000):
 tables=[buckets(x,i,mark[1] if j==mark[0] else -1,scale) for j,(x,i) in enumerate(zip(L,w))];Xs=[];count=0
 for a in C:
  O=[d.get(c,{}) for d,c in zip(tables,a)]
  if not all(O):continue
  rem=[sum(min(n for n,d in o) for o in O[j:]) for j in range(len(O)+1)]
  def go(j,N,p,P):
   nonlocal count
   if j==len(O):
    if N==4*scale and p==3*scale:
     n=math.prod(len(x) for x in P);count+=n
     if count>limit:raise OverflowError(count)
     Xs.append(np.array([np.concatenate(v) for v in it.product(*P)],dtype=np.int64))
    return
   for (n,d),A in O[j].items():
    if N+n+rem[j+1]<=4*scale:go(j+1,N+n,p+d,P+[A])
  go(0,0,0,[])
 X=np.concatenate(Xs) if Xs else np.empty((0,24),dtype=np.int64);assert len(X)==count
 return X

def coloring(A):
 n=len(A)
 if not n:return 0
 if n<=80:return omega(A)
 # Each color is an independent set. Two deterministic orders, take min.
 orderings=[np.arange(n),np.argsort(-np.sum(A,axis=1),kind='stable')]
 best=n
 for order in orderings:
  remaining=np.ones(n,dtype=bool);k=0
  while remaining.any():
   k+=1;available=remaining.copy()
   for i in order:
    if available[i]:remaining[i]=False;available[i]=False;available&=~A[i]
  best=min(best,k)
 return best

def test(L,C,scale,w,mark,expected):
 X=shell(L,C,scale,w,mark)
 assert len(X)==expected,(len(X),expected)
 W=np.concatenate([x.dom[i][2] for x,i in zip(L,w)]);XS={tuple(x) for x in X};assert len(XS)==len(X) and all(tuple(W-x) in XS for x in X);Y=2*X-W;idx=np.argmax(Y!=0,axis=1);Y=Y[Y[np.arange(len(Y)),idx]>0];assert 2*len(Y)==len(X)
 p=0;B=[];Cs=[];Ks=[]
 for j,(x,i) in enumerate(zip(L,w)):
  U=Y[:,p:p+x.n];p+=x.n;D=U@x.K*(scale//x.den);B.append(D)
  R=root_basis(x,i,mark[1] if j==mark[0] else -1)
  if len(R):
   M=sp.Matrix(R@x.C@R.T);det=int(M.det());assert 0<det<=2**24;inv=np.array(M.inv()*det).astype(np.int64);assert np.max(np.abs(inv))<=2**23
   coords=U@R.T@inv;Cs.append(coords);Ks.append(det*D-scale*coords@R)
  else:Ks.append(D)
 B=np.concatenate(B,axis=1);K=np.concatenate(Ks,axis=1);Cpos=np.concatenate(Cs,axis=1) if Cs else np.empty((len(Y),0),dtype=np.int64)
 assert np.all(np.sum(Y*B,axis=1)==10*scale)
 A0=np.abs(Y@B.T)==2*scale;A=A0.copy()
 for j in range(len(Y)):
  a=np.all(K==K[j],axis=1);b=np.all(K==-K[j],axis=1)
  A[j]&=~(a&(np.all(Cpos>=Cpos[j],axis=1)|np.all(Cpos<=Cpos[j],axis=1)))
  A[j]&=~(b&(np.all(Cpos>=-Cpos[j],axis=1)|np.all(Cpos<=-Cpos[j],axis=1)))
 assert np.array_equal(A,A.T)
 unseen=set(range(len(Y)));groups=[]
 while unseen:
  i=min(unseen);unseen.remove(i);q=[i]
  for j in q:
   nxt=unseen.intersection(np.flatnonzero(~A0[j]).tolist());q.extend(sorted(nxt));unseen.difference_update(nxt)
  groups.append(np.array(q))
 def bound(S):return sum(coloring(A[np.ix_(g[S[g]],g[S[g]])]) for g in groups)
 S=np.ones(len(Y),dtype=bool);base=bound(S);deleted=0
 if base>72:
  # use all simple-root pairing coordinates and other root directions as functionals
  fs=[Y[:,i] for i in range(24)];p=0
  for x in L:
   fs.extend((Y[:,p:p+x.n]@x.positive.T).T);p+=x.n
  change=True
  while change and bound(S)>72:
   change=False
   for f in fs:
    T=S&(f!=0)
    if T.any() and bound(T)<=13:S[T]=False;change=True;deleted+=1
  # If still unresolved, exact max clique in all cells up to size 150.
 final=bound(S)
 return base,final,deleted,max(map(len,groups),default=0)

def main(paths):
 data={}
 for path in paths:
  idx=None
  for line in open(path):
   if line.startswith('START '):idx=int(line.split()[1])
   elif line.startswith('REMAINING RECORDS '):data[idx]=ast.literal_eval(line[len('REMAINING RECORDS '):])
 cache={}
 for idx,records in sorted(data.items()):
  if not records:continue
  L=[]
  for t in TYPES[idx]:
   if t not in cache:cache[t]=Local(t)
   L.append(cache[t])
  scale=math.lcm(*(x.den for x in L));C=glue(L,scale);bad=[];mx=0;mbase=0;ndel=0
  print('START',idx,len(records),flush=True)
  for k,(w,j,m,n,b) in enumerate(records):
   try:a,z,d,s=test(L,C,scale,w,(j,m),n)
   except OverflowError:a=z=100000;d=s=0
   mx=max(mx,z);mbase=max(mbase,a);ndel+=d
   if z>72:bad.append((w,j,m,n,b,a,z,s));print('UNRESOLVED',idx,k,bad[-1],flush=True)
   if k%50==0:print('PROGRESS',idx,k,a,z,d,s,flush=True)
  print('RESULT',idx,'records',len(records),'max axes bound',mx,'remaining',len(bad),'max base',mbase,'deletions',ndel,flush=True)
  print('REMAINING GLOBAL',repr(bad),flush=True)
if __name__=='__main__':main(sys.argv[1:])
