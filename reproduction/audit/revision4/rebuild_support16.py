"""Regenerate all 21 compact rank-16 domains directly from the incoming frontier.
Usage: python3 rebuild_support16.py ROOT_PROFILES_DIRECTORY
Pure Python 3.9+, deterministic. No clique solver, optimization or orbital point coordinates.
Only the provisionally retained anchor pairs and optional comparison data are read.
"""
from fractions import Fraction as F
from itertools import product
from math import isqrt, lcm
from pathlib import Path
import hashlib,json,sys
def local(T):
    n=int(T[1:])
    if T[0]=='A':
        return tuple(F(2*k*(n+1-k),n+1) for k in range(n+1)),tuple(min(k,n+1-k) for k in range(n+1))
    assert T=='D4'
    return (F(0),F(2),F(2),F(2)),(0,2,2,2)
def rebuild(name,front):
    types=name.split('.')
    data=[local(t) for t in types]
    S=lcm(*(w.denominator for table,_ in data for w in table))
    wt=[tuple(int(S*w) for w in table) for table,_ in data]
    shapes={};V=[]
    groups=[[i for i,t in enumerate(types) if t==T] for T in dict.fromkeys(types)]
    for v in product(*(range(len(table)) for table in wt)):
        e=tuple(weights[k] for weights,k in zip(wt,v))
        B=5*S-sum(e)
        if B<0 or (B and sum(pair[1][k] for pair,k in zip(data,v))>3):continue
        sig=tuple(q for g in groups for q in sorted(e[i] for i in g))
        shapes.setdefault(sig,(v,B));V.append((v,B,sig))
    index={s:k for k,s in enumerate(sorted(shapes))}
    anchors={index[s]:ar for s,ar in shapes.items()}
    def valid(n,par):return n%(4*S)==(par%2)*S and (n==0 or n>=5*S)
    def norm(a,v,m,n):
        return sum(tab[(((m%2)*i)^((n%2)*j)) if t=='D4' else (m*i+n*j)%len(tab)]
                   for t,tab,i,j in zip(types,wt,a,v))
    def neg(v):
        return tuple(k if t=='D4' else -k%len(tab) for t,tab,k in zip(types,wt,v))
    out={}
    for anchor,p in sorted(front):
        a,A=anchors[anchor]
        assert A>0;W=set()
        for v,B,s in V:
            if B and (index[s],p) not in front:continue
            if not all(valid(norm(v,v,0,n)+n*n*B,n) for n in (2,3)):continue
            # Norms for all m*a+n*v listed below must have the required congruence and minimal size.
            pairs=((1,1),(-1,1),(2,1),(-2,1),(1,2),(-1,2),(3,1),(-3,1))
            checks=[(m,n,norm(a,v,m,n)+m*m*A+n*n*B) for m,n in pairs]
            lim=isqrt(4*A*B)
            for y in range(-lim,lim+1):
                K=4*A*B-y*y
                if K%p or isqrt(K//p)**2!=K//p:continue
                if not all(valid(base+m*n*y,m+n) for m,n,base in checks):continue
                b=isqrt(K//p)
                if b==0 and y<0:w=neg(v)+(-y,0)
                elif b==y==0:w=min(v,neg(v))+(0,0)
                else:w=v+(y,b)
                W.add(w)
        out[name+':'+str(anchor)+':'+str(p)]=sorted(W)
    return out
def main():
    root=Path(sys.argv[1])
    src=root/'frontier_current.json'
    print('frontier SHA256',hashlib.sha256(src.read_bytes()).hexdigest())
    front=json.loads(src.read_text())['rank16']
    result={}
    for name,f in sorted(front.items()):
        if f:result.update(rebuild(name,set(map(tuple,f))))
    total=sum(map(len,result.values()));assert len(result)==21 and total==2307
    canonical=json.dumps(result,sort_keys=True,separators=(',',':')).encode()
    print('terminal domains',len(result),'word occurrences',total,'label-set SHA256',hashlib.sha256(canonical).hexdigest())
    old=root/'orbit16.json'
    if old.exists():
        exp=json.loads(old.read_text())
        assert set(exp)==set(result)
        assert all([list(w) for w in W]==exp[key][2] for key,W in result.items())
        print('Every labelled word agrees with the archived 21 domains.')
if __name__=='__main__':main()
