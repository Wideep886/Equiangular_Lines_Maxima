"""Calculate parabolic midrank widths using the imbalance of root signs.
Run with joint23 and the revision-2 audit folder containing exceptional_shells.py.
"""
from collections import Counter
from functools import cache
from itertools import combinations
from pathlib import Path
import sys
sys.path.insert(0,sys.argv[2])
from exceptional_shells import blocks,codes,cases

@cache
def loc(n,tag,h,mark):
    # Sign convention makes e_i-e_j, e_i+e_j (i<j) positive.
    d=len(h)
    roots=[]
    for i,j in combinations(range(d),2):
        for s in ((-1,1) if (n==6 and d==4) else (-1,)):
            if h[i]+s*h[j]:continue
            if mark and (i<3 or j<3):continue # no root has equal coefficients in all 3 positions otherwise
            roots.append((i,j,s))
    hist=Counter()
    for E,v in blocks(n,d,tag,mark):
        if E>4*n*n:continue
        ip=sum(a*b for a,b in zip(h,v))
        values=(v[i]+s*v[j] for i,j,s in roots)
        score=sum((k>0)-(k<0) for k in values)
        hist[E,ip,score]+=1
    return tuple((e,i,s,m) for (e,i,s),m in hist.items())

def widths(n,w,mark,words):
    H=4*n*n; dot=3*n*n; hist=Counter()
    for word in words:
        parts=[loc(n,g,h,i==mark) for i,(g,h) in enumerate(zip(word,w))]
        m=[min((row[0] for row in p),default=H+1) for p in parts]
        if sum(m)>H:continue
        tails=[sum(m[i+1:]) for i in range(len(m))]
        freq={(0,0,0):1}
        for p,l in zip(parts,tails):
            out=Counter()
            for (E,B,K),a in freq.items():
                for e,b,k,c in p:
                    if e+E+l<=H:out[E+e,B+b,K+k]+=a*c
            freq=out
        hist.update({s:m for (e,b,s),m in freq.items() if (e,b)==(H,dot)})
    return hist

def main():
    records=[]
    for n in (6,5,4,3):
        words=codes(n); alive=[];table=[]
        for mark,k,w,total in cases(Path(sys.argv[1]),n):
            hist=widths(n,w,mark,words);assert sum(hist.values())==total
            b=hist[0]+hist[1];table.append(b)
            records.append((n,-1 if mark is None else mark,int(k),total,b))
            if b>=116:alive.append((mark,int(k),b))
        print(n,len(table),max(table),len(alive),alive,flush=True)
    if len(sys.argv)>3:
        Path(sys.argv[3]).write_text(''.join(' '.join(map(str,t))+'\n' for t in sorted(records)))

if __name__=='__main__':main()
