"""Per-record census of the eighteen hosts.
Usage: poset_records.py JOINT23 {direct,orbit,historical} OUT.
Only historical mode reads chains. It provides a diagnostic comparison
with the separate proof, not input to either midrank method.
"""
from pathlib import Path
from math import lcm
import json,gzip,sys
from poset_width import factors,glue,shapes,HOSTS
import poset_width,poset_direct
joint=Path(sys.argv[1]);method=sys.argv[2];out=Path(sys.argv[3]);assert method in ('direct','orbit','historical')
records=[]
for host in HOSTS:
    fs=factors(host.split('_'));code=glue(fs);L=lcm(*(f.d for f in fs))
    if method=='historical':
        path,=joint.glob('*/chains_'+host+'.json.gz')
        with gzip.open(path,'rt') as fp:
            for hw,n,links,*rest in json.load(fp):
                hw=tuple(hw);records.append((host,hw,-1,-1,n-len(links)))
                if rest:
                    records.extend((host,hw,i,j,m-len(edges)) for i,j,m,edges in rest[0])
        continue
    if method=='orbit':sh=shapes(fs,code,4)
    def width(w,mark):
        return poset_direct.bound(fs,code,w,mark,L) if method=='direct' else poset_width.bound(fs,w,sh,L,mark)
    blank=tuple(() for f in fs)
    for w in shapes(fs,code,6):
        h=sum(w,());b=width(w,blank);records.append((host,h,-1,-1,b))
        if b<116:continue
        off=0
        for nb,(f,u) in enumerate(zip(fs,w)):
            for i in range(f.n):
                for j in range(i):
                    if u[i] or u[j] or f.C[i][j]!=-1:continue
                    mark=list(blank);mark[nb]=(j,i)
                    records.append((host,h,off+i,off+j,width(w,tuple(mark))))
            off+=f.n
out.write_text(''.join(f'{host} {"".join(map(str,w))} {i} {j} {v}\n' for host,w,i,j,v in sorted(records)))
