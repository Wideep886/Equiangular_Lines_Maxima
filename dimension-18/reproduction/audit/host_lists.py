"""Independently enumerate ADE types and read raw Niemeier chain certificates (stdlib only)."""
from gzip import open as gzopen
from json import load
from pathlib import Path
import sys

def all_root_types():
    found = set()
    for h in range(2,47):
        blocks = [(h-1,f"A{h-1}")]
        if h%2 == 0 and h>=6:
            rank = 1+h//2
            blocks.append((rank,f"D{rank}"))
        if h in (12,18,30):
            rank = {12:6,18:7,30:8}[h]
            blocks.append((rank,f"E{rank}"))
        def add(total, start, taken):
            if total == 24:
                found.add(tuple(sorted(taken)))
            for i in range(start,len(blocks)):
                rank,name=blocks[i]
                if total+rank<=24:
                    add(total+rank,i,taken+[name])
        add(0,0,[])
    return found

def disjoint_paths(n,links):
    previous={}
    following={}
    for u,v in links:
        assert 0<=u<n and 0<=v<n and u!=v
        assert u not in following and v not in previous
        following[u]=v
        previous[v]=u
    roots=set(range(n))-set(previous)
    seen=set()
    for root in roots:
        v=root
        while v is not None:
            assert v not in seen
            seen.add(v)
            v=following.get(v)
    assert len(seen)==n
    return len(roots)

def main(joint):
    full=all_root_types()
    assert len(full)==23
    raw={}
    for name in ('thirteen','four-chains','marked-chains'):
        for path in sorted((joint/name).glob('chains_*.json.gz')):
            labels=tuple(sorted(path.name.removeprefix('chains_').removesuffix('.json.gz').split('_')))
            assert labels not in raw
            with gzopen(path,'rt') as f:
                contents=load(f)
            bounds=[]
            for record in contents:
                w,n,edges,*marked=record
                assert len(w)==24
                cap=disjoint_paths(n,edges)
                if marked and cap>=116:
                    refined=[disjoint_paths(m,subedges) for _,_,m,subedges in marked[0]]
                    cap=min(cap,max(refined,default=0))
                bounds.append(cap)
            raw[labels]=max(bounds)
    assert len(raw)==18 and max(raw.values())==112
    leftover=set(full)-set(raw)
    expected={('A1',)*24, ('A2',)*12, ('A3',)*8, ('A4',)*6, ('A5',)*4+('D4',)}
    assert set(raw).issubset(full) and leftover==expected
    for labels,cap in sorted(raw.items(),key=lambda pair: (pair[1],pair[0])):
        print(' '.join(labels),cap)
    print('Coxeter decompositions: 23. Raw chain types: 18. Largest oriented bound:',max(raw.values()))
    print('Complement: A1^24; A2^12, A3^8, A4^6, A5^4 D4. The Leech lattice is rootless.')
if __name__=='__main__':
    if not __debug__:
        raise RuntimeError("Enable assertions.")
    main(Path(sys.argv[1]))
