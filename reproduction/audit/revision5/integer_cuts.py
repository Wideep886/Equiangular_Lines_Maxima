"""Verify the five indecomposable tetrad cuts from a list of weight-four and weight-eight codewords."""
from collections import Counter
from itertools import combinations
from functools import lru_cache
from pathlib import Path
from math import comb
import json,sys

def subset(items,n):
    for t in combinations(items,n):yield sum(1<<x for x in t)

def permute(w,g):
    z=0
    while w:
        b=w&-w;w-=b;z|=1<<g[b.bit_length()-1]
    return z

def generators(s):
    m=int(s[1:])
    def swap(*pairs):
        v=list(range(m))
        for i,j in pairs:v[i],v[j]=v[j],v[i]
        return tuple(v)
    if s[0]=='d':
        # Swap adjacent pairs, or switch both of two pairs.
        return [t for i in range(0,m-2,2) for t in (swap((i,i+2),(i+1,i+3)),swap((i,i+1),(i+2,i+3)))]
    assert s=='e7'
    # Transvections on the seven nonzero three-bit columns of the simplex code.
    return [tuple((j^(((j>>a)&1)<<b))-1 for j in range(1,8)) for a,b in permutations2(3)]

def permutations2(m):
    return [(i,j) for i in range(m) for j in range(m) if i!=j]

@lru_cache(None)
def local(s):
    m=int(s[1:]); perms=generators(s)
    masks=[u for n in range(4) for u in subset(range(m),n)]
    maps=[{u:permute(u,g) for u in masks} for g in perms]
    lookup={};weights={};by={}
    for a in range(4):
        for b in range(4):
            by[a,b]=[]
            for u in subset(range(m),a):
                for v in subset([j for j in range(m) if not(u>>j&1)],b):
                    if (u,v) in lookup:continue
                    rep=(u,v);by[a,b].append(rep)
                    wave=[rep];lookup[rep]=rep;size=1
                    while wave:
                        x,y=wave.pop()
                        for f in maps:
                            nxt=(f[x],f[y])
                            if nxt not in lookup:lookup[nxt]=rep;wave.append(nxt);size+=1
                    weights[rep]=size
    return by,lookup,weights

def cut_types(name,outer,words):
    labels=name.split('_');dims=[int(s[1:]) for s in labels]
    starts=[sum(dims[:i]) for i in range(len(dims))]
    data=[local(s) for s in labels]
    generated=[]
    for s,m,start in zip(labels,dims,starts):
        for g in generators(s):
            generated.append(tuple(range(start))+tuple(x+start for x in g)+tuple(range(start+m,24)))
    for g in [*generated,*outer]:
        assert sorted(g)==list(range(24)) and {permute(w,g) for w in words}==words
    locweights={}
    def walk(i,aa,bb,mask,free,weight):
        if i==len(labels):
            if aa==0 and bb==0:locweights[mask,free]=weight
            return
        by,lookup,ws=data[i]
        for a in range(aa+1):
            for b in range(bb+1):
                for u,v in by[a,b]:
                    walk(i+1,aa-a,bb-b,mask | u<<starts[i],free | v<<starts[i],weight*ws[u,v])
    walk(0,3,3,0,0,1)
    # The global generators permute the blocks and normalize the local group.
    # Verify the latter directly on each generator before using local orbit representatives.
    outer_r=[]
    for g in outer:
        records=[]
        for s,m,start in zip(labels,dims,starts):
            tgt=next(k for k,b in enumerate(starts) if b<=g[start]<b+dims[k])
            h=tuple(x-starts[tgt] for x in g[start:start+m])
            assert sorted(h)==list(range(m)) and labels[tgt]==s
            # An exact test on all local cuts is enough: each local orbit has one image orbit.
            by,lookup,ws=data[tgt]
            masks={x for u,v in lookup for x in (u,v)}
            table={u:permute(u,h) for u in masks};orb={}
            for (u,v),rep in lookup.items():
                image=lookup[table[u],table[v]]
                if rep in orb:assert orb[rep]==image
                else:orb[rep]=image
            assert len(set(orb.values()))==len(ws)
            records.append((m,start,starts[tgt],orb))
        outer_r.append(records)
    def apply(uv,recs):
        u,v=uv;U=V=0
        for m,src,dst,lut in recs:
            mask=(1<<m)-1
            p,q=lut[(u>>src)&mask,(v>>src)&mask];U|=p<<dst;V|=q<<dst
        return U,V
    todo=set(locweights);orbs=[]
    for rep in sorted(locweights):
        if rep not in todo:continue
        todo.remove(rep);wave=[rep];weight=0
        while wave:
            now=wave.pop();weight+=locweights[now]
            for recs in outer_r:
                nxt=apply(now,recs)
                if nxt in todo:todo.remove(nxt);wave.append(nxt)
        orbs.append((rep,weight))
    assert not todo and sum(w for x,w in orbs)==comb(24,3)*comb(21,3)
    return orbs

def components(words):
    rest=set(words);parts=[]
    while rest:
        r=rest.pop();part=[r];wave=[r]
        while wave:
            s=wave.pop();found={t for t in rest if (s&t).bit_count()==3}
            rest-=found;part.extend(found);wave.extend(found)
        parts.append(part)
    return parts

def capacity(part):
    if len(part)==2:
        assert (part[0]&part[1]).bit_count()==3
        return 8
    if len(part)==4:
        common=part[0]
        for v in part:common &= v
        if common.bit_count()!=1:return 99
        parts=[v^common for v in part]
        pairs={a&b for a,b in combinations(parts,2)}-{0}
        if len(pairs)!=4 or any(x.bit_count()!=2 for x in pairs) or any(x&y for x,y in combinations(pairs,2)):return 99
        # The four edges on the four pair labels make a square.
        if not all(sum((p&w)==p for w in parts)==2 for p in pairs):return 99
        return 10
    if len(part)==12:
        used=0
        for w in part:used |= w
        deg={j:sum((w>>j)&1 for w in part) for j in range(24) if used>>j&1}
        M=sum(1<<j for j in deg if deg[j]==6)
        N=sum(1<<j for j in deg if deg[j]==4)
        if used!=M|N:return 99
        stars={w&M for w in part};pairs={w&N for w in part}
        if (len(stars),len(pairs))!=(4,3) or any(s.bit_count()!=3 for s in stars) or any(s.bit_count()!=2 for s in pairs):return 99
        if any(x&y for x,y in combinations(pairs,2)):return 99
        edges={a&b for a,b in combinations(stars,2)}
        if len(edges)!=6 or any(x.bit_count()!=1 for x in edges):return 99
        return 18
    return 99

def cutoff(supports):
    used=0
    for w in supports:
        if w.bit_count()==1:return 'singleton'
        assert w.bit_count()==5
        used |= w
    if len(supports)<12 or used.bit_count()<18:return 'capacity'
    ids=[i for i in range(24) if used>>i&1]
    if any(all(w&(1<<i|1<<j) for w in supports) for i,j in combinations(ids,2)):return 'transversal'
    ids.sort(key=lambda i:-sum((w>>i)&1 for w in supports))
    def violates(S):
        r=S.bit_count();vals=sorted((w&S).bit_count() for w in supports)
        # At most five lines on each support: the first eleven full bins and three more.
        low=5*sum(vals[:11])+3*vals[11];hi=5*sum(vals[-11:])+3*vals[-12]
        gap=max(0,18*low-290*r,290*r-18*hi)
        return gap*gap>1508*r*(18-r)
    mask=0
    for i in ids[:-1]:
        mask|=1<<i
        if violates(mask):return 'trace'
    parts=components(supports);bounds=sorted(map(capacity,parts))
    if sum(bounds)<=57:return tuple(bounds)
    raise AssertionError(('unresolved',supports))

def check_hand_models(words):
    octads=words['golay']
    assert len(octads)==759 and all(v.bit_count()==8 for v in octads)
    five_subsets={s for v in octads for s in subset([j for j in range(24) if v>>j&1],5)}
    assert len(five_subsets)==comb(24,5)==759*comb(8,5)
    for w in words['d24']:
        assert all((w>>(2*j)&3) in (0,3) for j in range(12))
    for name,block_sizes in [('3e8',(8,8,8)),('d16_e8',(16,8))]:
        offset=0; masks=[]
        for size in block_sizes:
            masks.append(((1<<size)-1)<<offset);offset+=size
        for w in words[name]:
            parts=[w&mask for mask in masks if w&mask]
            assert len(parts)==1 or (len(parts)==2 and all(v.bit_count()==4 and v in words[name] for v in parts))
    for j in [1,2]:
        for d in range(13-j,17-j):
            assert (18*58*(4-j)-290*d)**2>1508*d*(18-d)
    print('Four hand cases: Golay design, paired d24, and the two direct sums verified.')

def main():
    source=Path(sys.argv[1]); words={row[0]:set(map(int,row[2:])) for row in map(str.split,source.read_text().splitlines())}
    outer=json.loads(Path(__file__).with_name('outer.json').read_text())
    check_hand_models(words)
    assert set(words)==set(outer)|{'golay','3e8','d16_e8','d24'}
    total=0
    for name,lifts in outer.items():
        c=Counter();labelled=Counter();reps=cut_types(name,lifts,words[name]);total+=len(reps)
        for (one,zero),weight in reps:
            supports=[w^one for w in words[name] if w&(one|zero)==one]
            reason=cutoff(supports);c[reason]+=1;labelled[reason]+=weight
        print(name,len(reps),'{'+', '.join(f'{k}: {c[k]}' for k in sorted(c,key=str))+'}',
               'special cuts',[(str(k),labelled[k]) for k in sorted(labelled,key=str) if isinstance(k,tuple)],flush=True)
    print('All five codes:',total,'representatives of',5*comb(24,3)*comb(21,3),'labelled cuts.')
if __name__=='__main__':main()
