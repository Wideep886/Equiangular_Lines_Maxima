"""Build explicit positive-root chain covers for the four marked Niemeier models.

Usage: python3 generate_chains.py EL18_58_search/joint23 audit/revision2 output.json.gz
The generator uses the revision-2 local cosets and marking representatives.
The independent verifier does not require the matching algorithm.
"""
import argparse
from collections import defaultdict
from functools import cache
import gzip
from itertools import combinations
import json
from pathlib import Path
import sys


@cache
def local_candidates(scale, coset, wpart, marked):
    grouped = defaultdict(list)
    for squared, x in blocks(scale, len(wpart), coset, marked):
        if squared <= 4 * scale * scale:
            pairing = sum(t*h for t,h in zip(x, wpart))
            grouped[squared, pairing].append(x)
    return dict(grouped)


def shell(scale, mark, w, code):
    energy = 4 * scale * scale
    pairing = 3 * scale * scale
    points = []
    for c in code:
        groups = [local_candidates(scale, v, w[i], i == mark) for i, v in enumerate(c)]
        if any(not g for g in groups):
            continue
        n = len(c)
        possible = [{(0,0)} for _ in range(n+1)]
        for i in range(n-1, -1, -1):
            possible[i] = {(e+a,p+b) for e,p in groups[i] for a,b in possible[i+1]
                           if e+a<=energy}
        if (energy, pairing) not in possible[0]:
            continue
        def walk(i, e, p, prefix):
            if i == n:
                points.append(prefix)
                return
            for (a,b), vectors in groups[i].items():
                if (e-a, p-b) in possible[i+1]:
                    for v in vectors:
                        walk(i+1, e-a, p-b, prefix+v)
        walk(0,energy,pairing,())
    assert len(set(points)) == len(points)
    return sorted(points)


def positive_roots(scale, mark, w):
    roots = []
    start = 0
    for factor, h in enumerate(w):
        signs = (-1,1) if scale == 6 and len(h) == 4 else (-1,)
        for i,j in combinations(range(len(h)), 2):
            if factor == mark and i<3:
                continue
            for s in signs:
                if h[i]+s*h[j]==0:
                    roots.append((start+i,start+j,s))
        start+=len(h)
    return roots


def reflection_links(points, roots):
    idx = {x:i for i,x in enumerate(points)}
    neighbors = []
    for x in points:
        adj = []
        for a,b,s in roots:
            t = x[a]+s*x[b]
            if t < 0:
                y = list(x)
                y[a] -= t
                y[b] -= s*t
                adj.append(idx[tuple(y)])
        neighbors.append(sorted(adj))
    return neighbors


def match(neighbors):
    """Find vertex-disjoint incoming and outgoing links by augmenting paths."""
    left = {}
    right = {}
    for source in range(len(neighbors)):
        parent = {}
        seen = {source}
        queue = [source]
        end = None
        for i in queue:
            for j in neighbors[i]:
                if j in parent:
                    continue
                parent[j] = i
                if j not in right:
                    end = j
                    break
                k = right[j]
                if k not in seen:
                    seen.add(k)
                    queue.append(k)
            if end is not None:
                break
        while end is not None:
            i = parent[end]
            old_end = left.get(i)
            left[i],right[end] = end,i
            end = old_end
    return sorted(left.items())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('joint',type=Path)
    parser.add_argument('revision2',type=Path)
    parser.add_argument('output',type=Path)
    args=parser.parse_args()
    if not __debug__:
        raise SystemExit('Run with Python assertions enabled.')
    sys.path.insert(0,str(args.revision2.resolve()))
    global blocks
    from exceptional_shells import blocks, cases, codes
    records = []
    for scale in (3,4,5,6):
        chunk = []
        for mark, k, w, n in cases(args.joint,scale):
            pts = shell(scale,mark,w,codes(scale))
            assert len(pts)==n
            links = match(reflection_links(pts,positive_roots(scale,mark,w)))
            chunk.append(dict(scale=scale,mark=-1 if mark is None else int(mark),index=int(k),
                              points=n,links=links))
        cap = [r['points']-len(r['links']) for r in chunk]
        print(scale,len(chunk),sum(r['points'] for r in chunk),sum(len(r['links']) for r in chunk),
              max(cap),sum(v>=116 for v in cap),flush=True)
        records.extend(chunk)
    records.sort(key=lambda r:(r['scale'],r['mark'],r['index']))
    with args.output.open('wb') as file:
        with gzip.GzipFile(filename='',mode='wb',fileobj=file,mtime=0) as output:
            output.write(json.dumps(records,separators=(',',':')).encode('ascii'))

if __name__=='__main__':
    main()
