#!/usr/bin/env python3
"""Independent ADE type coverage check, not a test of certificate arithmetic."""
from collections import Counter
from pathlib import Path
import json

HERE = Path(__file__).resolve().parent
# Generate multisets independently by multiplicities for each irreducible type.
TYPES = [("A", i) for i in range(1, 20)] + [("D", i) for i in range(4, 20)] + [("E", i) for i in (6, 7, 8)]
def multisets(i, budget, current):
    if i == len(TYPES):
        if budget == 0:
            yield tuple(current)
        return
    t = TYPES[i]
    if budget == 0:
        yield tuple(current)
        return
    for k in range(budget // t[1] + 1):
        yield from multisets(i+1, budget-k*t[1], current+[t]*k)

def canonical(s):
    return tuple(sorted((a[0], int(a[1:])) for a in s.split('+')))

def read_types(filename, expected_count):
    ordinary=set(); branched={}; total=0
    with (HERE/filename).open() as f:
        for line in f:
            a=json.loads(line); t=canonical(a['type']); sc=a.get('square_class'); total+=1
            if sc is None:
                assert t not in ordinary and t not in branched
                ordinary.add(t)
            else:
                assert t not in ordinary and sc not in branched.setdefault(t,set())
                branched[t].add(sc)
    assert total==expected_count
    return ordinary|set(branched),branched

all_types=set().union(*(set(multisets(0,r,[])) for r in range(11,20)))
allowed={t for t in all_types if ('E',8) not in t and all(a!='D' or d<11 for a,d in t)}
mixed, branches=read_types('root19_all_mixed_certificates.jsonl',6571)
other, other_branches=read_types('root19_homogeneous_DE_recheck.jsonl',140)
special={tuple([('A',1)]*r) for r in (11,12,13)}|{tuple([('A',2)]*7)}
expected_mixed={t for t in allowed if len(set(t))>1 and any(a=='A' for a,d in t)}
assert len(all_types)==7328 and len(allowed)==6709
assert mixed==expected_mixed and len(mixed)==6565 and len(other)==140
assert not (mixed&other or mixed&special or other&special)
assert allowed==mixed|other|special
assert branches=={canonical('A1+A1+A3+A3+A3+A3+D4'):{1,2,3,5,6,10,14}}
assert not other_branches
print('rank  all_types  allowed  mixed  other  special')
for r in range(11,20):
    print(r,*(sum(sum(d for _,d in t)==r for t in S) for S in (all_types,allowed,mixed,other,special)))
print('PASS: 6709 = 6565 + 140 + 4; exact disjoint coverage.')
print('Metadata coverage only; run the other verifiers for the certificate arithmetic.')
