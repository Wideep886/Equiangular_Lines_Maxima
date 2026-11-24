"""Check which of the 15 exceptional centered shells are isometric, including their perpendicular roots."""
from ast import literal_eval
from collections import Counter
from itertools import combinations, product
from pathlib import Path
import json
import sys
import numpy as np

BASE = Path(__file__).resolve().parent
sys.path.insert(0, str(BASE / "revision2"))
from exceptional_shells import cases, codes
from odd_sections import signed_shell

RECORDS = {6: [32], 5: [26,27], 4: [0,1,5,6,10,14,15,16], 3: [1,2,3,4]}
PAIRS = [('A5_32','A2_1'),('A4_26','A3_10'),('A4_27','A2_3'),('A3_16','A2_2')]

def chart(joint, n, mark, h):
    scale = n
    if n == 6:
        sys.path.insert(0, str(joint/'a5'))
        from a5 import groups
        shell = {sum(choices, ()) for cell in groups(h, mark) for choices in product(*cell)}
    else:
        shell = signed_shell(n, codes(n), mark, h)
    w = sum(h, ())
    points = sorted(tuple(2*a-b for a,b in zip(x,w)) for x in shell)
    roots = []
    offset = 0
    for i, block in enumerate(h):
        for j,k in combinations(range(len(block)),2):
            if i == mark and (j<3 or k<3):
                continue
            for sign in ((-1,1) if n==6 and i==4 else (-1,)):
                if block[j] + sign * block[k] != 0:
                    continue
                v = [0]*len(w)
                v[offset+j] = n
                v[offset+k] = n*sign
                roots.extend((v, [-x for x in v]))
        offset += len(block)
    return scale, np.array(points, dtype=np.int64), np.array(roots, dtype=np.int64)

def gram(rows, cols, scale):
    product = rows @ cols.T
    assert (product % (scale*scale) == 0).all()
    return product // (scale*scale)

def independent(points):
    pivots = {}
    chosen = []
    for i,row in enumerate(points):
        row = row.copy() % 101
        for j in sorted(pivots):
            row = (row - row[j]*pivots[j]) % 101
        nz = np.flatnonzero(row)
        if len(nz):
            j = int(nz[0])
            pivots[j] = (row * pow(int(row[j]),-1,101)) % 101
            chosen.append(i)
    return chosen

def fingerprint(scale,points):
    G=gram(points,points,scale)
    assert np.all(np.diag(G)==10) and len(set(map(tuple,points)))==len(points)
    # tr(G^2), tr(G^3) and tr(G^4) via G^2. All entries fit in signed int64:
    # |G_ij|<=10, n<=576, so sum of squared entries of G^2 <= 10000 n^4.
    n=len(G)
    assert 10000*n**4 < 2**63
    G2=G@G
    return n,int(np.trace(G2)),int((G2*G).sum()),int((G2*G2).sum())

def main(checkpoint, export=None):
    tables={}
    for n,indices in RECORDS.items():
        for mark,k,h,count in cases(checkpoint,n):
            k=int(k)
            if k not in indices or (n==6 and mark!=4):
                continue
            label=f'A{n-1}_{k}'
            tab=chart(checkpoint,n,mark,h)
            assert len(tab[1])==count and label not in tables
            tables[label]=tab
            if export is not None:
                np.savez(export/label,Y=tab[1],s=n)
    assert len(tables)==15
    if export is not None:
        return
    signatures=Counter(fingerprint(n,p) for n,p,r in tables.values())
    assert sorted(signatures.values())==[1]*7+[2]*4
    images=json.loads((BASE/'section_isometries.json').read_text())
    for a,b in PAIRS:
        m,p,roots=tables[a];n,s,other=tables[b]
        ids=independent(p);left=p[ids];right=s[images[a+'='+b]]
        assert len(roots)==len(other) and len(ids)==21 and len(right)==21
        assert (gram(left,left,m)==gram(right,right,n)).all()
        assert sorted(map(tuple,gram(p,left,m)))==sorted(map(tuple,gram(s,right,n)))
        assert sorted(map(tuple,gram(roots,left,m)))==sorted(map(tuple,gram(other,right,n)))
        print(f'{a} = {b}: {len(p)} signed points, {len(roots)} roots, dimension 21')
    print('15 marked records: precisely 11 isometry types of centered point sets; all four maps preserve the perpendicular roots.')
if __name__=='__main__':
    main(Path(sys.argv[1]).resolve(), Path(sys.argv[2]) if len(sys.argv)==3 else None)
