"""Print terminal component-support multisets from the public archived word domains.
Reads orbit16.json. No search, class weights, witnesses or capacities are used.
A component name followed by a superscript is encoded by repeating the label.
"""
from collections import Counter
from pathlib import Path
import hashlib, json, sys
P = Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).with_name('orbit16.json')
data = json.loads(P.read_bytes())
print("input SHA256", hashlib.sha256(P.read_bytes()).hexdigest())
for k,(_,weights,words,_) in sorted(data.items()):
    types = k.split(':')[0].split('.')
    shapes = Counter()
    for word in words:
        assert len(word) == 2+len(types)
        # Triality makes every nonzero D4 class an orthogonal-axis orbit.
        shape = tuple(sorted(t+"["+str(min(j,int(t[1:])+1-j) if t.startswith('A') else 1)+"]"
                             for t,j in zip(types,word) if j))
        shapes[shape] += 1
    print(k, len(words), ' | '.join(' '.join(s)+':'+str(n) for s,n in sorted(shapes.items())))
from collections import defaultdict
from fractions import Fraction as F
print("Residual position counts and coordinates, two final root systems")
for key,(_,_,words,_) in sorted(data.items()):
    tags = key.split(':')[0].split('.')
    if (tags.count('A1'), tags.count('A3'), tags.count('A5'), tags.count('D4')) not in [(5,2,1,0),(6,2,0,1)]:
        continue
    coords=defaultdict(set)
    for w in words:
        proj=0
        for t,k in zip(tags,w):
            n=int(t[1:]); proj += F(2*k*(n+1-k),n+1) if t.startswith('A') else (2 if k else 0)
        coords[5-proj].add(tuple(w[-2:]))
    print(key, {str(h):sorted(z) for h,z in sorted(coords.items())})
