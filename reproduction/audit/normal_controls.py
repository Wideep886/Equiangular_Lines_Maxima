"""Deterministic satisfiable and unsatisfiable controls for the archived normal sieve."""
from itertools import combinations
from hashlib import sha256
from pathlib import Path
from subprocess import check_output
import sys

def make_rows(n):
    # 18 independent words followed by two-support columns; they span exactly a fixed 18-space.
    basis = [1 << i for i in range(18)]
    words = basis + [a ^ b for a,b in combinations(basis, 2)]
    return words[:n]

def make_input(d, words, conflict):
    # Add one incompatible point in each extra direction. The whole input spans exactly d.
    nbase=len(words)
    words=words+[1 << i for i in range(18,d)]
    n=len(words)
    all_bits=(1 << nbase)-1
    graph=[all_bits ^ (1 << i) for i in range(nbase)] + [0]*(n-nbase)
    if conflict:
        graph[0] &= ~(1 << 1)
        graph[1] &= ~1
    return ('\n'.join([f'{d} 1',str(n),' '.join(map(str,words)),' '.join(map(str,graph))])+'\n').encode()

def multistage_fixture(d,missing_edge):
    # Every nonzero equation on the 18-space deletes at least seven of the core words.
    nbase=58 if missing_edge else 57
    words=[]
    for i in range(1000):
        word=int.from_bytes(sha256(f"EL18-control-{i}".encode()).digest()[:3],"big") % (1<<18)
        if word and word not in words:
            words.append(word)
        if len(words)==nbase:
            break
    assert len(words)==nbase
    signatures=[sum(1<<i for i,w in enumerate(words) if w>>j & 1) for j in range(18)]
    minimum=nbase
    current=0
    for n in range(1,1<<18):
        current^=signatures[(n & -n).bit_length()-1]
        minimum=min(minimum,current.bit_count())
    assert minimum>=7
    quotient={20:[1,1,2], 21:[1,1,1,2,2,4,4]}[d]
    if missing_edge:
        quotient={20:[1,2,3],21:[1,2,3,4,5,6]}[d]
    all_words=words+[(v<<18)|i for i,v in enumerate(quotient)]
    n=len(all_words)
    full=(1<<n)-1
    graph=[full^(1<<i) for i in range(n)]
    if missing_edge:
        graph[0] &= ~2
        graph[1] &= ~1
    data=('\n'.join([f'{d} 1',str(n),' '.join(map(str,all_words)),' '.join(map(str,graph))])+'\n').encode()
    return minimum,data

def main(binary):
    for d in range(18, 22):
        words=make_rows(58)
        for conflict in (False,True):
            out=check_output([binary],input=make_input(d,words,conflict)).decode().splitlines()
            if conflict:
                assert out[-1]=='first-stage result 0',out
            elif d in (18,19):
                assert out[-1]=='first-stage result 1',out
            else:
                parts=out[-1].split()
                assert parts[:3]==['relaxed','kernel','witness'],out
                masks=list(map(int,parts[3:-2]))
                assert len(masks)==d-18 and parts[-2:]==['capacity','58']
                assert len({0}|{__import__('functools').reduce(int.__xor__,combo) for k in range(1,len(masks)+1) for combo in combinations(masks,k)})==1<<len(masks)
                assert all(not((w&a).bit_count()%2) for w in words for a in masks)
            print('ambient binary rank parameter',d,'cell capacity',58-int(conflict),out[-1],flush=True)
    for d in (20,21):
        for missing_edge in (False,True):
            minimum,data=multistage_fixture(d,missing_edge)
            out=check_output([binary],input=data).decode().splitlines()
            assert out[-1].startswith('no eligible kernel; anchors '),out
            fields=out[-1].split()
            assert int(fields[4])>=1,out
            if missing_edge:
                assert int(fields[7])>=1,out
                if d==21:
                    assert int(fields[-1])>=1,out
            print('Multistage negative control, dimension',d,'minimum core weight',minimum,
                  'core missing edge',int(missing_edge),'global capacity 57:',out[-1])
if __name__ == '__main__':
    if not __debug__:
        raise RuntimeError("Enable assertions.")
    main(str(Path(sys.argv[1]).resolve()))
