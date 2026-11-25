"""Check the six edge and four vertex blocks of A3^8 mark 1. Uses mark0.py."""
from mark0 import cells, rank

h = ((-2,-2,-2,6),(-2,-2,2,2),(-2,-2,2,2)) + ((0,0,0,0),)*4 + ((-2,-2,2,2),)
A=(1,2,7)
B=(3,4,5,6)

def main():
    big={};small={}
    for cell in cells(h):
        scaled=[tuple(tuple(2*z-w for z,w in zip(x,h[i]))for x in values)
                for i,values in enumerate(cell)]
        assert scaled[0]==((0,0,0,0),)
        fixed=tuple(sum(scaled[i][0][:2]) for i in A)
        assert not sum(fixed)
        dirs={}
        for i in A:
            if len(scaled[i])==2:
                a,b=scaled[i]
                assert all(x==-y for x,y in zip(a,b))
                coeff=tuple(j for j in range(4) if a[j])
                assert coeff in ((0,1),(2,3))
                dirs[i]=coeff[0]//2
            else:
                assert len(scaled[i])==1 and scaled[i][0][0]==scaled[i][0][1] and scaled[i][0][2]==scaled[i][0][3]
        facs=[i for i in B if len(scaled[i])>1]
        for i in B:
            vals=scaled[i]
            assert all(sum(col)==0 for col in zip(*vals))
            if len(vals)>1:
                v=vals[0]
                assert sorted(v.count(a) for a in set(v))==([1,3] if len(vals)==4 else [2,2])
        if any(fixed):
            if next(x for x in fixed if x)<0:continue
            assert len(dirs)==1 and len(facs)==2 and all(len(scaled[i])==4 for i in facs)
            assert tuple(facs) not in big
            big[tuple(facs)]=(next(iter(dirs.items())),fixed)
        else:
            assert len(dirs)==3 and len(facs)==1 and len(scaled[facs[0]])==6
            assert facs[0] not in small
            small[facs[0]]=set(dirs.items())
    assert len(big)==6 and len(small)==4 and len({a for a,c in big.values()})==6
    for i in B:
        assert {tag for edge,(tag,c) in big.items() if i not in edge}==small[i]
        assert rank([c for edge,(a,c) in big.items() if i not in edge])==2
    print('tetrahedron edge: binary factor, pair 0/1 and center coordinates')
    for edge,(a,c) in sorted(big.items()):
        print(''.join(map(str,edge)),a,c)
    print('vertex blocks: each uses the three opposite edges and its own octahedron')
    print('complementary triangle center ranks: 2 2 2 2')

if __name__=='__main__':
    main()
