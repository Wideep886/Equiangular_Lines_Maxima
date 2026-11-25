"""Finite margin checks for the elementary discriminant lemma. Python 3 stdlib."""
from itertools import combinations, product

def even_prime(rank):
    pairs = list(combinations(range(rank), 2))
    masks = list(range(1, 1 << rank))
    with_bits = [[i for i in range(rank) if (x >> i) & 1] for x in masks]
    pair_bits = [[j for j,(i,k) in enumerate(pairs) if i in bits and k in bits]
                 for bits in with_bits]
    total = anisotropic = all_even = rank_three_cyclic = 0
    for cross in product((0,2), repeat=len(pairs)):
        cross_values = [sum(cross[j] for j in active) for active in pair_bits]
        for diag in product(range(4), repeat=rank):
            total += 1
            Q = [(s + sum(diag[j] for j in active)) % 4
                 for s,active in zip(cross_values,with_bits)]
            if 0 in Q: continue
            anisotropic += 1
            if all(v % 2 == 0 for v in Q): all_even += 1
            if rank != 3: continue
            v = next(i for i,x in zip(masks,Q) if x % 2)
            Q = [0]+Q
            comp = {x for x in range(1 << rank) if (Q[x^v]-Q[x]-Q[v]) % 4 == 0}
            assert v not in comp and len(comp)==4 and len(comp|{x^v for x in comp})==8
            rank_three_cyclic += 1
    print(f'rank {rank}, tested {total}, anisotropic {anisotropic}, all even {all_even}, rank-three cyclic split {rank_three_cyclic}')
    if rank == 4: assert anisotropic == 0
    if rank == 3: assert rank_three_cyclic == anisotropic > 0 and all_even == 0

def odd_prime(p):
    # Scan every ternary quadratic form, including degenerate ones.
    line_vectors=[]
    for i in range(3):
        line_vectors.extend([(0,)*i+(1,)+x for x in product(range(p),repeat=2-i)])
    basis=[tuple(v[i]*v[j] % p for i in range(3) for j in range(i,3)) for v in line_vectors]
    count=0
    for coeffs in product(range(p),repeat=6):
        assert any(sum(a*b for a,b in zip(coeffs,b)) % p==0 for b in basis)
        count+=1
    print(f'odd prime {p}, tested {count} ternary forms, anisotropic 0')

for dimension in range(1,5): even_prime(dimension)
for prime in (3,5): odd_prime(prime)
