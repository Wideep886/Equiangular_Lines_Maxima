"""Check a saved Gram matrix using only integer Newton identities, independent of coordinates."""
import json
import sys
from pathlib import Path

def main(path):
    gram = json.loads(path.read_text())["gram"]
    n = len(gram)
    if n != 57 or any(len(row) != n for row in gram):
        raise ValueError("Gram shape differs from 57 by 57")
    for i, row in enumerate(gram):
        for j, value in enumerate(row):
            if type(value) is not int or value != gram[j][i]:
                raise ValueError("Gram is not a symmetric integer matrix")
            if value not in ((5,) if i == j else (-1, 1)):
                raise ValueError("Gram entry is invalid")
    # G = J - 2A + 4I, with A the graph of strictly negative entries. Use sparse sums.
    neg = [[j for j, val in enumerate(row) if val == -1] for row in gram]
    power = [[int(i == j) for j in range(n)] for i in range(n)]
    traces = []
    es = [1]
    for k in range(1, n+1):
        new = []
        for row in power:
            total = sum(row)
            new.append([total + 4*row[j] - 2*sum(row[l] for l in neg[j]) for j in range(n)])
        power = new
        traces.append(sum(power[i][i] for i in range(n)))
        num = sum((1 if i%2 else -1)*es[k-i]*traces[i-1] for i in range(1,k+1))
        e, rem = divmod(num,k)
        if rem:
            raise ArithmeticError("Newton identity returned a noninteger coefficient")
        es.append(e)
    if min(es[:19]) <= 0 or any(es[19:]):
        raise ValueError("Expected inertia not established")
    # det(tI+G) has positive coefficients in nonzero degrees and no positive real root.
    # G is symmetric, so the 18 nonzero eigenvalues are all positive.
    print("Shape 57 by 57, symmetric integral, diagonal 5, off diagonal +/-1")
    print("det(tI+G) = t^39 times monic degree 18; descending coefficients:")
    print(" ".join(map(str,es[:19])))
    # Integer factorisation checked by multiplication, independent of numerical roots.
    factored = [1]
    for term in ([1,18],)*8 + ([1,14],)*8 + ([1,29,192],):
        result = [0]*(len(factored)+len(term)-1)
        for i, a in enumerate(factored):
            for j, b in enumerate(term):
                result[i+j] += a*b
        factored = result
    if es[:19] != factored:
        raise ArithmeticError("Factorization mismatch")
    print("det(tI+G) = t^39 (t+14)^8 (t+18)^8 (t^2+29t+192)")
    print("Newton identities: strictly positive degree-18 factor coefficients, 39 trailing zeros.")
    print("Inertia (positive,negative,zero) = (18,0,39); 57 projective axes, angle arccos(1/5).")
if __name__ == "__main__":
    main(Path(sys.argv[1]))
