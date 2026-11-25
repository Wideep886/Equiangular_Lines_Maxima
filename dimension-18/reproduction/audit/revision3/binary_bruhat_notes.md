# Independent check of the parabolic Bruhat bound

Stanley, *Weyl groups, the hard Lefschetz theorem, and the Sperner property*, SIAM J. Algebraic Discrete Methods 1 (1980), Theorem 3.1 and its preceding rank-polynomial formula (2), apply to every parabolic quotient of a Weyl group. The paragraph following Theorem 3.1 includes the product/reducible case. The quotient is rank-symmetric and rank-unimodal and has the Sperner property. https://doi.org/10.1137/0601021

For a point x dominant on J its orbit is indexed by minimal representatives in W_J/W_K, with K its vanishing simple-root pairings in J. In a length-increasing step u→s_βu, both β and u⁻¹β are positive in this subsystem; ux−s_βux is the nonnegative multiple ⟨x,u⁻¹β∨⟩β. Summing in the pointed cone over any Bruhat chain proves order reversal into the root-cone order. Distinct coset endpoints have different images even if individual intermediate steps fix the point. A functional positive on J is strictly different at comparable distinct endpoints. Thus points with a prescribed functional value form an antichain and are bounded by the central rank number. The same functional is used separately on each orbit. Summing orbit bounds is a sufficient upper bound; it does not assert that there is a selection simultaneously attaining every orbit bound.

## Replay

From the directory of the three scripts set J to the extracted joint23 and B to an absolute temporary executable path.

```sh
c++ -O3 -std=c++17 binary_bruhat.cpp -o "$B"
python3 -B binary_bruhat.py "$J" "$B" --compare ./poset_width.py
```

The Python driver generates each complete norm-four slice in simple-root coordinates using the archived shell generator but does not invoke the archived chain verifier. The C++ program builds the simple-reflection graph, chooses each unique parabolic-dominant vertex and traverses the entire orbit by breadth-first search. The shortest graph distance from the dominant vertex to ux is the minimal number of simple reflections among all elements of the coset uW_K, hence gives its Bruhat rank. The program checks disjoint orbit coverage, rank symmetry and central maximality for every orbit. In comparison mode it compares each unmarked row and the complete set of pairs (mark, upper bound) with the polynomial implementation. Standard-output file binary_bruhat_checked.out contains the eighteen host summaries. The entry count is 1364 unmarked polarizations and 2646 marked refinements, of which all are matched. Neither numerical method reads the historical lists of chains. Both rely on the common Niemeier lattice model and the surrounding mathematical reductions.
