# Weyl orbit midrank bounds

Stanley's theorem gives a bound for any constant-level subset of a regular dominant functional on a Weyl orbit. Let J index a simple system, λ the dominant representative and K its stabilizing simple nodes. The bound is the central coefficient of the parabolic Poincaré polynomial P_J/P_K. In a marked model an equivalent count takes signed shell points x with N₊(x)−N₋(x) in {0,1}, where N₊ and N₋ count positive roots with positive and negative pairing with x in the subsystem perpendicular to both marks. These formulas bound the cardinality of the selected level, without requiring it to coincide with a middle rank.

All coordinates and polynomial coefficients are integers. The rank-24 ordinary-host comparison covers each of the 1,364 dominant polarizations and 2,646 additional adjacent-pair marks. Its 4,010 labelled records agree exactly among the ellipsoid enumeration, an expansion of component orbits, and the archived chain bounds. No chain record enters either of the first two methods. The ellipsoid enumeration imports the original ADE model descriptions, glue and dominant polarization generator. The comparison of exceptional hosts uses the marked models and polarization list of the revision-2 shell reconstruction. The independent second evaluation generates the Coxeter length by graph distance to the dominant point and multiplies local orbit polynomials.

Use Python 3.11+, with assertions on. The alternative expansion for ordinary hosts also needs NumPy. Here J is the extracted joint23 directory, and V2 contains exceptional_shells.py from the revision-2 manuscript supplement. Run from the directory of the new scripts.

```sh
python3 -B poset_direct.py "$J"
python3 -B poset_width.py "$J"
python3 -B poset_records.py "$J" direct direct.tbl
python3 -B poset_records.py "$J" orbit orbit.tbl
python3 -B poset_records.py "$J" historical archived.tbl
cmp direct.tbl orbit.tbl
cmp direct.tbl archived.tbl
python3 -B poset_extra.py "$J" "$V2" extra.tbl
python3 -B poset_extra_orbits.py "$J" "$V2" extra-orbits.tbl
cmp extra.tbl extra-orbits.tbl
```

Expected SHA-256 of the 4,010-row ordinary table: 237772fc57b07e04125aeac4663e153c186353d6503ebb532f2c2ccd2fc2f002. The five columns are host, concatenated 24 dominant Dynkin labels (each is a single digit), two global root nodes specifying the mark (−1,−1 for the unmarked run), and the signed bound.

Expected SHA-256 of the 216-row exceptional table: 81202120db24b98755b27c51eecf98944f5760f6ea8de5bbbc8ca5142c758a54. The columns are coordinate scale, marked factor (−1 for the shortened ternary Golay model), polarization index, signed candidate count and signed bound. There are 15 records with bound at least 116. These counts refer to auxiliary marks; they can overlap at the level of original line systems.

Reference: Richard P. Stanley, *Weyl groups, the hard Lefschetz theorem, and the Sperner property*, SIAM J. Algebraic Discrete Methods 1 (1980), 168–184, Theorem 3.1. https://doi.org/10.1137/0601021
