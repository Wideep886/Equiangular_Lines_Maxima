# Independent integer-shell check

The input is the original list of low-weight codewords in
`EL18_58_search/integer_frame/codes.txt`. Executing the commands below
from this directory checks all 24,227,280 ordered cuts and all 14,432
cuts remaining after the coordinate trace bounds. Python 3.10 or above
and a C++17 compiler are sufficient.

    c++ -O3 -std=c++17 recheck_cuts.cpp -o recheck_cuts
    ./recheck_cuts < /path/to/integer_frame/codes.txt > survivors.txt
    python3 recheck_fibres.py /path/to/integer_frame/codes.txt survivors.txt

The C++ program first visits a coordinate six-set and distributes all
weight-four and weight-eight codewords over its twenty ordered cuts. For
each admissible cut it checks the rank, support capacities and coordinate
trace inequalities. The nine inputs define doubly even self-dual codes.
The paired code of length 24 needs one generator of weight twelve in
addition to its low-weight words. The checker supplies and tests this
generator explicitly. The original representative generation and the
Pless--Sloane classification identify the nine isometry types.

The Python program classifies each surviving cut individually from its
coordinate incidences. It does not read orbit or component certificates.
It checks the two-coordinate transversal, the four-star by three-pair
decomposition, and the two smaller component types. For the small types
it runs a complete Bron--Kerbosch search on the signed supports.

The surviving counts in code type order as given in the input are
0, 0, 0, 0, 2880, 0, 0, 32 and 11520. The three sign capacities for a
single support, a two-support component and a four-support component are
5, 8 and 10. Sorted as raw UTF-8 records with a newline after each, the
14,432 surviving cuts have SHA-256
`b437a3e1a69dec4e70418da7cf51591710ce0c7278c5d679a65a62f33a8ac2b2`.
