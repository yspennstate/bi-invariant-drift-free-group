# Positive drift for a bi-invariant metric on a free group

Yitzchak Shmalo

This repository contains the paper and the computations for the following result. Let r >= 3, let mu be the uniform measure on the 2r letters of the free group F_r, and let L(g) be the cancelation norm: the least number of letters that must be deleted from a word representing g so that the remaining word represents the identity. Then d(g, h) = L(g^{-1} h) is a bi-invariant metric on F_r, every letter has length one, and the drift of the mu-random walk with respect to d is at least c_r > 0, where c_r is the root of

    -c log c - (1 - c) log(1 - c) = (1 - c) log sqrt(r/2)

in (0, 1/(1 + sqrt(2/r))). For example c_3 = 0.048, c_8 = 0.227 and c_32 = 1/2. This gives a negative answer to Problem 4 of M. Bjorklund, *Five remarks about random walks on groups* (arXiv:1406.0763), which asks whether the drift of a symmetric random walk vanishes for every bi-invariant metric with finite first moment.

The proof shows that L does not depend on the representing word, by a surgery on noncrossing matchings of inverse letters, and then bounds the probability that the random word of length n admits a noncrossing matching leaving at most cn letters unpaired: a fixed matching with m pairs succeeds with probability (2r)^{-m}, and there are at most binom(n, k) 4^m of them.

## Contents

- `paper/main.tex`, `paper/main.pdf`: the paper.
- `code/checks.py`: the dynamic program for L against the definition (all words of length at most 7 in F_2); reduction invariance, nondegeneracy, inverse symmetry, subadditivity and conjugation invariance on small words; the exact tail sum for r = 32 against 2^{-n/2}; the constants c_r; Monte Carlo estimates of E L(X_n)/n for r = 2, 3, 32. Output in `code/checks_out.txt` and `code/checks.json`.
- `lean/ExponentCheck.lean`: the arithmetic step k + 2m = n, 4k <= n implies 3n <= 8m used in the case r = 32, checked with Lean 4 and Mathlib.

## Running

    cd code
    python checks.py > checks_out.txt

The program uses only the Python standard library and runs in under a minute. The Lean file compiles inside any Lean 4 project with Mathlib (`lake env lean ExponentCheck.lean`).

## License

Paper: CC BY 4.0. Code: MIT.
