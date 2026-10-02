# File formats

Indices are zero based. GF(64) elements are integers 0 through 63 in the
polynomial basis for `u^6+u+1`, modulus 67. GF(4) digits 0, 1, 2, 3 embed as
0, 1, 58, 59. Polynomial coefficients are stored in increasing degree order.

## Public input

```text
WGPK1
m t n k modulus heldout d s
column_0 ... column_(n-1)
```

Here `m=3` is the extension degree over GF(4), `t` is the radical degree,
`D=n-4t-1`, and `s=d-1`. Bits `2a` and `2a+1` of `column_j` contain the GF(4)
digit of generator entry `(a,j)`. The supported inputs have `n<=64` and `k<=31`.

The retained columns must contain an information set of unit vectors. The
generator rank, field encoding, and interpolation bound are checked.

## Packed kernel

`WGK1` is followed by eight little-endian unsigned 32-bit integers:
`k, n, D, heldout, d, monomial_count, exported_count, full_kernel_dimension`.
Each exported row then occupies `ceil(monomial_count/4)` bytes. Coefficient `c`
is stored in bits `2(c mod 4)` and `2(c mod 4)+1` of byte `floor(c/4)`.
Coefficients are GF(4) digits.

Monomials are squarefree increasing index tuples of length `d`, in lexicographic
order. Recovery checks file length, parameters, GF(4) independence, all Hasse
constraints, and seed evaluations. The exported panel and full kernel dimension
are distinct quantities.

## Keys

The separate audit key has this format:

```text
WGSK1
m t n modulus
support_0 ... support_(n-1)
gamma_0 ... gamma_t
```

Recovered text keys use:

```text
WGKEY1
m t n modulus
support_0 ... support_(n-1)
lambda_0 ... lambda_(n-1)
gamma_0 ... gamma_t
```

The multiplier describes the ambient `GRS_(D+1)` code. The monic irreducible
radical describes the public code as both `Gamma(support,gamma^3)` and
`Gamma(support,gamma^4)`.

JSON keys contain `q`, `m`, `t`, `n`, `k`, `modulus`, `support`,
`grs_multiplier`, `goppa_radical`, and `goppa_exponent=3`. Verifiers recompute the
checks rather than trusting the saved `verification` object.

## Matrices

```text
GF64M1
rows columns
row_0_entries
...
```

- `branches.txt`: one tangent representative per row.
- `jet.txt`: coefficient of order `r` in row `r`.
- `quadrics.txt`: one quadratic form per row, with monomials `(a,b)`, `a<=b`,
  in lexicographic order, including squares.
- `tangents.txt`: one tangent representative per public column.
- `minor-matrix.txt`: square-root minor rows indexed by `(a,b)`, `a<b`.
- `direct-sections.txt`: numerator and denominator rows in quadratic monomial order.
- `reconstructed-curve.txt`: one coordinate polynomial per row.
- `coherence-series.txt`: scalar and parameter series used by the coefficient audit.
