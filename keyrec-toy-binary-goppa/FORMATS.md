# File formats

All positions and vector indices are zero based. Field elements are integers
0 through 63 in the polynomial basis for `x^6+x+1`, integer modulus 67 (`0x43`).
Polynomial coefficients are stored from lowest to highest degree.

## Public input: `BGPK1`

Whitespace-delimited text:

```text
BGPK1
m t n k modulus heldout d s
column_0 ... column_(n-1)
```

Bit `a` of `column_j` is generator entry `Y[a,j]`. The implementation validates
the target parameters `(m,t,n,k,d,s)=(6,6,64,28,5,4)` and generator rank. Retained
columns must contain all unit vectors. The held-out column is excluded from the
constraint matrix.

`prepare` converts this to the holdout operator format: a header
`k d number_of_points`, followed by `(column_mask, multiplicity)` rows. Unit
columns are omitted because their conditions hold automatically on this
squarefree polynomial space.

## Original audit key: `BGSK1`

```text
BGSK1
m t n modulus
support_0 ... support_(n-1)
g_0 ... g_t
```

The recovery command never reads this format. It is used by generation and the
separate original-secret audit.

## Recovered key: `BGKEY1` and JSON

```text
BGKEY1
m t n modulus
support_0 ... support_(n-1)
lambda_0 ... lambda_(n-1)
g_0 ... g_t
```

The multiplier places the public code in `GRS_(D+1)(support,lambda)`, with
`D=51`. The monic irreducible polynomial `g` defines a binary Goppa code exactly
equal to the public code. The two descriptions agree up to one nonzero global
multiplier scale.

The matching JSON file has `support`, `grs_multiplier`, `goppa_polynomial`, field
parameters, and a `verification` object. The independent verifier computes its
own checks and does not trust the saved verification flags.

## Intermediate matrices: `GF64M1`

```text
GF64M1
rows columns
row_0_entries
...
```

Use the field specified by the corresponding public input. Matrix meanings:

- `branches.txt`: one first coefficient per row, six rows.
- `jet.txt`: coefficient of order `r` in row `r`.
- `quadrics.txt`: one quadratic form per row; monomials `(a,b)` with `a<=b`
  in lexicographic order, including squares.
- `tangents.txt`: one tangent representative per public column.
- `minor-matrix.txt`: rows `(a,b)` with `a<b`, lexicographically ordered;
  entries are the square roots of the minors.
- `direct-sections.txt`: numerator and denominator coefficient rows, with
  the same quadratic monomial order.
- `reconstructed-curve.txt`: one coordinate polynomial per row.
- `coherence-series.txt`: the scalar series and parameter series fitted to
  the normalized local expansion of the reconstructed curve.

## Holdout kernel: `WHK1`

On little-endian hosts, the packed kernel format contains the four-byte magic `WHK1`, four unsigned 32-bit integers `(k,d,N,M)`,
then `M` rows of `ceil(N/8)` bytes. Bit `c%8` of byte `c/8` is coefficient `c`.
The degree-`d` squarefree monomials follow lexicographic order of increasing
index tuples. Here `N=binomial(28,5)=98,280`.

Recovery checks dimensions, complete file length, independence of the `M`
polynomials, held-out values, and the original Hasse constraints. This ties the
kernel to the supplied public matrix without relying on filenames or a cached
success report.
