# TII recovered holdout polynomials

These are the polynomial samples used in the completed TII-173, TII-249, and
TII-253 recoveries. Each sample is at original column **0**. The rows are
independent recovered polynomials, not a basis of the entire holdout kernel.
Polynomial coefficients lie in **GF(2)**; the extension field used in later
recovery is specified by `public.json:field_modulus`.

| Instance | Shortened code | Variables | Degree | Multiplicity away from holdout, including infinity | Polynomials | Monomials | Format |
|---|---|---:|---:|---:|---:|---:|---|
| [TII-173](tii-173/README.md) | `[91,35]` | 35 | 5 | 4 | 256 | 324,632 | WHK1, lexicographic |
| [TII-249](tii-249/README.md) | `[178,50]` | 50 | 5 | 4 | 280 | 2,118,760 | THK1, colex |
| [TII-253](tii-253/README.md) | `[118,46]` | 46 | 6 | 5 | 280 | 9,366,819 | THK1, colex |

`public.json` supplies the exact generator row basis, column maps, field
modulus, and constraint points. Variable `X_a` is row `a` of that generator.
Do not substitute a row-equivalent generator without transforming the forms.
`original_indices[j]` maps shortened column `j` to the original public matrix.
`heldout` and `information_set` use shortened indices; fields ending in
`_original` use original indices. All indices are zero based.

## Canonical packed storage

Both formats have a 20-byte header; all integers are unsigned little-endian.

| Byte offset | Size | Meaning |
|---:|---:|---|
| 0 | 4 | ASCII magic `WHK1` or `THK1` |
| 4 | 4 | Number of variables `k` |
| 8 | 4 | Homogeneous squarefree degree `d` |
| 12 | 4 | `N = binomial(k,d)` monomials |
| 16 | 4 | Number of saved polynomial rows `q` |

The payload contains `q` consecutive rows of `ceil(N/8)` bytes each.
Coefficient `c_j` is bit `j % 8` of byte `j // 8`, with the least significant
bit first. It multiplies the corresponding squarefree monomial
`X_a1 ... X_ad`. Unused bits in the last byte of each row are zero. The exact
file size is `20 + q*ceil(N/8)`; there is no row alignment or trailing data.

- **WHK1:** tuples `(a1,...,ad)` are in lexicographic order, exactly
  `itertools.combinations(range(k), d)`.
- **THK1:** masks with `d` set bits are in increasing numerical order (colex).
  For `a1 < ... < ad`, the zero-based rank is
  `sum(binomial(ai,i), i=1..d)`.

Changing the magic does not convert the monomial order. TII-173 retains the
WHK1 bytes read by its recovery and also includes the original THK1 sample
for the standalone CPU Hasse verifier. These are the same polynomial rows
with permuted monomial columns; the verifier's `export` command reproduces
the WHK1 file exactly. TII-249 and TII-253 retain their original THK1 bytes.
No polynomial coefficients were recomputed for this collection.

## Chunks and reconstruction

Files larger than 32 MiB are stored as consecutive raw chunks named
`kernel.thk1.part000`, `part001`, etc. Chunk boundaries have no mathematical
meaning and can fall inside a polynomial. Concatenating them in numerical
order restores the original file, including its single header.
`storage.json` gives every chunk's offset, byte length and SHA-256, the
complete file's SHA-256, and the hash of its public coordinate basis.
For TII-173, `alternate_encodings` also authenticates `kernel.thk1`.
The collection-wide `SHA256SUMS` also covers the metadata and documentation;
check it with `shasum -a 256 -c SHA256SUMS` from this directory.
Finder metadata (`.DS_Store`) and generated files are excluded from the manifest.

From this directory, verify the stored chunks without creating files:

```sh
python3 materialize.py tii-173
python3 materialize.py tii-249
python3 materialize.py tii-253
```

To reconstruct a saved input in a fresh directory:

```sh
python3 materialize.py tii-253 --output tii253-input
```

The output is `kernel.thk1` and `public.json`. TII-173 also writes `kernel.bin`
in WHK1 order for its recovery code. This is file conversion/copying only;
it does not run recovery. Each [implementation README](../tii-keyrec-implementation/README.md)
includes optional commands for auditing its saved kernel and intermediate evidence.

## Derived panel used by the implementations

TII-249 and TII-253 use a second representation during jet extraction:
`polynomials/kernel-panel.u64le`. It has no header and contains `N` rows in
colex monomial order, each consisting of **eight** little-endian `uint64`
words (64 bytes). Bit `b` of word `w` in row `j` is the coefficient of
monomial `j` in polynomial `64*w+b`. Only the first 280 polynomial columns
are populated; columns 280 through 511 are zero. This is a transpose of the
THK1 coefficient matrix plus zero padding, not a second set of polynomials.

`monomials-colex.u64` contains one little-endian 64-bit variable mask per
monomial, without a header. Bit `a` means variable `X_a` occurs. Both large
derived files are omitted from storage because they are exactly recoverable.
Their original byte sizes and SHA-256 hashes are in `storage.json`.
`polynomial-index.json` is the original index; its paths describe the restored
layout under `polynomials/`.

With NumPy (also included in the Sage environment), restore both derived
files and check them against their original hashes:

```sh
sage -python materialize.py tii-253 --output tii253-panel-input --panel
```

## Public constraints and verification

`operator.txt` starts with `TII_HOLDOUT_V1`, followed by labeled fields
`instance_id`, `k`, `degree`, `minimum`, `heldout_original`, `heldout_mask`,
and `points`. Each following point is `original_label binary_mask multiplicity`;
label `-1` denotes infinity. Bit `a` of a point mask is coordinate `a`.
There are no constraints at the held-out point. The native verifier checks
every Hasse order below each prescribed multiplicity, including orders omitted
from the reduced solver. `request.txt` for the larger instances is the exact
reduced GPU request; its `point MASK COUNT ORDERS...` entries retain orders
1 and 3 for TII-249 and 2 and 4 for TII-253.

`verification.json` is the saved independent membership/rank report.
Its `curve_identity_certified: false` describes the holdout-stage scope;
it does not supersede the later recovery and curve audits in
[tii-results](../tii-results/README.md). The implementation and direct commands
are in [tii-keyrec-implementation](../tii-keyrec-implementation/README.md).
