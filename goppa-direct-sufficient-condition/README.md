# Existence condition for binary direct recovery

This collection checks the sufficient condition for direct key recovery of
binary Goppa codes on generated keys whose secret curve is known. It covers
Classic McEliece keys shortened to the dimensions used by the distinguisher.

## The condition

Let `C = Γ(α, G)` be a binary square-free Goppa code of length `n` with
`deg G = t` over `K = GF(2^m)`. Then `C` lies in a GRS code of degree bound
`D = n − 2t − 1`. Let `f_1, ..., f_k` be the polynomials of degree at most `D`
that correspond to an F₂-basis of `C`, and put

```
W = span_K { f_a f_b : 1 ≤ a ≤ b ≤ k }.
```

Every product satisfies `G² | (f_a f_b)'`, so `dim W ≤ 2D − t + 1`. The
proposition on the existence of the forms (Section 7.4 of the revised paper)
states that equality implies that `(Z − α_i)^{2D}` and `(Z − α_i)^{2D−2}` lie
in `W` for every position `i`. Then the quadratic forms `R` and `S` used by the
binary direct recovery exist at every position.

The program generates a key, computes `dim W` from the secret curve, and stops
as soon as the rank reaches the upper bound `2D − t + 1`.

## Method

1. **Keys.** The Goppa polynomial is a random monic irreducible polynomial of
   degree `t`. The support consists of `n` random distinct elements of `K`. For
   `m = 12` and `m = 13`, `K` uses the Classic McEliece field polynomials.
2. **Shortening.** The distinguisher shortens the public code at `k − k_l`
   positions of an information set. Shortening a Goppa code gives the Goppa
   code with the same polynomial on the retained support. The Classic
   McEliece support order is random, so the retained support is a random
   subset of size `n = mt + k_l`. The shortened code therefore depends only on
   `(m, t, k_l)`; `mceliece6688128` and `mceliece8192128` give the same
   distribution. The program still draws different keys for them, because the
   full length enters the seed.
3. **Curve.** An F₂-basis of the code comes from the kernel of the binary
   parity-check matrix. It is systematic on an information set. The values of
   the curve polynomials are
   `f_a(α_j) = c_{a,j} Π'(α_j) / G(α_j)²` on the support, and
   `f_a(β) = Π(β) / G(β)² · Σ_j c_{a,j} / (β − α_j)` at further points `β`.
   Here `Π` is the support polynomial.
4. **Rank.** Each product is represented by its values at `N = 2D + 1`
   distinct points of `K`, which is injective on polynomials of degree at most
   `2D`. When `K` has fewer than `2D + 1` points, as for unshortened codes,
   the curve polynomials are interpolated and each product is stored by its
   coefficients. Products are added to a semi-echelon basis until the rank
   reaches `2D − t + 1` or every product has been used.
5. **Order of the products.** The `k` squares `f_a²` come first, then the
   products `f_a f_b`, `a < b`, each group in random order. The basis is
   systematic, so the unit vectors are public columns. Every quadratic form
   that vanishes on the curve vanishes at them, so its coefficients at
   `X_a²` are zero. Hence the products with `a < b` span a subspace of
   codimension `k` in `W`, and each square is needed. With this order the bound
   is reached after exactly `2D − t + 1` products.
6. **Membership.** At sampled positions `i`, the program tests whether
   `(Z − α_i)^{2D}` and `(Z − α_i)^{2D−2}` lie in `W`. As a control,
   `(Z − α_i)^{2D−1}` must not lie in `W`, since its derivative is not divisible
   by `G²`.

## Checks

Every key also passes these checks:

- each basis codeword satisfies `Σ_j c_j / (Z − α_j) ≡ 0 mod G²`, which makes
  `f_c / G²` a polynomial of degree at most `D`;
- interpolation of sampled values confirms `deg f_a ≤ D` and `G² | (f_a f_b)'`;
- the control power is never in `W`, and both target powers are in `W` whenever
  the bound is reached.

For small codes every product is used, so the upper bound `dim W ≤ 2D − t + 1`
is tested directly. The self-test compares the rank engine with plain dense
elimination, and recomputes each key in coefficient form to compare ranks and
memberships between the two representations. An independent SageMath script
recomputes small and unshortened keys from their coefficients and compares
ranks and memberships.

## Results

The batches ran on a 96-core AMD EPYC machine under Linux, one key per core.
The SageMath cross-check ran under macOS on an Apple M1 Pro. Rerunning a key on
either platform reproduces its JSON line exactly, apart from timings.
[results/summary.md](results/summary.md) has the complete tables, and
`results/*.jsonl.gz` has one JSON line per key.

All 165,400 keys pass every check listed above. The rank never exceeded
`2D − t + 1`, and the control power never lay in `W`. The cross-check
recomputed 67 keys, 20 of them deficient and 5 unshortened, and agreed on every
rank and on all 2,670 membership tests.

### Classic McEliece keys

The shortened dimensions `k_l` are those of the distinguisher, from the cost
table of the paper: GIAJS use `(d, s) = (8, 7)`, Weis uses `(d, s) = (7, 6)`.
Each row has 1,000 keys and 8 tested positions per key.

| Parameter set | Shortening | n | k | 2D − t + 1 | C(k+1, 2) | dim W = 2D − t + 1 | Positions with R and S | Median s per key |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| `mceliece348864` | Weis | 922 | 154 | 1523 | 11935 | 1000 / 1000 | 8000 / 8000 | 0.7 |
| `mceliece460896` | Weis | 1443 | 195 | 2405 | 19110 | 1000 / 1000 | 8000 / 8000 | 2.8 |
| `mceliece6688128` | Weis | 1888 | 224 | 3135 | 25200 | 1000 / 1000 | 8000 / 8000 | 6.2 |
| `mceliece6960119` | Weis | 1763 | 216 | 2930 | 23436 | 1000 / 1000 | 8000 / 8000 | 5.1 |
| `mceliece8192128` | Weis | 1888 | 224 | 3135 | 25200 | 1000 / 1000 | 8000 / 8000 | 6.2 |
| `mceliece348864` | GIAJS | 985 | 217 | 1649 | 23653 | 1000 / 1000 | 8000 / 8000 | 0.9 |
| `mceliece460896` | GIAJS | 1522 | 274 | 2563 | 37675 | 1000 / 1000 | 8000 / 8000 | 3.4 |
| `mceliece6688128` | GIAJS | 1979 | 315 | 3317 | 49770 | 1000 / 1000 | 8000 / 8000 | 7.2 |
| `mceliece6960119` | GIAJS | 1851 | 304 | 3106 | 46360 | 1000 / 1000 | 8000 / 8000 | 6.1 |
| `mceliece8192128` | GIAJS | 1979 | 315 | 3317 | 49770 | 1000 / 1000 | 8000 / 8000 | 6.1 |

Every key reached the bound, so the forms `R` and `S` exist at every position
of every key. The shortened code had dimension exactly `k_l` for every key.

### The counting threshold

The bound can be reached only if `C(k+1, 2) ≥ 2D − t + 1`. At the Classic
McEliece values of `(m, t)`, this needs `k ≥ 51` for `(12, 64)` and `k ≥ 75`
for `(13, 128)`. For 20 keys at each `k` from 47 to 56 and from 71 to 80,
`dim W = min(C(k+1, 2), 2D − t + 1)` held for all 400 keys. At these
parameters the count is the only obstruction observed. The distinguisher
dimensions are about 3 to 4 times the threshold for their `(m, t)`.

### Moderate shortened codes

These runs fill the range between the small codes and the Classic McEliece
shortenings. For `(m, t)` in `(9, 16)`, `(10, 20)`, `(10, 40)`, `(11, 32)`,
`(12, 32)` and `(13, 64)`, `k` runs from the counting threshold to 100 above
it, with 500 keys for each `k`.

| Outcome | Keys |
|---|---:|
| All moderate keys | 33,000 |
| `dim W = 2D − t + 1` | 32,993 |
| Deficient, one double point | 5 |
| Deficient at zero slack, one dependency more than predicted | 2 |

The five keys with a double point have `(m, t) = (9, 16)` and `k` from 22 to
25. The two zero-slack keys have `(m, t) = (10, 20)` and `k = 26`. No key with
`t ≥ 32` failed, and no key failed more than three steps above the threshold.
On every key that reached the bound, `R` and `S` existed at all tested
positions.

### Unshortened codes with the TII parameters

These codes have the length and Goppa degree of the TII challenge codes, but
random keys. Here `2D + 1` exceeds the field size, so the products are
computed in coefficient form.

| Parameters | (m, t) | n | k | 2D − t + 1 | C(k+1, 2) | dim W = 2D − t + 1 | Positions with R and S |
|---|---|---:|---:|---:|---:|---:|---:|
| TII-173 | (7, 8) | 96 | 40 | 151 | 820 | 1000 / 1000 | 8000 / 8000 |
| TII-249 | (8, 16) | 235 | 107 | 389 | 5778 | 1000 / 1000 | 8000 / 8000 |
| TII-253 | (8, 9) | 214 | 142 | 382 | 10153 | 1000 / 1000 | 8000 / 8000 |

### Small codes

For small codes every product is used. The runs cover `(m, t)` from `(6, 3)`
to `(11, 16)`, with `k` from a few values below the counting threshold to
values where no key fails, and 1,000 keys for each `k`.

| Outcome | Keys |
|---|---:|
| All small keys | 119,000 |
| Above the counting threshold | 98,008 |
| `dim W = 2D − t + 1` | 94,065 |
| Deficient, deficit equal to the number of conditions from double and base points | 3,917 |
| Deficient at zero slack, one dependency more than predicted | 26 |

Below the threshold, the products were independent for 20,988 of 20,992 keys.

**Double and base points.** The curve is the image of the projective line,
including infinity. The leading coefficient of `f_a` is the weight of the
codeword `c_a` modulo 2, so `F(∞)` is the parity column. Consider the generator
matrix of the basis extended by this column.

- **Double point.** Two equal nonzero columns `i` and `j` give
  `λ_i² U(α_i) = λ_j² U(α_j)` for every `U` in `W`, which is one linear
  condition. At infinity, the leading coefficient of `U` takes the place of
  `U(α_i)`.
- **Base point.** A zero column means that every `f_a` vanishes there. Then
  `U` and `U'` vanish there for every `U` in `W`, which is two conditions. If
  every codeword has even weight, the base point is at infinity.

In every deficient key with positive slack, the deficit equals the number of
these conditions. Among the 3,917 explained keys, 3,645 have special points
only on the support, 186 only at infinity, and 95 at both. A degenerate tangent
would add one more condition, but none occurred in the deficient keys, neither
on the support nor at infinity. Keys without special points reached the bound
whenever the slack was positive.

Both kinds of special points are visible in the public code: they are repeated
or zero columns of the public generator matrix extended by its parity column.
For uniformly random columns their probability is about `C(n+1, 2) / 2^k`,
which matches the observed rates in the summary. For the shortened Classic
McEliece codes above it is below `2^-135`.

**Zero slack.** When `C(k+1, 2) = 2D − t + 1`, every product must be
independent. The 26 remaining keys all have zero slack and exactly one extra
dependency. Among keys without special points, the rate of such failures is of
the same order as the probability that a random square matrix over `K` is
singular.

| (m, t) | k | Keys without special points | Failing | Rate | 1/(2^m − 1) |
|---|---:|---:|---:|---:|---:|
| (6, 3) | 8 | 297 | 11 | 3.7% | 1.6% |
| (7, 4) | 10 | 453 | 3 | 0.66% | 0.79% |
| (8, 5) | 12 | 718 | 3 | 0.42% | 0.39% |

On deficient keys the forms usually fail at every position. Only 138 of the
3,943 deficient keys had `R` and `S` at some tested position.

## Build and run

Requires a C++17 compiler and Python 3. No other libraries are needed.
SageMath 10.7 or later is optional and is used only for the cross-check.

```sh
make
make check
```

`make check` runs the self-test. Single keys are checked as follows:

```sh
./build/product_rank presets
./build/product_rank preset mceliece348864-weis --seeds 1-10 --positions 8
./build/product_rank run --m 6 --t 3 --k 8-12 --seeds 1-100 --full
./build/product_rank run --m 8 --t 16 --n 235 --seeds 1-10 --positions 8
```

Each key prints one JSON line. `--full` uses every product instead of stopping
at the bound. `--n` sets the length directly, and `--coeff` forces the
coefficient form, which is otherwise chosen when `2D + 1` exceeds the field
size. Deficient keys with `n ≤ 400` also report the special points of
the curve described above.

### Reproducing the results

`run.sh` builds the program, runs the self-test and every parameter set, and
writes `results/`. Its defaults are the settings of the results above. `JOBS`
sets the number of parallel processes and defaults to the number of CPUs. The
complete run takes about 16 CPU-hours, about ten minutes on 96 cores.

```sh
./run.sh
```

When SageMath is installed, `run.sh` also runs `crosscheck.sh`, which dumps
small keys and rechecks them independently. `summarize.py` turns the JSON lines
into [results/summary.md](results/summary.md). Keys are determined by their
seeds, so a rerun reproduces every JSON line apart from timings.
[results/build-info.json](results/build-info.json) records the compilers and
the source hash of the published results.

## Files

| File | Content |
|---|---|
| [src/product_rank.cpp](src/product_rank.cpp) | Key generation, curve values, rank engine, checks |
| [run.sh](run.sh) | Complete experiment |
| [crosscheck.py](crosscheck.py), [crosscheck.sh](crosscheck.sh) | Independent SageMath check |
| [summarize.py](summarize.py) | Summary tables |
| [results/](results/) | JSON lines per key, summary, cross-check, build information |
