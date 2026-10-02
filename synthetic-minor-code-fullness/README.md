# Fullness of the minor code

These sweeps test whether the minor code of a public code is a full GRS code
with the support of the public code, and whether the Sidelnikov–Shestakov
attack applied to it gives a support on which the public code embeds into a
GRS code of the expected dimension. The tangent planes are computed from the
secret key (hence *synthetic*), so the sweeps test the minor code alone.

The sweep covers 54 parameter sets in characteristic two: GRS subcodes,
alternant codes, wild Goppa codes and binary Goppa codes, including Classic
McEliece codes shortened to the dimensions used by the distinguisher. Each of
the 44 smaller sets has 100 independently seeded codes, and each of the 10
Classic McEliece sets has 20, for 4,600 codes in total.

## Construction

Let `C ⊆ GRS_{D+1}(α, λ)` be a public `[n, k]` code over the field `K` of the
ambient GRS code, with reduced generator matrix `Y` and columns `y_j`. The rows
of `Y` are `y_{a,j} = λ_j f_a(α_j)` with `deg f_a ≤ D`, and
`F = (f_1, …, f_k)`.

**Tangent planes.** With `Π(Z) = ∏_j (Z − α_j)` and `ν_j = (λ_j Π'(α_j))^{-1}`,
differentiating the Lagrange formula for `F` gives

```
ŷ_j = Σ_{l ≠ j} ν_l y_l / (α_j − α_l) ∈ F'(α_j) / Π'(α_j) + K·F(α_j),
```

so `V_j = span{y_j, ŷ_j}` is the tangent plane `span{F(α_j), F'(α_j)}`
whenever it is two-dimensional. Each representative is then replaced by
`a_j ŷ_j + b_j y_j` with random `a_j ≠ 0` and `b_j`. Positions where
`y_j, ŷ_j` span less than a plane are recorded as tangent failures.

**Minor code.** The minor matrix `M` has rows indexed by pairs `a < b` and
entries

```
m_{(a,b),j} = y_{a,j} ŷ_{b,j} − y_{b,j} ŷ_{a,j}.
```

The matrix `T` is obtained by taking the square root of every entry, and the
minor code `𝒯` is the row space of `T`.

**Upper bound.** Every row of `T` lies in `GRS_κ(α, η)`, where

- `κ = D` and `η_j = (λ_j a_j / Π'(α_j))^{1/2}` in general;
- `κ = D − t` and `η_j` multiplied by `G(α_j)` for binary square-free Goppa
  codes `Γ(α, G)` with `deg G = t`;
- `κ = D − t` and `η_j` multiplied by `g(α_j)` for wild Goppa codes
  `Γ(α, g^{q−1})` over `F_q` with `K = F_{q^2}` and `g` without roots in `K`.
  Here `C ⊆ GRS_{D−t+1}(α, (λ_j g(α_j))_j)`, which is checked for every code.

**Certificate of the rank.** The program forms the rows of a random set of
pairs and checks that each lies in `GRS_κ(α, η)`. If the rank of these rows is
`κ`, then the complete minor code, with all `C(k, 2)` rows, equals
`GRS_κ(α, η)`. Otherwise the set of pairs is doubled until the rank reaches
`κ` or all pairs are used. The recorded rank is therefore the dimension of the
complete minor code. The code is *full* if this dimension is `κ`, or if `𝒯`
equals a GRS code of its dimension with support `α`.

**Recovery.** Recovery reads only `Y`, `T` and `D`:

1. Sidelnikov–Shestakov on a reduced generator of the smaller of `𝒯` and its
   dual gives a support `α'`, with two coordinates fixed to 0 and 1 and every
   value of a third coordinate tried.
2. The multiplier of `𝒯` on `α'` is read from the Cauchy entries of its reduced
   generator and checked against every row.
3. A multiplier `λ'` of `C` solves `Y diag(z) H^T = 0`, where `H` is a
   parity-check matrix of `RS_{D+1}(α')`, and `λ'_j = z_j^{-1}`.

The embedding `C ⊆ GRS_{D+1}(α', λ')` is verified on all rows, once by the
recovery code and once independently (Sage for the GRS, alternant and wild
Goppa families; NumPy field arithmetic for binary Goppa codes). For shortened
codes the embedding covers the retained coordinates.

The secret support and multiplier are used to generate the codes and the
tangent planes, and in the containment checks with `η`. Recovery receives no
secret data.

## Code families

| Family | Sets | Code | `D` | Bound `κ` |
|---|---:|---|---|---|
| GRS subcodes | 7 | random `k`-dimensional subcode of `GRS_{D+1}(α, λ)` over `F_q`, `q = 32, …, 512` | given | `D` |
| Alternant | 11 | `A_r(α, λ) = GRS_{n−r}(α, λ) ∩ F_q^n`, `q = 2, 4, 8, 16`, `K = F_{q^m}` | `n − r − 1` | `D` |
| Wild Goppa | 15 | `Γ(α, g^{q−1}) = Γ(α, g^q)` over `F_q`, `q = 4, 8, 16`, `deg g = t` | `n − qt − 1` | `D − t` if `m = 2` and `g` irreducible of degree `t ≥ 2`, else `D` |
| Binary Goppa | 21 | `Γ(α, G)` over `F_2`, `K = F_{2^m}`, `deg G = t` | `n − 2t − 1` | `D − t` |

Supports are random subsets of `K` avoiding the roots of the Goppa polynomial;
multipliers are uniformly random nonzero elements. Goppa polynomials are
uniformly random monic irreducible polynomials, products of `t` distinct
random linear factors (`split`), or products of random irreducible
polynomials of distinct degrees (`separable`). The full parameter list is
printed by `sage -python src/families.py`.

**Binary Goppa codes and shortening.** The C++ program expands the parity
check `α_j^i / G(α_j)`, `0 ≤ i < t`, over `F_2`, reduces it, and shortens the
code at the first `ℓ` non-pivot positions, which belong to an information
set. The public matrix is the reduced generator of the shortened code. For
`m = 12, 13` the field polynomials are those of Classic McEliece. The codes
are generated by this program, not by the Classic McEliece key generator.

**Classic McEliece shortenings.** Each parameter set is shortened to the two
dimensions `k` used by the distinguisher with `(d, s) = (7, 6)` and
`(d, s) = (8, 7)`; the length after shortening is `n = k + mt`.

| Parameter set | `(m, n₀, t)` | `ℓ` | `n` | `k` | `D` | `D − t` | `C(k, 2)` |
|---|---|---:|---:|---:|---:|---:|---:|
| `mceliece348864` | `(12, 3488, 64)` | 2566 | 922 | 154 | 793 | 729 | 11781 |
| | | 2503 | 985 | 217 | 856 | 792 | 23436 |
| `mceliece460896` | `(13, 4608, 96)` | 3165 | 1443 | 195 | 1250 | 1154 | 18915 |
| | | 3086 | 1522 | 274 | 1329 | 1233 | 37401 |
| `mceliece6688128` | `(13, 6688, 128)` | 4800 | 1888 | 224 | 1631 | 1503 | 24976 |
| | | 4709 | 1979 | 315 | 1722 | 1594 | 49455 |
| `mceliece6960119` | `(13, 6960, 119)` | 5197 | 1763 | 216 | 1524 | 1405 | 23220 |
| | | 5109 | 1851 | 304 | 1612 | 1493 | 46056 |
| `mceliece8192128` | `(13, 8192, 128)` | 6304 | 1888 | 224 | 1631 | 1503 | 24976 |
| | | 6213 | 1979 | 315 | 1722 | 1594 | 49455 |

## Results

**4,599 of the 4,600 codes pass.** For each of them the complete minor code
is a full GRS code of dimension `κ` with support `α`, Sidelnikov–Shestakov
recovers a support `α'` from it, and the public code embeds into
`GRS_{D+1}(α', λ')`. Whenever all tangent planes were two-dimensional, the
minor code was full and recovery succeeded. No trial raised an error or timed
out.

| Family | Sets | Codes | Pass | Dimension of `𝒯` |
|---|---:|---:|---:|---|
| GRS subcodes | 7 | 700 | 700 | `D` |
| Alternant | 11 | 1,100 | 1,099 | `D` |
| Wild Goppa | 15 | 1,500 | 1,500 | `D − t` for the 5 sets with `m = 2` and irreducible `g`, `t ≥ 2`; `D` for the other 10 |
| Binary Goppa | 11 | 1,100 | 1,100 | `D − t` |
| Classic McEliece, shortened | 10 | 200 | 200 | `D − t` |

**The exception.** In trial 85 of `alternant-a3` (`n = 60`, `r = 8`,
`F_2 ⊂ F_64`, `k = 12`), the vectors `y_j` and `ŷ_j` span only a line at
position 40 (counted from 0). Column 40 of `T` is zero, and the complete minor
code, computed from all 66 pairs, has dimension 50 instead of 51. The trial is
recorded as a tangent failure.

**Classic McEliece.** All 200 shortened codes give a minor code of dimension
`D − t`. The last column is the median elapsed time per code, with 31 codes
processed in parallel:

| Set | `n` | `k` | `D − t` | Pass | Median s |
|---|---:|---:|---:|---:|---:|
| `mceliece348864-k154` | 922 | 154 | 729 | 20/20 | 3.0 |
| `mceliece348864-k217` | 985 | 217 | 792 | 20/20 | 3.6 |
| `mceliece460896-k195` | 1443 | 195 | 1154 | 20/20 | 17.7 |
| `mceliece460896-k274` | 1522 | 274 | 1233 | 20/20 | 20.7 |
| `mceliece6688128-k224` | 1888 | 224 | 1503 | 20/20 | 43.1 |
| `mceliece6688128-k315` | 1979 | 315 | 1594 | 20/20 | 53.2 |
| `mceliece6960119-k216` | 1763 | 216 | 1405 | 20/20 | 34.4 |
| `mceliece6960119-k304` | 1851 | 304 | 1493 | 20/20 | 42.7 |
| `mceliece8192128-k224` | 1888 | 224 | 1503 | 20/20 | 42.3 |
| `mceliece8192128-k315` | 1979 | 315 | 1594 | 20/20 | 53.2 |

[results/summary.md](results/summary.md) lists every set with its
parameters, ranks and failure counts. These are counts of passing codes, not
estimates of a probability: with 100/100 passes the two-sided Wilson 95%
interval for the per-set success rate starts at 96.3%, and with 20/20 at 83.9%.

### Data

| File | Content |
|---|---|
| [results/summary.md](results/summary.md), [results/summary.json](results/summary.json) | Per-set parameters, ranks, statuses, Wilson intervals, timings |
| [results/manifest.json](results/manifest.json) | Master seed, parameters, trial counts, SageMath version, source hashes |
| [results/raw-data.tar.gz](results/raw-data.tar.gz) | One JSON record per code, and a fixture with its generated key and public-matrix hash |
| [results/verification.json](results/verification.json) | Output of `verify_results.py`: coverage, seeds, fixture hashes, status flags |
| [results/validation.log](results/validation.log) | Output of `validate.py` |
| [results/SHA256SUMS](results/SHA256SUMS) | Checksums of these files |

To check the records:

```sh
mkdir -p runs/check && tar -xzf results/raw-data.tar.gz -C runs/check
python3 src/verify_results.py runs/check --require-complete
python3 src/summarize.py runs/check
```

## Building and running

**Requirements.** SageMath 10.7, which provides Python 3, NumPy and the M4RI
and M4RIE libraries with their headers; a C++14 compiler; a POSIX shell. The
results were produced with SageMath 10.7 from conda-forge, which can be
installed with

```sh
conda create -n sage -c conda-forge sage=10.7
conda activate sage
```

**Build.** From this directory, run

```sh
sh build.sh
```

This compiles [src/native.cpp](src/native.cpp) into `src/libminor_native.so`
against the M4RIE and M4RI libraries of the SageMath installation that `sage`
on `PATH` belongs to. Set `SAGE_PREFIX` to use another installation prefix and
`CXX` to choose the compiler.

**Validate.** The checks below take about 4 minutes on one core; each prints a
`PASS` line, as in [results/validation.log](results/validation.log).

```sh
sage -python src/validate.py
```

**Run.** The command

```sh
sage -python src/run.py --output runs/main --workers 8
```

runs all 54 sets with the master seed `minor-code-fullness`: 100 codes per set
and 20 per Classic McEliece set. The recorded run used about 2.7 thread-hours,
1.7 of them for the Classic McEliece sets; with 31 parallel workers it took
about 5 minutes. Each code is written to `runs/main/cases/<set>/<trial>.json`
and its fixture to `runs/main/fixtures/`. Repeating the command resumes the
missing trials; `run.py` refuses to mix different sources or settings in one
output directory. Further options:

| Option | Meaning |
|---|---|
| `--labels a,b,...` | Only these sets; `sage -python src/families.py` lists them |
| `--trials N`, `--mceliece-trials N`, `--trial-start N` | Trial numbers per set (default 100 and 20, from 0) |
| `--workers N` | Parallel processes (default: number of CPUs minus one) |
| `--master-seed S` | Another family of seeds |
| `--backend sage` | Sage reference for the minor code of the non-binary families |
| `--explicit-all-pairs` | Non-binary families: also check every pair against the parity check when the bound is attained |
| `--no-recovery` | Skip support recovery |

For example, a single code:

```sh
sage -python src/run.py --output runs/one --labels mceliece348864-k217 \
    --mceliece-trials 1 --trial-start 17 --workers 1
```

**Summarize and check.**

```sh
python3 src/summarize.py runs/main                     # writes summary.md and summary.json
python3 src/verify_results.py runs/main --require-complete
```

**Seeds.** SHA-256 of `master|label|trial|phase` gives a 128-bit seed for each
phase: the code, the tangent representatives, the sample of pairs and the
recovery. Sage uses the full seed; the C++ program seeds `std::mt19937_64`
with its low 64 bits. Binary Goppa codes therefore depend only on the seed,
since `std::mt19937_64` is specified by the C++ standard and the program uses
no library distributions. The other families use the random generator of
SageMath and need version 10.7 to regenerate the same codes. Failed trials are
recorded and never redrawn.

**Validation.** `src/validate.py` checks:

- table-based field arithmetic, the GRS parity-check formula and the multiplier
  recovery against Sage, including rejection of corrupted GRS matrices;
- irreducibility of the binary field polynomials and the NumPy field
  arithmetic against Sage;
- eight binary Goppa trials of the C++ program against an independent Sage
  reconstruction: the factorization of `G`, the code obtained from the
  reference construction with `G²` and its shortening, the tangent planes
  against derivatives of the interpolated polynomials `f_a`, the rank of the
  complete minor code over all pairs, and the recovered embedding;
- the C++ minor-code computation against Sage with an explicit parity
  certificate over all pairs, for GRS, alternant and wild Goppa codes;
- the tangent formula against interpolation for these three families.

## Files

| Path | Content |
|---|---|
| [build.sh](build.sh) | Builds the C++ library |
| [src/native.cpp](src/native.cpp) | C++/M4RIE: binary Goppa trials, minor matrices, rank, Sidelnikov–Shestakov, multipliers |
| [src/native.py](src/native.py) | Python interface to the C++ library |
| [src/families.py](src/families.py) | Parameter sets; GRS, alternant and wild Goppa codes in Sage |
| [src/minor.py](src/minor.py) | Sage reference: tangent planes, minor code, certificates, support recovery |
| [src/engine.py](src/engine.py) | One trial; records and fixtures |
| [src/gf2m.py](src/gf2m.py) | NumPy field arithmetic for the independent binary checks |
| [src/run.py](src/run.py) | Parallel, resumable runner |
| [src/summarize.py](src/summarize.py), [src/verify_results.py](src/verify_results.py) | Summary tables; coverage, seed and fixture checks |
| [src/validate.py](src/validate.py) | Validation against Sage |
| [results/](results/) | Summaries, manifest, raw records, verification and validation output |
