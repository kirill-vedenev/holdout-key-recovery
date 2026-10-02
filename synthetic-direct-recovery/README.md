# Synthetic direct recovery

This collection runs the direct key recovery of the revised paper (Section 7)
from one local jet of the secret curve. The jet is computed from the generated
key instead of a holdout kernel, so the test isolates the recovery step. It
covers binary Goppa codes, wild Goppa codes and generic alternant codes.

## Keys

Every key is a subfield subcode `C = GRS_{D+1}(α, λ) ∩ F_q^n` with support and
multiplier in `K = GF(q^m)`. The public generator matrix `Y` is systematic, and
its columns are `y_j = λ_j F(α_j)` for a curve `F` of degree at most `D`.

| Family | Code | Degree bound `D` | Multiplier `λ_j` |
|---|---|---|---|
| Binary Goppa, `q = 2` | `Γ(α, G)`, `G` irreducible of degree `t` | `n − 2t − 1` | `G(α_j)² / Π'(α_j)` |
| Wild Goppa, `q > 2` | `Γ(α, g^(q−1))`, `g` irreducible of degree `t` | `n − qt − 1` | `g(α_j)^q / Π'(α_j)` |
| Generic alternant | `GRS_{n−r}(α, λ) ∩ F_q^n`, random `λ` | `n − r − 1` | random |

Here `Π` is the support polynomial. For the Goppa families the degree bound
uses `Γ(α, G) = Γ(α, G²)` and `Γ(α, g^(q−1)) = Γ(α, g^q)`. The program checks
these bounds on every key by interpolating the curve through `D + 1` support
points and comparing the remaining values.

Odd characteristic is supported: `K` is any field of order at most `2^16`,
with Zech logarithms for addition.

## The jet

At a random public position `i`, the jet consists of the first `L`
coefficients of

```
u(T) · φ^σ(λ_i) · F^σ(φ^σ(α_i) + ψ(T)),
```

where `φ` is the Frobenius map `x ↦ x^q` and `σ` is a random branch. Both `u`
and `ψ` are random polynomials of degree at most 4 with `u(0) = 1`, `ψ(0) = 0`
and `ψ'(0) ≠ 0`, so the jet is known only up to reparametrization and scaling.
Its constant term is the public column `y_i`.

## Recovery

The recovery uses only `Y`, the position `i` and the jet `V(T)`.

1. **Forms.** For forms `h` of degree `δ`, the coefficients `[T^e] h(V(T))`
   are linear in the coefficients of `h`. The program solves for one `R` with
   `[T^e] R(V) = 0` for `e < δD` and `R(y_j0) = 1` at another public column,
   and for one `S` with `[T^e] S(V) = 0` for `e < δD − 1` and
   `[T^(δD−1)] S(V) = 1`. Then `R(F) = c (Z − α_i)^(δD)` and `S(F)` has a root
   of exact order `δD − 1` at `α_i`. The degree `δ` is the smallest with
   `C(k+δ−1, δ) ≥ δD + 2`. The jet order is `δD − 1`.
2. **Binary Goppa codes.** Adjacent orders are impossible here, so the program
   uses `δ = 2` with orders `2D` and `2D − 2`, and takes square roots of the
   quotients. The jet order is `2D − 1`.
3. **Support.** `β_j = R(y_j) / S(y_j)` for `j ≠ i`, and `β_i = 0`; a zero
   denominator means `∞`. These values are the image of the support under one
   Möbius map. If a value is `∞`, the map `β ↦ 1/(β − p)` with a free point `p`
   makes the support finite.
4. **Multiplier.** The program solves `Y diag(z) H_RS^T = 0`, where `H_RS` is a
   parity-check matrix of `RS_{D+1}(β)`. It requires a one-dimensional
   solution space and a solution without zero entries, and sets `λ'_j = 1/z_j`.
5. **Verification.** It checks `C ⊆ GRS_{D+1}(β, λ')` and that the subfield
   subcode `GRS_{D+1}(β, λ') ∩ F_q^n` has dimension `k`, so it equals `C`.

The products of the monomials with the jet are computed only for a subset of
monomials: all pure powers `X_a^δ` and random others, slightly more than the
number of conditions. The pure powers are needed because the generator matrix
is systematic. If a system is inconsistent, the subset grows up to all
monomials before the program reports that `R` or `S` does not exist.

**Audit.** With the secret key, the program checks that `β` is a Möbius image
of `φ^σ(α)` by comparing cross-ratios. A key counts as recovered only if both
the verification and the audit pass.

## Checks

- The self-test checks field arithmetic in 16 fields of characteristic 2, 3,
  5, 7, 11 and 13, subfields, Taylor shifts, interpolation, and the number of
  irreducible polynomials over GF(4) and GF(9). It also recovers 24 small keys
  of all families and confirms that the audit rejects a permuted support.
- **Control.** For binary Goppa codes, the same procedure with adjacent orders
  must fail, because no `S` exists. The control set runs it.
- **SageMath cross-check.** `crosscheck.sh` recovers small keys of every family
  and `crosscheck.py` rechecks them independently: containment in the recovered
  GRS code, the subfield-subcode dimension through the trace code of the dual,
  and the Möbius audit. Keys with a swapped support value or a changed
  multiplier entry are rejected.

## Results

The keys ran on a 96-core AMD EPYC machine under Linux, one process per key.
The SageMath cross-check ran under macOS on an Apple M1 Pro. Rerunning a key on
either platform reproduces its JSON line exactly, apart from timings.
[results/summary.md](results/summary.md) has the complete tables, and
`results/keys.jsonl.gz` has one JSON line per key.

Every key was recovered, 1,060 in total. The forms `R` and `S` existed for
every key, the recovered key passed the verification, and the audit confirmed
a Möbius image of the support on the branch of the jet. Times are for one core.

| Set | q | m | n | k | D | δ | Jet order | Keys recovered | Median s |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Binary Goppa | 2 | 8 | 100 | 36 | 83 | 2 | 165 | 100 / 100 | 0.01 |
| Binary Goppa | 2 | 10 | 400 | 200 | 359 | 2 | 717 | 20 / 20 | 1.3 |
| Binary Goppa, `mceliece348864` with Weis shortening | 2 | 12 | 922 | 154 | 793 | 2 | 1585 | 20 / 20 | 18 |
| Binary Goppa, `mceliece6688128` with Weis shortening | 2 | 13 | 1888 | 224 | 1631 | 2 | 3261 | 20 / 20 | 183 |
| Wild Goppa | 3 | 5 | 200 | 140 | 181 | 2 | 361 | 100 / 100 | 1.4 |
| Wild Goppa | 3 | 7 | 800 | 520 | 739 | 2 | 1477 | 20 / 20 | 112 |
| Wild Goppa | 4 | 4 | 200 | 128 | 175 | 2 | 349 | 100 / 100 | 0.2 |
| Wild Goppa | 4 | 5 | 600 | 420 | 551 | 2 | 1101 | 20 / 20 | 8.7 |
| Wild Goppa | 5 | 3 | 100 | 52 | 79 | 2 | 157 | 100 / 100 | 0.08 |
| Wild Goppa | 5 | 4 | 500 | 372 | 459 | 2 | 917 | 20 / 20 | 24 |
| Alternant | 2 | 8 | 200 | 72 | 183 | 2 | 365 | 100 / 100 | 0.1 |
| Alternant | 2 | 8 | 200 | 24 | 177 | 3 | 530 | 100 / 100 | 0.2 |
| Alternant | 2 | 11 | 1000 | 560 | 959 | 2 | 1917 | 20 / 20 | 39 |
| Alternant | 3 | 5 | 200 | 140 | 187 | 2 | 373 | 100 / 100 | 1.4 |
| Alternant | 3 | 7 | 800 | 520 | 759 | 2 | 1517 | 20 / 20 | 92 |
| Alternant | 4 | 4 | 200 | 152 | 187 | 2 | 373 | 100 / 100 | 0.2 |
| Alternant | 5 | 3 | 100 | 76 | 91 | 2 | 181 | 100 / 100 | 0.08 |

For binary Goppa codes the forms have orders `2D` and `2D − 2`, and the
support comes from square roots of the quotients. The other families use
adjacent orders. The alternant set with `k = 24` has too few quadrics,
`C(25, 2) < 2D + 2`, so it uses cubic forms.

**Control.** The binary Goppa control with adjacent orders, `(m, t, n) =
(10, 20, 400)`, failed on all 20 keys because no `S` exists, even with all
monomials. This is the obstruction that the square-root variant avoids.

**Cross-check.** SageMath rechecked 36 recovered keys, four of each small
parameter set, and all passed.

Odd-characteristic keys are slower, because addition goes through Zech
logarithms. Most of the time goes into the two linear systems for `R` and `S`,
whose size is about `δD`.

## Build and run

Requires a C++17 compiler and Python 3. No other libraries are needed.
SageMath 10.7 or later is optional and is used only for the cross-check.

```sh
make
make check
```

`make check` runs the self-test. Single keys are recovered as follows:

```sh
./build/synthetic_direct run --family binary-goppa --q 2 --m 12 --n 922 --t 64 --seeds 1-3
./build/synthetic_direct run --family wild-goppa --q 3 --m 5 --n 200 --t 6 --seeds 1-10
./build/synthetic_direct run --family alternant --q 4 --m 4 --n 200 --r 12 --seeds 1-10
```

Each key prints one JSON line. `--adjacent` runs the binary Goppa control,
`--delta` fixes the degree of the forms, and `--dump DIR` writes recovered keys
for the cross-check.

### Reproducing the results

`run.sh` builds the program, runs the self-test and all 1,080 keys, and writes
`results/`. `JOBS` sets the number of parallel processes and defaults to the
number of CPUs. The complete run takes about five CPU-hours, about four minutes
on 96 cores; the slowest single key takes about three minutes.

```sh
./run.sh
```

When SageMath is installed, `run.sh` also runs `crosscheck.sh`. `summarize.py`
turns the JSON lines into [results/summary.md](results/summary.md). Keys are
determined by their seeds, so a rerun reproduces every JSON line apart from
timings. [results/build-info.json](results/build-info.json) records the
compilers and the source hash of the published results.

## Files

| File | Content |
|---|---|
| [src/synthetic_direct.cpp](src/synthetic_direct.cpp) | Fields, keys, jets, direct recovery, verification, audit |
| [run.sh](run.sh) | Complete experiment |
| [crosscheck.sh](crosscheck.sh), [crosscheck.py](crosscheck.py) | Independent SageMath check of recovered keys |
| [summarize.py](summarize.py) | Summary tables |
| [results/](results/) | JSON lines per key, summary, cross-check, build information |
