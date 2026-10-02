# Toy binary Goppa key recovery

C++ implementation for unshortened binary Goppa parameters
`m=6, t=6, n=64, k=28`, over GF(64) with modulus `x^6+x+1`.
The ambient degree bound is `D=n-2t-1=51`.

The program computes a public holdout kernel, separates the six Frobenius
branches, continues a binary-normalized jet, and recovers an equivalent full key.

| Finishing method | Jet prefix | Output |
|---|---:|---|
| Quadratic bootstrap, aligned tangents, minor code, Sidelnikov–Shestakov | 39 | Full GRS `[64,45]` minor code and complete key |
| Direct extraction through a square-root quotient | 101 | Complete key from quadratic sections |

Each method returns 64 distinct finite support elements, a compatible GRS
multiplier, and a monic irreducible degree-six Goppa polynomial. The recovered
alternant and Goppa codes are checked for exact equality with the public binary
code. No shortening or support extension is needed.

## Build

Requires a C++17 compiler and M4RI. Extension-field computations use small
C++ table-based routines. The recovery executable does not require Sage or Python.

If M4RI is available through `pkg-config`:

```sh
make
make check
```

Otherwise supply the installation prefix containing `include/m4ri` and the library:

```sh
make M4RI_PREFIX=/path/to/m4ri
make M4RI_PREFIX=/path/to/m4ri check
```

Compiler, dependency, and source fingerprints are included in
[results/build-info.json](results/build-info.json).

The packed holdout matrix occupies about 1.32 GiB, with additional storage needed
by the solver. `TOY_MATRIX_GIB` sets the matrix-allocation cap, which defaults to
3 GiB; it does not bound total process memory.

## Run

Use a fresh output directory:

```sh
./run.sh instances/a/public.txt runs/demo both
```

The final argument is `minor`, `direct`, or `both`. This command computes a new
holdout kernel and verifies the recovered key. Minor-only mode computes a
normalized jet through order 40 and uses its order-39 prefix. Direct mode computes
through order 102 and uses only the order-101 prefix for its sections.

To run recovery from a saved kernel:

```sh
./build/toy-goppa recover instances/a/public.txt runs/demo/kernel.bin runs/replay both
```

Kernel independence and every retained Hasse condition are checked before the
polynomials are used.

## Input instances

[instances/a](instances/a/) and [instances/b](instances/b/) contain fixed test
inputs with the same parameters. Each has a public generator in `public.txt` and
an audit key in `audit/secret.txt`. Recovery reads only the public file. The
included seed and generator settings for instance B are in its `generation.json`.

Generate another key with an integer seed:

```sh
./build/toy-goppa generate 1 runs/new-key
./run.sh runs/new-key/public.txt runs/new-key-recovery both
```

Generation uses `std::mt19937_64`, an explicit support shuffle, and a degree-six
irreducibility test. It checks the public dimension and the binary square-free
Goppa identity. Failed hypotheses stop the run.

## Verify a saved key

The C++ verifier uses the public input alone:

```sh
./build/toy-goppa verify instances/a/public.txt results/a/recovery/direct-key.txt
```

The independent verifier uses Python's standard library and imports none of the
recovery implementation:

```sh
python3 verify.py instances/a/public.txt results/a/recovery/direct-key.json
python3 verify.py instances/a/public.txt results/a/recovery/minor-key.json
```

An optional audit compares the recovered key with the generated secret:

```sh
./build/toy-goppa audit instances/a/public.txt instances/a/audit/secret.txt results/a/recovery/direct-key.txt
```

It checks that one global Frobenius power and one Möbius transformation relate
the supports. Secret data enter only generation and this separate audit.

## Exact checks

1. **Holdout:** 98,280 squarefree degree-five monomials and multiplicity four
   at retained points. Orders 1 and 3 suffice for construction. All orders 0
   through 3 are verified, with unit-column conditions automatic.
2. **Branches:** gradient rank 21, nullity 7, and quadratic restriction rank 15.
   All six branches have continuation rank 26.
3. **Jets:** even coefficients are squares of their half-order coefficients.
   Odd coefficients satisfy the regular and next-even-order equations. A fresh
   full substitution verifies the completed series against every kernel polynomial.
4. **Minor finishing:** 406 quadratic unknowns, constraint rank 97, and 309
   quadrics; tangent rank 26 at all 64 positions; 378 minor rows spanning the
   full GRS `[64,45]` code.
5. **Direct finishing:** quadratic sections with contact orders 102 and 100;
   square roots of their public evaluation ratios recover the support.
6. **Certificates:** exact public-code equality, a reconstructed curve of degree
   at most 51, coefficient-level jet coherence, all branch planes, and global
   tangent alignment. The two recovered supports agree up to a Möbius map.

The initial ordinary root count is 252, below the degree bound 255. Verified
coherence through order 102 supplies 103 additional roots at the held-out point,
so `252+103=355>255`. Together with the Hasse checks and full substitution, this
certifies the returned curve identities using the recovered public representation.

## Results

| Instance | Kernel dimension | Holdout and verification | Recovery and verification | Outcome |
|---|---:|---:|---:|---|
| A | 552 | 186.6 s | 13.9 s | Both methods pass |
| B | 414 | 194.8 s | 13.1 s | Both methods pass |

Recovery timings include the kernel recheck, branch recovery, continuation,
both finishing methods, and public geometric certificates. Optional secret audits
are excluded. These are individual toy runs, not statistical success estimates.

[results/](results/) contains complete recovered keys, intermediate matrices,
and verification records. Minor-only and direct-only modes on instance A return
the same keys as combined mode. Mutated kernel and multiplier inputs are rejected.

See [FORMATS.md](FORMATS.md) for all encodings.
