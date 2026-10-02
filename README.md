# GRS-subcode key recovery experiments

Implementations and experiments for recovering an equivalent support and
multiplier from a public generator matrix. The public recovery pipeline computes
polynomial holdout kernels, separates local branches, and continues coherent
jets. Finishing uses quadratic bootstrapping and minor codes, or direct
extraction from high-contact sections.

This repository is the supplement to the paper

> Kirill Vedenev. *Extending Distinguishing to Key Recovery for Subfield
> Subcodes of GRS codes.* Cryptology ePrint Archive, Paper 2026/1747.
> <https://eprint.iacr.org/2026/1747>

See [Citation](#citation) for a BibTeX entry.

## Experiments

| Collection | Codes | Implementation | Finishing |
|---|---|---|---|
| [Toy key recovery](keyrec-toy/README.md) | GRS subcodes and alternant codes | SageMath notebooks | Minor and direct in characteristic two; direct in odd characteristic |
| [Toy binary Goppa recovery](keyrec-toy-binary-goppa/README.md) | Unshortened binary `[64,28]` Goppa codes, `m=6, t=6` | C++ and M4RI | Minor and binary square-root direct extraction |
| [Minor-code fullness](synthetic-minor-code-fullness/README.md) | GRS subcodes, alternant, wild Goppa and binary Goppa codes; 100 codes per smaller set and 20 per shortened Classic McEliece set (4,600 total) | SageMath, C++ and M4RIE | Minor, from tangent planes computed with the secret key |
| [Wild Goppa recovery over GF(4)](keyrec-toy-wild-goppa/README.md) | Wild Goppa codes over GF(4) with support in GF(64) | C++ and OpenMP-enabled M4RI | One holdout, deep jets, minor and direct finishing |
| [Direct-extraction existence condition](goppa-direct-sufficient-condition/README.md) | Binary Goppa codes: Classic McEliece shortened to the distinguisher dimension, and small codes | C++, SageMath cross-check | Checks `dim W = 2D−t+1` on the secret curve |
| [Synthetic direct recovery](synthetic-direct-recovery/README.md) | Binary Goppa (including shortened Classic McEliece), wild Goppa and generic alternant codes over q = 2, 3, 4, 5 | C++, SageMath cross-check | Direct, from a local jet computed with the secret key |

Each collection contains its implementation, run instructions, and saved results.
The [evidence summary](EVIDENCE.md) lists the checks and outcomes.

## TII challenge artifacts

The [TII experiment specification](TII-EXPERIMENTS.md) lists the exact
shortening sets, systematic positions, holdout and computational parameters,
measured hardware and stage runtimes, and the complete recovery pipeline for
each instance. All coordinate ranges there use explicit zero-based original
or shortened numbering.

Each [implementation README](tii-keyrec-implementation/README.md) gives a short
Python command for checking the saved key. Additional audits of intermediate
results are optional and documented with the corresponding implementation.

Completed recoveries for TII-173, TII-249, and TII-253 are organized into:

- [tii-results](tii-results/README.md): original public keys, recovered
  equivalent secret keys, explicit extension-field polynomials, and evidence.
- [tii-kernels](tii-kernels/README.md): the recovered holdout polynomial
  samples, coordinate maps, checksums, and complete storage specification.
- [tii-keyrec-implementation](tii-keyrec-implementation/README.md): the
  implementations used, grouped by instance, with direct build/run commands.

The saved keys recover all original coordinates. The source packages retain
the computational stages and omit the server orchestration. They were
assembled from existing artifacts. TII-253 also has a direct finishing run
completed on 2 October 2026 from the saved public jet and holdout sample;
both its direct and minor keys verify against the original public key.

## Run

For the notebooks, use a SageMath Jupyter kernel or run:

```sh
cd keyrec-toy
sage -python execute.py
```

For the binary Goppa experiments, install M4RI and a C++17 compiler, then run:

```sh
cd keyrec-toy-binary-goppa
make
make check
./run.sh instances/a/public.txt runs/demo both
```

The C++ collection accepts an explicit M4RI installation prefix when the library
is not available through `pkg-config`; see its build instructions.

For the existence condition of binary direct extraction, a C++17 compiler and
Python 3 suffice:

```sh
cd goppa-direct-sufficient-condition
make
make check
./build/product_rank preset mceliece348864-weis --seeds 1-10 --positions 8
```

The synthetic direct recovery has the same requirements:

```sh
cd synthetic-direct-recovery
make
make check
./build/synthetic_direct run --family wild-goppa --q 3 --m 5 --n 200 --t 6 --seeds 1-10
```

## Inputs and verification

The toy and TII recovery runs read public matrices and public parameters.
Generated secret data, where available, are used only for separate audits of
these runs. Three collections test individual steps using data derived from
generated keys:

- Minor-code fullness computes tangent planes from the secret key, then tests
  the minor code and support recovery using the public matrix and minor matrix.
- The direct-extraction existence condition is checked on products of the
  secret curve polynomials.
- Synthetic direct recovery computes one local jet from the secret key, then
  runs recovery using the public matrix and that jet.

Fixed input instances are included with the binary Goppa C++ collection, and
its generator can create additional instances.

The public recovery runs check kernel membership, branch and continuation
ranks, jet coherence, and the recovered algebraic representation. The toy,
TII, and synthetic direct recoveries for alternant and Goppa codes verify exact
equality with the public code. The minor-code fullness sweeps verify the minor
GRS code and containment of the public code in the recovered ambient GRS code;
for shortened codes these checks cover the retained coordinates. The binary
Goppa C++ collection also provides an independent verifier using Python's
standard library.

Saved results include parameters, ranks, timings, and verification outcomes.
The toy collections are fixed examples, not statistical estimates of success
probability. The TII collection records three specific completed challenge
recoveries; its saved measurements are not general performance claims.
Generated build products and temporary runs are excluded from version control.

## Citation

If you use this code or these results, please cite the paper.
[CITATION.cff](CITATION.cff) contains the same reference in machine-readable
form.

```bibtex
@misc{cryptoeprint:2026/1747,
  author       = {Kirill Vedenev},
  title        = {Extending Distinguishing to Key Recovery for Subfield
                  Subcodes of {GRS} codes},
  howpublished = {Cryptology {ePrint} Archive, Paper 2026/1747},
  year         = {2026},
  url          = {https://eprint.iacr.org/2026/1747}
}
```

## License

The code, documentation and results in this repository are released under the
[MIT License](LICENSE), with these exceptions:

- `tii-keyrec-implementation/*/holdout/third_party/saarinen/` contains code by
  Markku-Juhani O. Saarinen from
  [tii254-artifact](https://github.com/mjosaarinen/tii254-artifact), under his
  MIT License in the `LICENSE.txt` beside it.
- `tii-keyrec-implementation/*/holdout/patches/cado-coefficient-parallel.patch`
  modifies [CADO-NFS](https://gitlab.inria.fr/cado-nfs/cado-nfs) and is
  distributed under the GNU Lesser General Public License, version 2.1, as
  CADO-NFS itself; see `CADO-COPYING` beside it.

Some programs link against M4RI and M4RIE, which are licensed under the GNU
General Public License, version 2 or later, and some scripts run under
SageMath. These libraries are not included. Compiled programs that contain
them must be distributed under the terms of the GPL.
