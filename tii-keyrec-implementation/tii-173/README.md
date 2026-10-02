# TII-173 implementation

This is the implementation used to recover all 96 support positions and a
degree-eight Goppa polynomial. The holdout sample contains 256 degree-five
forms in 35 variables. Finishing uses an order-296 normalized jet, direct
extraction on the 91 retained columns, and restoration of five removed columns.
The field is GF(128), defined by `z^7+z+1` (`0x83`).

Run the commands below from this directory. Recovery has not been rerun while
preparing this distribution. Choose fresh output directories.

## Verify the saved key

Python 3 alone is sufficient; this does not run recovery:

```sh
python3 recovery/src/independent_verify_key.py \
  ../../tii-results/tii-173/recovered-sk.json \
  --public ../../tii-results/tii-173/pk_McEliece_173.txt \
  --report tii173-key-verification.json
```

## Recover from the saved holdout sample

Requirements: SageMath 10.7, NumPy, M4RI, GNU Make, and a C++17 compiler.
The native Makefile uses Sage's include/library paths; `SAGE` and
`SAGE_LOCAL` may be set explicitly for a nonstandard installation.

```sh
OUT="$PWD/runs/recovery"
python3 ../../tii-kernels/materialize.py ../../tii-kernels/tii-173 --output "$OUT"
make -C recovery/src libseries.so libincremental.so
export DOT_SAGE="$OUT/.sage"

sage -python recovery/src/extend_normalized.py "$OUT" --position 0 --depth 296
sage -python recovery/src/direct_jet_recovery.py "$OUT" --position 0 \
  --output "$OUT/direct-recovery.json"
sage -python recovery/src/deshorten.py "$OUT" \
  --public-key ../../tii-results/tii-173/pk_McEliece_173.txt
python3 recovery/src/independent_verify_key.py "$OUT/equivalent-key.json" \
  --public ../../tii-results/tii-173/pk_McEliece_173.txt \
  --report "$OUT/independent-full-key-verification.json"
sage -python recovery/src/verify_chain.py "$OUT"
```

The chain audit reconstructs the curve, checks every holdout composition
through degree 370, verifies the seven branches and normalized continuation,
and checks coefficient-level jet coherence. The new recovered key is written
to `runs/recovery/equivalent-key.json`.

## Compute a new holdout sample

Requirements: CUDA development toolkit, C++17 compiler, and two visible GPUs.
The saved run used two H200 NVL GPUs, CUDA 13.3, and the exact command parameters
below. The matrix occupies about 13.93 GiB in total; checkpointing requires
about 60 GB of free disk. The two-GPU elimination is part of the numerical
implementation and runs on one host.

```sh
make -C holdout cpu cuda CUDA_ARCH=90
KERNEL_OUT="$PWD/runs/holdout"
holdout/build/holdout-cuda solve ../../tii-kernels/tii-173/operator.txt \
  "$KERNEL_OUT" --count 256 --seed 17320260930 --devices 0,1 \
  --panel-bits 8 --max-gib 64 --checkpoint-every 2048
holdout/build/holdout-cpu verify ../../tii-kernels/tii-173/operator.txt \
  "$KERNEL_OUT/kernel-colex.bin" "$KERNEL_OUT/verification.json" --threads 16
holdout/build/holdout-cpu export "$KERNEL_OUT/kernel-colex.bin" "$KERNEL_OUT/kernel.bin"
cp ../../tii-kernels/tii-173/public.json "$KERNEL_OUT/public.json"
```

Require the verification report to accept 256 independent forms and every
original Hasse order before using this new sample. Set `OUT` to this new
directory and start the previous section at the native-library build.
`kernel-colex.bin` is THK1; export explicitly permutes it to WHK1.

The already prepared public input is included. To regenerate its shortening
and constraints from the original public key, without a kernel solve:

```sh
python3 holdout/scripts/prepare.py --output runs/prepared
```

This source is the original preparation implementation with relocated paths.
The [kernel format](../../tii-kernels/README.md),
[saved results](../../tii-results/tii-173/README.md), and
[source provenance](../PROVENANCE.json) specify the exact stored inputs.

## Recorded minor finishing

Both direct and minor finishings were tested for TII-173. The archived minor
run reused the exact public matrix and coherent jet from the direct run.
At bootstrap depth **58**, it recovered **489 quadrics** and a full
GRS **[91,66]** minor code, then recovered all 91 retained support
coordinates. Exact shortened public-code equality passed. Its internal
finishing time was **2.412 s**, excluding Sage startup and the shared jet.

The archived support comparison passes up to a Möbius transformation, with
Frobenius power 0. Full-key extension was not rerun separately. The evidence
is in [minor-finishing results](../../tii-results/tii-173/evidence/minor-finishing/).
The original bootstrap algorithm is in `recovery/minor/bootstrap.py`; the
bundled engine contains the unchanged mathematical helpers used by that run.

To check only the already saved support comparison, using SageMath:

```sh
sage -python recovery/minor/compare_finishings.py \
  --public ../../tii-kernels/tii-173/public.json \
  --minor ../../tii-results/tii-173/evidence/minor-finishing/bootstrap-public-data.json \
  --direct-key ../../tii-results/tii-173/recovered-sk.json \
  --report tii173-finishing-comparison.json
```

Use a fresh report path. The expected entry is `tii-173` with
`moebius_equivalent: true` and `frobenius_powers_matching: [0]`.
The saved full key is used only by this comparison, after recovery.

For a future reproduction of this finishing alone, starting from the saved
public jet, run in a fresh directory:

```sh
MINOR_OUT="$PWD/runs/minor-finishing"
mkdir -p runs
mkdir "$MINOR_OUT"
cp ../../tii-kernels/tii-173/public.json "$MINOR_OUT/public.json"
cp ../../tii-results/tii-173/evidence/deep-jets-0.json "$MINOR_OUT/deep-jets-0.json"
export DOT_SAGE="$MINOR_OUT/.sage"
sage -python recovery/minor/bootstrap.py "$MINOR_OUT" --position 0 --depth 58
```

This reads only `public.json` and `deep-jets-0.json` and writes
`bootstrap.json` and `bootstrap-public-data.json`. Packaging used the
archived outputs and did not rerun this finishing.

## Optional audits of saved results

The short saved-key check above is sufficient to validate the recovered
Goppa key. The commands below check additional paper claims about the
holdout sample and intermediate algebraic results. They use saved artifacts
and write only to a fresh temporary directory; no recovery pipeline is run.

Run from this `tii-173/` implementation directory. Requirements are Python
3.10+, a C++17 compiler, and SageMath 10.7 with NumPy for the deeper audits, plus M4RI for the native substitution helper.
Use normal Python/Sage mode without `-O`; these audits use assertions.

```sh
TII_REPO="$(cd ../.. && pwd)"
TII_CHECKS="$(mktemp -d "${TMPDIR:-/tmp}/tii-verification.XXXXXX")"
export PYTHONDONTWRITEBYTECODE=1
unset PYTHONOPTIMIZE
export DOT_SAGE="$TII_CHECKS/.sage"
mkdir "$TII_CHECKS/build"

(cd "$TII_REPO/tii-results" && shasum -a 256 -c SHA256SUMS)
(cd "$TII_REPO/tii-kernels" && shasum -a 256 -c SHA256SUMS)
(cd "$TII_REPO/tii-keyrec-implementation" && shasum -a 256 -c SHA256SUMS)
```


### Saved kernel

```sh
c++ -O3 -std=c++17 -pthread \
  "$TII_REPO/tii-keyrec-implementation/tii-173/holdout/src/cpu.cpp" \
  -o "$TII_CHECKS/build/holdout-cpu"
python3 "$TII_REPO/tii-kernels/materialize.py" "$TII_REPO/tii-kernels/tii-173" \
  --output "$TII_CHECKS/tii173"
"$TII_CHECKS/build/holdout-cpu" verify "$TII_REPO/tii-kernels/tii-173/operator.txt" \
  "$TII_CHECKS/tii173/kernel.thk1" "$TII_CHECKS/tii173-hasse.json" --threads 4
```

Expect `all_original_orders_verified: true` and `independent_rank: 256`.
The thread count is for this audit, not the original timing configuration.

Check that the archived THK1 and WHK1 encodings contain exactly the same forms:

```sh
"$TII_CHECKS/build/holdout-cpu" export "$TII_CHECKS/tii173/kernel.thk1" \
  "$TII_CHECKS/tii173-export.whk1"
cmp "$TII_CHECKS/tii173-export.whk1" "$TII_CHECKS/tii173/kernel.bin"
```

`cmp` should exit successfully without output.

### Saved jet and curve

The audit expects the key to be named `equivalent-key.json` beside
`public.json`, `kernel.bin`, and the saved jet/direct-extraction outputs.
Assemble those names in the already restored temporary directory. The native
library is built beside temporary copies of the original source because
this implementation loads the library from its source directory.

```sh
cp "$TII_REPO/tii-results/tii-173/recovered-sk.json" \
  "$TII_CHECKS/tii173/equivalent-key.json"
cp "$TII_REPO/tii-results/tii-173/evidence/deep-jets-0.json" \
  "$TII_REPO/tii-results/tii-173/evidence/direct-recovery.json" \
  "$TII_CHECKS/tii173/"

mkdir "$TII_CHECKS/tii173-code"
cp "$TII_REPO/tii-keyrec-implementation/tii-173/recovery/src/"*.py \
  "$TII_REPO/tii-keyrec-implementation/tii-173/recovery/src/native.cpp" \
  "$TII_REPO/tii-keyrec-implementation/tii-173/recovery/src/Makefile" \
  "$TII_CHECKS/tii173-code/"
make -C "$TII_CHECKS/tii173-code" libseries.so
sage -python "$TII_CHECKS/tii173-code/verify_chain.py" "$TII_CHECKS/tii173"
```

Sage's M4RI headers and library must be available to the supplied Makefile;
`SAGE` and `SAGE_LOCAL` can be passed explicitly to `make` for a nonstandard
installation. This target builds only the substitution library.

The new `tii173/chain-verification.json` should certify all 256 holdout
polynomials as exact curve identities through composition degree 370,
seven Frobenius branches, all continuation ranks equal to 33, and all
10,395 coordinate-jet coefficients coherent. It also checks support-only
extraction from the saved prefix through order 148. It does not compute a
new holdout sample, extend the jet, or restore a new full key.
