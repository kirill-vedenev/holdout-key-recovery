# TII-249 implementation

The saved 280 degree-five forms in 50 variables yield an order-290 coherent
jet, 178 retained support coordinates by direct extraction, and then all
235 original coordinates by IKZ extension. The final key has a monic
irreducible degree-16 Goppa polynomial over GF(256), modulus
`z^8+z^4+z^3+z^2+1` (`0x11d`).

Run these commands from this directory, using fresh output directories.
This layout has not been used for a new recovery run.

## Verify the saved key

This standard-library-only command checks the full key without recovery:

```sh
python3 recovery/full-key/independent_verify.py \
  ../../tii-results/tii-249/recovered-sk.json \
  --public ../../tii-results/tii-249/pk_McEliece_249.txt \
  --report tii249-key-verification.json
```

To also repeat the independent check of all 4,503 saved IKZ candidate tests,
add these arguments to that command:

```sh
--shortened-public ../../tii-kernels/tii-249/public.json \
--extension-evidence ../../tii-results/tii-249/evidence/full-key/positions
```

## Recover from the saved holdout sample

Requirements: SageMath 10.7 with NumPy and a C++17 compiler. This part needs
no GPU. The materializer restores the original THK1 bytes and the exact panel
and monomial map; it verifies their recorded hashes.

```sh
OUT="$PWD/runs/recovery"
sage -python ../../tii-kernels/materialize.py ../../tii-kernels/tii-249 \
  --output "$OUT" --panel
mkdir "$OUT/results" "$OUT/build"
c++ -O3 -std=c++17 recovery/jet-support/extract_jets.cpp -o "$OUT/build/extract_jets"
c++ -O3 -std=c++17 -shared -fPIC -pthread recovery/jet-support/series.cpp \
  -o "$OUT/build/libseries.so"
export DOT_SAGE="$OUT/.sage"
export TII_JET_SERIES_LIBRARY="$OUT/build/libseries.so"

sage -python recovery/jet-support/verify_input_encoding.py "$OUT/kernel.thk1" \
  "$OUT/polynomials/kernel-panel.u64le" "$OUT/results/input-encoding-verification.json"
"$OUT/build/extract_jets" "$OUT/polynomials/kernel-panel.u64le" \
  "$OUT/polynomials/monomials-colex.u64" 50 5 505067484167246 280 "$OUT/results"
sage -python recovery/jet-support/validate_native.py "$OUT/results"
sage -python recovery/jet-support/analyze_local.py "$OUT/public.json" "$OUT/results"
sage -python recovery/jet-support/extend_jet.py "$OUT/public.json" \
  "$OUT/polynomials/kernel-panel.u64le" "$OUT/results" --depth 290 --threads 4
sage -python recovery/jet-support/extract_support.py "$OUT/public.json" "$OUT/results"
python3 recovery/jet-support/independent_verify_support.py "$OUT/public.json" \
  "$OUT/results/direct-support.json" \
  --original-public-key ../../tii-results/tii-249/pk_McEliece_249.txt \
  --report "$OUT/results/independent-verification.json"
sage -python recovery/jet-support/audit_coherence.py "$OUT/public.json" "$OUT/results"

sage -python recovery/full-key/recover_full.py --public "$OUT/public.json" \
  --retained-support "$OUT/results/direct-support.json" \
  --original-public-key ../../tii-results/tii-249/pk_McEliece_249.txt \
  --output "$OUT/full-key"
python3 recovery/full-key/independent_verify.py "$OUT/full-key/equivalent-key.json" \
  --public ../../tii-results/tii-249/pk_McEliece_249.txt \
  --shortened-public "$OUT/public.json" --extension-evidence "$OUT/full-key/positions" \
  --report "$OUT/full-key/independent-verification.json"
```

The outputs are the recovered support, multipliers, Goppa polynomial,
candidate uniqueness evidence, and public-code equality checks. All inputs
to these computations are public data or newly computed intermediate values.

## Compute a new holdout sample

The original implementation has three stages. Build/setup and direct commands
are in [holdout/README.md](holdout/README.md). Once the final CPU verifier accepts
the new sample, its `kernel.thk1` and `polynomials/` already have the layout
used above. Copy the included `public.json` into that output, set `OUT` to its
path, and start at `mkdir "$OUT/results" "$OUT/build"`, skipping materialization
of the saved sample. The documented jet extraction expects 280 forms.

The deterministic public preparation source is also retained:

```sh
python3 holdout/scripts/prepare.py --output runs/prepared
```

See the [storage format](../../tii-kernels/README.md),
[saved results](../../tii-results/tii-249/README.md), and
[source provenance](../PROVENANCE.json).

## Recorded minor finishing

Both direct and minor finishings were tested for TII-249. The archived minor
run reused the exact public matrix and coherent jet from the direct run.
At bootstrap depth **113**, it recovered **1000 quadrics** and a full
GRS **[178,129]** minor code, then recovered all 178 retained support
coordinates. Exact shortened public-code equality passed. Its internal
finishing time was **7.087 s**, excluding Sage startup and the shared jet.

The archived support comparison passes up to a Möbius transformation, with
Frobenius power 0. Full-key extension was not rerun separately. The evidence
is in [minor-finishing results](../../tii-results/tii-249/evidence/minor-finishing/).
The original bootstrap algorithm is in `recovery/minor/bootstrap.py`; the
bundled engine contains the unchanged mathematical helpers used by that run.

To check only the already saved support comparison, using SageMath:

```sh
sage -python recovery/minor/compare_finishings.py \
  --public ../../tii-kernels/tii-249/public.json \
  --minor ../../tii-results/tii-249/evidence/minor-finishing/bootstrap-public-data.json \
  --direct-key ../../tii-results/tii-249/recovered-sk.json \
  --report tii249-finishing-comparison.json
```

Use a fresh report path. The expected entry is `tii-249` with
`moebius_equivalent: true` and `frobenius_powers_matching: [0]`.
The saved full key is used only by this comparison, after recovery.

For a future reproduction of this finishing alone, starting from the saved
public jet, run in a fresh directory:

```sh
MINOR_OUT="$PWD/runs/minor-finishing"
mkdir -p runs
mkdir "$MINOR_OUT"
cp ../../tii-kernels/tii-249/public.json "$MINOR_OUT/public.json"
cp ../../tii-results/tii-249/evidence/jet-support/deep-jets-0.json "$MINOR_OUT/deep-jets-0.json"
export DOT_SAGE="$MINOR_OUT/.sage"
sage -python recovery/minor/bootstrap.py "$MINOR_OUT" --position 0 --depth 113
```

This reads only `public.json` and `deep-jets-0.json` and writes
`bootstrap.json` and `bootstrap-public-data.json`. Packaging used the
archived outputs and did not rerun this finishing.

## Optional audits of saved results

The short saved-key check above is sufficient to validate the recovered
Goppa key. The commands below check additional paper claims about the
holdout sample and intermediate algebraic results. They use saved artifacts
and write only to a fresh temporary directory; no recovery pipeline is run.

Run from this `tii-249/` implementation directory. Requirements are Python
3.10+, a C++17 compiler, and SageMath 10.7 with NumPy for the deeper audits.
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
  "$TII_REPO/tii-keyrec-implementation/tii-249/holdout/src/cpu.cpp" \
  -o "$TII_CHECKS/build/holdout-cpu"
python3 "$TII_REPO/tii-kernels/materialize.py" "$TII_REPO/tii-kernels/tii-249" \
  --output "$TII_CHECKS/tii249"
"$TII_CHECKS/build/holdout-cpu" verify "$TII_REPO/tii-kernels/tii-249/operator.txt" \
  "$TII_CHECKS/tii249/kernel.thk1" "$TII_CHECKS/tii249-hasse.json" --threads 4
```

Expect `all_original_orders_verified: true` and `independent_rank: 280`.
The thread count is for this audit, not the original timing configuration.

### Saved support and coherent jet

Use the archived retained support and local branches together with the
restored public generator. The full-key and extension-candidate checks were
covered by the saved-key verification command above.

```sh
cp "$TII_REPO/tii-results/tii-249/evidence/jet-support/deep-jets-0.json" \
  "$TII_REPO/tii-results/tii-249/evidence/jet-support/direct-support.json" \
  "$TII_REPO/tii-results/tii-249/evidence/jet-support/local-branches.json" \
  "$TII_CHECKS/tii249/"

python3 "$TII_REPO/tii-keyrec-implementation/tii-249/recovery/jet-support/independent_verify_support.py" \
  "$TII_CHECKS/tii249/public.json" "$TII_CHECKS/tii249/direct-support.json" \
  --original-public-key "$TII_REPO/tii-results/tii-249/pk_McEliece_249.txt" \
  --report "$TII_CHECKS/tii249/support-verification.json"
sage -python "$TII_REPO/tii-keyrec-implementation/tii-249/recovery/jet-support/audit_coherence.py" \
  "$TII_CHECKS/tii249/public.json" "$TII_CHECKS/tii249"
```

Expect exact equality of the retained `[178,50]` public code with the recovered
alternant code, eight matching Frobenius branch planes, and all **14,550**
saved coordinate-jet coefficients fitting one scalar/parameter series.
The new coherence report is `tii249/coherence-audit.json`. This audit does
not claim a full degree-725 coefficient check of all quintic curve
compositions; the exact full-key certificate is the saved-key check above.
