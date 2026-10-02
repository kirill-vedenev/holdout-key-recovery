# TII-253 implementation

The saved 280 degree-six forms in 46 variables yield an order-82 coherent
jet. Quadratic bootstrapping and the full `[118,90]` minor GRS code recover
the 118 retained coordinates; IKZ extension restores all 214 original
coordinates and a degree-nine Goppa polynomial. The field is GF(256), modulus
`z^8+z^4+z^3+z^2+1` (`0x11d`).

Both finishings are now recorded. On 2 October 2026, direct finishing
continued the same public order-82 jet to order 198, recovered the 118
retained positions by quadratic sections, and restored a full 214-position
degree-nine Goppa key. The original minor-route key is preserved.

Run these commands from this directory, using fresh output directories.
The `recovery/direct/` sources are the ones executed for the new direct run.

## Verify the saved key

Python 3 alone suffices; this does not run recovery:

```sh
python3 recovery/src/independent_verify.py ../../tii-results/tii-253/recovered-sk.json \
  --public ../../tii-results/tii-253/pk_McEliece_253.txt \
  --report tii253-key-verification.json
```

To also check all 13,344 saved IKZ candidate tests, add:

```sh
--shortened-public ../../tii-kernels/tii-253/public.json \
--extension-evidence ../../tii-results/tii-253/evidence/full-support/positions
```

## Direct finishing

The new [direct key](../../tii-results/tii-253/direct-recovered-sk.json) passes
the same standalone check against the original public matrix:

```sh
python3 recovery/direct/independent_verify.py \
  ../../tii-results/tii-253/direct-recovered-sk.json \
  --public ../../tii-results/tii-253/pk_McEliece_253.txt \
  --report tii253-direct-key-verification.json
```

The [recorded result](../../tii-results/tii-253/evidence/direct-finishing/result-summary.json)
contains the parameters and timings. The direct-support system has 199 rows
and 1,081 quadratic unknowns, rank 190. Direct sections locate the public
infinity parameter at 7; the affine chart is `alpha=1/(beta+7)`. The existing
IKZ calculation then restores all 96 removed coordinates and extracts the
degree-nine Goppa polynomial. No recovered minor-route support or key is an
input to these stages.

To reproduce from the saved public jet prefix, use SageMath 10.7 with NumPy,
a C++17 compiler, and a fresh output directory. The native series cache needs
several GiB of RAM. No GPU or new holdout computation is needed.

```sh
DIRECT_OUT="$PWD/runs/direct-finishing"
sage -python ../../tii-kernels/materialize.py ../../tii-kernels/tii-253 \
  --output "$DIRECT_OUT" --panel
mkdir "$DIRECT_OUT/build"
cp ../../tii-results/tii-253/evidence/direct-finishing/local-branches.json \
  ../../tii-results/tii-253/evidence/direct-finishing/hasse-1.u64le \
  ../../tii-results/tii-253/evidence/direct-finishing/hasse-2.u64le \
  "$DIRECT_OUT/"
c++ -O3 -std=c++17 -shared -fPIC -pthread recovery/direct/series.cpp \
  -o "$DIRECT_OUT/build/libseries.so"
export DOT_SAGE="$DIRECT_OUT/.sage"
export TII_JET_SERIES_LIBRARY="$DIRECT_OUT/build/libseries.so"

sage -python recovery/direct/extend_jet.py "$DIRECT_OUT/public.json" \
  "$DIRECT_OUT/polynomials/kernel-panel.u64le" "$DIRECT_OUT" \
  --depth 198 --threads 4 --prefix ../../tii-results/tii-253/evidence/deep-jets-0.json
sage -python recovery/direct/extract_support.py "$DIRECT_OUT/public.json" "$DIRECT_OUT"
sage -python recovery/direct/recover_full.py --public "$DIRECT_OUT/public.json" \
  --retained-support "$DIRECT_OUT/direct-support.json" \
  --original-public-key ../../tii-results/tii-253/pk_McEliece_253.txt \
  --output "$DIRECT_OUT/full-key"
python3 recovery/direct/independent_verify.py "$DIRECT_OUT/full-key/equivalent-key.json" \
  --public ../../tii-results/tii-253/pk_McEliece_253.txt \
  --report "$DIRECT_OUT/full-key/verification.json"
```

The continuation verifies and preserves the order-82 prefix, then computes
orders 83–198. The final key is `full-key/equivalent-key.json`. The direct
and minor keys use different equivalent support coordinates; each verifies
against the original public key.

## Minor finishing from the saved holdout sample

Requirements: SageMath 10.7 with NumPy and a C++17 compiler. No GPU is needed
after acquiring the kernel. Allow space for the 328 MB THK1 file, 599 MB panel,
75 MB monomial map, local Hasse tensors, and later results.

First restore the exact saved coefficients and recompute the first-jet data:

```sh
OUT="$PWD/runs/recovery"
sage -python ../../tii-kernels/materialize.py ../../tii-kernels/tii-253 \
  --output "$OUT" --panel
mkdir "$OUT/results" "$OUT/build"
cp "$OUT/public.json" "$OUT/results/public.json"
c++ -O3 -std=c++17 recovery/first-jet/extract_jets.cpp -o "$OUT/build/extract_jets"
c++ -O3 -std=c++17 -shared -fPIC -pthread recovery/src/series.cpp \
  -o "$OUT/build/libseries.so"
export DOT_SAGE="$OUT/.sage"
export TII_JET_SERIES_LIBRARY="$OUT/build/libseries.so"

sage -python recovery/first-jet/verify_input_encoding.py "$OUT/kernel.thk1" \
  "$OUT/polynomials/kernel-panel.u64le" "$OUT/results/input-encoding-verification.json"
"$OUT/build/extract_jets" "$OUT/polynomials/kernel-panel.u64le" \
  "$OUT/polynomials/monomials-colex.u64" 46 6 42945927196204 280 "$OUT/results"
sage -python recovery/first-jet/check_criterion.py "$OUT/public.json" \
  "$OUT/polynomials/index.json" "$OUT/results"
```

The historical first-jet stage calls its branch record `branch_solver`;
the continuation stage reads the same record as `branches`. Apply the same
interface conversion used in the completed run:

```sh
python3 - "$OUT/results" <<'PY'
import json, sys
from pathlib import Path
out = Path(sys.argv[1])
record = json.loads((out / 'criterion.json').read_text())
record['branches'] = record['branch_solver']
(out / 'local-branches.json').write_text(json.dumps(record, indent=2) + '\n')
PY
```

Continue the jet, bootstrap and audit the curve, then restore all coordinates:

```sh
sage -python recovery/src/validate_native.py "$OUT/results"
sage -python recovery/src/extend_jet.py "$OUT/public.json" \
  "$OUT/polynomials/kernel-panel.u64le" "$OUT/results" --depth 82 --threads 4
sage -python recovery/src/bootstrap.py "$OUT/results" --position 0
python3 recovery/src/verify_minor.py "$OUT/results"
sage -python recovery/src/prepare_retained.py "$OUT/results"
sage -python recovery/src/verify_jet_substitution.py "$OUT/results" \
  "$OUT/polynomials/kernel-panel.u64le"
sage -python recovery/src/audit_coherence.py "$OUT/public.json" "$OUT/results"
sage -python recovery/src/audit_geometry.py "$OUT/results"
sage -python recovery/src/recover_full.py --public "$OUT/public.json" \
  --retained-support "$OUT/results/retained-support.json" \
  --original-public-key ../../tii-results/tii-253/pk_McEliece_253.txt \
  --output "$OUT/full-key"
python3 recovery/src/independent_verify.py "$OUT/full-key/equivalent-key.json" \
  --public ../../tii-results/tii-253/pk_McEliece_253.txt \
  --shortened-public "$OUT/public.json" --extension-evidence "$OUT/full-key/positions" \
  --report "$OUT/full-key/independent-verification.json"
```

The audits check every jet substitution coefficient, coherence, the recovered
quadratic ideal and minor code, and exact equality with the original public
code. The coherent jet and saved all-order Hasse constraints certify the
holdout identities by `590+83=673 > 594` on the reconstructed curve.

## Compute a new holdout sample

Build/setup and GPU-sequence → CPU-generator → GPU-reconstruction commands
are in [holdout/README.md](holdout/README.md). Once its CPU verifier accepts the
new 280-form sample, set `OUT` to that reconstruction directory, copy the
included `public.json` there, and start at `mkdir "$OUT/results" "$OUT/build"`.
Skip the materializer so the new sample is used.

To regenerate the public shortening and requests from the included original
public key, without solving the kernel:

```sh
python3 holdout/scripts/prepare.py --output runs/prepared
```

See the [storage format](../../tii-kernels/README.md),
[saved results](../../tii-results/tii-253/README.md), and
[source provenance](../PROVENANCE.json).

## Optional audits of saved results

The short saved-key check above is sufficient to validate the recovered
Goppa key. The commands below check additional paper claims about the
holdout sample and intermediate algebraic results. They use saved artifacts
and write only to a fresh temporary directory; no recovery pipeline is run.

Run from this `tii-253/` implementation directory. Requirements are Python
3.10+, a C++17 compiler, and SageMath 10.7 with NumPy for the deeper audits.
Allow about 2 GB of temporary storage for the restored kernel and panel.
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
  "$TII_REPO/tii-keyrec-implementation/tii-253/holdout/src/cpu.cpp" \
  -o "$TII_CHECKS/build/holdout-cpu"
python3 "$TII_REPO/tii-kernels/materialize.py" "$TII_REPO/tii-kernels/tii-253" \
  --output "$TII_CHECKS/tii253"
"$TII_CHECKS/build/holdout-cpu" verify "$TII_REPO/tii-kernels/tii-253/operator.txt" \
  "$TII_CHECKS/tii253/kernel.thk1" "$TII_CHECKS/tii253-hasse.json" --threads 4
```

Expect `all_original_orders_verified: true` and `independent_rank: 280`.
The thread count is for this audit, not the original timing configuration.

### Saved minor code, jet, and geometry

First assemble the saved evidence with `public.json`. That public file is
essential: the verifiers expect it in the same directory as the evidence.

```sh
cp "$TII_REPO/tii-results/tii-253/evidence/bootstrap-public-data.json" \
  "$TII_REPO/tii-results/tii-253/evidence/minor-matrix.u8" \
  "$TII_REPO/tii-results/tii-253/evidence/retained-support.json" \
  "$TII_REPO/tii-results/tii-253/evidence/local-branches.json" \
  "$TII_REPO/tii-results/tii-253/evidence/deep-jets-0.json" \
  "$TII_CHECKS/tii253/"

python3 "$TII_REPO/tii-keyrec-implementation/tii-253/recovery/src/verify_minor.py" \
  "$TII_CHECKS/tii253"
```

This Python-only check must compare all **122,130** saved minor entries,
obtain rank **90**, and prove exact GRS `[118,90]` equality. It writes new
`minor-independent-verification.json` and `minor-code.json` in the temporary
directory.

For the remaining audits, restore the polynomial panel in a separate fresh
directory and build the original substitution helper. The materializer
checks the restored panel and monomial map against their original hashes.

```sh
sage -python "$TII_REPO/tii-kernels/materialize.py" "$TII_REPO/tii-kernels/tii-253" \
  --output "$TII_CHECKS/tii253-panel" --panel
c++ -O3 -std=c++17 -shared -fPIC -pthread \
  "$TII_REPO/tii-keyrec-implementation/tii-253/recovery/src/series.cpp" \
  -o "$TII_CHECKS/build/libseries.so"
export TII_JET_SERIES_LIBRARY="$TII_CHECKS/build/libseries.so"

sage -python "$TII_REPO/tii-keyrec-implementation/tii-253/recovery/src/verify_jet_substitution.py" \
  "$TII_CHECKS/tii253" "$TII_CHECKS/tii253-panel/polynomials/kernel-panel.u64le"
sage -python "$TII_REPO/tii-keyrec-implementation/tii-253/recovery/src/audit_coherence.py" \
  "$TII_CHECKS/tii253/public.json" "$TII_CHECKS/tii253"
sage -python "$TII_REPO/tii-keyrec-implementation/tii-253/recovery/src/audit_geometry.py" \
  "$TII_CHECKS/tii253"
```

Run these in the listed order **after the successful Hasse check above**.
The geometry certificate combines those public constraints with the newly
verified jet substitution and coherence; do not substitute archived Boolean
report flags for these checks.

Expected results are **23,240** zero jet-substitution equalities, **3,818**
coherent coordinate coefficients, **891** global quadrics spanning the full
quadratic ideal, aligned tangent planes at all **118** positions, and the
holdout identity certificate **590 + 83 = 673 > 594**. New reports and the
updated infinity certificate in `retained-support.json` remain in the
temporary copy. The final full-key and all 13,344 extension-candidate checks
are covered by the saved-key check above.
