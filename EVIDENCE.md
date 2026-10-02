# Experimental evidence

The results below are included in this repository. The toy and TII recovery
runs start from public matrices and parameters; generated secret data, where
available, are used only for separate audits of those runs. Three collections
test individual steps using secret-derived data: tangent planes for minor-code
fullness, curve polynomials for the direct-extraction existence condition, and
one local jet for synthetic direct recovery.

## GRS subcodes and alternant codes

The [notebook collection](keyrec-toy/README.md) contains four fixed examples.
Each computes its holdout kernel, recovers a coherent jet, and extracts a support
and multiplier.

| Code | Base field | Ambient field | Jet order | Finishing | Result |
|---|---|---|---:|---|---|
| GRS subcode `[24,6]` | GF(32) | GF(32) | 13 | Minor and direct | Pass |
| Alternant `[60,20]` | GF(8) | GF(64) | 77 | Minor and direct | Pass |
| GRS subcode `[17,4]` | GF(19) | GF(19) | 14 | Direct | Pass |
| Alternant `[45,15]` | GF(11) | GF(121) | 57 | Direct | Pass |

All examples satisfy the ordinary interpolation bound `s(n-1)>dD`. Audits check
the curve identities, recovered branches, coefficient-level jet coherence, and
global Frobenius/Möbius equivalence of the supports. Both characteristic-two
examples verify exact minor-GRS equality and agreement of the finishing routes.
Both alternant examples verify exact equality of the recovered subfield code
with the public code. Corrupted multipliers are rejected.

The [notebook result records](keyrec-toy/results/README.md) contain the
parameters, ranks, timings, and individual checks.

## Binary Goppa codes

The [C++ collection](keyrec-toy-binary-goppa/README.md) contains two unshortened
instances with `(m,t,n,k)=(6,6,64,28)` and ambient degree bound `D=51`.

| Instance | Holdout dimension | Branches | Continuation rank | Minor code | Direct extraction |
|---|---:|---:|---:|---|---|
| [A](keyrec-toy-binary-goppa/results/a/recovery/recovery.json) | 552 | 6 | 26 | Full GRS `[64,45]` | Pass |
| [B](keyrec-toy-binary-goppa/results/b/recovery/recovery.json) | 414 | 6 | 26 | Full GRS `[64,45]` | Pass |

Both finishing routes recover all 64 support positions, a compatible multiplier,
and a monic irreducible degree-six Goppa polynomial. Independent public verifiers
check exact code equality. Separate audits check equivalence with the generated
secret key. Minor-only and direct-only runs on instance A return the same keys
as combined mode.

Minor finishing uses a jet prefix through order 39; binary direct finishing uses
order 101 and takes a square root of the quotient.

The initial ordinary root count, 252, is below the degree bound 255. After
recovery, exact coefficient-level coherence through order 102 and verified Hasse
conditions certify the holdout identities with `252+103=355>255` roots counted
with multiplicity. This certificate uses the reconstructed curve and public data.

The [validation record](keyrec-toy-binary-goppa/results/validation.json) includes
arithmetic checks, independent key verification, and rejection of corrupted
kernel and multiplier inputs. The results establish these fixed instances;
they do not assert that the rank hypotheses hold for every key.

## Minor-code fullness

The [minor-code fullness collection](synthetic-minor-code-fullness/README.md)
tests 54 parameter sets in characteristic two, using tangent planes computed
from generated secret keys. It checks the minor code and subsequent support
recovery; it does not compute the tangent planes from a public holdout kernel.
The 44 smaller sets have 100 codes each, and the 10 shortened Classic McEliece
sets have 20 codes each, for 4,600 trials.

| Family | Sets | Codes | Full minor GRS code and successful recovery |
|---|---:|---:|---:|
| GRS subcodes | 7 | 700 | 700 / 700 |
| Alternant | 11 | 1,100 | 1,099 / 1,100 |
| Wild Goppa | 15 | 1,500 | 1,500 / 1,500 |
| Binary Goppa, smaller sets | 11 | 1,100 | 1,100 / 1,100 |
| Classic McEliece, shortened | 10 | 200 | 200 / 200 |
| Total | 54 | 4,600 | 4,599 / 4,600 |

All 4,599 trials with two-dimensional tangent planes at every position yield
a full minor GRS code. Sidelnikov–Shestakov recovers an equivalent support,
and independent checks verify containment of the public code in the recovered
ambient GRS code. For shortened codes, these checks cover the retained
coordinates.

The remaining trial, `alternant-a3` trial 85, has a one-dimensional tangent
space at position 40 (zero-based). That column of the minor matrix is zero,
and the complete minor-code rank is 50 rather than 51. This tangent failure
is retained in the results; no trial raised an error or timed out. The
[per-set summary](synthetic-minor-code-fullness/results/summary.md) gives all
ranks, outcomes, and timings.

## Existence condition for binary direct extraction

The [existence-condition collection](goppa-direct-sufficient-condition/README.md)
generates binary Goppa keys and computes `W`, the span of the products
`f_a f_b` of the secret curve polynomials. Equality `dim W = 2D−t+1` implies
that the quadratic forms of the binary direct extraction exist at every
position. The rank computation stops once this upper bound is reached.

| Keys | Count | `dim W = 2D−t+1` | Positions with both forms |
|---|---:|---:|---:|
| Classic McEliece, shortened with the GIAJS and Weis dimensions | 10,000 | 10,000 | 80,000 / 80,000 |
| Classic McEliece `(m,t)`, `k` near the counting threshold | 400 | `min(C(k+1,2), 2D−t+1)` for all | — |
| Moderate shortened codes, `m` from 9 to 13 | 33,000 | 32,993 | all on keys at the bound |
| Unshortened codes with the TII challenge parameters | 3,000 | 3,000 | 24,000 / 24,000 |
| Small codes above the counting threshold | 98,008 | 94,065 | — |

In the small codes, every deficit at positive slack equals the number of linear
conditions from double and base points of the curve, including the point at
infinity. These points are repeated or zero columns of the public generator
matrix extended by its parity column. The remaining 26 deficient keys have zero
slack and one extra dependency each. The seven deficient moderate keys have
`t ≤ 20` and `k` at most three above the counting threshold: five have a double
point and two have zero slack. A SageMath cross-check recomputed 67 keys from
their coefficients and agreed on every rank and membership. The
[summary](goppa-direct-sufficient-condition/results/summary.md) contains all
tables.

## Synthetic direct recovery

The [synthetic collection](synthetic-direct-recovery/README.md) runs the direct
recovery from one local jet computed with the secret key, at a random position,
with a random reparametrization and on a random Frobenius branch. The recovery
uses only the public generator matrix and the jet, and every recovered key is
verified against the public code and audited against the secret support.

| Family | Fields | Sets | Keys recovered |
|---|---|---:|---:|
| Binary Goppa, square-root variant | GF(2^8) to GF(2^13), up to shortened `mceliece6688128` | 4 | 160 / 160 |
| Wild Goppa | q = 3, 4, 5 | 6 | 360 / 360 |
| Generic alternant, quadratic and cubic forms | q = 2, 3, 4, 5 | 7 | 540 / 540 |

The binary Goppa control with adjacent orders fails on all 20 keys because no
second form exists, as the paper predicts. A SageMath cross-check recomputed
36 recovered keys independently and rejected corrupted ones. The
[summary](synthetic-direct-recovery/results/summary.md) contains all tables.

## TII challenge recoveries

The [experiment specification](TII-EXPERIMENTS.md) records the shortening,
systematic coordinate sets, run settings, and full pipeline for each code.
Its [hardware and timing tables](TII-EXPERIMENTS.md#5-measured-hardware-and-runtimes)
give the measured durations and distinguish complete stages from component
timers and sums across separate runs.

The [TII results](tii-results/README.md) contain full equivalent binary Goppa
keys recovered from the original public matrices and single-position holdout
samples. No published secret challenge answers were inputs.

| Instance | Original code | Goppa degree | Saved holdout forms | Jet order | Full support |
|---|---|---:|---:|---:|---:|
| [TII-173](tii-results/tii-173/README.md) | `[96,40]` | 8 | 256 quintics in 35 variables | 296 | 96 |
| [TII-249](tii-results/tii-249/README.md) | `[235,107]` | 16 | 280 quintics in 50 variables | 290 | 235 |
| [TII-253](tii-results/tii-253/README.md) | `[214,142]` | 9 | 280 sextics in 46 variables | 82 | 214 |

TII-173 and TII-249 were also finished through quadratic bootstrapping and
minor codes on the same saved jets. Their minor ranks are 66 and 129, and
the retained supports match the direct results up to a Möbius transformation.
These additional runs did not repeat full-key extension. TII-253 now has
both full minor and direct recoveries: the direct run of 2 October extends
the shared public jet to order 198 and produces an independently verified
214-position degree-nine key. The [coverage table](TII-EXPERIMENTS.md#finishing-coverage)
and saved `evidence/minor-finishing/` outputs distinguish these scopes.

Each key retains its exact extension-field polynomial and original column
order. Saved independent verifiers establish equality of the reconstructed
binary Goppa parity-check row space with the original public check space,
distinct support, support avoidance, and irreducibility of the recovered
polynomial. TII-249 and TII-253 also retain the exhaustive coordinate-extension
candidate records. The larger-instance holdout reports certify sampled rank
and every original Hasse order; later geometry certificates are stored with
the recovery evidence.

[Kernel storage](tii-kernels/README.md) preserves the exact accepted samples.
[Source provenance and instructions](tii-keyrec-implementation/README.md)
identify the implementations used. The original artifacts were preserved;
the additional TII-253 direct run reused the holdout sample.

## Wild Goppa recovery over GF(4)

The [wild Goppa collection](keyrec-toy-wild-goppa/README.md) recovers complete
keys for GF(4) codes `[63,27]`, `[55,28]`, and `[64,28]`, with support in GF(64).
Every run computes one holdout kernel and continues one coherent branch jet.
Stationary quadratic rank three separates the three Frobenius branches; the
fixed continuation matrix has rank `k-2`.

Both finishing methods pass in all three cases. Bootstrapping yields full minor
GRS codes of dimensions 46, 42, and 47 respectively. Direct extraction uses
ordinary evaluation ratios, with jet orders 91, 83, and 93. The recovered radical
defines both a cube and a fourth-power Goppa description exactly equal to the
public GF(4) code. Independent verifiers check all six recovered keys.

The [result table](keyrec-toy-wild-goppa/results/README.md) gives kernel dimensions,
exported panel sizes, and measured runtimes. The large holdout uses OpenMP-enabled
M4RI; the public recovery and coefficient-level coherence checks use no generated
secret data.
