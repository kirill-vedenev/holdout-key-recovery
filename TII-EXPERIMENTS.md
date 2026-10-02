# TII experiments: shortening, coordinates, parameters, and recovery pipelines

This file describes the **completed TII-173, TII-249, and TII-253 runs** whose
keys, holdout samples, and implementations are distributed here. All three
used shortening before computing the holdout kernel, and all removed support
coordinates were subsequently recovered. The final keys cover the original,
unshortened public codes.

The [measured hardware and runtimes](#5-measured-hardware-and-runtimes)
separate GPU kernel work, CPU polynomial-generator work, and local key recovery.

Each [implementation README](tii-keyrec-implementation/README.md) includes a
short saved-key check and optional audits of the intermediate results.

**Every position in this file is zero based, and every range is inclusive.**
“Original position” means a column of `pk_McEliece_N.txt`. “Shortened position”
means a column of the saved `public.json:generator`. These two numberings must
not be interchanged. In all three runs, the held-out position was **original
0, shortened 0**, outside the systematic information set.

The tables were checked against the saved public matrices, coordinate maps,
source, and accepted run records. Original runs were preserved during
packaging. An additional TII-253 direct finishing was executed on 2 October
2026 using the existing holdout sample and public jet prefix.

## Finishing coverage

| Instance | Direct finishing | Minor-code finishing | Scope of the recorded comparison |
|---|---|---|---|
| TII-173 | Passed; support through jet order 148, curve through 296 | Passed; bootstrap depth 58, GRS `[91,66]` | Both recover equivalent support on all 91 retained positions |
| TII-249 | Passed; support through jet order 290 | Passed; bootstrap depth 113, GRS `[178,129]` | Both recover equivalent support on all 178 retained positions |
| TII-253 | Passed; support through jet order 198, followed by full-key extension | Passed; bootstrap depth 81 from an order-82 jet, GRS `[118,90]`, followed by full-key extension | Both full 214-position keys pass verification against the original public key |

The TII-173 and TII-249 minor finishings were run on 1 October 2026 using the
same saved public matrices and coherent jets as their direct runs. The
archived comparison verifies Möbius-equivalent retained supports, with
Frobenius power 0 in both cases. Full-key extension was not rerun separately
for the minor route. Thus both finishings were tested through retained-support
recovery, not as two separate complete executions from holdout to full key.

The original minor results, including quadrics, tangent representatives,
minor bases, recovered supports and multipliers, are included under
[TII-173 minor finishing](tii-results/tii-173/evidence/minor-finishing/) and
[TII-249 minor finishing](tii-results/tii-249/evidence/minor-finishing/).
Their input hashes match the already distributed public matrices and jets.
The original mathematical implementation and run commands accompany each
instance under `recovery/minor/`.

For TII-253, the direct run continued the shared order-82 public jet through
`2D=198`, with `D=99`, and independently recovered a full degree-nine Goppa
key. Its inputs did not include the minor-route support or recovered key.
The [direct key](tii-results/tii-253/direct-recovered-sk.json),
[public-key verification](tii-results/tii-253/direct-verification.json), and
[run summary](tii-results/tii-253/evidence/direct-finishing/result-summary.json)
are retained separately from the original minor results.

## 1. Public inputs and shortening

The public-key file contains a binary **parity-check** matrix `H`, followed
by the extension-field modulus. The starting code is `C = ker_GF(2)(H)`.
For a shortening set `S`, the code used in the holdout computation is

```text
C_short = { (c_j)_{j not in S} : c in C and c_j = 0 for every j in S }.
```

Thus the removed positions are constrained to zero before deletion; this is
shortening, not unrestricted puncturing of codewords. Preparation takes the
binary kernel of `H` restricted to the retained columns. Retained columns
stay in ascending original order. Row operations then put the chosen
information set into systematic form without changing that column order.

| Instance | Original `(m,t,n,k)` | Public `H` shape and rank | Removed original positions `S` | Number removed | Shortened `[n_s,k_s]` |
|---|---|---|---|---:|---|
| TII-173 | `(7,8,96,40)` | `56 × 96`, rank 56 | `56..60` | 5 | `[91,35]` |
| TII-249 | `(8,16,235,107)` | `128 × 235`, rank 128 | `128..184` | 57 | `[178,50]` |
| TII-253 | `(8,9,214,142)` | `72 × 214`, rank 72 | `72..167` | 96 | `[118,46]` |

Here `m` is the extension degree and `t` the Goppa-polynomial degree. The
extension fields and integer encodings are:

| Instance | Extension field | Defining polynomial over GF(2) | Integer modulus |
|---|---|---|---|
| TII-173 | GF(128) | `z^7 + z + 1` | `0x83` |
| TII-249, TII-253 | GF(256) | `z^8 + z^4 + z^3 + z^2 + 1` | `0x11d` |

An integer `a` represents `sum_i ((a >> i) & 1) z^i`. The saved field and
Goppa-polynomial coefficient lists are low degree first. Holdout polynomial
coefficients themselves are in GF(2); branch, jet, and support calculations
use the specified extension field.

The preparation code selects the first `k-k_s` free coordinates of the
original binary nullspace construction, excluding the held-out position.
For these inputs that rule gives exactly the three sets above. It then
selects an information set avoiding the held-out column and normalizes its
columns to the identity. These are fixed recorded choices, not a claim of
an optimal shortening search.

## 2. Systematic positions and the precise variable basis

Let `G` be the saved `k_s × n_s` generator. If `I` is the ordered information
set in the table below, **`G[:,I] = identity(k_s)`**. Polynomial variable
`X_a` refers to row `a` of this exact matrix. The holdout equations are
evaluated on its columns `p_j = G[:,j]`.

| Instance | Retained original positions, in order | Systematic positions, original numbering | Systematic positions, shortened numbering | Identity-column correspondence |
|---|---|---|---|---|
| TII-173 | `0..55`, then `61..95` | `61..95` | `56..90` | Original `61+a`, shortened `56+a`, is `e_a`, for `a=0..34` |
| TII-249 | `0..127`, then `185..234` | `185..234` | `128..177` | Original `185+a`, shortened `128+a`, is `e_a`, for `a=0..49` |
| TII-253 | `0..71`, then `168..213` | `168..213` | `72..117` | Original `168+a`, shortened `72+a`, is `e_a`, for `a=0..45` |

Equivalently, the saved matrices have the block form `G = [A | I]`, with
the retained nonsystematic original columns first. Their shapes are
`35 × 91`, `50 × 178`, and `46 × 118` respectively. The identity-column
equalities were checked directly in the saved matrices.

For conversion to **one-based** notation: the removed ranges are `57..61`,
`129..185`, and `73..168`; the systematic original ranges are `62..96`,
`186..235`, and `169..214`; the held-out position is 1.

The complete matrices and maps are the authoritative input records:
[TII-173 public.json](tii-kernels/tii-173/public.json),
[TII-249 public.json](tii-kernels/tii-249/public.json), and
[TII-253 public.json](tii-kernels/tii-253/public.json).
`original_indices[j]` maps shortened position `j` to its original label;
`shortened_original_indices` lists the removed labels; `information_set`
uses shortened labels, and `information_set_original` uses original labels.
Changing the generator row basis requires transforming the polynomial
coefficients as well.

## 3. Holdout parameters actually used

The unknown forms are homogeneous squarefree degree-`d` polynomials in
`k_s` variables. Multiplicity `s` is imposed at every other retained public
column and at the additional public point
`p_infinity = sum_j p_j` over GF(2). No equation is imposed at the held-out
column 0 during kernel construction. Infinity is labeled `-1` in
`operator.txt`; it is neither an original column nor a shortened coordinate.

| Parameter | TII-173 | TII-249 | TII-253 |
|---|---:|---:|---:|
| Variables `k_s` | 35 | 50 | 46 |
| Degree `d` | 5 | 5 | 6 |
| Multiplicity at retained columns other than 0 | 4 | 4 | 5 |
| Multiplicity at infinity | 4 | 4 | 5 |
| Multiplicity imposed at column 0 | 0 | 0 | 0 |
| Hasse orders retained in the reduced solver | 1, 3 | 1, 3 | 2, 4 |
| Hasse orders checked by independent verification | 0, 1, 2, 3 | 0, 1, 2, 3 | 0, 1, 2, 3, 4 |
| Active constraint points, including infinity | 56 | 128 | 72 |
| Reduced matrix rows | 368,480 | 2,515,200 | 11,823,840 |
| Monomials / matrix columns `binomial(k_s,d)` | 324,632 | 2,118,760 | 9,366,819 |
| Minimum sample size | 210 | 280 | 280 |
| **Actual saved independent forms** | **256** | **280** | **280** |
| Saved monomial ordering | WHK1, lexicographic; original THK1 also included for verification | THK1, colex | THK1, colex |
| Held-out point mask | 2,967,760,381 | 505,067,484,167,246 | 42,945,927,196,204 |
| Infinity point mask | 6,212,682,838 | 894,692,465,462,125 | 1,897,043,500,789 |

Bit `a` of a point mask is coordinate `a` in the saved row basis. The
systematic unit columns satisfy these squarefree Hasse constraints
identically, accounting for the difference between all constraint points
and the active points explicitly assembled. The reduced derivative-order
selection uses the exact Hasse-Euler relations; the independent verifier
checks the omitted orders as well.

The configured minimum is `max(10*binomial(m,2), 5*k_s)`, using the **shortened**
dimension. TII-173 rounded its minimum up to 256. For TII-249 and TII-253,
`public.json:default_elements = 512` is a preparation default; the accepted
reconstruction explicitly requested **280**. The 512-bit computational block
width below is also distinct from the number of output polynomials. Samples
from different held-out positions were not pooled.

The saved parameter choices have the following interpolation counts, with
shortened curve degree bound `D = n_s - 2t - 1`:

| Quantity | TII-173 | TII-249 | TII-253 |
|---|---:|---:|---:|
| `D` | 74 | 145 | 99 |
| Holdout composition bound `dD` | 370 | 725 | 594 |
| Imposed zeros counted with multiplicity `s*n_s` | 364 | 712 | 590 |
| Shortfall from the sufficient count `dD+1` | 7 | 14 | 5 |

The count includes `n_s-1` finite columns and infinity. In each case the
shortfall is at most `t`, as used by the configured heuristic parameter
criterion. **The initial ordinary interpolation count alone does not certify
curve identities.** Exact kernel membership, successful recovery, and the
later curve audits are recorded separately. The sample sizes are not claims
about the full kernel dimension. Only TII-173's dense elimination measured
the full rank/nullity: **312,551 / 12,081**.

## 4. Computational run settings

### TII-173: exact panel elimination

The accepted run used two H200 NVL GPUs on one host, devices `0,1`, CUDA 13.3,
and `sm_90` compilation. Its settings were `count=256`, `seed=17320260930`,
`panel-bits=8`, `max-gib=64`, and `checkpoint-every=2048` panels. The matrix
occupied 14,954,392,320 bytes, split by rows across the two GPUs.

After elimination, the solver chose free-coordinate assignments with an
identity minor to guarantee sample independence. The CPU verifier used
16 threads. The colex THK1 output was explicitly permuted into the WHK1
lexicographic order read by this instance's recovery code.

### TII-249 and TII-253: GPU → CPU → GPU

Both used the preserved square cuFFT worker and a 512-bit block (eight
64-bit words). The three stages were:

1. Generate `S_i = Z^T A^i Y` on a GPU, where `A=R*M` is the square operator
   formed from the public holdout matrix and deterministic Toeplitz compression.
2. Compute a polynomial matrix generator using the patched CADO-NFS
   `lingen_b64` on CPUs.
3. Reconstruct polynomial vectors on a GPU, check the literal uncompressed
   constraints, then independently verify all original Hasse orders on CPUs.

| Setting | TII-249 | TII-253 |
|---|---|---|
| GPU used per numerical stage | One H200 | One H200 |
| Initial-panel seed `seed_y` | `24920260930` | `25320260930` |
| Compressor seed | `6073466052459188813` | `6073466052459188813` |
| Panel PRNG | `splitmix64-indexed-v1` | `splitmix64-indexed-v1` |
| GPU `batch_words` | 2 | 2 |
| Sequence length passed to CADO | 8,705 terms | 36,846 terms |
| First sequence term | `S_0` | `S_0` |
| CADO parameters | `prime=2 m=512 n=512 rhs=none` | Same |
| CADO thread grid | `thr=4x8` | `thr=48x1` |
| CPU threads specified by the grid product | 32 | 48 |
| CADO memory setting | `max_ram=128` | `max_ram=128` |
| CADO tuning | `recursive:100000,notiming:1` | Same |
| Reconstruction orientation | `transpose=1` | `transpose=1` |
| Reconstruction output count | 280 | 280 |
| Independent Hasse-verifier threads | 32 | 32 |

The projection seed is derived as `seed_z = seed_y XOR 0x9e3779b97f4a7c15`.
CADO's `m=n=512` are block dimensions, unrelated to the code's extension
degree or length. The distributed CPU recipes explicitly set
`OMP_NUM_THREADS=32` and `48`, respectively, and use `OMP_DYNAMIC=FALSE`,
`OMP_PLACES=cores`, `OMP_PROC_BIND=spread`, and `OMP_MAX_ACTIVE_LEVELS=1`.
The recorded TII-253 CPU run also used NUMA interleaving; this scheduling
choice does not change the algebraic inputs.

The pinned CADO revision is `2d98fe176c11342cfa5d6c0d795a992cbe433c11`, with
the included coefficient-parallel patch. The accepted TII-249 sequence was
corrected from an earlier stream starting at `S_1` before entering CADO;
the included worker and current instructions produce `S_0` directly.
The accepted reconstructed generators contained 4,108 and 17,900 coefficient
matrices respectively; these are outputs, not requested input parameters.

All three CPU finishing stages used SageMath 10.7. TII-173's substitution
helpers use NumPy and M4RI; the larger instances use the included colex
series helper with **four threads**. Each continuation selected the first
returned branch (index 0), fixed even coefficients by componentwise squaring
`v_(2r) = v_r^2`, and set the two nonpivot coordinates of each odd solve to
zero. It used both the odd equations and all next-even equations. TII-173's
full-support-vector helper has fallback random seed **173**, used only if
its solution space has dimension greater than one. TII-253's bootstrap
fixes Sage's random seed to **3801** and allows up to **256**
Sidelnikov-Shestakov recovery attempts.

The additional TII-173/TII-249 minor runs also use Sage seed 3801, with
bootstrap depths 58 and 113 and support-recovery attempt limits 128 and 256.

## 5. Measured hardware and runtimes

These are measurements from the successful runs of 30 September and
1 October 2026, plus the TII-253 direct finishing of 2 October 2026.
An enclosing stage time includes its component timers.
Sums across separately executed stages exclude transfers, scheduling gaps,
development, validation experiments, and interrupted or superseded attempts.
They are not measurements of one uninterrupted end-to-end invocation.

### Hardware assigned to each stage

| Work | Recorded hardware and concurrency |
|---|---|
| TII-173 holdout elimination | **Two NVIDIA H200 NVL GPUs** on one host, CUDA 13.3, device architecture `sm_90`; CPU Hasse verification used 16 threads |
| TII-249 sequence generation and polynomial reconstruction | **One NVIDIA H200 per stage**; the accepted CPU generator used `thr=4x8` (32-thread grid), and Hasse verification used 32 threads |
| TII-253 sequence generation and polynomial reconstruction | **One NVIDIA H200 NVL per stage**, CUDA 13.3; Hasse verification used 32 threads |
| TII-253 CPU polynomial generator | **AMD EPYC-Genoa Processor**, 96 logical CPUs reported by the KVM guest; the successful CADO run used **48 threads**, `thr=48x1`, with NUMA interleaving |
| Jet, support, and full-key recovery | Local macOS ARM64, SageMath 10.7; TII-249 and TII-253 used four native series threads; these stages used no GPU |

The TII-253 CPU machine record reports two sockets, 24 cores per socket,
and two threads per core: **48 guest-reported cores / 96 logical CPUs**.
The recorded model is EPYC Genoa. “96 physical Ryzen cores” would not match
this record. The exact CPU model of the TII-249 CADO host and the local
Apple CPU model were not preserved in the available run summaries, so no
model-specific claim is made for them.

### Holdout computation: elapsed times and accepted outputs

| Instance / stage | Recorded elapsed time | Scope and result |
|---|---:|---|
| **TII-173 complete holdout run** | **175.71 s (2m 56s)** | Elimination, sample extraction, independent verification, and export; **256 independent quintics** |
| TII-173 GPU elimination | 172.54 s | Includes periodic checkpoints; component of the 175.71 s total |
| TII-173 CPU Hasse verification | 0.553 s | Component of the same total |
| TII-249 GPU sequence generation | **Approximately 27 min** | Earlier accepted-seed sequence; only this rounded timing survives in the contemporaneous summary |
| **TII-249 accepted solve/reconstruction pipeline** | **2,143.42 s (35m 43s)** | Sequence adaptation/checking, CPU generator, GPU reconstruction, and verification; **280 independent quintics**; excludes the earlier long sequence generation |
| TII-249 GPU polynomial reconstruction | 680.47 s (11m 20s) | Component of the 2,143.42 s pipeline |
| TII-249 CPU Hasse verification | 3.980 s | Component of the same pipeline |
| TII-253 GPU sequence generation | **24,326.59 s (6h 45m 27s)** | 36,846 terms for the seed used by the successful CPU run |
| TII-253 CPU polynomial generator | **14,132.62 s (3h 55m 33s)** | Successful 48-thread EPYC run; produces 17,900 coefficient matrices |
| **TII-253 GPU reconstruction and acceptance** | **11,580.99 s (3h 13m 1s)** | Reconstruction, output publication, and independent CPU verification; **280 independent sextics** |
| TII-253 numerical GPU worker | 11,558.89 s | Component of the reconstruction/acceptance stage |
| TII-253 CPU Hasse verification | 19.394 s | Component of the same stage |

For the complete successful holdout work, the recorded totals are:

- **TII-173: 175.71 seconds**, directly measured by the completed-run record.
- **TII-249: approximately 1h 03m**, from approximately 27 minutes of sequence
  generation plus 35m 43s for the accepted later pipeline. The exact standalone
  CADO duration is not separately recorded; subtracting the reconstruction
  time would also leave checks and other overhead, not just the CPU solver.
- **TII-253: 50,040.20 seconds, approximately 13h 54m**, the sum of the
  successful sequence, CPU-generator, and reconstruction/acceptance stages.
  This excludes the alternative-seed sequence and interrupted CPU attempts.

Every accepted sample passed independent rank and all-original-Hasse-order
checks. Their downstream pipelines recovered full equivalent keys on 96,
235, and 214 original positions, respectively.

### Local CPU finishing and audits

The kernel was already available for these measurements. Startup-inclusive
module wall times and internal mathematical timers are distinguished below.

| Instance / recorded work | Time | Timing boundary |
|---|---:|---|
| TII-173 normalized jet through order 296 | 67.41 s | Internal timer, including full substitution verification |
| TII-173 direct curve extraction core | 2.01 s | Internal core timer; excludes surrounding public-column checks |
| TII-173 deshortening, multipliers, and Goppa checks | 1.02 s | Internal timer |
| TII-173 full chain audit | 63.68 s | Separate audit; includes 32.10 s for all polynomial curve identities |
| TII-173 independent full-key verification | 0.006 s | Standard-library verifier's internal timer |
| TII-249 archived jet/support reproduction module | 89.27 s | Includes Sage startup, Hasse extraction, native validation, branch/jet/support stages and recorded audits; excludes native compilation |
| TII-249 archived full-key reproduction module | 18.51 s | Includes Sage startup, coordinate extension, full key, and independent verification |
| TII-173 additional minor finishing | 2.412 s | Original internal timer; quadrics, tangents, minor code and retained support; excludes Sage startup, the shared jet and full-key extension |
| TII-249 additional minor finishing | 7.087 s | Same scope; reuses the direct run's saved jet |
| TII-253 recorded successful frontend and finishing commands | 86.14 s | Sum of command wall times, including Sage startup, helper compilation, native validation and the listed finishing audits; starts from already prepared first-jet tensors and branches |
| TII-253 direct-run jet continuation to 198 | 83.26 s | Internal timer, including 22.05 s for cache setup and checking the saved order-82 prefix; excludes Sage startup |
| TII-253 direct support extraction | 2.663 s | Internal timer, including retained public-code checks; command wall time 7.47 s |
| TII-253 direct full-key extension and Goppa extraction | 13.754 s | Internal timer; command wall time 18.11 s |
| TII-253 direct full-key verification | 0.296 s | Standard-library verifier, including all extension candidates; command wall time 0.374 s |

TII-173 has no archived complete local wall timer; its internal timers should
not be presented as a complete startup-inclusive pipeline time. For TII-249,
the two already archived reproduction modules sum to **107.78 s (1m 48s)**.
Those reproductions were performed during the original experiments, not
during preparation of this repository. The original first-run internal
measurements are retained in the linked timing record as well.

For TII-253, the 86.14 s comprises 32.81 s of frontend commands and 53.33 s
of finishing commands. It excludes the earlier first-jet preparation. Separate
internal timers are 12.36 s for panel/THK1 encoding verification, 1.80 s for
the first-jet criterion, and 0.244 s for independent minor verification. No
complete standalone elapsed timer was retained for the initial Hasse-tensor
extraction, so these values do not establish a complete contiguous wall time
from kernel input through every check. The unsuccessful initial geometry
attempt and gaps between commands are excluded.

Exact values, timing boundaries, hardware metadata, and hashes of the source
records are saved in [TII-173 timings](tii-results/tii-173/timings.json),
[TII-249 timings](tii-results/tii-249/timings.json), and
[TII-253 timings](tii-results/tii-253/timings.json). These records contain
measurements and provenance, without the server orchestration or host addresses.

The new direct run used local macOS ARM64, SageMath 10.7 with Python 3.13.3
and NumPy 2.2.4, and four native series threads. It reused the original
holdout and initial jet. Its exact CPU model was unavailable to the process;
no whole-pipeline wall-time total is inferred from the mixed timer scopes.

## 6. Full pipeline: TII-173

1. **Prepare the public code.** Shorten original coordinates `56..60`, retain
   91 columns, and make original coordinates `61..95` systematic. Preserve
   original position 0 as the holdout.
2. **Compute and verify the holdout sample.** Use the degree-five,
   multiplicity-four constraints above. Exact GPU elimination returns 256
   independent forms. Check every original Hasse order on CPU, then convert
   colex coefficients to WHK1 lexicographic order.
3. **Recover a local branch and jet.** The gradient has rank 27 and its
   annihilator dimension is 8. The quotient by the public point yields seven
   branches. The fixed normalized continuation matrix has rank 33 and a
   two-dimensional point/tangent kernel. Continue through order **296**,
   retaining coefficients `0..296`.
4. **Extract the retained support and curve directly.** With `D=74`, quadratic
   sections use orders through `2D=148`; this implementation also reconstructs
   the curve using orders through `4D=296`. The section system has 630
   quadratic unknowns and 149 contact rows. Recover 91 distinct projective
   support values and verify exact equality with the shortened public code.
   One support value is infinite in this initial projective chart.
5. **Choose the Goppa affine chart and undo shortening.** The public infinity
   point has recovered parameter `beta_infinity=30`. Apply
   `alpha=1/(beta+30)`, making all retained coordinates finite. For each of
   the five removed positions, test all 37 unused field values against the
   public Vandermonde multiplier equations. Each has exactly one passing
   value: 185 candidate tests in total.
6. **Recover the full key.** Solve the `640 × 96` full multiplier system,
   whose nullity is 1. Interpolate the reciprocal dual multipliers, recover
   the monic degree-eight Goppa polynomial by an exact characteristic-two
   square root, and obtain all 96 support coordinates in original order.
7. **Audit the result.** The independent verifier proves equality with the
   original rank-56 binary check space. The chain audit reconstructs the
   degree-74 curve, checks all 371 coefficients of every one of the 256
   quintic compositions, verifies every branch and continuation kernel, and
   checks all 10,395 saved coordinate-jet coefficients for coherence. A
   separate support-only check uses only the prefix through order 148.

Sources and commands: [TII-173 implementation](tii-keyrec-implementation/tii-173/README.md).
Saved evidence: [direct extraction](tii-results/tii-173/evidence/direct-recovery.json),
[deshortening](tii-results/tii-173/evidence/deshortening-report.json),
[chain audit](tii-results/tii-173/evidence/chain-verification.json), and
[full key](tii-results/tii-173/recovered-sk.json).

## 7. Full pipeline: TII-249

1. **Prepare the public code.** Shorten original coordinates `128..184`;
   original `185..234` becomes the systematic block of the retained
   `[178,50]` code. Hold out original position 0.
2. **Compute and verify the holdout sample.** Run the three computational
   stages with degree five, multiplicity four, and the recorded seeds.
   Reconstruct 280 independent forms and verify every original Hasse order.
   Restore or read the exact colex panel and its monomial map.
3. **Recover a local branch and jet.** Extract the Hasse tensors. The
   gradient rank is 41, the annihilator dimension is 9, and the first-jet
   quotient yields eight branches. The continuation rank is 48, with only
   the two point/tangent freedoms. Continue through order **290**, saving
   coefficients `0..290`.
4. **Extract the retained support directly.** Here `D=145` and the chosen
   depth is `2D=290`. Solve the `291 × 1275` quadratic contact system for
   extreme sections and take square roots of their ratios to obtain 178
   distinct support values. The retained multiplier system is `1600 × 178`,
   rank 177, nullity 1. Verify exact shortened-code equality. A separate
   audit reconstructs the degree-145 curve and fits all 14,550 saved jet
   coefficients to one scalar series and one invertible parameter series.
5. **Choose the affine chart and restore 57 positions.** Evaluate the
   sections at public infinity, obtaining `beta_infinity=186`, and set
   `alpha=1/(beta+186)`. Factor the common retained multiplier system once.
   For each removed original coordinate `128..184`, test 78 unused finite
   values and infinity using all 32 moment equations: **79 candidates per
   position, 4,503 total**. Each position has one finite solution; every
   infinity candidate is rejected.
6. **Recover and verify the full key.** The full multiplier system is
   `3424 × 235`, rank 234, nullity 1. Interpolating reciprocal multipliers
   and taking the exact square root yields a monic irreducible degree-16
   Goppa polynomial. The independent verifier checks all 235 coordinates,
   the multiplier conventions, polynomial irreducibility and support
   avoidance, and equality with the original rank-128 binary check space.
   It can also recheck every saved coordinate-extension candidate.

This run used the direct-support finishing route. Its saved coherence audit
does not claim a full degree-725 coefficient audit of every quintic curve
composition; the full-key verification is an independent exact certificate
of the recovered decoding key.

Sources and commands: [TII-249 implementation](tii-keyrec-implementation/tii-249/README.md)
and [three-stage holdout](tii-keyrec-implementation/tii-249/holdout/README.md).
Saved evidence: [local branches](tii-results/tii-249/evidence/jet-support/local-branches.json),
[direct support](tii-results/tii-249/evidence/jet-support/direct-support.json),
[full recovery](tii-results/tii-249/evidence/full-key/recovery-report.json), and
[full key](tii-results/tii-249/recovered-sk.json).

## 8. Full pipeline: TII-253

1. **Prepare the public code.** Shorten original coordinates `72..167`;
   original `168..213` is the systematic block of the retained `[118,46]`
   code. Hold out original position 0.
2. **Compute and verify the holdout sample.** Run the three computational
   stages with degree six, multiplicity five, and the recorded seeds.
   Reconstruct 280 independent sextics. Verify the literal public constraints,
   all original Hasse orders, and the panel/THK1 coefficient correspondence.
3. **Recover a local branch and jet.** The gradient rank is 37, the
   annihilator dimension is 9, and the quotient gives eight branches. The
   normalized continuation rank is 44, leaving the point/tangent plane.
   Continue through order **82**, saving coefficients `0..82`. Full
   substitution later verifies all `280 × 83 = 23,240` field equalities.
4. **Bootstrap the quadratic ideal.** With `D=99`, use all 118 public points
   and **81 additional jet-contact orders**: `118+81=199 > 2D=198`.
   Thus the configured bootstrap depth is 81, even though the saved jet
   extends to 82 (the continuation works in odd/even pairs). Include all
   1,081 quadratic monomials, including squares. The constraint rank is 190
   and the quadratic ideal dimension is 891. All remaining saved contacts
   also pass. The observed first depth attaining rank 190 was 72; the run
   used 81, as required by its ordinary interpolation bound.
5. **Recover retained support through the minor code.** Align tangent planes
   at all 118 columns (rank 44, dimension 2), form all 1,035 coordinate-pair
   minors, and take their characteristic-two square roots. The resulting
   code is the full GRS `[118,90]` code. Sidelnikov-Shestakov recovery returns
   the 118 retained support values and multipliers. Independently verify
   the minor entries, exact GRS equality, and shortened public-code equality.
6. **Audit the curve and recover public infinity.** Reconstruct the
   degree-99 curve and verify all 3,818 saved jet coefficients are coherent.
   Check that all 891 quadrics are global identities spanning the quadratic
   ideal. Evaluating the curve at all 257 projective field points identifies
   public infinity uniquely at `beta_infinity=191`. The coherent jet gives
   held-out multiplicity 83; with the saved public Hasse constraints,
   `590+83=673 > 594` certifies all 280 holdout identities on this curve.
7. **Restore the full support and key.** Use `alpha=1/(beta+191)`. For each
   removed original coordinate `72..167`, test 138 unused finite values and
   infinity using all 18 moment equations: **139 candidates per position,
   13,344 total**. Each has exactly one finite solution. The full multiplier
   system is `2556 × 214`, rank 213, nullity 1. Recover the monic irreducible
   degree-nine Goppa polynomial and verify equality with the original
   rank-72 binary check space on all 214 coordinates. The independent
   verifier can recheck all extension candidates.

The original 1 October run used quadratic bootstrapping and minor-code finishing. The
first-jet record's `branch_solver` field is passed to continuation as
`branches`; the per-instance instructions include this interface conversion.

Sources and commands: [TII-253 implementation](tii-keyrec-implementation/tii-253/README.md)
and [three-stage holdout](tii-keyrec-implementation/tii-253/holdout/README.md).
Saved evidence: [first-jet criterion](tii-results/tii-253/evidence/criterion.json),
[bootstrap](tii-results/tii-253/evidence/bootstrap.json),
[geometry audit](tii-results/tii-253/evidence/geometry-verification.json),
[full recovery](tii-results/tii-253/evidence/full-support/recovery-report.json), and
[full key](tii-results/tii-253/recovered-sk.json).

### Additional direct finishing, 2 October 2026

The direct run preserves the saved coefficients through order 82 and computes
orders 83–198 using the same public branch and 280-form kernel. Quadratic
direct extraction solves a `199 × 1081` contact system of rank 190 and
recovers all 118 retained support positions. Its retained multiplier system
is `828 × 118`, rank 117, nullity 1.

The new sections locate public infinity at `beta=7`, giving the affine chart
`alpha=1/(beta+7)`. The original public parity check and this directly recovered
support then restore the 96 shortened positions and produce a full equivalent
degree-nine Goppa key. The final standalone verifier confirms exact equality
with the original public code on all 214 positions.

The [direct implementation](tii-keyrec-implementation/tii-253/README.md#direct-finishing)
includes the actual run instructions. The original minor-route key remains
`recovered-sk.json`; the additional key is `direct-recovered-sk.json`.

## 9. Reading the distributed artifacts

The input/output chain for every instance is:

```text
original public parity check + field definition
  -> shorten and normalize the public generator
  -> one-position holdout equations and verified polynomial sample
  -> first branch and coherent jet
  -> retained support (direct sections, or quadrics and minor GRS code)
  -> affine chart and restoration of every shortened coordinate
  -> full multipliers and Goppa polynomial
  -> independent verification against the original public parity check
```

[tii-results](tii-results/README.md) contains the original public inputs,
recovered full keys, and saved audits. [tii-kernels](tii-kernels/README.md)
defines coefficient packing, monomial order, chunk reassembly, and the
derived panel layout. [tii-keyrec-implementation](tii-keyrec-implementation/README.md)
contains the computational sources and direct commands, including dependency
versions and the CADO patch. Starting from saved kernels skips only their
recomputation. Unpacking coefficients does not perform a recovery stage.

Recovered keys and audit evidence are outputs, not inputs to a fresh
recovery. All final supports are finite and in original public column order;
`dual_grs_multipliers[0]=1` fixes the common multiplier scale. A common
coordinate change or Frobenius automorphism can yield a different equivalent
key, so public-code equality is the final correctness condition.
