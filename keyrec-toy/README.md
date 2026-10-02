# Toy end-to-end key recovery

Four SageMath Jupyter notebooks start with a freshly generated public matrix,
compute a holdout kernel, recover a coherent local jet, and return an equivalent
support and multiplier.

| Notebook | Public code | Ambient field | `D` | Holdout `(d,s)` | Finishing |
|---|---|---|---:|---|---|
| [01-grs-char2.ipynb](01-grs-char2.ipynb) | `[24,6]` over GF(32) | GF(32) | 7 | `(3,1)` | Minor and direct |
| [02-alternant-char2.ipynb](02-alternant-char2.ipynb) | `[60,20]` over GF(8) | GF(64) | 39 | `(3,2)` | Minor and direct |
| [03-grs-odd.ipynb](03-grs-odd.ipynb) | `[17,4]` over GF(19) | GF(19) | 5 | `(3,1)` | Direct only |
| [04-alternant-odd.ipynb](04-alternant-odd.ipynb) | `[45,15]` over GF(11) | GF(121) | 29 | `(3,2)` | Direct only |

All four use fixed generation seed 0. These are worked examples, not estimates
of success probability. In each case `s(n-1)>dD`, so curve vanishing follows from
the ordinary interpolation bound. The alternant examples exercise Frobenius
branch separation; no synthetic tangent or jet enters recovery.

## Read and run

The notebooks include their executed outputs. Open one with a SageMath kernel
and run all cells, or execute the complete collection from this directory:

```sh
sage -python execute.py
```

To execute one notebook:

```sh
sage -python execute.py 02-alternant-char2.ipynb
```

Tested with SageMath 10.7. The runner also uses `nbformat`, `nbclient`, and
`jupyter_client`, and discovers a registered SageMath kernel automatically;
`--kernel NAME` selects one explicitly. Every notebook gets a fresh kernel.
The runner saves outputs only after successful execution and writes a summary
with content hashes under [results/](results/).

All implementation code is in [reference.py](reference.py). The notebook working
directory must be this folder. A C++ compiler and a precomputed kernel archive are unnecessary
for these four examples.

## What is checked

1. **Public holdout.** Impose the retained Hasse conditions, verify the kernel
   residual, and check zero evaluation at the held-out point. Retained unit
   columns are handled exactly by the permitted monomial exponents.
2. **Local branch recovery.** Measure gradient rank and restricted stationary
   quadratic rank, and recover the projective branches by the quadratic quotient
   eigenvalue method. Select one branch using only these public equations.
3. **Coherent continuation.** Check that the fixed regular/next-order matrix has
   exactly the point/tangent plane as its kernel. Continue past saturation of the
   coefficient span and independently substitute the final series into every
   selected polynomial.
4. **Characteristic-two minor finishing.** Recover the full extension-field
   quadratic ideal from public values and a jet prefix, check tangent rank at
   every position, form all square-root minors, and verify exact equality with
   the recovered GRS code.
5. **Direct finishing.** Solve for the high-contact sections using exactly the
   sufficient jet prefix of order `delta*D-1`. Extract support ratios, handle
   infinity, and solve for a public-code multiplier. The odd GRS example uses
   cubic sections; the other three use quadrics.
6. **Independent audits.** After recovery, use the generated secret to check the
   polynomial identities, all branch planes, a coherent scalar/parameter series,
   and global Frobenius/Möbius equivalence of the recovered support. For the minor
   route, compare the recovered ideal with the full secret-curve ideal and check
   all tangent planes against one branch. The two characteristic-two finishing
   routes must agree up to a common Möbius transformation.

Both alternant cases additionally verify exact equality between the recovered
subfield code and the public code. A deliberately corrupted multiplier must fail
the public verification in every notebook.

The recovery functions never receive the generated secret. The secret is kept
in a separate notebook variable and passed only to the final audit functions.
The notebooks expose every stage, its ranks, and its input boundary. Complete
recovered supports and multipliers remain in `direct` and, where applicable,
`minor`; compact machine-readable records are written to `results/`.

## Scope and interpretation

These examples cover generic GRS subcodes and generic alternant codes.
Characteristic two alone does not imply the special binary-Goppa relation
`v_(2r)=v_r^2`: these notebooks use the general stationary continuation, and their
direct extraction uses adjacent orders and ordinary ratios. Binary Goppa recovery
is implemented in the [C++ collection](../keyrec-toy-binary-goppa/README.md).

In the `[24,6]` GRS example, the public points already determine the quadratic
ideal, so the bootstrap uses order zero. The `[60,20]` alternant example requires
order 19 and demonstrates selection of one branch by a recovered jet. Its direct
finishing uses order 77. The odd examples use orders 14 and 57 respectively.

No preliminary shortening is used in these four notebooks: every coordinate
of each displayed public code is recovered. They do not test restoration of
removed coordinates or performance at cryptographic sizes.
