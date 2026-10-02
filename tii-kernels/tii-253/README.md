# TII-253 holdout sample

**280 independent polynomials**, degree 6,
46 variables, 9,366,819 squarefree monomials;
held-out original column 0. Coefficients are binary. The file format is
**THK1**, using **colex** monomial order.
This is a sample from the kernel, not its full basis.

- `public.json`: exact variable basis, shortened public matrix, field modulus,
  original coordinate map, held-out point, and multiplicities.
- `operator.txt`: the complete public Hasse constraints.
- `verification.json`: saved independent membership and rank check.
- `storage.json`: dimensions, ordered chunks, original SHA-256, and byte sizes.

The original packed file has **327,838,860 bytes** and SHA-256
`8cf2d5e11f8de06411201e6af4968fcf1acfb2179a7c624487cf32b06e0de349`. It is stored in 10 file(s), listed
in `storage.json`. Concatenation preserves every original coefficient byte.
The [root storage specification](../README.md) defines the header, bit order,
monomial order, panel layout, and coordinate conventions.

From `tii-kernels/`, restore the input into a fresh directory:

```sh
python3 materialize.py tii-253 --output tii253-input
```

For the original jet implementation, also restore the exact panel and monomial map:

```sh
sage -python materialize.py tii-253 --output tii253-panel-input --panel
```

The [implementation](../../tii-keyrec-implementation/tii-253/README.md)
contains direct verification and recovery commands. The
[final key and later audits](../../tii-results/tii-253/README.md) are stored
separately. The holdout report's historical `curve_identity_certified` flag
describes that stage only.
