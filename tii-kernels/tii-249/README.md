# TII-249 holdout sample

**280 independent polynomials**, degree 5,
50 variables, 2,118,760 squarefree monomials;
held-out original column 0. Coefficients are binary. The file format is
**THK1**, using **colex** monomial order.
This is a sample from the kernel, not its full basis.

- `public.json`: exact variable basis, shortened public matrix, field modulus,
  original coordinate map, held-out point, and multiplicities.
- `operator.txt`: the complete public Hasse constraints.
- `verification.json`: saved independent membership and rank check.
- `storage.json`: dimensions, ordered chunks, original SHA-256, and byte sizes.

The original packed file has **74,156,620 bytes** and SHA-256
`4084644d4d079d8708aad0cf041296f5ddcd6003aa42e68c07ccccb63be74630`. It is stored in 3 file(s), listed
in `storage.json`. Concatenation preserves every original coefficient byte.
The [root storage specification](../README.md) defines the header, bit order,
monomial order, panel layout, and coordinate conventions.

From `tii-kernels/`, restore the input into a fresh directory:

```sh
python3 materialize.py tii-249 --output tii249-input
```

For the original jet implementation, also restore the exact panel and monomial map:

```sh
sage -python materialize.py tii-249 --output tii249-panel-input --panel
```

The [implementation](../../tii-keyrec-implementation/tii-249/README.md)
contains direct verification and recovery commands. The
[final key and later audits](../../tii-results/tii-249/README.md) are stored
separately. The holdout report's historical `curve_identity_certified` flag
describes that stage only.
