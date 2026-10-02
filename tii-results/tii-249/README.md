# TII-249 recovered key

All **235 original support positions** and a monic irreducible degree-16
Goppa polynomial were recovered for the public `[235,107]` binary code.
The original computation extended an order-290 jet and restored 57
coordinates removed by shortening. This is an equivalent decoding key;
the challenge author's original secret is not required.

- [Original public key](pk_McEliece_249.txt)
- [Recovered secret key](recovered-sk.json), copied unchanged from the completed run
- [Extension field](field.json): `z^8 + z^4 + z^3 + z^2 + 1`, integer `0x11d`
- [Saved independent verification](verification.json)
- [Recorded hardware and timings](timings.json)
- [Intermediate outputs and certificates](evidence/)
- [Archived minor-finishing results and support comparison](evidence/minor-finishing/)
- [Holdout sample](../../tii-kernels/tii-249/README.md)
- [Implementation and full run instructions](../../tii-keyrec-implementation/tii-249/README.md)

The arrays `support`, `dual_grs_multipliers`, and `ambient_grs_multipliers`
are in original column order. `g` is low degree first. The explicit field
polynomial is also retained in `recovered-sk.json:field_modulus` and as the
last line of the original public-key file. Full conventions are described in
the [results root](../README.md).

From this instance directory, independently verify the saved key with Python 3:

```sh
python3 ../../tii-keyrec-implementation/tii-249/recovery/full-key/independent_verify.py recovered-sk.json \
  --public pk_McEliece_249.txt --report tii249-verification.json
```

This checks the saved key; it does not recompute a kernel or run recovery.
The archived reports remain unchanged. Source-file hashes and any packaging
edits are recorded in [PROVENANCE.json](../../tii-keyrec-implementation/PROVENANCE.json).

Both finishings were tested through retained-support recovery on the same
saved jet. The minor result matches the direct key on retained positions up
to a Möbius transformation; full-key extension was not repeated for that route.
