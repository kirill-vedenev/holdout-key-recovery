# TII-253 recovered key

All **214 original support positions** and a monic irreducible degree-9
Goppa polynomial were recovered for the public `[214,142]` binary code.
The original computation extended an order-82 jet and restored 96
coordinates removed by shortening. This is an equivalent decoding key;
the challenge author's original secret is not required.

A direct finishing was added on 2 October 2026. It extends the shared public
jet from order 82 to 198 and independently produces a full equivalent key.
`recovered-sk.json` remains the original minor-route key;
`direct-recovered-sk.json` contains the new direct-route key.

- [Original public key](pk_McEliece_253.txt)
- [Recovered secret key](recovered-sk.json), copied unchanged from the completed run
- [Direct-finishing secret key](direct-recovered-sk.json) and [public-key verification](direct-verification.json)
- [Direct-finishing run results](evidence/direct-finishing/result-summary.json)
- [Extension field](field.json): `z^8 + z^4 + z^3 + z^2 + 1`, integer `0x11d`
- [Saved independent verification](verification.json)
- [Recorded hardware and timings](timings.json)
- [Intermediate outputs and certificates](evidence/)
- [Holdout sample](../../tii-kernels/tii-253/README.md)
- [Implementation and full run instructions](../../tii-keyrec-implementation/tii-253/README.md)

The arrays `support`, `dual_grs_multipliers`, and `ambient_grs_multipliers`
are in original column order. `g` is low degree first. The explicit field
polynomial is also retained in `recovered-sk.json:field_modulus` and as the
last line of the original public-key file. Full conventions are described in
the [results root](../README.md).

From this instance directory, independently verify the saved key with Python 3:

```sh
python3 ../../tii-keyrec-implementation/tii-253/recovery/src/independent_verify.py recovered-sk.json \
  --public pk_McEliece_253.txt --report tii253-verification.json
```

This checks the saved key; it does not recompute a kernel or run recovery.
The archived reports remain unchanged. Source-file hashes and any packaging
edits are recorded in [PROVENANCE.json](../../tii-keyrec-implementation/PROVENANCE.json).
