# TII challenge results

Completed public-only recoveries of equivalent binary Goppa decoding keys.
Each instance contains the unchanged original public input, the recovered
key, an explicit extension-field definition, and saved verification evidence.
The recovered keys need not equal the challenge author's original keys.

Each instance README gives a short Python command to verify its saved key.
Optional intermediate audits are documented in the corresponding
[implementation README](../tii-keyrec-implementation/README.md).

| Instance | Original code | Field | Field-defining polynomial | Goppa degree | Result |
|---|---|---|---|---:|---|
| [TII-173](tii-173/README.md) | `[96,40]` | GF(128) | `z^7 + z + 1` (`0x83`) | 8 | Full equivalent key |
| [TII-249](tii-249/README.md) | `[235,107]` | GF(256) | `z^8 + z^4 + z^3 + z^2 + 1` (`0x11d`) | 16 | Full equivalent key |
| [TII-253](tii-253/README.md) | `[214,142]` | GF(256) | `z^8 + z^4 + z^3 + z^2 + 1` (`0x11d`) | 9 | Full equivalent key |

Both direct and minor finishings were tested on the retained TII-173 and
TII-249 codes; their archived minor results are in `evidence/minor-finishing/`.
TII-253 has both a minor-route key and an independently produced direct-route
key (`direct-recovered-sk.json`), added on 2 October 2026. The shared full-key
extension was not repeated for the two additional TII-173/TII-249 minor finishings. See
[finishing coverage](../TII-EXPERIMENTS.md#finishing-coverage).

## Files and conventions

- `pk_McEliece_N.txt`: original binary parity-check rows, followed by the
  field-modulus coefficient list. The public code is the binary kernel of
  this matrix, not its row space.
- `recovered-sk.json`: byte-for-byte copy of the completed equivalent key.
  `support[j]` is the field element at original public column `j`;
  `g` lists the monic Goppa polynomial coefficients in increasing degree.
- `field.json`: explicit field polynomial, its coefficient list and integer
  encoding, original parameters, and coordinate convention.
- `verification.json`: the saved independent full-key verification report.
- `timings.json`: recorded hardware, elapsed stage times, timing boundaries,
  and the hashes of their original measurement records.
- `evidence/`: saved jets, intermediate algebraic outputs, and verification
  certificates. These are outputs, not inputs to fresh recovery.
- `source.json`: original public-input provenance. Historical scope labels
  describe the input-acquisition stage, not the final recovery status.

An integer `a` encodes `sum_i ((a >> i) & 1) z^i` modulo the field polynomial.
Thus both `field_modulus` and `g` are low-degree-first coefficient lists,
but entries of `field_modulus` are binary and entries of `g` are encoded
extension-field elements. Support and multiplier arrays retain the original
public column order. No permutation needs to be inferred.

`dual_grs_multipliers` stores `nu`, with check rows
`nu_j * alpha_j^r`, `r=0..2t-1`, normalized by `nu[0]=1`.
`ambient_grs_multipliers` stores `lambda`, with generator rows
`lambda_j * alpha_j^r`, `r=0..n-2t-1`, where
`lambda_j = 1/(nu_j * product_{i!=j}(alpha_j + alpha_i))`.
When present, `goppa_parity_multipliers` stores `1/g(alpha_j)` for the
degree-`t` Goppa check; this is a different convention from `nu`.
The pair `(support,g)` and its field definition suffice to describe the key.

The standard-library verifiers in the corresponding
[implementation folders](../tii-keyrec-implementation/README.md) check exact
equality of the reconstructed binary parity-check row space with the original
public matrix, distinct support, support avoidance, and irreducibility of `g`.
Each instance README gives its verification command.

The original holdout samples and storage specification are in
[tii-kernels](../tii-kernels/README.md). Checksums bind the distributed files;
`PROVENANCE.json` in the implementation root records their archived source
paths and hashes. Original archived runs were preserved; the additional
TII-253 direct finishing was executed on 2 October 2026 without recomputing
the holdout sample.
Some unchanged historical reports contain the original local input paths;
these are provenance labels, not runtime dependencies of this distribution.
From this directory, `shasum -a 256 -c SHA256SUMS` checks every distributed file.
