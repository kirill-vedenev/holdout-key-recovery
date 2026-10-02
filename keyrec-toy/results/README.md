# Recorded notebook executions

All four notebooks pass with SageMath 10.7. Every code cell completes without
an error; the saved outputs and JSON records contain the checks below.

| Example | Kernel dimension / used | Branches | Continuation rank | Jet order | Minor rank | Recovered positions |
|---|---:|---:|---:|---:|---:|---:|
| GRS, characteristic two | 34 / 34 | 1 | 4 | 13 | 7 | 24 |
| Alternant, characteristic two | 361 / 48 | 2 | 18 | 77 | 39 | 60 |
| GRS, odd characteristic | 4 / 4 | 1 | 2 | 14 | — | 17 |
| Alternant, odd characteristic | 21 / 21 | 2 | 13 | 57 | — | 45 |

The per-notebook JSON files contain the parameters, ranks, verification results,
and internal timings. `execution-summary.json` records execution times including
kernel startup, the runtime platform, and SHA-256 hashes of the notebooks and
shared implementation. The recorded execution times are approximately 4.0, 49.0, 3.7, and 12.6 seconds
respectively. Platform details are included in the execution summary.

All cases pass the independent secret-based coherence and support audits.
Both characteristic-two cases pass exact minor-GRS equality and agreement of
the two finishing routes. Both alternant cases pass exact public-subfield-code
equality. Every notebook rejects a deliberately corrupted multiplier.

These are four fixed examples, not a statistical sweep. Rerunning a notebook
updates its result record; rerun `execute.py` to refresh all notebook hashes in
the combined summary.
