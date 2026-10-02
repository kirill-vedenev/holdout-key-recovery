# Recorded recoveries

| Instance | Public code | Kernel dimension | Panel | Jet order | Minor rank | Kernel time | Recovery time |
|---|---|---:|---:|---:|---:|---:|---:|
| A | [63,27] | 5867 | 256 | 91 | 46 | 8.77 s | 3.29 s |
| B | [55,28] | 22223 | 256 | 83 | 42 | 279.83 s | 38.67 s |
| C | [64,28] | 7844 | 256 | 93 | 47 | 13.06 s | 2.38 s |

All three instances pass both finishing methods, independent public verification,
and comparison with the separately held generated key. Each successful run uses
one held-out position and one computed kernel. Recovery times include complete
panel verification and geometric certificates.

The recorded large-instance kernel uses an OpenMP-enabled M4RI build with a
16-thread limit. The two small-instance timings use the system M4RI build.
The driver defaults to the available CPU count; its thread limit is configurable.

Minor-only and direct-only checks on instance A reproduce the combined-mode keys.
These are fixed examples rather than a statistical sweep.

An additional OpenMP configuration check uses all 96 available threads. It
recomputes instance A in 2.12 seconds and produces a byte-identical kernel.
The driver defaults to the available CPU count and rejects a large allocation
when linked to single-threaded M4RI, unless explicitly overridden.
