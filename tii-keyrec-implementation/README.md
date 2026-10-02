# TII key-recovery implementations

The computational sources used for the completed TII-173, TII-249, and TII-253
recoveries, with explicit commands for each stage. Original results were
packaged without repeating their computations. A new TII-253 direct finishing
was run on 2 October 2026 from the saved public jet prefix and holdout sample.

Each instance README starts with a short standalone Python check of its saved
key. That check establishes key validity. Optional audits of kernels, jets,
and geometry are documented at the end of the corresponding README.

| Instance | Holdout computation | Finishing |
|---|---|---|
| [TII-173](tii-173/README.md) | Exact binary panel elimination on two GPUs | Direct and minor finishing on the same order-296 jet; full key via direct finishing and deshortening |
| [TII-249](tii-249/README.md) | GPU sequence → CPU generator → GPU reconstruction | Direct and minor finishing on the same order-290 jet; full key via direct finishing and IKZ extension |
| [TII-253](tii-253/README.md) | GPU sequence → CPU generator → GPU reconstruction | Minor finishing from an order-82 jet; direct finishing after continuation to 198; both produce full keys |

TII-173 and TII-249 also contain the recorded `recovery/minor/` implementation.
Those runs reached equivalent retained supports; full-key extension was not
repeated separately. TII-253 has separately completed direct and minor finishings. See the
[finishing coverage](../TII-EXPERIMENTS.md#finishing-coverage) for exact scope.

Each instance is self-contained at the source level. `holdout/` contains the
native computation and public preparation; `recovery/` contains jet/finishing
code and independent verifiers. The larger holdouts retain the three-stage
GPU/CPU/GPU interface. A reproducer chooses how to execute or schedule those
stages. There are no server launchers, transfers, supervisors, download scripts,
or multi-server orchestration in this collection.

The commands reference the sibling [tii-kernels](../tii-kernels/README.md)
and [tii-results](../tii-results/README.md) folders. Use fresh output folders.
Stored recovered keys and evidence are used only by the verification commands,
not by the recovery commands. The latter start from public inputs and either
the saved holdout sample or a newly computed sample.

## Dependencies and source fidelity

The completed CPU recoveries used SageMath 10.7 with NumPy and a C++17 compiler.
TII-173's native formal-substitution helpers also link M4RI from Sage. The
independent key verifiers use only Python's standard library. Holdout GPU code
requires a CUDA development toolkit; the original runs used H200 hardware.
TII-249 and TII-253 additionally use cuFFT and the pinned, patched CADO-NFS
`lingen_b64` described in their holdout instructions.

[PROVENANCE.json](PROVENANCE.json) maps every copied source and archived
artifact to its original path and SHA-256. Computational source is unchanged
except for documented path relocation, a fresh output argument for TII-253
preparation, and removal of an obsolete module search path. TII-253's public
preparation helper contains six unchanged functions extracted from the original
module; its unused command-line entry point is omitted. Build files and
documentation are adapted for this layout. The kernel materializer is new
format-only tooling; reconstructed bytes are checked against original hashes.

MIT notices for the reused Saarinen code are retained under each holdout's
`third_party/saarinen/`. The CADO patch, pinned revision, and upstream license
are retained under `holdout/patches/`. No precompiled machine-specific binaries
are distributed. Syntax, archive integrity, and independent saved-key checks
are separate from rerunning recovery; this relocated layout has not undergone
a new GPU-to-key end-to-end recovery run. The additional TII-253 direct
finishing in `recovery/direct/` was executed from the existing public inputs;
its source is copied from that completed run.

`PACKAGING-CHECKS.json` records the integrity and saved-key checks performed
for this layout, including the saved-artifact audits now documented in the
per-instance READMEs. No recovery pipeline was rerun. `SHA256SUMS`
covers all distributed files under this root, excluding Finder metadata,
caches, and generated build/run files;
run `shasum -a 256 -c SHA256SUMS` from this directory to check them. The results
and kernel collections each have their own checksum file.
