# TII-249 holdout: three computational stages

Run these commands from the parent `tii-249/` directory. These are the
actual numerical implementations used, with direct stage commands replacing
the server harness. Stages exchange ordinary files and can be scheduled
independently. No kernel computation was rerun while preparing this package.

## Build the GPU worker and CPU verifier

Requirements: C++17, GNU Make, a CUDA development toolkit with cuFFT, and an
H200-class GPU with sufficient memory. The saved computation used one H200.
Build for the target device with `CUDA_ARCH` (90 for H200):

```sh
make -C holdout cpu gpu CUDA_ARCH=90
```

`src/square_cufft_worker.cu` is byte-identical to the saved production source,
SHA-256 `9b17db24cc19985dd2fc4a935bc732ad2de968bfb51d11e95956638d21826307`.
It includes the preserved Saarinen worker under `third_party/saarinen/`.

## Build the CPU polynomial generator

Use CADO-NFS source at revision
`2d98fe176c11342cfa5d6c0d795a992cbe433c11` with its build dependencies installed
(CMake, a GNU C/C++ toolchain with OpenMP, and GMP development files).
`CADO_SRC` below must point to that source checkout. Apply the included patch
once to a clean checkout; `patches/cado.json` records its hash and purpose.
The patch parallelizes coefficient updates and avoids nested OpenMP teams.

```sh
CADO_SRC=/path/to/cado-nfs
patch -d "$CADO_SRC" -p1 < holdout/patches/cado-coefficient-parallel.patch
cmake -S "$CADO_SRC" -B "$CADO_SRC/build-tii" \
  -DCMAKE_BUILD_TYPE=Release -DWITH_MPI=OFF
cmake --build "$CADO_SRC/build-tii" --target lingen_b64 -j
LINGEN="$CADO_SRC/build-tii/linalg/bwc/lingen_b64"
```

The preserved CADO binary used by the completed run had SHA-256
`72af4da948a42e56be0717cdcaffa43f741990d54116695e56a3069992853630`.
A build on a different toolchain need not have the same executable hash.
The upstream license is retained as `patches/CADO-COPYING`.

## Stage 1: GPU sequence

Use the same seeds and block width for every stage. `request.txt` specifies
eight 64-bit words, hence a 512-bit block. The compressor seed is
`6073466052459188813`, and the starting panel seed is `24920260930`.

```sh
WORK="$PWD/runs/holdout"
mkdir -p "$WORK"
holdout/build/square-worker segment ../../tii-kernels/tii-249/request.txt \
  - 8705 "$WORK/sequence" 2 24920260930 6073466052459188813
```

The sequence is `S_i = Z^T A^i Y`, beginning with **S_0**. Each term is a
512-by-512 binary matrix, row-major with eight little-endian words per row,
32,768 bytes per term. `sequence.raw.bin` has no header and can be passed
directly to the generator. `state.bin` is an optional numerical continuation
state; no service or transfer harness is required.

The old TII-249 production stream initially began at S_1 and was corrected
before its accepted CPU calculation. The included worker emits S_0 directly;
the obsolete strip/prepend adapters are not part of this recipe.

## Stage 2: CPU generator

The settings below reproduce the recorded mathematical parameters and thread
configuration. Adapt available CPU resources as needed. `max_ram=128` is the
CADO memory setting in GB. Allow at least 5 GB of disk for these holdout
artifacts, plus the space needed for recovery outputs.

```sh
export OMP_NUM_THREADS=32
export OMP_DYNAMIC=FALSE
export OMP_PLACES=cores
export OMP_PROC_BIND=spread
export OMP_MAX_ACTIVE_LEVELS=1
"$LINGEN" prime=2 m=512 n=512 \
  "afile=$WORK/sequence/sequence.raw.bin" input_length=8705 \
  "ffile=$WORK/square.gen" rhs=none thr=4x8 max_ram=128 \
  tuning_thresholds=recursive:100000,notiming:1
```

`square.gen` is CADO's single-file polynomial matrix. The GPU reconstruction
must use **transpose=1**, matching this orientation. Preserve successful stage
outputs yourself if running the stages at different times or on different
machines. The original CPU basecase does not support mid-iteration restart.

## Stage 3: GPU reconstruction, then independent verification

```sh
holdout/build/square-worker recover ../../tii-kernels/tii-249/request.txt \
  "$WORK/square.gen" "$WORK/reconstruction" 2 24920260930 6073466052459188813 1 280
holdout/build/holdout-cpu verify ../../tii-kernels/tii-249/operator.txt \
  "$WORK/reconstruction/kernel.thk1" "$WORK/verification.json" --threads 32
cp ../../tii-kernels/tii-249/public.json "$WORK/reconstruction/public.json"
```

The worker checks the literal uncompressed constraints; the CPU verifier
independently checks **all original Hasse orders** and rank 280. Require the
verifier to pass before treating the sample as accepted. The output is
`kernel.thk1` plus the panel, monomial map, and index under `polynomials/`.
The [storage specification](../../../tii-kernels/README.md) defines these
encodings. The [parent instructions](../README.md) continue with jet recovery.
