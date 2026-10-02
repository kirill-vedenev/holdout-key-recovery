#!/bin/sh
# Build an OpenMP-enabled M4RI archive with GCC on Linux.
set -eu
root=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
prefix=${1:-"$root/deps/m4ri-omp"}
revision=35f06e132363d12592cdd94c98473d36e736cef9
mkdir -p "$prefix/src"
if [ ! -d "$prefix/src/m4ri-master" ]; then
    curl --fail --location --silent --show-error \
        "https://github.com/malb/m4ri/archive/$revision.zip" -o "$prefix/src/m4ri.zip"
    python3 - "$prefix/src" <<'PY'
from pathlib import Path
import sys, zipfile
directory = Path(sys.argv[1])
with zipfile.ZipFile(directory / 'm4ri.zip') as archive:
    names = {name.split('/')[0] for name in archive.namelist()}
    if len(names) != 1:
        raise ValueError('unexpected source archive layout')
    archive.extractall(directory)
(directory / names.pop()).rename(directory / 'm4ri-master')
PY
fi
config=
for candidate in "${M4RI_PREFIX:-/usr}/include/m4ri/m4ri_config.h" /usr/local/include/m4ri/m4ri_config.h; do
    if [ -f "$candidate" ]; then config=$candidate; break; fi
done
if [ -z "$config" ]; then
    echo "Install M4RI development headers or set M4RI_PREFIX" >&2
    exit 1
fi
cp "$config" "$prefix/src/m4ri-master/m4ri/m4ri_config.h"
python3 - "$prefix/src/m4ri-master/m4ri/m4ri_config.h" <<'PY'
from pathlib import Path
import re, sys
p = Path(sys.argv[1])
text, count = re.subn(r'(#define\s+__M4RI_HAVE_OPENMP\s+)\d+', r'\g<1>1', p.read_text())
if count != 1:
    raise ValueError('M4RI configuration lacks its OpenMP switch')
p.write_text(text)
PY
mkdir -p "$prefix/objects"
for source in "$prefix"/src/m4ri-master/m4ri/*.c; do
    case "$source" in */io.c) continue ;; esac
    object="$prefix/objects/$(basename "$source" .c).o"
    "${CC:-gcc}" -O3 -march=native -fopenmp -std=gnu11 \
        -I"$prefix/src/m4ri-master" -I"$prefix/src/m4ri-master/m4ri" \
        -c "$source" -o "$object"
done
ar rcs "$prefix/libm4ri_omp.a" "$prefix"/objects/*.o
echo "Built $prefix/libm4ri_omp.a"
echo "Build the programs with: make M4RI_OMP=$prefix"
