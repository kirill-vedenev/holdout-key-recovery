#!/bin/sh
# Build src/libminor_native.so against the M4RI and M4RIE libraries shipped
# with SageMath. SAGE_PREFIX overrides the installation prefix.
set -eu
cd "$(dirname "$0")"
prefix=${SAGE_PREFIX:-$(sage -python -c 'from sage.env import SAGE_LOCAL; print(SAGE_LOCAL)')}
${CXX:-c++} -O3 -std=c++14 -fPIC -shared -Wno-deprecated-register -Wno-expansion-to-defined \
    -I"$prefix/include" src/native.cpp -L"$prefix/lib" -Wl,-rpath,"$prefix/lib" \
    -lm4rie -lm4ri -o src/libminor_native.so
