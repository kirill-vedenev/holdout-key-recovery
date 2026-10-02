#!/bin/sh
set -eu
root=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
public=${1:-"$root/instances/a/public.txt"}
output=${2:-"$root/runs/demo"}
route=${3:-both}
OMP_NUM_THREADS=${WILD_THREADS:-${OMP_NUM_THREADS:-$(getconf _NPROCESSORS_ONLN)}}
export OMP_NUM_THREADS
if [ -e "$output" ]; then
    echo "Choose a fresh output directory: $output already exists" >&2
    exit 1
fi
if [ ! -x "$root/build/holdout" ] || [ ! -x "$root/build/wild-goppa" ]; then
    echo "Build the C++ programs with make first" >&2
    exit 1
fi
case "$route" in both|minor|direct) ;; *) echo "Route must be both, minor, or direct" >&2; exit 1 ;; esac
mkdir -p "$output"
echo "Computing one public holdout kernel over GF(4); thread limit: $OMP_NUM_THREADS"
"$root/build/holdout" "$public" "$output/kernel.bin" "${WILD_PANEL_SIZE:-256}" "${WILD_MATRIX_GIB:-16}" > "$output/holdout.json"
"$root/build/wild-goppa" recover "$public" "$output/kernel.bin" "$output/recovery" "$route"
for name in minor direct; do
    if [ -f "$output/recovery/$name-key.txt" ]; then
        "$root/build/wild-goppa" verify "$public" "$output/recovery/$name-key.txt" > "$output/$name-verification.json"
    fi
done
echo "Completed: $output/recovery/recovery.json"
