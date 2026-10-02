#!/bin/sh
# Run from a public key through a newly computed kernel and verified recovery.
set -eu
root=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
public=${1:-"$root/instances/a/public.txt"}
output=${2:-"$root/runs/demo"}
route=${3:-both}
if [ -e "$output" ]; then
    echo "Choose a fresh output directory: $output already exists" >&2
    exit 1
fi
if [ ! -x "$root/build/holdout" ] || [ ! -x "$root/build/toy-goppa" ]; then
    echo "Build the two C++ programs with make first" >&2
    exit 1
fi
case "$route" in both|minor|direct) ;; *) echo "Route must be both, minor, or direct" >&2; exit 1 ;; esac
mkdir -p "$output"
"$root/build/toy-goppa" prepare "$public" "$output/holdout.txt"
echo "Computing and verifying a fresh holdout kernel"
"$root/build/holdout" "$output/holdout.txt" "$output/kernel.bin" "${TOY_MATRIX_GIB:-3}" > "$output/holdout-report.json"
"$root/build/toy-goppa" recover "$public" "$output/kernel.bin" "$output/recovery" "$route"
for name in minor direct; do
    if [ -f "$output/recovery/$name-key.txt" ]; then
        "$root/build/toy-goppa" verify "$public" "$output/recovery/$name-key.txt" > "$output/$name-verification.json"
    fi
done
echo "Completed: $output/recovery/recovery.json"
