#!/usr/bin/env bash
# Dumps small keys with their C++ results and rechecks them in SageMath.
set -euo pipefail
cd "$(dirname "$0")"
BIN=./build/product_rank
make -s
rm -rf runs/crosscheck runs/crosscheck.jsonl
mkdir -p runs/crosscheck results
for spec in "6 3 8-12 1-6" "7 4 10-13 1-5" "8 8 15-17 1-4"; do
  read -r m t ks seeds <<< "$spec"
  "$BIN" run --m "$m" --t "$t" --k "$ks" --seeds "$seeds" --full --positions all \
    --cross-rank --dump runs/crosscheck --label "cc-m$m-t$t" --out runs/crosscheck.jsonl 2> /dev/null
done
# Unshortened keys with TII challenge parameters (coefficient form).
for spec in "7 8 96 1-3" "8 9 214 1-2"; do
  read -r m t n seeds <<< "$spec"
  "$BIN" run --m "$m" --t "$t" --n "$n" --seeds "$seeds" --positions 16 \
    --dump runs/crosscheck --label "cc-full-m$m-t$t" --out runs/crosscheck.jsonl 2> /dev/null
done
sage -python crosscheck.py runs/crosscheck results/crosscheck.json

# Build information for the cross-check.
SECTION=cross_check SRC=src/product_rank.cpp
python3 - "$SECTION" "$SRC" <<'PY'
import hashlib, json, os, platform, subprocess, sys
section, src = sys.argv[1], sys.argv[2]
path = "results/build-info.json"
info = json.load(open(path)) if os.path.exists(path) else {}
info.setdefault("source_sha256", {})[src] = hashlib.sha256(open(src, "rb").read()).hexdigest()
cxx = os.environ.get("CXX", "c++")
info[section] = {
    "platform": platform.system() + " " + platform.machine(),
    "compiler": subprocess.run([cxx, "--version"], capture_output=True, text=True).stdout.splitlines()[0],
    "flags": os.environ.get("CXXFLAGS", "Makefile default"),
}
if section == "cross_check":
    info[section]["sage"] = subprocess.run(["sage", "--version"], capture_output=True, text=True).stdout.strip()
with open(path, "w") as f:
    json.dump(info, f, indent=2)
    f.write("\n")
PY
