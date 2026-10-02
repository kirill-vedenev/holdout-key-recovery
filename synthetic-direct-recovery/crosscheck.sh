#!/usr/bin/env bash
# Recovers small keys of every family and rechecks them with SageMath.
set -euo pipefail
cd "$(dirname "$0")"
BIN=./build/synthetic_direct
make -s
rm -rf runs/crosscheck runs/crosscheck.jsonl
mkdir -p runs/crosscheck results
run() { "$BIN" run "$@" --dump runs/crosscheck --out runs/crosscheck.jsonl 2> /dev/null; }
run --family binary-goppa --q 2 --m 8 --n 100 --t 8 --seeds 1-4 --label cc-binary-goppa-q2
run --family wild-goppa --q 3 --m 5 --n 120 --t 4 --seeds 1-4 --label cc-wild-goppa-q3
run --family wild-goppa --q 4 --m 4 --n 120 --t 4 --seeds 1-4 --label cc-wild-goppa-q4
run --family wild-goppa --q 5 --m 3 --n 100 --t 4 --seeds 1-4 --label cc-wild-goppa-q5
run --family alternant --q 2 --m 8 --n 120 --r 10 --seeds 1-4 --label cc-alternant-q2
run --family alternant --q 2 --m 8 --n 120 --r 13 --seeds 1-4 --label cc-alternant-q2-cubic
run --family alternant --q 3 --m 5 --n 120 --r 8 --seeds 1-4 --label cc-alternant-q3
run --family alternant --q 4 --m 4 --n 120 --r 8 --seeds 1-4 --label cc-alternant-q4
run --family alternant --q 5 --m 3 --n 100 --r 6 --seeds 1-4 --label cc-alternant-q5
sage -python crosscheck.py runs/crosscheck results/crosscheck.json

# Build information for the cross-check.
SECTION=cross_check SRC=src/synthetic_direct.cpp
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
