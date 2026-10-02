#!/usr/bin/env bash
# Runs all parameter sets of this collection and writes results/.
#
# Environment variables:
#   JOBS    parallel processes (default: number of CPUs)
#   SCALE   multiplies the number of keys of every set (default 1)
#
# Each job is one process; jobs run in parallel.
set -euo pipefail
cd "$(dirname "$0")"

JOBS=${JOBS:-$(getconf _NPROCESSORS_ONLN)}
SCALE=${SCALE:-1}
BIN=./build/synthetic_direct

make -s
"$BIN" self-test

rm -rf runs/parts
mkdir -p runs/parts results
JOBLIST=runs/jobs.txt
: > "$JOBLIST"

# label family q m n t r keys keys-per-job [extra options]
SETS=(
  "binary-goppa-q2-m8-n100 binary-goppa 2 8 100 8 0 100 5"
  "binary-goppa-q2-m10-n400 binary-goppa 2 10 400 20 0 20 1"
  "binary-goppa-mceliece348864-weis binary-goppa 2 12 922 64 0 20 1"
  "binary-goppa-mceliece6688128-weis binary-goppa 2 13 1888 128 0 20 1"
  "control-binary-goppa-adjacent binary-goppa 2 10 400 20 0 20 1 --adjacent"
  "wild-goppa-q3-m5-n200 wild-goppa 3 5 200 6 0 100 5"
  "wild-goppa-q3-m7-n800 wild-goppa 3 7 800 20 0 20 1"
  "wild-goppa-q4-m4-n200 wild-goppa 4 4 200 6 0 100 5"
  "wild-goppa-q4-m5-n600 wild-goppa 4 5 600 12 0 20 1"
  "wild-goppa-q5-m3-n100 wild-goppa 5 3 100 4 0 100 5"
  "wild-goppa-q5-m4-n500 wild-goppa 5 4 500 8 0 20 1"
  "alternant-q2-m8-n200 alternant 2 8 200 0 16 100 5"
  "alternant-q2-m8-n200-cubic alternant 2 8 200 0 22 100 5"
  "alternant-q2-m11-n1000 alternant 2 11 1000 0 40 20 1"
  "alternant-q3-m5-n200 alternant 3 5 200 0 12 100 5"
  "alternant-q3-m7-n800 alternant 3 7 800 0 40 20 1"
  "alternant-q4-m4-n200 alternant 4 4 200 0 12 100 5"
  "alternant-q5-m3-n100 alternant 5 3 100 0 8 100 5"
)

for spec in "${SETS[@]}"; do
  read -r label family q m n t r keys chunk extra <<< "$spec"
  keys=$((keys * SCALE))
  lo=1
  while [ "$lo" -le "$keys" ]; do
    hi=$((lo + chunk - 1))
    [ "$hi" -gt "$keys" ] && hi=$keys
    echo "$BIN run --family $family --q $q --m $m --n $n --t $t --r $r --seeds $lo-$hi ${extra:-} --label $label --out runs/parts/$label--$lo.jsonl" >> "$JOBLIST"
    lo=$((hi + 1))
  done
done

echo "running $(wc -l < "$JOBLIST") jobs with $JOBS processes"
start=$(date +%s)
xargs -P "$JOBS" -I{} sh -c '{} 2>/dev/null' < "$JOBLIST"
echo "jobs finished in $(( $(date +%s) - start )) s"

cat runs/parts/*.jsonl > results/keys.jsonl
gzip -9 -f results/keys.jsonl

if command -v sage > /dev/null; then
  ./crosscheck.sh
fi

# Build information for the batch.
SECTION=batch SRC=src/synthetic_direct.cpp
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
python3 summarize.py results
