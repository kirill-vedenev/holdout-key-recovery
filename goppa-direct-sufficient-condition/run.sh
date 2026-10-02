#!/usr/bin/env bash
# Runs the experiments of this collection and writes results/.
#
# Environment variables:
#   JOBS        parallel processes (default: number of CPUs)
#   MCE_SEEDS   keys per Classic McEliece preset (default 1000)
#   TOY_SEEDS   keys per toy parameter set (default 1000)
#   EDGE_SEEDS  keys per parameter set near the counting threshold (default 20)
#   MOD_SEEDS   keys per moderate shortened parameter set (default 500)
#   FULL_SEEDS  keys per full-length parameter set (default 1000)
#   CHUNK       keys per job (default 50)
#
# Each job runs one process with one thread; jobs run in parallel.
set -euo pipefail
cd "$(dirname "$0")"

JOBS=${JOBS:-$(getconf _NPROCESSORS_ONLN)}
MCE_SEEDS=${MCE_SEEDS:-1000}
TOY_SEEDS=${TOY_SEEDS:-1000}
EDGE_SEEDS=${EDGE_SEEDS:-20}
MOD_SEEDS=${MOD_SEEDS:-500}
FULL_SEEDS=${FULL_SEEDS:-1000}
CHUNK=${CHUNK:-50}
BIN=./build/product_rank

make -s
"$BIN" self-test

rm -rf runs/parts
mkdir -p runs/parts results
JOBLIST=runs/jobs.txt
: > "$JOBLIST"

# add GROUP NAME SEEDS ARGS... : split the seeds into jobs of CHUNK keys.
add() {
  local group=$1 name=$2 total=$3
  shift 3
  local lo=1 hi
  while [ "$lo" -le "$total" ]; do
    hi=$((lo + CHUNK - 1))
    [ "$hi" -gt "$total" ] && hi=$total
    echo "$BIN $* --seeds $lo-$hi --threads 1 --out runs/parts/$group--$name--$lo.jsonl" >> "$JOBLIST"
    lo=$((hi + 1))
  done
}

# 1. Classic McEliece keys, shortened to the distinguisher dimension k_l.
for p in mceliece348864 mceliece460896 mceliece6688128 mceliece6960119 mceliece8192128; do
  for s in weis giajs; do
    add mceliece "$p-$s" "$MCE_SEEDS" preset "$p-$s" --positions 8
  done
done

# 2. The counting threshold C(k+1,2) >= 2D-t+1 at Classic McEliece (m, t).
for k in $(seq 47 56); do
  add edge "m12-t64-k$k" "$EDGE_SEEDS" run --m 12 --t 64 --k "$k" --label edge-m12-t64
done
for k in $(seq 71 80); do
  add edge "m13-t128-k$k" "$EDGE_SEEDS" run --m 13 --t 128 --k "$k" --label edge-m13-t128
done

# 3. Small codes: every product is used, so dim W <= 2D-t+1 is tested directly.
#    Format: m t k_min k_max (k_min is a few below the counting threshold).
for spec in "6 3 5 20" "7 4 7 22" "8 5 9 26" "8 8 12 30" "9 8 13 30" "10 10 16 32" "11 16 22 36"; do
  read -r m t k0 k1 <<< "$spec"
  for k in $(seq "$k0" "$k1"); do
    add toy "m$m-t$t-k$k" "$TOY_SEEDS" run --m "$m" --t "$t" --k "$k" --full --label "toy-m$m-t$t"
  done
done

# 4. Moderate shortened codes, between the small codes and Classic McEliece.
#    Format: m t followed by the k values (from the counting threshold upward).
for spec in "9 16 22 23 24 25 27 30 34 42 57 82 122" \
            "10 20 26 27 28 29 31 34 38 46 61 86 126" \
            "10 40 37 38 39 40 42 45 49 57 72 97 137" \
            "11 32 35 36 37 38 40 43 47 55 70 95 135" \
            "12 32 37 38 39 40 42 45 49 57 72 97 137" \
            "13 64 54 55 56 57 59 62 66 74 89 114 154"; do
  set -- $spec
  m=$1 t=$2
  shift 2
  for k in "$@"; do
    add moderate "m$m-t$t-k$k" "$MOD_SEEDS" run --m "$m" --t "$t" --k "$k" --positions 4 --diagnose --label "moderate-m$m-t$t"
  done
done

# 5. Unshortened codes with the TII challenge parameters. Here 2D+1 exceeds
#    the field size, so the products are computed in coefficient form.
for spec in "tii-173 7 8 96" "tii-249 8 16 235" "tii-253 8 9 214"; do
  read -r name m t n <<< "$spec"
  add full "$name" "$FULL_SEEDS" run --m "$m" --t "$t" --n "$n" --positions 8 --diagnose --label "full-$name"
done

echo "running $(wc -l < "$JOBLIST") jobs with $JOBS processes"
start=$(date +%s)
xargs -P "$JOBS" -I{} sh -c '{} 2>/dev/null' < "$JOBLIST"
echo "jobs finished in $(( $(date +%s) - start )) s"

for g in mceliece edge toy moderate full; do
  cat runs/parts/"$g"--*.jsonl > "results/$g.jsonl"
  gzip -9 -f "results/$g.jsonl"
done

# 4. Independent SageMath check of small dumped instances (needs SageMath).
if command -v sage > /dev/null; then
  ./crosscheck.sh
fi

# 5. Build information and summary tables.
# Build information for the batch.
SECTION=batch SRC=src/product_rank.cpp
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
