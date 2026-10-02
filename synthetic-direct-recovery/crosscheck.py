# Independent SageMath check of keys recovered by `synthetic_direct --dump DIR`.
#
# For each dump, with K = GF(p^e) built from the dumped modulus and
# F_q ⊂ K the subfield of order q, the script checks:
#   * the public generator Y has rank k and entries in F_q;
#   * C ⊆ GRS_{D+1}(alpha', lambda'), i.e. H Y^T = 0 for a parity-check
#     matrix H of the recovered GRS code;
#   * dim_{F_q}(GRS_{D+1}(alpha', lambda') ∩ F_q^n) = k, computed as
#     n - rank Tr_{K/F_q}(gamma H) (Delsarte), so the subfield subcode is C;
#   * the recovered support is a Moebius image of phi^sigma(alpha), where
#     sigma is the branch of the jet (audit with the secret support).
#
# Usage: sage -python crosscheck.py DIR [OUT.json]

import glob
import json
import os
import sys

from sage.all import GF, PolynomialRing, matrix, prod


def read_dump(path):
    d = {"Y": []}
    with open(path) as f:
        for line in f:
            parts = line.split()
            if not parts:
                continue
            vals = [int(x) for x in parts[1:]]
            if parts[0] == "Y":
                d["Y"].append(vals)
            else:
                d[parts[0]] = vals
    return d


def check(path):
    d = read_dump(path)
    p, e, q, m = d["p"][0], d["e"][0], d["q"][0], d["m"][0]
    n, k, D, sigma = d["n"][0], d["k"][0], d["D"][0], d["branch"][0]
    modulus = PolynomialRing(GF(p), "x")(d["modulus"])
    K = GF(p**e, "z", modulus=modulus)
    elt = K.from_integer
    Y = matrix(K, [[elt(x) for x in row] for row in d["Y"]])
    a = [elt(x) for x in d["alpha_recovered"]]
    lam = [elt(x) for x in d["lambda_recovered"]]
    secret = [elt(x) for x in d["alpha_secret"]]

    in_fq = all(x**q == x for x in Y.list())
    rank_ok = Y.rank() == k
    # Parity check of GRS_{D+1}(a, lam): rows a_j^l / (lam_j Pi'(a_j)).
    dpi = [prod(a[j] - a[l] for l in range(n) if l != j) for j in range(n)]
    w = [1 / (lam[j] * dpi[j]) for j in range(n)]
    H = matrix(K, [[w[j] * a[j]**l for j in range(n)] for l in range(n - D - 1)])
    contains = (H * Y.transpose()).is_zero()
    # Trace code of the dual over F_q.
    g = K.multiplicative_generator()
    basis = [g**i for i in range(m)]

    def trace(x):
        return sum(x**(q**i) for i in range(m))

    T = matrix(K, [[trace(b * H[l, j]) for j in range(n)] for b in basis for l in range(n - D - 1)])
    sub_dim = n - T.rank()
    # Moebius audit through cross-ratios.
    conj = [x**(q**sigma) for x in secret]

    def cr(z, j):
        return (z[2] - z[0]) * (z[j] - z[1]) / ((z[2] - z[1]) * (z[j] - z[0]))

    moebius = all(cr(conj, j) == cr(a, j) for j in range(3, n))
    ok = in_fq and rank_ok and contains and sub_dim == k and moebius and len(set(a)) == n
    return {"file": os.path.basename(path), "q": q, "m": m, "n": n, "k": k, "D": D,
            "entries_in_Fq": bool(in_fq), "rank_k": bool(rank_ok), "contains": bool(contains),
            "subfield_subcode_dim": int(sub_dim), "moebius_audit": bool(moebius), "ok": bool(ok)}


def main():
    if len(sys.argv) < 2:
        print("usage: sage -python crosscheck.py DIR [OUT.json]")
        sys.exit(2)
    files = sorted(glob.glob(os.path.join(sys.argv[1], "*.txt")))
    results = [check(f) for f in files]
    summary = {"keys": len(results), "all_ok": all(r["ok"] for r in results), "results": results}
    if len(sys.argv) > 2:
        with open(sys.argv[2], "w") as f:
            json.dump(summary, f, indent=1)
    print("crosscheck: %d recovered keys, %s" % (len(results), "ok" if summary["all_ok"] else "MISMATCH"))
    sys.exit(0 if summary["all_ok"] else 1)


main()
