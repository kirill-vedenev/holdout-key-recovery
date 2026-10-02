# Independent SageMath check of instances dumped by `product_rank --dump DIR`.
#
# For each dumped key (field modulus, Goppa polynomial, support), this script
# rebuilds the binary Goppa code from its own parity-check matrix, computes the
# curve polynomials f_c = (sum_j c_j Pi/(Z - alpha_j)) / G^2 in coefficient
# form, and compares with the C++ results:
#   * dim W = rank of all products f_a f_b, a <= b (no early stop);
#   * dim W <= 2D - t + 1 and G^2 | (f_a f_b)' for every product;
#   * membership of (Z - alpha)^{2D}, (Z - alpha)^{2D-2} and (Z - alpha)^{2D-1}
#     in W at the dumped positions;
#   * rank of the cross products f_a f_b, a < b, equals dim W - k.
#
# Usage: sage -python crosscheck.py DIR [OUT.json]

import glob
import json
import os
import sys

from sage.all import GF, PolynomialRing, matrix, prod, vector


def read_dump(path):
    d = {}
    with open(path) as f:
        for line in f:
            parts = line.split()
            if parts:
                d[parts[0]] = [int(x) for x in parts[1:]]
    return d


def check(path):
    d = read_dump(path)
    m, mod, t, n = d['m'][0], d['modulus'][0], d['t'][0], d['n'][0]
    modpoly = GF(2)['x']([(mod >> i) & 1 for i in range(m + 1)])
    F = GF(2**m, 'z', modulus=modpoly)
    P = PolynomialRing(F, 'Z')
    Z = P.gen()
    G = P([F.from_integer(c) for c in d['G']])
    alpha = [F.from_integer(a) for a in d['alpha']]
    assert len(alpha) == n and len(set(alpha)) == n
    assert all(G(a) != 0 for a in alpha)
    D = n - 2 * t - 1
    bound = 2 * D - t + 1

    rows = []
    for i in range(t):
        vals = [a**i / G(a) for a in alpha]
        for b in range(m):
            rows.append([(v.to_integer() >> b) & 1 for v in vals])
    code = matrix(GF(2), rows).right_kernel()
    k = code.dimension()

    Pi = prod(Z - a for a in alpha)
    G2 = G**2
    fs = []
    for c in code.basis():
        fc = sum(Pi // (Z - alpha[j]) for j in range(n) if c[j] != 0)
        q, r = fc.quo_rem(G2)
        assert r == 0, "G^2 does not divide f_c"
        assert q.degree() <= D
        fs.append(q)

    def vec(p):
        coeffs = p.list()
        return coeffs + [F(0)] * (2 * D + 1 - len(coeffs))

    prods = [fs[a] * fs[b] for a in range(k) for b in range(a, k)]
    in_S = all((p.derivative() % G2) == 0 for p in prods)
    M = matrix(F, [vec(p) for p in prods])
    rank = M.rank()
    cross = [vec(fs[a] * fs[b]) for a in range(k) for b in range(a + 1, k)]
    cross_rank = matrix(F, cross).rank() if cross else 0
    space = M.row_space()

    def member(e, a):
        return int(vector(F, vec((Z - a)**e)) in space)

    positions = d.get('positions', [])
    r_in = [member(2 * D, alpha[j]) for j in positions]
    s_in = [member(2 * D - 2, alpha[j]) for j in positions]
    neg_in = [member(2 * D - 1, alpha[j]) for j in positions]

    agree = (k == d['k'][0] and rank == d['rank'][0] and bound == d['bound'][0]
             and r_in == d.get('R_in_W', []) and s_in == d.get('S_in_W', [])
             and neg_in == d.get('neg_in_W', []))
    return {
        'file': os.path.basename(path), 'm': int(m), 't': int(t), 'n': int(n), 'k': int(k),
        'D': int(D), 'bound': int(bound), 'sage_rank': int(rank), 'cpp_rank': d['rank'][0],
        'cross_rank': int(cross_rank), 'cross_rank_is_rank_minus_k': bool(cross_rank == rank - k),
        'products_in_S': bool(in_S), 'rank_le_bound': bool(rank <= bound),
        'positions': len(positions), 'agree': bool(agree),
    }


def main():
    if len(sys.argv) < 2:
        print("usage: sage -python crosscheck.py DIR [OUT.json]")
        sys.exit(2)
    files = sorted(glob.glob(os.path.join(sys.argv[1], '*.txt')))
    results = [check(p) for p in files]
    ok = all(r['agree'] and r['products_in_S'] and r['rank_le_bound']
             and r['cross_rank_is_rank_minus_k'] for r in results)
    summary = {
        'instances': len(results),
        'deficient': int(sum(r['sage_rank'] < r['bound'] for r in results)),
        'all_agree': ok,
        'results': results,
    }
    if len(sys.argv) > 2:
        with open(sys.argv[2], 'w') as f:
            json.dump(summary, f, indent=1)
    print("crosscheck: %d instances (%d deficient), %s"
          % (summary['instances'], summary['deficient'], 'ok' if ok else 'MISMATCH'))
    sys.exit(0 if ok else 1)


main()
