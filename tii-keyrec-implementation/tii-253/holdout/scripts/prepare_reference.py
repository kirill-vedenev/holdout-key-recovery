"""Unchanged public-algebra helpers used by the TII-253 preparation stage."""
import json
import math


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":")).encode()



def rref(rows, n):
    a = list(rows)
    pivots = []
    for c in range(n):
        found = next((i for i in range(len(pivots), len(a)) if (a[i] >> c) & 1), None)
        if found is None:
            continue
        p = len(pivots)
        a[p], a[found] = a[found], a[p]
        for i in range(len(a)):
            if i != p and (a[i] >> c) & 1:
                a[i] ^= a[p]
        pivots.append(c)
    return a[: len(pivots)], pivots



def nullspace(rows, n):
    a, pivots = rref(rows, n)
    pivot_set = set(pivots)
    free = [j for j in range(n) if j not in pivot_set]
    basis = []
    for j in free:
        v = 1 << j
        for row, pivot in zip(a, pivots):
            v |= ((row >> j) & 1) << pivot
        basis.append(v)
    return basis, free



def restrict(v, keep):
    return sum(((v >> j) & 1) << i for i, j in enumerate(keep))



def retained_orders(d, multiplicity):
    keep = []
    for u in range(multiplicity - 1, -1, -1):
        if not any(math.comb(d - u, v - u) & 1 for v in keep):
            keep.append(u)
    return sorted(keep)



def kappa(k, d, p):
    coeffs = [
        sum((-1) ** j * math.comb(p, j) * math.comb(k, h - 2 * j) for j in range(min(p, h // 2) + 1))
        for h in range(d + 1)
    ]
    return (coeffs[-1] if min(coeffs) > 0 else 0), coeffs



