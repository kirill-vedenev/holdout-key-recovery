"""Public GRS recovery helpers extracted unchanged from the existing implementation."""
from sage.all import *
from functools import lru_cache
import numpy as np

@lru_cache(maxsize=8)
def binary_tables(E):
    """Polynomial-basis encodings, with deterministic primitive-element search."""
    q = int(E.order())
    assert E.characteristic() == 2 and q <= 65536
    factors = prime_divisors(q-1)
    primitive = next(E.from_integer(i) for i in range(1,q)
                     if all(E.from_integer(i)**((q-1)//p) != 1 for p in factors))
    powers = np.empty(q-1,dtype=np.uint16)
    x = E.one()
    for i in range(q-1):
        powers[i] = int(x.to_integer())
        x *= primitive
    logs = np.zeros(q,dtype=np.int64)
    logs[powers] = np.arange(q-1)
    return powers,logs


def evaluate_many(poly, points):
    """Horner evaluation with NumPy field tables, avoiding scalar Sage calls."""
    E = poly.base_ring()
    powers,logs = binary_tables(E)
    x = np.array([int(a.to_integer()) for a in points],dtype=np.uint16)
    y = np.zeros(len(points),dtype=np.uint16)
    for c in reversed(poly.list()):
        zero = (y==0)|(x==0)
        y = powers[(logs[y]+logs[x]) % len(powers)]
        y[zero] = 0
        y ^= int(c.to_integer())
    return [E.from_integer(int(v)) for v in y]


@lru_cache(maxsize=8)
def _support_derivatives(E, points):
    powers,logs = binary_tables(E)
    x = np.array([int(a.to_integer()) for a in points],dtype=np.uint16)
    assert len(np.unique(x)) == len(x), 'support is not distinct'
    # Diagonal entries are omitted: log[0]=0 contributes the exponent of 1.
    exponents = logs[x[:,None]^x[None,:]].sum(axis=1) % len(powers)
    return tuple(E.from_integer(int(v)) for v in powers[exponents])


def support_derivatives(E, points):
    """Values of the support polynomial's derivative at its distinct roots."""
    return _support_derivatives(E,tuple(points))


def grs_parity(E, alpha, dimension, multiplier=None):
    """Direct dual-GRS formula, avoiding an expensive generic right kernel."""
    n = len(alpha)
    if not 0 < dimension < n:
        raise ValueError('GRS dimension outside 1..n-1')
    dpi = support_derivatives(E,alpha)
    if multiplier is None:
        multiplier = [E.one()]*n
    row = vector(E, [1/(multiplier[j]*dpi[j]) for j in range(n)])
    av = vector(E, alpha)
    rows = [row]
    for _ in range(n-dimension-1):
        rows.append(rows[-1].pairwise_product(av))
    return matrix(E, rows)


def grs_multiplier_fast(E, G, alpha, dimension, probe_rows=8):
    n = len(alpha)
    if dimension < 1 or dimension >= n:
        return None
    H = grs_parity(E, alpha, dimension)
    Gr = G.echelon_form()[:G.rank()]
    used = min(Gr.nrows(), max(probe_rows, int(ceil((n+1)/H.nrows()))))
    while True:
        Gp = Gr if used == Gr.nrows() else random_matrix(E, used, Gr.nrows())*Gr
        blocks = [Gp.elementwise_product(matrix(E, [h]*used)) for h in H.rows()]
        sol = block_matrix(len(blocks), 1, blocks).right_kernel()
        if sol.dimension() <= 1 or used == Gr.nrows():
            break
        used = min(Gr.nrows(), 2*used)
    if not sol.dimension():
        return None
    candidates = [v for v in sol.basis() if 0 not in v]
    for _ in range(200):
        if len(candidates) >= 4:
            break
        v = sol.random_element()
        if 0 not in v:
            candidates.append(v)
    for b in candidates:
        if (Gr.elementwise_product(matrix(E, [b]*Gr.nrows()))*H.transpose()).is_zero():
            return vector(E, [1/x for x in b])
    return None


def grs_multiplier_systematic(E, G, alpha):
    """Recover a full GRS code's multiplier from its systematic Cauchy entries.

    Interpolation gives R[i,j] = v[j] L_i(alpha[j]) / v[pivot[i]].
    One non-pivot column fixes the pivot multipliers; one row fixes the rest.
    The final parity check certifies every row, including non-GRS inputs.
    """
    R = G.echelon_form()
    pivots = list(R.pivots())
    k, n = len(pivots), G.ncols()
    if not 1 < k < n:
        return None
    R = R[:k]
    pivot_set = set(pivots)
    rest = [j for j in range(n) if j not in pivot_set]
    ring = PolynomialRing(E, 'Z'); Z = ring.gen()
    poly = prod(Z-alpha[j] for j in pivots)
    dp = support_derivatives(E,[alpha[j] for j in pivots])
    values = evaluate_many(poly,[alpha[j] for j in rest])
    v = [E.zero()]*n
    for j,value in zip(rest,values):
        v[j] = R[0,j]*(alpha[j]-alpha[pivots[0]])*dp[0]/value
    j = rest[0]
    if not v[j] or any(not R[i,j] for i in range(k)):
        return None
    value = values[0]
    for i, p in enumerate(pivots):
        v[p] = v[j]*value/((alpha[j]-alpha[p])*dp[i]*R[i,j])
    if any(not x for x in v):
        return None
    if not (R*grs_parity(E,alpha,k,v).transpose()).is_zero():
        return None
    return vector(E,v)


def recover_GRS_support(G, x0, x1, xk, G1=None):
    n = G.ncols()
    k = G.nrows()
    if G1 is None:
        G1 = G.rref()
    x = [None] * n
    x[0], x[1], x[k] = x0, x1, xk
    seen = set([x0, x1, xk])
    i = 0; i_ = 1; j_ = k
    for j in range(k + 1, n):
        ratio = G1[i, j_] * G1[i_, j]
        if ratio == 0:          # a zero here means the code is not MDS, so not GRS
            return x, False
        gamma = (G1[i, j] * G1[i_, j_]) / ratio
        denom = (x[j_] - x[i]) - gamma * (x[j_] - x[i_])
        if denom == 0:
            return x, False
        xj = (x[i_] * (x[j_] - x[i]) - gamma * x[i] * (x[j_] - x[i_])) / denom
        if xj in seen:
            return x, False
        x[j] = xj
        seen.add(xj)
    i_ = 0; j = k; j_ = k + 1
    for i in range(2, k):
        ratio = G1[i, j_] * G1[i_, j]
        if ratio == 0:
            return x, False
        gamma = (G1[i, j] * G1[i_, j_]) / ratio
        denom = gamma * (x[j_] - x[i_]) - (x[j] - x[i_])
        if denom == 0:
            return x, False
        xi = (gamma * x[j] * (x[j_] - x[i_]) - (x[j] - x[i_]) * x[j_]) / denom
        if xi in seen:
            return x, False
        x[i] = xi
        seen.add(xi)
    return x, True

def sidelnikov_shestakov(F, G, tries=200):
    kappa = G.rank()
    if kappa < 3 or G.ncols() - kappa < 2:
        return None
    R = G.echelon_form()[:kappa]
    pivots = list(R.pivots())
    order = pivots + [j for j in range(G.ncols()) if j not in pivots]
    Gp = R.matrix_from_columns(order)
    G1 = Gp.rref()
    attempts = 0
    for xk in F:
        if xk == 0 or xk == 1:
            continue
        attempts += 1
        if attempts > tries:
            return None
        x, ok = recover_GRS_support(Gp, F(0), F(1), xk, G1)
        if ok:
            out = [None] * G.ncols()
            for idx, col in enumerate(order):
                out[col] = x[idx]
            return out
    return None
