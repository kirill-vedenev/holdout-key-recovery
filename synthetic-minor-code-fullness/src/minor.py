"""Sage reference implementation: tangent planes, minor codes, GRS checks and
support recovery over fields of characteristic two."""
from functools import lru_cache
import numpy as np
from sage.all import *


# Tangent planes from the secret key

def dual_weights(F, alpha, lam):
    """nu_j = 1 / (lam_j Pi'(alpha_j))."""
    R = PolynomialRing(F, 'Z')
    Z = R.gen()
    dPi = prod(Z - a for a in alpha).derivative()
    return [1 / (lam[j] * dPi(alpha[j])) for j in range(len(alpha))]


def tangent_matrix(F, Y, alpha, nu):
    """Column j is sum_{l != j} nu_l y_l / (alpha_j - alpha_l)."""
    n = len(alpha)
    Delta = matrix(F, n, n, lambda i, j: 0 if i == j else 1 / (alpha[i] - alpha[j]))
    return Y * diagonal_matrix(F, nu) * Delta.transpose()


def randomise_representatives(F, Y, Yhat, seed):
    """yhat_j <- a_j yhat_j + b_j y_j with a_j nonzero. Returns the new Yhat and a."""
    set_random_seed(seed)
    n = Y.ncols()
    a = []
    for _ in range(n):
        c = F.random_element()
        while c == 0:
            c = F.random_element()
        a.append(c)
    b = [F.random_element() for _ in range(n)]
    return Yhat * diagonal_matrix(F, a) + Y * diagonal_matrix(F, b), a


def tangent_defects(Y, Yhat):
    """Positions j where y_j and yhat_j do not span a plane."""
    return [i for i in range(Y.ncols())
            if matrix([Y.column(i), Yhat.column(i)]).rank() != 2]


def hidden_curve(F, Y, alpha, lam, D):
    """Interpolate the polynomials f_a with y_{a,j} = lam_j f_a(alpha_j)."""
    R = PolynomialRing(F, 'Z')
    f = [R.lagrange_polynomial([(alpha[j], Y[a][j] / lam[j]) for j in range(len(alpha))])
         for a in range(Y.nrows())]
    assert max(p.degree() for p in f) <= D, 'ambient degree exceeded'
    return f


def curve_matches_closed_form(F, Y, Yhat, curve, alpha):
    """Check span{y_j, yhat_j} = span{F(alpha_j), F'(alpha_j)} at every position."""
    for i in range(len(alpha)):
        honest = [vector(F, [p(alpha[i]) for p in curve]),
                  vector(F, [p.derivative()(alpha[i]) for p in curve])]
        got = [Y.column(i), Yhat.column(i)]
        if not (matrix(got).rank() == matrix(honest).rank()
                == matrix(got + honest).rank()):
            return False
    return True


# Minor matrix and minor code

def minor_matrix(Y, Yhat, pairs=None):
    """Rows (a,b), a < b, with entries y_{a,j} yhat_{b,j} - y_{b,j} yhat_{a,j}."""
    F = Y.base_ring()
    k = Y.nrows()
    if pairs is None:
        partners = {a: list(range(a + 1, k)) for a in range(k - 1)}
    else:
        partners = {}
        for a, b in pairs:
            partners.setdefault(a, []).append(b)
    blocks = []
    for a in sorted(partners):
        bs = sorted(partners[a])
        ya = matrix(F, [Y.row(a)] * len(bs))
        ha = matrix(F, [Yhat.row(a)] * len(bs))
        blocks.append(ya.elementwise_product(Yhat.matrix_from_rows(bs))
                      - ha.elementwise_product(Y.matrix_from_rows(bs)))
    return block_matrix(len(blocks), 1, blocks)


def square_root(M):
    F = M.base_ring()
    return M.apply_map(lambda x: x ** (F.order() // 2))


def minor_code(Y, Yhat, pairs=None):
    """The matrix T: entrywise square roots of the minor matrix."""
    return square_root(minor_matrix(Y, Yhat, pairs))


def sample_pairs(k, count, seed):
    """Random distinct pairs a < b; None means all binomial(k,2) pairs."""
    if binomial(k, 2) <= count:
        return None
    set_random_seed(seed)
    chosen = set()
    while len(chosen) < count:
        a = randint(0, k - 1)
        b = randint(0, k - 2)
        if b >= a:
            b += 1
        a, b = min(a, b), max(a, b)
        chosen.add((a, b))
    return sorted(chosen)


def row_budget(k, bound):
    return min(binomial(k, 2), max(600, 2 * bound + 50))


def all_minors_inside(Y, Yhat, H):
    """Check every minor row against every row h of H, without forming them.

    The (a,b) entry of A_h - A_h^T, with A_h = Y diag(h_j^2) Yhat^T, is the
    square of <T_(a,b), h>.
    """
    E = Y.base_ring()
    for h in H.rows():
        weights = vector(E, [x ** 2 for x in h])
        A = Y.elementwise_product(matrix(E, [weights] * Y.nrows())) * Yhat.transpose()
        if A != A.transpose():
            return False
    return True


# GRS codes

@lru_cache(maxsize=8)
def binary_tables(E):
    """Antilog and log tables of the field, from the first primitive element."""
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
    """Horner evaluation with the field tables."""
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
    """Values Pi'(alpha_j) of the support polynomial's derivative."""
    return _support_derivatives(E,tuple(points))


def grs_generator(F, support, dimension, mult=None):
    if mult is None:
        mult = [F(1)] * len(support)
    return matrix(F, dimension, len(support), lambda i, j: mult[j] * support[j] ** i)


def grs_parity(E, alpha, dimension, multiplier=None):
    """Parity-check matrix V_{n-dimension}(alpha, multiplier^perp) of GRS_dimension."""
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
    """A multiplier v with rowspace(G) inside GRS_dimension(alpha, v), or None.

    Solves G diag(z) H^T = 0 on dense row combinations, then checks every row.
    """
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
    """Multiplier of a full GRS code from its systematic Cauchy entries, or None.

    Interpolation gives R[i,j] = v[j] L_i(alpha[j]) / v[pivot[i]].
    One non-pivot column fixes the pivot multipliers; one row fixes the rest.
    The final parity check covers every row and rejects non-GRS inputs.
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


# Sidelnikov-Shestakov support recovery

def recover_GRS_support(G, x0, x1, xk, G1=None):
    """Support from the cross-ratios of a systematic generator, with the first
    two pivot coordinates and the first non-pivot coordinate fixed."""
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
    """A candidate support of the code spanned by G, or None."""
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


def recover_public_embedding(E, Y, T, D):
    """Recover a support from T and embed the public code Y in GRS_{D+1}.

    Inputs are public: the field, Y, the minor matrix T and the degree D.
    The smaller of the minor code and its dual is passed to Sidelnikov-Shestakov.
    """
    rank = T.rank()
    basis = T.echelon_form()[:rank]
    Gs = basis.right_kernel_matrix() if rank > T.ncols()-rank else basis
    recovered_support = sidelnikov_shestakov(E,Gs,tries=int(E.order()))
    support_ok = False
    embedding = False
    if recovered_support is not None:
        em = grs_multiplier_systematic(E,Gs,recovered_support)
        support_ok = em is not None
        if support_ok:
            lm = grs_multiplier_fast(E,Y,recovered_support,D+1)
            embedding = lm is not None and (Y*grs_parity(E,recovered_support,D+1,lm).transpose()).is_zero()
    return dict(support_recovered=bool(support_ok),embedding_verified=bool(embedding))
