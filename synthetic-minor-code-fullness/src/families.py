"""Parameter sets and random instances of the four code families."""
import json
from sage.all import *
from minor import dual_weights


GRS = [
    dict(id='A', q=32, n=20, k=5, D=9),
    dict(id='B', q=32, n=30, k=7, D=18),
    dict(id='C', q=64, n=40, k=8, D=24),
    dict(id='D', q=64, n=55, k=9, D=23),
    dict(id='E', q=128, n=100, k=11, D=45),
    dict(id='F', q=256, n=200, k=19, D=120),
    dict(id='G', q=512, n=300, k=23, D=200),
]
ALTERNANT = [
    dict(id='a1', q=2, m=5, n=28, r=3),
    dict(id='a2', q=2, m=6, n=50, r=3),
    dict(id='a3', q=2, m=6, n=60, r=8),
    dict(id='a4', q=2, m=7, n=100, r=5),
    dict(id='a5', q=2, m=8, n=200, r=8),
    dict(id='a6', q=2, m=8, n=220, r=12),
    dict(id='a7', q=4, m=2, n=14, r=3),
    dict(id='a8', q=4, m=3, n=50, r=5),
    dict(id='a9', q=4, m=4, n=180, r=6),
    dict(id='b1', q=8, m=2, n=55, r=5),
    dict(id='b2', q=16, m=2, n=200, r=6),
]
WILD = [
    dict(id='w1', q=4, m=3, t=1, g='irred', n=63),
    dict(id='w2', q=4, m=3, t=2, g='irred', n=64),
    dict(id='w3', q=4, m=3, t=2, g='split', n=62),
    dict(id='w4', q=4, m=3, t=3, g='irred', n=64),
    dict(id='w5', q=4, m=3, t=3, g='split', n=61),
    dict(id='w6', q=8, m=2, t=2, g='irred', n=64),
    dict(id='w7', q=8, m=2, t=2, g='irred', n=55),
    dict(id='w8', q=8, m=2, t=2, g='split', n=62),
    dict(id='w9', q=8, m=2, t=3, g='irred', n=64),
    dict(id='w10', q=8, m=3, t=1, g='irred', n=120),
    dict(id='w11', q=8, m=3, t=2, g='split', n=150),
    dict(id='w12', q=16, m=2, t=2, g='irred', n=256),
    dict(id='w13', q=16, m=2, t=2, g='irred', n=200),
    dict(id='w14', q=16, m=3, t=1, g='irred', n=120),
    dict(id='w15', q=16, m=3, t=2, g='irred', n=200),
]
BINARY = [
    dict(id='g1', m=4, n0=16, t=2, G='irred'),
    dict(id='g2', m=5, n0=28, t=2, G='irred'),
    dict(id='g3', m=6, n0=50, t=3, G='irred'),
    dict(id='g4', m=6, n0=55, t=3, G='split'),
    dict(id='g5', m=6, n0=60, t=5, G='separable'),
    dict(id='g6', m=7, n0=100, t=5, G='irred'),
    dict(id='g7', m=7, n0=110, t=8, G='split'),
    dict(id='g8', m=8, n0=200, t=8, G='irred'),
    dict(id='g9', m=8, n0=220, t=12, G='separable'),
    dict(id='s1', m=7, n0=128, t=5, G='irred', ell=28),
    dict(id='s2', m=8, n0=256, t=8, G='irred', ell=106),
]
# Field polynomials of GF(2^m) for the binary Goppa family, as integers whose
# bit i is the coefficient of z^i; m = 12 and 13 use the Classic McEliece fields.
MODULI = {4: 0x13, 5: 0x25, 6: 0x43, 7: 0x83, 8: 0x11d, 12: 0x1009, 13: 0x201b}
# Classic McEliece (m, n0, t), shortened to dimension k = n0 - mt - ell.
CLASSIC = {'mceliece348864': (12, 3488, 64), 'mceliece460896': (13, 4608, 96),
           'mceliece6688128': (13, 6688, 128), 'mceliece6960119': (13, 6960, 119),
           'mceliece8192128': (13, 8192, 128)}
SHORTENED_DIMENSIONS = {
    'mceliece348864': (154, 217), 'mceliece460896': (195, 274),
    'mceliece6688128': (224, 315), 'mceliece6960119': (216, 304),
    'mceliece8192128': (224, 315)}


def all_parameters():
    out = []
    for family, table in (('grs', GRS), ('alternant', ALTERNANT), ('wild', WILD), ('binary', BINARY)):
        out += [dict(P, family=family, label=family+'-'+P['id']) for P in table]
    for name, (m, n0, t) in CLASSIC.items():
        for k in SHORTENED_DIMENSIONS[name]:
            out.append(dict(id='%s-k%d' % (name, k), family='binary', m=m, n0=n0, t=t,
                            G='irred', ell=n0-m*t-k, label='%s-k%d' % (name, k)))
    return out


# GRS subcodes

def build_grs(P, seed):
    """Random k-dimensional subcode of GRS_{D+1}(alpha, lam) over F_q."""
    q, n, k, D = P["q"], P["n"], P["k"], P["D"]
    set_random_seed(seed)
    F = GF(q, 'a')
    points = list(F)
    shuffle(points)
    alpha = points[:n]
    lam = [F.random_element() for _ in range(n)]
    while 0 in lam:
        lam = [F.random_element() for _ in range(n)]
    M = random_matrix(F, k, D + 1)
    while M.rank() < k:
        M = random_matrix(F, k, D + 1)
    V = matrix(F, D + 1, n, lambda l, j: alpha[j] ** l)
    Yraw = M * V * diagonal_matrix(F, lam)
    return {"E": F, "Y": Yraw.echelon_form(), "alpha": alpha, "lam": lam,
            "nu": dual_weights(F, alpha, lam), "D": D}


# Alternant codes

def build_alternant(P, seed):
    """Kernel over F_q of the r x n matrix (beta_j alpha_j^i) over F_{q^m}."""
    q, m, n, r = P["q"], P["m"], P["n"], P["r"]
    set_random_seed(seed)
    E = GF(q ** m, 'b')
    S, emb = E.subfield(E.degree() // m, 'a', map=True)
    V, fromV, toV = E.vector_space(emb, map=True)
    points = list(E)
    shuffle(points)
    alpha = points[:n]
    beta = [E.random_element() for _ in range(n)]
    while 0 in beta:
        beta = [E.random_element() for _ in range(n)]
    H = matrix(E, r, n, lambda i, j: beta[j] * alpha[j] ** i)
    rows = []
    for i in range(r):
        cols = [toV(H[i][j]) for j in range(n)]
        for l in range(m):
            rows.append([cols[j][l] for j in range(n)])
    Ys = matrix(S, rows).right_kernel_matrix().echelon_form()
    Y = matrix(E, Ys.nrows(), n, lambda i, j: emb(Ys[i][j]))
    R = PolynomialRing(E, 'Z')
    Z = R.gen()
    dPi = prod(Z - a for a in alpha).derivative()
    lam = [1 / (beta[j] * dPi(alpha[j])) for j in range(n)]
    return {"E": E, "Y": Y, "alpha": alpha, "lam": lam, "nu": list(beta),
            "D": n - r - 1, "beta": beta}


# Wild Goppa codes

def goppa_subfield_subcode(E, S, toV, m, alpha, poly):
    n = len(alpha)
    r = poly.degree()
    H = matrix(E, r, n, lambda i, j: alpha[j] ** i / poly(alpha[j]))
    rows = []
    for i in range(r):
        cols = [toV(H[i][j]) for j in range(n)]
        for l in range(m):
            rows.append([cols[j][l] for j in range(n)])
    return matrix(S, rows).right_kernel_matrix().echelon_form()


def build_wild(P, seed):
    """Gamma(alpha, g^(q-1)) = Gamma(alpha, g^q) over F_q, with g square-free."""
    q, m, t, n = P["q"], P["m"], P["t"], P["n"]
    set_random_seed(seed)
    E = GF(q ** m, 'b')
    S, emb = E.subfield(E.degree() // m, 'a', map=True)
    V, fromV, toV = E.vector_space(emb, map=True)
    R = PolynomialRing(E, 'Z')
    Z = R.gen()
    if P["g"] == "split":
        pts = list(E)
        shuffle(pts)
        g = prod(Z - x for x in pts[:t])
    else:
        g = R.irreducible_element(t, algorithm='random')
    assert g.is_squarefree()
    avail = [a for a in E if g(a) != 0]
    shuffle(avail)
    alpha = avail[:n]
    Ys = goppa_subfield_subcode(E, S, toV, m, alpha, g ** q)
    assert Ys.row_space() == goppa_subfield_subcode(
        E, S, toV, m, alpha, g ** (q - 1)).row_space(), "the wild identity failed"
    Y = matrix(E, Ys.nrows(), len(alpha), lambda i, j: emb(Ys[i][j]))
    dPi = prod(Z - a for a in alpha).derivative()
    lam = [g(alpha[j]) ** q / dPi(alpha[j]) for j in range(len(alpha))]
    nu = [1 / g(alpha[j]) ** q for j in range(len(alpha))]
    return {"E": E, "Y": Y, "alpha": alpha, "lam": lam, "nu": nu,
            "D": len(alpha) - q * t - 1, "g": g}


# Binary Goppa codes (generated in C++; this construction is the reference)

def goppa_matrix(E, alpha, G):
    """Reference construction: kernel of the expanded check alpha_j^i / G(alpha_j)^2."""
    n = len(alpha)
    m = E.degree()
    V, fromV, toV = E.vector_space(map=True)
    ginv = [1 / G(a) ** 2 for a in alpha]
    rowsE = [vector(E, ginv)]
    av = vector(E, alpha)
    for _ in range(2 * G.degree() - 1):
        rowsE.append(rowsE[-1].pairwise_product(av))
    rows = []
    for rv in rowsE:
        cols = [toV(x) for x in rv]
        for l in range(m):
            rows.append([cols[j][l] for j in range(n)])
    return matrix(GF(2), rows).right_kernel_matrix().echelon_form()


BUILDERS = dict(grs=build_grs, alternant=build_alternant, wild=build_wild)


def build_instance(P, seed):
    return BUILDERS[P['family']](P, seed)


if __name__ == '__main__':
    print(json.dumps(all_parameters(), indent=2))
