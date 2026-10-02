"""NumPy arithmetic in GF(2^m), independent of Sage and of the C++ backend.

Elements are integers whose bit i is the coefficient of z^i modulo the field
polynomial. Used to recheck binary Goppa results.
"""
from functools import lru_cache
import numpy as np


def _multiply(a, b, modulus, m):
    r = 0
    while b:
        if b & 1:
            r ^= a
        b >>= 1
        a <<= 1
        if a >> m:
            a ^= modulus
    return r


@lru_cache(maxsize=8)
def tables(modulus):
    m = modulus.bit_length() - 1
    q = 1 << m
    for g in range(2, q):
        powers = [1]
        x = g
        while x != 1 and len(powers) < q:
            powers.append(x)
            x = _multiply(x, g, modulus, m)
        if len(powers) == q - 1:
            break
    else:
        raise ValueError('no primitive element')
    exp = np.array(powers + powers, dtype=np.int64)
    log = np.zeros(q, dtype=np.int64)
    log[np.array(powers)] = np.arange(q - 1)
    return m, exp, log


def multiply(modulus, a, b):
    m, exp, log = tables(modulus)
    a = np.asarray(a, dtype=np.int64); b = np.asarray(b, dtype=np.int64)
    out = exp[log[a] + log[b]]
    return np.where((a == 0) | (b == 0), 0, out)


def inverse(modulus, a):
    m, exp, log = tables(modulus)
    a = np.asarray(a, dtype=np.int64)
    assert np.all(a != 0), 'division by zero'
    return exp[(len(log) - 1 - log[a]) % (len(log) - 1)]


def evaluate(modulus, poly, points):
    y = np.zeros(len(points), dtype=np.int64)
    for c in reversed(poly):
        y = multiply(modulus, y, points) ^ c
    return y


def support_derivatives(modulus, alpha):
    """Pi'(alpha_j) = prod_{l != j} (alpha_j - alpha_l)."""
    m, exp, log = tables(modulus)
    a = np.asarray(alpha, dtype=np.int64)
    assert len(np.unique(a)) == len(a), 'support is not distinct'
    # Diagonal entries contribute log[0] = 0, the exponent of 1.
    return exp[log[a[:, None] ^ a[None, :]].sum(axis=1) % (len(log) - 1)]


def binary_in_grs(modulus, Y, alpha, dimension, multiplier):
    """Check that the binary rows of Y lie in GRS_dimension(alpha, multiplier).

    The parity check has rows alpha_j^r / (multiplier_j Pi'(alpha_j)) for
    0 <= r < n - dimension. Each bit plane of Y H^T is an integer product mod 2.
    """
    Y = np.asarray(Y, dtype=np.float64)
    alpha = np.asarray(alpha, dtype=np.int64)
    n = len(alpha)
    m = tables(modulus)[0]
    h = inverse(modulus, multiply(modulus, multiplier, support_derivatives(modulus, alpha)))
    rows = []
    for _ in range(n - dimension):
        rows.append(h)
        h = multiply(modulus, h, alpha)
    H = np.array(rows, dtype=np.int64)
    for b in range(m):
        counts = Y @ ((H >> b) & 1).astype(np.float64).T
        if np.any(np.fmod(counts, 2)):
            return False
    return True
