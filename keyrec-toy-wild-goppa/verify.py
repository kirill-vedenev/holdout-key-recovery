#!/usr/bin/env python3
"""Independent public verification of wild Goppa recovery over GF(4).

Uses only Python's standard library. No kernel, jet, original secret, or native
recovery component is read or imported.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import time


def multiply(a, b, modulus, size):
    value = 0
    while b:
        if b & 1:
            value ^= a
        a <<= 1
        b >>= 1
        if a & size:
            a ^= modulus
    return value


MUL4 = [[multiply(a, b, 7, 4) for b in range(4)] for a in range(4)]
MUL64 = [[multiply(a, b, 67, 64) for b in range(64)] for a in range(64)]
EMBED = [0, 1, 58, 59]
COMPONENTS = {}
for a in range(4):
    for b in range(4):
        for c in range(4):
            x = EMBED[a] ^ MUL64[EMBED[b]][2] ^ MUL64[EMBED[c]][4]
            if x in COMPONENTS:
                raise ValueError('subfield basis is not injective')
            COMPONENTS[x] = (a, b, c)


def power(a, exponent):
    result = 1
    while exponent:
        if exponent & 1:
            result = MUL64[result][a]
        a = MUL64[a][a]
        exponent >>= 1
    return result


def inverse(a):
    if not a:
        raise ValueError('zero denominator')
    return power(a, 62)


def rref4(rows, width):
    rows = [list(row) for row in rows]
    rank = 0
    for c in range(width):
        pivot = next((r for r in range(rank, len(rows)) if rows[r][c]), None)
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        scale = MUL4[rows[rank][c]][rows[rank][c]]
        rows[rank] = [MUL4[scale][x] for x in rows[rank]]
        for r in range(len(rows)):
            if r != rank and rows[r][c]:
                scale = rows[r][c]
                rows[r] = [a ^ MUL4[scale][b] for a, b in zip(rows[r], rows[rank])]
        rank += 1
    return rows[:rank]


def evaluate(poly, a):
    result = 0
    for c in reversed(poly):
        result = MUL64[result][a] ^ c
    return result


def trim(poly):
    poly = list(poly)
    while poly and not poly[-1]:
        poly.pop()
    return poly


def add(a, b):
    result = list(a)+[0]*max(0, len(b)-len(a))
    for i, c in enumerate(b):
        result[i] ^= c
    return trim(result)


def remainder(a, b):
    a, b = trim(a), trim(b)
    if not b:
        raise ValueError('zero polynomial divisor')
    while len(a) >= len(b):
        offset = len(a)-len(b)
        scale = MUL64[a[-1]][inverse(b[-1])]
        for i, c in enumerate(b):
            a[offset+i] ^= MUL64[scale][c]
        a = trim(a)
    return a


def gcd(a, b):
    while b:
        a, b = b, remainder(a, b)
    return [MUL64[c][inverse(a[-1])] for c in a] if a else []


def irreducible(poly):
    degree = len(poly)-1
    factors = []
    left = degree
    for p in range(2, degree+1):
        if left % p == 0:
            factors.append(p)
            while left % p == 0:
                left //= p
    h = [0, 1]
    for r in range(1, degree+1):
        for _ in range(6):
            squared = [0]*(2*len(h))
            for i, c in enumerate(h):
                squared[2*i] = MUL64[c][c]
            h = remainder(squared, poly)
        if any(r == degree//p for p in factors) and gcd(add(h, [0, 1]), poly) != [1]:
            return False
    return not add(h, [0, 1])


def verify(public_path, key_path):
    started = time.monotonic()
    data = public_path.read_text().split()
    if not data or data[0] != 'WGPK1':
        raise ValueError('invalid public format')
    m, t, n, k, modulus, heldout, degree, multiplicity = map(int, data[1:9])
    if m != 3 or modulus != 67 or len(data) != n+9 or not 0 <= heldout < n:
        raise ValueError('invalid public parameters')
    columns = list(map(int, data[9:]))
    if any(a < 0 or a >> (2*k) for a in columns):
        raise ValueError('invalid public column')
    public = [[(columns[j] >> (2*a)) & 3 for j in range(n)] for a in range(k)]
    if len(rref4(public, n)) != k:
        raise ValueError('public generator rank')
    key = json.loads(key_path.read_text())
    if any(key[name] != value for name, value in [('q', 4), ('m', m), ('t', t), ('n', n), ('k', k), ('modulus', modulus), ('goppa_exponent', 3)]):
        raise ValueError('key parameter mismatch')
    support, weights, gamma = key['support'], key['grs_multiplier'], key['goppa_radical']
    for values, length in [(support, n), (weights, n), (gamma, t+1)]:
        if len(values) != length or any(type(a) is not int or not 0 <= a < 64 for a in values):
            raise ValueError('invalid key field encoding')
    if len(set(support)) != n or not all(weights) or gamma[-1] != 1 or not irreducible(gamma):
        raise ValueError('invalid support, multiplier, or irreducible radical')
    gamma_values = [evaluate(gamma, a) for a in support]
    if not all(gamma_values):
        raise ValueError('radical has a support root')
    derivatives = []
    for i, a in enumerate(support):
        value = 1
        for j, b in enumerate(support):
            if i != j:
                value = MUL64[value][a ^ b]
        derivatives.append(value)

    def check_code(multipliers, dimension):
        rows = [[0]*n for _ in range(3*dimension)]
        for j, (a, value) in enumerate(zip(support, multipliers)):
            for r in range(dimension):
                for b, component in enumerate(COMPONENTS[value]):
                    rows[3*r+b][j] = component
                value = MUL64[value][a]
        rank = len(rref4(rows, n))
        if rank != n-k:
            raise ValueError('reconstructed parity rank differs from the public codimension')
        for g in public:
            for h in rows:
                residual = 0
                for a, b in zip(g, h):
                    residual ^= MUL4[a][b]
                if residual:
                    raise ValueError('reconstructed code differs from the public code')
        return rank

    dual = [inverse(MUL64[a][b]) for a, b in zip(weights, derivatives)]
    alternant_rank = check_code(dual, 4*t)
    cube_rank = check_code([inverse(power(a, 3)) for a in gamma_values], 3*t)
    fourth_rank = check_code([inverse(power(a, 4)) for a in gamma_values], 4*t)
    scales = [MUL64[MUL64[a][b]][inverse(power(c, 4))] for a, b, c in zip(weights, derivatives, gamma_values)]
    if len(set(scales)) != 1:
        raise ValueError('wild Goppa multiplier identity failed')
    return dict(status='passed', public_only=True, independent_of_recovery=True,
                support_count=n, distinct_finite_support=True, radical_degree=t,
                radical_irreducible=True, alternant_rank=alternant_rank,
                goppa_cube_rank=cube_rank, goppa_fourth_power_rank=fourth_rank,
                exact_public_code_equality=True, multiplier_fourth_power_identity=True,
                public_sha256=hashlib.sha256(public_path.read_bytes()).hexdigest(),
                recovered_key_sha256=hashlib.sha256(key_path.read_bytes()).hexdigest(),
                seconds=time.monotonic()-started)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('public', type=Path)
    parser.add_argument('key', type=Path)
    parser.add_argument('--report', type=Path)
    args = parser.parse_args()
    output = json.dumps(verify(args.public, args.key), indent=2)+'\n'
    if args.report:
        args.report.write_text(output)
    print(output, end='')


if __name__ == '__main__':
    try:
        main()
    except (ValueError, KeyError, OSError) as error:
        print(f'ERROR: {error}', file=sys.stderr)
        sys.exit(1)
