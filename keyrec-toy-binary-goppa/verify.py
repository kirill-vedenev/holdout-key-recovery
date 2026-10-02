#!/usr/bin/env python3
"""Independent standard-library verification of a recovered toy Goppa key.

No Sage, native recovery library, original secret, jets, or kernel are read.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import time


def binary_rref(rows, width):
    rows = list(rows)
    rank = 0
    for c in range(width):
        pivot = next((r for r in range(rank, len(rows)) if (rows[r] >> c) & 1), None)
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        for r in range(len(rows)):
            if r != rank and (rows[r] >> c) & 1:
                rows[r] ^= rows[rank]
        rank += 1
    return rows[:rank]


def verify(public_path, key_path):
    started = time.monotonic()
    words = public_path.read_text().split()
    if words[0] != 'BGPK1' or len(words) != 73:
        raise ValueError('invalid public input format')
    m, t, n, k, modulus, heldout, degree, multiplicity = map(int, words[1:9])
    if (m, t, n, k, degree, multiplicity) != (6, 6, 64, 28, 5, 4):
        raise ValueError('unexpected toy parameters')
    columns = list(map(int, words[9:]))
    if any(c < 0 or c >= 1 << k for c in columns):
        raise ValueError('invalid binary column')
    public = [sum(((columns[j] >> a) & 1) << j for j in range(n)) for a in range(k)]
    if len(binary_rref(public, n)) != k:
        raise ValueError('public generator rank')

    def mul(a, b):
        value = 0
        while b:
            if b & 1:
                value ^= a
            b >>= 1
            a <<= 1
            if a & (1 << m):
                a ^= modulus
        return value

    def power(a, exponent):
        value = 1
        while exponent:
            if exponent & 1:
                value = mul(value, a)
            a = mul(a, a)
            exponent >>= 1
        return value

    def inverse(a):
        if not a:
            raise ValueError('zero denominator')
        value = power(a, (1 << m)-2)
        if mul(value, a) != 1:
            raise ValueError('invalid field modulus')
        return value

    def evaluate(poly, x):
        value = 0
        for coefficient in reversed(poly):
            value = mul(value, x) ^ coefficient
        return value

    def trim(poly):
        poly = list(poly)
        while poly and not poly[-1]:
            poly.pop()
        return poly

    def add(a, b):
        result = list(a) + [0]*max(0, len(b)-len(a))
        for i, c in enumerate(b):
            result[i] ^= c
        return trim(result)

    def remainder(a, b):
        a, b = trim(a), trim(b)
        if not b:
            raise ValueError('zero polynomial divisor')
        while len(a) >= len(b):
            offset = len(a)-len(b)
            factor = mul(a[-1], inverse(b[-1]))
            for i, c in enumerate(b):
                a[offset+i] ^= mul(factor, c)
            a = trim(a)
        return a

    def gcd(a, b):
        while b:
            a, b = b, remainder(a, b)
        return [mul(c, inverse(a[-1])) for c in a] if a else []

    def square_mod(a, g):
        result = [0]*(2*len(a))
        for i, c in enumerate(a):
            result[2*i] = mul(c, c)
        return remainder(result, g)

    key = json.loads(key_path.read_text())
    if any(key[name] != value for name, value in [('m', m), ('t', t), ('n', n), ('modulus', modulus)]):
        raise ValueError('key parameter mismatch')
    support, weights, g = key['support'], key['grs_multiplier'], key['goppa_polynomial']
    for values, length in [(support, n), (weights, n), (g, t+1)]:
        if len(values) != length or any(type(x) is not int or not 0 <= x < 1 << m for x in values):
            raise ValueError('key field encoding or length')
    if len(set(support)) != n or not all(weights) or g[-1] != 1:
        raise ValueError('repeated support, zero multiplier, or nonmonic polynomial')
    h = [0, 1]
    for r in range(1, t+1):
        for _ in range(m):
            h = square_mod(h, g)
        if r in (2, 3) and gcd(add(h, [0, 1]), g) != [1]:
            raise ValueError('Goppa polynomial is reducible')
    if add(h, [0, 1]):
        raise ValueError('Goppa polynomial fails Frobenius irreducibility test')

    derivative_values = []
    for i, a in enumerate(support):
        value = 1
        for j, b in enumerate(support):
            if i != j:
                value = mul(value, a ^ b)
        derivative_values.append(value)
    g_values = [evaluate(g, a) for a in support]
    if not all(g_values):
        raise ValueError('Goppa polynomial has a support root')

    def check_code(multipliers, dimension):
        rows = [0]*(m*dimension)
        for j, (a, value) in enumerate(zip(support, multipliers)):
            for r in range(dimension):
                for bit in range(m):
                    rows[r*m+bit] |= ((value >> bit) & 1) << j
                value = mul(value, a)
        rank = len(binary_rref(rows, n))
        if rank != n-k or any((a & b).bit_count() & 1 for a in public for b in rows):
            raise ValueError('reconstructed code differs from the original public code')
        return rank

    dual = [inverse(mul(a, b)) for a, b in zip(weights, derivative_values)]
    alternant_rank = check_code(dual, 2*t)
    goppa_rank = check_code([inverse(a) for a in g_values], t)
    scales = [mul(mul(a, b), inverse(mul(c, c))) for a, b, c in zip(weights, derivative_values, g_values)]
    if len(set(scales)) != 1:
        raise ValueError('multiplier square identity failed')
    return dict(status='passed', public_only=True, independent_of_recovery=True,
                support_count=n, distinct_finite_support=True, goppa_degree=t,
                goppa_irreducible=True, alternant_binary_rank=alternant_rank,
                goppa_binary_rank=goppa_rank, exact_public_code_equality=True,
                multiplier_square_identity=True,
                public_sha256=hashlib.sha256(public_path.read_bytes()).hexdigest(),
                recovered_key_sha256=hashlib.sha256(key_path.read_bytes()).hexdigest(),
                seconds=time.monotonic()-started)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('public', type=Path)
    parser.add_argument('key', type=Path)
    parser.add_argument('--report', type=Path)
    args = parser.parse_args()
    result = verify(args.public, args.key)
    text = json.dumps(result, indent=2)+'\n'
    if args.report:
        args.report.write_text(text)
    print(text, end='')


if __name__ == '__main__':
    try:
        main()
    except (ValueError, KeyError, OSError) as error:
        print(f'ERROR: {error}', file=sys.stderr)
        sys.exit(1)
