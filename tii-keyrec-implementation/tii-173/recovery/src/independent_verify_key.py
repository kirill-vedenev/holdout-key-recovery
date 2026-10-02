#!/usr/bin/env python3
"""Independent, standard-library-only verification of a full TII-173 Goppa key.

No recovery modules or secret reference key are imported. Field elements use
the public polynomial-basis integer encoding; g is low coefficient first.
"""
import argparse
import hashlib
import json
from pathlib import Path
import time

PUBLIC_SHA256 = '8f40c2e315e7f7dbd41cfdb3ee3e5aafe9abf0b06693d5225d0953916cbe79bb'
PUBLIC = Path(__file__).resolve().parents[4] / 'tii-results' / 'tii-173' / 'pk_McEliece_173.txt'
MODULUS = 0x83


def mul(a, b):
    out = 0
    while b:
        if b & 1:
            out ^= a
        b >>= 1
        a <<= 1
        if a & 128:
            a ^= MODULUS
    return out


def power(a, n):
    out = 1
    while n:
        if n & 1:
            out = mul(out, a)
        a = mul(a, a)
        n >>= 1
    return out


def inverse(a):
    if not a:
        raise ValueError('division by zero')
    return power(a, 126)


def trim(p):
    p = list(p)
    while p and not p[-1]:
        p.pop()
    return p


def poly_add(a, b):
    c = list(a) + [0] * max(0, len(b)-len(a))
    for i, x in enumerate(b):
        c[i] ^= x
    return trim(c)


def remainder(a, b):
    a, b = trim(a), trim(b)
    if not b:
        raise ValueError('zero polynomial divisor')
    lead_inverse = inverse(b[-1])
    while a and len(a) >= len(b):
        shift = len(a)-len(b)
        coefficient = mul(a[-1], lead_inverse)
        for i, x in enumerate(b):
            a[shift+i] ^= mul(coefficient, x)
        a = trim(a)
    return a


def gcd(a, b):
    a, b = trim(a), trim(b)
    while b:
        a, b = b, remainder(a, b)
    return [mul(x, inverse(a[-1])) for x in a] if a else []


def square_mod(a, modulus):
    result = [0] * (2*len(a))
    for i, x in enumerate(a):
        result[2*i] = mul(x, x)
    return remainder(result, modulus)


def irreducible_degree_eight(g):
    # Rabin's exact criterion: 8 has the sole prime divisor 2.
    x = [0, 1]
    h = x
    for degree in range(1, 9):
        for _ in range(7):
            h = square_mod(h, g)
        if degree == 4 and gcd(poly_add(h, x), g) != [1]:
            return False
    return not poly_add(h, x)


def evaluate(poly, x):
    value = 0
    for a in reversed(poly):
        value = mul(value, x) ^ a
    return value


def rref(rows, n):
    rows = list(rows)
    rank = 0
    for col in range(n):
        pivot = next((j for j in range(rank, len(rows)) if rows[j] >> col & 1), None)
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        for j in range(len(rows)):
            if j != rank and rows[j] >> col & 1:
                rows[j] ^= rows[rank]
        rank += 1
    return tuple(rows[:rank])


def verify(candidate_path, public_path=PUBLIC):
    start = time.monotonic()
    candidate_path, public_path = Path(candidate_path), Path(public_path)
    candidate_raw, public_raw = candidate_path.read_bytes(), public_path.read_bytes()
    if hashlib.sha256(public_raw).hexdigest() != PUBLIC_SHA256:
        raise ValueError('public input differs from the authenticated original TII-173 file')
    lines = public_raw.decode().splitlines()
    if len(lines) != 57:
        raise ValueError('unexpected public matrix format')
    original = []
    for line in lines[:-1]:
        bits = line.strip()[1:-1].split()
        if len(bits) != 96 or any(x not in ('0', '1') for x in bits):
            raise ValueError('invalid original binary parity-check row')
        original.append(sum(int(x) << j for j, x in enumerate(bits)))
    field_modulus = json.loads(lines[-1])
    key = json.loads(candidate_raw)
    if key['field_modulus'] != field_modulus:
        raise ValueError('candidate uses a different field encoding')
    support, g = key['support'], key['g']
    if len(support) != 96 or len(g) != 9:
        raise ValueError('expected all 96 original support positions and a degree-8 polynomial')
    if any(type(x) is not int or not 0 <= x < 128 for x in support+g):
        raise ValueError('invalid canonical field element')
    if len(set(support)) != 96:
        raise ValueError('support elements are not distinct')
    if g[-1] != 1 or not irreducible_degree_eight(g):
        raise ValueError('the polynomial must be monic and irreducible of degree 8')
    recovered = [0] * 56
    for j, x in enumerate(support):
        gx = evaluate(g, x)
        if not gx:
            raise ValueError('Goppa polynomial vanishes on the support')
        v = inverse(gx)
        for degree in range(8):
            for bit in range(7):
                recovered[degree*7+bit] |= ((v >> bit) & 1) << j
            v = mul(v, x)
    actual_space, expected_space = rref(recovered, 96), rref(original, 96)
    if len(expected_space) != 56 or actual_space != expected_space:
        raise ValueError('reconstructed Goppa parity-check row space differs from the original public code')
    return dict(schema='tii173-independent-full-key-verification-v1',
                verified=True, original_public_sha256=PUBLIC_SHA256,
                candidate_sha256=hashlib.sha256(candidate_raw).hexdigest(),
                support_count=96, support_distinct=True, goppa_degree=8,
                polynomial_monic=True, polynomial_irreducible=True,
                polynomial_nonzero_on_support=True, field_modulus=field_modulus,
                reconstructed_binary_rank=len(actual_space), original_binary_rank=len(expected_space),
                original_public_row_space_equal=True, verified_equivalent_goppa_key=True,
                secret_reference_used=False, recovery_modules_imported=False,
                seconds=time.monotonic()-start)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('candidate', type=Path)
    ap.add_argument('--public', type=Path, default=PUBLIC)
    ap.add_argument('--report', type=Path)
    args = ap.parse_args()
    try:
        report = verify(args.candidate, args.public)
    except (ValueError, KeyError, OSError) as exc:
        ap.exit(1, f'Independent verification failed: {exc}\n')
    if args.report:
        with args.report.open('x') as out:
            json.dump(report, out, indent=2)
            out.write('\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
