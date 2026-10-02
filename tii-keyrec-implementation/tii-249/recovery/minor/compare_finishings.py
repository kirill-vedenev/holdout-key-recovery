"""Compare a saved minor-finishing support with the saved direct full key.

The field and Moebius calculations are retained from the archived comparison.
Explicit paths replace the original fixed two-instance directory layout.
"""
from sage.all import *
from pathlib import Path
import argparse
import json

def field(modulus):
    m = len(modulus) - 1
    return GF(2**m, 'a', modulus=PolynomialRing(GF(2), 'x')(modulus))


def moebius_from(src, dst):
    # solve dst = (a s + b)/(c s + d) from three points; return (a,b,c,d)
    E = src[0].parent()
    M = matrix(E, [[s, 1, -s*d, -d] for s, d in zip(src, dst)])
    K = M.right_kernel()
    return K.basis()[0] if K.dimension() == 1 else None



def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--public', type=Path, required=True)
    parser.add_argument('--minor', type=Path, required=True)
    parser.add_argument('--direct-key', type=Path, required=True)
    parser.add_argument('--report', type=Path, required=True)
    args = parser.parse_args()
    pub = json.loads(args.public.read_text())
    r2 = json.loads(args.minor.read_text())
    key = json.loads(args.direct_key.read_text())
    E = field(pub['field_modulus']); m = E.degree(); dec = E.from_integer
    assert pub['field_modulus'] == key['field_modulus'] == r2['field_modulus']
    orig = pub['original_indices']
    s2 = [dec(x) for x in r2['recovered_support']]
    s1 = [dec(key['support'][o]) for o in orig]
    assert len(s2) == len(s1) == len(orig)
    found = []
    for sigma in range(m):
        c1 = [x**(2**sigma) for x in s1]
        coeffs = moebius_from(c1[:3], s2[:3])
        if coeffs is None: continue
        a, b, c, d = coeffs
        if all(c*x + d != 0 and (a*x + b)/(c*x + d) == y for x, y in zip(c1, s2)):
            found.append(sigma)
    result = dict(retained_positions=len(orig), route2_support_distinct=len(set(s2)) == len(s2),
                  frobenius_powers_matching=found, moebius_equivalent=bool(found))
    if not result['route2_support_distinct'] or not result['moebius_equivalent']:
        raise ValueError('Saved finishing routes do not have equivalent distinct supports')
    report = {pub['challenge'].lower(): result}
    with args.report.open('x') as stream:
        json.dump(report, stream, indent=2)
        stream.write('\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
