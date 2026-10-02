"""Direct characteristic-two recovery from one coherent reparametrized jet.

recover_curve takes only the coefficient field, degree bound, and jet vectors.
It neither constructs a curve ideal nor reads public columns or a hidden key.
"""
from sage.all import *
from itertools import combinations_with_replacement
from pathlib import Path
import argparse
import json
import time


def recover_curve(E, D, jets, form_degree=2):
    start = time.monotonic()
    q = int(form_degree)
    if q < 2 or q & (q-1):
        raise ValueError('form degree must be a positive power of two, at least two')
    if E.characteristic() != 2 or D < 1 or len(jets) < 2*q*D+1:
        raise ValueError('requires characteristic two, D >= 1, and orders 0..2qD')
    k = len(jets[0]); depth = 2*q*D
    R = PolynomialRing(E, 'U'); U = R.gen()
    pairs = list(combinations_with_replacement(range(k), q))
    arcs = [R([row[a] for row in jets[:depth+1]]) for a in range(k)]
    products = [prod(arcs[a] for a in monomial) for monomial in pairs]
    contact = matrix(E, q*D+1, len(pairs), lambda r,j: products[j][r])
    targets = matrix(E, q*D+1, 2,
                     lambda r,j: int(r == (q*D if j == 0 else q*(D-1))))
    try:
        solutions = contact.solve_right(targets)
    except ValueError:
        return None, dict(status='extreme_sections_absent', D=D, k=k,form_degree=q,
                          product_rank=contact.rank(), seconds=time.monotonic()-start)
    assert contact*solutions == targets
    selected = [sum((solutions[j,b]*products[j] for j in range(len(pairs))), R.zero())
         for b in range(2)]
    if any(f[r] for f in selected for r in range(depth+1) if r%q):
        raise ValueError('selected sections fail the power-series consequence')
    exponent = 2**((-(q.bit_length()-1)) % E.degree())
    root = lambda x: x**exponent
    roots = [R([root(f[q*r]) for r in range(2*D+1)]) for f in selected]
    assert [f.valuation() for f in roots] == [D,D-1]
    PS = PowerSeriesRing(E, 'u', default_prec=D+1); u = PS.gen()
    high_unit = PS(roots[0] // U**D).add_bigoh(D+1)
    lower_unit = PS(roots[1] // U**(D-1)).add_bigoh(D+1)
    parameter = (u*high_unit/lower_unit).add_bigoh(D+1)
    scale = (lower_unit**D / high_unit**(D-1)).add_bigoh(D+1)
    powers = [parameter**r for r in range(D+1)]
    change = matrix(E,D+1,D+1,lambda r,j:powers[j][r])
    normalized = [PS(f)/scale for f in arcs]
    values = matrix(E,D+1,k,lambda r,j:normalized[j][r])
    coefficients = change.solve_right(values)
    curve = [R(list(coefficients.column(j))) for j in range(k)]
    assert max(f.degree() for f in curve) <= D
    # Return the polynomial curve and the two explicit evaluation forms.
    result = dict(curve=curve, forms=solutions.transpose(), pairs=pairs,form_degree=q,
                  parameter=parameter, scale=scale)
    report = dict(status='recovered',D=D,k=k,form_degree=q,jet_order_for_sections=q*D,
                  jet_order_for_curve=2*q*D,form_unknowns=len(pairs),
                  contact_rows=q*D+1,product_rank=contact.rank(),
                  perfect_power_series_check=True,seconds=time.monotonic()-start)
    if q==2:
        odd=contact.matrix_from_rows(range(1,2*D+1,2))
        report['square_section_dimension']=contact.rank()-odd.rank()
    return result,report


def support_from_sections(E, P, seed_position, forms, pairs,form_degree=2):
    exponent=2**((-(int(form_degree).bit_length()-1)) % E.degree())
    root = lambda x:x**exponent
    support = []
    for i,p in enumerate(P.columns()):
        if i == seed_position:
            support.append(E.zero())
            continue
        monomials = vector(E,[prod(p[a] for a in monomial) for monomial in pairs])
        high,lower = forms*monomials
        if not high:
            raise ValueError('high section vanishes away from the seed')
        support.append(root(high/lower) if lower else None)
    if len(set(support)) != len(support):
        raise ValueError('recovered support is not distinct')
    return support


def validate_public(E,P,D,curve,support,subfield=True):
    multipliers = []
    for p,z in zip(P.columns(),support):
        value = vector(E,[f(z) if z is not None else f[D] for f in curve])
        anchor = next(j for j,x in enumerate(value) if x)
        lam = p[anchor]/value[anchor]
        assert lam and lam*value == p
        multipliers.append(lam)
    report = dict(all_public_columns_matched=True,distinct_support=True,
                  infinite_support_entries=support.count(None))
    if subfield:
        ambient = matrix(E,D+1,P.ncols(),lambda r,j:
            multipliers[j]*(int(r==D) if support[j] is None else support[j]**r))
        H = ambient.right_kernel_matrix()
        m = E.degree()
        binary_check = matrix(GF(2),H.nrows()*m,H.ncols(),lambda r,j:
                              (int(H[r//m,j].to_integer()) >> (r%m)) & 1)
        public_binary = matrix(GF(2),P)
        equal = binary_check.right_kernel() == public_binary.row_space()
        assert equal
        report.update(ambient_grs_dimension=ambient.rank(),
                      recovered_binary_subfield_code_equals_public=True)
    return multipliers,report


def encode_result(result,support,multipliers,report):
    enc = lambda x:int(x.to_integer())
    return dict(report=report,curve=[[enc(x) for x in f.list()] for f in result['curve']],
                forms=[[enc(x) for x in row] for row in result['forms']],
                support=[None if x is None else enc(x) for x in support],
                multipliers=[enc(x) for x in multipliers])


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('directory',type=Path)
    ap.add_argument('--jets',type=Path)
    ap.add_argument('--position',type=int,default=28)
    ap.add_argument('--form-degree',type=int,default=2)
    ap.add_argument('--output',type=Path,required=True)
    args = ap.parse_args()
    public_file = args.directory/'public.json'
    jet_file = args.jets or args.directory/f'deep-jets-{args.position}.json'
    public = json.loads(public_file.read_text()); saved = json.loads(jet_file.read_text())
    m,n,t,k = [public['parameters'][x] for x in ['m','n','t','k']]; D=n-2*t-1
    E = GF(2**m,'a',modulus=PolynomialRing(GF(2),'x')(public['field_modulus']))
    jets = [vector(E,[E.from_integer(x) for x in row]) for row in saved['jets']]
    result,report = recover_curve(E,D,jets,args.form_degree)
    report.update(public_only=True,instance=args.directory.name,position=args.position,
                  input_files=[str(public_file.resolve()),str(jet_file.resolve())],
                  parameters=dict(m=m,n=n,t=t,k=k,D=D))
    if result is None:
        args.output.write_text(json.dumps(dict(report=report),indent=2)+'\n')
    else:
        # The columns enter only after the polynomial curve has been recovered.
        P = matrix(E,public['generator'])
        support = support_from_sections(E,P,args.position,result['forms'],result['pairs'],args.form_degree)
        multipliers,checks = validate_public(E,P,D,result['curve'],support)
        report.update(checks)
        args.output.write_text(json.dumps(encode_result(result,support,multipliers,report),indent=2)+'\n')
    print(json.dumps(report,default=int),flush=True)


if __name__ == '__main__':
    main()
