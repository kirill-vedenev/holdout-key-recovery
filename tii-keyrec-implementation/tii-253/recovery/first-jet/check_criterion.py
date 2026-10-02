"""Measure the public first-jet and continuation rank conditions on saved forms."""
from sage.all import GF, PolynomialRing, matrix, vector
from pathlib import Path
from itertools import combinations, combinations_with_replacement
from math import comb
import argparse
import hashlib
import json
import time
from gf2_io import panel_matrix
from public_local import public_branches, quadratic_restriction, polar_matrix


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def binary_rank(rows):
    """Independent integer-XOR rank, separate from Sage's matrix elimination."""
    pivots = {}
    for row in rows:
        value = sum(int(x) << i for i, x in enumerate(row))
        while value:
            i = value.bit_length() - 1
            if i not in pivots:
                pivots[i] = value
                break
            value ^= pivots[i]
    return len(pivots)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('public', type=Path)
    ap.add_argument('index', type=Path)
    ap.add_argument('results', type=Path)
    args = ap.parse_args()
    start = time.monotonic()
    data = json.loads(args.public.read_text())
    index = json.loads(args.index.read_text())
    k, m, d = (data['parameters'][x] for x in ('k', 'm', 'degree'))
    count = index['count']
    words = (count + 63) // 64
    assert (k, m, d, count) == (46, 8, 6, 280)
    assert data['public_sha256'] == index['public_sha256']
    tensors = [panel_matrix(args.results / f'hasse-{u}.u64le', count, words) for u in range(3)]
    assert [a.ncols() for a in tensors] == [comb(k, u) for u in range(3)]
    V, J, Q = tensors
    p = vector(GF(2), [(data['heldout_mask'] >> a) & 1 for a in range(k)])
    assert p == matrix(GF(2), data['generator']).column(data['heldout'][0])
    pairs = sorted(combinations(range(k), 2), key=lambda pair: sum(comb(a, j + 1) for j, a in enumerate(pair)))
    E = GF(2**m, 'a', modulus=PolynomialRing(GF(2), 'x')(data['field_modulus']))
    W = J.right_kernel()
    S = J.left_kernel().basis_matrix()
    B = W.basis_matrix()
    R = quadratic_restriction(Q, B, pairs)
    R0 = S * R
    assert J.rank() == binary_rank(J.rows())
    assert R.rank() == binary_rank(R.rows())
    assert R0.rank() == binary_rank(R0.rows())
    assert not J * p
    assert not Q * vector(GF(2), [p[a] * p[b] for a, b in pairs])
    # The point coordinate must disappear in restriction, including all polar terms.
    point_polar = polar_matrix(Q, p, pairs)
    assert not point_polar * B.transpose()
    report = dict(
        challenge=data['challenge'], instance_id=data['instance_id'], public_only=True,
        heldout_original=data['heldout_original'], heldout_shortened=data['heldout'][0],
        heldout_weight=sum(int(x) for x in p), parameters=data['parameters'], kernel_polynomials=count,
        degree=d, field_modulus=data['field_modulus'],
        heldout_value_rank=int(V.rank()), gradient_rank=int(J.rank()),
        expected_gradient_rank=k-m-1, gradient_annihilator_dimension=int(W.dimension()),
        quotient_dimension=int(W.dimension())-1, stationary_polynomial_dimension=int(S.nrows()),
        restricted_quadratic_rank=int(R.rank()), restricted_stationary_quadratic_rank=int(R0.rank()),
        expected_restricted_quadratic_rank=comb(m, 2), independent_binary_rank_checks=True,
        maximal_gradient_rank=bool(J.rank()==k-m-1),
        binary_rank_criterion_holds=bool(not V and J.rank()==k-m-1 and R.rank()==comb(m, 2)),
        stationary_rank_criterion_holds=bool(not V and J.rank()==k-m-1 and R0.rank()==comb(m, 2)),
        gradient_annihilator_basis=[[int(x) for x in v] for v in B.rows()],
        restricted_quadratic_monomials=[list(uv) for uv in combinations_with_replacement(range(W.dimension()), 2)],
        restricted_quadratic_row_basis=[[int(x) for x in v] for v in R.row_space().basis()],
        restricted_stationary_quadratic_row_basis=[[int(x) for x in v] for v in R0.row_space().basis()],
        curve_identity_certified=False,
        curve_identity_scope='Local rank conditions are measured exactly; identifying branches with the hidden curve remains conditional on the forms being curve identities.',
    )
    print(json.dumps({key:value for key,value in report.items() if 'basis' not in key}), flush=True)
    branches, detail = public_branches(J, Q, p, E, pairs)
    report['branch_solver'] = detail
    report['branch_vectors'] = [[int(x.to_integer()) for x in v] for v in branches]
    report['continuation_checks'] = []
    JE = J.change_ring(E)
    SE = S.change_ring(E)
    for branch_id, v in enumerate(branches):
        L = polar_matrix(Q, v, pairs)
        A = JE.stack(L)
        A0 = JE.stack(SE * L)
        plane = matrix(E, [p.change_ring(E), v]).row_space()
        report['continuation_checks'].append(dict(
            branch=branch_id, binary_rank=int(A.rank()), binary_nullity=k-int(A.rank()),
            binary_kernel_equals_point_tangent_plane=bool(A.right_kernel()==plane),
            stationary_rank=int(A0.rank()), stationary_nullity=k-int(A0.rank()),
            stationary_kernel_equals_point_tangent_plane=bool(A0.right_kernel()==plane)))
    if branches:
        planes = [matrix(E, [p.change_ring(E), v]).row_space() for v in branches]
        frobenius_permutation = []
        for v in branches:
            square_plane = matrix(E, [p.change_ring(E), vector(E, [x**2 for x in v])]).row_space()
            matches = [i for i, plane in enumerate(planes) if plane == square_plane]
            assert len(matches) == 1
            frobenius_permutation.append(matches[0])
        report['frobenius_branch_permutation'] = frobenius_permutation
        orbit = []; current = 0
        while current not in orbit:
            orbit.append(current); current = frobenius_permutation[current]
        report['one_full_frobenius_orbit'] = bool(current==0 and len(orbit)==m and len(branches)==m)
        report['point_and_branches_span_gradient_kernel'] = bool(matrix(E, [p.change_ring(E)] + branches).row_space()==B.change_ring(E).row_space())
    report['all_binary_continuation_checks_pass'] = bool(len(branches)==m and all(x['binary_rank']==k-2 and x['binary_kernel_equals_point_tangent_plane'] for x in report['continuation_checks']))
    report['all_stationary_continuation_checks_pass'] = bool(len(branches)==m and all(x['stationary_rank']==k-2 and x['stationary_kernel_equals_point_tangent_plane'] for x in report['continuation_checks']))
    D = data['parameters']['n']-2*data['parameters']['t']-1
    forced = sum(data['multiplicities']) + data['infinity_multiplicity']
    observed_extra = int(not V and data['multiplicities'][data['heldout'][0]]==0)
    shortfall = d*D+1-forced
    report['zero_count'] = dict(curve_degree_bound=D,
        composition_degree_bound=d*D, forced_zeros=forced, strict_sufficient_count=d*D+1,
        strict_zero_count_certificate=bool(forced>d*D), weis_shortfall=shortfall,
        goppa_degree=data['parameters']['t'],
        weis_hidden_constraint_heuristic_applies=bool(shortfall<=data['parameters']['t']),
        observed_heldout_value_additional_zero=observed_extra,
        zeros_including_observed_heldout_value=forced+observed_extra,
        strict_certificate_including_heldout_value=bool(forced+observed_extra>d*D))
    report['input_sha256'] = {str(path):sha(path) for path in [args.public,args.index] + [args.results / f'hasse-{u}.u64le' for u in range(3)]}
    report['seconds'] = time.monotonic()-start
    (args.results/'criterion.json').write_text(json.dumps(report, indent=2, default=int)+'\n')
    print(json.dumps({key:value for key,value in report.items() if 'basis' not in key and key!='branch_vectors'}, default=int), flush=True)


if __name__ == '__main__':
    main()
