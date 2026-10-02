"""Public-only IKZ-type support extension, multipliers, and Goppa key recovery.

Input: the recovered shortened parametrization and the original public matrix.
No published private key, supplied support, or challenge answer is read.
"""
from sage.all import *
from pathlib import Path
from itertools import combinations_with_replacement
import argparse, hashlib, json, random, time


def encode(x):
    return int(x.to_integer())


def dump(path, value):
    path.write_text(json.dumps(value, indent=2, default=int) + '\n')


def expand_binary(A):
    m = A.base_ring().degree()
    return matrix(GF(2), A.nrows()*m, A.ncols(),
                  lambda r,c: (encode(A[r//m,c]) >> (r % m)) & 1)


def multiplier_matrix(G, alpha, codimension):
    E = alpha[0].parent()
    return matrix(E, G.nrows()*codimension, len(alpha),
                  lambda row,j: G[row//codimension,j]*alpha[j]**(row % codimension))


def full_support_vector(B, seed=173):
    """Return a full-support kernel vector, or an exact impossibility witness.

    A subspace over GF(128) with length at most 96 contains a full-support
    vector iff none of its coordinate functionals vanishes identically:
    fewer than |E| proper hyperplanes cannot cover the subspace.
    Randomization here only finds the vector; it does not decide rejection.
    """
    E = B.base_ring()
    if not B.nrows():
        return None
    if any(not B.column(j) for j in range(B.ncols())):
        return None
    if B.nrows() == 1:
        return B.row(0)
    if B.ncols() >= E.order():
        raise ValueError('full-support existence certificate requires length < field size')
    rng = random.Random(seed)
    for _ in range(10000):
        v = vector(E, [E.from_integer(rng.randrange(E.order())) for _ in range(B.nrows())])*B
        if all(v):
            return v
    raise RuntimeError('full-support vector is certified to exist, but deterministic sampling failed')


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('directory', type=Path)
    ap.add_argument('--public-key', type=Path, required=True)
    args = ap.parse_args()
    start = time.monotonic(); folder = args.directory
    public = json.loads((folder/'public.json').read_text())
    direct = json.loads((folder/'direct-recovery.json').read_text())
    params = public['original_parameters']; n,m,t = [params[s] for s in ['n','m','t']]
    E = GF(2**m, 'a', modulus=PolynomialRing(GF(2),'x')(public['field_modulus']))
    dec = E.from_integer
    R = PolynomialRing(E,'X'); X=R.gen()
    P = matrix(E,public['generator']); k=P.nrows(); D=P.ncols()-2*t-1
    keep = public['original_indices']; removed=public['shortened_original_indices']
    lines = [s.strip() for s in args.public_key.read_text().splitlines() if s.strip()]
    if hashlib.sha256(args.public_key.read_bytes()).hexdigest() != public['public_sha256']:
        raise ValueError('public input SHA-256 mismatch')
    H = matrix(GF(2), [[int(x) for x in line[1:-1].split()] for line in lines[:-1]])
    assert H.ncols()==n and H.rank()==n-params['k']
    assert json.loads(lines[-1])==public['field_modulus']
    Gfull=H.right_kernel_matrix()
    assert H.matrix_from_columns(keep).right_kernel()==P.change_ring(GF(2)).row_space()

    # Direct recovery leaves a PGL(2) freedom. The public sum of columns is
    # the original point at infinity. Sending its recovered parameter to
    # infinity gives a finite affine support and the Goppa polynomial chart.
    forms=matrix(E,[[dec(x) for x in row] for row in direct['forms']])
    pairs=list(combinations_with_replacement(range(k),2))
    pinf=sum(P.columns(),vector(E,k))
    assert pinf==vector(E,[(public['infinity_mask']>>j)&1 for j in range(k)])
    high,lower=forms*vector(E,[pinf[a]*pinf[b] for a,b in pairs])
    assert high
    beta_infinity=(high/lower)**(2**(m-1)) if lower else None
    beta=[None if a is None else dec(a) for a in direct['support']]
    assert beta_infinity not in beta
    def to_affine(z):
        if beta_infinity is None:
            assert z is not None
            return z
        return E.zero() if z is None else 1/(z+beta_infinity)
    alpha=[to_affine(a) for a in beta]
    assert len(set(alpha))==len(alpha) and len(alpha)==len(keep)
    chart_report=dict(public_infinity_parameter=None if beta_infinity is None else encode(beta_infinity),
                      map='identity' if beta_infinity is None else 'x=1/(beta+beta_infinity)',
                      all_retained_support_finite=True,
                      retained_support=[encode(a) for a in alpha])
    dump(folder/'affine-chart.json',chart_report)

    # Each candidate is accepted exactly by the existence of a full-support
    # solution to G_j diag(nu) V_(2t)(alpha,w)^T=0.
    support=[None]*n
    for j,a in zip(keep,alpha): support[j]=a
    candidates=[a for a in E if a not in set(alpha)]
    histories=[]
    for j in removed:
        positions=keep+[j]
        Gj=H.matrix_from_columns(positions).right_kernel_matrix().change_ring(E)
        assert Gj.nrows()==P.nrows()+1
        base=multiplier_matrix(Gj.matrix_from_columns(range(len(keep))),alpha,2*t)
        base_rank=base.rank()
        # In this instance base has full column rank. Factoring once makes
        # exhaustive support testing cheap and still checks every equation.
        if base_rank!=len(keep):
            raise ValueError('unexpected non-unique restricted multiplier: use a general kernel solve')
        rows=base.transpose().pivots()
        inv=base.matrix_from_rows(rows).inverse()
        accepted=[]; tests=[]
        for w in candidates:
            column=vector(E,[Gj[r//(2*t),-1]*w**(r%(2*t)) for r in range(Gj.nrows()*2*t)])
            v=inv*vector(E,[column[r] for r in rows])
            consistent=base*v==column
            full_support=consistent and all(v)
            tests.append(dict(value=encode(w),kernel_dimension=int(consistent),
                              has_full_support_vector=bool(full_support)))
            if full_support: accepted.append(w)
        assert len(accepted)==1, (j, accepted)
        support[j]=accepted[0]
        history=dict(original_position=j,candidates_tested=len(tests),
                     known_multiplier_matrix_rank=base_rank,
                     retained_support_value=encode(accepted[0]),unique_candidate=True,tests=tests)
        histories.append(history)
        print(json.dumps(dict(event='deshortened_position',position=j,
                              support_value=encode(accepted[0]),candidates_tested=len(tests),
                              seconds=time.monotonic()-start)),flush=True)
    assert all(a is not None for a in support) and len(set(support))==n

    A=multiplier_matrix(Gfull.change_ring(E),support,2*t)
    M=A.right_kernel_matrix()
    nu=full_support_vector(M)
    assert nu is not None and not A*nu
    nu=nu/nu[0]
    Pi=prod(X-a for a in support); Pip=Pi.derivative()
    lam=vector(E,[1/(nu[j]*Pip(support[j])) for j in range(n)])
    parity=matrix(E,2*t,n,lambda r,j:nu[j]*support[j]**r)
    ambient=matrix(E,n-2*t,n,lambda r,j:lam[j]*support[j]**r)
    assert not Gfull.change_ring(E)*parity.transpose()
    assert not ambient*parity.transpose()
    assert ambient.rank()==n-2*t and parity.rank()==2*t
    assert expand_binary(parity).row_space()==H.row_space()

    # In the recovered original affine chart reciprocal parity multipliers
    # are a scalar multiple of g(x)^2. Interpolate, then square-root exactly.
    square=R.lagrange_polynomial([(a,1/b) for a,b in zip(support,nu)])
    assert square.degree()==2*t and not square.derivative()
    g=R([square[2*r]**(2**(m-1)) for r in range(t+1)]).monic()
    assert g.degree()==t and gcd(g,g.derivative())==1
    assert all(g(a) for a in support)
    scalar=square.leading_coefficient()
    assert square==scalar*g*g
    gp=matrix(E,t,n,lambda r,j:support[j]**r/g(support[j]))
    assert expand_binary(gp).row_space()==H.row_space()
    assert not Gfull.change_ring(E)*gp.transpose()
    report=dict(public_only=True,original_parameters=params,
                input_files=['public.json','direct-recovery.json',str(args.public_key.resolve())],
                field_modulus=public['field_modulus'],all_96_support_values_distinct=True,
                all_96_support_values_finite=True,removed_positions=removed,
                support_extension='IKZ-type exhaustive Vandermonde multiplier test',
                support_extension_histories=histories,full_multiplier_space_dimension=M.nrows(),
                final_multiplier_system_shape=[A.nrows(),A.ncols()],
                recovered_ambient_grs_dimension=ambient.rank(),original_code_dimension=Gfull.nrows(),
                ambient_grs_contains_original_public_code=True,
                ambient_binary_subfield_code_equals_original_public_code=True,
                reciprocal_multiplier_polynomial_degree=square.degree(),
                reciprocal_multiplier_polynomial_is_square=True,
                goppa_degree=g.degree(),goppa_squarefree=True,goppa_irreducible=g.is_irreducible(),
                goppa_nonzero_on_full_support=True,goppa_parity_binary_rank=expand_binary(gp).rank(),
                goppa_public_row_space_exactly_equal=True,
                recovered_full_equivalent_goppa_decoding_key=True,seconds=time.monotonic()-start)
    dump(folder/'deshortening-report.json',report)
    dump(folder/'equivalent-key.json',dict(challenge='TII-173',public_only_recovery=True,
        field_modulus=public['field_modulus'],support=[encode(x) for x in support],
        g=[encode(x) for x in g.list()],
        dual_grs_multipliers=[encode(x) for x in nu],
        ambient_grs_multipliers=[encode(x) for x in lam],
        public_key_sha256=public['public_sha256'],
        note='Equivalent full degree-8 Goppa key in original column order; not the published secret key.'))
    print(json.dumps({k:v for k,v in report.items() if k!='support_extension_histories'},default=int),flush=True)


if __name__=='__main__': main()
