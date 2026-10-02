"""Audit the entire public recovery chain using only its recovered key.

The curve is rebuilt from the original public shortened generator and the
newly recovered Goppa key. No supplied/published secret is an input.
"""
from sage.all import *
from pathlib import Path
from itertools import combinations, combinations_with_replacement
import argparse, hashlib, json, time
from core import read_kernel, jets_at_point, products, apply_binary
from public_local import public_branches, polar_matrix
from direct_jet_recovery import support_from_sections


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('directory',type=Path)
    args=ap.parse_args(); folder=args.directory; start=time.monotonic()
    public=json.loads((folder/'public.json').read_text())
    recovered=json.loads((folder/'deep-jets-0.json').read_text())
    direct=json.loads((folder/'direct-recovery.json').read_text())
    key=json.loads((folder/'equivalent-key.json').read_text())
    E=GF(2**public['parameters']['m'],'a',modulus=PolynomialRing(GF(2),'x')(public['field_modulus']))
    dec=E.from_integer; enc=lambda x:int(x.to_integer())
    R=PolynomialRing(E,'U');U=R.gen()
    P=matrix(E,public['generator']);k=P.nrows();i=public['heldout'][0]
    alpha=[dec(key['support'][j]) for j in public['original_indices']]
    g=R([dec(x) for x in key['g']]);Pi=prod(U-a for a in alpha)
    partial=[Pi//(U-a) for a in alpha]
    raw=[sum((P[a,j]*partial[j] for j in range(P.ncols())),R.zero()) for a in range(k)]
    F=[]
    for f in raw:
        q,r=f.quo_rem(g*g)
        assert not r
        F.append(q)
    D=P.ncols()-2*g.degree()-1
    assert max(f.degree() for f in F)==D and gcd(F).degree()==0
    for j,a in enumerate(alpha):
        assert g(a)**2/Pi.derivative()(a)*vector(E,[f(a) for f in F])==P.column(j)
    print(json.dumps(dict(event='recovered_key_curve_rebuilt',degree=D),default=int),flush=True)

    # A complete coefficient check, not evaluation merely on the support.
    kk,d,K=read_kernel(folder/'kernel.bin');assert kk==k
    composition=apply_binary(K,products([f.list() for f in F],d,E,d*D),E)
    assert not composition
    identity_seconds=time.monotonic()-start
    print(json.dumps(dict(event='all_kernel_curve_identities_verified',forms=K.nrows(),
                          composition_degree=d*D,seconds=identity_seconds),default=int),flush=True)

    jets=[vector(E,[dec(x) for x in row]) for row in recovered['jets']]
    depth=recovered['depth'];p=P.column(i)
    coordinate=next(a for a,x in enumerate(p) if x)
    scalar=p[coordinate]/F[coordinate](alpha[i])
    local=[scalar*f(alpha[i]+U) for f in F]
    C=matrix(E,k,D+1,lambda a,r:local[a][r])
    first=C.matrix_from_columns([0,1]);assert first.rank()==2
    assert jets[1] in first.column_space()

    # Compare the exhaustive public branch locus to the seven Frobenius
    # tangent planes of the independently rebuilt curve.
    local_maps=K*jets_at_point(k,d,p).transpose()
    J=local_maps.matrix_from_columns(range(1,k+1))
    Q=local_maps.matrix_from_columns(range(k+1,local_maps.ncols()))
    pairs=list(combinations(range(k),2))
    branches,branch_report=public_branches(J,Q,vector(GF(2),[enc(x) for x in p]),E,pairs)
    found=[matrix(E,[p,v]).row_space() for v in branches]
    expected=[matrix(E,[p,vector(E,[x**(2**s) for x in first.column(1)])]).row_space()
              for s in range(E.degree())]
    assert len(found)==E.degree() and all(w in expected for w in found) and all(w in found for w in expected)
    continuation_ranks=[]
    for v in branches:
        A=J.change_ring(E).stack(polar_matrix(Q,v,pairs))
        assert A.right_kernel()==matrix(E,[p,v]).row_space()
        continuation_ranks.append(A.rank())

    # Fit one scalar and one invertible parameter series to all jet orders.
    pivot_rows=first.transpose().pivots(); inverse=first.matrix_from_rows(pivot_rows).inverse()
    scale=[E.one()]+[E.zero()]*depth; parameter=[E.zero()]*(depth+1)
    for r in range(1,depth+1):
        T=PowerSeriesRing(E,'z',default_prec=r+1)
        lam=T(scale[:r+1],prec=r+1);psi=T(parameter[:r+1],prec=r+1)
        power=T.one();coefficients=[]
        for j in range(D+1):
            coefficients.append((lam*power)[r]);power*=psi
        residual=jets[r]-C*vector(E,coefficients)
        correction=inverse*vector(E,[residual[a] for a in pivot_rows])
        assert first*correction==residual, ('incoherent recovered jet',r)
        scale[r],parameter[r]=correction
    assert parameter[1] and scale[0]
    T=PowerSeriesRing(E,'z',default_prec=depth+1)
    lam=T(scale,prec=depth+1);psi=T(parameter,prec=depth+1)
    powers=[T.one()]
    for j in range(D):powers.append(powers[-1]*psi)
    for a in range(k):
        composed=lam*sum((C[a,j]*powers[j] for j in range(D+1)),T.zero())
        assert all(composed[r]==jets[r][a] for r in range(depth+1))

    # Discard orders > 2D and independently repeat support-only extraction.
    arcs=[R([v[a] for v in jets[:2*D+1]]) for a in range(k)]
    qpairs=list(combinations_with_replacement(range(k),2))
    prods=[arcs[a]*arcs[b] for a,b in qpairs]
    contact=matrix(E,2*D+1,len(qpairs),lambda r,j:prods[j][r])
    target=matrix(E,2*D+1,2,lambda r,j:int(r==(2*D if j==0 else 2*D-2)))
    forms=contact.solve_right(target).transpose()
    support=support_from_sections(E,P,i,forms,qpairs)
    saved=[None if x is None else dec(x) for x in direct['support']]
    assert support==saved

    report=dict(public_only=True,published_secret_key_read=False,
        input_files=['public.json','kernel.bin','deep-jets-0.json','direct-recovery.json','equivalent-key.json'],
        kernel_sha256=hashlib.sha256((folder/'kernel.bin').read_bytes()).hexdigest(),
        recovered_key_sha256=hashlib.sha256((folder/'equivalent-key.json').read_bytes()).hexdigest(),
        curve_rebuilt_from_public_generator_and_recovered_goppa_key=True,
        curve_degree=D,common_coordinate_factor_degree=0,all_public_curve_columns_verified=True,
        kernel_forms_checked=K.nrows(),composition_degree_bound=d*D,
        coefficients_checked_per_kernel_form=d*D+1,
        all_kernel_forms_exact_polynomial_curve_identities=True,
        branch_report=branch_report,all_public_branches_equal_recovered_key_frobenius_branches=True,
        all_branch_continuation_ranks=continuation_ranks,
        continuation_kernel_equals_tangent_plane_on_every_branch=True,
        jet_depth=depth,one_coherent_scalar_and_reparametrization=True,
        exact_jet_coefficient_equalities=k*(depth+1),
        scalar_constant_nonzero=True,parameter_linear_coefficient_nonzero=True,
        support_only_jet_order_used=2*D,later_coefficients_discarded_for_support_check=True,
        support_only_result_matches_full_curve_recovery=True,
        curve_identity_seconds=identity_seconds,seconds=time.monotonic()-start,
        scalar_series=[enc(x) for x in scale],parameter_series=[enc(x) for x in parameter])
    (folder/'chain-verification.json').write_text(json.dumps(report,indent=2,default=int)+'\n')
    print(json.dumps({a:b for a,b in report.items() if a not in ('scalar_series','parameter_series')},default=int),flush=True)


if __name__=='__main__':main()
