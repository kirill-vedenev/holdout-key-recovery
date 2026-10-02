"""One recovered E-valued jet -> global quadrics -> aligned tangents -> GRS.

Public inputs only: public.json and deep-jets-POSITION.json. In particular,
the original holdout basis and secret.json are not read by this program.
"""
from sage.all import *
from pathlib import Path
from itertools import combinations, combinations_with_replacement
from collections import Counter
import argparse,json,sys,time
HERE=Path(__file__).resolve().parent
from engine import sidelnikov_shestakov,grs_multiplier_systematic,grs_multiplier_fast,grs_parity


def encode_matrix(A):
    return [[int(x.to_integer()) for x in row] for row in A]


def main():
    ap=argparse.ArgumentParser(); ap.add_argument('directory'); ap.add_argument('--position',type=int,default=28)
    ap.add_argument('--depth',type=int); args=ap.parse_args(); source=Path(args.directory); start=time.monotonic()
    public=json.loads((source/'public.json').read_text())
    recovered=json.loads((source/f'deep-jets-{args.position}.json').read_text())
    parameters=public['parameters']; m=parameters['m']; n=parameters['n']; t=parameters['t']; D=n-2*t-1
    E=GF(2**m,'a',modulus=PolynomialRing(GF(2),'x')(public['field_modulus']))
    decode=lambda x:E.from_integer(x)
    P=matrix(E,public['generator']); k=P.nrows(); set_random_seed(3801)
    jets=[vector(E,[decode(x) for x in row]) for row in recovered['jets']]
    assert jets[0]==P.column(args.position)
    depth=max(0,2*D-n+1) if args.depth is None else args.depth
    if depth<0 or depth>=len(jets): raise ValueError('requested depth not available in saved coherent jets')
    quadratic_monomials=list(combinations_with_replacement(range(k),2))
    R=PolynomialRing(E,'U')
    arcs=[R([jet[a] for jet in jets]) for a in range(k)]
    local_products=[arcs[a]*arcs[b] for a,b in quadratic_monomials]
    values=matrix(E,[[p[a]*p[b] for a,b in quadratic_monomials] for p in P.columns()])
    contacts=matrix(E,len(jets)-1,len(quadratic_monomials),lambda r,c:local_products[c][r+1])
    M=values.stack(contacts[:depth]); Q=M.right_kernel_matrix()
    assert not M*Q.transpose()
    all_saved_contact=not bool(Q*contacts.transpose())
    profile=[dict(depth=r,rank=values.stack(contacts[:r]).rank()) for r in range(depth+1)]
    tangent_spaces=[]; representatives=[]; gradients=[]
    for i,p in enumerate(P.columns()):
        B=matrix(E,len(quadratic_monomials),k)
        for row,(a,b) in enumerate(quadratic_monomials):
            if a!=b: B[row,a]=p[b]; B[row,b]=p[a]
        G=Q*B; W=G.right_kernel(); tangent_spaces.append(W); gradients.append(G.rank())
        assert p in W
        if W.dimension()==2:
            representatives.append(next(v for v in W.basis() if matrix(E,[p,v]).rank()==2))
    report=dict(public_only=True,input_files=['public.json',f'deep-jets-{args.position}.json'],
        parameters=dict(m=m,n=n,t=t,k=k,D=D),seed_position=args.position,
        seed_depth_used=depth,seed_depth_available=len(jets)-1,
        ordinary_interpolation_bound_satisfied=n+depth>2*D,
        coefficient_field_order=E.order(),square_monomials_included=True,
        quadratic_unknowns=len(quadratic_monomials),constraint_rows=M.nrows(),constraint_rank=M.rank(),
        quadratic_ideal_dimension=Q.nrows(),all_saved_jet_contacts_satisfied=all_saved_contact,
        prefix_rank_profile=profile,first_depth_attaining_final_rank=next(r['depth'] for r in profile if r['rank']==M.rank()),
        gradient_rank_counts=dict(Counter(int(x) for x in gradients)),
        tangent_dimension_counts=dict(Counter(int(W.dimension()) for W in tangent_spaces)),
        bootstrap_seconds=time.monotonic()-start)
    payload=dict(public_only=True,field_modulus=public['field_modulus'],quadratic_monomials=quadratic_monomials,
                 quadratic_basis=encode_matrix(Q))
    if len(representatives)==n:
        tick=time.monotonic(); V=matrix(E,representatives).transpose()
        pairs=list(combinations(range(k),2))
        minors=matrix(E,len(pairs),n,lambda r,i:
            (P[pairs[r][0],i]*V[pairs[r][1],i]+P[pairs[r][1],i]*V[pairs[r][0],i])**(2**(m-1)))
        rank=minors.rank(); T=minors.row_space().basis_matrix(); predicted=D-t
        report['minor']=dict(rows=minors.nrows(),length=n,rank=rank,predicted_binary_goppa_dimension=predicted,
                             all_pairs_used=True,construction_seconds=time.monotonic()-tick)
        payload.update(tangent_representatives=encode_matrix(V),minor_basis=encode_matrix(T))
        tick=time.monotonic(); small=T.right_kernel_matrix() if rank>n-rank else T
        support=sidelnikov_shestakov(E,small,tries=E.order())
        recovered_grs=False; embedding=False
        if support is not None:
            multiplier=grs_multiplier_systematic(E,T,support)
            if multiplier is not None:
                grs=matrix(E,rank,n,lambda a,i:multiplier[i]*support[i]**a)
                recovered_grs=T.row_space()==grs.row_space()
                assert recovered_grs
                payload.update(recovered_support=[int(x.to_integer()) for x in support],
                               recovered_minor_multiplier=[int(x.to_integer()) for x in multiplier])
                lam=grs_multiplier_fast(E,P,support,D+1)
                if lam is not None:
                    H=grs_parity(E,support,D+1,lam)
                    embedding=not bool(P*H.transpose())
                    _,_,toV=E.vector_space(map=True)
                    HB=matrix(GF(2),[[toV(x)[bit] for x in row] for row in H for bit in range(m)])
                    subfield=HB.right_kernel()
                    binary_public=matrix(GF(2),public['generator']).row_space()
                    report['recovered_subfield_code_dimension']=subfield.dimension()
                    report['recovered_subfield_code_equals_public_code']=subfield==binary_public
                    payload['recovered_ambient_multiplier']=[int(x.to_integer()) for x in lam]
        report['public_grs_recovery']=dict(support_recovered=support is not None,
            exact_grs_row_space_verified=recovered_grs,ambient_public_embedding_verified=embedding,
            seconds=time.monotonic()-tick)
    else:
        report['minor']=dict(status='not_constructed',reason='some gradient kernels are not two-dimensional')
    report['total_seconds']=time.monotonic()-start
    (source/'bootstrap.json').write_text(json.dumps(report,indent=2,default=int)+'\n')
    (source/'bootstrap-public-data.json').write_text(json.dumps(payload,default=int)+'\n')
    print(json.dumps({a:b for a,b in report.items() if a!='prefix_rank_profile'},default=int),flush=True)


if __name__=='__main__':main()
