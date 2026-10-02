"""Certify global quadrics, aligned tangents, and the public point at infinity."""
from sage.all import *
from pathlib import Path
import argparse,json,time

ap=argparse.ArgumentParser();ap.add_argument('results',type=Path);args=ap.parse_args()
out=args.results;start=time.monotonic();public=json.loads((out/'public.json').read_text())
support=json.loads((out/'retained-support.json').read_text());payload=json.loads((out/'bootstrap-public-data.json').read_text())
coherence=json.loads((out/'coherence-audit.json').read_text());jets=json.loads((out/'deep-jets-0.json').read_text())
substitution=json.loads((out/'full-jet-substitution-verification.json').read_text())
m,n,t,k=[public['parameters'][a] for a in ['m','n','t','k']];D=n-2*t-1
E=GF(2**m,'a',modulus=PolynomialRing(GF(2),'x')(public['field_modulus']));dec=E.from_integer;enc=lambda a:int(a.to_integer())
R=PolynomialRing(E,'x');x=R.gen();F=[R([dec(a) for a in f]) for f in coherence['curve']]
P=matrix(E,public['generator']);alpha=[dec(a) for a in support['support']]
Q=matrix(E,[[dec(a) for a in row] for row in payload['quadratic_basis']]);pairs=payload['quadratic_monomials']
products=[F[a]*F[b] for a,b in pairs];coeff=matrix(E,2*D+1,len(pairs),lambda r,j:products[j][r])
assert not Q*coeff.transpose()
assert Q.nrows()+coeff.rank()==len(pairs) and Q.rank()==Q.nrows()
V=matrix(E,[[dec(a) for a in row] for row in payload['tangent_representatives']])
for j in range(n):
    tangent=vector(E,[f.derivative()(alpha[j]) for f in F])
    assert matrix(E,[P.column(j),tangent]).rank()==2
    assert matrix(E,[P.column(j),tangent]).row_space()==matrix(E,[P.column(j),V.column(j)]).row_space()
pinfinity=sum(P.columns(),vector(E,k));assert pinfinity==vector(E,[(public['infinity_mask']>>a)&1 for a in range(k)])
anchor=next(a for a,z in enumerate(pinfinity) if z);matches=[]
for parameter in [dec(a) for a in range(E.order())]+[None]:
    value=vector(E,[f[D] if parameter is None else f(parameter) for f in F])
    if value[anchor] and value==value[anchor]/pinfinity[anchor]*pinfinity:matches.append(parameter)
assert len(matches)==1 and matches[0] not in alpha
infinity=dict(unique_projective_preimage=True,projective_candidates_checked=E.order()+1,
    parameter=None if matches[0] is None else enc(matches[0]),public_point_mask=public['infinity_mask'],
    method='exact projective evaluation of the reconstructed degree-99 curve')
support['public_infinity']=infinity;(out/'retained-support.json').write_text(json.dumps(support,indent=2,default=int)+'\n')
forced=sum(public['multiplicities'])+public['infinity_multiplicity'];total=forced+jets['depth']+1
assert coherence['report']['coherent_scalar_and_parameter_series_verified']
assert substitution['all_residuals_zero'] and substitution['jet_depth']==jets['depth']
assert total>public['parameters']['degree']*D
report=dict(status='passed',public_only=True,global_quadric_identities_checked=Q.nrows(),
    quadratic_identity_coefficient_orders=2*D+1,full_quadratic_ideal_certified=True,
    quadratic_product_rank=coeff.rank(),aligned_tangent_planes_verified=n,
    public_infinity=infinity,coherent_seed_jet_depth=jets['depth'],
    holdout_curve_identity_certificate=dict(certified=True,composition_degree_bound=public['parameters']['degree']*D,
        preexisting_forced_zeros=forced,heldout_multiplicity_from_verified_coherent_jet=jets['depth']+1,
        total_zeros_counted_with_multiplicity=total,polynomials_certified=jets['kernel_dimension'],
        justification='Full jet substitution and one coherent invertible local parametrization add held-out multiplicity to the already verified public Hasse constraints.'),
    seconds=time.monotonic()-start)
(out/'geometry-verification.json').write_text(json.dumps(report,indent=2,default=int)+'\n')
print(json.dumps(report,default=int),flush=True)
