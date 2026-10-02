"""Public algebraic reconstruction and an independent coherent-jet audit."""
from sage.all import *
from pathlib import Path
import argparse,json,time


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('public',type=Path);ap.add_argument('result',type=Path)
    args=ap.parse_args();begin=time.monotonic();public=json.loads(args.public.read_text());candidate=json.loads((args.result/'retained-support.json').read_text());saved=json.loads((args.result/'deep-jets-0.json').read_text());branches=json.loads((args.result/'local-branches.json').read_text())
    m,n,t,k=[public['parameters'][x] for x in ['m','n','t','k']];D=n-2*t-1
    E=GF(2**m,'a',modulus=PolynomialRing(GF(2),'x')(public['field_modulus']));dec=E.from_integer;enc=lambda x:int(x.to_integer())
    if any(a is None for a in candidate['support']):raise ValueError('this independent affine curve audit requires finite support')
    alpha=[dec(a) for a in candidate['support']];nu=[dec(x) for x in candidate['dual_multipliers']];P=matrix(E,public['generator'])
    R=PolynomialRing(E,'x');x=R.gen();Pi=prod(x-a for a in alpha);partial=[Pi//(x-a) for a in alpha]
    F=[sum((P[a,j]*nu[j]*partial[j] for j in range(n)),R.zero()) for a in range(k)]
    assert max(f.degree() for f in F)<=D and gcd(F).degree()==0
    multipliers=[1/(nu[j]*Pi.derivative()(alpha[j])) for j in range(n)]
    assert all(multipliers[j]*vector(E,[f(alpha[j]) for f in F])==P.column(j) for j in range(n))
    seed=public['heldout'][0];local=[f(alpha[seed]+x) for f in F]
    p=P.column(seed);tangent=vector(E,[f[1] for f in local]);expected=[matrix(E,[p,vector(E,[a**(2**s) for a in tangent])]).row_space() for s in range(m)]
    found=[matrix(E,[p,vector(E,[dec(a) for a in v])]).row_space() for v in branches['branch_vectors']]
    assert len(set(str(V.basis_matrix()) for V in expected))==m and len(found)==m and all(V in expected for V in found)
    permutation=[expected.index(V) for V in found]
    jets=[vector(E,[dec(a) for a in row]) for row in saved['jets']];depth=saved['depth'];anchor=next(a for a,z in enumerate(p) if z)
    other=next(a for a in range(k) if local[a][1]*local[anchor][0]+local[a][0]*local[anchor][1])
    PS=PowerSeriesRing(E,'u',default_prec=depth+1)
    arcs=[PS([v[a] for v in jets],prec=depth+1) for a in range(k)];ratio=(arcs[other]/arcs[anchor]).add_bigoh(depth+1)
    psi=PS.zero();prec=1
    # Newton inversion of a ratio with nonzero derivative determines one
    # common parameter series; the anchor determines the common scalar.
    while prec<depth+1:
        prec=min(2*prec,depth+1);psi=PS(psi.list(),prec=prec);ratio_part=PS(ratio.list()[:prec],prec=prec)
        residual=local[other](psi)-ratio_part*local[anchor](psi)
        derivative=local[other].derivative()(psi)-ratio_part*local[anchor].derivative()(psi)
        assert derivative[0];psi=(psi-residual/derivative).add_bigoh(prec)
    assert psi[0]==0 and psi[1]
    powers=[PS.one().add_bigoh(depth+1)]
    for r in range(D):powers.append((powers[-1]*psi).add_bigoh(depth+1))
    composed=[sum((f[r]*powers[r] for r in range(D+1) if f[r]),PS.zero()) for f in local]
    scale=(arcs[anchor]/composed[anchor]).add_bigoh(depth+1);assert scale[0]
    for a in range(k):assert all((scale*composed[a])[r]==jets[r][a] for r in range(depth+1)),('incoherent jet',a)
    report=dict(public_only=True,published_secret_key_read=False,curve_derived_from_public_generator_and_extracted_support=True,
        degree_bound=D,actual_curve_degree=max(f.degree() for f in F),common_coordinate_factor_degree=0,
        all_public_columns_matched=True,all_eight_branch_planes_equal_curve_frobenius_planes=True,
        branch_to_frobenius_power=permutation,coherent_scalar_and_parameter_series_verified=True,
        exact_jet_coefficients_checked=k*(depth+1),jet_depth=depth,parameter_constant_zero=True,
        parameter_linear_nonzero=True,scalar_constant_nonzero=True,seconds=time.monotonic()-begin)
    result=dict(report=report,curve=[[enc(a) for a in f.list()] for f in F],ambient_multipliers=[enc(a) for a in multipliers],
        parameter_series=[enc(psi[r]) for r in range(depth+1)],scalar_series=[enc(scale[r]) for r in range(depth+1)],field_modulus=public['field_modulus'])
    (args.result/'coherence-audit.json').write_text(json.dumps(result,indent=2,default=int)+'\n');print(json.dumps(report,default=int),flush=True)

if __name__=='__main__':main()
