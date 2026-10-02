"""Public-only IKZ-type extension of the recovered TII-249 support.

Factor the common retained-coordinate multiplier block once. Each restored
coordinate adds one generator row, so its candidate test reduces exactly to
32 moments. Every passing candidate is then checked against the entire public
one-coordinate extension. No known support or challenge answer is an input.
"""
from sage.all import *
from pathlib import Path
import argparse,hashlib,json,time

PUBLIC_SHA256='63ef408ba6f954272e43a0563dff373064a70c24601a84611a00d977d558b598'


def enc(a):return int(a.to_integer())
def dump(path,value):path.write_text(json.dumps(value,indent=2,default=int)+'\n')


def binary_expansion(A):
    m=A.base_ring().degree()
    return matrix(GF(2),A.nrows()*m,A.ncols(),lambda r,j:(enc(A[r//m,j])>>(r%m))&1)


def multiplier_matrix(G,alpha,r):
    E=G.base_ring();V=matrix(E,r,len(alpha),lambda a,j:alpha[j]**a)
    return matrix(E,G.nrows()*r,len(alpha),lambda a,j:G[a//r,j]*V[a%r,j])


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--public',type=Path,required=True);ap.add_argument('--retained-support',type=Path,required=True)
    ap.add_argument('--original-public-key',type=Path,required=True);ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args();begin=time.monotonic();out=args.output;out.mkdir(parents=True,exist_ok=True);(out/'positions').mkdir(exist_ok=True)
    if (out/'equivalent-key.json').exists():ap.error('completed output already exists; choose a fresh folder')
    public=json.loads(args.public.read_text());direct=json.loads(args.retained_support.read_text());params=public['original_parameters']
    n,m,t,k_full=[params[x] for x in ['n','m','t','k']];r=2*t;keep=public['original_indices'];removed=public['shortened_original_indices']
    assert n==235 and m==8 and t==16 and k_full==107 and direct['original_indices']==keep
    raw=args.original_public_key.read_bytes();assert hashlib.sha256(raw).hexdigest()==PUBLIC_SHA256==public['public_sha256']
    lines=[line.strip() for line in raw.decode().splitlines() if line.strip()]
    H=matrix(GF(2),[[int(x) for x in line[1:-1].split()] for line in lines[:-1]])
    assert H.ncols()==n and H.rank()==128 and json.loads(lines[-1])==public['field_modulus']
    E=GF(2**m,'a',modulus=PolynomialRing(GF(2),'x')(public['field_modulus']));dec=E.from_integer
    R=PolynomialRing(E,'X');X=R.gen();P=matrix(GF(2),public['generator']);k=P.nrows();Gfull=H.right_kernel_matrix()
    assert k==50 and Gfull.nrows()==107 and H.matrix_from_columns(keep).right_kernel()==P.row_space()
    assert len(keep)==178 and len(removed)==57 and sorted(keep+removed)==list(range(n))
    beta=[None if x is None else dec(x) for x in direct['support']]
    forms=matrix(E,[[dec(x) for x in row] for row in direct['forms']]);pairs=direct['quadratic_pairs']
    pinfinity=sum(P.change_ring(E).columns(),vector(E,k))
    assert pinfinity==vector(E,[(public['infinity_mask']>>a)&1 for a in range(k)])
    high,lower=forms*vector(E,[pinfinity[a]*pinfinity[b] for a,b in pairs])
    assert high,'public infinity section degenerates'
    beta_infinity=(high/lower)**(2**(m-1)) if lower else None
    assert beta_infinity not in beta
    def affine(z):
        if beta_infinity is None:assert z is not None;return z
        return E.zero() if z is None else 1/(z+beta_infinity)
    alpha=[affine(z) for z in beta];assert len(set(alpha))==len(alpha)
    chart=dict(public_only=True,public_point_at_infinity_mask=public['infinity_mask'],
        quadratic_section_values=[enc(high),enc(lower)],recovered_beta_at_public_infinity=None if beta_infinity is None else enc(beta_infinity),
        transformation='identity' if beta_infinity is None else 'alpha=1/(beta+beta_infinity)',
        original_infinity_maps_to_projective_infinity=True,all_retained_support_finite=True,
        retained_original_indices=keep,retained_support=[enc(a) for a in alpha])
    dump(out/'affine-chart.json',chart)

    # Every extension contains the same shortened code with zero at the new
    # coordinate. Its equations force the old multipliers into this 1D space.
    fixed=multiplier_matrix(P.change_ring(E),alpha,r);fixed_basis=fixed.right_kernel_matrix()
    assert fixed_basis.nrows()==1 and all(fixed_basis.row(0))
    nu_retained=fixed_basis.row(0)/fixed_basis[0,0];assert not fixed*nu_retained
    fixed_rank=fixed.rank();assert fixed_rank==len(keep)-1
    powers=matrix(E,r,len(keep),lambda a,j:alpha[j]**a)
    finite_candidates=[dec(a) for a in range(E.order()) if dec(a) not in set(alpha)]
    projective_candidates=finite_candidates+[None]
    support=[None]*n
    for j,a in zip(keep,alpha):support[j]=a
    histories=[];extensions_begin=time.monotonic()
    print(json.dumps(dict(event='affine_chart_and_fixed_block_ready',beta_infinity=None if beta_infinity is None else enc(beta_infinity),
        fixed_shape=[fixed.nrows(),fixed.ncols()],fixed_rank=fixed_rank,fixed_nullity=1,candidates=len(projective_candidates),seconds=time.monotonic()-begin)),flush=True)
    for j in removed:
        positions=keep+[j];Gj=H.matrix_from_columns(positions).right_kernel_matrix();assert Gj.nrows()==k+1
        embedded=P.augment(matrix(GF(2),k,1));extra=next(row for row in Gj.rows() if row[-1])
        assert embedded.stack(matrix(GF(2),[extra])).row_space()==Gj.row_space()
        extra=extra.change_ring(E)/E(extra[-1]);moments=powers*vector(E,[extra[a]*nu_retained[a] for a in range(len(keep))])
        assert moments,'the extension failed to remove the retained multiplier freedom'
        pivot=next(a for a,x in enumerate(moments) if x);accepted=[];tests=[]
        for w in projective_candidates:
            target=vector(E,[int(a==r-1) if w is None else w**a for a in range(r)])
            scalar=target[pivot]/moments[pivot];residual=scalar*moments-target;consistent=not residual
            full_support=consistent and bool(scalar)
            tests.append(dict(value=None if w is None else enc(w),projective_infinity=w is None,
                kernel_dimension=int(consistent),has_full_support_multiplier=bool(full_support),
                retained_multiplier_scalar=enc(scalar),moment_residual=[enc(x) for x in residual],
                nonzero_moment_equations=sum(bool(x) for x in residual)))
            if full_support:accepted.append(w)
        assert len(accepted)==1,(j,'non-unique extension',accepted)
        assert accepted[0] is not None,(j,'public affine support contains infinity')
        support[j]=accepted[0];winner=next(test for test in tests if test['has_full_support_multiplier'])
        candidate_nu=vector(E,list(dec(winner['retained_multiplier_scalar'])*nu_retained)+[E.one()])
        full_equations=multiplier_matrix(Gj.change_ring(E),alpha+[accepted[0]],r)
        assert all(candidate_nu) and not full_equations*candidate_nu
        history=dict(original_position=j,extension_positions=positions,extension_generator_dimension=Gj.nrows(),
            common_known_block_rank=fixed_rank,common_known_block_nullity=1,known_multiplier_matrix_rank=fixed_rank+1,
            known_block_rank_certificate='the added row has nonzero image on the one-dimensional retained nullspace',
            fixed_block_factorizations=0,common_fixed_block_factorizations_for_entire_run=1,
            moment_pivot=pivot,moments=[enc(x) for x in moments],finite_candidates_tested=len(finite_candidates),
            projective_infinity_tested=True,total_candidates_tested=len(tests),accepted_value=enc(accepted[0]),unique_candidate=True,
            all_winner_public_equations_verified=True,winner_full_multiplier=[enc(x) for x in candidate_nu],tests=tests)
        dump(out/'positions'/f'position-{j}.json',history)
        histories.append(dict(original_position=j,support_value=enc(accepted[0]),unique_candidate=True,
            candidate_count=len(tests),infinity_rejected=True,evidence=f'positions/position-{j}.json'))
        print(json.dumps(dict(event='deshortened_position',position=j,support_value=enc(accepted[0]),candidates=len(tests),seconds=time.monotonic()-begin)),flush=True)
    assert all(a is not None for a in support) and len(set(support))==n
    extension_seconds=time.monotonic()-extensions_begin

    full_begin=time.monotonic();A=multiplier_matrix(Gfull.change_ring(E),support,r);space=A.right_kernel_matrix()
    assert space.nrows()==1 and all(space.row(0));nu=space.row(0)/space[0,0];assert not A*nu
    Pi=prod(X-a for a in support);Pi_prime=Pi.derivative();lam=vector(E,[1/(nu[j]*Pi_prime(support[j])) for j in range(n)])
    parity=matrix(E,r,n,lambda a,j:nu[j]*support[j]**a)
    ambient=matrix(E,n-r,n,lambda a,j:lam[j]*support[j]**a)
    assert not Gfull.change_ring(E)*parity.transpose() and not ambient*parity.transpose()
    assert parity.rank()==r and ambient.rank()==n-r and binary_expansion(parity).row_space()==H.row_space()
    multiplier_report=dict(full_multiplier_system_shape=[A.nrows(),A.ncols()],rank=A.rank(),nullity=space.nrows(),
        all_multipliers_nonzero=True,normalization='dual_grs_multipliers[0]=1',
        dual_convention='H32[a,j]=nu[j]*alpha[j]^a for a=0..31',
        primal_convention='C203[a,j]=lambda[j]*alpha[j]^a for a=0..202',
        dual_primal_relation='lambda[j]=1/(nu[j]*prod_{l!=j}(alpha[j]+alpha[l]))',
        unavoidable_common_scalar=True,ambient_grs_dimension=ambient.rank(),binary_parity_rank=binary_expansion(parity).rank(),
        exact_original_public_row_space_equality=True,seconds=time.monotonic()-full_begin)
    dump(out/'full-multipliers.json',dict(report=multiplier_report,support=[enc(x) for x in support],
        dual_grs_multipliers=[enc(x) for x in nu],ambient_grs_multipliers=[enc(x) for x in lam],
        field_modulus=public['field_modulus'],original_indices=list(range(n))))

    goppa_begin=time.monotonic();square=R.lagrange_polynomial([(a,1/b) for a,b in zip(support,nu)])
    square_ok=square.degree()==r and not square.derivative();g=None;goppa_report=dict(reciprocal_multiplier_polynomial_degree=square.degree(),
        reciprocal_multiplier_polynomial_is_square=bool(square_ok))
    if square_ok:
        g=R([square[2*a]**(2**(m-1)) for a in range(t+1)]).monic();assert g.degree()==t
        scalar=square.leading_coefficient();assert square==scalar*g*g and all(g(a) for a in support)
        gp=matrix(E,t,n,lambda a,j:support[j]**a/g(support[j]));gp_binary=binary_expansion(gp)
        assert gp_binary.row_space()==H.row_space() and not Gfull.change_ring(E)*gp.transpose()
        goppa_report.update(goppa_degree=g.degree(),goppa_monic=True,goppa_squarefree=gcd(g,g.derivative()).degree()==0,
            goppa_irreducible=g.is_irreducible(),goppa_nonzero_on_full_support=True,goppa_binary_parity_rank=gp_binary.rank(),
            goppa_public_row_space_exactly_equal=True,reciprocal_multiplier_scale=enc(scalar),
            reciprocal_multiplier_identity='1/nu[j]=reciprocal_multiplier_scale*g(alpha[j])^2',
            seconds=time.monotonic()-goppa_begin)
        assert goppa_report['goppa_irreducible'],'equivalent Goppa polynomial is reducible'
    else:goppa_report.update(status='specific_multiplier_obstruction',obstruction='reciprocal multipliers do not interpolate to a degree-32 characteristic-two square')
    key=dict(challenge='TII-249',parameters=params,public_only_recovery=True,field_modulus=public['field_modulus'],
        original_indices=list(range(n)),support=[enc(x) for x in support],dual_grs_multipliers=[enc(x) for x in nu],
        ambient_grs_multipliers=[enc(x) for x in lam],multiplier_conventions=multiplier_report,
        public_key_sha256=PUBLIC_SHA256,affine_chart=chart,reciprocal_multiplier_polynomial=[enc(x) for x in square.list()])
    if g is not None:key.update(g=[enc(x) for x in g.list()],goppa_parity_multipliers=[enc(1/g(a)) for a in support],
        reciprocal_multiplier_scale=enc(square.leading_coefficient()))
    dump(out/'equivalent-key.json',key)
    report=dict(status='success',public_only=True,local_only=True,published_secret_key_read=False,
        original_parameters=params,retained_support_count=len(keep),removed_positions=removed,restored_support_count=len(removed),
        full_support_count=n,full_support_distinct=True,full_support_finite=True,
        support_extension='IKZ-type exhaustive projective Vandermonde multiplier test',common_block_factorizations=1,
        per_position_histories=histories,total_candidates_tested=len(removed)*len(projective_candidates),
        all_57_candidates_unique=True,all_projective_infinity_candidates_rejected=True,
        original_public_key_sha256=PUBLIC_SHA256,input_retained_support_sha256=hashlib.sha256(args.retained_support.read_bytes()).hexdigest(),
        multiplier_report=multiplier_report,goppa_report=goppa_report,
        phase_seconds=dict(support_extensions=extension_seconds,full_multipliers=multiplier_report['seconds'],goppa=goppa_report.get('seconds',0)),
        seconds=time.monotonic()-begin)
    dump(out/'recovery-report.json',report)
    print(json.dumps({a:b for a,b in report.items() if a!='per_position_histories'},default=int),flush=True)

if __name__=='__main__':main()
