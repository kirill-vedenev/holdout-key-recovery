"""Readable SageMath implementation shared by the toy recovery notebooks.

Recovery functions take public matrices and recovered data only. Instance
generation and independent secret-based audits are separate functions below.
"""
from itertools import combinations, combinations_with_replacement, product
from collections import Counter
from functools import lru_cache
from time import perf_counter
import random

from sage.all import (
    GF, PolynomialRing, MatrixSpace, matrix, vector, binomial, prod,
    set_random_seed, random_matrix,
)


def mat(F, *args, **kwargs):
    """Use the generic backend for odd-characteristic extension fields."""
    A = matrix(F, *args, **kwargs)
    if F.characteristic() != 2 and F.degree() > 1:
        A = MatrixSpace(F, A.nrows(), A.ncols(), implementation='generic')(A.list())
    return A


def over(A, F):
    return mat(F, A.nrows(), A.ncols(), A.list())


def kernel(A):
    K = over(A.right_kernel_matrix(), A.base_ring())
    assert (A * K.transpose()).is_zero(), 'kernel residual'
    return K


def same_span(A, B):
    """Exact rank comparison without a backend-changing vector-space constructor."""
    A, B = over(A, A.base_ring()), over(B, A.base_ring())
    return A.rank() == B.rank() == A.stack(B).rank()


def grs_matrix(E, support, dimension, multiplier=None):
    weights = multiplier if multiplier is not None else [E.one()] * len(support)
    return mat(E, dimension, len(support),
               lambda r, j: weights[j] * support[j] ** r)


def generate_grs(q, n, k, D, seed):
    """Return (public generator, audit-only secret)."""
    E = GF(q, 'a')
    rng = random.Random(int(seed))
    support = rng.sample(list(E), n)
    weights = [rng.choice([a for a in E if a]) for _ in range(n)]
    set_random_seed(seed)
    coefficients = random_matrix(E, k, D + 1)
    while coefficients.rank() != k:
        coefficients = random_matrix(E, k, D + 1)
    Y = (coefficients * grs_matrix(E, support, D + 1, weights)).echelon_form()
    return Y, dict(field=E, support=support, multiplier=weights, D=D)


def generate_binary_goppa(m, n, t, seed):
    E = GF(2 ** m, 'a')
    R = PolynomialRing(E, 'Z')
    set_random_seed(seed)
    g = R.irreducible_element(t)
    rng = random.Random(int(seed))
    support = rng.sample([a for a in E if g(a)], n)
    checks = grs_matrix(E, support, 2*t, [1/g(a)**2 for a in support])
    _, _, to_vector = E.vector_space(map=True)
    rows = [[to_vector(checks[r, j])[b] for j in range(n)]
            for r in range(2*t) for b in range(m)]
    Y = kernel(mat(GF(2), rows)).echelon_form()
    derivative = prod(R.gen()-a for a in support).derivative()
    weights = [g(a)**2 / derivative(a) for a in support]
    return Y, dict(field=E, support=support, multiplier=weights,
                   D=n-2*t-1, goppa_polynomial=g, t=t)


def generate_alternant(q, m, n, r, seed):
    E = GF(q**m, 'a')
    F, embedding = E.subfield(E.degree()//m, 'b', map=True)
    if not E.has_coerce_map_from(F):
        E.register_coercion(embedding)
    _, _, to_vector = E.vector_space(embedding, map=True)
    rng = random.Random(int(seed))
    support = rng.sample(list(E),n)
    beta = [rng.choice([a for a in E if a]) for _ in range(n)]
    checks = grs_matrix(E,support,r,beta)
    rows = [[to_vector(checks[a,j])[b] for j in range(n)] for a in range(r) for b in range(m)]
    Y = kernel(mat(F,rows)).echelon_form()
    R = PolynomialRing(E,'Z')
    derivative = prod(R.gen()-a for a in support).derivative()
    weights = [1/(beta[j]*derivative(a)) for j,a in enumerate(support)]
    return Y,dict(field=E,support=support,multiplier=weights,D=n-r-1)


@lru_cache(None)
def monomials(k, degree, cap=None):
    return tuple(I for I in combinations_with_replacement(range(k), degree)
                 if cap is None or max(Counter(I).values(), default=0) <= cap)


def hasse_row(I, point, mons):
    E = point.base_ring()
    deriv = Counter(I)
    result = []
    for J in mons:
        powers = Counter(J)
        if any(powers[a] < b for a, b in deriv.items()):
            result.append(E.zero())
        else:
            result.append(prod((E(binomial(b, deriv[a])) * point[a] ** (b-deriv[a])
                                for a, b in powers.items()), E.one()))
    return result


def holdout_kernel(Y, position, degree, multiplicity, infinity=False):
    """Full literal Hasse constraints; binary assembly exploits squarefreeness."""
    started = perf_counter()
    F = Y.base_ring()
    k, n = Y.dimensions()
    points = [Y.column(j) for j in range(n) if j != position]
    if infinity:
        points.append(sum(Y.columns(), vector(F, k)))
    unit_rows = {next(a for a in range(k) if p[a]) for p in points
                 if sum(bool(x) for x in p) == 1 and next(x for x in p if x) == 1}
    cap = int(F.order())-1
    if len(unit_rows) == k and degree > multiplicity:
        # Exactly impose all retained unit-column conditions by the basis choice.
        cap = min(cap, degree-multiplicity)
        points = [p for p in points if sum(bool(x) for x in p) != 1]
    mons = monomials(k, degree, cap)
    indices = [I for r in range(multiplicity)
               for I in monomials(k, r, int(F.order())-1)]
    if F.order() == 2:
        index = {I: r for r, I in enumerate(indices)}
        A = mat(F, len(points)*len(indices), len(mons))
        for block, p in enumerate(points):
            offset = block*len(indices)
            for c, I in enumerate(mons):
                zeros = tuple(j for j in I if not p[j])
                ones = tuple(j for j in I if p[j])
                for r in range(multiplicity-len(zeros)):
                    for extra in combinations(ones, r):
                        A[offset+index[tuple(sorted(zeros+extra))], c] = 1
    else:
        A = mat(F, [hasse_row(I, p, mons) for p in points for I in indices])
    K = kernel(A)
    values = K * vector(F, hasse_row((), Y.column(position), mons))
    report = dict(matrix_rows=A.nrows(), monomials=len(mons),
                  kernel_dimension=K.nrows(), heldout_values_zero=not bool(values),
                  seconds=perf_counter()-started)
    return K, mons, report


def local_maps(K, mons, point):
    F, k = point.base_ring(), len(point)
    J = K * mat(F, [hasse_row((a,), point, mons) for a in range(k)]).transpose()
    pairs = monomials(k, 2)
    Q = K * mat(F, [hasse_row(I, point, mons) for I in pairs]).transpose()
    return J, Q, pairs


def quadratic_values(v, pairs):
    return vector(v.base_ring(), [v[a]*v[b] for a, b in pairs])


def polar_matrix(Q, v, pairs):
    E, k = v.base_ring(), len(v)
    P = mat(E, len(pairs), k)
    for r, (a, b) in enumerate(pairs):
        P[r, a] += v[b]
        P[r, b] += v[a]
    return over(Q, E)*P


def separate_branches(J, Q, pairs, point, E, extension_degree, binary_goppa=False):
    """Recover projective branches by the quadratic quotient eigenvalue method."""
    F, k = J.base_ring(), J.ncols()
    W = kernel(J)
    assert W.nrows() == extension_degree+1, 'gradient nullity differs from m+1'
    basis = [point]
    for row in W:
        if mat(F, basis+[row]).rank() > len(basis):
            basis.append(row)
    B = mat(F, basis[1:])
    stationary = kernel(J.transpose())
    QQ = Q if binary_goppa else stationary*Q
    b = B.nrows()
    quotient_pairs = monomials(b, 2)
    substitution = mat(F, len(pairs), len(quotient_pairs),
        lambda r, c: B[quotient_pairs[c][0], pairs[r][0]] * B[quotient_pairs[c][1], pairs[r][1]]
        if quotient_pairs[c][0] == quotient_pairs[c][1] else
        B[quotient_pairs[c][0], pairs[r][0]] * B[quotient_pairs[c][1], pairs[r][1]] +
        B[quotient_pairs[c][1], pairs[r][0]] * B[quotient_pairs[c][0], pairs[r][1]])
    restrictions = QQ*substitution
    assert restrictions.rank() == binomial(extension_degree, 2), 'quadratic rank criterion'
    report = dict(gradient_rank=J.rank(), gradient_nullity=W.nrows(),
                  quadratic_rank=restrictions.rank(), quotient_dimension=b,
                  quadratic_family='all forms' if binary_goppa else 'stationary forms')
    if b == 1:
        branches = [vector(E, B[0])]
    else:
        quotient = kernel(restrictions).transpose()
        index = {I: a for a, I in enumerate(quotient_pairs)}
        multiplications = []
        for a in range(b):
            columns = [quotient[index[tuple(sorted((a, j)))] ] for j in range(b)]
            multiplications.append(mat(F, columns).transpose())
        A1 = multiplications[0]
        assert A1.is_invertible(), 'first projective chart misses a branch'
        rng = random.Random(0)
        candidates = multiplications[1:] + [sum(
            (rng.choice(list(F))*A for A in multiplications), mat(F, b, b))
            for _ in range(32)]
        branches = []
        for A2 in candidates:
            characteristic = (A2*A1.inverse()).charpoly().change_ring(E)
            roots = characteristic.roots(multiplicities=False)
            if len(roots) != b:
                continue
            for theta in roots:
                left = kernel((over(A2,E)-theta*over(A1,E)).transpose())
                assert left.nrows() == 1
                coords = left[0]*over(A1,E)
                branches.append(coords*over(B,E))
            break
        assert len(branches) == extension_degree, 'no separating eigenvalue pencil found'
    for v in branches:
        assert not over(J,E)*v
        assert not over(QQ,E)*quadratic_values(v,pairs)
        assert mat(E,[point,v]).rank() == 2
    report['branches'] = len(branches)
    return branches, stationary, report


def series_products(jets, mons, depth):
    E = jets[0].base_ring()
    R = PolynomialRing(E, 'T')
    coordinates = [R([v[a] for v in jets]) for a in range(len(jets[0]))]
    cache = {(): R.one()}
    def get(I):
        if I not in cache:
            cache[I] = (get(I[:-1])*coordinates[I[-1]]).truncate(depth+1)
        return cache[I]
    polys = [get(I) for I in mons]
    return mat(E, len(mons), depth+1, lambda r,c: polys[r][c])


def fixed_solver(A):
    columns = A.pivots()
    independent = A.matrix_from_columns(columns)
    rows = independent.transpose().pivots()
    inverse = independent.matrix_from_rows(rows).inverse()
    def solve(rhs):
        v = vector(A.base_ring(), A.ncols())
        values = inverse*vector(A.base_ring(), [rhs[r] for r in rows])
        for c, value in zip(columns,values):
            v[c] = value
        assert A*v == rhs, 'inconsistent continuation equations'
        return v
    return solve


def continue_jet(K, mons, J, Q, pairs, stationary, point, v1, depth, binary_goppa=False):
    started = perf_counter()
    E, k = v1.base_ring(), len(v1)
    JE, KE, SE = over(J,E), over(K,E), over(stationary,E)
    L = polar_matrix(Q,v1,pairs)
    A = JE.stack(L if binary_goppa else SE*L)
    assert A.rank() == k-2, 'continuation rank condition'
    assert (A*mat(E,[point,v1]).transpose()).is_zero()
    solve = fixed_solver(A)
    jets = [vector(E,point),v1]
    if binary_goppa:
        jets.append(vector(E,[a*a for a in v1]))
        for r in range(3, depth+1, 2):
            even = vector(E,[a*a for a in jets[(r+1)//2]])
            residual = KE*series_products(jets,mons,r+1)
            rhs = -vector(E,list(residual.column(r))+list(residual.column(r+1)+JE*even))
            jets.extend([solve(rhs),even])
        assert all(jets[2*r] == vector(E,[a*a for a in jets[r]])
                   for r in range(1,len(jets)//2))
    else:
        for r in range(2,depth+1):
            residual = KE*series_products(jets,mons,r+1)
            rhs = -vector(E,list(residual.column(r))+list(SE*residual.column(r+1)))
            jets.append(solve(rhs))
    # Fresh full substitution checks the final series, independently of each solve.
    assert (KE*series_products(jets,mons,len(jets)-1)).is_zero()
    first_full = next((r for r in range(len(jets)) if mat(E,jets[:r+1]).rank()==k), None)
    return jets, dict(continuation_rank=A.rank(), ambiguity_dimension=k-A.rank(),
                     requested_depth=depth, recovered_depth=len(jets)-1,
                     first_full_span_order=first_full, full_substitution=True,
                     binary_normalization=binary_goppa, seconds=perf_counter()-started)


def evaluate_forms(H, mons, Y, E):
    values = mat(E,len(mons),Y.ncols(),
                 lambda r,j: prod((E(Y[a,j]) for a in mons[r]),E.one()))
    return over(H,E)*values


def bootstrap(Y, E, jets, D, form_degree=2):
    depth = max(0,form_degree*D-Y.ncols()+1)
    assert len(jets)>depth
    mons = monomials(Y.nrows(),form_degree)
    values = mat(E,len(mons),Y.ncols(),
                 lambda r,j: prod((E(Y[a,j]) for a in mons[r]),E.one()))
    contact = series_products(jets[:depth+1],mons,depth)
    A = values.augment(contact.matrix_from_columns(range(1,depth+1))).transpose()
    H = kernel(A)
    tangents, ranks = [], []
    for j in range(Y.ncols()):
        p = vector(E,Y.column(j))
        J = H*mat(E,[hasse_row((a,),p,mons) for a in range(Y.nrows())]).transpose()
        W = kernel(J)
        assert W.nrows()==2, 'bootstrap tangent rank failure at column '+str(j)
        tangents.append(next(v for v in W if mat(E,[p,v]).rank()==2))
        ranks.append(J.rank())
    return dict(forms=H,monomials=mons,tangents=mat(E,tangents).transpose(),
                report=dict(degree=form_degree,jet_order=depth,constraint_rank=A.rank(),
                            ideal_dimension=H.nrows(),gradient_ranks=sorted(set(ranks)),
                            positions_verified=Y.ncols()))


def public_multiplier(Y, E, support, D):
    V = grs_matrix(E,support,D+1)
    H = kernel(V)
    YE = over(Y,E)
    A = mat(E,[g.pairwise_product(h) for g in YE for h in H])
    B = kernel(A)
    assert B.nrows(), 'no compatible multiplier'
    candidates = list(B.rows())
    rng = random.Random(0)
    candidates += [vector(E,[rng.choice(list(E)) for _ in range(B.nrows())])*B
                   for _ in range(200)]
    z = next((v for v in candidates if all(v)),None)
    assert z is not None, 'nonvanishing multiplier search exhausted'
    weights = [1/a for a in z]
    recovered = grs_matrix(E,support,D+1,weights)
    assert (YE*kernel(recovered).transpose()).is_zero()
    return weights, dict(multiplier_nullity=B.nrows(),public_grs_containment=True)


def finite_support(E, projective):
    if None not in projective:
        return projective
    pole = next(a for a in E if a not in projective)
    return [E.zero() if a is None else 1/(a-pole) for a in projective]


def direct_finish(Y,E,jets,position,D,form_degree=2,binary_goppa=False):
    L = form_degree*D
    depth = L-1
    assert len(jets)>depth
    mons = monomials(Y.nrows(),form_degree)
    contact = series_products(jets[:depth+1],mons,depth).transpose()
    WR = kernel(contact)
    RV = evaluate_forms(WR,mons,Y,E)
    r_index = next((a for a in range(WR.nrows()) if RV[a]),None)
    assert r_index is not None, 'required numerator section is absent'
    R = WR[r_index]
    target = L-2 if binary_goppa else L-1
    WS = kernel(contact.matrix_from_rows(range(target)))
    leading = WS*contact.row(target)
    s_index = next((a for a,c in enumerate(leading) if c),None)
    assert s_index is not None, 'required denominator section is absent'
    S = WS[s_index]/leading[s_index]
    values = evaluate_forms(mat(E,[R,S]),mons,Y,E)
    support = []
    for j in range(Y.ncols()):
        if j==position:
            support.append(E.zero())
        elif not values[1,j]:
            assert values[0,j]
            support.append(None)
        else:
            value=values[0,j]/values[1,j]
            support.append(value**(E.order()//2) if binary_goppa else value)
    assert len(set(support))==Y.ncols()
    support=finite_support(E,support)
    weights,report=public_multiplier(Y,E,support,D)
    report.update(form_degree=form_degree,jet_order=depth,numerator_space=WR.nrows(),
                  denominator_space=WS.nrows(),square_root=binary_goppa,
                  recovered_positions=len(support),distinct_support=True)
    return dict(support=support,multiplier=weights,sections=mat(E,[R,S]),
                monomials=mons,report=report)


def sidelnikov_shestakov(T):
    E,n=T.base_ring(),T.ncols()
    G=T.row_space().basis_matrix()
    if n-G.nrows()<G.nrows():
        G=kernel(G)
    G=G.echelon_form()
    k=G.nrows()
    assert 2<=k<=n-2
    pivots=list(G.pivots())
    order=pivots+[j for j in range(n) if j not in pivots]
    M=G.matrix_from_columns(order)
    assert all(M[a,j] for a in range(k) for j in range(k,n))
    for b in E:
        if b in (0,1):
            continue
        alpha=[None]*n
        alpha[0],alpha[1],alpha[k]=E.zero(),E.one(),b
        try:
            for j in range(k+1,n):
                ratio=M[0,j]*M[1,k]/(M[1,j]*M[0,k])
                alpha[j]=b/(b-ratio*(b-1))
            c=alpha[k+1]
            for a in range(2,k):
                ratio=M[a,k]*M[0,k+1]/(M[0,k]*M[a,k+1])
                alpha[a]=b*c*(ratio-1)/(ratio*c-b)
        except ZeroDivisionError:
            continue
        if len(set(alpha))!=n:
            continue
        support=[E.zero()]*n
        for a,j in enumerate(order):
            support[j]=alpha[a]
        try:
            weights,_=public_multiplier(G,E,support,k-1)
        except AssertionError:
            continue
        assert same_span(grs_matrix(E,support,k,weights),G)
        return support
    raise AssertionError('no finite GRS representation found')


def minor_finish(Y,E,bootstrap_data,D,expected_dimension):
    assert E.characteristic()==2
    tangent=bootstrap_data['tangents']
    T=mat(E,[[ (E(Y[a,j])*tangent[b,j]-E(Y[b,j])*tangent[a,j])**(E.order()//2)
              for j in range(Y.ncols())] for a,b in combinations(range(Y.nrows()),2)])
    assert T.rank()==expected_dimension, 'minor fullness condition'
    support=sidelnikov_shestakov(T)
    minor_weights,_=public_multiplier(T,E,support,expected_dimension-1)
    assert same_span(grs_matrix(E,support,expected_dimension,minor_weights),T)
    weights,report=public_multiplier(Y,E,support,D)
    report.update(minor_rows=T.nrows(),minor_rank=T.rank(),exact_grs_equality=True,
                  recovered_positions=len(support))
    return dict(support=support,multiplier=weights,minor=T,report=report)


def equivalent_support(original,recovered,E,q):
    """Independent audit: one Frobenius power followed by one Mobius map."""
    m=0
    while q**m<E.cardinality():
        m+=1
    assert q**m==E.cardinality()
    for sigma in range(m):
        alpha=[a**(q**sigma) for a in original]
        rows=[[alpha[j],1,-recovered[j]*alpha[j],-recovered[j]] for j in range(3)]
        B=kernel(mat(E,rows))
        if B.nrows()!=1:
            continue
        a,b,c,d=B[0]
        if a*d-b*c and all(c*x+d and (a*x+b)/(c*x+d)==y for x,y in zip(alpha,recovered)):
            return dict(frobenius_power=sigma,mobius_equivalence=True)
    raise AssertionError('recovered support is not globally equivalent')


def audit_key(Y,result,D):
    """Verify the recovered supercode and exact subfield equality when applicable."""
    F=Y.base_ring()
    E=result['support'][0].parent()
    recovered=grs_matrix(E,result['support'],D+1,result['multiplier'])
    H=kernel(recovered)
    assert (over(Y,E)*H.transpose()).is_zero()
    report=dict(distinct_support=len(set(result['support']))==Y.ncols(),
                public_grs_containment=True)
    if F.order()!=E.order():
        embedding=E.coerce_map_from(F)
        _,_,to_vector=E.vector_space(embedding,map=True)
        m=E.degree()//F.degree()
        rows=[[to_vector(h[j])[b] for j in range(Y.ncols())]
              for h in H for b in range(m)]
        recovered_public=kernel(mat(F,rows))
        assert same_span(recovered_public,Y)
        report['exact_public_subfield_code_equality']=True
    return report


def audit_geometry(Y, secret, K, mons, position, jets, q, bootstrap_data=None, branches=None):
    """Secret-based audit, called only after public recovery has finished.

    Fit one scalar series and one parameter series to every recovered coefficient.
    This checks more than equality of successive derivative spans.
    """
    E,D=secret['field'],secret['D']
    R=PolynomialRing(E,'Z')
    F=[R.lagrange_polynomial([(a,E(Y[r,j])/secret['multiplier'][j])
                             for j,a in enumerate(secret['support'])])
       for r in range(Y.nrows())]
    assert max(f.degree() for f in F)<=D
    polynomial_cache={():R.one()}
    def product_polynomial(I):
        if I not in polynomial_cache:
            polynomial_cache[I]=product_polynomial(I[:-1])*F[I[-1]]
        return polynomial_cache[I]
    pulled=[product_polynomial(I) for I in mons]
    coefficient_matrix=mat(E,len(mons),max(f.degree()+1 for f in pulled),
                           lambda r,c:pulled[r][c])
    assert (over(K,E)*coefficient_matrix).is_zero()
    seed=secret['support'][position]
    scale=secret['multiplier'][position]
    P=[scale*f(R.gen()+seed) for f in F]
    extension_degree=0
    while q**extension_degree<E.order():
        extension_degree+=1
    if branches is not None:
        expected=[mat(E,[jets[0],vector(E,[f[1]**(q**s) for f in P])])
                  for s in range(extension_degree)]
        found=[mat(E,[jets[0],v]) for v in branches]
        assert len(found)==len(expected) and all(any(same_span(V,W) for W in expected) for V in found)
        assert all(any(same_span(V,W) for W in found) for V in expected)
    sigma=next(s for s in range(extension_degree)
               if mat(E,[jets[0],jets[1],vector(E,[f[1]**(q**s) for f in P])]).rank()==2)
    P=[R([c**(q**sigma) for c in f]) for f in P]
    basis=mat(E,[jets[0],vector(E,[f[1] for f in P])]).transpose()
    solve=fixed_solver(basis)
    scalar,parameter=R.one(),R.zero()
    depth=len(jets)-1
    def substitute_truncated(f,parameter,precision):
        out=R.zero()
        for c in reversed(f.list()):
            out=(out*parameter+c).truncate(precision)
        return out
    for r in range(1,depth+1):
        predicted=vector(E,[(scalar*substitute_truncated(f,parameter,r+1))[r] for f in P])
        a,b=solve(jets[r]-predicted)
        scalar+=a*R.gen()**r
        parameter+=b*R.gen()**r
    assert scalar[0]==1 and parameter[0]==0 and parameter[1]
    actual=[(scalar*substitute_truncated(f,parameter,depth+1)).truncate(depth+1) for f in P]
    assert all(actual[a][r]==jets[r][a] for a in range(Y.nrows()) for r in range(depth+1))
    report=dict(exact_holdout_curve_identities=True,audited_polynomials=K.nrows(),coherent_jet=True,
                coherent_coordinate_coefficients=Y.nrows()*(depth+1),
                audited_frobenius_branch=sigma)
    if branches is not None:
        report['all_branches_match_secret_audit']=True
    if bootstrap_data is not None:
        selected_curve=[R([c**(q**sigma) for c in f]) for f in F]
        bm=bootstrap_data['monomials']
        pulled=[prod((selected_curve[a] for a in I),R.one()) for I in bm]
        coefficient_matrix=mat(E,len(bm),max(f.degree()+1 for f in pulled),
                               lambda r,c:pulled[r][c])
        true_ideal=kernel(coefficient_matrix.transpose())
        assert same_span(bootstrap_data['forms'],true_ideal)
        for j,alpha in enumerate(secret['support']):
            alpha=alpha**(q**sigma)
            expected=mat(E,[Y.column(j),vector(E,[f.derivative()(alpha) for f in selected_curve])])
            found=mat(E,[Y.column(j),bootstrap_data['tangents'].column(j)])
            assert same_span(expected,found)
        report.update(exact_full_branch_ideal=True,all_tangents_on_same_branch=True)
    return report


def reject_tampered_key(Y,result,D):
    """Negative control: corrupt one multiplier and require the public check to fail."""
    E=result['support'][0].parent()
    altered=list(result['multiplier'])
    altered[0]*=next(a for a in E if a not in (0,1))
    bad_grs=grs_matrix(E,result['support'],D+1,altered)
    rejected=not (over(Y,E)*kernel(bad_grs).transpose()).is_zero()
    assert rejected, 'the negative control was unexpectedly accepted'
    return True
