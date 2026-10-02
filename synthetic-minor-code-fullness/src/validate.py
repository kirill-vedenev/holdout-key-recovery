"""Independent checks of the backends on small instances.

Run with: sage -python src/validate.py
"""
import tempfile
from copy import copy
from math import comb
import numpy as np
from sage.all import *
import gf2m
import native
from engine import run_case, seed_for
from families import CLASSIC, MODULI, SHORTENED_DIMENSIONS, all_parameters, build_instance, goppa_matrix
from minor import (curve_matches_closed_form, evaluate_many, grs_generator,
                   grs_multiplier_fast, grs_multiplier_systematic, grs_parity, hidden_curve,
                   minor_code, support_derivatives, tangent_defects, tangent_matrix)

PARAMS = {P['label']: P for P in all_parameters()}
set_random_seed(219)

# Field arithmetic and GRS helpers against Sage.
for degree in (2,4,6,8,12,13):
    E=GF(2**degree,'v');R=PolynomialRing(E,'Z');Z=R.gen()
    points=list(E)[:min(37,E.order())]
    polynomial=R.random_element(19)
    assert evaluate_many(polynomial,points)==[polynomial(x) for x in points]
    derivative=prod(Z-x for x in points).derivative()
    assert list(support_derivatives(E,points))==[derivative(x) for x in points]
print('PASS table evaluation and support derivatives against Sage polynomials')
for order in (16,64,256):
    E=GF(order,'v'); points=list(E)[:12]
    mult=[E.random_element() or E.one() for _ in points]
    for d in (3,7,10):
        H=grs_parity(E,points,d,mult)
        G=grs_generator(E,points,d,mult)
        assert H.rank()==len(points)-d and (G*H.transpose()).is_zero()
        assert H.row_space()==codes.GeneralizedReedSolomonCode(points,d,mult).parity_check_matrix().row_space()
        change=random_matrix(E,d,d)
        while not change.is_invertible():
            change=random_matrix(E,d,d)
        recovered=grs_multiplier_systematic(E,change*G,points)
        assert recovered is not None
        assert (G*grs_parity(E,points,d,recovered).transpose()).is_zero()
        damaged=copy(G.echelon_form())
        damaged[0,-1]+=1
        assert grs_multiplier_systematic(E,damaged,points) is None
print('PASS dual-GRS parity formula, systematic multiplier recovery, rejection of corrupted GRS matrices')

# NumPy field arithmetic against Sage.
for m,modulus in MODULI.items():
    poly=PolynomialRing(GF(2),'z')([(modulus>>i)&1 for i in range(m+1)])
    assert poly.is_irreducible() and poly.degree()==m
    E=GF(2**m,'z',modulus=poly)
    a=[E.random_element() for _ in range(50)]; b=[E.random_element() or E.one() for _ in range(50)]
    ai=np.array([x.to_integer() for x in a]); bi=np.array([x.to_integer() for x in b])
    assert list(gf2m.multiply(modulus,ai,bi))==[(x*y).to_integer() for x,y in zip(a,b)]
    assert list(gf2m.inverse(modulus,bi))==[(1/y).to_integer() for y in b]
    points=list(dict.fromkeys(a))
    R=PolynomialRing(E,'Z');Z=R.gen();f=R.random_element(9)
    pi=np.array([x.to_integer() for x in points])
    assert list(gf2m.evaluate(modulus,[c.to_integer() for c in f],pi))==[f(x).to_integer() for x in points]
    derivative=prod(Z-x for x in points).derivative()
    assert list(gf2m.support_derivatives(modulus,pi))==[derivative(x).to_integer() for x in points]
print('PASS field polynomials irreducible; NumPy GF(2^m) arithmetic against Sage')


def sage_field(m):
    poly=PolynomialRing(GF(2),'z')([(MODULI[m]>>i)&1 for i in range(m+1)])
    return GF(2**m,'z',modulus=poly)


# C++ binary Goppa trials against an independent Sage reconstruction.
SPLITS={5:[2,3],12:[5,7]}
for label in ('binary-g1','binary-g3','binary-g4','binary-g5','binary-g7','binary-g9','binary-s1','binary-s2'):
    P=PARAMS[label]
    m,t,n0,ell=P['m'],P['t'],P['n0'],P.get('ell',0)
    seeds=[seed_for('validation',label,0,phase)&((1<<64)-1) for phase in ('instance','representatives','minor_sample','recovery')]
    r=native.goppa_trial(MODULI[m],n0,t,ell,P['G'],seeds,True,True)
    again=native.goppa_trial(MODULI[m],n0,t,ell,P['G'],seeds,True,True)
    assert np.array_equal(r['Y'],again['Y']) and np.array_equal(r['Yhat'],again['Yhat']) and r['minor_rank']==again['minor_rank']
    E=sage_field(m);R=PolynomialRing(E,'Z');Z=R.gen()
    G=R([E.from_integer(int(c)) for c in r['goppa_polynomial']])
    assert G.degree()==t and G.is_monic() and G.is_squarefree()
    factors=sorted(f.degree() for f,e in G.factor())
    assert factors==({'irred':[t],'split':[1]*t,'separable':SPLITS.get(t)}[P['G']])
    alpha0=[E.from_integer(int(x)) for x in r['alpha0']]
    assert len(set(alpha0))==n0 and all(G(a)!=0 for a in alpha0)
    keep=list(r['keep']); dropped=[j for j in range(n0) if j not in set(keep)]
    full=goppa_matrix(E,alpha0,G)
    k0=full.nrows()
    if dropped:
        assert full.matrix_from_columns(dropped).rank()==ell     # information positions
        shortened=(full.matrix_from_columns(dropped).left_kernel_matrix()*full).matrix_from_columns(keep)
    else:
        shortened=full
    Yb=matrix(GF(2),r['Y'].tolist())
    assert Yb==shortened.echelon_form() and Yb.nrows()==k0-ell==r['k'] and k0==r['k0']
    alpha=[alpha0[j] for j in keep]; n=len(alpha); D=n-2*t-1
    Y=Yb.change_ring(E)
    Yhat=matrix(E,[[E.from_integer(int(x)) for x in row] for row in r['Yhat']])
    dPi=prod(Z-a for a in alpha).derivative()
    lam=[G(a)**2/dPi(a) for a in alpha]
    curve=hidden_curve(E,Y,alpha,lam,D)
    assert curve_matches_closed_form(E,Y,Yhat,curve,alpha), label
    assert tangent_defects(Y,Yhat)==list(r['defects'])
    T=minor_code(Y,Yhat)
    rank=T.rank()
    assert rank==r['minor_rank']==D-t, (label,rank,r['minor_rank'])
    assert grs_multiplier_fast(E,T,alpha,rank) is not None
    ar=[E.from_integer(int(x)) for x in r['recovered_support']]
    lm=[E.from_integer(int(x)) for x in r['recovered_multiplier']]
    assert (Y*grs_parity(E,ar,D+1,lm).transpose()).is_zero()
    print('PASS C++ binary Goppa trial against Sage:',label,'[%d,%d] -> [%d,%d], minor rank %d over all %d pairs'
          %(n0,k0,n,r['k'],rank,comb(r['k'],2)))

# The NumPy embedding check rejects a wrong multiplier.
modulus=MODULI[P['m']]
Y=r['Y'];alpha=r['recovered_support'].astype(np.int64);lm=r['recovered_multiplier'].astype(np.int64)
assert gf2m.binary_in_grs(modulus,Y,alpha,r['D']+1,lm)
lm[0]=gf2m.multiply(modulus,lm[0],3)
assert not gf2m.binary_in_grs(modulus,Y,alpha,r['D']+1,lm)
print('PASS NumPy embedding check accepts the recovered multiplier and rejects a changed one')

# Tangent formula against interpolation for the Sage families.
for label in ('grs-B','alternant-a8','wild-w2','wild-w6'):
    inst=build_instance(PARAMS[label],230)
    E,Y,alpha,lam,D=inst['E'],inst['Y'],inst['alpha'],inst['lam'],inst['D']
    curve=hidden_curve(E,Y,alpha,lam,D)
    assert curve_matches_closed_form(E,Y,tangent_matrix(E,Y,alpha,inst['nu']),curve,alpha), label
print('PASS tangent formula against interpolated curves (GRS, alternant, wild Goppa)')

# C++ minor backend against the Sage backend with explicit all-pairs certificates.
with tempfile.TemporaryDirectory() as tmp:
    for label in ('grs-A','grs-B','alternant-a1','alternant-a8','wild-w2','wild-w6','wild-w8'):
        P=PARAMS[label]
        a=run_case(P,0,'validation',tmp,False)
        b=run_case(P,0,'validation',tmp,True,backend='sage')
        for key in ('status','fixture_sha256','minor_rank','full_grs','embedding_verified'):
            assert a.get(key)==b.get(key),(label,key,a,b)
        assert a['status']=='pass',(label,a)
        print('PASS C++ against Sage explicit all-pairs certificate:',label,'minor rank',a['minor_rank'])

# The native verifier rejects a repeated support element.
inst=build_instance(PARAMS['grs-B'],5)
E,Y=inst['E'],inst['Y']
Yhat=tangent_matrix(E,Y,inst['alpha'],inst['nu'])
handle=native.Minor(E,Y,Yhat,None)
bad=list(inst['alpha']); bad[1]=bad[0]
try:
    handle.verify(bad,[E.one()]*len(bad),inst['D'])
    raise AssertionError('duplicate support accepted')
except RuntimeError as exc:
    assert 'duplicate support' in str(exc)
handle.close()
print('PASS native verifier rejects a repeated support element')

# Classic McEliece shortenings.
for name,(m,n0,t) in CLASSIC.items():
    for k in SHORTENED_DIMENSIONS[name]:
        P=PARAMS['%s-k%d'%(name,k)]
        n=k+m*t
        assert P['ell']==n0-n and comb(k,2)>n-3*t-1
        print('PASS %s shortened by %d positions to n=%d, k=%d, D=%d, D-t=%d'%(name,P['ell'],n,k,n-2*t-1,n-3*t-1))
