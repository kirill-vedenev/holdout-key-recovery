"""Only public branch selection and local polar maps, extracted from audit.py."""
from sage.all import *
from itertools import combinations, combinations_with_replacement
from pathlib import Path
import argparse, json, time
from gf2_io import binary_times_field


def quadratic_restriction(Q, basis, pairs):
    E=basis.base_ring(); b=basis.nrows()
    columns=[]
    for u,v in combinations_with_replacement(range(b),2):
        columns.append([basis[u,a]*basis[u,c] if u==v else
                        basis[u,a]*basis[v,c]+basis[v,a]*basis[u,c] for a,c in pairs])
    return binary_times_field(Q,matrix(E,columns).transpose())


def public_branches(G,Q,p,E,pairs):
    W=G.right_kernel(); basis=[p]
    for row in W.basis():
        if row not in span(basis,GF(2)): basis.append(row)
    B=matrix(GF(2),basis[1:]); b=B.nrows()
    if b>8: return [],dict(status='quotient_dimension_exceeds_enumeration_limit',dimension=b)
    R=PolynomialRing(E,b,'z',order='lex'); z=R.gens()
    monomials=[z[u]*z[v] for u,v in combinations_with_replacement(range(b),2)]
    restrictions=quadratic_restriction(Q,B,pairs).row_space().basis_matrix()
    equations=[sum((E(c)*mon for c,mon in zip(row,monomials)),R.zero()) for row in restrictions]
    solutions=[]; dimensions=[]
    for chart in range(b):
        I=R.ideal(equations+[z[a] for a in range(chart)]+[z[chart]-1])
        dim=I.dimension(); dimensions.append(dim)
        if dim>0:
            return [],dict(status='positive_dimensional_branch_locus',chart=chart,dimension=dim)
        if dim==0:
            for sol in I.variety(ring=E):
                coordinates=vector(E,[sol[x] for x in z])
                v=coordinates*B.change_ring(E)
                assert not G.change_ring(E)*v
                assert not binary_times_field(Q,vector(E,[v[a]*v[c] for a,c in pairs]))
                solutions.append(v)
    return solutions,dict(status='complete',quotient_dimension=b,quadratic_rank=restrictions.nrows(),
                          projective_chart_dimensions=dimensions,branches=len(solutions))


def polar_matrix(Q,v,pairs):
    E=v.base_ring(); k=len(v)
    polar=matrix(E,len(pairs),k)
    for row,(a,b) in enumerate(pairs): polar[row,a],polar[row,b]=v[b],v[a]
    return binary_times_field(Q,polar)
