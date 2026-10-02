"""Public binary-normalized continuation, following the paper's recurrence.

Reads only public.json and kernel.bin. Fixes v_(2r)=v_r^2, then solves the
regular odd-order and ALL next-even-order equations with the fixed [J; L].
"""
from sage.all import *
from pathlib import Path
from itertools import combinations
import argparse,json,time
import numpy as np
from core import read_kernel,jets_at_point,apply_binary,products
from incremental import SeriesCache
from public_local import public_branches,polar_matrix


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('directory',type=Path)
    ap.add_argument('--position',type=int)
    ap.add_argument('--depth',type=int,default=512)
    args=ap.parse_args()
    if args.depth<2 or args.depth%2:
        ap.error('depth must be an even integer at least two')
    folder=args.directory;started=time.monotonic()
    data=json.loads((folder/'public.json').read_text())
    k,d,K=read_kernel(folder/'kernel.bin')
    P=matrix(GF(2),data['generator'])
    i=data['heldout'][0] if args.position is None else args.position
    E=GF(2**data['parameters']['m'],'a',modulus=PolynomialRing(GF(2),'x')(data['field_modulus']))
    p=P.column(i);pairs=list(combinations(range(k),2))
    local=K*jets_at_point(k,d,p).transpose()
    if not K.nrows() or local.column(0):raise ValueError('empty or rejected holdout')
    G=local.matrix_from_columns(range(1,k+1))
    Q=local.matrix_from_columns(range(k+1,local.ncols()))
    branches,branch_report=public_branches(G,Q,p,E,pairs)
    if not branches:raise ValueError(f'first-jet recovery failed: {branch_report}')
    v1=branches[0];J=G.change_ring(E);L=polar_matrix(Q,v1,pairs)
    A=J.stack(L);plane=matrix(E,[p,v1]).row_space()
    if A.right_kernel()!=plane:
        raise ValueError('normalized continuation does not have exactly the two allowed freedoms')
    # Factor the fixed solve once. The two non-pivot coordinates are set to zero.
    columns=A.pivots();independent=A.matrix_from_columns(columns)
    rows=independent.transpose().pivots()
    inverse=independent.matrix_from_rows(rows).inverse()
    def solve(rhs):
        coefficients=inverse*vector(E,[rhs[r] for r in rows])
        answer=vector(E,k)
        for j,c in zip(columns,coefficients):answer[j]=c
        assert A*answer==rhs
        return answer
    jets=[p.change_ring(E),v1,vector(E,[x*x for x in v1])]
    cache=SeriesCache(k,d,E,args.depth)
    try:
        for r,v in enumerate(jets):cache.set(r,v)
        for r in range(3):
            assert not apply_binary(K,cache.evaluate(r).reshape((-1,1)),E)
        history=[];first_full=None
        for r in range(3,args.depth,2):
            even=vector(E,[x*x for x in jets[(r+1)//2]])
            constant=apply_binary(K,np.column_stack((cache.evaluate(r),cache.evaluate(r+1))),E)
            rhs=-vector(E,list(constant.column(0))+list(constant.column(1)+J*even))
            odd=solve(rhs)
            jets.extend([odd,even]);cache.set(r,odd);cache.set(r+1,even)
            actual=apply_binary(K,np.column_stack((cache.evaluate(r),cache.evaluate(r+1))),E)
            assert not actual
            odd_rank=matrix(E,jets[:-1]).rank();even_rank=matrix(E,jets).rank()
            if first_full is None:
                if odd_rank==k:first_full=r
                elif even_rank==k:first_full=r+1
            history.append(dict(odd_order=r,even_order=r+1,
                                odd_flag_dimension=odd_rank,even_flag_dimension=even_rank))
            if (r+1)%64==0:
                print(json.dumps(dict(event='normalized_extended',order=r+1,
                    flag_dimension=even_rank,seconds=time.monotonic()-started)),flush=True)
    finally:
        cache.close()
    recursion_seconds=time.monotonic()-started
    assert len(jets)==args.depth+1
    assert all(jets[2*r]==vector(E,[x*x for x in jets[r]]) for r in range(1,args.depth//2+1))
    # Fresh multiplication, independent of the incremental product cache.
    composed=apply_binary(K,products([list(row) for row in zip(*jets)],d,E,args.depth),E)
    assert not composed
    report=dict(public_only=True,input_files=['public.json','kernel.bin'],
        binary_normalization=True,even_coefficients_fixed_by_squaring=True,
        all_next_even_order_equations_used=True,stationary_projection_used=False,
        position=i,depth=args.depth,kernel_dimension=K.nrows(),branches=branch_report,
        fixed_matrix='[J; L]',fixed_matrix_factorizations=1,fixed_continuation_rank=A.rank(),
        affine_nullity=k-A.rank(),ambiguity_equals_first_plane=True,
        odd_linear_solves=(args.depth-2)//2,even_square_relations_checked=args.depth//2,
        all_even_square_relations_verified=True,all_recurrence_equations_verified=True,
        independent_full_substitution_verified=True,first_full_flag_order=first_full,
        field_modulus=data['field_modulus'],history=history,
        jets=[[int(x.to_integer()) for x in v] for v in jets],
        recursion_seconds=recursion_seconds,seconds=time.monotonic()-started)
    (folder/f'deep-jets-{i}.json').write_text(json.dumps(report,indent=2,default=int)+'\n')
    print(json.dumps({a:b for a,b in report.items() if a not in ('history','jets')},default=int),flush=True)


if __name__=='__main__':main()
