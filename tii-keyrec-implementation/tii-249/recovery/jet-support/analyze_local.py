"""Analyze exact local Hasse tensors, select branches using public algebra only."""
from sage.all import *
from pathlib import Path
from math import comb
from itertools import combinations
import argparse,json,time
from gf2_io import panel_matrix
from public_local import public_branches,polar_matrix


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('public',type=Path);ap.add_argument('result',type=Path)
    args=ap.parse_args();begin=time.monotonic()
    data=json.loads(args.public.read_text());k=data['parameters']['k'];d=data['parameters']['degree'];count=280;words=5
    E=GF(2**data['parameters']['m'],'a',modulus=PolynomialRing(GF(2),'x')(data['field_modulus']))
    tensors=[panel_matrix(args.result/f'hasse-{u}.u64le',count,words) for u in range(d)]
    assert [M.ncols() for M in tensors]==[comb(k,u) for u in range(d)]
    assert not tensors[0]
    report=dict(public_only=True,point_original=data['heldout_original'],point_shortened=data['heldout'][0],
                variables=k,kernel_polynomials=count,monomial_order='colex',
                hasse_orders=[dict(order=u,coordinates=comb(k,u),rank=int(M.rank()),nonzero_forms=int(sum(bool(row) for row in M.rows()))) for u,M in enumerate(tensors)],
                degree_five_rank=count,degree_five_tensor='input kernel-panel.u64le, no translation needed',
                repeated_variable_hasse_coefficients='all zero in the original squarefree basis',
                field_modulus=data['field_modulus'])
    G,Q=tensors[1:3]
    p=vector(GF(2),[(data['heldout_mask']>>a)&1 for a in range(k)])
    assert p==matrix(GF(2),data['generator']).column(data['heldout'][0])
    pairs=sorted(combinations(range(k),2),key=lambda pair:sum(comb(a,j+1) for j,a in enumerate(pair)))
    W=G.right_kernel();report['gradient_annihilator_dimension']=W.dimension()
    report['gradient_annihilator_basis']=[[int(x) for x in v] for v in W.basis()]
    report['gradient_annihilates_point']=not G*p
    print(json.dumps({a:b for a,b in report.items() if a!='gradient_annihilator_basis'},default=int),flush=True)
    (args.result/'local-ranks.json').write_text(json.dumps(report,indent=2,default=int)+'\n')
    branches,branch_report=public_branches(G,Q,p,E,pairs)
    report['branches']=branch_report
    report['branch_vectors']=[[int(x.to_integer()) for x in v] for v in branches]
    report['continuation_checks']=[]
    for i,v in enumerate(branches):
        J=G.change_ring(E);L=polar_matrix(Q,v,pairs);A=J.stack(L)
        report['continuation_checks'].append(dict(branch=i,rank=A.rank(),nullity=k-A.rank(),
            kernel_equals_point_tangent_plane=A.right_kernel()==matrix(E,[p.change_ring(E),v]).row_space()))
    report['seconds']=time.monotonic()-begin
    (args.result/'local-branches.json').write_text(json.dumps(report,indent=2,default=int)+'\n')
    print(json.dumps({a:b for a,b in report.items() if a not in ['gradient_annihilator_basis','branch_vectors']},default=int),flush=True)

if __name__=='__main__':main()
