"""Direct quadratic-section support extraction using only orders 0..2D.

The recovered projective support is checked by an independent public alternant
system; no curve coefficients or previously known support are input.
"""
from sage.all import *
from itertools import combinations_with_replacement
from pathlib import Path
import argparse,json,time
from direct_jet_recovery import support_from_sections


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('public',type=Path);ap.add_argument('result',type=Path)
    args=ap.parse_args();begin=time.monotonic();public=json.loads(args.public.read_text());saved=json.loads((args.result/'deep-jets-0.json').read_text())
    m,n,t,k=[public['parameters'][x] for x in ['m','n','t','k']];D=n-2*t-1
    E=GF(2**m,'a',modulus=PolynomialRing(GF(2),'x')(public['field_modulus']));enc=lambda x:int(x.to_integer())
    jets=[vector(E,[E.from_integer(x) for x in row]) for row in saved['jets'][:2*D+1]]
    if len(jets)!=2*D+1:raise ValueError('insufficient jet depth')
    R=PolynomialRing(E,'u');arcs=[R([row[a] for row in jets]) for a in range(k)]
    pairs=list(combinations_with_replacement(range(k),2));prods=[arcs[a]*arcs[b] for a,b in pairs]
    contact=matrix(E,2*D+1,len(pairs),lambda r,j:prods[j][r])
    target=matrix(E,2*D+1,2,lambda r,j:int(r==(2*D if j==0 else 2*D-2)))
    report=dict(public_only=True,public_file=str(args.public.resolve()),jet_file=str((args.result/'deep-jets-0.json').resolve()),
        published_secret_key_read=False,known_support_read=False,m=m,n=n,t=t,k=k,D=D,jet_order_used=2*D,form_degree=2,
        contact_rows=2*D+1,quadratic_form_unknowns=len(pairs),product_rank=contact.rank(),
        odd_coefficient_rank=contact.matrix_from_rows(range(1,2*D+1,2)).rank())
    try:solutions=contact.solve_right(target)
    except ValueError:
        report.update(status='extreme_sections_absent',seconds=time.monotonic()-begin)
        (args.result/'direct-support.json').write_text(json.dumps(dict(report=report),indent=2,default=int)+'\n');print(json.dumps(report,default=int));return
    assert contact*solutions==target
    forms=solutions.transpose();P=matrix(E,public['generator'])
    support=support_from_sections(E,P,public['heldout'][0],forms,pairs,2)
    assert len(support)==n and len(set(support))==n
    report.update(extreme_sections_found=True,section_targets_exact=True,distinct_projective_support=True,
        infinite_support_entries=support.count(None),section_seconds=time.monotonic()-begin)
    print(json.dumps(dict(event='direct_support_extracted',**report),default=int),flush=True)
    r=2*t
    V=matrix(E,r,n,lambda a,j:int(a==r-1) if support[j] is None else support[j]**a)
    equations=matrix(E,k*r,n,lambda a,j:P[a//r,j]*V[a%r,j]);space=equations.right_kernel()
    report['public_multiplier_system']=dict(rows=k*r,columns=n,rank=equations.rank(),nullity=space.dimension())
    if space.dimension()!=1:
        report.update(status='multiplier_space_requires_further_resolution',seconds=time.monotonic()-begin)
        result=dict(report=report,support=[None if x is None else enc(x) for x in support],forms=[[enc(x) for x in row] for row in forms],
            multiplier_space_basis=[[enc(x) for x in row] for row in space.basis()],original_indices=public['original_indices'])
    else:
        nu=space.basis()[0]
        if any(not x for x in nu):raise ValueError('public multiplier solution has zero entries')
        H=V*diagonal_matrix(E,nu);assert not P*H.transpose()
        binary=matrix(GF(2),r*m,n,lambda a,j:(enc(H[a//m,j])>>(a%m))&1)
        G=matrix(GF(2),public['generator']);equal=binary.right_kernel()==G.row_space();assert equal
        report.update(status='verified',public_alternant_equations_verified=True,all_multipliers_nonzero=True,
            recovered_binary_subfield_code_equals_public=True,binary_parity_rank=binary.rank(),public_binary_dimension=G.rank(),
            ambient_grs_dimension=n-H.rank(),seconds=time.monotonic()-begin,
            ambiguity='PGL(2,GF(256)) coordinate change, global multiplier scaling, and a common Frobenius automorphism')
        result=dict(report=report,support=[None if x is None else enc(x) for x in support],
            dual_multipliers=[enc(x) for x in nu],forms=[[enc(x) for x in row] for row in forms],quadratic_pairs=[list(pair) for pair in pairs],
            original_indices=public['original_indices'],field_modulus=public['field_modulus'])
    (args.result/'direct-support.json').write_text(json.dumps(result,indent=2,default=int)+'\n');print(json.dumps(report,default=int),flush=True)

if __name__=='__main__':main()
