#!/usr/bin/env python3
"""Recompute every minor and verify its exact GRS row space with Python alone."""
import argparse,hashlib,json,time
from itertools import combinations
from pathlib import Path
from independent_verify import TABLE,power,inverse,field_rank

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('results',type=Path);args=ap.parse_args()
    out=args.results;begin=time.monotonic()
    public=json.loads((out/'public.json').read_text());payload=json.loads((out/'bootstrap-public-data.json').read_text())
    n,k,m,t=[public['parameters'][x] for x in ['n','k','m','t']]
    assert m==8 and sum(a<<i for i,a in enumerate(public['field_modulus']))==285
    G=public['generator'];V=payload['tangent_representatives'];basis=payload['minor_basis']
    alpha=payload['recovered_support'];eta=payload['recovered_minor_multiplier'];pairs=list(combinations(range(k),2))
    assert len(V)==k and all(len(row)==n for row in V)
    assert len(alpha)==len(eta)==n and len(set(alpha))==n and all(eta)
    assert all(0<=x<256 for row in V+basis for x in row)
    roots=[power(a,128) for a in range(256)]
    minors=[bytes(roots[(V[b][j] if G[a][j] else 0)^(V[a][j] if G[b][j] else 0)] for j in range(n)) for a,b in pairs]
    packed=b''.join(minors);assert packed==(out/'minor-matrix.u8').read_bytes()
    rank=field_rank([int.from_bytes(row,'little') for row in minors],n)
    assert rank==n-3*t-1==90
    assert field_rank([int.from_bytes(bytes(row),'little') for row in basis],n)==rank
    nu=[]
    for j,a in enumerate(alpha):
        dp=1
        for i,b in enumerate(alpha):
            if i!=j:dp=TABLE[dp][a^b]
        nu.append(inverse(TABLE[eta[j]][dp]))
    checks=[[TABLE[nu[j]][power(alpha[j],r)] for j in range(n)] for r in range(n-rank)]
    products_checked=0
    for row in minors+[bytes(row) for row in basis]:
        for h in checks:
            z=0
            for a,b in zip(row,h):z^=TABLE[a][b]
            assert not z;products_checked+=1
    report=dict(status='passed',standard_library_only=True,sage_imported=False,all_square_root_minors_recomputed=True,
        all_pairs=len(pairs),minor_entries_compared=len(packed),minor_rank=rank,length=n,
        saved_basis_rank=rank,distinct_support=n,all_grs_multipliers_nonzero=True,
        exact_minor_code_equals_recovered_grs=True,parity_inner_products_checked=products_checked,
        minor_matrix_sha256=hashlib.sha256(packed).hexdigest(),seconds=time.monotonic()-begin)
    (out/'minor-independent-verification.json').write_text(json.dumps(report,indent=2)+'\n')
    artifact=dict(challenge=public['challenge'],field_order=256,field_modulus=public['field_modulus'],
        length=n,dimension=rank,original_indices=public['original_indices'],basis=basis,
        support=alpha,grs_multipliers=eta,full_matrix_file='minor-matrix.u8',
        full_matrix_shape=[len(pairs),n],full_matrix_encoding='row-major, one polynomial-basis GF(256) integer byte per entry',
        minor_row_pairs=[list(pair) for pair in pairs],exact_grs_row_space_verified=True)
    (out/'minor-code.json').write_text(json.dumps(artifact,indent=2)+'\n')
    print(json.dumps(report),flush=True)

if __name__=='__main__':main()
