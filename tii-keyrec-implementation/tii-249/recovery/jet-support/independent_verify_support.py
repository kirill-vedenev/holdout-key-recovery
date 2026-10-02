#!/usr/bin/env python3
"""Standard-library-only verification of the extracted projective support.

Reads only the shortened public generator and newly recovered candidate. No
Sage or recovery module is imported. Inclusion plus exact dimensions proves
equality of the binary subfield code and public row space.
"""
import argparse,hashlib,json,time
from pathlib import Path


def binary_rank(rows):
    pivots={}
    for row in rows:
        while row:
            p=row.bit_length()-1
            if p in pivots:row^=pivots[p]
            else:pivots[p]=row;break
    return len(pivots)


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('public',type=Path);ap.add_argument('candidate',type=Path);ap.add_argument('--report',type=Path,required=True)
    ap.add_argument('--original-public-key',type=Path)
    args=ap.parse_args();begin=time.monotonic();public=json.loads(args.public.read_text());saved=json.loads(args.candidate.read_text())
    m,n,t,k=[public['parameters'][x] for x in ['m','n','t','k']];q=1<<m;modulus=sum(x<<i for i,x in enumerate(public['field_modulus']))
    def mul(a,b):
        out=0
        while b:
            if b&1:out^=a
            b>>=1;a<<=1
            if a&q:a^=modulus
        return out
    table=[[mul(a,b) for b in range(q)] for a in range(q)]
    def power(a,r):
        out=1
        while r:
            if r&1:out=table[out][a]
            a=table[a][a];r>>=1
        return out
    def inverse(a):
        assert a;return power(a,q-2)
    support=saved['support'];nu=saved['dual_multipliers'];G=public['generator'];pairs=saved['quadratic_pairs'];forms=saved['forms']
    assert len(G)==k and all(len(row)==n and all(x in (0,1) for x in row) for row in G)
    assert len(support)==len(nu)==n and len(set(support))==n
    assert all(x is None or isinstance(x,int) and 0<=x<q for x in support)
    assert all(isinstance(x,int) and 0<x<q for x in nu)
    assert saved['original_indices']==public['original_indices']
    seed=public['heldout'][0];ratios=[]
    for j in range(n):
        values=[0,0]
        for row,form in enumerate(forms):
            for a,(u,v) in enumerate(pairs):
                if G[u][j] and G[v][j]:values[row]^=form[a]
        high,lower=values
        if j==seed:assert high==lower==0;alpha=0
        else:
            assert high;alpha=power(table[high][inverse(lower)],q//2) if lower else None
        assert alpha==support[j];ratios.append(alpha)
    H=[[table[nu[j]][int(r==2*t-1) if support[j] is None else power(support[j],r)] for j in range(n)] for r in range(2*t)]
    for row in G:
        for check in H:
            value=0
            for x,y in zip(row,check):
                if x:value^=y
            assert value==0
    parity=[sum(((H[r][j]>>bit)&1)<<j for j in range(n)) for r in range(2*t) for bit in range(m)]
    rank_H=binary_rank(parity);rank_G=binary_rank([sum(x<<j for j,x in enumerate(row)) for row in G]);assert rank_G==k and rank_H==n-k
    original_checks={}
    if args.original_public_key:
        lines=[line.strip() for line in args.original_public_key.read_text().splitlines() if line.strip()]
        full_rows=[];original_n=public['original_parameters']['n'];keep=public['original_indices']
        assert len(set(keep))==n and all(0<=j<original_n for j in keep)
        for line in lines[:-1]:
            assert line.startswith('[') and line.endswith(']')
            bits=[int(x) for x in line[1:-1].split()];assert len(bits)==original_n and all(x in (0,1) for x in bits)
            full_rows.append(sum(bits[j]<<a for a,j in enumerate(keep)))
        assert json.loads(lines[-1])==public['field_modulus']
        g_rows=[sum(x<<j for j,x in enumerate(row)) for row in G]
        assert all((h&g).bit_count()%2==0 for h in full_rows for g in g_rows)
        assert binary_rank(full_rows)==n-k
        assert keep[seed]==public['heldout_original'] and sum(G[a][seed]<<a for a in range(k))==public['heldout_mask']
        original_sha=hashlib.sha256(args.original_public_key.read_bytes()).hexdigest();assert original_sha==public['public_sha256']
        original_checks=dict(original_public_key_sha256=original_sha,
            shortened_generator_equals_original_public_code_shortened_at_saved_positions=True,
            original_public_parity_check_restriction_rank=binary_rank(full_rows))
    report=dict(status='passed',standard_library_only=True,independent_of_sage_and_recovery_modules=True,
        public_sha256=hashlib.sha256(args.public.read_bytes()).hexdigest(),candidate_sha256=hashlib.sha256(args.candidate.read_bytes()).hexdigest(),
        public_only=True,published_secret_key_read=False,field_modulus_integer=modulus,positions=n,
        distinct_projective_support=True,all_dual_multipliers_nonzero=True,quadratic_section_ratios_reproduce_saved_support=True,
        shortened_to_original_order_verified=True,all_public_alternant_equations_zero=True,
        binary_parity_rank=rank_H,public_generator_rank=rank_G,recovered_binary_subfield_code_equals_public=True,seconds=time.monotonic()-begin,**original_checks)
    args.report.write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report),flush=True)

if __name__=='__main__':main()
