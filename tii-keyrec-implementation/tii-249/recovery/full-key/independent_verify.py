#!/usr/bin/env python3
"""Independently verify the full TII-249 support, multipliers, and Goppa key.

Only Python's standard library is used. No recovery modules or secret reference
key are imported. Field elements use the original public polynomial basis.
"""
import argparse,hashlib,json,time
from pathlib import Path

PUBLIC_SHA256='63ef408ba6f954272e43a0563dff373064a70c24601a84611a00d977d558b598'
N,M,T,K,MODULUS=235,8,16,107,285


def mul(a,b):
    result=0
    while b:
        if b&1:result^=a
        b>>=1;a<<=1
        if a&256:a^=MODULUS
    return result


TABLE=[[mul(a,b) for b in range(256)] for a in range(256)]


def power(a,r):
    result=1
    while r:
        if r&1:result=TABLE[result][a]
        a=TABLE[a][a];r>>=1
    return result


def inverse(a):
    if not a:raise ValueError('division by zero')
    return power(a,254)


def evaluate(p,a):
    value=0
    for c in reversed(p):value=TABLE[value][a]^c
    return value


def trim(p):
    p=list(p)
    while p and not p[-1]:p.pop()
    return p


def poly_add(a,b):
    result=list(a)+[0]*max(0,len(b)-len(a))
    for i,c in enumerate(b):result[i]^=c
    return trim(result)


def remainder(a,b):
    a,b=trim(a),trim(b)
    if not b:raise ValueError('zero polynomial divisor')
    leading_inverse=inverse(b[-1])
    while a and len(a)>=len(b):
        shift=len(a)-len(b);c=TABLE[a[-1]][leading_inverse]
        for i,v in enumerate(b):a[shift+i]^=TABLE[c][v]
        a=trim(a)
    return a


def gcd(a,b):
    a,b=trim(a),trim(b)
    while b:a,b=b,remainder(a,b)
    return [TABLE[c][inverse(a[-1])] for c in a] if a else []


def square_mod(a,g):
    result=[0]*(2*len(a))
    for i,c in enumerate(a):result[2*i]=TABLE[c][c]
    return remainder(result,g)


def irreducible_degree_16(g):
    # Exact Rabin criterion: the sole prime divisor of degree 16 is 2.
    x=[0,1];h=x
    for r in range(1,17):
        for _ in range(M):h=square_mod(h,g)
        if r==8 and gcd(poly_add(h,x),g)!=[1]:return False
    return not poly_add(h,x)


def rref(rows,n):
    rows=list(rows);rank=0
    for c in range(n):
        pivot=next((j for j in range(rank,len(rows)) if (rows[j]>>c)&1),None)
        if pivot is None:continue
        rows[rank],rows[pivot]=rows[pivot],rows[rank]
        for j in range(len(rows)):
            if j!=rank and ((rows[j]>>c)&1):rows[j]^=rows[rank]
        rank+=1
    return tuple(rows[:rank])


def binary_check(alpha,multipliers,dimension):
    rows=[0]*(M*dimension)
    for j,(a,nu) in enumerate(zip(alpha,multipliers)):
        value=nu
        for degree in range(dimension):
            for bit in range(M):rows[degree*M+bit]|=((value>>bit)&1)<<j
            value=TABLE[value][a]
    return rows


def field_rank(packed_rows,width):
    """Exact GF(256) elimination using one byte per packed field coefficient."""
    rows=list(packed_rows);rank=0
    for c in range(width):
        shift=8*c;pivot=next((j for j in range(rank,len(rows)) if (rows[j]>>shift)&255),None)
        if pivot is None:continue
        rows[rank],rows[pivot]=rows[pivot],rows[rank]
        coefficient=(rows[rank]>>shift)&255
        normalized=rows[rank].to_bytes(width,'little').translate(bytes(TABLE[inverse(coefficient)]))
        rows[rank]=int.from_bytes(normalized,'little');multiples={1:rows[rank]}
        for j in range(rank+1,len(rows)):
            coefficient=(rows[j]>>shift)&255
            if not coefficient:continue
            if coefficient not in multiples:multiples[coefficient]=int.from_bytes(normalized.translate(bytes(TABLE[coefficient])),'little')
            rows[j]^=multiples[coefficient]
        rank+=1
        if rank==width:break
    return rank


def verify_extensions(key,original,shortened_public,evidence_dir):
    short=json.loads(shortened_public.read_text());keep=short['original_indices'];removed=short['shortened_original_indices'];G=short['generator']
    if len(keep)!=178 or len(removed)!=57 or sorted(keep+removed)!=list(range(N)):raise ValueError('invalid shortening map')
    if len(G)!=50 or any(len(row)!=178 or any(x not in (0,1) for x in row) for row in G):raise ValueError('invalid shortened binary generator')
    restricted_known=[sum(((row>>b)&1)<<a for a,b in enumerate(keep)) for row in original]
    packed_G=[sum(x<<a for a,x in enumerate(row)) for row in G]
    if len(rref(restricted_known,178))!=128 or len(rref(packed_G,178))!=50:raise ValueError('invalid shortened public dimensions')
    if any((h&g).bit_count()%2 for h in restricted_known for g in packed_G):raise ValueError('shortened generator is outside the original public code')
    alpha=[key['support'][j] for j in keep];nu=[key['dual_grs_multipliers'][j] for j in keep]
    if alpha!=key['affine_chart']['retained_support']:raise ValueError('affine retained support mismatch')
    powers=[[power(a,r) for a in alpha] for r in range(2*T)];common=[]
    for row in G:
        for r in range(2*T):common.append(int.from_bytes(bytes(powers[r][j] if row[j] else 0 for j in range(len(keep))),'little'))
    common_rank=field_rank(common,len(keep))
    if common_rank!=177:raise ValueError('retained multiplier block is not one-dimensional')
    for row in common:
        value=0;coefficients=row.to_bytes(len(keep),'little')
        for a,b in zip(coefficients,nu):value^=TABLE[a][b]
        if value:raise ValueError('retained multipliers do not satisfy every shortened equation')
    candidates=[a for a in range(256) if a not in set(alpha)]+[None];total=0
    for j in removed:
        record=json.loads((evidence_dir/f'position-{j}.json').read_text());positions=keep+[j]
        restricted=[sum(((row>>b)&1)<<a for a,b in enumerate(positions)) for row in original];reduced=rref(restricted,len(positions))
        if len(reduced)!=128:raise ValueError('one-position extension parity rank changed')
        c=1<<len(keep)
        for row in reduced:
            pivot=(row&-row).bit_length()-1
            if pivot>=len(keep):raise ValueError('new extension coordinate is not free')
            if (row>>len(keep))&1:c|=1<<pivot
        if any((row&c).bit_count()%2 for row in restricted):raise ValueError('independent extension codeword is outside the public code')
        moments=[]
        for r in range(2*T):
            value=0
            for a in range(len(keep)):
                if (c>>a)&1:value^=TABLE[nu[a]][powers[r][a]]
            moments.append(value)
        if moments!=record['moments'] or not moments:raise ValueError('saved extension moments differ from public recomputation')
        pivot=next(a for a,x in enumerate(moments) if x);accepted=[]
        if [test['value'] for test in record['tests']]!=candidates:raise ValueError('candidate enumeration is incomplete')
        for test,w in zip(record['tests'],candidates):
            target=[int(r==2*T-1) if w is None else power(w,r) for r in range(2*T)]
            scalar=TABLE[target[pivot]][inverse(moments[pivot])]
            residual=[TABLE[scalar][a]^b for a,b in zip(moments,target)];consistent=not any(residual);full=consistent and bool(scalar)
            if residual!=test['moment_residual'] or test['retained_multiplier_scalar']!=scalar:raise ValueError('candidate residual differs from saved evidence')
            if test['kernel_dimension']!=int(consistent) or test['has_full_support_multiplier']!=full:raise ValueError('incorrect candidate outcome')
            if full:accepted.append(w)
            total+=1
        if accepted!=[key['support'][j]] or record['accepted_value']!=key['support'][j]:raise ValueError('extension candidate is not uniquely the recovered original coordinate')
        if moments[0]!=key['dual_grs_multipliers'][j]:raise ValueError('zeroth moment does not recover the missing multiplier')
        if moments[1]!=TABLE[moments[0]][key['support'][j]]:raise ValueError('first moment does not recover missing support')
        scalar=inverse(moments[0]);winner=[TABLE[scalar][a] for a in nu]+[1]
        if winner!=record['winner_full_multiplier']:raise ValueError('winning extension multiplier vector differs')
    return dict(independent_candidate_reenumeration=True,common_retained_multiplier_rank=common_rank,
        shortened_generator_checked_against_original_public_key=True,
        restored_positions_verified=len(removed),projective_candidates_verified=total,
        all_restored_positions_have_exactly_one_candidate=True,all_projective_infinity_candidates_rejected=True,
        missing_multipliers_equal_zeroth_public_moments=True,missing_support_equals_first_over_zeroth_public_moment=True)


def verify(candidate,public,shortened_public=None,evidence_dir=None):
    begin=time.monotonic();original_raw=public.read_bytes();candidate_raw=candidate.read_bytes()
    if hashlib.sha256(original_raw).hexdigest()!=PUBLIC_SHA256:raise ValueError('original public-key SHA-256 mismatch')
    lines=[line.strip() for line in original_raw.decode().splitlines() if line.strip()]
    if len(lines)!=129:raise ValueError('expected 128 public rows plus field modulus')
    original=[]
    for line in lines[:-1]:
        if not line.startswith('[') or not line.endswith(']'):raise ValueError('invalid public row')
        bits=[int(x) for x in line[1:-1].split()]
        if len(bits)!=N or any(x not in (0,1) for x in bits):raise ValueError('invalid public row width')
        original.append(sum(x<<j for j,x in enumerate(bits)))
    expected=rref(original,N)
    if len(expected)!=128:raise ValueError('original parity rank is not 128')
    key=json.loads(candidate_raw);modulus=json.loads(lines[-1])
    if key['field_modulus']!=modulus or sum(x<<i for i,x in enumerate(modulus))!=MODULUS:raise ValueError('field encoding mismatch')
    alpha,nu,lam=key['support'],key['dual_grs_multipliers'],key['ambient_grs_multipliers']
    if key['original_indices']!=list(range(N)):raise ValueError('original coordinate ordering was lost')
    if any(len(v)!=N for v in [alpha,nu,lam]) or len(set(alpha))!=N:raise ValueError('invalid full support length or duplicates')
    if any(type(x) is not int or not 0<=x<256 for x in alpha+nu+lam):raise ValueError('invalid field element encoding')
    if any(not x for x in nu+lam) or nu[0]!=1:raise ValueError('invalid zero multiplier or normalization')
    actual_alternant=rref(binary_check(alpha,nu,2*T),N)
    if actual_alternant!=expected:raise ValueError('full alternant binary row space differs from the original public code')
    derivative_values=[]
    for i,a in enumerate(alpha):
        derivative=1
        for j,b in enumerate(alpha):
            if i!=j:derivative=TABLE[derivative][a^b]
        if TABLE[lam[i]][TABLE[nu[i]][derivative]]!=1:raise ValueError('primal/dual multiplier convention mismatch')
        derivative_values.append(derivative)
    # Check every possible exponent in the ambient/dual orthogonality product.
    weighted=[TABLE[x][y] for x,y in zip(lam,nu)];powers=[1]*N
    for degree in range(N-1):
        value=0
        for j in range(N):value^=TABLE[weighted[j]][powers[j]];powers[j]=TABLE[powers[j]][alpha[j]]
        if value:raise ValueError('ambient GRS and dual check fail orthogonality')
    report=dict(schema='tii249-independent-full-key-verification-v1',status='passed',public_only=True,
        standard_library_only=True,recovery_modules_imported=False,secret_reference_used=False,
        original_public_key_sha256=PUBLIC_SHA256,candidate_sha256=hashlib.sha256(candidate_raw).hexdigest(),
        original_order_verified=True,support_count=N,support_finite=True,support_distinct=True,
        all_dual_and_primal_multipliers_nonzero=True,dual_multiplier_at_original_zero=1,
        primal_dual_derivative_relation_verified=True,ambient_dual_orthogonality_exponents_checked=N-1,
        ambient_grs_dimension=N-2*T,ambient_dimension_certificate='nonzero-scaled Vandermonde evaluation at 235 distinct field elements',
        reconstructed_alternant_binary_rank=len(actual_alternant),original_binary_rank=len(expected),
        alternant_public_row_space_exactly_equal=True,original_binary_code_dimension=K)
    if 'g' in key:
        g=key['g']
        if len(g)!=T+1 or any(type(x) is not int or not 0<=x<256 for x in g) or g[-1]!=1:raise ValueError('expected a monic degree-16 Goppa polynomial')
        if not irreducible_degree_16(g):raise ValueError('Goppa polynomial is reducible')
        values=[evaluate(g,a) for a in alpha]
        if any(not x for x in values):raise ValueError('Goppa polynomial vanishes on support')
        eta=[inverse(x) for x in values]
        if key.get('goppa_parity_multipliers')!=eta:raise ValueError('Goppa parity multiplier mismatch')
        actual_goppa=rref(binary_check(alpha,eta,T),N)
        if actual_goppa!=expected:raise ValueError('Goppa binary row space differs from original public code')
        scale=key['reciprocal_multiplier_scale'];square=key['reciprocal_multiplier_polynomial']
        if len(square)!=2*T+1 or not scale or square[-1]!=scale:raise ValueError('invalid reciprocal multiplier square')
        target=[0]*(2*T+1)
        for i,c in enumerate(g):target[2*i]=TABLE[scale][TABLE[c][c]]
        if square!=target:raise ValueError('reciprocal multiplier polynomial is not the claimed square')
        for a,b,gx in zip(alpha,nu,values):
            if evaluate(square,a)!=inverse(b) or TABLE[b][TABLE[scale][TABLE[gx][gx]]]!=1:raise ValueError('reciprocal multiplier/Goppa identity fails')
        report.update(goppa_degree=T,goppa_monic=True,goppa_irreducible=True,goppa_nonzero_on_support=True,
            reconstructed_goppa_binary_rank=len(actual_goppa),goppa_public_row_space_exactly_equal=True,
            reciprocal_multiplier_square_identity_verified=True,verified_full_equivalent_goppa_key=True)
    else:report.update(verified_full_support_and_multipliers=True,goppa_polynomial_absent=True)
    if shortened_public is not None or evidence_dir is not None:
        if shortened_public is None or evidence_dir is None:raise ValueError('supply both shortened public input and extension evidence')
        report.update(verify_extensions(key,original,shortened_public,evidence_dir))
    report['seconds']=time.monotonic()-begin;return report


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('candidate',type=Path);ap.add_argument('--public',type=Path,required=True);ap.add_argument('--report',type=Path,required=True)
    ap.add_argument('--shortened-public',type=Path);ap.add_argument('--extension-evidence',type=Path)
    args=ap.parse_args();report=verify(args.candidate,args.public,args.shortened_public,args.extension_evidence);args.report.write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2),flush=True)

if __name__=='__main__':main()
