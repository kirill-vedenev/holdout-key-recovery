"""Convert the minor-code support to an explicitly verified alternant certificate."""
from sage.all import *
from pathlib import Path
import argparse, json

ap=argparse.ArgumentParser();ap.add_argument('results',type=Path);args=ap.parse_args()
out=args.results;public=json.loads((out/'public.json').read_text());saved=json.loads((out/'bootstrap-public-data.json').read_text())
m,n,t,k=[public['parameters'][x] for x in ['m','n','t','k']]
E=GF(2**m,'a',modulus=PolynomialRing(GF(2),'x')(public['field_modulus']));dec=E.from_integer;enc=lambda a:int(a.to_integer())
alpha=[dec(a) for a in saved['recovered_support']];lam=[dec(a) for a in saved['recovered_ambient_multiplier']]
assert len(alpha)==n and len(set(alpha))==n and all(lam)
R=PolynomialRing(E,'x');x=R.gen();Pi=prod(x-a for a in alpha)
nu=[1/(lam[j]*Pi.derivative()(alpha[j])) for j in range(n)]
G=matrix(GF(2),public['generator']);H=matrix(E,2*t,n,lambda r,j:nu[j]*alpha[j]**r)
assert not G.change_ring(E)*H.transpose()
B=matrix(GF(2),2*t*m,n,lambda r,j:(enc(H[r//m,j])>>(r%m))&1)
assert B.right_kernel()==G.row_space()
data=dict(challenge=public['challenge'],public_only=True,support=[enc(a) for a in alpha],
    dual_multipliers=[enc(a) for a in nu],ambient_multipliers=[enc(a) for a in lam],
    original_indices=public['original_indices'],field_modulus=public['field_modulus'],
    verification=dict(distinct_support=n,binary_parity_rank=B.rank(),binary_generator_rank=G.rank(),
        recovered_binary_subfield_code_equals_public=True,route='jet -> extension-field quadrics -> square-root minor code -> Sidelnikov-Shestakov'))
(out/'retained-support.json').write_text(json.dumps(data,indent=2,default=int)+'\n')
print(json.dumps(data['verification'],default=int),flush=True)
