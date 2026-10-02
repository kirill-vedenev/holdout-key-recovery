"""Recompute every saved jet coefficient by full convolution, without patch updates."""
from sage.all import *
from pathlib import Path
from extend_jet import Cache
import argparse,json,time

ap=argparse.ArgumentParser();ap.add_argument('results',type=Path);ap.add_argument('panel',type=Path);args=ap.parse_args()
out=args.results;public=json.loads((out/'public.json').read_text());saved=json.loads((out/'deep-jets-0.json').read_text())
k,d,m=[public['parameters'][x] for x in ['k','degree','m']];count=saved['kernel_dimension'];depth=saved['depth']
E=GF(2**m,'a',modulus=PolynomialRing(GF(2),'x')(public['field_modulus']))
jets=[vector(E,[E.from_integer(a) for a in v]) for v in saved['jets']]
assert len(jets)==depth+1
start=time.monotonic();cache=Cache(k,d,E,depth,count,4,args.panel)
try:
    for r,v in enumerate(jets):cache.set(r,v)
    for r in range(depth+1):
        residual=cache.evaluate(r)
        assert not residual,('nonzero coefficient in complete substitution',r)
finally:cache.close()
assert all(jets[2*r]==vector(E,[a*a for a in jets[r]]) for r in range(1,depth//2+1))
report=dict(status='passed',full_convolution_recomputed=True,incremental_patch_updates_used=False,
    polynomial_count=count,jet_depth=depth,coefficient_vectors_checked=depth+1,
    extension_field_equalities_checked=count*(depth+1),all_residuals_zero=True,
    all_even_square_relations_verified=True,seconds=time.monotonic()-start)
(out/'full-jet-substitution-verification.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report),flush=True)
