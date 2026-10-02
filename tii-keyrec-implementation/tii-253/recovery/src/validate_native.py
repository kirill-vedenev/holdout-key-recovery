"""Independent exact small-instance checks of the new colex series cache."""
from sage.all import *
from pathlib import Path
from itertools import combinations
import argparse,json,random,tempfile
import numpy as np
from extend_jet import Cache


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('result',type=Path);args=ap.parse_args()
    rng=random.Random(25320261001);E=GF(256,'a',modulus=PolynomialRing(GF(2),'x')([1,0,1,1,1,0,0,0,1]));R=PolynomialRing(E,'u')
    k,d,count,depth=8,6,11,18;mons=sorted(combinations(range(k),d),key=lambda mon:sum(1<<a for a in mon))
    panel=np.asarray([rng.randrange(1<<count) for _ in mons],dtype='<u8')
    with tempfile.TemporaryDirectory(dir=args.result) as temp:
        path=Path(temp)/'panel.u64le';panel.tofile(path);cache=Cache(k,d,E,depth,count,3,path)
        arcs=[[E.zero()]*(depth+1) for _ in range(k)];checks=0
        def expected(r):
            products=[prod(R(arcs[a]) for a in mon)[r] for mon in mons]
            return vector(E,[sum((value for c,value in enumerate(products) if (int(panel[c])>>j)&1),E.zero()) for j in range(count)])
        try:
            for r in range(3):
                v=vector(E,[E.from_integer(rng.randrange(256)) for _ in range(k)]);cache.set(r,v)
                for a in range(k):arcs[a][r]=v[a]
            for r in range(3):assert cache.evaluate(r)==expected(r);checks+=1
            for r in range(3,depth,2):
                before=cache.pair(r);assert before==[expected(r),expected(r+1)];checks+=2
                odd=vector(E,[E.from_integer(rng.randrange(256)) for _ in range(k)])
                even=vector(E,[E.from_integer(rng.randrange(256)) for _ in range(k)])
                after=cache.patch(r,odd,even)
                for a in range(k):arcs[a][r],arcs[a][r+1]=odd[a],even[a]
                assert after==[expected(r),expected(r+1)];checks+=2
                assert cache.pair(r)==after;checks+=2
        finally:cache.close()
    report=dict(status='passed',independent_sage_polynomial_multiplication=True,k=k,degree=d,polynomials=count,
        depth=depth,checked_coefficient_vectors=checks,checked_field_equalities=checks*count,arbitrary_non_normalized_odd_and_even_corrections=True)
    (args.result/'native-validation.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report),flush=True)

if __name__=='__main__':main()
