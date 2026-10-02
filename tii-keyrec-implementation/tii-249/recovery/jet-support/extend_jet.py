"""Exact binary-normalized continuation from the public first-jet branch."""
from sage.all import *
from pathlib import Path
from itertools import combinations
from math import comb
import argparse,json,time,ctypes as C,os
import numpy as np
from gf2_io import panel_matrix
from public_local import polar_matrix


class Cache:
    def __init__(self,k,d,E,depth,count,threads,panel):
        self.count=count;self.lib=C.CDLL(os.environ.get('TII_JET_SERIES_LIBRARY',str(Path(__file__).with_name('libseries.so'))))
        self.U8=np.ctypeslib.ndpointer(dtype=np.uint8,flags='C_CONTIGUOUS')
        self.lib.cache_create.argtypes=[C.c_int]*7+[C.c_char_p];self.lib.cache_create.restype=C.c_void_p
        self.lib.cache_free.argtypes=[C.c_void_p]
        self.lib.cache_set.argtypes=[C.c_void_p,C.c_int,self.U8]
        self.lib.cache_eval.argtypes=[C.c_void_p,C.c_int,self.U8]
        self.lib.cache_pair.argtypes=[C.c_void_p,C.c_int,self.U8]
        self.lib.cache_patch.argtypes=[C.c_void_p,C.c_int,self.U8,self.U8,self.U8]
        modulus=sum(int(x)<<i for i,x in enumerate(E.modulus()))
        self.handle=self.lib.cache_create(k,d,E.degree(),modulus,depth,count,threads,str(panel.resolve()).encode())
        if not self.handle:raise ValueError('formal substitution cache allocation failed')
        self.elements=[E.from_integer(a) for a in range(E.order())]
    def close(self):self.lib.cache_free(self.handle);self.handle=None
    def array(self,v):return np.asarray([int(x.to_integer()) for x in v],dtype=np.uint8)
    def set(self,r,v):self.lib.cache_set(self.handle,r,self.array(v))
    def evaluate(self,r):
        out=np.zeros(self.count,dtype=np.uint8);self.lib.cache_eval(self.handle,r,out)
        return vector(self.elements[0].parent(),[self.elements[c] for c in out])
    def pair(self,r):
        out=np.zeros(2*self.count,dtype=np.uint8);self.lib.cache_pair(self.handle,r,out)
        return [vector(self.elements[0].parent(),[self.elements[c] for c in row]) for row in out.reshape(2,-1)]
    def patch(self,r,odd,even):
        out=np.zeros(2*self.count,dtype=np.uint8);self.lib.cache_patch(self.handle,r,self.array(odd),self.array(even),out)
        return [vector(self.elements[0].parent(),[self.elements[c] for c in row]) for row in out.reshape(2,-1)]


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('public',type=Path);ap.add_argument('panel',type=Path);ap.add_argument('result',type=Path)
    ap.add_argument('--depth',type=int,default=290);ap.add_argument('--threads',type=int,default=4)
    args=ap.parse_args();begin=time.monotonic();data=json.loads(args.public.read_text())
    saved=json.loads((args.result/'local-branches.json').read_text());k=data['parameters']['k'];d=data['parameters']['degree'];count=saved['kernel_polynomials']
    E=GF(2**data['parameters']['m'],'a',modulus=PolynomialRing(GF(2),'x')(data['field_modulus']))
    p=vector(E,[(data['heldout_mask']>>a)&1 for a in range(k)]);v1=vector(E,[E.from_integer(x) for x in saved['branch_vectors'][0]])
    G=panel_matrix(args.result/'hasse-1.u64le',count,5);Q=panel_matrix(args.result/'hasse-2.u64le',count,5)
    pairs=sorted(combinations(range(k),2),key=lambda pair:sum(comb(a,j+1) for j,a in enumerate(pair)))
    J=G.change_ring(E);L=polar_matrix(Q,v1,pairs);A=J.stack(L)
    assert A.right_kernel()==matrix(E,[p,v1]).row_space()
    columns=A.pivots();independent=A.matrix_from_columns(columns);rows=independent.transpose().pivots();inverse=independent.matrix_from_rows(rows).inverse()
    def solve(rhs):
        coefficients=inverse*vector(E,[rhs[r] for r in rows]);answer=vector(E,k)
        for j,c in zip(columns,coefficients):answer[j]=c
        assert A*answer==rhs
        return answer
    jets=[p,v1,vector(E,[x*x for x in v1])];history=[];cache=Cache(k,d,E,args.depth,count,args.threads,args.panel)
    first_full=None
    try:
        for r,v in enumerate(jets):cache.set(r,v)
        for r in range(3):assert not cache.evaluate(r)
        for r in range(3,args.depth,2):
            even=vector(E,[x*x for x in jets[(r+1)//2]])
            c0,c1=cache.pair(r);odd=solve(-vector(E,list(c0)+list(c1+J*even)))
            actual=cache.patch(r,odd,even);assert all(not v for v in actual)
            # Independently recompute both complete convolutions at a few orders;
            # this checks the incremental correction identity against direct sums.
            if r in (3,17,args.depth-1):assert all(not v for v in cache.pair(r))
            jets.extend([odd,even]);flag=matrix(E,jets).rank()
            if first_full is None and flag==k:first_full=r+1
            history.append(dict(odd_order=r,even_order=r+1,flag_dimension=flag,exact_residual_zero=True))
            if (r+1)%16==0:
                print(json.dumps(dict(event='jet_extended',order=r+1,flag_dimension=flag,seconds=time.monotonic()-begin)),flush=True)
    finally:cache.close()
    assert len(jets)==args.depth+1
    assert all(jets[2*r]==vector(E,[x*x for x in jets[r]]) for r in range(1,args.depth//2+1))
    report=dict(public_only=True,position=data['heldout'][0],position_original=data['heldout_original'],depth=args.depth,
        binary_normalization=True,even_coefficients_fixed_by_squaring=True,all_next_even_order_equations_used=True,
        monomial_order='colex',kernel_dimension=count,branches=saved['branches'],branch_index=0,
        fixed_continuation_rank=A.rank(),affine_nullity=k-A.rank(),ambiguity_equals_first_plane=True,
        all_recurrence_equations_verified=True,all_even_square_relations_verified=True,
        native_direct_recomputation_orders=[0,1,2,3,4,17,18,args.depth-1,args.depth],first_full_flag_even_order=first_full,
        jets=[[int(x.to_integer()) for x in v] for v in jets],history=history,
        field_modulus=data['field_modulus'],seconds=time.monotonic()-begin)
    (args.result/f'deep-jets-{data["heldout"][0]}.json').write_text(json.dumps(report,indent=2,default=int)+'\n')
    print(json.dumps({a:b for a,b in report.items() if a not in ['jets','history']},default=int),flush=True)

if __name__=='__main__':main()
