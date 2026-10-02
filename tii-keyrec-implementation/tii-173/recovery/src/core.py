"""Packed binary kernel IO and fast formal substitution; no hidden code data."""
from pathlib import Path
from math import comb
import ctypes as C
import struct, sys
import numpy as np
from sage.all import GF, matrix
HERE=Path(__file__).resolve().parent
from packed_gf2 import packed_matrix
lib=C.CDLL(str(HERE/'libseries.so'))
U8=np.ctypeslib.ndpointer(dtype=np.uint8,flags='C_CONTIGUOUS')
lib.series_products.argtypes=[C.c_int]*5+[U8,U8]
lib.series_products.restype=C.c_int
lib.local_rows.argtypes=[C.c_int,C.c_int,C.c_uint64,U8]
lib.local_rows.restype=C.c_int


def read_kernel(path):
    with open(path,'rb') as f:
        if f.read(4)!=b'WHK1': raise ValueError('invalid kernel format')
        k,d,N,dim=struct.unpack('<4I',f.read(16))
        if N!=comb(k,d): raise ValueError('invalid kernel monomial count')
        width=(N+7)//8
        rows=[]
        for _ in range(dim):
            row=f.read(width)
            if len(row)!=width: raise ValueError('truncated kernel')
            rows.append(int.from_bytes(row,'little'))
        if f.read(1): raise ValueError('trailing kernel data')
    return k,d,packed_matrix(rows,N)


def jets_at_point(k,d,point):
    N=comb(k,d); width=(N+7)//8
    data=np.zeros((1+k+comb(k,2),width),dtype=np.uint8)
    mask=sum(int(x)<<a for a,x in enumerate(point))
    if lib.local_rows(k,d,mask,data): raise ValueError('invalid local jet request')
    return packed_matrix([int.from_bytes(row.tobytes(),'little') for row in data],N)


def products(arcs,d,field,depth):
    k=len(arcs); m=field.degree()
    if not 1<=m<=8: raise ValueError('formal substitution supports extension degrees 1..8')
    modulus=sum(int(x)<<i for i,x in enumerate(field.modulus()))
    data=np.zeros((k,depth+1),dtype=np.uint8)
    for a,arc in enumerate(arcs):
        for r,c in enumerate(arc[:depth+1]): data[a,r]=int(field(c).to_integer())
    result=np.zeros((comb(k,d),depth+1),dtype=np.uint8)
    if lib.series_products(k,d,m,modulus,depth,data,result):
        raise ValueError('native substitution supports fields of degree 1..8')
    return result


def apply_binary(K,coefficients,field):
    """Apply binary polynomial coefficients to field-valued functionals."""
    n,r=coefficients.shape
    assert K.ncols()==n
    m=field.degree(); rows=[]
    for j in range(r):
        for bit in range(m):
            packed=np.packbits((coefficients[:,j]>>bit)&1,bitorder='little')
            rows.append(int.from_bytes(packed.tobytes(),'little'))
    B=K*packed_matrix(rows,n).transpose()
    elements=[field.from_integer(a) for a in range(field.order())]
    return matrix(field,K.nrows(),r,
                  lambda a,j:elements[sum(int(B[a,j*m+bit])<<bit for bit in range(m))])


def binary_times_field(B, X):
    """Multiply without converting a potentially huge binary matrix to GF(2^m)."""
    E=X.base_ring(); is_vector=not hasattr(X,'ncols')
    if E.degree()>8: raise ValueError('packed extension-field products support degrees at most 8')
    M=matrix(E,len(X),1,list(X)) if is_vector else X
    if E==GF(2): result=B*M
    else:
        values=np.fromiter((int(x.to_integer()) for x in M.list()),dtype=np.uint8).reshape(M.nrows(),M.ncols())
        result=apply_binary(B,values,E)
    return result.column(0) if is_vector else result
