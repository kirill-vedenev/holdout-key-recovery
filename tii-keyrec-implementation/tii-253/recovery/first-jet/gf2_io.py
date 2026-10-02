"""Exact small packed binary and extension-field linear algebra."""
from sage.all import GF,matrix
from packed_gf2 import packed_matrix
import numpy as np


def binary_times_field(B,X):
    E=X.base_ring(); is_vector=not hasattr(X,'ncols')
    M=matrix(E,len(X),1,list(X)) if is_vector else X
    if E==GF(2): result=B*M
    else:
        values=np.fromiter((int(x.to_integer()) for x in M.list()),dtype=np.uint8).reshape(M.nrows(),M.ncols())
        rows=[]
        for j in range(M.ncols()):
            for bit in range(E.degree()):
                rows.append(int.from_bytes(np.packbits((values[:,j]>>bit)&1,bitorder='little').tobytes(),'little'))
        R=B*packed_matrix(rows,M.nrows()).transpose()
        elements=[E.from_integer(a) for a in range(E.order())]
        result=matrix(E,B.nrows(),M.ncols(),lambda a,j:elements[sum(int(R[a,j*E.degree()+bit])<<bit for bit in range(E.degree()))])
    return result.column(0) if is_vector else result


def panel_matrix(path,width,words):
    raw=np.fromfile(path,dtype='<u8').reshape(-1,words)
    return packed_matrix([int.from_bytes(row.tobytes(),'little') for row in raw],width).transpose()
