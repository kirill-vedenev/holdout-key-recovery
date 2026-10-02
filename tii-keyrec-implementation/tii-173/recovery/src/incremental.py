"""An incremental formal substitution cache, independent of the jet solver."""
import ctypes as C
from math import comb
import numpy as np
from core import HERE,U8
lib=C.CDLL(str(HERE/'libincremental.so'))
lib.series_cache_create.argtypes=[C.c_int]*5; lib.series_cache_create.restype=C.c_void_p
lib.series_cache_free.argtypes=[C.c_void_p]
lib.series_cache_set.argtypes=[C.c_void_p,C.c_int,U8]; lib.series_cache_set.restype=C.c_int
lib.series_cache_eval.argtypes=[C.c_void_p,C.c_int,U8]; lib.series_cache_eval.restype=C.c_int


class SeriesCache:
    def __init__(self,k,d,E,depth):
        self.k=k; self.N=comb(k,d); self.E=E
        modulus=sum(int(x)<<i for i,x in enumerate(E.modulus()))
        self.handle=lib.series_cache_create(k,d,E.degree(),modulus,depth)
        if not self.handle: raise ValueError('could not allocate formal series cache')
    def close(self):
        if self.handle: lib.series_cache_free(self.handle); self.handle=None
    def __del__(self): self.close()
    def set(self,r,jet):
        a=np.asarray([int(self.E(x).to_integer()) for x in jet],dtype=np.uint8)
        if len(a)!=self.k or lib.series_cache_set(self.handle,r,a): raise ValueError('invalid jet')
    def evaluate(self,r):
        a=np.zeros(self.N,dtype=np.uint8)
        if lib.series_cache_eval(self.handle,r,a): raise ValueError('invalid coefficient order')
        return a
