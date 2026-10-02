"""Array interface to the C++/M4RIE backend (libminor_native.so)."""
import ctypes as C
from pathlib import Path
import numpy as np

LIBRARY=Path(__file__).with_name('libminor_native.so')
lib=C.CDLL(str(LIBRARY))
U8=np.ctypeslib.ndpointer(dtype=np.uint8,flags='C_CONTIGUOUS')
U16=np.ctypeslib.ndpointer(dtype=np.uint16,flags='C_CONTIGUOUS')
I32=np.ctypeslib.ndpointer(dtype=np.int32,flags='C_CONTIGUOUS')
U64=np.ctypeslib.ndpointer(dtype=np.uint64,flags='C_CONTIGUOUS')
F64=np.ctypeslib.ndpointer(dtype=np.float64,flags='C_CONTIGUOUS')
lib.minor_error.restype=C.c_char_p
lib.minor_create.argtypes=[C.c_uint,C.c_int,C.c_int,U16,U16,C.c_int,I32]
lib.minor_create.restype=C.c_void_p
lib.minor_free.argtypes=[C.c_void_p]
lib.minor_rank.argtypes=[C.c_void_p]
lib.minor_times.argtypes=[C.c_void_p,F64]
lib.minor_verify.argtypes=[C.c_void_p,U16,U16,C.c_int]
lib.minor_recover.argtypes=[C.c_void_p,C.c_int,C.c_uint64,U16,U16,F64]
lib.goppa_trial.argtypes=[C.c_uint,C.c_int,C.c_int,C.c_int,C.c_int,U64,C.c_int,
                          I32,F64,U16,I32,U16,U8,I32,U16,U16,C.c_void_p]

GOPPA_KINDS={'irred':0,'split':1,'separable':2}
GOPPA_OUTPUTS=('n','k','k0','D','predicted_bound','parity_check_rank','public_ambient_contains',
               'nominal_degree_attained','squarefree','tangent_defect_count','minor_rank',
               'minor_rows','available_minor_rows','complete_rank','eta_containment',
               'full_grs','recovery_code','irreducibility_candidates')
GOPPA_TIMES=('key','code','ambient','tangents','minor_construction','minor_rref',
             'certificate','support','minor_multiplier','public_embedding')


def encoded(values):
    return np.fromiter((int(x.to_integer()) for x in values),dtype=np.uint16)


def encoded_matrix(matrix):
    return encoded(x for row in matrix.rows() for x in row)


def error():
    return RuntimeError(lib.minor_error().decode())


class Minor:
    """Minor matrix of a Sage public matrix Y and tangent representatives Yhat."""
    def __init__(self,E,Y,Yhat,pairs):
        self.handle=None
        self.E=E;self.n=Y.ncols()
        modulus=sum(int(c)<<i for i,c in enumerate(E.modulus()))
        if pairs is None:
            pairs=[(a,b) for a in range(Y.nrows()) for b in range(a+1,Y.nrows())]
        self.rows=len(pairs)
        self.handle=lib.minor_create(modulus,Y.nrows(),self.n,encoded_matrix(Y),encoded_matrix(Yhat),
                                     self.rows,np.asarray(pairs,dtype=np.int32))
        if not self.handle:
            raise error()
        self.rank=lib.minor_rank(self.handle)
        timing=np.zeros(2)
        lib.minor_times(self.handle,timing)
        self.times=dict(construction=float(timing[0]),rref=float(timing[1]))

    def close(self):
        if self.handle:
            lib.minor_free(self.handle);self.handle=None

    def __del__(self):
        if lib is not None:
            self.close()

    def verify(self,support,multiplier,dimension):
        value=lib.minor_verify(self.handle,encoded(support),encoded(multiplier),dimension)
        if value<0:
            raise error()
        return bool(value)

    def recover(self,D,seed):
        # Only the public state built in __init__ and D are used.
        support=np.zeros(self.n,dtype=np.uint16);multiplier=np.zeros_like(support)
        timing=np.zeros(3)
        value=lib.minor_recover(self.handle,D,seed&((1<<64)-1),support,multiplier,timing)
        if value<0:
            raise error()
        out=dict(support_recovered=value>=1,embedding_verified=value==2,
                 native_recovery_seconds=dict(zip(('support','minor_multiplier','public_embedding'),map(float,timing))))
        if value==2:
            out['recovered_support']=[self.E.from_integer(int(x)) for x in support]
            out['recovered_multiplier']=[self.E.from_integer(int(x)) for x in multiplier]
        return out


def goppa_trial(modulus,n0,t,ell,kind,seeds,recover=True,tangents=False):
    """Run one binary Goppa trial in C++; seeds are four 64-bit integers for the
    key, the tangent representatives, the minor sample and the recovery."""
    n=n0-ell
    out=np.zeros(len(GOPPA_OUTPUTS),dtype=np.int32);times=np.zeros(len(GOPPA_TIMES))
    alpha0=np.zeros(n0,dtype=np.uint16);keep=np.zeros(n0,dtype=np.int32)
    poly=np.zeros(t+1,dtype=np.uint16);Y=np.zeros(n*n,dtype=np.uint8)
    defects=np.zeros(n,dtype=np.int32)
    support=np.zeros(n,dtype=np.uint16);multiplier=np.zeros(n,dtype=np.uint16)
    Yhat=np.zeros(n*n,dtype=np.uint16) if tangents else None
    value=lib.goppa_trial(modulus,n0,t,ell,GOPPA_KINDS[kind],np.asarray(seeds,dtype=np.uint64),int(recover),
                          out,times,alpha0,keep,poly,Y,defects,support,multiplier,
                          None if Yhat is None else Yhat.ctypes.data)
    if value<0:
        raise error()
    r=dict(zip(GOPPA_OUTPUTS,map(int,out)))
    k=r['k']
    r.update(times=dict(zip(GOPPA_TIMES,map(float,times))),alpha0=alpha0,keep=keep[:n],
             goppa_polynomial=poly,Y=Y[:k*n].reshape(k,n),defects=defects[:r['tangent_defect_count']])
    if Yhat is not None:
        r['Yhat']=Yhat[:k*n].reshape(k,n)
    if r['recovery_code']==2:
        r['recovered_support']=support;r['recovered_multiplier']=multiplier
    return r
