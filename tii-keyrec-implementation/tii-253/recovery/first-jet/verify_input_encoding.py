"""Check every coefficient of the supplied panel against the verified THK1 file."""
from pathlib import Path
from math import comb
import argparse,hashlib,json,struct,time
import numpy as np


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('kernel',type=Path);ap.add_argument('panel',type=Path);ap.add_argument('report',type=Path)
    args=ap.parse_args();begin=time.monotonic()
    with args.kernel.open('rb') as f:header=f.read(20)
    assert header[:4]==b'THK1';k,d,N,count=struct.unpack('<4I',header[4:]);assert N==comb(k,d)
    stride=(N+7)//8;assert args.kernel.stat().st_size==20+stride*count
    panel_words=args.panel.stat().st_size//(N*8);assert args.panel.stat().st_size==N*panel_words*8
    packed=np.memmap(args.kernel,dtype=np.uint8,mode='r',offset=20,shape=(count,stride))
    panel=np.memmap(args.panel,dtype='<u8',mode='r',shape=(N,panel_words))
    for start in range(0,N,65536):
        end=min(start+65536,N)
        original=np.unpackbits(np.ascontiguousarray(packed[:,start//8:(end+7)//8]),axis=1,bitorder='little')[:,:end-start].transpose()
        recovered=np.unpackbits(np.ascontiguousarray(panel[start:end]).view(np.uint8).reshape(end-start,panel_words*8),axis=1,bitorder='little')
        assert np.array_equal(original,recovered[:,:count]),('panel/THK1 mismatch',start)
        assert not recovered[:,count:].any()
    def digest(p):
        h=hashlib.sha256()
        with p.open('rb') as f:
            for b in iter(lambda:f.read(1<<20),b''):h.update(b)
        return h.hexdigest()
    report=dict(status='passed',verified_thk1_equals_input_panel=True,every_coefficient_checked=True,unused_panel_columns_zero=True,
        k=k,degree=d,monomials=N,polynomials=count,exact_binary_coefficients_checked=N*count,
        kernel_sha256=digest(args.kernel),panel_sha256=digest(args.panel),seconds=time.monotonic()-begin)
    args.report.write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report),flush=True)

if __name__=='__main__':main()
