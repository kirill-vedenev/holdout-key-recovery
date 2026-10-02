#!/usr/bin/env python3
"""Public TII-253 dimensions and a short-run operator request; no full solve."""
import hashlib
import importlib.util
import json
import math
from pathlib import Path

import argparse
SOURCE_DIR = Path(__file__).resolve().parent
PUBLIC_DATA = SOURCE_DIR.parents[3] / 'tii-results' / 'tii-253'
reference = SOURCE_DIR / 'prepare_reference.py'
module_spec = importlib.util.spec_from_file_location('reference', reference)
reference_module = importlib.util.module_from_spec(module_spec)
module_spec.loader.exec_module(reference_module)
rref = reference_module.rref
nullspace = reference_module.nullspace
restrict = reference_module.restrict


def choose(n, r):
    return math.comb(n, r) if 0 <= r <= n else 0


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output', type=Path, required=True)
    args = ap.parse_args()
    HERE = args.output.resolve()
    HERE.mkdir(parents=True, exist_ok=False)
    source = json.loads((PUBLIC_DATA / 'source.json').read_text())
    raw = (PUBLIC_DATA / 'pk_McEliece_253.txt').read_bytes()
    assert hashlib.sha256(raw).hexdigest() == source['sha256']
    lines = raw.decode().strip().splitlines()
    assert len(lines) == 73
    rows = []
    for line in lines[:-1]:
        bits = line.strip()[1:-1].split()
        assert len(bits) == 214 and set(bits) <= {'0', '1'}
        rows.append(sum(int(bit) << j for j, bit in enumerate(bits)))
    assert len(rref(rows, 214)[1]) == 72
    k, degree, multiplicity, heldout = 46, 6, 5, 0
    _, free = nullspace(rows, 214)
    removed = [j for j in free if j != heldout][:142-k]
    keep = [j for j in range(214) if j not in removed]
    hs = [restrict(row, keep) for row in rows]
    generator, free = nullspace(hs, len(keep))
    assert len(generator) == k and len(keep) == 118
    held_local = keep.index(heldout)
    order = [j for j in free if j != held_local]
    order += [j for j in range(len(keep)) if j not in order and j != held_local]
    _, pivots = rref([restrict(row, order) for row in generator], len(order))
    info = [order[j] for j in pivots]
    augmented, pivots = rref([restrict(row, info) | (row << k) for row in generator], k+len(keep))
    assert pivots == list(range(k))
    generator = [row >> k for row in augmented]
    for i, row in enumerate(generator):
        assert restrict(row, info) == 1 << i
        assert all((row & h).bit_count() % 2 == 0 for h in hs)
    columns = [sum(((row >> j) & 1) << i for i, row in enumerate(generator)) for j in range(len(keep))]
    infinity = 0
    for point in columns:
        infinity ^= point
    assert columns[held_local].bit_count() >= degree
    points = [dict(original=keep[j], mask=point, multiplicity=multiplicity)
              for j, point in enumerate(columns) if j != held_local]
    points.append(dict(original=-1, mask=infinity, multiplicity=multiplicity))
    active = [p for p in points if p['mask'].bit_count() >= degree-multiplicity+1]
    levels = reference_module.retained_orders(degree, multiplicity)
    assert levels == [2, 4]
    nmon = math.comb(k, degree)
    nrows = len(active) * sum(math.comb(k, u) for u in levels)
    nnz = sum(choose(p['mask'].bit_count(), a)*choose(k-p['mask'].bit_count(), u-a)
              *choose(p['mask'].bit_count()-a, degree-u)
              for p in active for u in levels for a in range(u+1))
    predicted, coefficients = reference_module.kappa(k, degree, 72)
    public = dict(challenge='TII-253', parameters=dict(m=8,t=9,n=118,k=46,degree=6),
                  original_parameters=dict(m=8,t=9,n=214,k=142),
                  public_sha256=source['sha256'], field_modulus=json.loads(lines[-1]),
                  original_indices=keep,shortened_original_indices=removed,
                  information_set=info,information_set_original=[keep[j] for j in info],
                  heldout=[held_local],heldout_original=heldout,heldout_mask=columns[held_local],
                  minimum_elements=280,default_elements=512,points=points,
                  infinity_mask=infinity,infinity_multiplicity=multiplicity,
                  multiplicities=[0 if j==held_local else multiplicity for j in range(len(keep))],
                  generator=[[(row >> j)&1 for j in range(len(keep))] for row in generator])
    public['instance_id'] = hashlib.sha256(reference_module.canonical(public)).hexdigest()
    (HERE/'public.json').write_text(json.dumps(public,indent=2)+'\n')
    op = ['TII_HOLDOUT_V1',f'instance_id {public["instance_id"]}',f'k {k}',f'degree {degree}',
          'minimum 280','heldout_original 0',f'heldout_mask {columns[held_local]}',f'points {len(points)}']
    op += [f'{p["original"]} {p["mask"]} {multiplicity}' for p in points]
    (HERE/'operator.txt').write_text('\n'.join(op)+'\n')
    request = ['mceliecex-holdout-cuda-request-v1',f'public_sha256 {source["sha256"]}',
               f'k {k}',f'degree {degree}','words 8','repetitions 1','seed 25320260930',
               'emit_output 0',f'expected_columns {nmon}',f'expected_rows {nrows}',
               f'expected_nonzeros {nnz}',f'point_count {len(active)}']
    request += [f'point {p["mask"]} 2 2 4' for p in active]
    (HERE/'request.txt').write_text('\n'.join(request)+'\n')
    report = dict(status='public-input-preflight',production_started=False,
                  parameters=public['parameters'],original_parameters=public['original_parameters'],
                  heldout_original=0,heldout_weight=columns[held_local].bit_count(),
                  monomials=nmon,reduced_rows=nrows,nonzeros=nnz,active_points=len(active),
                  retained_orders=levels,minimum_elements=280,block_bits=512,
                  predicted_kernel_dimension=predicted,prediction_verified=False,
                  kappa_coefficients=coefficients,csr_bytes=(nrows+1)*8+nnz*4,
                  panel_bytes=nmon*64,fft_length=1<<(nrows-2).bit_length(),
                  planned_sequence_terms=2*math.ceil(nmon/512)+256,
                  composition_degree_bound=degree*(len(keep)-18-1),
                  forced_zero_count=multiplicity*len(keep),
                  weis_count_shortfall=degree*(len(keep)-18-1)+1-multiplicity*len(keep))
    (HERE/'dimensions.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__ == '__main__':
    main()
