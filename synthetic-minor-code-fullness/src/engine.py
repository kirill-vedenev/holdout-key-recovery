"""One seeded trial: build a code, form its minor code from the tangent planes,
certify the rank, test fullness, recover a support and verify the embedding.
Failures are recorded, never replaced by new draws."""
import os
for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[key] = '1'
import gzip
import hashlib
import json
import resource
import time
from pathlib import Path
import numpy as np
from sage.all import *
import gf2m
from families import MODULI, build_instance
from minor import (all_minors_inside, grs_multiplier_fast, grs_parity, minor_code,
                   randomise_representatives, recover_public_embedding, row_budget,
                   sample_pairs, support_derivatives, tangent_defects, tangent_matrix)


def seed_for(master, label, trial, phase):
    text = '%s|%s|%s|%s' % (master, label, trial, phase)
    return int.from_bytes(hashlib.sha256(text.encode()).digest()[:16], 'big')


def encode_field(x):
    try:
        return int(x.to_integer())
    except AttributeError:
        return int(x)


def generator_sha256(rows):
    data = np.asarray(rows, dtype='<u2').tobytes()
    return hashlib.sha256(data).hexdigest()


def dump_fixture(path, obj):
    raw = json.dumps(obj, sort_keys=True, separators=(',',':')).encode()
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(gzip.compress(raw, mtime=0))
    return hashlib.sha256(raw).hexdigest()


def sage_fixture(inst, P):
    E, Y = inst['E'], inst['Y']
    obj = dict(field_order=int(E.order()), field_modulus=[int(x) for x in E.modulus()],
               parameters=P, support=[encode_field(a) for a in inst['alpha']])
    if P['family'] == 'wild':
        obj['goppa_polynomial'] = [encode_field(x) for x in inst['g']]
    else:
        obj['multiplier'] = [encode_field(a) for a in inst['lam']]
        if 'beta' in inst:
            obj['beta'] = [encode_field(a) for a in inst['beta']]
        else:
            obj['public_generator'] = [[encode_field(x) for x in row] for row in Y.rows()]
    obj['public_generator_sha256'] = generator_sha256([[encode_field(x) for x in row] for row in Y.rows()])
    return obj


def run_case(P, trial, master, output, explicit_all_pairs=False, recover=True, backend='cpp'):
    start = time.monotonic()
    result = dict(label=P['label'], family=P['family'], trial=trial, master_seed=master,
                  parameters=P, status='error', seeds={}, backend=backend)
    try:
        if P['family'] == 'binary':
            if backend != 'cpp':
                raise ValueError('binary Goppa trials run in C++; the Sage reference is in validate.py')
            run_binary(P, trial, master, output, recover, result)
        else:
            run_sage(P, trial, master, output, explicit_all_pairs, recover, backend, result)
    except Exception as exc:
        import traceback
        result.update(status='error', error_type=type(exc).__name__, error=str(exc),
                      traceback=traceback.format_exc())
    result['wall_seconds'] = time.monotonic()-start
    result['peak_rss_kib'] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return result


def set_status(result, defects, full, recover):
    if defects:
        result['status'] = 'tangent_assumption_failure'
    elif not full:
        result['status'] = 'fullness_failure'
    elif recover and not result['embedding_verified']:
        result['status'] = 'recovery_failure'
    else:
        result['status'] = 'pass'


def run_binary(P, trial, master, output, recover, result):
    """Binary Goppa trial in C++, rechecked with NumPy field arithmetic."""
    seeds = {}
    for phase in ('instance', 'representatives', 'minor_sample', 'recovery'):
        seeds[phase] = seed_for(master, P['label'], trial, phase)
        result['seeds'][phase] = str(seeds[phase])
    from native import goppa_trial
    modulus = MODULI[P['m']]
    r = goppa_trial(modulus, P['n0'], P['t'], P.get('ell', 0), P['G'],
                    [s & ((1 << 64)-1) for s in seeds.values()], recover)
    Y, keep = r['Y'], r['keep']
    alpha = r['alpha0'][keep].astype(np.int64)
    poly = [int(c) for c in r['goppa_polynomial']]
    fixture = dict(field_order=1 << P['m'], field_modulus=[(modulus >> i) & 1 for i in range(P['m']+1)],
                   parameters=P, support=[int(a) for a in alpha], goppa_polynomial=poly,
                   original_support=[int(a) for a in r['alpha0']], kept_positions=[int(j) for j in keep],
                   public_generator_sha256=generator_sha256(Y))
    result['fixture_sha256'] = dump_fixture(Path(output)/'fixtures'/P['label']/('%03d.json.gz' % trial), fixture)
    result['public_generator_sha256'] = fixture['public_generator_sha256']
    for key in ('n', 'k', 'k0', 'D', 'predicted_bound', 'parity_check_rank', 'minor_rank', 'minor_rows',
                'available_minor_rows', 'irreducibility_candidates'):
        result[key] = r[key]
    result['n0'] = P['n0']
    for key in ('public_ambient_contains', 'nominal_degree_attained', 'squarefree', 'eta_containment',
                'complete_rank', 'full_grs'):
        result[key] = bool(r[key])
    result['upper_bound_reason'] = 'binary square-free Goppa factor'
    result['tangent_defects'] = [int(j) for j in r['defects']]
    result['eligible'] = not r['tangent_defect_count']
    result['certificate'] = ('upper bound attained' if r['minor_rank'] == r['predicted_bound']
                             else 'all minor rows materialized')
    result['dimension_matches_bound'] = r['minor_rank'] == r['predicted_bound']
    result['times_seconds'] = r['times']
    assert result['squarefree'], 'Goppa polynomial is not square-free'
    assert result['public_ambient_contains'], 'public ambient containment failed'

    # Independent recheck of the ambient containment and of the recovered embedding.
    t0 = time.monotonic()
    galpha = gf2m.evaluate(modulus, poly, alpha)
    assert np.all(galpha != 0) and len(set(alpha.tolist())) == len(alpha)
    lam = gf2m.multiply(modulus, gf2m.multiply(modulus, galpha, galpha),
                        gf2m.inverse(modulus, gf2m.support_derivatives(modulus, alpha)))
    assert gf2m.binary_in_grs(modulus, Y, alpha, r['D']+1, lam), 'NumPy ambient check failed'
    if r['recovery_code'] == 2:
        ar, lm = r['recovered_support'].astype(np.int64), r['recovered_multiplier'].astype(np.int64)
        assert gf2m.binary_in_grs(modulus, Y, ar, r['D']+1, lm), 'independent embedding verification failed'
        result['recovered_support'] = [int(x) for x in ar]
        result['recovered_multiplier'] = [int(x) for x in lm]
    result['times_seconds']['numpy_checks'] = time.monotonic()-t0
    if recover and r['full_grs']:
        result['support_recovered'] = r['recovery_code'] >= 1
        result['embedding_verified'] = r['recovery_code'] == 2
    else:
        result.update(support_recovered=None, embedding_verified=None)
    set_status(result, r['tangent_defect_count'], r['full_grs'], recover)


def run_sage(P, trial, master, output, explicit_all_pairs, recover, backend, result):
    """GRS, alternant and wild Goppa trials: instances in Sage, minor code in C++ or Sage."""
    times = {}
    compiled = None
    def phase(name):
        s = seed_for(master, P['label'], trial, name)
        result['seeds'][name] = str(s)
        set_random_seed(s)
        return s
    def stamp(name, t):
        times[name] = time.monotonic()-t
    try:
        t0=time.monotonic(); s=phase('instance')
        inst = build_instance(P,s)
        E,Y,alpha,lam,D=inst['E'],inst['Y'],inst['alpha'],inst['lam'],inst['D']
        n,k=Y.ncols(),Y.nrows()
        result.update(n=int(n), k=int(k), D=int(D))
        obj=sage_fixture(inst,P)
        result['fixture_sha256']=dump_fixture(Path(output)/'fixtures'/P['label']/('%03d.json.gz'%trial),obj)
        result['public_generator_sha256']=obj['public_generator_sha256']
        stamp('build',t0)

        # Membership of every minor row in the ambient code below is a theorem.
        t0=time.monotonic()
        H=grs_parity(E,alpha,D+1,lam)
        ambient=(Y*H.transpose()).is_zero()
        result['public_ambient_contains']=bool(ambient)
        assert ambient, 'public ambient containment failed'
        nominal_degree_attained=not (Y*grs_parity(E,alpha,D,lam).transpose()).is_zero()
        result['nominal_degree_attained']=bool(nominal_degree_attained)
        bound=D
        factor=[E.one()]*n
        if P['family']=='wild' and P['m']==2 and P['g']=='irred' and P['t']>1:
            g=inst['g']; factor=[g(x) for x in alpha]
            bound=D-P['t']
            stronger=vector(E,[lam[i]*factor[i] for i in range(n)])
            member=(Y*grs_parity(E,alpha,bound+1,stronger).transpose()).is_zero()
            assert member, 'membership in GRS_{D-t+1}(alpha, lambda g(alpha)) failed'
            result['quadratic_identity_membership']=bool(member)
            result['upper_bound_reason']='m = 2, g without roots in F_{q^2}'
        else:
            result['upper_bound_reason']='characteristic two'
        result['predicted_bound']=int(bound)
        stamp('ambient',t0)

        t0=time.monotonic(); s=phase('representatives')
        Yhat,a=randomise_representatives(E,Y,tangent_matrix(E,Y,alpha,inst['nu']),s)
        defects=tangent_defects(Y,Yhat)
        result['tangent_defects']=[int(i) for i in defects]
        result['eligible']=not defects
        dpi=support_derivatives(E,alpha)
        assert all(lam[i]*inst['nu'][i]*dpi[i]==1 for i in range(n))
        eta=vector(E,[(lam[i]*a[i]/dpi[i])**(E.order()//2)*factor[i] for i in range(n)])
        He=grs_parity(E,alpha,bound,eta)
        stamp('tangents',t0)

        t0=time.monotonic(); s=phase('minor_sample')
        available=int(binomial(k,2)); budget=int(row_budget(k,bound))
        last_rank=-1; complete=False
        while True:
            pairs=sample_pairs(k,budget,s)
            if backend=='cpp':
                from native import Minor
                if compiled is not None:
                    compiled.close()
                compiled=Minor(E,Y,Yhat,pairs)
                rank=compiled.rank
                assert compiled.verify(alpha,eta,bound), 'sample violates the predicted containment'
                rows=compiled.rows
            else:
                T=minor_code(Y,Yhat,pairs)
                rank=int(T.rank());rows=T.nrows()
                assert (T*He.transpose()).is_zero(), 'sample violates the predicted containment'
            if rank==bound:
                # Every minor row lies in the bound's GRS code; the sample spans it.
                complete=True
                method='upper bound attained'
                if explicit_all_pairs:
                    complete=all_minors_inside(Y,Yhat,He)
                    assert complete, 'explicit all-pairs certificate failed'
                    method='explicit all-pairs parity certificate'
                break
            if backend=='cpp':
                T=minor_code(Y,Yhat,pairs)
                assert T.rank()==rank, 'C++/Sage rank mismatch'
            if budget==available:
                complete=True; method='all minor rows materialized'; break
            if rank==last_rank:
                Hspan=T.right_kernel_matrix()
                if all_minors_inside(Y,Yhat,Hspan):
                    complete=True; method='explicit all-pairs span certificate'; break
            last_rank=rank; budget=min(available,2*budget)
        result.update(minor_rank=rank,minor_rows=int(rows),available_minor_rows=available,
                      complete_rank=bool(complete),certificate=method,dimension_matches_bound=rank==bound)
        if compiled is not None:
            result['native_minor_seconds']=compiled.times
        phase('fullness')
        if rank==bound:
            full=True
        elif rank>0 and not defects:
            full=grs_multiplier_fast(E,T,alpha,rank) is not None
        else:
            full=False
        result['full_grs']=bool(full)
        stamp('minor_and_certificate',t0)

        t0=time.monotonic(); phase('recovery')
        if recover and full:
            if compiled is not None:
                recovered=compiled.recover(D,int(result['seeds']['recovery']))
                ar=recovered.pop('recovered_support',None)
                lm=recovered.pop('recovered_multiplier',None)
                if ar is not None:
                    assert (Y*grs_parity(E,ar,D+1,lm).transpose()).is_zero(), 'independent embedding verification failed'
                    recovered['recovered_support']=[encode_field(x) for x in ar]
                    recovered['recovered_multiplier']=[encode_field(x) for x in lm]
                result.update(recovered)
            else:
                result.update(recover_public_embedding(E,Y,T,D))
        else:
            result.update(support_recovered=None,embedding_verified=None)
        stamp('recovery',t0)
        set_status(result, defects, full, recover)
    finally:
        result['times_seconds']=times
        if compiled is not None:
            compiled.close()
