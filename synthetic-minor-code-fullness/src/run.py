#!/usr/bin/env python3
"""Run or resume seeded trials in parallel; save every outcome atomically."""
import argparse
import concurrent.futures
import hashlib
import json
import multiprocessing
import os
from pathlib import Path
import platform
import signal
import time

for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[name]='1'
SRC=Path(__file__).resolve().parent


def atomic_json(path, value):
    path=Path(path); path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_suffix(path.suffix+'.tmp')
    tmp.write_text(json.dumps(value,sort_keys=True,indent=2)+'\n')
    os.replace(tmp,path)


def source_hashes():
    paths=sorted(p for p in SRC.iterdir() if p.suffix in ('.py','.cpp'))
    return {'src/'+p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}


def work(args):
    from engine import run_case
    P,trial,master,output,explicit,recover,timeout,backend=args
    def alarm(_sig,_frame):
        raise TimeoutError('per-instance timeout after %d seconds'%timeout)
    signal.signal(signal.SIGALRM,alarm)
    signal.alarm(timeout)
    try:
        result=run_case(P,trial,master,output,explicit,recover,backend)
        if result.get('error_type')=='TimeoutError':
            result['status']='timeout'
        atomic_json(Path(output)/'cases'/P['label']/('%03d.json'%trial),result)
        return result
    finally:
        signal.alarm(0)


def size(P):
    return P.get('n0',P.get('n',0))-P.get('ell',0)


def main():
    from families import all_parameters
    ap=argparse.ArgumentParser()
    ap.add_argument('--output',required=True)
    ap.add_argument('--trials',type=int,default=100)
    ap.add_argument('--mceliece-trials',type=int,default=20,help='trials for the Classic McEliece sets')
    ap.add_argument('--trial-start',type=int,default=0)
    ap.add_argument('--workers',type=int,default=max(1,(os.cpu_count() or 2)-1))
    ap.add_argument('--master-seed',default='minor-code-fullness')
    ap.add_argument('--labels',help='comma-separated labels; default all sets')
    ap.add_argument('--explicit-all-pairs',action='store_true')
    ap.add_argument('--no-recovery',action='store_true')
    ap.add_argument('--backend',choices=('cpp','sage'),default='cpp')
    ap.add_argument('--timeout',type=int,default=3600)
    ap.add_argument('--retry-status',default='',help='comma-separated saved statuses to rerun with the SAME seed')
    a=ap.parse_args()
    if a.trials<1 or a.mceliece_trials<1 or a.workers<1 or a.trial_start<0:
        ap.error('trials and workers must be positive')
    params=all_parameters()
    if a.labels:
        labels=set(a.labels.split(',')); params=[P for P in params if P['label'] in labels]
        if {P['label'] for P in params}!=labels:
            ap.error('unknown labels: '+str(labels-{P['label'] for P in params}))
    out=Path(a.output).resolve(); out.mkdir(parents=True,exist_ok=True)
    import engine
    native_hash=None
    if a.backend=='cpp':
        from native import LIBRARY
        native_hash=hashlib.sha256(LIBRARY.read_bytes()).hexdigest()
    from sage.env import SAGE_VERSION
    trials={P['label']:a.mceliece_trials if P['label'].startswith('mceliece') else a.trials for P in params}
    spec=dict(master_seed=a.master_seed,trials=trials,trial_start=a.trial_start,parameters=params,
              explicit_all_pairs=a.explicit_all_pairs,recovery=not a.no_recovery,
              sources=source_hashes(),sage_version=SAGE_VERSION,backend=a.backend)
    fingerprint=hashlib.sha256(json.dumps(spec,sort_keys=True).encode()).hexdigest()
    path=out/'manifest.json'
    if path.exists():
        old=json.loads(path.read_text())
        if old['spec_sha256']!=fingerprint:
            raise SystemExit('Refusing to mix different sources/settings: choose another output directory.')
    else:
        atomic_json(path,dict(spec=spec,spec_sha256=fingerprint,platform=platform.platform(),
                             python=platform.python_version(),cpus=os.cpu_count(),workers=a.workers,
                             timeout_seconds=a.timeout,native_binary_sha256=native_hash))
    retry=set(a.retry_status.split(','))-{''}
    tasks=[]; saved=0
    # Start larger jobs first to avoid a long serial tail.
    ordered=sorted(params,key=size,reverse=True)
    for trial in range(a.trial_start,a.trial_start+max(trials.values())):
        for P in ordered:
            if trial>=a.trial_start+trials[P['label']]:
                continue
            dest=out/'cases'/P['label']/('%03d.json'%trial)
            if dest.exists() and json.loads(dest.read_text())['status'] not in retry:
                saved+=1; continue
            tasks.append((P,trial,a.master_seed,str(out),a.explicit_all_pairs,not a.no_recovery,a.timeout,a.backend))
    tasks.sort(key=lambda task:size(task[0]),reverse=True)
    total=sum(trials.values())
    print(json.dumps(dict(event='start',total=total,saved=saved,pending=len(tasks),workers=a.workers)),flush=True)
    started=time.monotonic(); counts={}; done=saved
    # Fork after importing Sage, sharing read-only runtime pages.
    context=multiprocessing.get_context('fork')
    with concurrent.futures.ProcessPoolExecutor(max_workers=a.workers,mp_context=context) as pool:
        futures={pool.submit(work,task):(task[0]['label'],task[1]) for task in tasks}
        for future in concurrent.futures.as_completed(futures):
            label,trial=futures[future]
            result=future.result()
            done+=1; status=result['status']; counts[status]=counts.get(status,0)+1
            elapsed=time.monotonic()-started
            progress=dict(done=done,total=total,counts_this_invocation=counts,elapsed_seconds=elapsed,
                          last=dict(label=label,trial=trial,status=status,seconds=result['wall_seconds']))
            atomic_json(out/'progress.json',progress)
            print(json.dumps(progress),flush=True)
    print(json.dumps(dict(event='complete',done=done,total=total,elapsed_seconds=time.monotonic()-started)),flush=True)


if __name__=='__main__':
    main()
