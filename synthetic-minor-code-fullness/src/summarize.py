#!/usr/bin/env python3
"""Summarize all recorded trials, including failed assumptions and errors."""
import argparse
import collections
import json
import math
import statistics
from pathlib import Path


def wilson(success,total):
    if not total:
        return [None,None]
    z=1.959963984540054; p=success/total; den=1+z*z/total
    centre=(p+z*z/(2*total))/den
    radius=z*math.sqrt(p*(1-p)/total+z*z/(4*total*total))/den
    return [max(0,centre-radius),min(1,centre+radius)]


def field(P):
    if P['family']=='grs':
        return 'GF(%d)'%P['q']
    if P['family']=='binary':
        return 'GF(2^%d)'%P['m']
    return 'GF(%d) in GF(%d^%d)'%(P['q'],P['q'],P['m'])


def histogram(values):
    return dict(sorted(collections.Counter(values).items()))


def compact(h):
    return ', '.join('%s: %d'%(k,v) for k,v in h.items()) if len(h)>1 else (str(next(iter(h))) if h else '—')


def summarize(root):
    manifest=json.loads((root/'manifest.json').read_text())
    spec=manifest['spec']; records=[]; rows=[]
    for P in spec['parameters']:
        cases=[json.loads(p.read_text()) for p in sorted((root/'cases'/P['label']).glob('*.json'))]
        records.extend(cases)
        status=collections.Counter(c['status'] for c in cases)
        assessed=sum(c.get('eligible') is True and c.get('complete_rank') is True and 'full_grs' in c for c in cases)
        full=sum(c.get('eligible') is True and c.get('full_grs') is True for c in cases)
        times=[c['wall_seconds'] for c in cases]
        row=dict(label=P['label'],family=P['family'],field=field(P),
                 original_length=P.get('n0'),shortened_positions=P.get('ell',0) if P['family']=='binary' else None,
                 length=histogram(c['n'] for c in cases if 'n' in c),
                 dimension=histogram(c['k'] for c in cases if 'k' in c),
                 D=histogram(c['D'] for c in cases if 'D' in c),
                 predicted_bound=histogram(c['predicted_bound'] for c in cases if 'predicted_bound' in c),
                 minor_rank=histogram(c['minor_rank'] for c in cases if 'minor_rank' in c),
                 expected=spec['trials'][P['label']],recorded=len(cases),passed=status['pass'],
                 fullness_assessed=assessed,full_grs=full,
                 tangent_failures=status['tangent_assumption_failure'],
                 fullness_failures=status['fullness_failure'],
                 recovery_failures=status['recovery_failure'],errors=status['error'],timeouts=status['timeout'],
                 embedding_verified=sum(c.get('embedding_verified') is True for c in cases),
                 unique_public_generators=len({c['public_generator_sha256'] for c in cases if 'public_generator_sha256' in c}),
                 conditional_fullness_wilson95=wilson(full,assessed),
                 median_wall_seconds=statistics.median(times) if times else None,
                 max_wall_seconds=max(times) if times else None)
        rows.append(row)
    summary=dict(spec_sha256=manifest['spec_sha256'],expected=sum(spec['trials'].values()),
                 recorded=len(records),status_counts=dict(sorted(collections.Counter(c['status'] for c in records).items())),
                 sets=rows)
    (root/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    lines=['# Minor-code fullness: all trials','',
           '%d of %d prescribed trials recorded. Statuses: %s.'%(
               len(records),summary['expected'],', '.join('%s %d'%kv for kv in summary['status_counts'].items())),'',
           'Columns: `n0`/`ell` original length and shortened positions (binary Goppa), `n`, `k`, `D` of the '
           'code after shortening, `bound` the upper bound on the minor-code dimension, `rank` the certified '
           'dimension of the complete minor code. Several values are shown as `value: count`. '
           'Failures: tangent planes / fullness / recovery / errors and timeouts.','',
           '| Set | Field | n0 | ell | n | k | D | bound | rank | Pass | Failures | Median s |',
           '|---|---|---:|---:|---:|---:|---:|---:|---|---:|---|---:|']
    for r in rows:
        lines.append('| %s | %s | %s | %s | %s | %s | %s | %s | %s | %d/%d | %d / %d / %d / %d | %s |'%(
            r['label'],r['field'],r['original_length'] or '—',r['shortened_positions'] if r['shortened_positions'] is not None else '—',
            compact(r['length']),compact(r['dimension']),compact(r['D']),compact(r['predicted_bound']),
            compact(r['minor_rank']),r['passed'],r['recorded'],r['tangent_failures'],r['fullness_failures'],
            r['recovery_failures'],r['errors']+r['timeouts'],
            '%.3f'%r['median_wall_seconds'] if r['median_wall_seconds'] is not None else '—'))
    lines+=['','Per-set Wilson 95% intervals for fullness, conditional on two-dimensional tangent planes '
            'and a certified rank, are in `summary.json`. Times are elapsed seconds per trial under '
            'parallel load, including generation and all checks.']
    (root/'summary.md').write_text('\n'.join(lines)+'\n')
    print(json.dumps({k:v for k,v in summary.items() if k!='sets'},indent=2))
    return summary


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('directory',type=Path)
    summarize(p.parse_args().directory)
