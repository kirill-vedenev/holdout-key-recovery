#!/usr/bin/env python3
"""Check result coverage, trial seeds, fixture hashes and status flags."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path


def seed_for(master,label,trial,phase):
    return int.from_bytes(hashlib.sha256(('%s|%s|%s|%s'%(master,label,trial,phase)).encode()).digest()[:16],'big')


def verify(root,require_complete=False):
    manifest=json.loads((root/'manifest.json').read_text());spec=manifest['spec']
    sha=hashlib.sha256(json.dumps(spec,sort_keys=True).encode()).hexdigest()
    assert sha==manifest['spec_sha256']
    start=spec.get('trial_start',0)
    expected={(P['label'],i) for P in spec['parameters'] for i in range(start,start+spec['trials'][P['label']])}
    seen=set();errors=[]; fixtures=0
    for p in sorted((root/'cases').glob('*/*.json')):
        r=json.loads(p.read_text());key=(r['label'],r['trial'])
        assert key in expected and key not in seen,(p,key)
        seen.add(key)
        assert p.parent.name==r['label'] and p.stem=='%03d'%r['trial']
        assert r['master_seed']==spec['master_seed']
        for phase,seed in r['seeds'].items():
            assert int(seed)==seed_for(spec['master_seed'],*key,phase),(p,phase)
        if 'fixture_sha256' in r:
            path=root/'fixtures'/r['label']/('%03d.json.gz'%r['trial'])
            data=gzip.decompress(path.read_bytes())
            assert hashlib.sha256(data).hexdigest()==r['fixture_sha256'],path
            obj=json.loads(data)
            assert obj['public_generator_sha256']==r['public_generator_sha256']
            assert obj['parameters']==r['parameters']
            fixtures+=1
        if r['status']=='pass':
            assert r['eligible'] and r['complete_rank'] and r['full_grs'] and r['public_ambient_contains'],p
            if spec['recovery']:
                assert r['support_recovered'] and r['embedding_verified'],p
        if r['status'] in ('error','timeout'):
            errors.append(key)
    report=dict(expected=len(expected),recorded=len(seen),missing=len(expected-seen),verified_fixtures=fixtures,
                execution_errors=errors,manifest_sha256=sha)
    print(json.dumps(report,indent=2))
    if require_complete:
        assert seen==expected, 'missing prescribed trials'
        assert not errors, 'execution errors need investigation'
    return report


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('directory',type=Path);p.add_argument('--require-complete',action='store_true')
    a=p.parse_args();verify(a.directory,a.require_complete)
