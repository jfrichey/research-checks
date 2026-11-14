#!/usr/bin/env python3
"""Restore the exact saved certificate inputs; refuse to replace other data."""
from pathlib import Path
import gzip,hashlib,json
ROOT=Path(__file__).resolve().parent
for name,item in json.loads((ROOT/'certificate-inputs.json').read_text()).items():
    target=ROOT/name
    assert target.parent==ROOT and target.suffix=='.json'
    with gzip.open(ROOT/(name+'.gz'),'rb') as f:
        data=f.read(item['bytes']+1)
    if len(data)!=item['bytes'] or hashlib.sha256(data).hexdigest()!=item['sha256']:
        raise ValueError('Certificate identity mismatch: '+name)
    if target.exists() and target.read_bytes()!=data:
        raise FileExistsError('Refusing to replace different data: '+name)
    target.write_bytes(data)
    print('Verified and restored '+name)
