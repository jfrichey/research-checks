"""Recompute and certify the final enclosure (several GB and many minutes)."""
import argparse
from fractions import Fraction as F
import gzip
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import time

sys.dont_write_bytecode = True
from certify_bounds import slope_interval, outward


def validate(results, claim):
    k = claim['width']
    lower = (slope_interval(results[f'square_free_{k}'])[0]
             - slope_interval(results[f'square_free_{k-1}'])[1])
    upper = (slope_interval(results[f'square_wired_{k}'])[1]
             - slope_interval(results[f'square_wired_{k-1}'])[0])
    assert F(claim['lower']) <= lower <= upper <= F(claim['upper'])
    diameter = F(claim['upper']) - F(claim['lower'])
    assert diameter == F(claim['diameter']) < F(claim['strict_width_bound'])
    return lower, upper, diameter


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--certificates-dir', type=Path,
                   help='Optional directory of compressed certificate-*-*.txt.gz files; skips floating iteration.')
    args = p.parse_args()
    root = Path(__file__).resolve().parents[1]
    claim = json.loads((root/'data/enclosure.json').read_text())
    results = {}
    with tempfile.TemporaryDirectory(prefix='verify-cluster-enclosure-') as tmp:
        work = Path(tmp)
        exe = work/'certify_layers'
        subprocess.run(['g++', '-O3', '-std=c++17', str(root/'code/certify_layers.cpp'),
                        '-o', str(exe)], check=True)
        for k in (claim['width']-1, claim['width']):
            for wired in (False, True):
                kind = 'wired' if wired else 'free'
                start = time.monotonic()
                print(f'Checking width {k}, {kind} boundary ...', flush=True)
                cmd = [str(exe), str(k), str(int(wired))]
                if args.certificates_dir:
                    name = f'certificate-{kind}-{k}.txt'
                    cert = work/name
                    with gzip.open(args.certificates_dir/(name+'.gz'), 'rb') as src, cert.open('wb') as dst:
                        shutil.copyfileobj(src, dst)
                    cmd += ['--verify', str(cert)]
                cp = subprocess.run(cmd, capture_output=True, text=True, check=True)
                results[f'square_{kind}_{k}'] = json.loads(cp.stdout)
                print(f'Exact residual certified ({time.monotonic()-start:.1f} seconds).', flush=True)
    lower, upper, diameter = validate(results, claim)
    print('Computed enclosure:', outward(lower, 10), outward(upper, 10, up=True))
    print(f"VERIFIED: {claim['lower']} <= lambda <= {claim['upper']}")
    print(f"Window width {float(diameter):.8g} < {claim['strict_width_bound']}")


if __name__ == '__main__':
    main()
