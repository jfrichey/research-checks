#!/usr/bin/env python3
"""Run explicitly listed research checks and report their stated scope."""
from pathlib import Path
import argparse
import hashlib
import json
import os
import subprocess
import sys

ROOT = Path(__file__).resolve().parent


def local_path(relative):
    path = (ROOT / relative).resolve()
    if not path.is_relative_to(ROOT) or not path.is_file():
        raise ValueError('Invalid companion path: ' + relative)
    return path


def verify_files():
    records = json.loads((ROOT / 'files.sha256.json').read_text())
    for relative, expected in records.items():
        actual = hashlib.sha256(local_path(relative).read_bytes()).hexdigest()
        if actual != expected:
            raise ValueError('File checksum differs from this release: ' + relative)
    print(f'PASS: {len(records)} companion file checksums', flush=True)


def main():
    if not __debug__:
        raise SystemExit('Run without -O: the mathematical scripts use assertions.')
    catalog = json.loads((ROOT / 'checks.json').read_text())
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('project', nargs='?', default='all', choices=['all', *catalog['projects']])
    parser.add_argument('--list', action='store_true')
    parser.add_argument('--verify-only', action='store_true')
    args = parser.parse_args()
    verify_files()
    if args.list:
        for name, project in catalog['projects'].items():
            print(name + ': ' + project['kind'] + '\n  ' + project['scope'])
        return
    if args.verify_only:
        return
    environment = dict(os.environ)
    environment.pop('PYTHONOPTIMIZE', None)
    environment['PYTHONDONTWRITEBYTECODE'] = '1'
    names = catalog['projects'] if args.project == 'all' else [args.project]
    for name in names:
        project = catalog['projects'][name]
        print('\n' + name + ': ' + project['kind'], flush=True)
        for command in project['commands']:
            print('Running ' + ' '.join(command), flush=True)
            subprocess.run([sys.executable, '-B', '-u', str(local_path(command[0])), *command[1:]],
                           cwd=ROOT, env=environment, check=True,
                           timeout=project.get('timeout_seconds_per_command', 900))
        print('PASS: ' + name + '\nScope: ' + project['limits'], flush=True)


if __name__ == '__main__':
    main()
