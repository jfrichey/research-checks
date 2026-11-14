#!/usr/bin/env python3
"""Build all project proofs and audit the axioms of every public theorem."""
import pathlib
import re
import subprocess

ROOT = pathlib.Path(__file__).resolve().parent
subprocess.run([str(ROOT / 'lake.sh'), 'build', 'Balloons'], cwd=ROOT, check=True)

names = []
for path in sorted((ROOT / 'Balloons').glob('*.lean')):
    source = path.read_text()
    if re.search(r'\b(?:sorry|admit|native_decide)\b|^\s*axiom\s', source, re.M):
        raise SystemExit(f'Unverified proof construct in {path.name}')
    names.extend('Balloons.' + name for name in re.findall(
        r'^(?:@\[[^\n]*\]\s*)?(?:theorem|lemma)\s+([\w.]+)', source, re.M))
if not names or len(names) != len(set(names)):
    raise SystemExit('Empty or duplicate theorem audit list')

audit = ROOT / '.lake' / 'axiom-audit.lean'
audit.write_text('module\nimport Balloons\n\n' +
                 '\n'.join('#print axioms ' + name for name in names) + '\n')
result = subprocess.run([str(ROOT / 'lake.sh'), 'env', 'lean', '-j1', str(audit)],
                        cwd=ROOT, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
(ROOT / '.lake' / 'axiom-audit.txt').write_text(result.stdout)
if result.returncode:
    print(result.stdout)
    raise SystemExit(result.returncode)
allowed = {'propext', 'Classical.choice', 'Quot.sound'}
for block in re.findall(r'depends on axioms:\s*\[([^]]*)\]', result.stdout, re.S):
    axioms = {a.strip() for a in block.split(',') if a.strip()}
    if not axioms <= allowed:
        raise SystemExit(f'Unexpected axioms: {sorted(axioms - allowed)}')
for name in names:
    if "'" + name + "'" not in result.stdout:
        raise SystemExit(f'Missing axiom report: {name}')
print(f'Checked {len(names)} public theorems; only standard Lean axioms occur.')
print('Full axiom report: .lake/axiom-audit.txt')
