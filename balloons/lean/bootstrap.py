#!/usr/bin/env python3
"""Install the pinned Linux Lean toolchain and mathlib sources without Git mutations.

Everything installed stays in this directory; downloaded archives stay in /tmp.
The recorded SHA256 hashes must match before any archive is extracted.
"""
import concurrent.futures
import hashlib
import json
import pathlib
import shutil
import subprocess
import tarfile

ROOT = pathlib.Path(__file__).resolve().parent
DOWNLOADS = pathlib.Path('/tmp/balloons-lean-bootstrap')
PINS = json.loads((ROOT / 'dependency-pins.json').read_text())


def archive(url, name, digest):
    DOWNLOADS.mkdir(parents=True, exist_ok=True)
    path = DOWNLOADS / name
    if not path.exists():
        temporary = path.with_suffix(path.suffix + '.part')
        subprocess.run(['curl', '-fLsS', '--retry', '2', '-o', str(temporary), url], check=True)
        temporary.replace(path)
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    if h.hexdigest() != digest:
        raise RuntimeError(f'SHA256 mismatch: {path}')
    return path


def unpack_source(url, name, digest, target):
    if (target / '.balloons-archive-sha256').exists():
        if (target / '.balloons-archive-sha256').read_text().strip() == digest:
            return
    path = archive(url, name, digest)
    target.mkdir(parents=True, exist_ok=True)
    with tarfile.open(path) as tf:
        for member in tf.getmembers():
            parts = member.name.split('/', 1)
            if len(parts) != 2 or not parts[1]:
                continue
            member.name = parts[1]
            tf.extract(member, target, filter='data')
    (target / '.balloons-archive-sha256').write_text(digest + '\n')


def main():
    version = PINS['lean'].removeprefix('v')
    runtime = ROOT / '.toolchain' / f'lean-{version}-linux'
    if not (runtime / 'bin' / 'lean').exists():
        path = archive(
            f'https://github.com/leanprover/lean4/releases/download/v{version}/lean-{version}-linux.tar.zst',
            f'lean-{version}-linux.tar.zst', PINS['leanArchiveSha256'])
        runtime.parent.mkdir(parents=True, exist_ok=True)
        subprocess.run(['tar', '--zstd', '-xf', str(path), '-C', str(runtime.parent)], check=True)
    packages = ROOT / '.lake' / 'packages'
    mathlib = PINS['mathlib']
    unpack_source(
        'https://codeload.github.com/leanprover-community/mathlib4/tar.gz/refs/tags/' + mathlib['tag'],
        'mathlib-' + mathlib['tag'] + '.tar.gz', mathlib['archiveSha256'], packages / 'mathlib')

    def fetch(p):
        repo = p['url'].removeprefix('https://github.com/')
        unpack_source('https://codeload.github.com/' + repo + '/tar.gz/' + p['rev'],
                      p['name'] + '.tar.gz', p['archiveSha256'], packages / p['name'])

    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
        list(pool.map(fetch, PINS['packages']))
    # The upstream cache unpacks dependency artifacts relative to mathlib.
    # Point those paths at the same archive-backed packages used by this project.
    nested = packages / 'mathlib' / '.lake' / 'packages'
    nested.mkdir(parents=True, exist_ok=True)
    for p in PINS['packages']:
        link = nested / p['name']
        if link.exists() and not link.is_symlink():
            shutil.copytree(link, packages / p['name'], dirs_exist_ok=True)
            shutil.rmtree(link)
        if not link.is_symlink():
            link.symlink_to(pathlib.Path('..') / '..' / '..' / p['name'], target_is_directory=True)
    entries = [{'type': 'path', 'name': 'mathlib', 'dir': '.lake/packages/mathlib',
                'scope': '', 'inherited': False, 'configFile': 'lakefile.lean',
                'manifestFile': 'lake-manifest.json'}]
    entries.extend({'type': 'path', 'name': p['name'], 'dir': '.lake/packages/' + p['name'],
                    'scope': p['scope'], 'inherited': True, 'configFile': p['configFile'],
                    'manifestFile': p['manifestFile']} for p in PINS['packages'])
    manifest = {'version': '1.3.0', 'packagesDir': '.lake/packages', 'packages': entries,
                'name': 'balloons', 'lakeDir': '.lake'}
    (ROOT / 'lake-manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    subprocess.run([str(runtime / 'bin' / 'lean'), '--version'], check=True)
    print('Pinned toolchain and sources ready. Use ./lake.sh for this project.')


if __name__ == '__main__':
    main()
