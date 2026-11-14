"""Exact checks of rectangular insertion into a fixed periodic background.

Research code developed by ChatGPT from Jacob Richey's research direction.
The all-size statement is a working proof in entropy-ordering.md. These checks
audit the face construction and literal boundary occurrences, using integers.
"""

import argparse
from collections import Counter
from itertools import product
import json
from pathlib import Path
import random


EVIDENCE = Path(__file__).with_name('rectangular-filling-evidence.json')


def background(shape, word, q):
    sites = list(product(*(range(n) for n in shape)))
    assert len(word) == len(sites) and all(n >= 2 for n in shape)
    faces = []
    for axis, n in enumerate(shape):
        for endpoint in (0, n-1):
            values = {value for site, value in zip(sites, word) if site[axis] == endpoint}
            if len(values) == 1:
                faces.append((axis, endpoint, next(iter(values))))
    colors = {c for _, _, c in faces}
    for c in range(q):
        if c not in colors:
            return ('constant', c)
    assert q == 2
    a = next(f for f in faces if f[2] == 0)
    b = next(f for f in faces if f[2] == 1)
    assert a[0] == b[0] and a[1] != b[1]
    return ('stripes', next(axis for axis in range(len(shape)) if axis != a[0]))


def background_value(scheme, site):
    kind, value = scheme
    return value if kind == 'constant' else site[value] % 2


def verify_faces(shape, word, q):
    scheme = background(shape, word, q)
    sites = list(product(*(range(n) for n in shape)))
    for axis, n in enumerate(shape):
        for endpoint in (0, n-1):
            face = [(site, value) for site, value in zip(sites, word) if site[axis] == endpoint]
            for phase in (0, 1):
                assert any(value != background_value(scheme, tuple(t+phase for t in site))
                           for site, value in face), (shape, word, scheme, axis, endpoint, phase)
    return scheme


def allowed(word, patch):
    h, w = len(word), len(word[0])
    a, b = len(patch), len(patch[0])
    return not any(all(patch[r+i][s+j] == word[i][j]
                       for i in range(h) for j in range(w))
                   for r in range(a-h+1) for s in range(b-w+1))


def verify_insertion(word, patches, q):
    h, w = len(word), len(word[0])
    scheme = background((h, w), tuple(c for row in word for c in row), q)
    assignments = {}
    anchors = set()
    for origin, patch in patches:
        assert allowed(word, patch)
        a, b = len(patch), len(patch[0])
        for i in range(a):
            for j in range(b):
                site = (origin[0]+i, origin[1]+j)
                assert site not in assignments
                assignments[site] = patch[i][j]
        anchors.update((origin[0]+r, origin[1]+s)
                       for r in range(1-h, a) for s in range(1-w, b))
    for r, s in anchors:
        assert any(assignments.get((r+i, s+j), background_value(scheme, (r+i, s+j)))
                   != word[i][j] for i in range(h) for j in range(w)), (word, patches, r, s)


def face_checks():
    records = []
    for q, shapes in [(2, [(2, 2), (2, 3), (2, 4), (3, 3), (3, 4), (4, 4),
                           (2, 2, 2), (2, 2, 3)]),
                      (3, [(2, 2), (2, 3), (3, 3), (2, 2, 2)])]:
        for shape in shapes:
            volume = 1
            for n in shape:
                volume *= n
            counts = Counter()
            for word in product(range(q), repeat=volume):
                counts[verify_faces(shape, word, q)[0]] += 1
            records.append(dict(q=q, shape=list(shape), patterns=q**volume, schemes=dict(counts)))
    return records


def exhaustive_insertions():
    records = []
    for h, w, a, b in [(2, 2, 3, 3), (2, 3, 3, 4), (3, 2, 4, 3)]:
        checked = 0
        for data in product(range(2), repeat=h*w):
            word = tuple(tuple(data[i*w:(i+1)*w]) for i in range(h))
            for data in product(range(2), repeat=a*b):
                patch = tuple(tuple(data[i*b:(i+1)*b]) for i in range(a))
                if allowed(word, patch):
                    verify_insertion(word, [((0, 0), patch)], 2)
                    checked += 1
        records.append(dict(pattern_shape=[h, w], patch_shape=[a, b], exact_insertions=checked))
    return records


def random_insertions():
    rng = random.Random(20261113)
    counts = Counter()
    for _ in range(3000):
        q = rng.randrange(2, 5)
        h, w = rng.randrange(2, 6), rng.randrange(2, 6)
        word = tuple(tuple(rng.randrange(q) for _ in range(w)) for _ in range(h))
        patches = []
        for __ in range(2):
            a, b = rng.randrange(1, 9), rng.randrange(1, 9)
            while True:
                patch = tuple(tuple(rng.randrange(q) for _ in range(b)) for _ in range(a))
                if allowed(word, patch):
                    break
            patches.append(patch)
        r, s = rng.randrange(-10, 11), rng.randrange(-10, 11)
        if rng.randrange(2):
            other = (r+len(patches[0])+h-1, s+rng.randrange(-3, 4))
        else:
            other = (r+rng.randrange(-3, 4), s+len(patches[0][0])+w-1)
        verify_insertion(word, [((r, s), patches[0]), (other, patches[1])], q)
        counts[str(q)] += 1
    return dict(seed=20261113, exact_two_patch_insertions=dict(counts))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--rebuild', action='store_true')
    args = parser.parse_args()
    result = dict(status='Exact face and insertion checks for the working rectangular '
                  'background theorem; no arbitrary-exterior filling assertion',
                  face_audits=face_checks(), exhaustive_insertions=exhaustive_insertions(),
                  random_insertions=random_insertions())
    if args.rebuild:
        EVIDENCE.write_text(json.dumps(result, indent=2)+'\n')
    else:
        assert result == json.loads(EVIDENCE.read_text())
    print(json.dumps(result, indent=2), flush=True)


if __name__ == '__main__':
    main()
