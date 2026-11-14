"""Exact agreement-alphabet spectra and all-alphabet line certificates.

Research code developed by ChatGPT from Jacob Richey's research direction.
Only standard-library integer calculations are used. --rebuild searches and
saves; the default reconstructs the agreement closure system, checks minimum
alphabets and independently replays every literal line certificate. The
all-volume consequences use the analytic interval and strictification proofs.
"""

import argparse
from collections import Counter
from functools import lru_cache
from itertools import product
import json
from math import gcd
from pathlib import Path

from ordering_checks import agreement, displacements, rectangle
from merger_reachability import binary_masks

ROOT = Path(__file__).parent
EVIDENCE = ROOT / 'alphabet-sweep-evidence.json'
SMALL = [(2, 2), (2, 3), (2, 4), (3, 3), (2, 5), (2, 6), (3, 4)]
PLANAR = [(h, w) for h in range(2, 6) for w in range(h, 25//h+1)
          if h*w <= 25]
SPECTRA = PLANAR + [(2, 2, 2), (2, 2, 3)]


def geometry(dims):
    sites = tuple(product(*(range(a) for a in dims)))
    index = {z: i for i, z in enumerate(sites)}
    ds = sorted({tuple(a-b for a, b in zip(x, y))
                 for x in sites for y in sites if x > y})
    pairs = [[(i, index[z]) for i, x in enumerate(sites)
              if (z := tuple(a-b for a, b in zip(x, d))) in index] for d in ds]
    return sites, ds, pairs


def mask_of(pairs, word):
    return sum(1 << k for k, edges in enumerate(pairs)
               if all(word[a] == word[b] for a, b in edges))


def close(n, pairs, mask):
    """Join every equality required by mask, then add all implied periods."""
    parent = list(range(n))

    def root(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    for k, edges in enumerate(pairs):
        if mask >> k & 1:
            for a, b in edges:
                parent[root(a)] = root(b)
    ids, labels = {}, []
    for i in range(n):
        r = root(i)
        if r not in ids:
            ids[r] = len(ids)
        labels.append(ids[r])
    return mask_of(pairs, labels), tuple(labels)


def all_classes(n, pairs):
    """Every closed G is reached by adding its displacements one at a time."""
    g, labels = close(n, pairs, 0)
    classes, queue = {g: labels}, [g]
    for g in queue:
        for k in range(len(pairs)):
            if not g >> k & 1:
                h, labels = close(n, pairs, g | (1 << k))
                if h not in classes:
                    classes[h] = labels
                    queue.append(h)
    return classes


def clauses_for(labels, pairs, g):
    return sorted({tuple(sorted({tuple(sorted((labels[a], labels[b])))
                                 for a, b in edges if labels[a] != labels[b]}))
                   for k, edges in enumerate(pairs) if not g >> k & 1})


def realize(labels, pairs, g, q):
    """Exhaust component colorings modulo alphabet permutations, with pruning."""
    clauses = clauses_for(labels, pairs, g)
    assert all(clauses)
    c = max(labels)+1
    score = Counter()
    for clause in clauses:
        for a, b in clause:
            score[a] += 10000//len(clause)**2
            score[b] += 10000//len(clause)**2
    order = sorted(range(c), key=lambda x: (-score[x], x))
    inverse = {a: i for i, a in enumerate(order)}
    at = [[] for _ in range(c)]
    for clause in clauses:
        cl = [(inverse[a], inverse[b]) for a, b in clause]
        at[max(max(a, b) for a, b in cl)].append(cl)
    word = [0]*c

    def visit(i, top):
        if i == c:
            return tuple(word[inverse[a]] for a in labels)
        for x in range(min(top+2, q)):
            word[i] = x
            if all(any(word[a] != word[b] for a, b in cl) for cl in at[i]):
                result = visit(i+1, max(top, x))
                if result is not None:
                    return result
        return None

    return visit(0, -1)


def spectrum(dims):
    sites, ds, pairs = geometry(dims)
    classes = all_classes(len(sites), pairs)
    records = []
    for g, labels in sorted(classes.items()):
        for q in range(1, max(labels)+2):
            word = realize(labels, pairs, g, q)
            if word is not None:
                break
        assert mask_of(pairs, word) == g
        records.append([g, q, ''.join(map(str, word))])
    return dict(shape=list(dims), classes=records,
                minimum_alphabets=dict(sorted(Counter(r[1] for r in records).items())))


def bell_audit(dims):
    """Independent exhaustive equality partitions for supports of size <=9."""
    sites, _, pairs = geometry(dims)
    n, word, minimum = len(sites), [0]*len(sites), {}
    count = 0

    def visit(i, top):
        nonlocal count
        if i == n:
            count += 1
            g = mask_of(pairs, word)
            minimum[g] = min(top+1, minimum.get(g, n+1))
            return
        for x in range(top+2):
            word[i] = x
            visit(i+1, max(top, x))

    visit(1, 0)
    return count, minimum


def covers(classes):
    result = []
    for g in classes:
        above = sorted((h for h in classes if h != g and g & h == g),
                       key=lambda h: (h.bit_count(), h))
        minimal = []
        for h in above:
            if not any(z & h == z for z in minimal):
                minimal.append(h)
                result.append((g, h))
    return result


def search_lines(h, w):
    sites, _, pairs = geometry((h, w))
    word, profiles = [0]*len(sites), {}

    def visit(i, top):
        if i == len(sites):
            conflicts = tuple(sum(1 << a for a, b in edges if word[a] != word[b])
                              for edges in pairs)
            g = sum(1 << k for k, c in enumerate(conflicts) if not c)
            profiles.setdefault(g, []).append((tuple(word), conflicts))
            return
        for x in range(min(top+2, 3)):
            word[i] = x
            visit(i+1, max(top, x))

    visit(1, 0)
    assert set(profiles) == set(all_classes(len(sites), pairs))
    normals = [(1, 0), (-1, 0), (0, 1), (0, -1)]
    normals += [(a, b) for a in range(1-w, w) for b in range(1-h, h)
                if gcd(a, b) == 1 and (a, b) not in normals]
    keys = [[a*r+b*c for r, c in sites] for a, b in normals]

    @lru_cache(None)
    def allowed(old, new):
        assert old and new
        if old & new:
            return (1 << len(normals))-1
        oi = [i for i in range(h*w) if old >> i & 1]
        ni = [i for i in range(h*w) if new >> i & 1]
        return sum((min(key[i] for i in ni) <= max(key[i] for i in oi)) << k
                   for k, key in enumerate(keys))

    generators = []
    for g, target in covers(sorted(profiles)):
        certificate = None
        for u, old in profiles[g]:
            for v, new in profiles[target]:
                choices = (1 << len(normals))-1
                for a, b in zip(old, new):
                    if b:
                        choices &= allowed(a, b)
                    if not choices:
                        break
                if choices:
                    normal = normals[(choices & -choices).bit_length()-1]
                    certificate = [g, target, ''.join(map(str, u)),
                                   ''.join(map(str, v)), *normal]
                    break
            if certificate:
                break
        assert certificate is not None, (h, w, g, target)
        generators.append(certificate)
    return dict(shape=[h, w], generators=generators,
                proper_inclusions=sum(g != k and g & k == g
                                      for g in profiles for k in profiles))


def independent_binary_classes(dims):
    h, w = dims
    if h*w <= 16:
        return set(binary_masks(rectangle(h, w)))
    for name, key in [('switching-evidence.json', 'systems'),
                      ('line-sweep-evidence.json', 'systems')]:
        data = json.loads((ROOT/name).read_text())
        for record in data[key]:
            if record['shape'] == [h, w]:
                if 'classes' in record:
                    return set(record['classes'])
                # Every class participates in a saved generator or is constant.
                edges = record['generators']+record['row_generators']
                return {g for edge in edges for g in edge[:2]}
    raise AssertionError(('No independent binary enumeration', dims))


def cover_merger_audit(data):
    """Split one old forced component off an upper-class representative."""
    spectra = {tuple(r['shape']): r for r in data['spectra']}
    counts, covers_checked = Counter(), 0
    for record in data['lines']:
        dims = tuple(record['shape'])
        sites, _, pairs = geometry(dims)
        words = {g: tuple(map(int, w)) for g, q, w in spectra[dims]['classes']}
        for g, h, *_ in record['generators']:
            _, labels = close(len(sites), pairs, g)
            target = words[h]
            fresh = max(target)+1
            successes = 0
            for component in set(labels):
                colors = {target[i] for i, c in enumerate(labels) if c == component}
                assert len(colors) == 1
                original = next(iter(colors))
                source = tuple(fresh if c == component else target[i]
                               for i, c in enumerate(labels))
                actual = agreement(sites, source)
                assert g & actual == g and actual & h == actual
                if actual != h:
                    assert actual == g
                    assert tuple(original if c == fresh else c for c in source) == target
                    successes += 1
                    counts[fresh+1] += 1
            assert successes
            covers_checked += 1
    return dict(covers=covers_checked, valid_splits_by_alphabet=dict(sorted(counts.items())))


def alphabet_bound_audit(data):
    classes, mergers = 0, 0
    for record in data['spectra']:
        _, ds, pairs = geometry(record['shape'])
        for g, q, letters in record['classes']:
            assert q*(q-1)//2 <= len(ds)-g.bit_count()
            word = tuple(map(int, letters))
            used = 0
            for a in range(q):
                for b in range(a+1, q):
                    merged = tuple(a if c == b else c for c in word)
                    target = mask_of(pairs, merged)
                    assert g & target == g and target != g
                    new = target ^ g
                    assert not new & used
                    used |= new
                    mergers += 1
            classes += 1
    return dict(classes=classes, distinct_pair_mergers=mergers)


def lower_covers(classes, target):
    if not target:
        return []
    result, g = [], (target-1) & target
    while True:
        if g in classes and not any(g & h == g for h in result):
            result.append(g)
        if not g:
            return result
        g = (g-1) & target


def search_critical_lines(data, q):
    """Only upper classes needing every ambient symbol escape cover mergers."""
    results = []
    for spectrum_record in data['spectra']:
        if len(spectrum_record['shape']) != 2:
            continue
        h, w = spectrum_record['shape']
        sites, _, pairs = geometry((h, w))
        words = {g: tuple(map(int, letters))
                 for g, r, letters in spectrum_record['classes'] if r <= q}
        targets = [g for g, r, _ in spectrum_record['classes'] if r == q]
        normals = [(1, 0), (-1, 0), (0, 1), (0, -1)]
        normals += [(a, b) for a in range(1-w, w) for b in range(1-h, h)
                    if gcd(a, b) == 1 and (a, b) not in normals]
        keys = [[a*r+b*c for r, c in sites] for a, b in normals]

        @lru_cache(None)
        def conflicts(word):
            return tuple(tuple(i for i, j in edges if word[i] != word[j])
                         for edges in pairs)

        def test(g, target, u, v):
            old, new = conflicts(u), conflicts(v)
            for normal, key in zip(normals, keys):
                if all(not b or min(key[i] for i in b) <= max(key[i] for i in a)
                       for a, b in zip(old, new)):
                    return [g, target, ''.join(map(str, u)), ''.join(map(str, v)), *normal]
            return None

        generators = []
        for target in targets:
            v = words[target]
            for g in lower_covers(words, target):
                certificate = test(g, target, words[g], v)
                if certificate is None:
                    _, labels = close(len(sites), pairs, g)
                    guesses = set()
                    for component in set(labels):
                        for color in range(q):
                            u = tuple(color if c == component else v[i]
                                      for i, c in enumerate(labels))
                            if mask_of(pairs, u) == g:
                                guesses.add(u)
                    for u in sorted(guesses):
                        certificate = test(g, target, u, v)
                        if certificate:
                            break
                assert certificate is not None, (h, w, q, g, target)
                generators.append(certificate)
        results.append(dict(shape=[h, w], q=q, generators=generators))
        print('CRITICAL LINES BUILT', (h, w), q, len(generators), flush=True)
    return results


def verify_critical_lines(data):
    spectra = {tuple(r['shape']): r for r in data['spectra']}
    totals = Counter()
    for record in data['critical_lines']:
        dims, q = tuple(record['shape']), record['q']
        shape = rectangle(*dims)
        classes = {g for g, r, w in spectra[dims]['classes'] if r <= q}
        targets = [g for g, r, w in spectra[dims]['classes'] if r == q]
        expected = {(g, h) for h in targets for g in lower_covers(classes, h)}
        assert {(r[0], r[1]) for r in record['generators']} == expected
        assert len(record['generators']) == len(expected)

        @lru_cache(None)
        def literal(letters):
            word = tuple(map(int, letters))
            assert max(word) < q
            labels = dict(zip(shape, word))
            mismatches = tuple(tuple((x, y) for x, y in shape
                              if (x-dx, y-dy) in labels
                              and labels[x, y] != labels[x-dx, y-dy])
                              for dx, dy in displacements(shape))
            return agreement(shape, word), mismatches

        for g, target, u, v, a, b in record['generators']:
            assert g != target and g & target == g and gcd(a, b) == 1
            # Independent cover audit by all subsets of the gained directions.
            gain = target ^ g
            sub = (gain-1) & gain
            while sub:
                assert (g | sub) not in classes
                sub = (sub-1) & gain
            gu, old = literal(u)
            gv, new = literal(v)
            assert gu == g and gv == target
            key = lambda z: a*z[0]+b*z[1]
            assert all(not after or before and min(map(key, after)) <= max(map(key, before))
                       for before, after in zip(old, new))
        totals[q] += len(record['generators'])
        print('CRITICAL LINES VERIFIED', dims, q, len(record['generators']), flush=True)
    assert len(data['critical_lines']) == 2*len(PLANAR)
    assert {(tuple(r['shape']), r['q']) for r in data['critical_lines']} == {
        (dims, q) for dims in PLANAR for q in (3, 4)}
    return dict(sorted(totals.items()))


def antipodal_five_audit():
    """Independent coloring-partition audit for one five-dimensional cube G."""
    from itertools import combinations
    n, full = 16, (1 << 16)-1
    adjacency = [sum(1 << ((x ^ (1 << k)) if k < 4 else x ^ 15)
                     for k in range(5)) for x in range(n)]
    independent = [mask for mask in range(1, 1 << n)
                   if all(not (adjacency[i] & mask)
                          for i in range(n) if mask >> i & 1)]
    assert max(m.bit_count() for m in independent) == 5
    containing = [[m for m in independent if m >> i & 1] for i in range(n)]
    _, ds, pairs = geometry((2,)*5)
    body = sum(1 << i for i, d in enumerate(ds) if all(abs(x) == 1 for x in d))
    faces = []
    for a, b in combinations(range(5), 2):
        for fixed in range(32):
            if fixed & ((1 << a) | (1 << b)):
                continue
            face = [fixed, fixed ^ (1 << a), fixed ^ (1 << b),
                    fixed ^ (1 << a) ^ (1 << b)]
            faces.append([min(z, 31-z) for z in face])
    parts, histogram, survivors = [], Counter(), Counter()

    def visit(left, remaining):
        if not left:
            colors = [next(c for c, m in enumerate(parts) if m >> i & 1)
                      for i in range(n)]
            bad = sum(colors[z[0]] == colors[z[3]] and colors[z[1]] == colors[z[2]]
                      for z in faces)
            histogram[len(parts), bad] += 1
            if not bad:
                word = tuple(colors[min(i, 31-i)] for i in range(32))
                actual = mask_of(pairs, word)
                assert actual & body == body and actual != body
                extra = actual ^ body
                profile = tuple(sum(1 for j, d in enumerate(ds)
                                    if extra >> j & 1 and d.count(0) == k)
                                for k in range(1, 5))
                survivors[profile] += 1
            return
        if not remaining or left.bit_count() > 5*remaining:
            return
        i = (left & -left).bit_length()-1
        for block in containing[i]:
            if block & left == block:
                parts.append(block)
                visit(left ^ block, remaining-1)
                parts.pop()

    visit(full, 4)
    word = (0,1,1,2,1,2,3,0,1,2,3,4,4,0,2,3,
            3,2,0,4,4,3,2,1,0,3,2,1,2,1,1,0)
    assert mask_of(pairs, word) == body and len(set(word)) == 5
    return dict(independent_sets=len(independent),
                independent_size_histogram=dict(sorted(Counter(m.bit_count() for m in independent).items())),
                proper_coloring_counts=[[q, bad, count] for (q, bad), count in sorted(histogram.items())],
                square_free_survivors=[[list(profile), count] for profile, count in sorted(survivors.items())],
                five_symbol_word=list(word))


def verify(data):
    counts = {}
    for record in data['spectra']:
        dims = tuple(record['shape'])
        rebuilt = spectrum(dims)
        assert json.loads(json.dumps(rebuilt)) == record
        minimum = {g: q for g, q, _ in record['classes']}
        if len(dims) == 2:
            assert {g for g, q in minimum.items() if q <= 2} == independent_binary_classes(dims)
        if len(geometry(dims)[0]) <= 9:
            bell, direct = bell_audit(dims)
            assert direct == minimum
            counts['x'.join(map(str, dims))] = bell
        print('SPECTRUM VERIFIED', dims, len(minimum), record['minimum_alphabets'], flush=True)
    assert counts == data['bell_partition_counts']
    spectra = {tuple(r['shape']): r for r in data['spectra']}
    for record in data['lines']:
        dims = tuple(record['shape'])
        shape = rectangle(*dims)
        classes = {row[0] for row in spectra[dims]['classes']}
        assert {(r[0], r[1]) for r in record['generators']} == set(covers(classes))
        assert len(record['generators']) == len(covers(classes))
        for g, target, uw, vw, a, b in record['generators']:
            u, v = tuple(map(int, uw)), tuple(map(int, vw))
            assert max(u+v) <= 2 and gcd(a, b) == 1
            assert agreement(shape, u) == g and agreement(shape, v) == target
            old, new = dict(zip(shape, u)), dict(zip(shape, v))
            for dx, dy in displacements(shape):
                before = [a*x+b*y for x, y in shape if (x-dx, y-dy) in old
                          and old[x, y] != old[x-dx, y-dy]]
                after = [a*x+b*y for x, y in shape if (x-dx, y-dy) in new
                         and new[x, y] != new[x-dx, y-dy]]
                assert not after or before and min(after) <= max(before)
        assert record['proper_inclusions'] == sum(g != h and g & h == g
                                                   for g in classes for h in classes)
        print('LINES VERIFIED', dims, len(record['generators']),
              record['proper_inclusions'], flush=True)
    # Explicit lower bounds: forced components have all six K4 inequalities.
    for dims, g in [((3, 7), 2686976544), ((2, 2, 2), 5200)]:
        sites, ds, pairs = geometry(dims)
        actual, labels = close(len(sites), pairs, g)
        assert actual == g and max(labels) == 3
        singleton_edges = {cl[0] for cl in clauses_for(labels, pairs, g) if len(cl) == 1}
        assert singleton_edges == {(a, b) for a in range(4) for b in range(a+1, 4)}
    print('K4 LOWER BOUNDS VERIFIED', flush=True)
    assert json.loads(json.dumps(cover_merger_audit(data))) == data['cover_merger_audit']
    print('COVER MERGERS VERIFIED', data['cover_merger_audit'], flush=True)
    assert alphabet_bound_audit(data) == data['alphabet_bound_audit']
    print('ALPHABET BOUND VERIFIED', data['alphabet_bound_audit'], flush=True)
    totals = verify_critical_lines(data)
    assert json.loads(json.dumps(totals)) == data['critical_line_totals']
    assert json.loads(json.dumps(antipodal_five_audit())) == data['antipodal_five_audit']
    print('FIVE-DIMENSIONAL ALPHABET OBSTRUCTION VERIFIED', flush=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--rebuild', action='store_true')
    args = parser.parse_args()
    if args.rebuild:
        spectra = []
        for dims in SPECTRA:
            record = spectrum(dims)
            spectra.append(record)
            print('SPECTRUM BUILT', dims, len(record['classes']), record['minimum_alphabets'], flush=True)
        lines = []
        for dims in SMALL:
            record = search_lines(*dims)
            lines.append(record)
            print('LINES BUILT', dims, len(record['generators']), record['proper_inclusions'], flush=True)
        data = dict(status='Exact finite structural certificates; all-volume implications use working proofs.',
                    spectra=spectra, lines=lines,
                    bell_partition_counts={'x'.join(map(str, dims)): bell_audit(dims)[0]
                                           for dims in SPECTRA if len(geometry(dims)[0]) <= 9})
        data['cover_merger_audit'] = cover_merger_audit(data)
        data['alphabet_bound_audit'] = alphabet_bound_audit(data)
        data['critical_lines'] = search_critical_lines(data, 3)+search_critical_lines(data, 4)
        data['critical_line_totals'] = dict(Counter({q: sum(len(r['generators'])
            for r in data['critical_lines'] if r['q'] == q) for q in (3, 4)}))
        data['antipodal_five_audit'] = antipodal_five_audit()
        data = json.loads(json.dumps(data))
        verify(data)
        EVIDENCE.write_text(json.dumps(data, indent=2)+'\n')
    else:
        verify(json.loads(EVIDENCE.read_text()))


if __name__ == '__main__':
    main()
