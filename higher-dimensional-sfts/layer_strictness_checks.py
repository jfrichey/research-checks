"""Exact local audits for strictification with a shared independent layer.

Research code developed by ChatGPT from Jacob Richey's research direction.
The working strictification proof combines these local identities with
the rectangular background, universal collar and all-anchor comparison proofs.
The uniform conditional amplification is proved separately in the note;
these finite checks alone do not establish the limiting entropy conclusion.
Standard library only. --rebuild writes evidence; default recomputes it.
"""

import argparse
from collections import Counter
from functools import lru_cache
from itertools import permutations, product
import json
from pathlib import Path

from ordering_checks import agreement, displacements, rectangle, finite_moments, distribution
from rectangular_filling_checks import background, background_value, verify_faces
from uniform_filling_checks import verify as verify_collar
import alphabet_sweep_checks as alphabet

EVIDENCE = Path(__file__).with_name('layer-strictness-evidence.json')


def reflect(word, h, w, d):
    return tuple(word[(h-1-r if d[0] < 0 else r)*w+(w-1-c if d[1] < 0 else c)]
                 for r in range(h) for c in range(w))


@lru_cache(None)
def marker_patch(h, w, word, k, dr, dc, audit=False):
    """After reflections dr,dc>=0. All marker occurrences lie inside U."""
    assert 0 <= dr < h and 0 <= dc < w and (dr or dc)
    sites = rectangle(h, w)
    labels = dict(zip(sites, word))
    union = {}
    for ar, ac in [(0, 0), (dr, dc)]:
        for (r, c), letter in labels.items():
            z = r+ar, c+ac
            assert z not in union or union[z] == letter
            union[z] = letter
    patch = dict(union)
    upper = next(a for a in range(k) if a != labels[0, w-1])
    lower = next(a for a in range(k) if a != labels[h-1, 0])
    for r, c in rectangle(h+dr, w+dc):
        if (r, c) in patch:
            continue
        if r < dr and c >= w:
            patch[r, c] = upper
        else:
            assert r >= h and c < dc
            patch[r, c] = lower
    scheme = background((h, w), word, k)

    def value(z):
        return patch.get(z, background_value(scheme, z))

    active = []
    windows = 0
    for ar, ac in product(range(1-h, h+dr), range(1-w, w+dc)):
        occurs = all(value((ar+r, ac+c)) == letter for (r, c), letter in labels.items())
        if occurs:
            assert all((ar+r, ac+c) in union for r, c in sites)
            active.append((ar, ac))
        windows += 1
    assert (0, 0) in active and (dr, dc) in active
    if audit:
        verify_faces((h, w), word, k)
    # R0=bounding box of U+(S-S). Only this marker restriction is pinned.
    r0 = {z: value(z) for z in product(range(1-h, 2*h+dr-1),
                                      range(1-w, 2*w+dc-1))}
    return dict(active=tuple(active), union=tuple(sorted(union)), marker=r0,
                inner=(3*h+dr-2, 3*w+dc-2), checked_windows=windows)


def local_identity(h, w, first_old, first_new, marker, r=2, k=2):
    shape = rectangle(h, w)
    ga, gb, gm = [agreement(shape, word) for word in (first_old, first_new, marker)]
    assert ga & gb == ga
    old, new = ga & gm, gb & gm
    assert old != new and old & new == old
    ds = displacements(shape)
    d = max((d for i, d in enumerate(ds) if (new ^ old) >> i & 1),
            key=lambda d: (h-abs(d[0]))*(w-abs(d[1])))
    aa, bb, mm = [reflect(word, h, w, d) for word in (first_old, first_new, marker)]
    dr, dc = abs(d[0]), abs(d[1])
    patch = marker_patch(h, w, mm, k, dr, dc)
    anchors = patch['active']
    laws = [distribution(finite_moments(shape, anchors, agreement(shape, word), r))
            for word in (aa, bb)]
    diff = [b-a for a, b in zip(*laws)]
    assert diff[:3] == [1, -2, 1] and not any(diff[3:])
    return dict(displacement=list(d), active_anchors=[list(z) for z in anchors],
                local_volume=len(patch['union']), old_avoidance= laws[0][0],
                new_avoidance=laws[1][0], marker_inner=list(patch['inner']))


def marker_audit():
    records = []
    for h, w, k in [(2, 2, 2), (2, 3, 2), (2, 4, 2), (3, 3, 2),
                    (2, 2, 3), (2, 3, 3)]:
        shape = rectangle(h, w)
        counts = Counter()
        for word in product(range(k), repeat=h*w):
            g = agreement(shape, word)
            for i, d in enumerate(displacements(shape)):
                if not g >> i & 1:
                    continue
                mm = reflect(word, h, w, d)
                patch = marker_patch(h, w, mm, k, abs(d[0]), abs(d[1]), True)
                counts['compatible_pairs'] += 1
                counts['translated_windows'] += patch['checked_windows']
                counts['extra_contained_occurrences'] += len(patch['active'])-2
        records.append(dict(shape=[h, w], k=k, **counts))
    return records


def small_layer_audit(data):
    records = []
    for h, w in [(2, 2), (2, 3), (2, 4), (3, 3)]:
        spectrum = next(s for s in data['spectra'] if s['shape'] == [h, w])
        words = {g: tuple(map(int, word)) for g, q, word in spectrum['classes'] if q <= 2}
        checked = 0
        for ga, first_old in words.items():
            for gb, first_new in words.items():
                if ga == gb or ga & gb != ga:
                    continue
                for gm, marker in words.items():
                    if ga & gm == gb & gm:
                        continue
                    local_identity(h, w, first_old, first_new, marker)
                    checked += 1
        records.append(dict(shape=[h, w], exact_layer_identities=checked))
    return records


def cover_certificates(data):
    spectra = {tuple(s['shape']): s for s in data['spectra']}
    records = []
    for lines in data['critical_lines']:
        if lines['q'] != 4 or not lines['generators']:
            continue
        h, w = lines['shape']
        shape = rectangle(h, w)
        spectrum = spectra[h, w]
        binary = {g: tuple(map(int, word)) for g, q, word in spectrum['classes'] if q <= 2}
        certificates = []
        for g, target, _, vw, *_ in lines['generators']:
            word = tuple(map(int, vw))
            bits = [tuple(c & 1 for c in word), tuple(c >> 1 for c in word),
                    tuple((c & 1) ^ (c >> 1) for c in word)]
            masks = [agreement(shape, bit) for bit in bits]
            certificate = None
            for i, j in permutations(range(3), 2):
                new_first_g, marker_g = masks[i], masks[j]
                assert new_first_g & marker_g == target
                first_g = next((a for a in binary
                                if a & new_first_g == a and a & marker_g == g), None)
                if first_g is None:
                    continue
                old_first, new_first, marker = binary[first_g], bits[i], bits[j]
                # Independent literal four-symbol agreement check.
                old_joint = tuple(2*a+m for a, m in zip(old_first, marker))
                new_joint = tuple(2*b+m for b, m in zip(new_first, marker))
                assert agreement(shape, old_joint) == g
                assert agreement(shape, new_joint) == target
                identity = local_identity(h, w, old_first, new_first, marker)
                joint_rows = [list(old_joint[i*w:(i+1)*w]) for i in range(h)]
                windows, patches = verify_collar(joint_rows, 4, tuple(identity['marker_inner']))
                certificate = dict(g=g, h=target, first_G=first_g,
                    new_first_G=new_first_g, marker_G=marker_g,
                    old_first=''.join(map(str, old_first)),
                    new_first=''.join(map(str, new_first)), marker=''.join(map(str, marker)),
                    local_identity=identity, universal_collar_windows=windows,
                    universal_collar_patches=patches)
                break
            assert certificate is not None, (h, w, g, target)
            certificates.append(certificate)
        records.append(dict(shape=[h, w], covers=certificates))
        print('LAYER COVERS VERIFIED', (h, w), len(certificates), flush=True)
    assert sum(len(r['covers']) for r in records) == 74
    return records


def audit():
    data = json.loads(alphabet.EVIDENCE.read_text())
    return dict(status='Exact finite local and factorization certificates; limiting entropy uses the separate working proof.',
                marker_isolation=marker_audit(), small_layer_identities=small_layer_audit(data),
                critical_four_covers=cover_certificates(data))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--rebuild', action='store_true')
    args = parser.parse_args()
    result = audit()
    if args.rebuild:
        EVIDENCE.write_text(json.dumps(result, indent=2)+'\n')
    else:
        assert result == json.loads(EVIDENCE.read_text())
    print('LAYER STRICTIFICATION VERIFIED', flush=True)


if __name__ == '__main__':
    main()
