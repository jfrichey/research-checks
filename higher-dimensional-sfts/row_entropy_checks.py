"""Exact audits for hard-endpoint amplification along rectangular row sweeps.

Research implementation developed by ChatGPT from Jacob Richey's direction.
Checks universal collars for inhomogeneous prescriptions, adjacent-stage
collars equal off the switched row, and the local resolvent witnesses.
The limiting entropy assertion requires the separate analytic proof.
"""
import argparse
from collections import Counter
from itertools import product
import json
from pathlib import Path
import random

from ordering_checks import agreement, displacements, rectangle, distribution

ROOT = Path(__file__).parent
EVIDENCE = ROOT / 'row-entropy-evidence.json'


def face_colors(word):
    faces = [word[0], word[-1], tuple(row[0] for row in word),
             tuple(row[-1] for row in word)]
    return {face[0] for face in faces if len(set(face)) == 1}


def background_colors(old, new, q):
    before = [c for c in range(q) if c not in face_colors(old)]
    after = [c for c in range(q) if c not in face_colors(new)]
    forbidden = {(new[0][j], old[-1][j]) for j in (0, -1)}
    return next((a, b) for a in after for b in before if (a, b) not in forbidden)


def guard(old, new, q, inner, cut):
    """Interior origin is (gap,gap); global rows <=cut use new requests."""
    h, w = len(old), len(old[0])
    gap = 4*max(h, w)
    assert gap-1 <= cut < gap+inner[0]
    a, b = background_colors(old, new, q)
    height, width = inner[0]+2*gap, inner[1]+2*gap
    out = [[a if r <= cut else b for _ in range(width)] for r in range(height)]
    for r in range(gap, gap+inner[0]):
        for c in range(gap, gap+inner[1]):
            out[r][c] = None
    for vertical, horizontal in product((False, True), repeat=2):
        word, color = (old, b) if vertical else (new, a)
        word = [list(row[::-1] if horizontal else row)
                for row in (word[::-1] if vertical else word)]
        if word[-1][-1] != color:
            continue
        col, row = [line[-1] for line in word], word[-1]
        corner = next(c for c in range(q) if c != color and
                      all(col != [c]*k+[color]*(h-k) for k in range(1, h)))
        bridge = next(c for c in range(q) if c not in (row[0], row[-1]))
        sequence = [corner]*w+[bridge]*w+[color]*w
        for r in range(h):
            for c, value in enumerate(sequence):
                out[height-1-r if vertical else r][width-1-c if horizontal else c] = value
    return out, gap


def verify_pins(old, new, pins, cut, exempt=lambda r, c: False):
    height, width = len(pins), len(pins[0])
    h, w = len(old), len(old[0]); checked = 0
    for r in range(1-h, height):
        for c in range(1-w, width):
            if exempt(r, c):
                continue
            assert any(0 <= r+i < height and 0 <= c+j < width and
                       pins[r+i][c+j] is not None and
                       pins[r+i][c+j] != (new if r+i <= cut else old)[i][j]
                       for i in range(h) for j in range(w)), (old, new, cut, r, c)
            checked += 1
    return checked


def verify_guard(old, new, q, inner, relative_cut):
    gap = 4*max(len(old), len(old[0])); cut = gap+relative_cut
    pins, gap = guard(old, new, q, inner, cut)
    return verify_pins(old, new, pins, cut,
        lambda r, c: gap <= r and r+len(old) <= gap+inner[0] and
                     gap <= c and c+len(old[0]) <= gap+inner[1])


def face_audit():
    counts = Counter()
    for flat_u in product(range(3), repeat=4):
        old = [flat_u[:2], flat_u[2:]]
        for flat_v in product(range(3), repeat=4):
            new = [flat_v[:2], flat_v[2:]]
            a, b = background_colors(old, new, 3)
            for j in (0, 1):
                for k in range(3):
                    face = [new[i][j] if i < k else old[i][j] for i in range(2)]
                    assert face != [a]*k+[b]*(2-k)
            counts['pattern_pairs'] += 1
            counts['blocked_windows'] += verify_guard(old, new, 3, (2, 3), 0)
    return dict(counts)


def normalized_patterns(shape, u, v, normal):
    """Move an axis sweep to ascending physical rows."""
    h, w = shape; u, v = tuple(map(int, u)), tuple(map(int, v))
    if normal[0]:
        index = [((r if normal[0] > 0 else h-1-r)*w+c)
                 for r in range(h) for c in range(w)]
        dims = h, w
    else:
        index = [(r*w+(c if normal[1] > 0 else w-1-c))
                 for c in range(w) for r in range(h)]
        dims = w, h
    h, w = dims
    words = [tuple(word[i] for i in index) for word in (u, v)]
    return tuple(tuple(tuple(word[r*w:(r+1)*w]) for r in range(h)) for word in words)


def conditional_law(supports, words, n, q):
    """Clique expansion using literal requests, independent of G masks."""
    moments = [0]*(len(words)+1)
    for subset in range(1 << len(words)):
        letters = {}; compatible = True
        for i, (support, word) in enumerate(zip(supports, words)):
            if not subset >> i & 1:
                continue
            for x, value in zip(support, word):
                if x in letters and letters[x] != value:
                    compatible = False; break
                letters[x] = value
            if not compatible:
                break
        if compatible:
            moments[subset.bit_count()] += q**(n-len(letters))
    return distribution(moments)


def local_witness(shape, u, v, normal, q=3):
    old, new = normalized_patterns(shape, u, v, normal)
    h, w = len(old), len(old[0]); sites = rectangle(h, w)
    old_g, new_g = [agreement(sites, sum(word, ())) for word in (old, new)]
    assert old_g != new_g and old_g & new_g == old_g
    dr, dc = max((d for k, d in enumerate(displacements(sites)) if (new_g ^ old_g) >> k & 1),
                 key=lambda d: (h-abs(d[0]))*(w-abs(d[1])))
    if dc < 0:
        old, new = [tuple(row[::-1] for row in word) for word in (old, new)]
        dc = -dc
    assert dr >= 0 and dc >= 0
    row = max(r for r, c in sites if r >= dr and c >= dc and
              old[r][c] != old[r-dr][c-dc])
    union = {}
    for ar, ac in ((0, 0), (dr, dc)):
        for r, c in sites:
            z = ar+r, ac+c
            value = (new if z[0] <= row else old)[r][c]
            assert z not in union or union[z] == value
            union[z] = value
    inner = h+dr, w+dc
    patch = [[None]*inner[1] for _ in range(inner[0])]
    upper = next(c for c in range(q) if c != new[0][-1])
    lower = next(c for c in range(q) if c != old[-1][0])
    for r, c in product(range(inner[0]), range(inner[1])):
        patch[r][c] = union[r, c] if (r, c) in union else upper if r < dr else lower
    active, supports, before, after = [], [], [], []
    for ar in range(dr+1):
        for ac in range(dc+1):
            if all(patch[ar+r][ac+c] == (new if ar+r < row else old)[r][c]
                   for r, c in sites if ar+r != row):
                assert ar <= row < ar+h
                assert all((ar+r, ac+c) in union for r, c in sites)
                active.append((ar, ac)); supports.append(tuple(range(ac, ac+w)))
                before.append(old[row-ar]); after.append(new[row-ar])
    assert (0, 0) in active and (dr, dc) in active
    laws = [conditional_law(supports, words, inner[1], q) for words in (before, after)]
    difference = [b-a for a, b in zip(*laws)]
    assert difference[:3] == [1, -2, 1] and not any(difference[3:])
    direct = q**inner[1] <= 25000
    if direct:
        for words, expected in zip((before, after), laws):
            observed = [0]*len(expected)
            for letters in product(range(q), repeat=inner[1]):
                matches = sum(all(letters[x] == value for x, value in zip(s, word))
                              for s, word in zip(supports, words))
                observed[matches] += 1
            assert observed == expected
    avoiding = next(letters for letters in product(range(q), repeat=inner[1])
                    if not any(all(letters[x] == value for x, value in zip(s, word))
                               for s, word in zip(supports, after)))
    gap = 4*max(h, w)
    pins_before, _ = guard(old, new, q, inner, gap+row-1)
    pins_after, _ = guard(old, new, q, inner, gap+row)
    assert all(a == b for r, (aa, bb) in enumerate(zip(pins_before, pins_after))
               if r != gap+row for a, b in zip(aa, bb))
    checks = 0
    for pins, cut in ((pins_before, gap+row-1), (pins_after, gap+row)):
        checks += verify_pins(old, new, pins, cut,
            lambda r, c: gap <= r and r+h <= gap+inner[0] and
                         gap <= c and c+w <= gap+inner[1])
        for r, c in product(range(inner[0]), range(inner[1])):
            pins[gap+r][gap+c] = patch[r][c]
    # New pair extension has exactly its two designated occurrences.
    checks += verify_pins(old, new, pins_after, gap+row,
                          lambda r, c: (r, c) in ((gap, gap), (gap+dr, gap+dc)))
    # Old suffix witness: impose endpoint j and exclude every later event.
    for c in range(inner[1]):
        pins_before[gap+row][gap+c] = old[row-dr][c-dc] if dc <= c < dc+w else 0
    j_key = gap+dc, gap+dr
    checks += verify_pins(old, new, pins_before, gap+row-1,
                          lambda r, c: (c, r) <= j_key)
    # A new-avoiding filling gives the conditional probability bound.
    for c, value in enumerate(avoiding):
        pins_after[gap+row][gap+c] = value
    checks += verify_pins(old, new, pins_after, gap+row)
    return dict(shape=list(shape), g=old_g, h=new_g, displacement=[dr, dc],
                switched_row=row, active_anchors=[list(z) for z in active],
                old_conditional_avoidance=laws[0][0], new_conditional_avoidance=laws[1][0],
                literal_row_audit=direct, collar_shape=[len(pins_after), len(pins_after[0])],
                pinned_window_checks=checks)


def audit():
    seed = 20261020; rng = random.Random(seed)
    faces = face_audit()
    print('ROW BACKGROUND / COLLAR AUDIT', faces, flush=True)
    random_windows = 0
    for _ in range(250):
        q = rng.randrange(3, 6); h, w = rng.randrange(2, 7), rng.randrange(2, 7)
        old, new = [tuple(tuple(rng.randrange(q) for _ in range(w)) for _ in range(h))
                    for _ in range(2)]
        inner = rng.randrange(1, 8), rng.randrange(1, 8)
        random_windows += verify_guard(old, new, q, inner, rng.randrange(inner[0]))
    searches = json.loads((ROOT/'ternary-layer-evidence.json').read_text())['systems']
    remaining = {(tuple(r['shape']), g, h) for r in searches for g, h, *_ in r['failures']}
    data = json.loads((ROOT/'alphabet-sweep-evidence.json').read_text())
    witnesses = []
    for record in data['critical_lines']:
        if record['q'] != 3:
            continue
        for g, h, u, v, a, b in record['generators']:
            if (tuple(record['shape']), g, h) in remaining:
                assert (a == 0) != (b == 0)
                witness = local_witness(record['shape'], u, v, (a, b))
                witness['original_g'], witness['original_h'] = g, h
                witnesses.append(witness)
                print('ROW ENTROPY WITNESS', record['shape'], g, h,
                      witness['old_conditional_avoidance'], witness['new_conditional_avoidance'], flush=True)
    assert len(witnesses) == len(remaining) == 32
    return dict(status='Exact finite collar and local-row audits; limiting entropy uses the analytic proof.',
                exhaustive=faces, random=dict(seed=seed, pattern_pairs=250, blocked_windows=random_windows),
                critical_ternary_witnesses=witnesses)


def main():
    parser = argparse.ArgumentParser(); parser.add_argument('--rebuild', action='store_true')
    args = parser.parse_args(); result = audit()
    if args.rebuild:
        EVIDENCE.write_text(json.dumps(result, indent=2)+'\n')
    else:
        assert result == json.loads(EVIDENCE.read_text())
    print('ROW ENTROPY AUDITS VERIFIED', flush=True)


if __name__ == '__main__':
    main()
