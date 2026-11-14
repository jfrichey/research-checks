"""Exact universal-collar checks for one forbidden rectangle, q >= 3.

Research code developed by ChatGPT from Jacob Richey's research direction.
The all-size working proof is in entropy-ordering.md. Each tested translated
window must have a mismatch at a pinned collar site; unpinned letters inside
and outside are universally quantified, not sampled.
"""

import argparse
from collections import Counter
from itertools import product
import json
from pathlib import Path
import random


EVIDENCE = Path(__file__).with_name("uniform-filling-evidence.json")


def guard(word, q, inner):
    h, width = len(word), len(word[0])
    assert q >= 3 and h >= 2 and width >= 2 and min(inner) >= 1
    assert all(len(row) == width for row in word)
    assert all(0 <= value < q for row in word for value in row)
    faces = [word[0], word[-1], [row[0] for row in word],
             [row[-1] for row in word]]
    face_colors = {face[0] for face in faces if len(set(face)) == 1}
    assert len(face_colors) <= 2
    c = next(value for value in range(q) if value not in face_colors)
    gap = 4 * max(h, width)
    height, breadth = inner[0] + 2 * gap, inner[1] + 2 * gap
    out = [[c] * breadth for _ in range(height)]
    for i in range(gap, gap + inner[0]):
        for j in range(gap, gap + inner[1]):
            out[i][j] = None

    patches = []
    for vertical, horizontal in product((False, True), repeat=2):
        reflected = [list(row[::-1] if horizontal else row)
                     for row in (word[::-1] if vertical else word)]
        if reflected[-1][-1] != c:
            continue
        col, row = [line[-1] for line in reflected], reflected[-1]
        a = next(a for a in range(q) if a != c and
                 all(col != [a] * k + [c] * (h-k) for k in range(1, h)))
        bridge = next(b for b in range(q) if b not in (row[0], row[-1]))
        sequence = [a] * width + [bridge] * width + [c] * width
        for i in range(h):
            for j, value in enumerate(sequence):
                out[height-1-i if vertical else i][breadth-1-j if horizontal else j] = value
        patches.append(dict(vertical=vertical, horizontal=horizontal,
                            corner=a, bridge=bridge))
    return out, gap, c, patches


def verify(word, q, inner=(2, 3)):
    out, gap, c, patches = guard(word, q, inner)
    height, breadth = len(out), len(out[0])
    h, width = len(word), len(word[0])
    checked = 0
    for r in range(1-h, height):
        for s in range(1-width, breadth):
            if (gap <= r and r+h <= gap+inner[0] and
                    gap <= s and s+width <= gap+inner[1]):
                continue
            witness = any(0 <= r+i < height and 0 <= s+j < breadth and
                          out[r+i][s+j] is not None and
                          out[r+i][s+j] != word[i][j]
                          for i, j in product(range(h), range(width)))
            assert witness, dict(word=word, q=q, inner=inner, origin=(r, s),
                                 background=c, patches=patches)
            checked += 1
    return checked, len(patches)


def audit():
    records, totals = [], Counter()
    for h, width, q in ((2, 2, 3), (2, 3, 3), (3, 2, 3), (2, 2, 4)):
        count = Counter()
        for letters in product(range(q), repeat=h*width):
            word = [list(letters[i*width:(i+1)*width]) for i in range(h)]
            windows, patches = verify(word, q)
            count.update(patterns=1, blocked_windows=windows, patched_corners=patches)
        records.append(dict(shape=[h, width], q=q, inner=[2, 3], **count))
        totals.update(count)

    seed, rng = 20261120, random.Random(20261120)
    random_count = Counter()
    for _ in range(1000):
        q = rng.randrange(3, 7)
        h, width = rng.randrange(2, 7), rng.randrange(2, 7)
        word = [[rng.randrange(q) for _ in range(width)] for _ in range(h)]
        inner = rng.randrange(1, 8), rng.randrange(1, 8)
        windows, patches = verify(word, q, inner)
        random_count.update(patterns=1, blocked_windows=windows, patched_corners=patches)
    totals.update(random_count)
    return dict(status="Exact universal pinned-letter obstruction checks; "
                "the all-size UFP conclusion uses the working proof",
                exhaustive=records, random=dict(seed=seed, **random_count),
                totals=dict(totals))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--rebuild", action="store_true")
    args = parser.parse_args()
    result = audit()
    if args.rebuild:
        EVIDENCE.write_text(json.dumps(result, indent=2) + "\n")
    else:
        assert result == json.loads(EVIDENCE.read_text())
    print(json.dumps(result, indent=2), flush=True)


if __name__ == "__main__":
    main()
