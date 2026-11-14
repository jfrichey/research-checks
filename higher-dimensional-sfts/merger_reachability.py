"""Exact reachability by equal-G relabeling and ternary letter mergers.

An edge is an actual merger of a three-symbol pattern. Thus a path proves
convex order at q=3 on EVERY finite anchor set via the merging theorem; this
is not a check of a single finite-volume law. Failure to find a path is only
a limitation of this proof method. The algorithm exhausts alphabet partitions,
not labeled words. Standard library only.
"""

import argparse
import json
from pathlib import Path

from ordering_checks import agreement, displacements, rectangle


def binary_masks(shape):
    """Agreement masks for all binary words, with a bit per support site."""
    ids = {z: i for i, z in enumerate(shape)}
    specs = []
    for dx, dy in displacements(shape):
        by_shift = {}
        for i, (x, y) in enumerate(shape):
            j = ids.get((x + dx, y + dy))
            if j is not None:
                assert j > i  # shape is lexicographically ordered
                by_shift[j-i] = by_shift.get(j-i, 0) | (1 << i)
        specs.append(tuple(by_shift.items()))
    masks = []
    for word in range(1 << len(shape)):
        masks.append(sum(1 << k for k, spec in enumerate(specs)
                         if all(((word ^ (word >> shift)) & sites) == 0
                                for shift, sites in spec)))
    return masks


def scan(shape):
    shape = tuple(sorted(shape))
    n = len(shape)
    masks = binary_masks(shape)
    full = (1 << n) - 1
    edges = {}
    representatives = {}
    witnesses = {}

    def record(a, b):
        c = full ^ a ^ b
        ga, gb, gc = masks[a], masks[b], masks[c]
        g = ga & gb
        assert g == ga & gc == gb & gc
        if g not in edges:
            edges[g] = set()
            representatives[g] = (a, b)
        for target, remaining in [(ga, a), (gb, b), (gc, c)]:
            if target != g:
                assert g & target == g
                edges[g].add(target)
                witnesses.setdefault((g, target), (a, b, remaining))

    # Restricted-growth partitions; first site belongs to symbol class zero.
    def visit(i, a, b, used):
        if i == n:
            record(a, b)
            return
        bit = 1 << i
        visit(i+1, a | bit, b, used)
        visit(i+1, a, b | bit, max(used, 2))
        if used >= 2:
            visit(i+1, a, b, 3)

    visit(1, 1, 0, 1)
    binary = set(masks)
    binary_words = {}
    for word, g in enumerate(masks):
        binary_words.setdefault(g, word)
    assert binary <= set(edges)
    reach = {g: set(targets) for g, targets in edges.items()}
    for g in sorted(reach, key=int.bit_count, reverse=True):
        for target in edges[g]:
            reach[g].update(reach[target])

    # Audit every stored edge against the original, independent comparator.
    for (g, target), (a, b, remaining) in witnesses.items():
        word = tuple(0 if a >> i & 1 else 1 if b >> i & 1 else 2
                     for i in range(n))
        merged = tuple((remaining >> i) & 1 for i in range(n))
        assert agreement(shape, word) == g
        assert agreement(shape, merged) == target

    pairs = [(a, b) for a in edges for b in edges if a != b and a & b == a]
    missing = [(a, b) for a, b in pairs if b not in reach[a]]
    binary_pairs = [(a, b) for a, b in pairs if a in binary and b in binary]
    binary_missing = [(a, b) for a, b in missing if a in binary and b in binary]
    direct_missing = sum(b not in edges[a] for a, b in binary_pairs)

    def binary_word_for(g):
        return [(binary_words[g] >> i) & 1 for i in range(n)]

    return dict(shape=shape, binary_classes=len(binary), ternary_classes=len(edges),
                merger_edges=len(witnesses), inclusion_pairs=len(pairs),
                unreached_inclusions=len(missing),
                binary_inclusion_pairs=len(binary_pairs),
                binary_inclusions_without_direct_edge=direct_missing,
                unreached_binary_inclusions=len(binary_missing),
                unreached_binary_examples=[dict(mask_u=a, mask_v=b,
                    u=binary_word_for(a), v=binary_word_for(b))
                    for a, b in binary_missing[:3]])


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--large', action='store_true', help='also enumerate 3x5 and 4x4')
    args = p.parse_args()
    shapes = [(2, 2), (2, 3), (2, 4), (3, 3), (2, 5), (3, 4), (2, 6)]
    if args.large:
        shapes += [(3, 5), (4, 4)]
    results = []
    for h, w in shapes:
        result = scan(rectangle(h, w))
        results.append(result)
        print((h, w), 'classes', result['ternary_classes'], 'unreached binary',
              result['unreached_binary_inclusions'], 'unreached total',
              result['unreached_inclusions'], flush=True)
    obstructions = [((0, 0), (0, 1), (1, 0), (1, 1), (1, 2), (1, 3)),
                    tuple((0, j) for j in (0, 1, 2, 5, 7))]
    for support in obstructions:
        results.append(scan(support))
        assert results[-1]['unreached_binary_inclusions'] > 0
    output = dict(status='Exact finite pattern-class enumeration; paths give '
                        'all-anchor comparisons using the merging theorem. '
                        'Unreached pairs do not refute entropy or convex order.',
                  records=results)
    dest = Path(__file__).with_name('merger-evidence.json')
    dest.write_text(json.dumps(output, indent=2) + '\n')


if __name__ == '__main__':
    main()
