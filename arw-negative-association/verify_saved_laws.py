"""Recheck all event-gap polynomials of the 44 exported four-site laws.

This validates normalization, coefficientwise mass positivity, graph coverage,
and all disjoint increasing-event inequalities. It does not recompute the
absorption solutions; four_site.py provides that separate, longer calculation.
"""
from pathlib import Path
import json
from itertools import combinations

from exact import sp, increasing_events
from four_site import graph_types


def main():
    data = json.loads(Path(__file__).with_name('four_site_certificates.json').read_text())
    graphs = {str(key): neighbors for key, neighbors in graph_types()}
    if set(graphs) != set(data['cases']) or len(graphs) != 44:
        raise ValueError('Four-site graph coverage does not match')
    t = sp.Symbol('t')
    total = 0
    for key, saved in data['cases'].items():
        if tuple(map(tuple, saved['neighbors'])) != graphs[key]:
            raise ValueError('Graph does not match its label: ' + key)
        den = sp.Poly(sp.sympify(saved['denominator']), t)
        weights = [sp.Poly(sp.sympify(x), t) for x in saved['numerators']]
        zero = sp.Poly(0, t)
        if len(weights) != 16 or sum(weights, zero) != den:
            raise ValueError('Masses do not normalize: ' + key)
        if any(any(c < 0 for c in p.all_coeffs()) or p.is_zero for p in [den, *weights]):
            raise ValueError('Mass positivity check failed: ' + key)
        count = 0
        for am, bm in combinations(range(1, 16), 2):
            if am & bm:
                continue
            aa = tuple(i for i in range(4) if am >> i & 1)
            bb = tuple(i for i in range(4) if bm >> i & 1)
            for ea in increasing_events(aa):
                ia = [sum(((s >> i) & 1) << j for j, i in enumerate(aa)) in ea for s in range(16)]
                pa = sum((w for w, keep in zip(weights, ia) if keep), zero)
                for eb in increasing_events(bb):
                    ib = [sum(((s >> i) & 1) << j for j, i in enumerate(bb)) in eb for s in range(16)]
                    pb = sum((w for w, keep in zip(weights, ib) if keep), zero)
                    joint = sum((w for w, a, b in zip(weights, ia, ib) if a and b), zero)
                    if any(c < 0 for c in (pa * pb - den * joint).all_coeffs()):
                        raise ValueError('Negative event-gap coefficient: ' + key)
                    count += 1
        if count != 174:
            raise ValueError('Wrong event count: ' + key)
        total += count
    print(f'PASS: {len(graphs)} normalized positive laws; {total} event-gap polynomials')
    print('Scope: exported laws checked; absorption solutions not recomputed by this command.')


if __name__ == '__main__':
    main()
