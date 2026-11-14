#!/usr/bin/env python3
"""Rational finite ARW stabilization and negative-dependence diagnostics.

States use -1 for a sleeper, 0 for empty, and k>0 for k active particles.
The sink is indexed by n. Failed sleep instructions are skipped at sites
with two or more particles, which does not change absorption probabilities.
All arithmetic and residual checks are exact. Intended for purposeful small
examples, not an uncontrolled graph sweep; see report.md.

Initial implementation developed with ChatGPT by OpenAI from Jacob Richey's
research direction. Scientific validation is documented in report.md.
"""

from collections import defaultdict
from functools import lru_cache
from itertools import combinations

import sympy as sp
from sympy.polys.matrices import DomainMatrix


def graph_kernel(neighbors):
    """Uniform choice from each row; n denotes the sink."""
    return tuple(tuple((y, sp.Rational(1, len(row))) for y in row)
                 for row in neighbors)


def stable(mask, n):
    return tuple(-1 if mask >> i & 1 else 0 for i in range(n))


def add(state, x):
    out = list(state)
    out[x] = 2 if out[x] == -1 else out[x] + 1
    return tuple(out)


def stable_mask(state):
    return sum(1 << i for i, x in enumerate(state) if x == -1)


def transitions(state, kernel, sleep, reverse=False):
    n = len(state)
    order = range(n - 1, -1, -1) if reverse else range(n)
    x = next(i for i in order if state[i] > 0)
    out = defaultdict(lambda: sp.S.Zero)
    if state[x] == 1:
        asleep = list(state)
        asleep[x] = -1
        out[tuple(asleep)] += sleep[x]
        jump = 1 - sleep[x]
    else:
        jump = sp.S.One
    for y, prob in kernel[x]:
        nxt = list(state)
        nxt[x] -= 1
        if y < n:
            nxt[y] = 2 if nxt[y] == -1 else nxt[y] + 1
        out[tuple(nxt)] += jump * prob
    assert sp.expand(sum(out.values()) - 1) == 0
    return {s: p for s, p in out.items() if p != 0}


def absorption(starts, kernel, sleep, reverse=False):
    """Return absorption matrix for listed starts, checking full equations."""
    n = len(kernel)
    sleep = tuple(map(sp.Rational, sleep))
    for row in kernel:
        assert sum(p for _, p in row) == 1
        assert all(0 <= y <= n and p >= 0 for y, p in row)
    assert all(0 < p < 1 for p in sleep)
    transient = []
    index = {}

    def discover(s):
        if any(v > 0 for v in s) and s not in index:
            index[s] = len(transient)
            transient.append(s)

    for s in starts:
        discover(s)
    rows = []
    for state in transient:
        row = transitions(state, kernel, sleep, reverse)
        rows.append(row)
        for nxt in row:
            discover(nxt)
    m = len(transient)
    left = sp.MutableSparseMatrix.eye(m)
    right = sp.MutableSparseMatrix.zeros(m, 1 << n)
    for i, row in enumerate(rows):
        for nxt, prob in row.items():
            if nxt in index:
                left[i, index[nxt]] -= prob
            else:
                right[i, stable_mask(nxt)] += prob
    a = DomainMatrix.from_Matrix(left).to_field()
    b = DomainMatrix.from_Matrix(right).convert_to(a.domain)
    sol = a.inv().matmul(b).to_Matrix()
    assert left * sol == right
    assert all(sum(sol.row(i)) == 1 for i in range(m))
    assert all(v >= 0 for v in sol)
    result = []
    for state in starts:
        if state in index:
            result.append(list(sol.row(index[state])) )
        else:
            result.append([sp.S.One if mask == stable_mask(state) else sp.S.Zero
                           for mask in range(1 << n)])
    return sp.Matrix(result), m


def stationary_law(kernel, sleep, verify_additions=True):
    n = len(kernel)
    starts = [(1,) * n]
    if verify_additions:
        starts += [add(stable(mask, n), x)
                   for x in range(n) for mask in range(1 << n)]
    laws, m = absorption(starts, kernel, sleep)
    mu = laws[0, :]
    if verify_additions:
        for x in range(n):
            transition = laws[1 + x * (1 << n):1 + (x + 1) * (1 << n), :]
            assert mu * transition == mu
    # An independently chosen legal toppling order must agree by abelianity.
    other, other_m = absorption([(1,) * n], kernel, sleep, reverse=True)
    assert other == mu
    return tuple(mu), (m, other_m)


@lru_cache(None)
def increasing_events(coords):
    """All nonconstant increasing Boolean functions, encoded by true masks."""
    k = len(coords)
    events = []
    for bits in range(1, (1 << (1 << k)) - 1):
        if all(not (bits >> mask & 1) or all(
                bits >> (mask | (1 << j)) & 1 for j in range(k))
               for mask in range(1 << k)):
            events.append(tuple(mask for mask in range(1 << k) if bits >> mask & 1))
    return tuple(events)


def na_gaps(mu):
    """Yield all disjoint increasing-event pairs and E[f]E[g]-E[fg]."""
    n = len(mu).bit_length() - 1
    if n > 5:
        raise ValueError("Exhaustive event enumeration is limited to at most five sites")
    for a_mask in range(1, 1 << n):
        a = tuple(i for i in range(n) if a_mask >> i & 1)
        for b_mask in range(a_mask + 1, 1 << n):
            if a_mask & b_mask:
                continue
            b = tuple(i for i in range(n) if b_mask >> i & 1)
            for ea in increasing_events(a):
                ia = [sum(((s >> i) & 1) << j for j, i in enumerate(a)) in ea
                      for s in range(1 << n)]
                pa = sum(p for p, keep in zip(mu, ia) if keep)
                for eb in increasing_events(b):
                    ib = [sum(((s >> i) & 1) << j for j, i in enumerate(b)) in eb
                          for s in range(1 << n)]
                    pb = sum(p for p, keep in zip(mu, ib) if keep)
                    pab = sum(p for p, x, y in zip(mu, ia, ib) if x and y)
                    yield (a, ea, b, eb), pa * pb - pab


def na_check(mu):
    """Exhaustive increasing-event test, sufficient for NA on this finite cube."""
    tests = 0
    largest = None
    witness = None
    for pair, gap in na_gaps(mu):
        tests += 1
        if largest is None or -gap > largest:
            largest = -gap
            witness = pair
    return tests, largest, witness


def polynomial(mu):
    n = len(mu).bit_length() - 1
    z = sp.symbols(f"z0:{n}")
    return z, sum(p * sp.prod(z[i] for i in range(n) if s >> i & 1)
                  for s, p in enumerate(mu))


def rayleigh_differences(mu):
    z, poly = polynomial(mu)
    return {(i, j): sp.factor(sp.diff(poly, z[i]) * sp.diff(poly, z[j])
                             - poly * sp.diff(poly, z[i], z[j]))
            for i, j in combinations(range(len(z)), 2)}


def main():
    cases = {
        "dag3": ((1,), (2,), (3,)),
        "path2": ((1, 2), (0, 2)),
        "path3": ((1, 3), (0, 2), (1, 3)),
        "triangle3": ((1, 2, 3), (0, 2, 3), (0, 1, 3)),
        "path4": ((1, 4), (0, 2), (1, 3), (2, 4)),
        "star4": ((1, 2, 3), (0, 4), (0, 4), (0, 4)),
        "complete4": ((1, 2, 3, 4), (0, 2, 3, 4), (0, 1, 3, 4), (0, 1, 2, 4)),
    }
    for name, neighbors in cases.items():
        n = len(neighbors)
        mu, sizes = stationary_law(graph_kernel(neighbors), (sp.Rational(1, 2),) * n)
        comparison = na_check(mu)
        assert comparison[1] <= 0
        differences = rayleigh_differences(mu)
        z, _ = polynomial(mu)
        for pair, difference in differences.items():
            remaining = [z[i] for i in range(n) if i not in pair]
            coefficients = sp.Poly(difference, *remaining).coeffs() if remaining else [difference]
            assert all(c >= 0 for c in coefficients)
        if name == "dag3":
            assert mu == (sp.Rational(1, 8),) * 8
        print(name, "states", sizes, "mu", mu, flush=True)
        print("NA", comparison, flush=True)
        print("Rayleigh differences", differences, flush=True)


if __name__ == "__main__":
    main()
