"""One valuation example and small fixed combinatorial checks; no search.

AI-written research code, 2026-10-19. Standard library only. Run with
python3 -B. No output files. Uses RN (115)--(117) for an exact rational
identity and independently checks raw Taylor coefficients and one
rational parameter point under both toppling priorities. A separate
small combinatorial audit checks source insertion, trace, contraction,
circuit leak locations, and the distinction between rank and guaranteed
IDLA occupancy. One further three-site graph checks possible filling,
exact support, and minimum occupancy by reachability matching. Three
particles on the fixed two-arm tree test output-mark exchange and the
failure of the Lorentzian Hessian condition.
Finite checks do not replace the general proofs.
"""
from collections import defaultdict
from fractions import Fraction as F
from itertools import product
from math import prod
import sys
sys.dont_write_bytecode = True
import verify_response_renewal as rn

ZERO = (0, 0, 0)
ONE = {ZERO: F(1)}


def add(*polys):
    out = defaultdict(F)
    for poly in polys:
        for mon, coefficient in poly.items():
            out[mon] += coefficient
    return {mon: coefficient for mon, coefficient in out.items() if coefficient}


def scale(poly, coefficient):
    return {mon: value*coefficient for mon, value in poly.items()
            if value*coefficient}


def mul(*polys):
    out = ONE
    for poly in polys:
        new = defaultdict(F)
        for left, a in out.items():
            for right, b in poly.items():
                new[tuple(x+y for x, y in zip(left, right))] += a*b
        out = {mon: value for mon, value in new.items() if value}
    return out


def divide_monomial(poly, powers):
    assert all(all(a >= b for a, b in zip(mon, powers)) for mon in poly)
    return {tuple(a-b for a, b in zip(mon, powers)): coefficient
            for mon, coefficient in poly.items()}


def value(poly, point):
    return sum(coefficient*prod(x**power for x, power in zip(point, mon))
               for mon, coefficient in poly.items())


def minimal_support(poly):
    return {mon: coefficient for mon, coefficient in poly.items()
            if not any(other != mon and all(a <= b for a, b in zip(other, mon))
                       for other in poly)}


def determinant(matrix):
    if len(matrix) == 1:
        return matrix[0][0]
    return add(*(scale(mul(matrix[0][j], determinant(
        [[row[k] for k in range(len(matrix)) if k != j]
         for row in matrix[1:]])), (-1)**j) for j in range(len(matrix))))


def adjugate(matrix):
    n = len(matrix)
    return [[scale(determinant([[matrix[a][b] for b in range(n) if b != i]
                                for a in range(n) if a != j]), (-1)**(i+j))
             for j in range(n)] for i in range(n)]


def cycle_rational_gain():
    """Target 2; lower {0,2}; upper {0,1,2}; return/leak c=d=1/2.

    M=diag(p)(I-Rcal)diag(1/p) in RN (115). Pair masses beta,
    gamma obey M diag(p)(beta-gamma)=P*(d*t/L-e_1).
    The three killed-return differences in bracket come from RN (117).
    """
    p = [{tuple(int(i == j) for i in range(3)): F(1)} for j in range(3)]
    q = [add(ONE, scale(a, -1)) for a in p]
    c = d = F(1, 2)
    P = mul(*p)
    matrix = [[ONE, scale(q[0], -1), scale(mul(q[0], p[1]), -1)],
              [scale(mul(q[1], p[2]), -c), ONE, scale(q[1], -1)],
              [scale(q[2], -c), scale(mul(q[2], p[0]), -c), ONE]]
    det, adj = determinant(matrix), adjugate(matrix)
    for i in range(3):
        for j in range(3):
            assert add(*(mul(matrix[i][k], adj[k][j]) for k in range(3))) == (
                det if i == j else {})
    L = add(scale(ONE, d), scale(P, c))
    t = [mul(q[0], p[1], p[2]), mul(q[1], p[2]), q[2]]
    T = [add(scale(add(*(mul(adj[i][j], t[j]) for j in range(3))), d),
             scale(mul(L, adj[i][1]), -1)) for i in range(3)]
    bracket = add(
        mul(add(p[0], p[1], scale(mul(p[0], p[1]), -2)),
            divide_monomial(T[1], (0, 1, 0))),
        mul(q[0], p[1], divide_monomial(T[0], (1, 0, 0))),
        mul(add(p[0], p[1], scale(mul(p[0], p[1]), -1)), det))
    numerator = scale(mul(q[2], P, bracket), c)
    denominator = mul(add(ONE, scale(mul(*q), -c)), L, det)
    # A compact numerator identity, with s=p0+p1 and t=p0*p1.
    s, t = add(p[0], p[1]), mul(p[0], p[1])
    B = add(s, scale(mul(s, s), -1), mul(add(ONE, p[2]), t),
            mul(add(scale(ONE, 4), scale(p[2], -1)), s, t),
            scale(mul(add(scale(ONE, 3), scale(p[2], -1)), t, t), -1))
    assert numerator == scale(mul(P, q[2], B), F(1, 8))
    assert denominator[ZERO] == F(1, 8)
    assert minimal_support(numerator) == {(2, 1, 1): F(1, 8), (1, 2, 1): F(1, 8)}
    return numerator, denominator


KERNEL = [[F(0), F(1), F(0), F(0)],
          [F(0), F(0), F(1), F(0)],
          [F(1, 2), F(0), F(0), F(1, 2)]]
LOWER, UPPER = (1, 0, 1), (1, 1, 1)


def raw_taylor(degree=5, terminal_support=None):
    """Differentiate the maintained raw absorption equations at zero."""
    if terminal_support is not None:
        assert 2 in terminal_support
    def reward(end):
        if terminal_support is None:
            return end[2] == -1
        return frozenset(i for i, a in enumerate(end) if a == -1) == terminal_support
    states = [LOWER, UPPER]
    pos = {state: i for i, state in enumerate(states)}
    for state in states:
        for end, weight in rn.arw_step(KERNEL, [F(1, 2)]*3, state).items():
            if weight and any(a > 0 for a in end) and end not in pos:
                pos[end] = len(states)
                states.append(end)
    rows, changes, exits = [], [], []
    for i, state in enumerate(states):
        z = next(j for j, a in enumerate(state) if a > 0)
        base = rn.arw_step(KERNEL, [F(0)]*3, state)
        unit = [F(0)]*3
        unit[z] = F(1)
        at_one = rn.arw_step(KERNEL, unit, state)
        derivative = {end: at_one.get(end, F(0))-base.get(end, F(0))
                      for end in base.keys() | at_one.keys()}
        row = {i: F(1)}
        for end, weight in base.items():
            if end in pos:
                row[pos[end]] = row.get(pos[end], F(0))-weight
            else:
                assert not reward(end) or weight == 0
        rows.append(row)
        changes.append((z, {pos[end]: weight for end, weight in derivative.items()
                            if end in pos and weight}))
        exits.append(sum(weight for end, weight in derivative.items()
                         if end not in pos and reward(end)))
    coefficients, gain = {}, {}
    monomials = sorted((mon for mon in product(range(degree+1), repeat=3)
                        if sum(mon) <= degree), key=lambda mon: (sum(mon), mon))
    for mon in monomials:
        rhs = []
        for i, (z, change) in enumerate(changes):
            if not mon[z]:
                rhs.append([F(0)])
                continue
            predecessor = tuple(a-int(j == z) for j, a in enumerate(mon))
            previous = coefficients[predecessor]
            rhs.append([sum(weight*previous[j] for j, weight in change.items())
                        + (exits[i] if sum(mon) == 1 else 0)])
        solution = [row[0] for row in rn.old.solve_rows(rows, rhs)]
        coefficients[mon] = solution
        gain[mon] = solution[pos[UPPER]]-solution[pos[LOWER]]
    return gain, len(states)


def numeric_audit(numerator, denominator):
    p = [F(1, 3), F(2, 5), F(3, 7)]
    expected = value(numerator, p)/value(denominator, p)
    for reverse in (False, True):
        indices = list(reversed(range(3))) if reverse else list(range(3))
        kernel = [[KERNEL[i][j] for j in indices]+[KERNEL[i][3]] for i in indices]
        lower, upper = (tuple(state[i] for i in indices) for state in (LOWER, UPPER))
        target, pp = indices.index(2), [p[i] for i in indices]
        laws, _ = rn.stopped([lower, upper], lambda state: rn.arw_step(kernel, pp, state),
                             lambda state: all(a <= 0 for a in state))
        gain = sum(weight for end, weight in laws[upper].items() if end[target] == -1)
        gain -= sum(weight for end, weight in laws[lower].items() if end[target] == -1)
        assert gain == expected
    return expected


def subsets(ground):
    ground = tuple(sorted(ground))
    return [frozenset(i for i, bit in zip(ground, bits) if bit)
            for bits in product((0, 1), repeat=len(ground))]


def dominator_regions(kernel):
    n = len(kernel)
    def reaches_sink(start, forbidden):
        if start == forbidden:
            return False
        seen, queue = {start}, [start]
        for z in queue:
            if z == n:
                return True
            for y, weight in enumerate(kernel[z]):
                if weight and y != forbidden and y not in seen:
                    seen.add(y)
                    queue.append(y)
        return False
    return {d: frozenset(a for a in range(n) if not reaches_sink(a, d))
            for d in range(n)}


def independent_sets(regions, eta, ground=None):
    if ground is None:
        ground = range(len(eta))
    return {C for C in subsets(ground)
            if all(len(C & B) <= sum(eta[a] for a in B) for B in regions.values())}


def circuits(independent, ground):
    return {C for C in subsets(ground) if C not in independent
            and all(C-{a} in independent for a in C)}


def trace_kernel(kernel, retained):
    retained = tuple(sorted(retained))
    n = len(kernel)
    laws, _ = rn.stopped(range(n+1),
        lambda z: {y: weight for y, weight in enumerate(kernel[z]) if weight},
        lambda z: z == n or z in retained)
    trace = [[sum(weight*laws[y].get(z, F(0)) for y, weight in enumerate(kernel[x]))
              for z in retained+(n,)] for x in retained]
    assert all(sum(row) == 1 for row in trace)
    return trace


def idla_law(kernel, capacity, sources, initial=frozenset()):
    """Exact sequential IDLA, using killed-walk absorption for each particle."""
    n = len(kernel)
    law = {frozenset(initial): F(1)}
    for source in sources:
        new = defaultdict(F)
        for occupied, probability in law.items():
            vacancies = capacity-occupied
            hits, _ = rn.stopped([source],
                lambda z: {y: weight for y, weight in enumerate(kernel[z]) if weight},
                lambda z: z == n or z in vacancies)
            for end, weight in hits[source].items():
                new[occupied if end == n else occupied | {end}] += probability*weight
        law = dict(new)
    assert sum(law.values()) == 1
    return law


def matroid_neighborhood_audit():
    # A directed path: one source insertion replaces all pair circuits
    # by one triple circuit. This checks the circuit-ideal intersection.
    path = [[F(0), F(1), F(0), F(0)],
            [F(0), F(0), F(1), F(0)],
            [F(0), F(0), F(0), F(1)]]
    regions = dominator_regions(path)
    eta = (1, 0, 0)
    lower = independent_sets(regions, eta)
    upper = independent_sets(regions, (1, 1, 0))
    D = {d for d, B in regions.items() if 1 in B}
    assert upper == lower | {I | {d} for I in lower for d in D}
    assert circuits(lower, range(3)) == {frozenset(C) for C in ((0, 1), (0, 2), (1, 2))}
    assert circuits(upper, range(3)) == {frozenset(range(3))}
    for C in subsets(range(3)):
        assert (C not in upper) == all(C-{d} not in lower for d in D)
    W = (0, 2)
    traced = trace_kernel(path, W)
    traced_independent = independent_sets(dominator_regions(traced), (1, 0))
    assert {frozenset(W[i] for i in C) for C in traced_independent} == {
        C for C in lower if C <= set(W)}

    # A cyclic fork with an incomparable fourth vertex. Its nontrivial
    # circuit is {a,b,r}; trace has one leak at r, or two with v adjoined.
    fork = [[F(0), F(0), F(1), F(0), F(0)],
            [F(0), F(0), F(1), F(0), F(0)],
            [F(1, 4), F(1, 4), F(0), F(1, 4), F(1, 4)],
            [F(0), F(0), F(1, 2), F(0), F(1, 2)]]
    regions = dominator_regions(fork)
    eta, R = (1, 1, 0, 0), frozenset({2})
    independent = independent_sets(regions, eta)
    cs = circuits(independent, range(4))
    assert cs == {frozenset({3}), frozenset({0, 1, 2})}
    for C in cs:
        violated = [d for d, B in regions.items() if len(C & B) > sum(eta[a] for a in B)]
        d = min(violated, key=lambda z: len(regions[z]))
        assert d in C and eta[d] == 0 and C <= regions[d]
        children = [e for e, B in regions.items() if B < regions[d]
                    and not any(B < BB < regions[d] for BB in regions.values())]
        assert all(len(C & regions[e]) == sum(eta[a] for a in regions[e]) for e in children)
        Kc = trace_kernel(fork, C)
        assert {z for z, row in zip(sorted(C), Kc) if row[-1]} == {d}
        for x in (0, 3):
            added = tuple(a+int(i == x) for i, a in enumerate(eta))
            assert (C in independent_sets(regions, added)) == (x in regions[d])
    assert {i for i, row in enumerate(fork) if row[-1]} == {2, 3}
    Kc = trace_kernel(fork, {0, 1, 2})
    assert Kc[2][-1] == F(3, 8)
    with_root = independent_sets(regions, (1, 1, 1, 0))
    assert {C-R for C in with_root if R <= C} == {C for C in independent if not C & R}
    contracted = {C-R for C in independent if R <= C}
    assert contracted == {frozenset(), frozenset({0}), frozenset({1})}
    assert all(contracted != independent_sets(regions, residual, {0, 1, 3})
               for residual in ((1, 0, 0, 0), (0, 1, 0, 0)))

    # Two alternative length-two branches: identical sure-filling
    # matroids for k=1,2,3 at x, but minimum occupancies 1,2,2.
    branches = [[F(0), F(1, 2), F(0), F(1, 2), F(0), F(0)],
                [F(0), F(0), F(1), F(0), F(0), F(0)],
                [F(0), F(0), F(0), F(0), F(0), F(1)],
                [F(0), F(0), F(0), F(0), F(1), F(0)],
                [F(0), F(0), F(0), F(0), F(0), F(1)]]
    regions = dominator_regions(branches)
    C = frozenset({1, 2, 3, 4})
    count_minima = []
    for k in (1, 2, 3):
        assert independent_sets(regions, (k, 0, 0, 0, 0)) == {frozenset(), frozenset({0})}
        law = idla_law(branches, C, (0,)*k)
        count_minima.append(min(map(len, law)))
        if k == 2:
            assert law == {frozenset({1, 2}): F(1, 4), frozenset({1, 3}): F(1, 2),
                           frozenset({3, 4}): F(1, 4)}
    assert count_minima == [1, 2, 2]
    traced = trace_kernel(branches, C)
    tr = dominator_regions(traced)
    entrance_intersection = independent_sets(tr, (1, 0, 0, 0)) & independent_sets(tr, (0, 0, 1, 0))
    assert entrance_intersection == {frozenset()}
    return dict(fixed_graphs=3, branch_minimum_occupancies=count_minima,
                fork_trace_leak=str(Kc[2][-1]), source_insertion='passed',
                restriction_and_contraction='passed', circuit_anchors='passed')


def reachable_avoiding(kernel, source, forbidden):
    """Reachability includes time zero; a deleted source has no neighbors."""
    if source in forbidden:
        return frozenset()
    sink = len(kernel)
    seen, queue = {source}, [source]
    for z in queue:
        if z == sink:
            continue
        for y, weight in enumerate(kernel[z]):
            if weight and y not in forbidden and y not in seen:
                seen.add(y)
                queue.append(y)
    return frozenset(seen)


def has_left_saturating_matching(neighbors):
    """Ordinary bipartite augmenting paths; right slots are integers."""
    owner = {}
    def augment(i, seen):
        for slot in sorted(neighbors[i]):
            if slot in seen:
                continue
            seen.add(slot)
            if slot not in owner or augment(owner[slot], seen):
                owner[slot] = i
                return True
        return False
    return all(augment(i, set()) for i in range(len(neighbors)))


def exact_support_matching(kernel, capacity, sources, occupied):
    if not occupied <= capacity or len(occupied) > len(sources):
        return False
    sink = len(kernel)
    sink_slots = set(range(sink, sink+len(sources)-len(occupied)))
    reaches = [reachable_avoiding(kernel, a, capacity-occupied) for a in sources]
    neighbors = [(R & occupied) | (sink_slots if sink in R else set())
                 for R in reaches]
    return has_left_saturating_matching(neighbors)


def contained_support_matching(kernel, capacity, sources, allowed):
    """An output contained in allowed needs only sink-blocked tokens matched."""
    assert allowed <= capacity
    sink = len(kernel)
    reaches = [reachable_avoiding(kernel, a, capacity-allowed) for a in sources]
    return has_left_saturating_matching([R & allowed for R in reaches if sink not in R])


def possible_support_audit():
    # One graph only: x->a, a->b or sink equally, b->x. Two tokens at x.
    # Its cycle tests genuine absorption; x in a forbidden set fails at time zero.
    kernel = [[F(0), F(1), F(0), F(0)],
              [F(0), F(0), F(1, 2), F(1, 2)],
              [F(1), F(0), F(0), F(0)]]
    sources = (0, 0)
    reaches = [reachable_avoiding(kernel, a, frozenset()) for a in sources]
    comparisons = 0
    for C in subsets(range(3)):
        law = idla_law(kernel, C, sources)
        supports = {B for B, weight in law.items() if weight > 0}
        full_matching = has_left_saturating_matching(
            [{i for i, R in enumerate(reaches) if c in R} for c in sorted(C)])
        assert (C in supports) == full_matching == (len(C) <= 2)
        candidates = []
        for B in subsets(C):
            assert (B in supports) == exact_support_matching(kernel, C, sources, B)
            contained = contained_support_matching(kernel, C, sources, B)
            assert contained == any(A <= B for A in supports)
            if contained:
                candidates.append(len(B))
            comparisons += 1
        assert min(candidates) == min(map(len, supports))
    law = idla_law(kernel, frozenset({1, 2}), sources)
    assert law == {frozenset({1}): F(1, 2), frozenset({1, 2}): F(1, 2)}
    return dict(fixed_graphs=1, capacity_support_pairs=comparisons,
                full_filling='passed', exact_support='passed',
                minimum_occupancy='passed', forbidden_source='passed')


def output_mark_audit():
    # Exactly the requested tree: x->a or c; a->b->s; c->d->s.
    kernel = [[F(0), F(1, 2), F(0), F(1, 2), F(0), F(0)],
              [F(0), F(0), F(1), F(0), F(0), F(0)],
              [F(0), F(0), F(0), F(0), F(0), F(1)],
              [F(0), F(0), F(0), F(0), F(1), F(0)],
              [F(0), F(0), F(0), F(0), F(0), F(1)]]
    law = idla_law(kernel, frozenset(range(5)), (0, 0, 0))
    expected = {frozenset({0, 1, 2}): F(1, 4),
                frozenset({0, 1, 3}): F(1, 2),
                frozenset({0, 3, 4}): F(1, 4)}
    assert law == expected and all(len(B) == 3 for B in law)
    left, right = frozenset({0, 1, 2}), frozenset({0, 3, 4})
    assert all((left-{1}) | {j} not in law for j in right-left)

    # Rooted feasible sets consist of x and a prefix of each arm.
    feasible = {frozenset()} | {
        frozenset({0}) | frozenset((1, 2)[:i]) | frozenset((3, 4)[:j])
        for i in range(3) for j in range(3)}
    assert set(law) == {B for B in feasible if len(B) == 3}

    # Hessian of d/dz_x of the output polynomial, in coordinates a,b,c,d.
    hessian = [[F(0) for _ in range(4)] for _ in range(4)]
    for B, weight in law.items():
        i, j = sorted(B-{0})
        hessian[i-1][j-1] += weight
        hessian[j-1][i-1] += weight
    directions = ((1, 1, 0, 0), (0, -2, 1, 1))
    restriction = [[sum(u[i]*hessian[i][j]*v[j]
                        for i in range(4) for j in range(4))
                    for v in directions] for u in directions]
    assert restriction == [[F(1, 2), F(0)], [F(0), F(1, 2)]]
    return dict(particles=3, nonlosing_probability='1',
                terminal_law={str(tuple(sorted(B))): str(p) for B, p in law.items()},
                basis_exchange='fails as claimed', antimatroid_layer='passed',
                hessian_positive_subspace_dimension=2)


def main():
    numerator, denominator = cycle_rational_gain()
    raw, states = raw_taylor()
    for mon in raw:
        product_coefficient = sum(
            weight*raw.get(tuple(a-b for a, b in zip(mon, term)), F(0))
            for term, weight in denominator.items()
            if all(a >= b for a, b in zip(mon, term)))
        assert product_coefficient == numerator.get(mon, F(0)), mon
    assert {mon: weight for mon, weight in raw.items() if weight and sum(mon) <= 4} == {
        (2, 1, 1): F(1), (1, 2, 1): F(1)}
    assert raw[3, 1, 1] == -2  # Leading-face positivity is not coefficient positivity.
    print('VERIFIED exact numerator minimal exponents:', minimal_support(numerator))
    print('VERIFIED raw Taylor coefficients:', len(raw), 'on', states, 'transient states')
    print('VERIFIED rational point under both priorities:', numeric_audit(numerator, denominator))
    print('VERIFIED signed higher coefficient [p0^3 p1 p2] Delta =', raw[3, 1, 1])
    support_coefficients = {}
    for support in ((2,), (0, 2), (1, 2), (0, 1, 2)):
        table, _ = raw_taylor(3, frozenset(support))
        support_coefficients[support] = table[1, 1, 1]
    assert support_coefficients == {(2,): -1, (0, 2): -1, (1, 2): 0, (0, 1, 2): 2}
    print('VERIFIED squarefree cubic cancellation on supports containing 2:', support_coefficients)
    print('VERIFIED matroid neighborhood:', matroid_neighborhood_audit())
    print('VERIFIED possible IDLA support:', possible_support_audit())
    print('VERIFIED output-mark obstruction:', output_mark_audit())


if __name__ == '__main__':
    main()
