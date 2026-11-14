"""Focused exact addition-algebra, polytope, and duality certificates.

Standard library only. AI-written exploratory certificate, not a sweep.
Reuses a single-leak cycle, the asymmetric two-site current kernel, and
a two-site path. Checks algebra/Jordan formulas, volume normalizations,
Poisson cell integrals, complementary-sleep reversal links, and ordered
bridge books, including one independent raw release. General arguments
belong in response-renewal.md.
"""

from fractions import Fraction as F
from itertools import product
from collections import defaultdict
import json

from verify_creator_potentials import idla_addition_matrices
from verify_response_renewal import multiply, inverse


def rank(matrix):
    a = [list(row) for row in matrix]
    pivot_row = 0
    for col in range(len(a[0])):
        pivot = next((i for i in range(pivot_row, len(a)) if a[i][col]), None)
        if pivot is None:
            continue
        a[pivot_row], a[pivot] = a[pivot], a[pivot_row]
        divisor = a[pivot_row][col]
        a[pivot_row] = [x / divisor for x in a[pivot_row]]
        for i in range(len(a)):
            if i != pivot_row and a[i][col]:
                factor = a[i][col]
                a[i] = [x - factor * y for x, y in zip(a[i], a[pivot_row])]
        pivot_row += 1
        if pivot_row == len(a):
            break
    return pivot_row


def supports(n):
    return tuple(frozenset(i for i, b in enumerate(bits) if b)
                 for bits in product((0, 1), repeat=n))


def sink_core(K, B):
    n = len(K)
    reachable = {i for i in B if K[i][n]}
    while True:
        new = {i for i in B if any(K[i][j] for j in reachable)} | reachable
        if new == reachable:
            return frozenset(new)
        reachable = new


def determinant(matrix):
    if not matrix:
        return F(1)
    return sum((-1)**j * entry * determinant(
        [[value for k, value in enumerate(row) if k != j]
         for row in matrix[1:]]) for j, entry in enumerate(matrix[0]))


def homogeneous(K):
    sets = supports(len(K))
    additions = idla_addition_matrices(K, sets)
    return sets, [[[value if i != j else F(0)
                    for j, value in enumerate(row)] for i, row in enumerate(a)]
                 for a in additions]


def top_product(sets, factors):
    row = [F(not C) for C in sets]
    for factor in factors:
        row = multiply([row], factor)[0]
    return row[-1]


def capacity_additions(K, capacities):
    """Killed-at-first-loss IDLA matrices, retaining counts up to capacity."""
    states = tuple(product(*(range(d+1) for d in capacities)))
    pos = {state: i for i, state in enumerate(states)}
    n, dim = len(K), len(states)
    additions = [[[F(0)]*dim for _ in states] for _ in range(n)]
    for index, state in enumerate(states):
        full = [i for i in range(n) if state[i] == capacities[i]]
        green = inverse([[F(i == j)-K[i][j] for j in full] for i in full])
        for x in range(n):
            if x not in full:
                end = list(state)
                end[x] += 1
                additions[x][index][pos[tuple(end)]] = F(1)
            else:
                for y in range(n):
                    if y in full:
                        continue
                    chance = sum(g*K[z][y] for g, z in zip(green[full.index(x)], full))
                    end = list(state)
                    end[y] += 1
                    additions[x][index][pos[tuple(end)]] += chance
    return states, additions


def capacity_audit():
    # Four inputs on the same asymmetric two-site kernel, with a capacity-two
    # loss buffer. Independently integrate its four-dimensional polytope by
    # simplex fibers over the original two-dimensional polygon.
    K = [[F(1, 5), F(2, 5), F(2, 5)],
         [F(1, 2), F(0), F(1, 2)]]
    augmented = [row+[F(0)] for row in K]+[[F(0), F(0), F(0), F(1)]]
    capacities = (1, 1, 2)
    states, U = capacity_additions(augmented, capacities)
    dim = len(states)
    identity = [[F(i == j) for j in range(dim)] for i in range(dim)]
    for x in range(3):
        power = identity
        for _ in range(capacities[x]):
            power = multiply(power, U[x])
        routed = [[sum(augmented[x][y]*U[y][i][j] for y in range(3))
                   for j in range(dim)] for i in range(dim)]
        assert multiply(power, U[x]) == multiply(power, routed)
        for y in range(3):
            assert multiply(U[x], U[y]) == multiply(U[y], U[x])
    b = [F(2, 3), F(1, 3)]
    ell = [[sum(b[x]*U[x][i][j] for x in range(2))
            for j in range(dim)] for i in range(dim)]
    row = [F(i == 0) for i in range(dim)]
    for _ in range(4):
        row = multiply([row], ell)[0]
    capacity_full = row[-1]
    original_sets = supports(2)
    original = idla_addition_matrices(K, original_sets)
    driven = [[sum(b[x]*original[x][i][j] for x in range(2))
               for j in range(4)] for i in range(4)]
    row = [F(i == 0) for i in range(4)]
    for _ in range(4):
        row = multiply([row], driven)[0]
    assert capacity_full == row[-1]

    M = [[F(i == j)-K[i][j] for j in range(2)] for i in range(2)]
    green_T = inverse([[M[j][i] for j in range(2)] for i in range(2)])
    full = [sum(row[j]*b[j] for j in range(2)) for row in green_T]
    polygon = [[F(0), F(0)], [b[0]/M[0][0], F(0)],
               full, [F(0), b[1]/M[1][1]]]
    volume = F(0)
    for a, c in ((polygon[1], polygon[2]), (polygon[2], polygon[3])):
        area = (a[0]*c[1]-a[1]*c[0])/2
        ha, hc = (sum(K[i][2]*v[i] for i in range(2)) for v in (a, c))
        volume += area*(ha*ha+ha*hc+hc*hc)/12
    assert capacity_full == 24*determinant(M)*volume
    probabilities = []
    for at_zero in range(5):
        row = [F(i == 0) for i in range(4)]
        for x in [0]*at_zero+[1]*(4-at_zero):
            row = multiply([row], original[x])[0]
        probabilities.append(row[-1])
    assert all(probabilities[i]**2 >= probabilities[i-1]*probabilities[i+1]
               for i in range(1, 4))
    return {"capacities": capacities, "states": dim,
            "full_filling_probability": str(capacity_full),
            "independent_4d_volume": str(volume),
            "fixed_input_filling_probabilities": list(map(str, probabilities))}


def poisson_cell_audit(current_weight=F(1)):
    """Compare exact exponential expressions from birth ODEs and polytope integrals."""
    K = [[F(1, 5), F(2, 5), F(2, 5)],
         [F(1, 2), F(0), F(1, 2)]]
    b = [F(2, 3), F(1, 3)]
    total = sum(b)
    sets = supports(2)
    pos = {B: i for i, B in enumerate(sets)}
    D = idla_addition_matrices(K, sets)
    diagonal = {B: -total+current_weight*sum(
        b[y]*D[y][pos[B]][pos[B]] for y in range(2)) for B in sets}
    laws = {frozenset(): {-total: F(1)}}
    for C in sorted(sets, key=len):
        if not C:
            continue
        law = defaultdict(F)
        for B in laws:
            rate = sum(b[y]*D[y][pos[B]][pos[C]] for y in range(2))
            if not rate:
                continue
            for exponent, coefficient in laws[B].items():
                assert exponent != diagonal[C]
                factor = rate*coefficient/(exponent-diagonal[C])
                law[exponent] += factor
                law[diagonal[C]] -= factor
        laws[C] = {e: c for e, c in law.items() if c}

    M = [[F(i == j)-K[i][j] for j in range(2)] for i in range(2)]
    integrated = {frozenset(): {-total: F(1)}}
    for x in range(2):
        factor = M[x][x]/(current_weight*K[x][2])
        integrated[frozenset({x})] = {
            -total+current_weight*K[x][2]*b[x]/M[x][x]: factor, -total: -factor}
    green_T = inverse([[M[j][i] for j in range(2)] for i in range(2)])
    full = [sum(row[j]*b[j] for j in range(2)) for row in green_T]
    polygon = [[F(0), F(0)], [b[0]/M[0][0], F(0)],
               full, [F(0), b[1]/M[1][1]]]
    integral = defaultdict(F)
    for a, c in ((polygon[1], polygon[2]), (polygon[2], polygon[3])):
        jacobian = a[0]*c[1]-a[1]*c[0]
        ha, hc = (current_weight*sum(K[i][2]*v[i] for i in range(2)) for v in (a, c))
        assert ha and hc and ha != hc
        # Integral of exp(ha*t+hc*s) over t,s>=0, t+s<=1.
        factor = jacobian*determinant(M)
        integral[-total+ha] += factor/(ha*(ha-hc))
        integral[-total+hc] -= factor/(hc*(ha-hc))
        integral[-total] += factor/(ha*hc)
    integrated[frozenset({0, 1})] = {e: c for e, c in integral.items() if c}
    assert integrated == laws
    return {"exact_support_cells": 4, "exit_weight": str(current_weight),
            "full_law_exponential_terms": {str(e): str(c) for e, c in laws[frozenset({0, 1})].items()},
            "method": "symbolic exponential coefficients; no quadrature or rounding"}


def cycle_all_loss_check(sets, pairing, p):
    """All-loss routing formula on the same single-leak three-cycle."""
    from math import prod

    full = frozenset(range(3))
    a, leak = (F(1), F(1), F(1, 2)), (F(0), F(0), F(1, 2))
    mu = [prod(p[i] if i in B else 1-p[i] for i in range(3)) for B in sets]
    kappa, d = prod(a), 1-prod(a)
    moments = (kappa/d, kappa/d**2, kappa*(1+kappa)/d**3,
               kappa*(1+4*kappa+kappa*kappa)/d**4)
    assert moments == (1, 2, 6, 26)
    R = [[mu[j]*pairing[i][sets.index(full-B)]/d
          for j, B in enumerate(sets)] for i in range(len(sets))]

    def routing(C, B):
        k = len(C)-len(B)
        if k < 0:
            return F(0)
        answer = F(0)
        for e in product(range(k+1), repeat=3):
            if sum(e) != k:
                continue
            factor = prod(leak[i]**e[i] for i in range(3))
            if not factor:
                continue
            h, running = [], 0
            for i in range(3):
                running += int(i in C)-int(i in B)-e[i]
                h.append(running)
            assert running == 0
            low = min(h)
            h = [z-low for z in h]
            assert all(h[i]-h[i-1] == int(i in C)-int(i in B)-e[i]
                       for i in range(3))
            factor *= prod(a[i]**h[i] for i in range(3))
            polynomial = [F(1)]
            for i in range(3):
                for r in range(1, e[i]+1):
                    nxt = [F(0)]*(len(polynomial)+1)
                    for j, coefficient in enumerate(polynomial):
                        nxt[j] += coefficient*F(h[i]+r, r)
                        nxt[j+1] += coefficient*F(1, r)
                    polynomial = nxt
            answer += factor*sum(c*m for c, m in zip(polynomial, moments))
        return answer

    calculated = [[F(i == j)+routing(C, B)*mu[j] for j, B in enumerate(sets)]
                  for i, C in enumerate(sets)]
    assert calculated == R
    assert R[-1][0] == F(3, 10)
    Rinv = inverse(R)
    f = [sum(row) for row in Rinv]
    assert min(f) > 0
    assert all(Rinv[i][j] <= 0 for i in range(8) for j in range(8) if i != j)
    checks = 0
    for v in range(3):
        O = [[F(i == j and v in C) for j in range(8)] for i, C in enumerate(sets)]
        G = multiply(multiply(Rinv, O), R)
        for i, C in enumerate(sets):
            for j in range(8):
                assert (2*int(v in C)-1)*(G[i][j]-O[i][j]) >= 0
                checks += 1
    return {"routing_entries": 64, "target_operator_entries": checks,
            "full_to_empty_remainder": str(R[-1][0])}


def duality_audit():
    from verify_feedback_cone import binary_stabilization
    import verify_response_renewal as rn
    K = [[F(0), F(1), F(0), F(0)],
         [F(0), F(0), F(1), F(0)],
         [F(1, 2), F(0), F(0), F(1, 2)]]
    p = [F(1, 3), F(2, 5), F(3, 5)]
    complementary = [1-x for x in p]
    sets, P, E, D = binary_stabilization(K, p)
    _, P2, E2, _ = binary_stabilization(K, complementary)
    dim = len(sets)
    gram = []
    for C in sets:
        entries = []
        for B in sets:
            row = [F(T == C) for T in sets]
            for i in B:
                row = multiply([row], D[i])[0]
            entries.append(row[-1])
        gram.append(entries)
    transpose = lambda a: list(map(list, zip(*a)))
    H = multiply(multiply(E, gram), transpose(E2))
    cycle_loss_checks = cycle_all_loss_check(sets, H, p)
    assert rank(H) == dim and min(map(min, H)) >= 0
    pi, pi2 = P[-1], P2[-1]
    assert min(pi) > 0 and min(pi2) > 0
    density = {}
    for j, A in enumerate(sets):
        mass = F(1)
        for i in range(3):
            mass *= p[i] if i in A else 1-p[i]
        density[A] = pi[j]/mass
    density_checks = 0
    for A in sets:
        for i in range(3):
            if i not in A:
                assert density[A|{i}] >= density[A]
                density_checks += 1
    assert multiply([pi], H)[0] == [F(1)]*dim
    assert multiply(H, [[x] for x in pi2]) == [[F(1)]]*dim
    link = [[value*pi2[j] for j, value in enumerate(row)] for row in H]
    assert all(sum(row) == 1 for row in link)
    assert link[0] == [F(i == dim-1) for i in range(dim)]
    assert multiply([pi], link)[0] == pi2
    for x in range(3):
        A = multiply(multiply(E, D[x]), P)
        A2 = multiply(multiply(E2, D[x]), P2)
        reversed_A2 = [[A2[j][i]*pi2[j]/pi2[i] for j in range(dim)]
                       for i in range(dim)]
        assert all(sum(row) == 1 for row in reversed_A2)
        assert multiply(A, H) == multiply(H, transpose(A2))
        assert multiply(A, link) == multiply(link, reversed_A2)
    assert all(not H[i][j] for i, C in enumerate(sets) for j, B in enumerate(sets)
               if len(C)+len(B) < 3)

    # Exact target obstruction in the same matrix: after complementing
    # the dual state, undoing the link need not give a positive or
    # inclusion-increasing target payoff.
    full = frozenset(range(3))
    vacancy_link = [[row[sets.index(full-C)] for C in sets] for row in link]
    target = [[F(0 in C)] for C in sets]
    payoff = multiply(inverse(vacancy_link), target)
    assert multiply(vacancy_link, payoff) == target
    assert payoff[sets.index(frozenset({1}))][0] == -F(29, 598)
    assert payoff[sets.index(frozenset({2}))][0] == -F(29, 598)
    assert payoff[sets.index(frozenset({1, 2}))][0] == -F(65493339, 194915225)
    assert payoff[sets.index(frozenset({1, 2}))][0] < payoff[sets.index(frozenset({1}))][0]
    assert all((link[i][j] > 0) == (len(C)+len(B) >= 3)
               for i, C in enumerate(sets) for j, B in enumerate(sets))
    transient = [[F(i == j)-vacancy_link[i+1][j+1]
                  for j in range(dim-1)] for i in range(dim-1)]
    green = inverse(transient)
    assert min(map(min, green)) >= 0
    potential = [[F(0)]] + multiply(green, target[1:])
    evolved_potential = multiply(vacancy_link, potential)
    assert [[a[0]-b[0]] for a, b in zip(potential, evolved_potential)] == target
    singleton = potential[sets.index(frozenset({0}))][0]
    pair = potential[sets.index(frozenset({0, 1}))][0]
    assert singleton == F(189979, 74354)
    assert pair == F(40836185505115661, 21542872302627536) < singleton

    # Independent raw realization of one overlap entry: sleepers {0,2},
    # two forced jumps from overlap site 1, both landing at 2.
    t_sleep = [F(1, 3), F(1, 2), F(2, 5)]
    initial = (-1, 0, 3)
    laws, transient = rn.stopped([initial], lambda state: rn.arw_step(K, t_sleep, state),
                                lambda state: all(x <= 0 for x in state))
    raw = laws[initial].get((-1, -1, -1), F(0))
    _, _, Et, _ = binary_stabilization(K, t_sleep)
    a, b = sets.index(frozenset({0, 1})), sets.index(frozenset({1, 2}))
    assert raw == F(3, 32) and Et[-1][-1] == 8
    assert H[a][b] == raw*Et[-1][-1] == F(3, 4)
    return {"complementary_sleep_vectors": [list(map(str, p)), list(map(str, complementary))],
            "addition_dualities": 3, "stochastic_reversal_links": 3,
            "stationary_density_cover_comparisons": density_checks,
            "cycle_all_loss_routing": cycle_loss_checks,
            "target_vacancy_payoff_at_singleton_1": str(payoff[sets.index(frozenset({1}))][0]),
            "iterated_link_target_potential_inclusion_difference": str(pair-singleton),
            "raw_release_states": transient, "raw_release_probability": str(raw),
            "release_pairing_entry": str(H[a][b])}


def full_minus_one_product_audit():
    """J112: all-source product-average identities, including mandatory portals."""
    from math import prod
    from verify_feedback_cone import binary_stabilization

    def data(K, p):
        n = len(K)
        sets, P, E, D = binary_stabilization(K, p)
        _, _, Eq, _ = binary_stabilization(K, [1-t for t in p])
        full, dim = frozenset(range(n)), len(sets)
        pos = {C: j for j, C in enumerate(sets)}
        gram = []
        for C in sets:
            entries = []
            for B in sets:
                row = [F(T == C) for T in sets]
                for i in B:
                    row = multiply([row], D[i])[0]
                entries.append(row[-1])
            gram.append(entries)
        H = multiply(multiply(E, gram), list(map(list, zip(*Eq))))
        d = determinant([[F(i == j)-K[i][j] for j in range(n)] for i in range(n)])
        mass = [prod(p[i] if i in B else 1-p[i] for i in range(n)) for B in sets]
        R = [[mass[j]*H[a][pos[full-B]]/d for j, B in enumerate(sets)]
             for a in range(dim)]
        Rinv = inverse(R)
        f = [sum(row) for row in Rinv]
        gamma = [[sum(row[j] for j, B in enumerate(sets) if v in B)
                  for v in range(n)] for row in Rinv]
        return sets, pos, P, R, f, gamma, d

    cycle = [[F(0), F(1), F(0), F(0), F(0)],
             [F(0), F(0), F(1), F(0), F(0)],
             [F(0), F(0), F(0), F(1), F(0)],
             [F(1, 2), F(0), F(0), F(0), F(1, 2)]]
    bidirected = [[F(int(j in ((i-1) % 4, (i+1) % 4)), 3) for j in range(4)]
                  + [F(1, 3)] for i in range(4)]
    p = [F(1, 3), F(2, 5), F(3, 7), F(4, 9)]
    coefficient_checks = target_checks = balance_checks = 0
    for K, sleeps in ((cycle, p), (bidirected, p), (bidirected, [F(1, 2)]*4)):
        sets, pos, P, R, f, gamma, d = data(K, sleeps)
        full = frozenset(range(4))
        for x in range(4):
            others = sorted(full-{x})
            h = [F(0)]*4
            green = inverse([[F(i == j)-K[i][j] for j in others] for i in others])
            for i, row in zip(others, green):
                h[i] = sum(g*K[j][4] for g, j in zip(row, others))
            sigma = K[x][4]+sum(K[x][i]*h[i] for i in others)
            Z = prod(sleeps[i]+(1-sleeps[i])*h[i] for i in others)
            theta = {i: sleeps[i]/(sleeps[i]+(1-sleeps[i])*h[i]) for i in others}
            delta = [a-b for a, b in zip(P[-1], P[pos[full-{x}]])]
            formal = multiply([delta], R)[0]
            average = [F(0)]*4
            balance = F(0)
            for B in sets:
                if x in B:
                    continue
                weight = prod(theta[i] if i in B else 1-theta[i] for i in others)
                factor = sleeps[x]*Z*weight/d
                assert formal[pos[B|{x}]] == factor
                assert formal[pos[B]] == -sigma*factor
                coefficient_checks += 2
                balance += weight*(f[pos[B|{x}]]-sigma*f[pos[B]])
                for v in range(4):
                    average[v] += weight*(gamma[pos[B|{x}]][v]-sigma*gamma[pos[B]][v])
            assert balance == 0
            balance_checks += 1
            for v in range(4):
                actual = sum(delta[j] for j, B in enumerate(sets) if v in B)
                assert actual == sleeps[x]*Z*average[v]/d
                target_checks += 1
            if K == bidirected and sleeps == [F(1, 2)]*4 and x == 0:
                assert sigma == F(5, 7) and theta[2] == F(7, 12)
                background = frozenset({1, 3})
                bracket = sum(w*(gamma[pos[background|extra|{x}]][2]
                                 -sigma*gamma[pos[background|extra]][2])
                              for extra, w in ((frozenset(), 1-theta[2]),
                                               (frozenset({2}), theta[2])))
                assert bracket == F(-87744420438699, 6964520196491312)
                gain = sum(delta[j] for j, B in enumerate(sets) if 2 in B)
                assert gain == F(6235498643, 991955589872)
    two = [[F(0), F(1), F(0)], [F(1, 2), F(0), F(1, 2)]]
    sets, pos, _, _, _, gamma, _ = data(two, [F(1, 2)]*2)
    assert gamma[pos[frozenset()]][1] == 0
    assert gamma[pos[frozenset({0})]][1] == F(-8, 49)
    return {"formal_coefficients": coefficient_checks, "target_identities": target_checks,
            "normalization_identities": balance_checks,
            "target_averaged_obstruction": str(bracket), "actual_increment": str(gain)}


def bridge_book_audit():
    """Sixteen decorated books on one existing three-edge current table."""
    from collections import deque
    K = [[F(1, 5), F(2, 5), F(2, 5)],
         [F(1, 2), F(0), F(1, 2)]]
    p = [F(1, 3), F(2, 5)]
    all_sites = frozenset({0, 1})

    def normal_form(book):
        arrivals, emitted, residents = [0, 0], [0, 0], [False, False]
        queue = deque(range(2))
        while queue:
            i = queue.popleft()
            arrivals[i] += 1
            bit = book[i][arrivals[i]-1][0] if arrivals[i] <= len(book[i]) else 1
            assert arrivals[i] <= len(book[i])+1
            outgoing = int(residents[i]) + int(not bit)
            residents[i] = bool(bit)
            for _ in range(outgoing):
                assert emitted[i] < len(book[i])
                queue.append(book[i][emitted[i]][1])
                emitted[i] += 1
        assert all(residents)
        return (tuple(row[:k] for row, k in zip(book, emitted)),
                tuple(row[k:] for row, k in zip(book, emitted)))

    def balanced_cuts(book):
        for cut in product(*(range(len(row)+1) for row in book)):
            incoming = [0, 0]
            for row, k in zip(book, cut):
                for _, j in row[:k]:
                    incoming[j] += 1
            if tuple(incoming) == cut:
                yield (tuple(row[:k] for row, k in zip(book, cut)),
                       tuple(row[k:] for row, k in zip(book, cut)))

    def first_ones(book, sites):
        return all(not book[i] or book[i][0][0] for i in sites)

    def is_bridge(book):
        return normal_form(book)[0] == book

    def greatest_suffix(book, sites):
        candidates = [tail for _, tail in balanced_cuts(book) if first_ones(tail, sites)]
        maximum_lengths = tuple(max(len(tail[i]) for tail in candidates) for i in range(2))
        maxima = [tail for tail in candidates if tuple(map(len, tail)) == maximum_lengths]
        assert len(maxima) == 1
        return maxima[0]

    def graphical_suffix(book, sites):
        positions, letters = [0, 0], []
        while any(positions[i] < len(book[i]) for i in range(2)):
            i = next(i for i in range(2) if positions[i] < len(book[i]))
            route = []
            while i not in route:
                route.append(i)
                assert positions[i] < len(book[i])
                i = book[i][positions[i]][1]
            cycle = route[route.index(i):]
            letters.append({j: book[j][positions[j]] for j in cycle})
            for j in cycle:
                positions[j] += 1
        cemetery = len(letters)
        reverse_edges = [set() for _ in range(cemetery+1)]
        previous = [None, None]
        for current, letter in enumerate(letters):
            for i, (bit, _) in letter.items():
                earlier = previous[i]
                if earlier is not None:
                    reverse_edges[current].add(earlier)  # Earlier -> current.
                if i in sites and not bit:
                    destination = cemetery if earlier is None else earlier
                    reverse_edges[destination].add(current)
                previous[i] = current
        doomed, pending = {cemetery}, [cemetery]
        while pending:
            for earlier in reverse_edges[pending.pop()]:
                if earlier not in doomed:
                    doomed.add(earlier)
                    pending.append(earlier)
        return tuple(tuple(letter[i] for j, letter in enumerate(letters)
                           if j not in doomed and i in letter) for i in range(2))

    def weight(book):
        value = F(1)
        for i, row in enumerate(book):
            for bit, j in row:
                value *= K[i][j] * (p[i] if bit else 1-p[i])
        return value

    cases = 0
    mixed_weights = {order: F(0) for order in ((0, 1), (1, 0))}
    for order in mixed_weights:
        for bits in product((0, 1), repeat=3):
            book = (tuple(zip(bits[:2], order)), ((bits[2], 0),))
            a, b = normal_form(book)
            decompositions = [(head, tail) for head, tail in balanced_cuts(book)
                              if is_bridge(head) and first_ones(tail, all_sites)]
            assert decompositions == [(a, b)]
            assert weight(book) == weight(a)*weight(b)
            assert greatest_suffix(book, all_sites) == b
            projections = {S: greatest_suffix(a, S) for S in supports(2)}
            for S, tail in projections.items():
                assert graphical_suffix(a, S) == tail
                assert graphical_suffix(book, S) == greatest_suffix(book, S)
                assert is_bridge(tail) and first_ones(tail, S)
                for smaller, larger_book in projections.items():
                    if smaller <= S:
                        assert greatest_suffix(larger_book, S) == tail
                        assert all(not row or larger[-len(row):] == row
                                   for row, larger in zip(tail, larger_book))
            if a == book and first_ones(book, {0}):
                mixed_weights[order] += weight(book)
            cases += 1
    total = sum(mixed_weights.values())
    conditional = [mixed_weights[order]/total for order in mixed_weights]
    assert conditional == [F(3, 5), F(2, 5)]
    return {"decorated_balanced_books": cases, "suffix_projections": 4*cases,
            "heap_reachability_checks": 8*cases,
            "mixed_bridge_arrow_order_probabilities": list(map(str, conditional))}


def stable_upper_book_audit():
    """One partial-terminal table and the one-site remainder identity."""
    from collections import deque
    from math import comb
    from verify_feedback_cone import binary_stabilization
    p = [F(1, 3), F(2, 5)]
    K = [[F(0), F(1), F(0)], [F(1, 2), F(0), F(1, 2)]]
    proposed = frozenset({0})

    def weight(book):
        answer = F(1)
        for i, row in enumerate(book):
            for bit, dest in row:
                answer *= K[i][dest] * (p[i] if bit else 1-p[i])
        return answer

    def probability_product(sites):
        answer = F(1)
        for i in sites:
            answer *= p[i]
        return answer

    decompositions, outcomes = set(), defaultdict(int)
    witness = (((0, 1), (0, 1)), ((1, 0), (0, 0)))
    for bits in product((0, 1), repeat=3):
        book = (((bits[0], 1), (bits[1], 1)), ((bits[2], 0), (0, 0)))
        arrivals, emitted, resident = [0, 0], [0, 0], [0, 0]
        queue = deque([0])
        while queue:
            i = queue.popleft()
            arrivals[i] += 1
            assert arrivals[i] <= len(book[i]) + int(i in proposed)
            bit = book[i][arrivals[i]-1][0] if arrivals[i] <= len(book[i]) else 1
            outgoing = resident[i] + 1-bit
            resident[i] = bit
            for _ in range(outgoing):
                assert emitted[i] < len(book[i])
                dest = book[i][emitted[i]][1]
                emitted[i] += 1
                if dest < 2:
                    queue.append(dest)
        actual = frozenset(i for i in range(2) if resident[i])
        head = tuple(row[:k] for row, k in zip(book, emitted))
        tail = tuple(row[k:] for row, k in zip(book, emitted))
        incoming = [0, 0]
        for row in tail:
            for _, dest in row:
                incoming[dest] += 1
        assert all(int(i in actual)+incoming[i]-len(tail[i]) == int(i in proposed)
                   for i in range(2))
        assert all(not tail[i] or tail[i][0][0] == 1 for i in actual)
        assert all(not tail[i] or tail[i][-1][0] == 0 for i in range(2) if i not in proposed)
        lhs = probability_product(proposed)*weight(book)
        rhs = (probability_product(actual)*weight(head)
               *probability_product(proposed)/probability_product(actual)*weight(tail))
        assert lhs == rhs
        decompositions.add((head, actual, tail))
        outcomes[tuple(actual)] += 1
        if book == witness:
            assert actual == frozenset({1}) and emitted == [1, 0]
            assert tail == (((0, 1),), ((1, 0), (0, 0)))
    assert len(decompositions) == 8

    # Direct calculation, including a sink arrow mark and an internal mark.
    one_site_checks = 0
    sleep = F(2, 5)
    awake = 1-sleep
    for internal, exit_weight in ((F(1, 4), F(1, 2)), (F(2, 5), F(1, 3))):
        Q = [[(1-sleep*internal)/(1-internal), F(0)],
             [awake*exit_weight/(1-internal)**2, sleep/(1-internal)]]
        P = [[F(1), F(0)],
             [awake*exit_weight/(1-awake*internal), sleep/(1-awake*internal)]]
        R = [[(1-sleep*internal)/(1-internal), F(0)],
             [awake*exit_weight*internal*(2-internal)/(1-internal)**2,
              (1-awake*internal)/(1-internal)]]
        assert multiply(P, R) == Q
        one_site_checks += 1

    # Sum every remainder book on the same two-site single-leak cycle.
    # At site 1 let m be the emitted count and k the forced sink count;
    # its arrow words have binom(m,k) orders. For m>=3 the bit factor
    # is constant, so the infinite tail is an exact geometric derivative.
    sets, P, E, D = binary_stabilization(K, p)
    q = [1-x for x in p]
    _, Pq, Eq, _ = binary_stabilization(K, q)

    def remainders(sleeps):
        result = [[F(0) for _ in sets] for _ in sets]
        a, b, leak = K[0][1], K[1][0], K[1][2]
        t = a*b
        for ci, C in enumerate(sets):
            for bi, B in enumerate(sets):
                k = len(C)-len(B)
                if k < 0:
                    continue
                shift = int(1 in B)-int(1 in C)
                mu_B, boundary = F(1), F(1)
                for i in range(2):
                    mu_B *= sleeps[i] if i in B else 1-sleeps[i]
                    boundary *= (sleeps[i] if i in B else F(1))/(sleeps[i] if i in C else F(1))
                tail = t**k/(1-t)**(k+1)-sum(F(comb(m, k))*t**m for m in range(k, 3))
                value = mu_B*a**shift*b**(-k)*leak**k*tail
                for m in range(3):
                    m0 = m+shift
                    if m < k or m0 < 0:
                        continue
                    bit_weight = F(1)
                    for i, length in enumerate((m0, m)):
                        if not length:
                            continue
                        if length == 1 and i in C and i not in B:
                            bit_weight = F(0)
                            break
                        if i in C:
                            bit_weight *= sleeps[i]
                        if i not in B:
                            bit_weight *= 1-sleeps[i]
                    value += boundary*bit_weight*comb(m, k)*a**m0*b**(m-k)*leak**k
                result[ci][bi] = value
        return result

    R, Rq = remainders(p), remainders(q)
    gram = []
    for C in sets:
        entries = []
        for B in sets:
            row = [F(T == C) for T in sets]
            for i in B:
                row = multiply([row], D[i])[0]
            entries.append(row[-1])
        gram.append(entries)
    H = multiply(multiply(E, gram), list(map(list, zip(*Eq))))
    det = 1-K[0][1]*K[1][0]
    all_sites = frozenset({0, 1})
    mu = []
    factors = []
    for B in sets:
        mass = F(1)
        for i in range(2):
            mass *= p[i] if i in B else q[i]
        mu.append(mass)
        factors.append(det*Pq[-1][sets.index(all_sites-B)]/mass)
    for ci, C in enumerate(sets):
        for bi, B in enumerate(sets):
            assert R[ci][bi] == mu[bi]*H[ci][sets.index(all_sites-B)]/det, (C, B, R[ci][bi], mu[bi]*H[ci][sets.index(all_sites-B)]/det)
            assert mu[ci]*R[ci][bi] == mu[bi]*Rq[sets.index(all_sites-B)][sets.index(all_sites-C)]
        assert sum(R[ci][bi]*factors[bi] for bi in range(len(sets))) == 1
    schur_checks = 0
    for ci, C in enumerate(sets):
        for di, Dcut in enumerate(sets):
            first_cut_mass = R[ci][di]/R[di][di]
            assert first_cut_mass == 1 if ci == di else first_cut_mass <= 1-factors[ci]
            for bi, B in enumerate(sets):
                avoidance = R[ci][bi]-first_cut_mass*R[di][bi]
                assert avoidance >= 0
                if ci == bi and ci != di:
                    assert avoidance >= 1
                schur_checks += 1
    return {"partial_terminal_books": 8,
            "actual_terminal_counts": {str(k): v for k, v in outcomes.items()},
            "witness_proposed_terminal": [0], "witness_actual_terminal": [1],
            "one_site_positive_factorizations": one_site_checks,
            "two_site_remainder_entries": 16,
            "two_site_reciprocity_entries": 16,
            "single_endpoint_Schur_checks": schur_checks,
            "stochastic_column_normalization": list(map(str, factors))}


def _endpoint_fiber_cuts(book, C, B, D):
    """Finite IDLA prefix, maximal cycle heap, then constrained ideals."""
    from collections import deque
    n = len(book)
    q = [0]*n
    occupied = set()
    queue = deque(i for i in range(n) for _ in range(int(i in C)+int(i not in D)))
    while queue:
        i = queue.popleft()
        if i not in occupied:
            occupied.add(i)
            continue
        if q[i] == len(book[i]):
            return set()
        dest = book[i][q[i]][1]
        q[i] += 1
        if dest < n:
            queue.append(dest)
    if len(occupied) != n:
        return set()

    tail = [row[k:] for row, k in zip(book, q)]
    pos = [0]*n
    letters = []
    while True:
        cycle = None
        for start in range(n):
            route, seen = [], {}
            i = start
            while i < n and pos[i] < len(tail[i]) and i not in seen:
                seen[i] = len(route)
                route.append(i)
                i = tail[i][pos[i]][1]
            if i in seen:
                cycle = route[seen[i]:]
                break
        if cycle is None:
            break
        letters.append({i: tail[i][pos[i]] for i in cycle})
        for i in cycle:
            pos[i] += 1

    source, dead = len(letters), len(letters)+1
    edges = [set() for _ in range(len(letters)+2)]
    for i in range(n):
        chain = [j for j, letter in enumerate(letters) if i in letter]
        for earlier, later in zip(chain, chain[1:]):
            edges[later].add(earlier)
        for k in range(len(chain)+1):
            arrival = q[i]+k+int(i in D)
            if arrival == 0:
                bit = 0
            elif arrival <= len(book[i]):
                bit = book[i][arrival-1][0]
            else:
                assert i in B and arrival == len(book[i])+1
                bit = 1
            if bit == int(i in D):
                continue
            if k == 0:
                edges[source].add(chain[0] if chain else dead)
            elif k == len(chain):
                edges[chain[-1]].add(dead)
            else:
                edges[chain[k-1]].add(chain[k])

    # Exhaust only the tiny extracted heaps in these two certificates.
    cuts = set()
    for mask in product((False, True), repeat=len(letters)):
        chosen = {source} | {j for j, yes in enumerate(mask) if yes}
        if any(not edges[j] <= chosen for j in chosen):
            continue
        cuts.add(tuple(q[i]+sum(i in letters[j] for j in range(len(letters)) if mask[j])
                       for i in range(n)))
    return cuts



def stable_cut_lattice_audit():
    """Exhaust cuts of the twelve-record M3 and eight-record N5 examples."""
    arrows = ((3, 5), (4, 3), (5, 4), (0, 1), (1, 2), (2, 0))
    initial = (1, 1, 1, 0, 0, 0)
    zero, top = (0,)*6, (2,)*6

    def cut_map(m, words=arrows, occupied=initial):
        incoming = [0]*len(words)
        for row, length in zip(words, m):
            for j in row[:length]:
                incoming[j] += 1
        arrivals = [occupied[i]+incoming[i] for i in range(len(words))]
        return tuple(n-(0 if n == 0 else (1, 0, 1)[n-1]) for n in arrivals)

    fixed = {m for m in product(range(3), repeat=6) if cut_map(m) == m}
    atoms = ((2, 2, 0, 2, 0, 0), (0, 2, 2, 0, 2, 0), (2, 0, 2, 0, 0, 2))
    assert fixed == {zero, top, *atoms}
    for i, a in enumerate(atoms):
        for b in atoms[i+1:]:
            assert cut_map(tuple(min(x, y) for x, y in zip(a, b))) == zero
            assert cut_map(tuple(max(x, y) for x, y in zip(a, b))) == top
    n5_arrows = ((0, 3), (1, 0), (2, 3), (1, 2))
    n5_initial = (1, 1, 1, 0)
    n5_map = lambda m: cut_map(m, n5_arrows, n5_initial)
    n5_fixed = {m for m in product(range(3), repeat=4) if n5_map(m) == m}
    u, v, w = (2, 0, 0, 0), (2, 2, 0, 0), (0, 0, 2, 0)
    assert n5_fixed == {(0,)*4, u, v, w, (2,)*4}
    assert n5_map(tuple(min(x, y) for x, y in zip(v, w))) == (0,)*4
    proposed_join = tuple(max(x, y) for x, y in zip(u, w))
    assert n5_map(n5_map(proposed_join)) == (2,)*4

    def mobius_and_fibers(cuts, words, occupied):
        mu = {}
        endpoints = {}
        for m in sorted(cuts, key=sum):
            mu[m] = (1 if not any(m) else
                     -sum(value for x, value in mu.items()
                          if all(a <= b for a, b in zip(x, m))))
            incoming = [0]*len(words)
            for row, length in zip(words, m):
                for j in row[:length]:
                    incoming[j] += 1
            endpoints[m] = tuple(occupied[i]+incoming[i]-m[i] for i in range(len(words)))
        for a in cuts:
            for b in cuts:
                if endpoints[a] == endpoints[b]:
                    low = tuple(min(x, y) for x, y in zip(a, b))
                    high = tuple(max(x, y) for x, y in zip(a, b))
                    assert low in cuts and high in cuts
                    assert endpoints[low] == endpoints[high] == endpoints[a]
        return mu[(2,)*len(words)], endpoints

    m3_mu, m3_endpoints = mobius_and_fibers(fixed, arrows, initial)
    n5_mu, n5_endpoints = mobius_and_fibers(n5_fixed, n5_arrows, n5_initial)
    assert (m3_mu, n5_mu) == (2, 1)
    endpoint_queries = 0
    for words, occupied, endpoints in ((arrows, initial, m3_endpoints),
                                        (n5_arrows, n5_initial, n5_endpoints)):
        C = frozenset(i for i, yes in enumerate(occupied) if yes)
        book = tuple(tuple(zip((1, 0), row)) for row in words)
        for marks in product((False, True), repeat=len(words)):
            D = frozenset(i for i, yes in enumerate(marks) if yes)
            expected = {m for m, endpoint in endpoints.items()
                        if frozenset(i for i, yes in enumerate(endpoint) if yes) == D}
            assert _endpoint_fiber_cuts(book, C, C, D) == expected
            endpoint_queries += 1
    return {"endpoint_fiber_queries": endpoint_queries,
            "M3": {"cuts_checked": 3**6, "fixed_cuts": len(fixed),
                   "middle_cuts": [list(a) for a in atoms], "mobius_top": m3_mu},
            "N5": {"cuts_checked": 3**4, "fixed_cuts": len(n5_fixed),
                   "middle_cuts": [list(a) for a in (u, v, w)], "mobius_top": n5_mu}}


def target_cut_cancellation_audit():
    """A four-site, one-particle obstruction to a fixed-book target sign."""
    words = ((1, 2, 3), (0, 1), (0, 2), ())
    bits = ((1, 0, 0), (1, 0), (1, 0), ())
    C, B = (1, 0, 0, 0), (0, 0, 0, 1)
    def cut_data(arrow_words):
        endpoints = {}
        candidates = 0
        for m in product(*(range(len(row)+1) for row in arrow_words)):
            candidates += 1
            arrivals = list(C)
            for row, length in zip(arrow_words, m):
                for j in row[:length]:
                    arrivals[j] += 1
            assert all(arrivals[i] <= len(arrow_words[i])+B[i] for i in range(4))
            resident = tuple(0 if a == 0 else bits[i][a-1]
                             if a <= len(bits[i]) else B[i]
                             for i, a in enumerate(arrivals))
            if tuple(a-b for a, b in zip(arrivals, resident)) == m:
                endpoints[m] = resident
        mu = {}
        for m in sorted(endpoints, key=sum):
            mu[m] = 1 if not any(m) else -sum(
                value for u, value in mu.items()
                if all(x <= y for x, y in zip(u, m)))
        coefficients = tuple(sum(mu[m]*D[v] for m, D in endpoints.items())
                             for v in range(4))
        return endpoints, mu, coefficients, candidates

    endpoints, mu, target_coefficients, candidates = cut_data(words)
    assert endpoints == {(0, 0, 0, 0): C,
                         (2, 2, 0, 0): (0, 0, 1, 0),
                         (2, 0, 2, 0): (0, 1, 0, 0),
                         (3, 2, 2, 0): B}
    assert target_coefficients == (1, -1, -1, 1)
    assert C[3] == 0 and target_coefficients[3] > 0
    # The only other admissible ordering of these paired-record multisets
    # exchanges the two root zero records and cancels the bad target term.
    swapped = ((1, 3, 2), *words[1:])
    partner, partner_mu, partner_coefficients, partner_count = cut_data(swapped)
    assert partner == {(0, 0, 0, 0): C, (2, 2, 0, 0): B, (3, 2, 2, 0): B}
    assert partner_mu[(3, 2, 2, 0)] == 0
    assert partner_coefficients == (1, 0, 0, -1)
    assert target_coefficients[3]+partner_coefficients[3] == 0

    # Summing all books restores the favorable sign in this one example.
    # Use the independently proved complementary-pairing formula for R.
    from verify_feedback_cone import binary_stabilization
    K = [[F(0), F(1, 4), F(1, 4), F(1, 4), F(1, 4)],
         [F(1, 3), F(1, 3), F(0), F(0), F(1, 3)],
         [F(1, 3), F(0), F(1, 3), F(0), F(1, 3)],
         [F(0), F(0), F(0), F(0), F(1)]]
    sets, _, E, D = binary_stabilization(K, [F(1, 2)]*4)
    gram = []
    for S in sets:
        entries = []
        for T in sets:
            row = [F(A == S) for A in sets]
            for i in T:
                row = multiply([row], D[i])[0]
            entries.append(row[-1])
        gram.append(entries)
    H = multiply(multiply(E, gram), [list(row) for row in zip(*E)])
    d = determinant([[F(i == j)-K[i][j] for j in range(4)] for i in range(4)])
    singletons = [sets.index(frozenset({i})) for i in range(4)]
    R1 = [[H[i][sets.index(frozenset(range(4))-sets[j])]/(16*d)
           for j in singletons] for i in singletons]
    assert R1[3] == [F(0), F(0), F(0), F(383, 224)]
    averaged = inverse(R1)[0][3]*R1[3][3]
    assert averaged == -F(5625, 238336)
    return {"cut_candidates": candidates, "fixed_cuts": len(endpoints),
            "mobius_top": mu[(3, 2, 2, 0)],
            "target_coefficients": target_coefficients,
            "partner_cut_candidates": partner_count,
            "partner_target_coefficients": partner_coefficients,
            "averaged_target_entry": str(averaged)}


def target_inventory_cancellation_audit():
    """One fixed three-particle inventory, including all 192 local orders."""
    from itertools import permutations
    base = ((3, 5, 6), (4, 3), (5, 4), (0, 1), (1, 2), (2, 0), ())
    bits = ((1, 0, 0), (1, 0), (1, 0), (1, 0), (1, 0), (1, 0), ())
    C, B = (1, 1, 1, 0, 0, 0, 0), (0, 1, 1, 0, 0, 0, 1)
    # The displayed bit words force these emission-coordinate possibilities.
    candidates = tuple(product((0, 2, 3), *((0, 2),)*5, (0,)))
    histogram = defaultdict(int)
    for words in product(*(tuple(permutations(row)) for row in base)):
        endpoints = {}
        for m in candidates:
            arrivals = list(C)
            for row, length in zip(words, m):
                for j in row[:length]:
                    arrivals[j] += 1
            assert all(a <= len(words[i])+B[i] for i, a in enumerate(arrivals))
            resident = tuple(0 if a == 0 else bits[i][a-1]
                             if a <= len(bits[i]) else B[i]
                             for i, a in enumerate(arrivals))
            if tuple(a-b for a, b in zip(arrivals, resident)) == m:
                endpoints[m] = resident
        mu = {}
        for m in sorted(endpoints, key=sum):
            mu[m] = 1 if not any(m) else -sum(
                value for u, value in mu.items()
                if all(x <= y for x, y in zip(u, m)))
        coefficient = sum(mu[m]*D[6] for m, D in endpoints.items())
        assert coefficient == (2 if words[0][-1] == 6 else 0)
        histogram[coefficient] += 1
    assert dict(histogram) == {2: 64, 0: 128}
    return {"local_orders": sum(histogram.values()),
            "target_coefficient_histogram": dict(histogram),
            "grouped_target_coefficient": sum(k*v for k, v in histogram.items())}


def one_particle_mobius_audit():
    """A five-site one-particle M3 with a two-record root."""
    words = ((1, 2), (3, 4, 0, 0), (1, 2), (1, 3), (1, 4))
    bits = ((1, 0), (1, 1, 0, 0), (1, 0), (1, 0), (1, 0))
    C = (1, 0, 0, 0, 0)
    endpoints = {}
    candidates = 0
    for m in product(*(range(len(row)+1) for row in words)):
        candidates += 1
        arrivals = list(C)
        for row, length in zip(words, m):
            for j in row[:length]:
                arrivals[j] += 1
        resident = tuple(0 if a == 0 else bits[i][a-1]
                         if a <= len(bits[i]) else C[i]
                         for i, a in enumerate(arrivals))
        if tuple(a-b for a, b in zip(arrivals, resident)) == m:
            endpoints[m] = resident
    top = (2, 4, 2, 2, 2)
    atoms = ((2, 3, 0, 2, 2), (2, 3, 2, 0, 2), (2, 3, 2, 2, 0))
    assert set(endpoints) == {(0,)*5, top, *atoms}
    assert [endpoints[a].index(1) for a in atoms] == [2, 3, 4]
    mu = {}
    for m in sorted(endpoints, key=sum):
        mu[m] = 1 if not any(m) else -sum(
            value for u, value in mu.items()
            if all(x <= y for x, y in zip(u, m)))
    assert mu[top] == 2
    return {"cut_candidates": candidates, "stable_cuts": len(endpoints),
            "transporting_atoms": len(atoms), "mobius_top": mu[top]}


def greatest_cut_field_audit():
    """Two finite arrow fields, all bits and binary inputs; no truncation."""
    # The next arrow at both sites is a sink arrow, so complete cycle
    # popping after draining gives the exact universal bound (2,2).
    # Reversal changes the drain/cycle split but not the greatest drain.
    arrow_fields = (((1, 2), (0, 2)), ((2, 1), (2, 0)))
    inputs = tuple(product((0, 1), repeat=2))
    bit_fields = tuple(product((0, 1), repeat=4))
    candidate_cuts = tuple(product(range(3), repeat=2))
    le = lambda a, b: all(x <= y for x, y in zip(a, b))
    field_checks = adjunction_checks = sleep_checks = 0
    drain_splits = []
    for arrows in arrow_fields:
        def next_arrow(i, m):
            return arrows[i][m] if m < 2 else 2

        mass, drain = [1, 1], [0, 0]
        while any(mass):
            i = next(i for i in range(2) if mass[i])
            j = next_arrow(i, drain[i])
            mass[i] -= 1
            drain[i] += 1
            if j < 2:
                mass[j] += 1
        bound = drain[:]
        while next_arrow(0, bound[0]) == 1 and next_arrow(1, bound[1]) == 0:
            bound = [x+1 for x in bound]
        assert tuple(bound) == (2, 2)
        assert all(next_arrow(i, bound[i]) == 2 for i in range(2))
        drain_splits.append((tuple(drain), tuple(bound[i]-drain[i] for i in range(2))))

        def arrivals(eta, m):
            result = list(eta)
            for i in range(2):
                for j in arrows[i][:m[i]]:
                    if j < 2:
                        result[j] += 1
            return result

        data = {}
        for bits in bit_fields:
            def buffer_map(eta, m):
                a = arrivals(eta, m)
                assert all(x <= 2 for x in a)
                return tuple(x-(bits[2*i+x-1] if x else 0)
                             for i, x in enumerate(a))

            for eta in inputs:
                cuts = tuple(m for m in candidate_cuts if buffer_map(eta, m) == m)
                lower = next(m for m in cuts if all(le(m, n) for n in cuts))
                upper = next(m for m in cuts if all(le(n, m) for n in cuts))
                for start, expected in (((0, 0), lower), ((2, 2), upper)):
                    m = start
                    for _ in range(5):
                        n = buffer_map(eta, m)
                        if n == m:
                            break
                        m = n
                    assert m == expected
                # Independent arrival-buffer execution checks the lower cut.
                queue = [i for i, amount in enumerate(eta) for _ in range(amount)]
                received, emitted, resident = [0, 0], [0, 0], [0, 0]
                while queue:
                    i = queue.pop(0)
                    received[i] += 1
                    bit = bits[2*i+received[i]-1]
                    count = resident[i]+1-bit
                    resident[i] = bit
                    for _ in range(count):
                        j = arrows[i][emitted[i]]
                        emitted[i] += 1
                        if j < 2:
                            queue.append(j)
                assert tuple(emitted) == lower
                endpoints = {m: tuple(a-x for a, x in zip(arrivals(eta, m), m))
                             for m in cuts}
                data[bits, eta] = cuts, lower, upper, endpoints
                field_checks += 1

        for bits in bit_fields:
            for eta in inputs:
                cuts, lower, upper, endpoints = data[bits, eta]
                for larger in inputs:
                    if not le(eta, larger):
                        continue
                    new_cuts = data[bits, larger][0]
                    alpha = {m: next(n for n in new_cuts if le(m, n)
                                      and all(le(n, z) for z in new_cuts if le(m, z)))
                             for m in cuts}
                    beta = {n: next(m for m in cuts if le(m, n)
                                     and all(le(z, m) for z in cuts if le(z, n)))
                            for n in new_cuts}
                    for m in cuts:
                        for n in new_cuts:
                            assert le(alpha[m], n) == le(m, beta[n]) == le(m, n)
                            adjunction_checks += 1
                for record, bit in enumerate(bits):
                    if bit:
                        continue
                    changed = bits[:record]+(1,)+bits[record+1:]
                    _, new_lower, new_upper, new_endpoints = data[changed, eta]
                    assert le(new_lower, lower) and le(new_upper, upper)
                    site = record//2
                    if endpoints[upper][site]:
                        assert new_upper == upper and new_endpoints[upper] == endpoints[upper]
                    if not new_endpoints[new_lower][site]:
                        assert new_lower == lower and new_endpoints[lower] == endpoints[lower]
                    sleep_checks += 1
    assert drain_splits == [((2, 2), (0, 0)), ((1, 1), (1, 1))]
    return {"complete_bit_input_fields": field_checks,
            "adjunction_cut_pairs": adjunction_checks,
            "single_bit_sleep_comparisons": sleep_checks,
            "drain_cycle_splits": drain_splits}


def envelope_inventory_audit():
    """One four-order envelope inventory and its coin-preserving involution."""
    delta = (4, 1, 2)
    data, counts, sectors = {}, [], defaultdict(lambda: [0, 0, 0])
    for t in range(4):
        root = tuple(2 if j == t else 3 for j in range(4))
        arrows = (root, (0,), (0, 0))
        occupied = [[0, 0, 0], [0, 0, 0]]
        loss = gain = 0
        for field in product((0, 1), repeat=7):
            bits = (field[:4], field[4:5], field[5:])
            ends = []
            for ci, eta in enumerate(((1, 0, 0), (1, 1, 0))):
                m = (0, 0, 0)
                for _ in range(8):
                    arrivals = list(eta)
                    for row, length in zip(arrows, m):
                        for dest in row[:length]:
                            if dest < 3:
                                arrivals[dest] += 1
                    assert all(a <= d for a, d in zip(arrivals, delta))
                    endpoint = tuple(bits[i][a-1] if a else 0
                                     for i, a in enumerate(arrivals))
                    new = tuple(a-b for a, b in zip(arrivals, endpoint))
                    assert all(x <= y <= z for x, y, z in zip(m, new, delta))
                    if new == m:
                        break
                    m = new
                else:
                    raise AssertionError("finite envelope failed to stabilize")
                ends.append(endpoint)
                for v in range(3):
                    occupied[ci][v] += endpoint[v]
            difference = tuple(b-a for a, b in zip(*ends))
            assert difference[1] >= 0 and difference[2] >= 0
            data[t, field] = difference
            loss += difference[0] < 0
            gain += difference[0] > 0
            key = (sum(field[:4]), field[4], sum(field[5:]))
            for v in range(3):
                sectors[key][v] += difference[v]
        counts.append((occupied[0][0], occupied[1][0], loss, gain))

    def partner(t, field):
        x1, x2, x3, x4, y, z, z2 = field
        if y:
            return t, field
        if not z and t == 0:
            return 1, (x3, x1, x2, x4, y, z, z2)
        if not z and t == 1:
            return 0, (x2, x3, x1, x4, y, z, z2)
        return t, (x2, x1, x3, x4, y, z, z2)

    for (t, field), difference in data.items():
        other = partner(t, field)
        assert partner(*other) == (t, field)
        assert data[other][0] == -difference[0]
        assert (sum(field[:4]), field[4], sum(field[5:])) == (
            sum(other[1][:4]), other[1][4], sum(other[1][5:]))
    assert counts == [(80, 72, 20, 12), (64, 72, 12, 20),
                      (64, 64, 16, 16), (64, 64, 16, 16)]
    assert len(sectors) == 30
    assert all(s[0] == 0 and min(s) >= 0 for s in sectors.values())
    return {"decorated_envelopes": len(data), "coin_count_sectors": len(sectors),
            "root_occupation_loss_gain_counts": counts}


def stationary_lattice_audit():
    """Two fixed stationary laws; evidence only for log-submodularity."""
    from itertools import combinations
    from verify_feedback_cone import binary_stabilization

    cases = [
        ([[0, 1, 0, 0], [0, 0, 1, 0], [F(1, 2), 0, 0, F(1, 2)]],
         [F(1, 3), F(2, 5), F(3, 5)], (0, 6, 0)),
        ([[0, F(1, 4), F(1, 4), F(1, 4), F(1, 4)],
          [F(1, 3), F(1, 3), 0, 0, F(1, 3)],
          [F(1, 3), 0, F(1, 3), 0, F(1, 3)], [0, 0, 0, 0, 1]],
         [F(1, 2)]*4, (0, 12, 12))]
    results = []
    for raw_K, p, expected in cases:
        K = [[F(x) for x in row] for row in raw_K]
        n, V = len(p), frozenset(range(len(p)))
        sets, P, E, _ = binary_stabilization(K, p)
        pi = dict(zip(sets, P[sets.index(V)]))
        assert sum(pi.values()) == 1 and min(pi.values()) > 0
        assert all(sum(pi[B]*E[b][a] for b, B in enumerate(sets)) == F(A == V)
                   for a, A in enumerate(sets))
        minors = []
        for i, j in combinations(range(n), 2):
            for A in sets:
                if i not in A and j not in A:
                    minors.append(pi[A]*pi[A | {i, j}]-pi[A | {i}]*pi[A | {j}])
        signs = (sum(d > 0 for d in minors), sum(d < 0 for d in minors),
                 sum(d == 0 for d in minors))
        assert signs == expected
        if n == 3:
            mass = sum(pi[A] for A in sets if 0 not in A)
            assert [pi[A]/mass for A in (frozenset(), frozenset({1}),
                                         frozenset({2}), frozenset({1, 2}))] == [
                                             F(w, 688) for w in (75, 112, 225, 276)]
        else:
            assert all(pi[A] == pi[A | {3}] for A in sets if 3 not in A)
        results.append({"vertices": n, "positive_negative_zero_minors": signs})
    return results


def cyclic_envelope_audit():
    """One 768-word inventory: physical signs, burning, and phase counts."""
    from math import comb

    lengths = (3, 2, 2)
    orders = [(tuple(1 if j == pos else 3 for j in range(3)), row1, (0, 0))
              for pos in range(3) for row1 in ((3, 2), (2, 3))]
    population, counts = defaultdict(int), [defaultdict(int) for _ in range(3)]
    pairs, recurrent, histograms = {}, {}, [defaultdict(int) for _ in range(3)]
    rotation_groups, drain_images = defaultdict(dict), {}
    phase_checks = rotor_checks = 0

    def arrivals(arrows, emitted, eta):
        result = list(eta)
        for row, used in zip(arrows, emitted):
            for dest in row[:used]:
                if dest < 3:
                    result[dest] += 1
        return tuple(result)

    def cut_map(arrows, bits, eta, emitted):
        n = arrivals(arrows, emitted, eta)
        assert all(0 <= x <= length for x, length in zip(n, lengths))
        return tuple(x - (bits[i][x-1] if x else 0) for i, x in enumerate(n))

    def fixed_cut(arrows, bits, eta, greatest=False):
        emitted = lengths if greatest else (0, 0, 0)
        for _ in range(sum(lengths) + 1):
            nxt = cut_map(arrows, bits, eta, emitted)
            assert all((b <= a if greatest else a <= b)
                       for a, b in zip(emitted, nxt))
            if nxt == emitted:
                n = arrivals(arrows, emitted, eta)
                return emitted, tuple(a-b for a, b in zip(n, emitted))
            emitted = nxt
        raise AssertionError("envelope iteration did not terminate")

    def buffer_add(arrows, bits, initial, letters, reset=False):
        phase, queue = list(initial), list(letters)
        resident = [0 if reset else bits[i][(a-1) % lengths[i]]
                    for i, a in enumerate(initial)]
        processed, emitted = [0, 0, 0], [0, 0, 0]
        for _ in range(sum(lengths) + 1):
            if not any(queue):
                return tuple(phase)
            i = next(i for i in range(3) if queue[i])
            queue[i] -= 1
            a, length = phase[i], lengths[i]
            old, new = resident[i], bits[i][a]
            head = (a-old) % length
            out = 1 + old - new
            for offset in range(out):
                dest = arrows[i][(head+offset) % length]
                if dest < 3:
                    queue[dest] += 1
            phase[i] = (a+1) % length
            resident[i] = new
            processed[i] += 1
            emitted[i] += out
            assert processed[i] <= length and emitted[i] <= length
        raise AssertionError("burn exceeded one-period upper cut")

    def rotor_drain(arrows, initial_heads, particles):
        heads, queue = list(initial_heads), list(particles)
        for _ in range(sum(lengths)+1):
            if not any(queue):
                return tuple(heads)
            i = next(i for i in range(3) if queue[i])
            queue[i] -= 1
            dest = arrows[i][heads[i]]
            heads[i] = (heads[i]+1) % lengths[i]
            if dest < 3:
                queue[dest] += 1
        raise AssertionError("rotor drain exceeded its full-period upper bound")

    def drain_residents(arrows, bits, phase):
        residents = tuple(bits[i][(a-1) % lengths[i]] for i, a in enumerate(phase))
        heads = tuple((a-c) % m for a, c, m in zip(phase, residents, lengths))
        return rotor_drain(arrows, heads, residents)

    def is_last_forest(arrows, heads):
        last = tuple(arrows[i][(h-1) % lengths[i]] for i, h in enumerate(heads))
        for start in range(3):
            seen, current = set(), start
            while current < 3:
                if current in seen:
                    return False
                seen.add(current)
                current = last[current]
        return True

    for arrows in orders:
        for field in product((0, 1), repeat=7):
            bits = (field[:3], field[3:5], field[5:])
            k = tuple(map(sum, bits))
            population[k] += 1
            _, low = fixed_cut(arrows, bits, (1, 0, 0))
            _, high = fixed_cut(arrows, bits, (1, 1, 0))
            paired_word = tuple(tuple(zip(a, b)) for a, b in zip(arrows, bits))
            canonical = tuple(min(row[e:]+row[:e] for e in range(len(row)))
                              for row in paired_word)
            rotation_groups[canonical][paired_word] = tuple(h-l for l, h in zip(low, high))
            for v in range(3):
                delta = high[v] - low[v]
                counts[v][k] += delta
                histograms[v][delta] += 1
            rec_count, rec_images = 0, set()
            for phase in product(*(range(m) for m in lengths)):
                is_rec = buffer_add(arrows, bits, phase, (1, 1, 1)) == phase
                rec_count += is_rec
                phase_checks += 1
                rotor = drain_residents(arrows, bits, phase)
                if is_rec:
                    assert is_last_forest(arrows, rotor)
                    rec_images.add(rotor)
                for x in range(3):
                    letter = tuple(int(i == x) for i in range(3))
                    changed = buffer_add(arrows, bits, phase, letter)
                    assert drain_residents(arrows, bits, changed) == rotor_drain(
                        arrows, rotor, letter)
                    rotor_checks += 1
                C = tuple(bits[i][(a-1) % lengths[i]] for i, a in enumerate(phase))
                shifts = tuple((a-c) % m for a, c, m in zip(phase, C, lengths))
                rotated_arrows = tuple(row[e:]+row[:e] for row, e in zip(arrows, shifts))
                rotated_bits = tuple(row[e:]+row[:e] for row, e in zip(bits, shifts))
                assert all(row[0] == 1 if c else row[-1] == 0
                           for row, c in zip(rotated_bits, C))
                key = (rotated_arrows, rotated_bits, C)
                pairs[key] = pairs.get(key, 0) + 1
                if key in recurrent:
                    assert recurrent[key] == is_rec
                recurrent[key] = is_rec
            assert rec_count == 10  # det(diag(3,2,2)-J^T)
            assert len(rec_images) == 10
            assert rec_images == {h for h in product(*(range(m) for m in lengths))
                                  if is_last_forest(arrows, h)}
            for h in rec_images:
                preimage = buffer_add(arrows, bits, h, (1, 1, 1), reset=True)
                assert buffer_add(arrows, bits, preimage, (1, 1, 1)) == preimage
                assert drain_residents(arrows, bits, preimage) == h
            if is_last_forest(arrows, (0, 0, 0)):
                used, C = fixed_cut(arrows, bits, (1, 1, 1))
                shifted_arrows = tuple(row[e:]+row[:e] for row, raw in zip(arrows, used)
                                       for e in [raw % len(row)])
                shifted_bits = tuple(row[e:]+row[:e] for row, raw in zip(bits, used)
                                     for e in [raw % len(row)])
                key = (shifted_arrows, shifted_bits, C)
                assert key not in drain_images
                drain_images[key] = (arrows, bits)

    def add(*polys):
        result = defaultdict(int)
        for poly in polys:
            for exponent, coefficient in poly.items():
                result[exponent] += coefficient
        return dict(result)

    def mul(*polys):
        result = {(0, 0, 0): 1}
        for poly in polys:
            nxt = defaultdict(int)
            for a, ca in result.items():
                for b, cb in poly.items():
                    nxt[tuple(x+y for x, y in zip(a, b))] += ca*cb
            result = dict(nxt)
        return result

    one, two, three = ({(0, 0, 0): n} for n in (1, 2, 3))
    x, y, z = ({(1, 0, 0): 1}, {(0, 1, 0): 1}, {(0, 0, 1): 1})
    ox, oy, oz = add(one, x), add(one, y), add(one, z)
    expected = [mul(x, ox, y, oy, oz),
                mul(ox, ox, y, oz, add(ox, mul(two, add(two, mul(three, x)), oy, oz))),
                mul(ox, ox, oy, z, oz, add(three, mul(three, x), y))]
    assert len(population) == 36 and sum(population.values()) == 768
    for k, number in population.items():
        assert number == 6 * comb(3, k[0]) * comb(2, k[1]) * comb(2, k[2])
    for actual, predicted in zip(counts, expected):
        assert {k: c for k, c in actual.items() if c} == predicted
        assert min(actual.values()) >= 0
    assert [sum(c.values()) for c in counts] == [8, 336, 112]
    assert [h[-1] for h in histograms] == [24, 32, 0]

    assert len(pairs) == 768 and set(pairs.values()) == {12}
    assert sum(recurrent.values()) == 640
    assert set(drain_images) == {key for key, value in recurrent.items() if value}
    assert len(rotation_groups) == 96
    assert sum(len(g) == 6 for g in rotation_groups.values()) == 64
    assert sum(len(g) == 12 for g in rotation_groups.values()) == 32
    for representative, members in rotation_groups.items():
        orbit = {tuple(row[e:]+row[:e] for row, e in zip(representative, shifts))
                 for shifts in product(*(range(m) for m in lengths))}
        assert orbit == set(members)
        assert all(sum(d[v] for d in members.values()) >= 0 for v in range(3))
    rec_sectors, cut_checks = defaultdict(int), 0
    for (arrows, bits, C), is_rec in recurrent.items():
        rec_sectors[tuple(map(sum, bits))] += is_rec
        if is_rec:
            heads = drain_residents(arrows, bits, C)
            recovered = tuple(tuple(row[e:]+row[:e] for row, e in zip(rows, heads))
                              for rows in (arrows, bits))
            assert recovered == drain_images[arrows, bits, C]
        reverse_arrows = tuple(row[::-1] for row in arrows)
        reverse_bits = tuple(tuple(1-b for b in row[::-1]) for row in bits)
        S = tuple(1-c for c in C)
        greatest, _ = fixed_cut(reverse_arrows, reverse_bits, S, greatest=True)
        assert is_rec == (greatest == (0, 0, 0))
        for used in product(*(range(m+1) for m in lengths)):
            n = arrivals(arrows, used, tuple(1+c for c in C))
            B = tuple(a-u for a, u in zip(n, used))
            stable = all(b in (0, 1) and bits[i][(n[i]-1) % lengths[i]] == b
                         for i, b in enumerate(B))
            other = tuple(m-u for m, u in zip(lengths, used))
            reverse_stable = cut_map(reverse_arrows, reverse_bits, S, other) == other
            assert stable == reverse_stable
            if stable:
                rn = arrivals(reverse_arrows, other, S)
                assert tuple(a-u for a, u in zip(rn, other)) == tuple(1-b for b in B)
            cut_checks += 1
    assert all(12*rec_sectors[k] == 10*number for k, number in population.items())
    return {"decorated_words": 768, "coin_sectors": 36,
            "physical_signed_totals": [8, 336, 112],
            "phase_burning_tests": phase_checks, "recurrent_phases_per_word": 10,
            "compatible_pairs": len(pairs), "recurrent_pairs": 640,
            "cut_reversal_checks": cut_checks,
            "buffer_rotor_intertwining_checks": rotor_checks,
            "recurrent_forest_bijections": 768,
            "one_full_input_inverse_checks": 7680,
            "drain_recurrent_pair_bijection": 640,
            "nonnegative_physical_rotation_orbits": 96}


def no_cut_ideal_audit():
    """A positive-probability field whose no-cut ideal is not a matroid."""
    arrows = ((1, 3, 1), (2, 3, 2, 3, 2), (0, 1, 0, 1))
    lengths = tuple(map(len, arrows))
    good, nonzero, checks = set(), {}, 0
    for eta in product((0, 1), repeat=3):
        C = frozenset(i for i, value in enumerate(eta) if value)
        cuts = []
        for m in product(*(range(0, length+1, 2) for length in lengths)):
            checks += 1
            arrivals = list(eta)
            for row, length in zip(arrows, m):
                for dest in row[:length]:
                    if dest < 3:
                        arrivals[dest] += 1
            assert all(a <= length for a, length in zip(arrivals, lengths))
            if tuple(a-a % 2 for a in arrivals) == m:
                cuts.append(m)
        assert (0, 0, 0) in cuts
        positive = [m for m in cuts if any(m)]
        if positive:
            assert len(positive) == 1
            nonzero[C] = positive[0]
        else:
            good.add(C)
    assert good == {frozenset(), frozenset({0}), frozenset({1}),
                    frozenset({2}), frozenset({0, 1})}
    assert nonzero == {frozenset({0, 2}): (2, 2, 2),
                       frozenset({1, 2}): (0, 2, 2),
                       frozenset({0, 1, 2}): (2, 2, 2)}
    assert frozenset({2}) in good and frozenset({0, 1}) in good
    assert all(frozenset({2, i}) not in good for i in (0, 1))
    bits = ((1, 0, 1), (1, 0, 1, 0, 1), (1, 0, 1, 0))
    base = tuple(tuple(zip(a, b)) for a, b in zip(arrows, bits))
    orbit = {tuple(row[e:]+row[:e] for row, e in zip(base, shifts))
             for shifts in product(*(range(m) for m in lengths))}
    assert len(orbit) == 30
    totals = {eta: [0, 0, 0] for eta in product((0, 1), repeat=3)}
    for word in orbit:
        for eta in totals:
            used = (0, 0, 0)
            for _ in range(sum(lengths)+1):
                incoming = list(eta)
                for row, u in zip(word, used):
                    for dest, bit in row[:u]:
                        if dest < 3:
                            incoming[dest] += 1
                assert all(0 <= a <= m for a, m in zip(incoming, lengths))
                endpoint = tuple(word[i][a-1][1] if a else 0
                                 for i, a in enumerate(incoming))
                nxt = tuple(a-b for a, b in zip(incoming, endpoint))
                assert all(a <= b for a, b in zip(used, nxt))
                if used == nxt:
                    for v in range(3):
                        totals[eta][v] += endpoint[v]
                    break
                used = nxt
            else:
                raise AssertionError("rotated nonmatroid book did not stabilize")
    for (c0, c1, c2), total in totals.items():
        assert total == [20*c0, 18*c1+9*(1-c1)*c2, (15+9*c1)*c2]
    covers = 0
    for eta, low in totals.items():
        for x in range(3):
            if eta[x] == 0:
                upper = tuple(1 if i == x else eta[i] for i in range(3))
                assert all(a <= b for a, b in zip(low, totals[upper]))
                covers += 1
    assert covers == 12
    return {"even_cut_candidates": checks,
            "no_cut_ideal": sorted(map(sorted, good)),
            "matroid_exchange": "fails", "paired_rotation_orbit": 30,
            "physical_stabilizations": 240, "nonnegative_cover_target_sums": 36}


def rotation_obstruction_audit():
    """J110: paired rotations fail, but all orders repair one cycle family."""
    from itertools import combinations, permutations
    from math import comb

    def run(arrows, bits, C):
        n = len(arrows)
        lengths = tuple(map(len, arrows))
        # Independent legal arrival-buffer execution.
        queue = [int(i in C) for i in range(n)]
        arrivals, emitted, resident = [0]*n, [0]*n, [0]*n
        for _ in range(sum(lengths)+1):
            if not any(queue):
                break
            i = next(j for j in range(n) if queue[j])
            queue[i] -= 1
            arrivals[i] += 1
            assert arrivals[i] <= lengths[i]
            new = bits[i][arrivals[i]-1]
            out = 1+resident[i]-new
            for _ in range(out):
                assert emitted[i] < lengths[i]
                dest = arrows[i][emitted[i]]
                emitted[i] += 1
                if dest < n:
                    queue[dest] += 1
            resident[i] = new
        assert not any(queue)
        actual = tuple(emitted), tuple(resident)
        # Monotone least-fixed-point iteration has no legal event schedule.
        cut = (0,)*n
        for _ in range(sum(lengths)+1):
            incoming = [int(i in C) for i in range(n)]
            for row, used in zip(arrows, cut):
                for dest in row[:used]:
                    if dest < n:
                        incoming[dest] += 1
            assert all(x <= m for x, m in zip(incoming, lengths))
            endpoint = tuple(bits[i][x-1] if x else 0
                             for i, x in enumerate(incoming))
            nxt = tuple(x-b for x, b in zip(incoming, endpoint))
            assert all(a <= b for a, b in zip(cut, nxt))
            if nxt == cut:
                assert actual == (cut, endpoint)
                return actual
            cut = nxt
        raise AssertionError("unit-divergence upper bound did not terminate")

    def rotate(row, j):
        return row[j:]+row[:j]

    arrows, bits = ((0, 2, 2), (0,)), ((0, 1, 0), (0,))
    rows = []
    for j in range(3):
        a, b = (rotate(arrows[0], j), arrows[1]), (rotate(bits[0], j), bits[1])
        rows.append((run(a, b, {0}), run(a, b, {0, 1})))
    assert rows == [(((1, 0), (1, 0)), ((3, 1), (0, 0))),
                    (((0, 0), (1, 0)), ((2, 1), (0, 0))),
                    (((1, 0), (0, 0)), ((2, 1), (1, 0)))]
    totals = [0, 0]
    for t, b in product(range(3), repeat=2):
        a = (tuple(0 if j == t else 2 for j in range(3)), (0,))
        beta = (tuple(int(j == b) for j in range(3)), (0,))
        for index, C in enumerate(({0}, {0, 1})):
            totals[index] += run(a, beta, C)[1][0]
    assert totals == [4, 4]

    n, L = 4, 5
    arrows = ((1,)+(n,)*n,)+tuple(((i+1) % n,)*(i+1) for i in range(1, n))
    tails = tuple((0,)*(i+1) for i in range(1, n))
    sets = supports(n)
    orbit_checks = order_checks = 0
    for d in (1, 2, 3):
        row = tuple(int(j == d) for j in range(L))
        totals = {C: 0 for C in sets}
        for shift in range(L):
            a = (rotate(arrows[0], shift),)+arrows[1:]
            b = (rotate(row, shift),)+tails
            for C in sets:
                endpoint = run(a, b, C)[1]
                assert not any(endpoint[1:])
                totals[C] += endpoint[0]
                orbit_checks += 1
        assert all(totals[C] == (1+int(len(C) == d) if C else 0) for C in sets)
    for k in range(L+1):
        totals = {C: 0 for C in sets}
        for t in range(L):
            a = (tuple(1 if j == t else n for j in range(L)),)+arrows[1:]
            for ones in combinations(range(L), k):
                b = (tuple(int(j in ones) for j in range(L)),)+tails
                for C in sets:
                    endpoint = run(a, b, C)[1]
                    assert not any(endpoint[1:])
                    totals[C] += endpoint[0]
                    order_checks += 1
        count = L*comb(L, k)
        prediction = F(k, L)+F(k*(L-k), L*L*(L-1))
        assert all(F(totals[C], count) == (prediction if C else 0) for C in sets)
    # J111: conditioning rotations on the recurrent rotor heads also fails.
    arrows = ((4, 1, 3), (2, 0, 4), (3, 1, 4), (0, 2, 4))
    bits = ((1, 1, 1), (0, 1, 0), (1, 0, 0), (1, 1, 1))

    def forest(a):
        for start in range(n):
            seen, vertex = set(), start
            while vertex < n and vertex not in seen:
                seen.add(vertex)
                vertex = a[vertex][-1]
            if vertex < n:
                return False
        return True

    totals = [[0]*n for _ in range(n+1)]
    count = 0
    for shifts in product(range(3), repeat=n):
        a = tuple(rotate(row, j) for row, j in zip(arrows, shifts))
        b = tuple(rotate(row, j) for row, j in zip(bits, shifts))
        if not forest(a):
            continue
        count += 1
        for missing in range(-1, n):
            endpoint = run(a, b, set(range(n))-{missing})[1]
            totals[missing+1] = [x+y for x, y in zip(totals[missing+1], endpoint)]
    assert count == 45
    assert totals == [[45, 20, 21, 45], [7, 21, 21, 45], [45, 4, 21, 45],
                      [45, 17, 8, 45], [45, 21, 21, 9]]

    def pool(options):
        count = upper = lower = 0
        for rows in product(*options):
            a = tuple(tuple(pair[0] for pair in row) for row in rows)
            if not forest(a):
                continue
            b = tuple(tuple(pair[1] for pair in row) for row in rows)
            count += 1
            upper += run(a, b, {0, 1, 2, 3})[1][1]
            lower += run(a, b, {1, 2, 3})[1][1]
        return count, upper, lower

    paired = pool([tuple(set(permutations(tuple(zip(a, b)))))
                   for a, b in zip(arrows, bits)])
    independent = pool([tuple(tuple(zip(aa, bb))
                              for aa in set(permutations(a))
                              for bb in set(permutations(b)))
                        for a, b in zip(arrows, bits)])
    assert paired == (720, 332, 326)
    assert independent == (6480, 3004, 2512)
    return {"rotation_increment": "-1/3", "minimal_all_order_totals": [4, 4],
            "four_cycle_rotation_runs": orbit_checks,
            "four_cycle_all_order_runs": order_checks,
            "recurrent_rotation_increment": "-1/45",
            "recurrent_paired_pool": paired, "recurrent_all_order_pool": independent}


def cycle_remainder_geometry_audit():
    """Four-cycle height determinants, and two SCCs with one-way routing."""
    from math import prod

    def blocks(K, p):
        sets, X = homogeneous(K)
        size = len(sets)
        W = [[[X[i][a][b]-(1-p[i])*sum(K[i][j]*X[j][a][b]
                                              for j in range(len(K)))
               for b in range(size)] for a in range(size)] for i in range(len(K))]
        rows = []
        for C in sets:
            row = [F(not S) for S in sets]
            for i in C:
                row = multiply([row], W[i])[0]
            rows.append(row)
        return sets, [[[rows[i][j] for j, B in enumerate(sets) if len(B) == k]
                       for i, C in enumerate(sets) if len(C) == k]
                      for k in range(len(K)+1)]

    forward = (F(1, 2), F(2, 3), F(3, 4), F(4, 5))
    K = [[forward[i] if j == (i+1) % 4 else F(0) for j in range(4)]
         + [1-forward[i]] for i in range(4)]
    p = (F(1, 3), F(2, 5), F(3, 7), F(4, 9))
    sets, Bp = blocks(K, p)
    _, Bq = blocks(K, [1-x for x in p])
    _, Bzero = blocks(K, [F(0)]*4)
    kappa, d = prod(forward), 1-prod(forward)
    assert [determinant(B) for B in Bzero] == [F(1), d, d, d, d]
    sector_determinants = []
    target_sign_entries = 0
    for k in range(5):
        labels = [C for C in sets if len(C) == k]
        H = []
        for C in labels:
            row = []
            for B in labels:
                height, running = [], 0
                for i in range(4):
                    running += int(i in C)-int(i in B)
                    height.append(running)
                low = min(height)
                height = [h-low for h in height]
                assert all(height[i]-height[i-1] == int(i in C)-int(i in B)
                           for i in range(4))
                row.append(prod(a**h for a, h in zip(forward, height)))
            H.append(row)
        assert determinant(H) == d**(len(labels)-1)
        masses = [prod(p[i] if i in B else 1-p[i] for i in range(4)) for B in labels]
        Rk = [[F(i == j)+kappa/d*H[i][j]*masses[j]
               for j in range(len(labels))] for i in range(len(labels))]
        sector = determinant(Rk)
        assert sector == determinant(Bp[k])*determinant(Bq[4-k])/d
        sector_determinants.append(str(sector))
        if 0 < k < 4:
            assert inverse(H) == [[value/d for value in row] for row in Bzero[k]]
            Rinv = inverse(Rk)
            assert all(Rinv[i][j] <= 0 for i in range(len(labels))
                       for j in range(len(labels)) if i != j)
            O = [[F(i == j and 0 in C) for j in range(len(labels))]
                 for i, C in enumerate(labels)]
            target = multiply(multiply(Rinv, O), Rk)
            for i, C in enumerate(labels):
                for j in range(len(labels)):
                    assert (2*int(0 in C)-1)*(target[i][j]-O[i][j]) >= 0
                    target_sign_entries += 1

    # A directed edge between the SCCs affects off-diagonal blocks only.
    coupled = [[F(0), F(1, 2), F(0), F(0), F(1, 2)],
               [F(1, 4), F(0), F(1, 4), F(0), F(1, 2)],
               [F(0), F(0), F(0), F(1, 2), F(1, 2)],
               [F(0), F(0), F(1, 2), F(0), F(1, 2)]]
    _, global_blocks = blocks(coupled, p)
    factors = []
    for S in ((0, 1), (2, 3)):
        local = [[coupled[i][j] for j in S]+[1-sum(coupled[i][j] for j in S)]
                 for i in S]
        _, local_blocks = blocks(local, [p[i] for i in S])
        factors.append(prod(determinant(B) for B in local_blocks))
    assert prod(determinant(B) for B in global_blocks) == prod(a**4 for a in factors)
    return {"cycle_particle_sectors": 5, "cycle_A_at_zero": str(d**4),
            "cycle_sector_determinants": sector_determinants,
            "cycle_target_sign_entries": target_sign_entries,
            "SCC_factorization": "PASS"}


def polytope_audit():
    # Same single-leak cycle as audit(), now its homogeneous killed algebra.
    K = [[F(0), F(1), F(0), F(0)],
         [F(0), F(0), F(1), F(0)],
         [F(1, 2), F(0), F(0), F(1, 2)]]
    sets, U = homogeneous(K)
    n, dim = len(K), len(sets)
    M = [[F(i == j)-K[i][j] for j in range(n)] for i in range(n)]
    V = [[[sum(M[x][y]*U[y][i][j] for y in range(n))
           for j in range(dim)] for i in range(dim)] for x in range(n)]
    det = determinant(M)
    for C in sets:
        principal = [[M[i][j] for j in sorted(C)] for i in sorted(C)]
        assert top_product(sets, [U[i] if i in C else V[i]
                                  for i in range(n)]) == det/determinant(principal)
    # Nine Alexandrov--Fenchel inequalities for three-particle input laws.
    for a in range(n):
        for i in range(n):
            for j in range(i+1, n):
                mixed = top_product(sets, [U[a], U[i], U[j]])
                left = top_product(sets, [U[a], U[i], U[i]])
                right = top_product(sets, [U[a], U[j], U[j]])
                assert mixed*mixed >= left*right

    # The asymmetric two-site kernel from the current-tilt certificate,
    # with its two exit labels combined. Independent polygon area check.
    K2 = [[F(1, 5), F(2, 5), F(2, 5)],
          [F(1, 2), F(0), F(1, 2)]]
    sets2, U2 = homogeneous(K2)
    M2 = [[F(i == j)-K2[i][j] for j in range(2)] for i in range(2)]
    b = [F(2, 3), F(5, 4)]
    MT_inverse = inverse([[M2[j][i] for j in range(2)] for i in range(2)])
    full = [sum(row[j]*b[j] for j in range(2)) for row in MT_inverse]
    vertices = [[F(0), F(0)], [b[0]/M2[0][0], F(0)],
                full, [F(0), b[1]/M2[1][1]]]
    area = sum(x[0]*y[1]-x[1]*y[0]
               for x, y in zip(vertices, vertices[1:]+vertices[:1]))/2
    ell = [[sum(b[x]*U2[x][i][j] for x in range(2))
            for j in range(4)] for i in range(4)]
    top = top_product(sets2, [ell, ell])
    assert top == 2*determinant(M2)*area

    # Two-site path: an inverse sleep class can have negative square.
    path = [[F(0), F(1), F(0)], [F(0), F(0), F(1)]]
    paths, path_U = homogeneous(path)
    for odds in (F(1, 2), F(1), F(2)):
        w = [[(1+odds)*path_U[0][i][j]-odds*path_U[1][i][j]
              for j in range(4)] for i in range(4)]
        assert top_product(paths, [w, w]) == 1-odds*odds
    return {"facet_products": 8, "filling_AF_inequalities": 9,
            "independent_polygon_area": str(area), "top_square": str(top),
            "path_inverse_class_square": "1-r^2"}


def audit():
    K = [[F(0), F(1), F(0), F(0)],
         [F(0), F(0), F(1), F(0)],
         [F(1, 2), F(0), F(0), F(1, 2)]]
    n = len(K)
    sets = supports(n)
    pos = {B: i for i, B in enumerate(sets)}
    m = len(sets)
    identity = [[F(i == j) for j in range(m)] for i in range(m)]
    D = idla_addition_matrices(K, sets)
    for x in range(n):
        J = [[sum(K[x][y] * D[y][i][j] for y in range(n))
              + K[x][n] * identity[i][j] for j in range(m)] for i in range(m)]
        assert multiply(D[x], D[x]) == multiply(D[x], J)
        for y in range(n):
            assert multiply(D[x], D[y]) == multiply(D[y], D[x])
    for B in sets:
        row = [F(not C) for C in sets]
        for x in sorted(B):
            row = multiply([row], D[x])[0]
        assert row == [F(C == B) for C in sets]

    fibers = defaultdict(list)
    for B in sets:
        R = sink_core(K, B)
        fibers[R].append(B)
        assert tuple(D[x][pos[B]][pos[B]] for x in range(n)) == tuple(
            D[x][pos[R]][pos[R]] for x in range(n))
    expected = {frozenset(): 4, frozenset({2}): 2,
                frozenset({1, 2}): 1, frozenset(range(n)): 1}
    assert {R: len(fiber) for R, fiber in fibers.items()} == expected
    character_table = []
    for R, fiber in fibers.items():
        X = frozenset(range(n)) - R
        frontier = {x for x in X if K[x][n] + sum(K[x][y] for y in R)}
        Y = X - frontier
        assert set(fiber) == {B for B in sets if R <= B <= R | Y}
        character_table.append({"R": sorted(R), "deep_vertices": sorted(Y),
                                "multiplicity": len(fiber),
                                "character": [str(D[x][pos[R]][pos[R]]) for x in range(n)]})
    feasible = set(fibers)
    assert all(R | S in feasible for R in feasible for S in feasible)
    assert all(not R or any(R - {x} in feasible for x in R) for R in feasible)

    rho = [F(1, 2), F(1, 3), F(1, 6)]
    driven = [[sum(rho[x] * D[x][i][j] for x in range(n))
               for j in range(m)] for i in range(m)]
    checks = {F(0): [2, 3, 4, 4], F(1, 12): [1, 2, 2],
              F(1, 4): [1, 1], F(1): [1, 1]}
    for eigenvalue, expected_nullities in checks.items():
        shifted = [[x - eigenvalue * identity[i][j] for j, x in enumerate(row)]
                   for i, row in enumerate(driven)]
        power = identity
        for nullity in expected_nullities:
            power = multiply(power, shifted)
            assert m - rank(power) == nullity

    # The local factor at the zero character is the nonlosing IDLA
    # algebra on the two-site path: u0^2=u0*u1 and u1^2=0.
    local_K = [[F(0), F(1), F(0)], [F(0), F(0), F(1)]]
    local_sets = supports(2)
    local_D = idla_addition_matrices(local_K, local_sets)
    killed = [[[w if i != j else F(0) for j, w in enumerate(row)]
               for i, row in enumerate(a)] for a in local_D]
    zero = [[F(0)] * 4 for _ in range(4)]
    assert multiply(killed[0], killed[0]) == multiply(killed[0], killed[1])
    assert multiply(killed[1], killed[1]) == zero
    assert multiply(multiply(killed[0], killed[0]), killed[0]) == zero
    assert multiply(killed[0], killed[0]) != zero
    return {"status": "PASS", "characters": character_table,
            "driving_law": list(map(str, rho)),
            "Jordan_blocks": {"0": [3, 1], "1/12": [2], "1/4": [1], "1": [1]},
            "local_zero_character_relations": ["u0^2=u0*u1", "u1^2=0"],
            "polytope_checks": polytope_audit(),
            "capacity_checks": capacity_audit(),
            "poisson_cells": [poisson_cell_audit(), poisson_cell_audit(F(2))],
            "complementary_duality": duality_audit(),
            "full_minus_one_product": full_minus_one_product_audit(),
            "ordered_bridge_books": bridge_book_audit(),
            "stable_upper_books": stable_upper_book_audit(),
            "stable_cut_lattice": stable_cut_lattice_audit(),
            "target_cut_cancellation": target_cut_cancellation_audit(),
            "target_inventory_cancellation": target_inventory_cancellation_audit(),
            "one_particle_mobius": one_particle_mobius_audit(),
            "greatest_cut_fields": greatest_cut_field_audit(),
            "envelope_inventory": envelope_inventory_audit(),
            "cyclic_envelope": cyclic_envelope_audit(),
            "stationary_lattice": stationary_lattice_audit(),
            "no_cut_ideal": no_cut_ideal_audit(),
            "rotation_obstruction": rotation_obstruction_audit(),
            "cycle_remainder_geometry": cycle_remainder_geometry_audit()}


if __name__ == "__main__":
    print(json.dumps(audit(), indent=2))
