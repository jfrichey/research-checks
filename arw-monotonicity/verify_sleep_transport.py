"""Exact audits of the Markov flow obtained by varying sleep odds.

The complete-block inverse gives the derivatives. Independent raw ARW
absorption gives the resident-release generator and waking transport.
Two fixed motifs only; no parameter search. Standard library only.
AI-written research code; general statements are proved in RN section 29.
"""

from collections import defaultdict
from fractions import Fraction as F
from itertools import product
from math import prod
import json
import sys

sys.dont_write_bytecode = True
import verify_response_renewal as rn
from verify_feedback_cone import binary_stabilization


def subtract(A, B):
    return [[a-b for a, b in zip(u, v)] for u, v in zip(A, B)]


def add(A, B):
    return [[a+b for a, b in zip(u, v)] for u, v in zip(A, B)]


def negate(A):
    return [[-a for a in row] for row in A]


def at_odds(K, r):
    return binary_stabilization(K, [1/(1+x) for x in r])


def raw_release_generators(K, r, supports):
    n = len(K)
    p = [1/(1+x) for x in r]
    starts = {}
    for C in supports:
        for x in C:
            law = defaultdict(F)
            for y, weight in enumerate(K[x]):
                if not weight:
                    continue
                state = [-int(i in C and i != x) for i in range(n)]
                if y < n:
                    state[y] = 2 if state[y] == -1 else 1
                law[tuple(state)] += weight
            starts[C, x] = law
    initial = {state for law in starts.values() for state in law}
    laws, count = rn.stopped(initial, lambda state: rn.arw_step(K, p, state),
                            lambda state: all(x <= 0 for x in state))
    terminal = {B: tuple(-int(i in B) for i in range(n)) for B in supports}
    G = []
    for x in range(n):
        matrix = []
        for C in supports:
            if x not in C:
                matrix.append([F(0)]*len(supports))
                continue
            released = [sum(w*laws[state].get(terminal[B], F(0))
                            for state, w in starts[C, x].items()) for B in supports]
            assert sum(released) == 1
            matrix.append([p[x]*(w-F(C == B)) for B, w in zip(supports, released)])
        G.append(matrix)
    return G, count


def raw_waking_transport(K, r, s, supports):
    """Wake each current sleeper i with probability 1-r_i/s_i, then stabilize."""
    assert all(0 <= a <= b for a, b in zip(r, s))
    n = len(K)
    p = [1/(1+x) for x in s]
    mixtures = {}
    for C in supports:
        law = defaultdict(F)
        for bits in product((0, 1), repeat=len(C)):
            active = {i for i, b in zip(sorted(C), bits) if b}
            weight = prod((1-r[i]/s[i] if i in active else r[i]/s[i])
                          if s[i] else F(i not in active) for i in C)
            if weight:
                state = tuple(1 if i in active else -int(i in C) for i in range(n))
                law[state] += weight
        assert sum(law.values()) == 1
        mixtures[C] = law
    starts = {state for law in mixtures.values() for state in law}
    laws, count = rn.stopped(starts, lambda state: rn.arw_step(K, p, state),
                            lambda state: all(x <= 0 for x in state))
    terminal = {B: tuple(-int(i in B) for i in range(n)) for B in supports}
    T = [[sum(w*laws[state].get(terminal[B], F(0))
              for state, w in mixtures[C].items()) for B in supports] for C in supports]
    return T, count


def flow_audit(name, K, r):
    n = len(K)
    supports, P, E, D = at_odds(K, r)
    m = len(supports)
    zero = [[F(0)]*m for _ in range(m)]
    identity = [[F(i == j) for j in range(m)] for i in range(m)]
    Gi, raw_states = raw_release_generators(K, r, supports)
    Eplus = []
    for x in range(n):
        new = list(r)
        new[x] += 1
        Eplus.append(at_odds(K, new)[2])
    dE = [subtract(ep, E) for ep in Eplus]
    G = [negate(rn.multiply(de, P)) for de in dE]
    assert G == Gi
    A = [rn.multiply(rn.multiply(E, d), P) for d in D]
    for x, generator in enumerate(G):
        assert all(sum(row) == 0 for row in generator)
        for i, C in enumerate(supports):
            for j, B in enumerate(supports):
                assert i == j or generator[i][j] >= 0
                assert len(B) <= len(C) or generator[i][j] == 0
                assert x in C or generator[i][j] == 0
        # Compare the infinitesimal Lax equation with differentiating E D_z P.
        dP = negate(rn.multiply(rn.multiply(P, dE[x]), P))
        for z in range(n):
            dA = add(rn.multiply(rn.multiply(dE[x], D[z]), P),
                     rn.multiply(rn.multiply(E, D[z]), dP))
            assert dA == subtract(rn.multiply(A[z], G[x]), rn.multiply(G[x], A[z]))
    curvature_checks = 0
    for x in range(n):
        for y in range(x):
            new = list(r)
            new[x] += 1
            new[y] += 1
            mixed = add(subtract(subtract(at_odds(K, new)[2], Eplus[x]), Eplus[y]), E)
            dxy = add(negate(rn.multiply(mixed, P)),
                      rn.multiply(rn.multiply(rn.multiply(dE[y], P), dE[x]), P))
            dyx = add(negate(rn.multiply(mixed, P)),
                      rn.multiply(rn.multiply(rn.multiply(dE[x], P), dE[y]), P))
            comm = subtract(rn.multiply(G[x], G[y]), rn.multiply(G[y], G[x]))
            assert add(subtract(dxy, dyx), comm) == zero
            curvature_checks += m*m

    s = [x+1 for x in r]
    paths = [list(range(n)), list(reversed(range(n)))]
    path_states = 0
    direct, count = raw_waking_transport(K, r, s, supports)
    path_states += count
    assert direct == rn.multiply(E, at_odds(K, s)[1])
    for order in paths:
        current = list(r)
        composed = identity
        for x in order:
            following = list(current)
            following[x] = s[x]
            step, count = raw_waking_transport(K, current, following, supports)
            path_states += count
            composed = rn.multiply(composed, step)
            current = following
        assert composed == direct
    from_zero, count = raw_waking_transport(K, [F(0)]*n, r, supports)
    assert from_zero == P
    path_states += count
    result = dict(case=name, raw_release_states=raw_states,
                  raw_transport_states=path_states,
                  release_generator_entries=n*m*m,
                  curvature_entries=curvature_checks,
                  physical_path_comparisons=2,
                  raw_zero_odds_stabilization_entries=m*m)
    if name == 'two-half':
        comm = subtract(rn.multiply(G[0], G[1]), rn.multiply(G[1], G[0]))
        value = comm[supports.index(frozenset((0,)))][supports.index(frozenset())]
        assert value == F(2, 45)
        result['frozen_commutator_0_to_empty'] = str(value)
        small = supports.index(frozenset((0,)))
        large = supports.index(frozenset((0, 1)))
        drift = sum(G[1][large][j]-G[1][small][j]
                    for j, B in enumerate(supports) if 0 in B)
        assert drift == F(-22, 195)
        result['site_flow_nonattractive_drift'] = str(drift)
    return result


def inverse_odds_series(K, degree):
    """Coefficients of the finite polynomial E_(epsilon,...,epsilon)."""
    n = len(K)
    supports, _, _, D = at_odds(K, [F(0)]*n)
    m = len(supports)
    coefficient = [[[F(0)]*m for _ in range(m)] for _ in range(degree+1)]
    for i, C in enumerate(supports):
        row = [[F(not B) for B in supports]]+[[F(0)]*m for _ in range(degree)]
        for x in sorted(C):
            jump = [[K[x][n]*F(j == k)+sum(K[x][y]*D[y][j][k] for y in range(n))
                     for k in range(m)] for j in range(m)]
            change = subtract(D[x], jump)
            row = [[a+b for a, b in zip(
                rn.multiply([row[ell]], D[x])[0],
                rn.multiply([row[ell-1]], change)[0] if ell else [F(0)]*m)]
                for ell in range(degree+1)]
        for ell in range(degree+1):
            coefficient[ell][i] = row[ell]
    assert coefficient[0] == [[F(i == j) for j in range(m)] for i in range(m)]
    return supports, coefficient


def logarithm_counterexample():
    """Exact analytic certificate: stochastic quarter-powers need not be monotone."""
    K = [[F(0), F(1), F(0), F(0)],
         [F(0), F(0), F(1), F(0)],
         [F(1, 2), F(0), F(0), F(1, 2)]]
    degree = 4
    supports, E = inverse_odds_series(K, degree)
    m = len(supports)
    zero = lambda: [[F(0)]*m for _ in range(m)]
    H = [zero()]+E[1:]
    logP = [zero() for _ in range(degree+1)]
    power = [E[0]]+[zero() for _ in range(degree)]
    for k in range(1, degree+1):
        result = [zero() for _ in range(degree+1)]
        for d in range(1, degree+1):
            for j in range(1, d+1):
                result[d] = add(result[d], rn.multiply(power[d-j], H[j]))
        power = result
        for d in range(1, degree+1):
            logP[d] = add(logP[d], [[F((-1)**k, k)*q for q in row] for row in power[d]])

    # Independent raw absorption checks the whole inverse polynomial at one point.
    epsilon = F(1, 2)
    rawP, count = raw_waking_transport(K, [F(0)]*3, [epsilon]*3, supports)
    exactE = [[sum(epsilon**d*E[d][i][j] for d in range(degree+1))
               for j in range(m)] for i in range(m)]
    assert rn.multiply(exactE, rawP) == E[0]

    leading = defaultdict(list)
    for i, A in enumerate(supports):
        for j, B in enumerate(supports):
            if i == j:
                continue
            if len(B) > len(A):
                assert all(logP[d][i][j] == 0 for d in range(degree+1))
                continue
            d, value = next((d, logP[d][i][j]) for d in range(1, degree+1)
                            if logP[d][i][j])
            assert value > 0
            leading[d].append(value)
    assert {d: len(v) for d, v in leading.items()} == {1: 16, 2: 13, 3: 5}
    assert {d: min(v) for d, v in leading.items()} == {1: F(1, 2), 2: F(1, 8), 3: F(1, 24)}
    assert all(sum(row) == 0 for matrix in logP for row in matrix)
    a = supports.index(frozenset((2,)))
    b = supports.index(frozenset((1, 2)))
    score = lambda M: sum(M[b][j]-M[a][j] for j, C in enumerate(supports) if 2 in C)
    assert [score(logP[d]) for d in range(1, 5)] == [F(0), F(-1, 4), F(1, 6), F(-1, 24)]
    assert score(rn.multiply(logP[1], logP[1]))/2 == F(3, 4)
    t = F(1, 4)
    quarter = t*score(logP[2])+t*t*score(rn.multiply(logP[1], logP[1]))/2
    assert quarter == F(-1, 64)
    return dict(raw_inverse_audit_states=count,
                positive_log_leading_counts={d: len(v) for d, v in leading.items()},
                minimum_leading_coefficients={d: str(min(v)) for d, v in leading.items()},
                log_drift_coefficients=[str(score(logP[d])) for d in range(1, 5)],
                fractional_time_quadratic='t*(3*t-1)/4',
                quarter_time_quadratic=str(quarter),
                scope='analytic example for sufficiently small positive uniform jump odds')


def release_resolvent_audit():
    """J113: raw release, geometric cancellation, and forest coefficients."""
    K = [[F(0), F(1, 2), F(1, 2)], [F(1, 2), F(0), F(1, 2)]]
    sets, P0, E0, D = at_odds(K, [F(0), F(1)])
    pos = {C: j for j, C in enumerate(sets)}
    dim = len(sets)
    Id = [[F(a == b) for b in range(dim)] for a in range(dim)]
    # At p_0=1 use the continuous boundary law: collapse failed sleeps
    # when two particles meet into their compulsory next jump.
    sleeps = [F(1), F(1, 2)]

    def step(state):
        i = next(j for j, count in enumerate(state) if count > 0)
        out = defaultdict(F)
        jump_weight = F(1)
        if state[i] == 1:
            following = list(state)
            following[i] = -1
            out[tuple(following)] += sleeps[i]
            jump_weight = 1-sleeps[i]
        if jump_weight:
            for dest, weight in enumerate(K[i]):
                if not weight:
                    continue
                following = list(state)
                following[i] -= 1
                if dest < 2:
                    following[dest] = 2 if following[dest] == -1 else following[dest]+1
                out[tuple(following)] += jump_weight*weight
        return dict(out)

    initial = {}
    for C in sets:
        if 0 not in C:
            continue
        law = defaultdict(F)
        for dest, weight in enumerate(K[0]):
            if not weight:
                continue
            state = [0, -int(1 in C)]
            if dest < 2:
                state[dest] = 2 if state[dest] == -1 else 1
            law[tuple(state)] += weight
        initial[C] = law
    laws, states = rn.stopped({state for law in initial.values() for state in law},
                              step, lambda state: all(x <= 0 for x in state))
    W = []
    for C in sets:
        if 0 not in C:
            W.append(Id[pos[C]])
        else:
            W.append([sum(weight*laws[state].get(tuple(-int(j in B) for j in range(2)), F(0))
                          for state, weight in initial[C].items()) for B in sets])
    assert W == [[F(1), F(0), F(0), F(0)],
                 [F(0), F(1), F(0), F(0)],
                 [F(5, 8), F(1, 4), F(1, 8), F(0)],
                 [F(1, 16), F(5, 8), F(19, 112), F(1, 7)]]
    L = subtract(Id, W)
    A = rn.multiply(rn.multiply(E0, D[0]), P0)
    remove = [[F(B == C-{0}) for B in sets] for C in sets]
    projection = rn.multiply(A, remove)
    assert rn.multiply(W, A) == A
    assert rn.multiply(W, projection) == projection
    assert rn.multiply(projection, projection) == projection
    assert all(projection[a][b] == 0 for a in range(dim)
               for b, B in enumerate(sets) if 0 in B)
    for a, C in enumerate(sets):
        if 0 not in C:
            assert projection[a] == Id[a]

    low, high = pos[frozenset({1})], pos[frozenset({0, 1})]
    delta = [[a-b for a, b in zip(P0[high], P0[low])]]
    target = [[F(1 in C)] for C in sets]
    one = rn.multiply(rn.multiply(P0, W), target)
    assert one[low][0] == F(9, 16) and one[high][0] == F(107, 196)
    power, cumulative = Id, F(0)
    for k in range(13):
        gain = rn.multiply(rn.multiply(delta, power), target)[0][0]
        assert gain == F(3, 2)*F(1, 8)**k-F(10, 7)*F(1, 7)**k
        assert gain >= 0 if k == 0 else gain < 0
        cumulative += gain
        assert cumulative >= F(1, 21)
        power = rn.multiply(power, W)

    coefficients, forests = [F(1)], [Id]
    for k in (1, 2, 3):
        LH = rn.multiply(L, forests[-1])
        dk = sum(LH[j][j] for j in range(dim))/k
        Hk = subtract([[dk*x for x in row] for row in Id], LH)
        coefficients.append(dk)
        forests.append(Hk)
        assert min(map(min, Hk)) >= 0
        assert all(sum(row) == dk for row in Hk)
    assert coefficients == [F(1), F(97, 56), F(3, 4), F(0)]
    assert not any(x for row in forests[-1] for x in row)
    assert forests[2] == [[F(3, 4)*x for x in row] for row in projection]
    numerators = [rn.multiply(rn.multiply(delta, H), target)[0][0] for H in forests]
    assert numerators == [F(1, 14), F(1, 28), F(0), F(0)]
    for s in (F(1, 3), F(1), F(7)):
        resolvent = rn.inverse(add(Id, [[s*x for x in row] for row in L]))
        actual = at_odds(K, [s, F(1)])[1]
        assert actual == rn.multiply(P0, resolvent)
        gain = rn.multiply(rn.multiply(delta, resolvent), target)[0][0]
        assert gain == 2*(s+2)/((6*s+7)*(7*s+8))
        den = sum(c*s**k for k, c in enumerate(coefficients))
        adj = [[sum(H[a][b]*s**k for k, H in enumerate(forests))
                for b in range(dim)] for a in range(dim)]
        assert adj == [[den*x for x in row] for row in resolvent]
    return {"raw_release_states": states, "fixed_time_increment": "-13/784",
            "cumulative_limit": "1/21", "resolvent_matrix_checks": 3,
            "determinant_coefficients": list(map(str, coefficients[:-1])),
            "increment_numerator": list(map(str, numerators[:-1]))}


def integrated_release_audit():
    """J114: two-site formulas and a complete three-site generating identity."""
    def release(K, p, i=0):
        sets, P, E, D = binary_stabilization(K, p)
        dim, n = len(sets), len(K)
        pos = {C: j for j, C in enumerate(sets)}
        Jd = [[K[i][n]*F(a == b)+sum(K[i][j]*D[j][a][b] for j in range(n))
               for b in range(dim)] for a in range(dim)]
        J = rn.multiply(rn.multiply(E, Jd), P)
        W = [J[pos[C-{i}]] if i in C else [F(C == B) for B in sets] for C in sets]
        return sets, pos, P, W

    checks = 0
    for a, b, p in ((F(1, 2), F(1, 2), F(1, 2)),
                    (F(1), F(1, 2), F(1, 2)),
                    (F(1, 5), F(9, 10), F(3, 10)),
                    (F(0), F(1), F(7, 10)),
                    (F(1), F(0), F(2, 5))):
        K = [[F(0), a, 1-a], [b, F(0), 1-b]]
        sets, pos, P, W = release(K, [F(1), p])
        rho = a*b*(1-p)
        tau = p/(1-rho)
        sigma = a*b*tau
        c = a*b*(1-tau+(1-b)*(1-p))
        e = 1-a+a*(1-b)*p
        assert W[pos[frozenset({0})]] == [1-rho-a*p, a*p, rho, F(0)]
        assert W[-1] == [a*(1-b)**2*(1-p), e, c, sigma]
        delta = [[x-y for x, y in zip(P[-1], P[pos[frozenset({0})]])]]
        assert delta[0] == [F(0), F(0), -tau, tau]
        power = [[F(i == j) for j in range(4)] for i in range(4)]
        prefix = [F(0), F(0)]
        # Positive series factors in J114, compared with actual matrix powers.
        geometric = [sum(rho**j*sigma**(k-j) for j in range(k+1)) for k in range(25)]
        h1 = [F(1)] + [1-a]*24
        h2 = [F(1)] + [1-a*b]*24
        numerator_series = [sum(h1[j]*h2[k-j] for j in range(k+1)) for k in range(25)]
        for M in range(25):
            row = rn.multiply(delta, power)[0]
            for v in range(2):
                prefix[v] += sum(row[j] for j, B in enumerate(sets) if v in B)
            predicted1 = tau*sum(numerator_series[j]*geometric[M-j] for j in range(M+1))
            predicted0 = tau*a*b*(1-b*(1-p))*sum(geometric[:M])
            assert prefix == [predicted0, predicted1]
            assert min(prefix) >= 0
            power = rn.multiply(power, W)
            checks += 2

    K = [[F(0), F(1), F(0), F(0)],
         [F(0), F(0), F(1), F(0)],
         [F(1, 2), F(0), F(0), F(1, 2)]]
    sets, pos, P, W = release(K, [F(1), F(1, 2), F(1, 2)])
    transient = [j for j, B in enumerate(sets) if 0 in B]
    absorbing = [j for j, B in enumerate(sets) if 0 not in B]
    Q = [[W[a][b] for b in transient] for a in transient]
    assert Q == [[F(1, 8), F(0), F(0), F(0)],
                 [F(15, 112), F(1, 7), F(1, 28), F(0)],
                 [F(15, 112), F(1, 7), F(9, 56), F(0)],
                 [F(139, 1120), F(13, 70), F(117, 560), F(1, 5)]]
    delta = [[P[-1][j]-P[pos[frozenset({0, 2})]][j] for j in transient]]
    target = [[F(2 in sets[j])] for j in transient]
    absorbed = [[sum(W[j][k] for k in absorbing if 2 in sets[k])] for j in transient]
    Id = [[F(a == b) for b in range(4)] for a in range(4)]
    H, determinant_coefficients = [Id], [F(1)]
    minusQ = negate(Q)
    for k in range(1, 5):
        product = rn.multiply(minusQ, H[-1])
        coefficient = sum(product[j][j] for j in range(4))/k
        determinant_coefficients.append(coefficient)
        H.append(subtract([[coefficient*x for x in row] for row in Id], product))
    assert not any(x for row in H[-1] for x in row)
    numerator = [F(0)]*5
    for k, matrix in enumerate(H[:-1]):
        row = rn.multiply(delta, matrix)
        numerator[k] += rn.multiply(row, target)[0][0]
        numerator[k+1] += rn.multiply(row, subtract(absorbed, target))[0][0]
    assert numerator == [F(16, 280), F(-31, 280), F(16, 280), F(-1, 280), F(0)]

    def polymul(a, b):
        result = [F(0)]*(len(a)+len(b)-1)
        for i, x in enumerate(a):
            for j, y in enumerate(b):
                result[i+j] += x*y
        return result

    assert determinant_coefficients == polymul(
        polymul([F(1), -F(1, 8)], [F(1), -F(17, 56), F(1, 56)]),
        [F(1), -F(1, 5)])
    return {"two_site_parameter_cases": 5, "exact_prefix_checks": checks,
            "three_site_numerator": list(map(str, numerator)),
            "three_site_denominator": list(map(str, determinant_coefficients))}


def boundary_interpolation_audit():
    """Positive boundary factorization and an independent raw ARW audit."""
    import verify_release_forest as rf

    def add(*matrices):
        return [[sum(entries, F(0)) for entries in zip(*rows)]
                for rows in zip(*matrices)]

    def mat(*matrices):
        result = matrices[0]
        for right in matrices[1:]:
            result = rf.multiply(result, right)
        return result

    def blocks(a, b, c, d):
        return [x+y for x,y in zip(a,b)] + [x+y for x,y in zip(c,d)]

    K = [[F(1, 2), F(1, 2), 0, 0],
         [0, 0, 1, 0],
         [F(1, 2), 0, 0, F(1, 2)]]
    i, x, s = 0, 1, F(2)
    alg = rf.AdditionAlgebra(K)
    data = rf.ReleaseForest(alg, (0, 1, 1), i, degrees=False)
    ids = data.vacant + data.present
    m, I, Z = data.m, rf.identity(data.m), rf.zeros(data.m)
    A = [mat(data.E0, D, data.P0) for D in alg.D]
    Ax = rf.block(A[x], ids, ids)
    Ai = rf.block(A[i], ids, ids)
    U = rf.block(A[x], data.vacant, data.vacant)
    V = rf.block(A[x], data.vacant, data.present)
    At = rf.block(A[x], data.present, data.present)
    B, Q, M, T = data.B, data.Q, data.A, data.T
    N = add(I, rf.scale(s, M))
    R = rf.inverse(N)
    L = rf.scale(s, mat(R, B))
    G = add(R, L)
    S = blocks(I, Z, rf.scale(-s, B), N)
    Sinv = blocks(I, Z, L, R)
    Ji = rf.block(mat(data.E0, alg.J[i], data.P0), ids, ids)
    assert Ai == blocks(Z, I, Z, T)
    assert Ax == blocks(U, V, Z, At)
    assert Ji == blocks(B, Q, Z, T)
    assert mat(S, Sinv) == rf.identity(2*m)
    assert mat(Ax, Ji) == mat(Ji, Ax)
    assert mat(Ax, Ai) == mat(Ai, Ax)
    assert At == add(U, mat(V, T))
    assert mat(At, T) == mat(T, At)
    assert B == mat(M, T)
    assert mat(U, B) == mat(B, U)
    assert add(mat(U, Q), mat(V, T)) == add(mat(B, V), mat(Q, At))
    assert mat(M, At) == add(mat(U, M), mat(B, V))
    assert mat(At, R) == add(mat(R, U),
        mat(R, add(mat(V, T), rf.scale(s, mat(B, V))), R))
    predicted_i = blocks(L, R, mat(T, L), mat(T, R))
    predicted_x = blocks(add(U, mat(V, L)), mat(V, R),
                         mat(V, T, L), add(U, mat(V, T, R)))
    assert predicted_i == mat(S, Ai, Sinv)
    assert predicted_x == mat(S, Ax, Sinv)
    assert rf.stochastic(predicted_i) and rf.stochastic(predicted_x)
    assert rf.stochastic(G)
    assert all(q >= 0 for matrix in (R, L) for row in matrix for q in row)
    embedding = [left + right for left, right in zip(L, R)]
    assert rf.stochastic(embedding)
    assert mat(At, embedding) == mat(embedding, predicted_x)
    assert mat(T, embedding) == mat(embedding, predicted_i)
    assert G == add(T, mat(R, rf.subtract(I, T)))
    assert mat(rf.subtract(At, I), G) == add(
        mat(T, rf.subtract(At, I)),
        mat(rf.subtract(At, I), R, rf.subtract(I, T)))

    # The signed partner is used ONLY as an algebraic identity.
    Up = rf.subtract(At, mat(T, V))
    assert Up == mat(rf.inverse(M), U, M)
    negatives = [(a, b, q) for a, row in enumerate(Up)
                 for b, q in enumerate(row) if q < 0]
    assert negatives
    assert T[2] == T[3]  # In particular T is singular.

    # Independently solve physical addition from all sleeping supports.
    starts = {}
    for z in (i, x):
        for a, support in enumerate(alg.sets):
            state = [-int(j in support) for j in range(alg.n)]
            state[z] = 2 if state[z] == -1 else 1
            starts[z, a] = tuple(state)
    active = [tuple(int(j in support) for j in range(alg.n))
              for support in alg.sets]
    raw0, count0 = rf.raw_absorption(
        list(starts.values()) + active, alg.K, [F(1), F(1, 2), F(1, 2)], alg.sets)
    raws, counts = rf.raw_absorption(
        list(starts.values()) + active, alg.K, [F(1, 3), F(1, 2), F(1, 2)], alg.sets)
    for z, prediction in ((i, predicted_i), (x, predicted_x)):
        raw_boundary = [raw0[starts[z, a]] for a in range(len(alg.sets))]
        assert raw_boundary == A[z]
        raw_matrix = [raws[starts[z, a]] for a in range(len(alg.sets))]
        assert rf.block(raw_matrix, ids, ids) == prediction
    assert [raw0[state] for state in active] == data.P0
    ordered_P0 = rf.block(data.P0, ids, ids)
    Ps = mat(ordered_P0, Sinv)
    assert Ps == rf.block([raws[state] for state in active], ids, ids)

    # Only one occupation comparison: C={0,2}, added x=1, target i=0.
    # Evaluate the signed forcing vector on every outside state, not events.
    ones = [[F(1)] for _ in range(m)]
    forcing = mat(rf.subtract(mat(B, V, R), mat(V, T, R, B)), ones)
    root_gain_vector = mat(rf.subtract(At, I), R, ones)
    assert root_gain_vector == rf.scale(s, mat(R, forcing))
    C = frozenset({0, 2})
    nu = [data.P0[alg.pos[C]][j] for j in data.present]
    gain = mat([nu], root_gain_vector)[0][0]
    raw_gain = sum(raws[active[alg.pos[C | {x}]]][j] -
                   raws[active[alg.pos[C]]][j] for j in data.present)
    assert gain == raw_gain
    # Explicitly expose where noncommutation matters.
    commutators = {
        "TR-RT": rf.subtract(mat(T, R), mat(R, T)),
        "BV-VB": rf.subtract(mat(B, V), mat(V, B)),
    }
    for commutator in commutators.values():
        assert any(q for row in commutator for q in row)
    # Certify exactly where replacing VT by TV would change Ax.
    assert At != add(U, mat(T, V))

    assert (count0, counts) == (44, 44)
    assert gain == F(355, 5157)
    assert Up[1][0] == -F(19, 700)
    return dict(raw_states=[count0, counts], physical_gain=str(gain),
                raw_addition_entries=256, raw_binary_entries=128,
                signed_partner_entry=str(Up[1][0]), identities=True)


def numerator_shortcut_audit():
    """J116: independent reachable-block certificate, not a release failure."""
    import verify_release_forest as rf
    K = [[0,F(1,2),F(1,2),0,0], [0,0,0,1,0],
         [F(1,2),0,0,0,F(1,2)], [F(1,2),0,F(1,2),0,0]]
    alg = rf.AdditionAlgebra(K)
    data = rf.ReleaseForest(alg, [0,4,F(1,4),1], 0, degrees=False)
    low, high = [alg.pos[frozenset(C)] for C in ({0,3},{0,1,3})]
    mu = [a-b for a,b in zip(data.P0[high], data.P0[low])]
    reachable = {j for j in data.present if data.P0[low][j] or data.P0[high][j]}
    while True:
        enlarged = reachable | {k for j in reachable for k in data.present
                                if data.W[j][k]}
        if enlarged == reachable:
            break
        reachable = enlarged
    h = [F(2 in B) for B in alg.sets]
    A0 = sum(x*y for x,y in zip(mu,h))
    A1 = A0 + sum(x*y for x,y in zip(rf.multiply([mu],data.W)[0],h))
    trace = sum(data.W[j][j] for j in reachable)
    assert len(reachable) == 7
    assert A0 == F(48869,121282)
    assert A1 == F(9220242269537,17651188228800)
    assert trace == F(14775529,9702560)
    assert A1-trace*A0 == -F(13280929,145538400)
    return dict(reachable_states=7, A0=str(A0), A1=str(A1),
                negative_numerator_coefficient=str(A1-trace*A0))


def main():
    cases = [
        ('two-half', [[F(0), F(1, 2), F(1, 2)],
                      [F(1, 2), F(0), F(1, 2)]], [F(1)]*2),
        ('single-leak-three', [[F(0), F(1), F(0), F(0)],
                               [F(0), F(0), F(1), F(0)],
                               [F(1, 2), F(0), F(0), F(1, 2)]],
         [F(1, 2), F(1), F(2)]),
    ]
    return dict(parameter_flow=[flow_audit(*case) for case in cases],
                logarithm=logarithm_counterexample(),
                release_resolvent=release_resolvent_audit(),
                integrated_release=integrated_release_audit(),
                boundary_interpolation=boundary_interpolation_audit(),
                numerator_shortcut=numerator_shortcut_audit())


if __name__ == '__main__':
    print(json.dumps(main(), indent=2))
