"""Exact source-flow and Hadamard-potential audits for response-renewal §§13–23.

Standard library only. These audit proved identities and scoped obstructions;
the occupation and stopping-flow signs remain open. AI-written research code.
No files are written by this script.
"""

from fractions import Fraction as F
from itertools import combinations, combinations_with_replacement, product
import verify_contact_renewal as cr
import verify_response_renewal as rn


def first_sleep(K, p, disabled=(), taboo=()):
    n = len(p)
    outside = [z for z in range(n) if z not in taboo]
    coins = [F(0) if z in disabled else p[z] for z in range(n)]
    G = rn.inverse([[F(i == j)-(1-coins[i])*K[i][j] for j in outside]
                    for i in outside])
    S = [[F(0)]*n for _ in range(n)]
    for i, t in enumerate(outside):
        for j, u in enumerate(outside):
            S[t][u] = G[i][j]*coins[u]
    return S


def pair_law(S, two_same, supports, a, b):
    n = len(S)
    pos = {C: j for j, C in enumerate(supports)}
    e = [1-sum(row) for row in S]
    M = [S[a][u]*S[b][u]/S[u][u] for u in range(n)]
    row = [sum(q*law[j] for q, law in zip(M, two_same))
           for j in range(len(supports))]
    row[pos[()]] += e[a]*e[b]
    for u in range(n):
        row[pos[(u,)]] += e[a]*S[b][u]+e[b]*S[a][u]-M[u]*e[u]
    for u, z in combinations(range(n), 2):
        row[pos[(u, z)]] += (S[a][u]*S[b][z]+S[a][z]*S[b][u]
                             -M[u]*S[u][z]-M[z]*S[z][u])
    assert all(q >= 0 for q in row) and sum(row) == 1
    return row


def source_flows(K, p, data):
    supports, states, C, T, _, two_same, _, triples, _ = data
    n = len(p)
    idx = {state: i for i, state in enumerate(states)}
    S = first_sleep(K, p)
    h = [[S[t][w]/S[w][w] for w in range(n)] for t in range(n)]
    D, L, _ = cr.creator_flow(K, p, data)
    J = [[q-F(u == w)*h[t][w] for u, q in enumerate(row)]
         for (t, w), row in zip(states, L)]
    IminusD = [[F(u == v)-D[u][v] for v in range(n)] for u in range(n)]
    flows = {}
    pairs = {}
    predicted = {}
    for a, b in combinations_with_replacement(range(n), 2):
        lower = pair_law(S, two_same, supports, a, b)
        pairs[a, b] = lower
        initial = tuple(int(z == a)+int(z == b) for z in range(n))
        predicted[initial] = [sum(q for B, q in zip(supports, lower) if v in B)
                              for v in range(n)]
        for w in range(n):
            if w in (a, b):
                continue
            A = [[S[t][u]-h[t][w]*S[w][u] for u in range(n)]
                 for t in range(n)]
            weights = [F(0)]*len(states)
            for u in range(n):
                if u == w:
                    continue
                collision = A[a][u]*A[b][u]/A[u][u]
                cross = (h[a][w]*A[b][u]+h[b][w]*A[a][u]
                         -collision*h[u][w]
                         +h[a][w]*h[b][w]*S[w][u]/(1-S[w][w]))
                assert collision >= 0 and cross >= 0
                weights[idx[u, w]] += collision
                weights[idx[w, u]] += cross
                traced = (S[a][u]*S[b][u]/S[u][u]
                          +S[a][w]*S[b][w]/S[w][w]*C[w][u]/(1-C[w][w]))
                assert cross*h[w][u] == traced-collision
            assert sum(weights) <= 1
            first = h[a][w]+h[b][w]-h[a][w]*h[b][w]
            flow = [F(u == w)*first+sum(q*row[u] for q, row in zip(weights, J))
                    for u in range(n)]
            assert all(q >= 0 for q in flow)
            flows[a, b, w] = flow
            gain = [F(v == w)-sum(flow[u]*IminusD[u][v] for u in range(n))
                    for v in range(n)]
            initial = [int(z == a)+int(z == b) for z in range(n)]
            initial[w] = -1
            predicted[tuple(initial)] = [q+r for q, r in zip(predicted[tuple(
                int(z == a)+int(z == b) for z in range(n))], gain)]
    laws, _ = rn.stopped(predicted, lambda state: rn.arw_step(K, p, state),
                         lambda state: all(q <= 0 for q in state))
    checked = 0
    for initial, row in predicted.items():
        for v, q in enumerate(row):
            assert q == sum(prob for state, prob in laws[initial].items() if state[v] == -1)
            checked += 1
    source_checks = 0
    for (a, b), lower in pairs.items():
        rho = first_sleep(K, p, disabled={a, b})
        for x in range(n):
            flow = [sum(rho[x][w]*flows[a, b, w][u] for w in range(n) if w not in (a, b))
                    for u in range(n)]
            upper = triples[tuple(sorted((a, b, x)))]
            for v in range(n):
                gain = sum(q-r for B, q, r in zip(supports, upper, lower) if v in B)
                assert gain == rho[x][v]-sum(flow[u]*IminusD[u][v] for u in range(n))
                source_checks += 1
    return S, pairs, dict(raw_source_flow_entries=checked, add_one_source_entries=source_checks)


def hadamard_audit(K, p, data, S, pairs):
    supports, states, _, T, reward, two_same, _, triples, _ = data
    n = len(p)
    idx = {state: i for i, state in enumerate(states)}
    inv_square = rn.inverse([[q*q for q in row] for row in S])
    assert all(inv_square[i][j] <= 0 for i in range(n) for j in range(n) if i != j)
    four_point_checks = 0
    residual_bad = []
    killed_bad = []
    for a, b in combinations(range(n), 2):
        lam = [sum(S[a][j]*S[b][j]*inv_square[j][t] for j in range(n))
               for t in range(n)]
        assert all(q >= 0 for q in lam) and sum(lam) <= 1
        assert lam[a] >= S[b][a]/(2*S[a][a])
        assert lam[b] >= S[a][b]/(2*S[b][b])
        for u, z in combinations(range(n), 2):
            defect = (S[a][u]*S[b][z]+S[a][z]*S[b][u]
                      -2*sum(lam[t]*S[t][u]*S[t][z] for t in range(n)))
            assert defect >= 0
            four_point_checks += 1
        residual = [q-sum(lam[t]*two_same[t][j] for t in range(n))
                    for j, q in enumerate(pairs[a, b])]
        assert all(q >= 0 for q in residual)
        assert sum(residual) == 1-sum(lam)
        for w in range(n):
            if w in (a, b):
                continue
            ts = [t for t in range(n) if t != w]
            A = [[S[t][u]-S[t][w]*S[w][u]/S[w][w] for u in ts] for t in ts]
            inv = rn.inverse([[q*q for q in row] for row in A])
            killed = [sum(A[ts.index(a)][j]*A[ts.index(b)][j]*inv[j][t]
                          for j in range(len(ts))) for t in range(len(ts))]
            for t, q in zip(ts, killed):
                if q > lam[t]:
                    killed_bad.append((a, b, w, t, str(q-lam[t])))
        for x in range(n):
            for v in range(n):
                gain = sum(q-r for B, q, r in zip(
                    supports, triples[tuple(sorted((a, b, x)))], pairs[a, b]) if v in B)
                diagonal = [sum(q-r for B, q, r in zip(
                    supports, triples[tuple(sorted((t, t, x)))], two_same[t]) if v in B)
                            for t in range(n)]
                residual_gain = gain-sum(q*r for q, r in zip(lam, diagonal))
                if residual_gain < 0:
                    residual_bad.append((a, b, x, v, str(residual_gain)))
    rate_checks = 0
    for w in range(n):
        ts = [t for t in range(n) if t != w]
        A = [[S[t][u]-S[t][w]*S[w][u]/S[w][w] for u in ts] for t in ts]
        invA2 = rn.inverse([[q*q for q in row] for row in A])
        for j, u in enumerate(ts):
            for k, z in enumerate(ts):
                coeff = [sum(q*A[r][j]*A[r][k] for r, q in enumerate(row))
                         for row in invA2]
                assert all(q >= 0 for q in coeff)
                if j != k:
                    assert coeff[j] >= A[j][k]/(2*A[j][j])
        CW = [[T[idx[t, w]][idx[u, w]] for u in ts] for t in ts]
        cross = [[T[idx[t, w]][idx[w, u]] for u in ts] for t in ts]
        inv = rn.inverse(CW)
        pre_cross = rn.multiply(inv, cross)
        pre_reward = rn.multiply(inv, [reward[idx[t, w]] for t in ts])
        b = [S[t][w]/S[w][w] for t in ts]
        pre_rate = [sum(q*z*(1-z) for q, z in zip(row, b)) for row in inv]
        assert all(q >= 0 for row in pre_cross+pre_reward for q in row)
        assert all(q >= 0 for q in pre_rate)
        assert all(sum(row) >= 0 for row in pre_reward)
        rate_checks += len(ts)*(len(ts)+len(supports)+1)
    return dict(four_point_entries=four_point_checks, preconditioned_entries=rate_checks,
                negative_residual_response=residual_bad[:1],
                killed_weight_increase=killed_bad[:1])


def power_algebra_audit(S):
    """Degree-three audits of the general Holder/product identities."""
    from collections import defaultdict
    from itertools import combinations_with_replacement
    n = len(S)

    def mul(P, Q):
        result = defaultdict(F)
        for a, x in P.items():
            for b, y in Q.items():
                result[tuple(sorted(a+b))] += x*y
        return dict(result)

    def combine(P, coefficient, Q):
        result = defaultdict(F, P)
        for a, x in Q.items():
            result[a] += coefficient*x
        return {a: x for a, x in result.items() if x}

    linear = [{(): 1-sum(row), **{(u,): q for u, q in enumerate(row)}}
              for row in S]
    weights, remainder, independent = {}, {}, {}
    checks = 0
    for m in (2, 3):
        inv = rn.inverse([[q**m for q in row] for row in S])
        all_same = {}
        for t in range(n):
            P = {(): F(1)}
            for _ in range(m):
                P = mul(P, linear[t])
            all_same[t] = P
        for roots in combinations_with_replacement(range(n), m):
            values = [F(1)]*n
            P = {(): F(1)}
            for a in roots:
                values = [q*S[a][j] for j, q in enumerate(values)]
                P = mul(P, linear[a])
            lam = [sum(values[j]*inv[j][t] for j in range(n)) for t in range(n)]
            assert all(q >= 0 for q in lam) and sum(lam) <= 1
            R = P
            for t, q in enumerate(lam):
                R = combine(R, -q, all_same[t])
            assert all(q >= 0 for q in R.values())
            assert sum(R.values()) == 1-sum(lam)
            assert all(R.get((u,)*m, F(0)) == 0 for u in range(n))
            weights[roots], remainder[roots], independent[roots] = lam, R, P
            checks += len(R)+n
    associative = 0
    for a in range(n):
        for b in range(n):
            ab = tuple(sorted((a, b)))
            for c in range(n):
                abc = tuple(sorted((a, b, c)))
                composed = [sum(weights[ab][t]*weights[tuple(sorted((t, t, c)))][u]
                                for t in range(n)) for u in range(n)]
                assert composed == weights[abc]
                R = mul(remainder[ab], linear[c])
                for t, q in enumerate(weights[ab]):
                    R = combine(R, q, remainder[tuple(sorted((t, t, c)))])
                assert {a: x for a, x in R.items() if x} == remainder[abc]
                associative += n
    return dict(positive_polynomial_entries=checks, associative_entries=associative)


def reference_creator_audit(K, p, data):
    supports, _, _, _, _, two_same, _, triples, _ = data
    n = len(p)
    reference = first_sleep(K, [max(q, F(1, 2)) for q in p])
    actual = first_sleep(K, p)
    checked = 0
    for m in (2, 3):
        inv = rn.inverse([[q**m for q in row] for row in reference])
        diagonal = two_same if m == 2 else [triples[(t, t, t)] for t in range(n)]
        for roots in combinations_with_replacement(range(n), m):
            products = [F(1)]*n
            for a in roots:
                products = [q*reference[a][u] for u, q in enumerate(products)]
            lam = [sum(products[j]*inv[j][t] for j in range(n)) for t in range(n)]
            assert all(q >= 0 for q in lam) and sum(lam) <= 1
            law = (pair_law(actual, two_same, supports, *roots)
                   if m == 2 else triples[roots])
            residual = [q-sum(lam[t]*diagonal[t][j] for t in range(n))
                        for j, q in enumerate(law)]
            assert all(q >= 0 for q in residual)
            assert sum(residual) == 1-sum(lam)
            checked += len(residual)
    return dict(reference_creator_law_entries=checked)


def toppling_laws(Q, alpha, bound=3):
    """Capacity-one sandpile; unstable pairs emit one or two Q-walkers."""
    from collections import defaultdict
    n = len(Q)
    initials = set()
    for roots in combinations_with_replacement(range(n+1), bound):
        initials.add(tuple(roots.count(z) for z in range(n)))

    def step(state):
        z = next(z for z, count in enumerate(state) if count >= 2)
        out = defaultdict(F)
        if alpha[z]:
            for y, prob in enumerate(Q[z]):
                if prob:
                    nxt = list(state)
                    nxt[z] -= 1
                    if y < n:
                        nxt[y] += 1
                    out[tuple(nxt)] += alpha[z]*prob
        if alpha[z] < 1:
            for y, prob in enumerate(Q[z]):
                for w, prob2 in enumerate(Q[z]):
                    if prob*prob2:
                        nxt = list(state)
                        nxt[z] -= 2
                        if y < n:
                            nxt[y] += 1
                        if w < n:
                            nxt[w] += 1
                        out[tuple(nxt)] += (1-alpha[z])*prob*prob2
        return dict(out)

    laws, _ = rn.stopped(initials, step, lambda state: max(state) <= 1)
    return laws


def deposited_law(S, roots, laws, supports):
    from collections import defaultdict
    n = len(S)
    mix = {(0,)*n: F(1)}
    for a in roots:
        nxt = defaultdict(F)
        for state, mass in mix.items():
            for y, prob in enumerate(S[a]+[1-sum(S[a])]):
                if prob:
                    state2 = list(state)
                    if y < n:
                        state2[y] += 1
                    nxt[tuple(state2)] += mass*prob
        mix = dict(nxt)
    final = defaultdict(F)
    for state, mass in mix.items():
        for end, prob in laws[state].items():
            final[tuple(z for z, count in enumerate(end) if count)] += mass*prob
    return [final[C] for C in supports]


def toppling_audit(K, p):
    reference = [max(q, F(1, 2)) for q in p]
    n = len(p)
    S = first_sleep(K, reference)
    e = [1-sum(row) for row in S]
    Q = [[sum(K[z][y]*S[y][u] for y in range(n)) for u in range(n)]
         + [K[z][n]+sum(K[z][y]*e[y] for y in range(n))]
         for z in range(n)]
    alpha = [(2*q-1)/q for q in reference]
    assert all(sum(row) == 1 and min(row) >= 0 for row in Q)
    assert all(0 <= q <= 1 for q in alpha)
    for z in range(n):
        for u in range(n):
            assert S[z][u] == reference[z]*F(z == u)+(1-reference[z])*Q[z][u]
    laws = toppling_laws(Q, alpha)
    data = cr.contact_model(K, reference)
    supports, _, _, _, _, two_same, _, triples, _ = data
    checked = 0
    for m in (2, 3):
        for roots in combinations_with_replacement(range(n), m):
            expected = (pair_law(S, two_same, supports, *roots) if m == 2
                        else triples[roots])
            assert deposited_law(S, roots, laws, supports) == expected
            checked += len(supports)
    return dict(toppling_law_entries=checked)


def toppling_and_partition_obstructions():
    Q = [[F(0), F(1, 3), F(1, 2), F(1, 6)],
         [F(0), F(1, 3), F(1, 3), F(1, 3)],
         [F(1, 2), F(0), F(0), F(1, 2)]]
    n = len(Q)
    S = [[(F(i == j)+Q[i][j])/2 for j in range(n)] for i in range(n)]
    supports = [C for k in range(n+1) for C in combinations(range(n), k)]
    laws = toppling_laws(Q, [F(0)]*n)
    lower = deposited_law(S, (0, 2), laws, supports)
    upper = deposited_law(S, (0, 1, 2), laws, supports)
    gain = sum(hi-lo for C, lo, hi in zip(supports, lower, upper) if 2 in C)
    assert gain == -F(2273, 434166)
    assert rn.inverse(S)[1][0] == F(6, 19)

    # Deleting labels cannot give all these pair/full-merge marginals.
    S = [[F(8, 15), F(2, 15)], [F(2, 15), F(8, 15)]]
    masses = {}
    for roots in ((0, 0), (0, 1), (0, 0, 1)):
        m = len(roots)
        inv = rn.inverse([[q**m for q in row] for row in S])
        products = [F(1)]*2
        for a in roots:
            products = [q*S[a][u] for u, q in enumerate(products)]
        masses[roots] = sum(products[j]*inv[j][t] for j in range(2) for t in range(2))
    assert masses == {(0, 0): F(1), (0, 1): F(8, 17), (0, 0, 1): F(4, 13)}
    assert masses[0, 0, 1] < masses[0, 1]
    return dict(general_sandpile_gain=str(gain),
                pair_merge=str(masses[0, 1]), triple_merge=str(masses[0, 0, 1]))


def determinant(matrix):
    """Exact determinant, including the empty principal minor."""
    A = [list(row) for row in matrix]
    value = F(1)
    for i in range(len(A)):
        pivot = next((j for j in range(i, len(A)) if A[j][i]), None)
        if pivot is None:
            return F(0)
        if pivot != i:
            A[i], A[pivot] = A[pivot], A[i]
            value = -value
        entry = A[i][i]
        value *= entry
        for j in range(i+1, len(A)):
            ratio = A[j][i]/entry
            for k in range(i+1, len(A)):
                A[j][k] -= ratio*A[i][k]
    return value


def subsets(vertices):
    vertices = tuple(vertices)
    return [C for k in range(len(vertices)+1) for C in combinations(vertices, k)]


def internal_block_expansion(K, p, roots):
    """Expand a complete sleep block after cancelling internal root traffic."""
    from collections import defaultdict
    n = len(p)
    R = tuple(sorted(roots))
    outside = [y for y in range(n+1) if y not in R]
    odds = [(1-q)/q for q in p]

    def minor(C):
        return determinant([[F(i == j)-K[i][j] for j in C] for i in C])

    result = defaultdict(F)
    for B in subsets(R):
        released = [z for z in R if z not in B]
        for C in subsets(B):
            G = rn.inverse([[F(i == j)-K[i][j] for j in C] for i in C]) if C else []
            factor = F((-1)**len(released))*minor(B)/minor(C)
            for z in R:
                if z not in C:
                    factor *= odds[z]
            poly = {B: factor}
            for z in released:
                nxt = defaultdict(F)
                for y in outside:
                    weight = K[z][y]+sum(K[z][a]*G[i][j]*K[b][y]
                                         for i, a in enumerate(C) for j, b in enumerate(C))
                    assert weight >= 0
                    for active, prob in poly.items():
                        new = active if y == n else tuple(sorted(active+(y,)))
                        nxt[new] += prob*weight
                poly = dict(nxt)
            for active, prob in poly.items():
                result[active] += prob
    return {active: prob for active, prob in result.items() if prob}


def internal_block_audit(K, p, data, S):
    supports, _, _, _, _, two_same, _, triples, _ = data
    n = len(p)
    cache = {(): [F(C == ()) for C in supports]}
    for z in range(n):
        cache[z,] = [1-sum(S[z]) if not C else S[z][C[0]] if len(C) == 1
                     else F(0) for C in supports]
    for a, b in combinations_with_replacement(range(n), 2):
        cache[a, b] = pair_law(S, two_same, supports, a, b)
    cache.update(triples)
    checks = 0
    for R in subsets(range(n)):
        if not 1 <= len(R) <= 3:
            continue
        coefficients = internal_block_expansion(K, p, R)
        law = [sum(weight*cache[active][j] for active, weight in coefficients.items())
               for j in range(len(supports))]
        assert law == [F(C == R) for C in supports]
        for active, weight in coefficients.items():
            retained = set(active) & set(R)
            assert weight*(-1)**(len(R)-len(retained)) > 0
        checks += len(supports)
    return dict(internal_block_entries=checks)


def full_inverse_audit():
    checks = 0
    for case in rn.CASES:
        K = [[F(q) for q in row] for row in case["K"]]
        p = list(map(F, case["p"]))
        n = len(p)
        expansion = internal_block_expansion(K, p, range(n))
        row = [expansion.get(tuple(z for z, bit in enumerate(C) if bit), F(0))
               for C in case["subsets"]]
        assert row == list(map(F, case["inverse_matrix"][-1]))
        assert F(case["stabilization_matrix"][-1][-1])*row[-1] == 1
        checks += len(row)
    return dict(full_inverse_entries=checks)

def inverse_diagonal_audit():
    checks = 0
    for case in rn.CASES:
        K = [[F(q) for q in row] for row in case["K"]]
        p = list(map(F, case["p"]))
        odds = [(1-q)/q for q in p]
        for i, bits in enumerate(case["subsets"]):
            R = tuple(z for z, bit in enumerate(bits) if bit)
            dR = determinant([[F(a == b)-K[a][b] for b in R] for a in R])
            total = F(0)
            for C in subsets(R):
                weight = F(1)
                for z in R:
                    if z not in C:
                        weight *= odds[z]
                dC = determinant([[F(a == b)-K[a][b] for b in C] for a in C])
                total += weight/dC
            assert dR*total == F(case["inverse_matrix"][i][i])
            checks += 1
    return dict(inverse_diagonal_entries=checks)


def causal_charge_obstruction():
    K = [[F(0), F(1), F(0), F(0)],
         [F(0), F(0), F(1), F(0)],
         [F(0), F(0), F(0), F(1)]]
    checked = 0
    for p in ([F(1, 2)]*3, [F(2, 5), F(3, 7), F(5, 11)]):
        q = [1-z for z in p]
        S = first_sleep(K, p)
        h = [[S[t][w]/S[w][w] for w in range(3)] for t in range(3)]
        data = cr.contact_model(K, p)
        _, states, _, T, _, _, _, _, _ = data
        D, _, _ = cr.creator_flow(K, p, data)
        assert D[2] == [F(0)]*3
        c = states.index((1, 2))
        taboo = first_sleep(K, p, taboo={2})
        contact_entry = taboo[0][1]*taboo[1][1]/taboo[1][1]
        reward = h[1][2]*(1-h[1][2])
        f = h[0][2]+h[1][2]-h[0][2]*h[1][2]
        assert contact_entry == q[0]*p[1] and T[c][c] == p[1]
        assert reward == q[1]*p[1] and f == q[1]*(1+q[0]*p[1])
        direct = p[2]*f
        contact = p[2]*contact_entry/(1-T[c][c])*reward
        shortage = direct*contact_entry+contact-p[2]*contact_entry
        assert contact == p[2]*q[0]*p[1]**2
        assert shortage == p[2]*q[0]**2*q[1]*p[1]**2 > 0
        assert direct+contact == p[2]*(1-p[0]*p[1]) < p[2]
        assert contact <= p[2]*contact_entry  # Complete-path singleton cut.
        initials = [(1, 1, 0), (1, 1, -1), (1, 1, 1)]
        laws, _ = rn.stopped(initials, lambda state: rn.arw_step(K, p, state),
                             lambda state: all(z <= 0 for z in state))
        marginal = {state: [sum(prob for end, prob in law.items() if end[v] == -1)
                            for v in range(3)] for state, law in laws.items()}
        assert [b-a for a, b in zip(marginal[initials[0]], marginal[initials[1]])] == [F(0), F(0), p[0]*p[1]]
        assert [b-a for a, b in zip(marginal[initials[0]], marginal[initials[2]])] == [F(0), F(0), p[0]*p[1]*p[2]]
        if p == [F(1, 2)]*3:
            assert direct*contact_entry+contact == F(9, 64)
            assert p[2]*contact_entry == F(8, 64)
        checked += 1
    return dict(causal_charge_kernels=checked, half_sleep_cut="9/64 > 8/64")

def frozen_step(K,p,arrows,state):
    counts, frozen = state
    n=len(p)
    z=next((i for i,q in enumerate(counts) if q>0),None)
    if z is None:return {}
    out={}
    if counts[z]==1:
        new=list(counts);new[z]=-1
        out[tuple(new),frozen]=p[z]
    jump=1 if counts[z]>1 else 1-p[z]
    destinations=[(arrows[z],F(1))] if frozen>>z&1 else [(y,q) for y,q in enumerate(K[z]) if q]
    for y,weight in destinations:
        new=list(counts);new[z]-=1
        if y<n:new[y]=2 if new[y]==-1 else new[y]+1
        flag=frozen&~(1<<z)
        key=tuple(new),flag
        out[key]=out.get(key,F(0))+jump*weight
    return out

def pair_forest_law(K,p,R,supports):
    n=len(p);a,b=R
    assert K[a][a]==K[b][b]==0
    odds=[(1-q)/q for q in p]
    cycle=K[a][b]*K[b][a]
    Z=1+(1-cycle)*(odds[a]+odds[b]+odds[a]*odds[b])
    result={B:F(B==R)/Z for B in supports}
    total=F(1)/Z
    for ya,yb in product(range(n+1),repeat=2):
        mass=K[a][ya]*K[b][yb]
        if not mass or (ya==b and yb==a):continue
        arrows={a:ya,b:yb}
        initials={}
        for mask in (1,2,3):
            counts=[0]*n;frozen=0;weight=mass/Z
            for j,z in enumerate(R):
                if mask>>j&1:weight*=odds[z]
                else:counts[z]=-1;frozen|=1<<z
            for j,z in enumerate(R):
                if mask>>j&1 and arrows[z]<n:
                    y=arrows[z]
                    counts[y]=2 if counts[y]==-1 else counts[y]+1
            state=tuple(counts),frozen
            initials[state]=initials.get(state,F(0))+weight
        laws,_=rn.stopped(initials,lambda state:frozen_step(K,p,arrows,state),
                          lambda state:all(q<=0 for q in state[0]))
        for initial,weight in initials.items():
            total+=weight
            for (counts,_),prob in laws[initial].items():
                B=tuple(z for z,q in enumerate(counts) if q==-1)
                result[B]+=weight*prob
    assert total==1
    return [result[B] for B in supports]


def pair_forest_audit(K, p, data, S):
    supports, _, _, _, _, two_same, _, _, _ = data
    checks = 0
    for R in combinations(range(len(p)), 2):
        if any(K[z][z] for z in R):
            continue
        expected = pair_law(S, two_same, supports, *R)
        assert pair_forest_law(K, p, R, supports) == expected
        checks += len(supports)
    return dict(pair_forest_law_entries=checks)



def inverse_diagonal_set_audit(K, p):
    """All higher mixed set differences of log Z have nonpositive sign."""
    n = len(p)
    odds = [(1-q)/q for q in p]
    sets = [frozenset(C) for C in subsets(range(n))]
    minors = {C: determinant([[F(i == j)-K[i][j] for j in sorted(C)]
                             for i in sorted(C)]) for C in sets}
    Z = {}
    for R in sets:
        value = F(0)
        for C in sets:
            if C <= R:
                weight = F(1)
                for i in R-C:
                    weight *= odds[i]
                value += weight/minors[C]
        Z[R] = minors[R]*value
    checked = 0
    for A in sets:
        for B in sets:
            if A & B or len(B) < 2:
                continue
            positive = negative = F(1)
            for C in sets:
                if C <= B:
                    if (len(B)-len(C)) % 2:
                        negative *= Z[A | C]
                    else:
                        positive *= Z[A | C]
            assert positive <= negative, (A, B, positive, negative)
            checked += 1
    return dict(log_Z_mixed_differences=checked)



def loss_marker_audit(K, p):
    """Conditional marker correlations, using exact intersection probabilities."""
    from math import prod
    n = len(p)
    V = frozenset(range(n))
    sets = [frozenset(C) for C in subsets(V)]
    odds = [(1-q)/q for q in p]
    d = {C: determinant([[F(i == j)-K[i][j] for j in sorted(C)]
                         for i in sorted(C)]) for C in sets}
    Z = {R: d[R]*sum(prod(odds[i] for i in R-C)/d[C]
                     for C in sets if C <= R) for R in sets}
    a = [odds[i]*K[i][n] for i in range(n)]
    moments = {S: prod(a[i] for i in S)*Z[V-S]/Z[V] for S in sets}
    masses = {S: sum((-1)**len(D-S)*moments[D] for D in sets if S <= D)
              for S in sets}
    assert min(masses.values()) >= 0 and sum(masses.values()) == 1
    checks = 0
    for assignment in product((-1, 0, 1), repeat=n):
        yes = {i for i, state in enumerate(assignment) if state == 1}
        no = {i for i, state in enumerate(assignment) if state == 0}
        free = [i for i, state in enumerate(assignment) if state == -1]
        rows = [(S, q) for S, q in masses.items() if yes <= S and not no & S]
        total = sum(q for _, q in rows)
        if not total:
            continue
        for i, j in combinations(free, 2):
            pi = sum(q for S, q in rows if i in S)
            pj = sum(q for S, q in rows if j in S)
            pij = sum(q for S, q in rows if i in S and j in S)
            assert total*pij <= pi*pj
            checks += 1
    return dict(conditional_marker_pairs=checks)


def marker_occupation_obstruction():
    K = [[F(0), F(1, 2), F(1, 2)], [F(1, 2), F(0), F(1, 2)]]
    p = [F(1, 2)]*2
    S = first_sleep(K, p)
    data = cr.contact_model(K, p)
    supports, _, _, _, _, two_same, _, _, _ = data
    law = pair_law(S, two_same, supports, 0, 1)
    q1 = sum(q for C, q in zip(supports, law) if 1 in C)
    assert q1 == F(22, 39)
    marker_prob = F(4, 13)
    covariance = marker_prob*(S[1][1]-q1)
    assert covariance == -F(8, 845)
    # Given C={1}, two of the three equiprobable forests have F_0=s.
    # On that event root 1 is left sleeping; otherwise 0->1 wakes it,
    # its frozen forest arrow is 1->s, and one fresh active remains at 1.
    conditional_covariance = F(2, 3)*F(1, 3)*(1-S[1][1])
    assert conditional_covariance == F(14, 135)
    assert S[1][0] == F(2, 15)
    return dict(marker_occupation_cov=str(covariance),
                given_latent_subset_cov=str(conditional_covariance))

def terminal_pair_flow_audit():
    """Direct absorption audit of the terminal-pair reset and stopping rule."""
    K = [[F(1, 6), F(1, 3), F(1, 3), F(0), F(1, 6)],
         [F(1, 4), F(1, 4), F(0), F(1, 4), F(1, 4)],
         [F(0), F(0), F(1, 4), F(1, 2), F(1, 4)],
         [F(0), F(0), F(1, 3), F(1, 3), F(1, 3)]]
    p = [F(2, 5), F(3, 7), F(4, 9), F(5, 11)]
    n, a, w = 4, 2, 3
    data = cr.contact_model(K, p)
    D, canonical_L, _ = cr.creator_flow(K, p, data)
    G = rn.inverse([[F(i == j)-D[i][j] for j in range(n)] for i in range(n)])
    states = data[1]
    common = [(2, 0, 0, 0), (1, 1, 0, 0), (3, 0, 0, 0)]
    initial = common+[state[:w]+(-1,) for state in common]
    laws, _ = rn.stopped(initial, lambda state: rn.arw_step(K, p, state),
                         lambda state: all(q <= 0 for q in state))
    def visit_step(state):
        counts, visited = state
        return {(end, visited or end[w] != 0): prob
                for end, prob in rn.arw_step(K, p, counts).items()}
    visit_laws, _ = rn.stopped([(state, False) for state in common], visit_step,
                               lambda state: all(q <= 0 for q in state[0]))
    def partial_step(state):
        if any(state[i] > 0 for i in (0, 1)):
            return rn.arw_step(K, p, state)
        z = next((i for i in (a, w) if state[i] > 1), None)
        if z is None:
            return {}
        out = {}
        for y, prob in enumerate(K[z]):
            if not prob:
                continue
            end = list(state)
            end[z] -= 1
            if y < n:
                end[y] += 1
            key = tuple(end)
            out[key] = out.get(key, F(0))+prob
        return out
    idla_laws, _ = rn.stopped(common, partial_step,
        lambda state: all(state[i] <= 0 for i in (0, 1))
                      and all(state[i] <= 1 for i in (a, w)))
    checked = 0
    for state in common:
        f = sum(q for (_, visited), q in visit_laws[state, False].items() if visited)
        zeta = sum(q for end, q in idla_laws[state].items()
                   if end[a] == end[w] == 1)
        assert 0 <= zeta <= f <= 1
        expected = [(1-f)*F(v == w)+(f-zeta)*D[w][v] for v in range(n)]
        upper = state[:w]+(-1,)
        actual = [sum(q for end, q in laws[upper].items() if end[v] == -1)
                  -sum(q for end, q in laws[state].items() if end[v] == -1)
                  for v in range(n)]
        assert actual == expected and all(q >= 0 for q in expected)
        if state[0] == 2 and sum(state) == 2:
            row = canonical_L[states.index((0, w))]
            assert row == [zeta*G[w][v]+(f-zeta)*F(v == w) for v in range(n)]
        checked += 1
    return dict(terminal_pair_inputs=checked, arbitrary_common_multiplicity=3)

def loop_thinning_coefficients(bound, p):
    """Formal coefficients with the common K-edge monomial removed."""
    from math import factorial, prod
    n = len(p)
    tables = []
    for counts in product(*[range(k+1) for k in bound]):
        outgoing = [sum(counts[i*n+j] for j in range(n)) for i in range(n)]
        incoming = [sum(counts[i*n+j] for i in range(n)) for j in range(n)]
        if outgoing != incoming:
            continue
        number = F(prod(factorial(k) for k in outgoing),
                   prod(factorial(k) for k in counts))
        sleep_weight = prod(p[i] for i, k in enumerate(outgoing) if k)
        tables.append((sum(counts), counts, number*sleep_weight, number))
    tables.sort()
    g, f, c, ell = {}, {}, {}, {}
    for degree, counts, g_value, f_value in tables:
        g[counts], f[counts] = g_value, f_value
        if not degree:
            continue
        c_value, ell_value = g_value, f_value
        for part, coefficient in c.items():
            if all(x <= y for x, y in zip(part, counts)):
                rest = tuple(y-x for x, y in zip(part, counts))
                c_value -= F(sum(part), degree)*coefficient*g.get(rest, 0)
                ell_value -= F(sum(part), degree)*ell[part]*f.get(rest, 0)
        c[counts], ell[counts] = c_value, ell_value
    return c, ell


def closed_words_with_counts(counts, n):
    """Rooted closed vertex words, with each vertex word counted once."""
    result = []
    def extend(start, path, remaining):
        z = path[-1]
        if not any(remaining):
            if z == start:
                result.append(path)
            return
        for y in range(n):
            edge = z*n+y
            if remaining[edge]:
                rest = list(remaining)
                rest[edge] -= 1
                extend(start, path+[y], tuple(rest))
    for start in range(n):
        extend(start, [start], counts)
    return result


def loop_thinning_coefficient_audit():
    examples = [
        (2, (0, 2, 2, 0), F(7, 16)),
        (3, (0, 1, 1, 1, 0, 0, 1, 0, 0), F(3, 16)),
        (3, (0, 1, 1, 1, 0, 1, 1, 1, 0), F(53, 192)),
        (3, (1, 1, 0, 0, 1, 1, 1, 0, 0), F(9, 32)),
        (3, (0, 2, 0, 1, 0, 1, 1, 0, 0), F(7, 32)),
    ]
    checked = 0
    for n, counts, ratio in examples:
        c, ell = loop_thinning_coefficients(counts, [F(1, 2)]*n)
        larger, _ = loop_thinning_coefficients(counts, [F(2, 3)]*n)
        for part in c:
            assert 0 <= c[part] <= larger[part] <= ell[part]
            checked += 1
        assert c[counts]/ell[counts] == ratio
        # Independently count the trace-log coefficient, including repeated tours.
        assert ell[counts] == F(len(closed_words_with_counts(counts, n)), sum(counts))
    return dict(thinned_loop_coefficients=checked, closed_word_classes=len(examples))


def periodic_full_stabilizes(arrows, bits):
    """Closed periodic buffer book: exact termination by repeated-state detection."""
    n = len(arrows)
    phase = [0]*n
    residents = [-1]*n
    pending = [int(bool(row)) for row in arrows]
    seen = set()
    while any(pending):
        state = tuple(phase), tuple(residents), tuple(pending)
        if state in seen:
            return False
        seen.add(state)
        z = next(i for i, q in enumerate(pending) if q)
        pending[z] -= 1
        k = phase[z]
        phase[z] = (k+1) % len(arrows[z])
        if residents[z] >= 0:
            pending[residents[z]] += 1
            residents[z] = -1
        if bits[z][k]:
            residents[z] = arrows[z][k]
        else:
            pending[arrows[z][k]] += 1
    return True


def periodic_loop_audit():
    examples = [
        (2, (0, 2, 2, 0), [F(1, 2)]*2),
        (3, (0, 1, 1, 1, 0, 0, 1, 0, 0), [F(1, 2)]*3),
        (3, (0, 1, 1, 1, 0, 1, 1, 1, 0), [F(1, 2)]*3),
        (3, (1, 1, 0, 0, 1, 1, 1, 0, 0), [F(1, 2)]*3),
        (3, (0, 2, 0, 1, 0, 1, 1, 0, 0), [F(1, 2)]*3),
        (3, (0, 1, 1, 1, 0, 1, 1, 1, 0), [F(1, 3), F(2, 5), F(3, 7)]),
        (3, (0, 2, 2, 2, 0, 0, 2, 0, 0), [F(1, 3), F(2, 5), F(3, 7)]),
    ]
    checks = 0
    probabilities = []
    for n, counts, p in examples:
        c, ell = loop_thinning_coefficients(counts, p)
        words = closed_words_with_counts(counts, n)
        total = F(0)
        for word in words:
            arrows = [[] for _ in range(n)]
            for z, y in zip(word, word[1:]):
                arrows[z].append(y)
            for raw in product((0, 1), repeat=sum(counts)):
                bits, offset, weight = [], 0, F(1)
                for z, row in enumerate(arrows):
                    local = raw[offset:offset+len(row)]
                    offset += len(row)
                    bits.append(local)
                    for bit in local:
                        weight *= p[z] if bit else 1-p[z]
                if periodic_full_stabilizes(arrows, bits):
                    total += weight
                checks += 1
        probability = total/len(words)
        assert probability == c[counts]/ell[counts]
        probabilities.append(str(probability))
    assert probabilities[-2:] == ['1504/11025', '6152/33075']
    return dict(periodic_decorated_words=checks, periodic_loop_classes=len(examples),
                heterogeneous_acceptance=probabilities[-2:])

def common_retained_loop_obstruction():
    """The common retained-loop sector cannot be fixed during comparison."""
    from collections import defaultdict
    K = [[F(0), F(1, 2), F(0), F(1, 2)],
         [F(1, 2), F(0), F(0), F(1, 2)],
         [F(1), F(0), F(0), F(0)]]
    p = [F(1, 2)]*3
    values = []
    for m in (1, 2):
        upper = F(0)
        for forest in ((1, 3, 0), (3, 0, 0), (3, 3, 0)):
            words = [[1]*m+[forest[0]], [0]*m+[forest[1]], [0]]
            def step(state):
                counts, indices = state
                z = next(i for i, q in enumerate(counts) if q > 0)
                coin = F(1) if z in (0, 1) and indices[z] == 0 else p[z]
                out = defaultdict(F)
                if counts[z] == 1:
                    end = list(counts)
                    end[z] = -1
                    out[tuple(end), indices] += coin
                jump = F(1) if counts[z] > 1 else 1-coin
                destinations = ([(words[z][indices[z]], F(1))]
                                if indices[z] < len(words[z]) else enumerate(K[z]))
                for y, weight in destinations:
                    if not weight or not jump:
                        continue
                    end = list(counts)
                    end[z] -= 1
                    if y < 3:
                        end[y] = 2 if end[y] == -1 else end[y]+1
                    changed = list(indices)
                    changed[z] = min(indices[z]+1, len(words[z]))
                    out[tuple(end), tuple(changed)] += jump*weight
                return out
            initial = ((1, 1, 1), (0, 0, 0))
            laws, _ = rn.stopped([initial], step,
                                 lambda state: all(q <= 0 for q in state[0]))
            upper += sum(q for (state, _), q in laws[initial].items()
                         if state[0] == -1)/3
        assert upper == F(61, 78)
        values.append(str(upper-1))
    return dict(common_loop_conditional_gains=values, actual_gain="0")

def promotion_latent_audit(K, p):
    from math import prod
    n = len(p)
    sets = [frozenset(C) for C in subsets(range(n))]
    odds = [(1-q)/q for q in p]
    d = {C: determinant([[F(i == j)-K[i][j] for j in sorted(C)]
                         for i in sorted(C)]) for C in sets}
    H = {R: sum(prod(odds[i] for i in R-C)/d[C] for C in sets if C <= R)
         for R in sets}
    Z = {R: d[R]*H[R] for R in sets}
    checked = conditional = 0
    for T in sets:
        for x in T:
            R = T-{x}
            sigma = d[T]/d[R]
            accepted = {}
            for C in sets:
                if not C <= R:
                    continue
                lower = prod(odds[i] for i in R-C)/(d[C]*H[R])
                accepted[C] = lower*(1-p[x])*sigma
                accepted[C | {x}] = lower*p[x]*sigma*d[C]/d[C | {x}]
            total = sum(accepted.values())
            assert total == p[x]*Z[T]/Z[R] and sigma <= total <= 1
            for C, mass in accepted.items():
                expected = prod(odds[i] for i in T-C)/(d[C]*H[T])
                assert mass/total == expected
                checked += 1
            for S in sets:
                if not S <= R:
                    continue
                old = {C: prod(odds[i] for i in R-C)/d[C]
                       for C in sets if S <= C <= R}
                modified = {C: prod((F(0) if i in S else odds[i]) for i in R-C)/d[C]
                            for C in sets if C <= R}
                normalizer = sum(old.values())
                assert sum(modified.values()) == normalizer
                assert all(modified[C] == old.get(C, 0) for C in modified)
                conditional += len(modified)
    return dict(promotion_subset_weights=checked, promotion_refresh_weights=conditional)

def audit(K, p):
    data = cr.contact_model(K, p)
    S, pairs, result = source_flows(K, p, data)
    result.update(hadamard_audit(K, p, data, S, pairs))
    result.update(power_algebra_audit(S))
    result.update(reference_creator_audit(K, p, data))
    result.update(toppling_audit(K, p))
    result.update(internal_block_audit(K, p, data, S))
    result.update(pair_forest_audit(K, p, data, S))
    result.update(inverse_diagonal_set_audit(K, p))
    result.update(loss_marker_audit(K, p))
    result.update(promotion_latent_audit(K, p))
    return result



def creator_exit_forest_obstruction():
    K = [[F(0), F(1, 2), F(1, 2)], [F(1, 2), F(0), F(1, 2)]]
    p = [F(1, 2)]*2
    data = cr.contact_model(K, p)
    supports, _, _, _, _, two_same, _, _, _ = data
    S = first_sleep(K, p)
    law = pair_law(S, two_same, supports, 0, 1)
    stationary = sum(prob for C, prob in zip(supports, law) if 0 in C)
    assert stationary == F(22, 39) and S[1][0] == F(2, 15)
    lower = F(3, 13)*(S[1][0]+1+S[1][0])+F(4, 13)
    assert lower == F(3, 5)
    gain = (stationary-lower)/2
    assert gain == -F(7, 390)
    # Independent source 2 sleeps or sends an active to the stationary pair.
    K3 = [[F(0), F(1, 2), F(0), F(1, 2)],
          [F(1, 2), F(0), F(0), F(1, 2)],
          [F(1), F(0), F(0), F(0)]]
    initials = [(1, 1, 0), (1, 1, 1)]
    laws, _ = rn.stopped(initials, lambda state: rn.arw_step(K3, [F(1, 2)]*3, state),
                         lambda state: all(q <= 0 for q in state))
    margins = [sum(prob for C, prob in laws[state].items() if C[0] == -1)
               for state in initials]
    assert margins == [stationary, stationary]
    return dict(fixed_creator_forest_gain=str(gain), actual_gain="0")



def surplus_reset_normalization_audit():
    """Exact walk rejection normalization and labelled exit law."""
    from collections import defaultdict
    from math import prod
    examples = [
        ([[0, F(1,2), 0, F(1,2)],
          [F(1,2), 0, 0, F(1,2)], [0, 0, 0, 1]],
         [F(1,2)]*3, frozenset((0,1))),
        ([[F(1,4), F(1,4), F(1,4), F(1,4)],
          [F(1,3), F(1,3), 0, F(1,3)], [0, 0, 0, 1]],
         [F(1,3),F(2,5),F(1,2)], frozenset((0,1))),
        ([[0,F(1,3),F(1,3),F(1,3)],
          [F(1,4),0,F(1,4),F(1,2)],
          [F(1,5),F(1,5),F(1,5),F(2,5)]],
         [F(1,3),F(2,5),F(3,7)], frozenset((0,1,2))),
    ]
    checked = flags = 0
    for raw, p, R in examples:
        K = [[F(q) for q in row] for row in raw]
        n = len(p)
        sets = [frozenset(C) for C in subsets(R)]
        d = {C: determinant([[F(i == j)-K[i][j] for j in sorted(C)]
                             for i in sorted(C)]) for C in sets}
        def phi(U,A):
            return d[U]*sum(prod(p[i] if i in B else 1-p[i] for i in A)
                            /d[(U-A)|B] for B in sets if B <= A)
        normalizer = phi(R,R)
        U0 = sorted(R)
        G = rn.inverse([[F(i == j)-K[i][j] for j in U0] for i in U0])
        def step(state):
            z, le, counts = state
            newcounts = list(counts)
            newcounts[z] = min(2,newcounts[z]+1)
            out = defaultdict(F)
            for y,q in enumerate(K[z]):
                if not q: continue
                newle = le
                if y in R:
                    newle = le[:le.index(y)+1] if y in le else le+(y,)
                out[y,newle,tuple(newcounts)] += q
            return out
        starts = [(z,(z,),tuple(0 for _ in range(n))) for z in sorted(R)]
        laws,_ = rn.stopped(starts,step,lambda state:state[0] not in R)
        for start in starts:
            exit_weights = defaultdict(F)
            for (exit_label, le, counts), mass in laws[start].items():
                P = frozenset(le)
                L = frozenset(i for i in R if counts[i] > int(i in P))
                A = R-(P|L)
                U = R-P
                ph = phi(U,A)
                accept = prod(p[i] for i in L)*ph
                assert 0 <= accept <= 1
                exit_weights[exit_label] += mass*accept
                flag_weights = [prod(p[i] if i in B else 1-p[i] for i in A)
                                *d[U]/(d[(U-A)|B]*ph) for B in sets if B <= A]
                assert sum(flag_weights) == 1 and all(q >= 0 for q in flag_weights)
                if not A or (len(R)==2 and all(K[i][i]==0 for i in R)):
                    assert ph == 1
                flags += len(flag_weights)
                checked += 1
            assert sum(exit_weights.values()) == normalizer
            for y in range(n+1):
                if y in R:continue
                ordinary = sum(q*K[i][y] for i,q in zip(U0,G[U0.index(start[0])]))
                assert exit_weights[y]/normalizer == ordinary
    assert F(1,2)/(1-F(1,2)*F(1,2)) == F(2,3)
    return dict(surplus_walk_terminal_classes=checked,
                residual_flag_weights=flags, unvisited_first_bit="2/3")



def idla_addition_matrices(K, sets):
    n, m = len(K), len(sets)
    pos = {A: i for i, A in enumerate(sets)}
    D = [[[F(0)]*m for _ in range(m)] for _ in range(n)]
    for a, A in enumerate(sets):
        U = sorted(A)
        G = rn.inverse([[F(i == j)-K[i][j] for j in U] for i in U])
        for x in range(n):
            if x not in A:
                D[x][a][pos[A | {x}]] = F(1)
                continue
            row = G[U.index(x)]
            D[x][a][a] = sum(q*K[z][n] for z, q in zip(U, row))
            for y in range(n):
                if y not in A:
                    D[x][a][pos[A | {y}]] = sum(q*K[z][y] for z, q in zip(U, row))
    return D


def wake_field_coefficients(K, p, sets, P):
    """Return exact coefficients; their nonnegativity is not assumed."""
    from math import prod
    n, m = len(p), len(sets)
    V = frozenset(range(n))
    D = idla_addition_matrices(K, sets)
    Fs = [[[D[x][i][j]-sum(K[x][y]*D[y][i][j] for y in range(n))
            -K[x][n]*F(i == j) for j in range(m)] for i in range(m)]
          for x in range(n)]
    moments = [[sum(q for B, q in zip(sets, row) if T <= B) for T in sets]
               for row in P]
    coefficient = {}
    for C in sets:
        row = [F(A == V-C) for A in sets]
        for x in sorted(C):
            row = rn.multiply([row], Fs[x])[0]
        coefficient[C] = rn.multiply([row], moments)[0]
    for j in range(m):
        assert sum(prod((1-p[i])/p[i] for i in C)*coefficient[C][j]
                   for C in sets) == 1
    return coefficient, D, Fs


def full_wake_archived_audit():
    from math import prod
    checks = complete = degenerate = one_target = joint_target = 0
    bad = []
    for ci, case in enumerate(rn.CASES):
        K = [[F(q) for q in row] for row in case["K"]]
        p = list(map(F, case["p"]))
        n = len(p)
        V = frozenset(range(n))
        sets = [frozenset(i for i, b in enumerate(bits) if b)
                for bits in case["subsets"]]
        P = [[F(q) for q in row] for row in case["stabilization_matrix"]]
        coeff, D, Fs = wake_field_coefficients(K, p, sets, P)
        dv = determinant([[F(i == j)-K[i][j] for j in range(n)] for i in range(n)])
        for C in sets:
            R = sorted(V-C)
            dr = determinant([[F(i == j)-K[i][j] for j in R] for i in R])
            assert coeff[C][sets.index(V)] == dv/dr*P[sets.index(V)][sets.index(V)]
        m = len(sets)
        I = [[F(i == j) for j in range(m)] for i in range(m)]
        total = [[F(0)]*m for _ in range(m)]
        nondegenerate = True
        for C in sets:
            R, U = V-C, sorted(V-C)
            G = rn.inverse([[F(i == j)-K[i][j] for j in U] for i in U])
            h = [F(0)]*n
            for i, row in zip(U, G):
                h[i] = sum(q*K[j][n] for j, q in zip(U, row))
            escape = {x: K[x][n]+sum(K[x][y]*h[y] for y in U) for x in C}
            lam = (-1)**len(C)*prod(h[x] for x in R)*prod(escape.values())
            mat = I
            for x in sorted(R):
                mat = rn.multiply(mat, D[x])
            for x in sorted(C):
                mat = rn.multiply(mat, Fs[x])
            for z in range(n):
                assert rn.multiply(mat, D[z]) == [[h[z]*a for a in row] for row in mat]
                checks += m*m
            assert rn.multiply(mat, mat) == [[lam*a for a in row] for row in mat]
            if lam:
                total = [[q+a/lam for q, a in zip(row1, row2)]
                         for row1, row2 in zip(total, mat)]
            else:
                nondegenerate = False
                degenerate += 1
            for T, a in zip(sets, coeff[C]):
                joint_target += 1
                one_target += int(len(T) == 1)
                if a < 0:
                    bad.append((ci, sorted(C), sorted(T), str(a)))
        if nondegenerate:
            assert total == I
            complete += 1
    return dict(character_eigen_entries=checks, complete_domains=complete,
                zero_lambda_sectors=degenerate, full_wake_one_target=one_target,
                full_wake_joint_target=joint_target, negative_field_coefficients=bad[:2])


def full_wake_sandpile_obstruction():
    Q = [[F(0), F(1,3), F(1,2), F(1,6)],
         [F(0), F(1,3), F(1,3), F(1,3)],
         [F(1,2), F(0), F(0), F(1,2)]]
    n = 3
    S = [[(F(i == j)+Q[i][j])/2 for j in range(n)] for i in range(n)]
    supports = [tuple(C) for C in subsets(range(n))]
    laws = toppling_laws(Q, [F(0)]*n)
    active = deposited_law(S, range(n), laws, supports)
    mixedS = [list(row) for row in S]
    mixedS[1] = [F(i == 1) for i in range(n)]
    mixed = deposited_law(mixedS, range(n), laws, supports)
    coefficient = sum(hi-lo for B, lo, hi in zip(supports, active, mixed) if 2 in B)
    assert coefficient == -F(9092,4124577)
    return dict(general_sandpile_wake_coefficient=str(coefficient))


def full_wake_screen():
    """Optional finite search and an independent raw mixed-state audit."""
    from random import Random
    from math import prod
    rng = Random(2026101612)
    checked = 0
    for ci in range(30):
        n = 3 if ci < 25 else 4
        sets = [frozenset(C) for C in subsets(range(n))]
        K = []
        for i in range(n):
            row = [0 if j == i or rng.randrange(4) == 0 else rng.randrange(1,7)
                   for j in range(n)]+[rng.randrange(1,4)]
            K.append([F(a,sum(row)) for a in row])
        p = [rng.choice([F(1,20),F(1,3),F(1,2),F(4,5),F(19,20)]) for i in range(n)]
        starts = [tuple(int(i in A) for i in range(n)) for A in sets]
        laws, _ = rn.stopped(starts, lambda state: rn.arw_step(K,p,state),
                             lambda state: all(a <= 0 for a in state))
        P = [[laws[start].get(tuple(-int(i in B) for i in range(n)),F(0))
              for B in sets] for start in starts]
        coeff, _, _ = wake_field_coefficients(K,p,sets,P)
        bad = [(sorted(C), sorted(T), str(a)) for C in sets
               for T,a in zip(sets,coeff[C]) if len(T) == 1 and a < 0]
        checked += n*len(sets)
        if bad:
            return dict(seed=2026101612, case=ci, negative=bad,
                        K=[[str(a) for a in row] for row in K], p=list(map(str,p)))

    rng = Random(2026101615)
    n = 5
    sets = [frozenset(C) for C in subsets(range(n))]
    K = []
    for i in range(n):
        row = [0 if j == i or rng.randrange(3) == 0 else rng.randrange(1,4)
               for j in range(n)]+[3]
        K.append([F(a,sum(row)) for a in row])
    p = [rng.choice([F(1,3),F(2,5),F(4,5)]) for _ in range(n)]
    binary = {A: tuple(int(i in A) for i in range(n)) for A in sets}
    mixed = {S: tuple(-1 if i in S else 1 for i in range(n)) for S in sets}
    laws, _ = rn.stopped(list(binary.values())+list(mixed.values()),
                         lambda state: rn.arw_step(K,p,state),
                         lambda state: all(a <= 0 for a in state))
    P = [[laws[binary[A]].get(tuple(-int(i in B) for i in range(n)),F(0))
          for B in sets] for A in sets]
    coeff, _, _ = wake_field_coefficients(K,p,sets,P)
    verified = 0
    bad = []
    for v in range(n):
        j = sets.index(frozenset((v,)))
        mixedq = {S: sum(q for end,q in laws[mixed[S]].items() if end[v] == -1)
                  for S in sets}
        for C in sets:
            direct = sum((-1)**len(C-S)*mixedq[S] for S in sets if S <= C)
            assert direct == prod((1-p[i])/p[i] for i in C)*coeff[C][j]
            if direct < 0:
                bad.append((sorted(C),v,str(direct)))
            verified += 1
    return dict(screen_one_target_entries=checked, screen_cases=30,
                direct_mixed_state_entries=verified, negative=bad[:2])

def singleton_field_audit():
    """First-contact vector cancellation and the no-visit Doob transform."""
    from math import prod
    from collections import defaultdict
    vector_entries = doob_entries = diagonal_checks = 0
    for case in rn.CASES:
        K = [[F(q) for q in row] for row in case["K"]]
        p = list(map(F,case["p"]))
        n = len(p)
        V = frozenset(range(n))
        sets = [frozenset(i for i,b in enumerate(bits) if b)
                for bits in case["subsets"]]
        P = [[F(q) for q in row] for row in case["stabilization_matrix"]]
        pos = {A:i for i,A in enumerate(sets)}
        pi = P[pos[V]]
        for x in range(n):
            R = V-{x}
            U = sorted(R)
            G = rn.inverse([[F(i == j)-K[i][j] for j in U] for i in U])
            h = {i:sum(q*K[j][n] for j,q in zip(U,row)) for i,row in zip(U,G)}
            sigma = K[x][n]+sum(K[x][i]*h[i] for i in U)
            mu = P[pos[R]]
            start = tuple(int(i in R) for i in range(n))
            killed,_ = rn.stopped([start],lambda state:rn.arw_step(K,p,state),
                                  lambda state:state[x] > 0 or all(a <= 0 for a in state))
            nu = defaultdict(F)
            for state,q in killed[start].items():
                if state[x] == 0:
                    nu[frozenset(i for i,a in enumerate(state) if a == -1)] += q
            starts = {B:tuple(1 if i == x else -int(i in B) for i in range(n))
                      for B in sets if B <= R}
            laws,_ = rn.stopped(starts.values(),lambda state:rn.arw_step(K,p,state),
                                lambda state:all(a <= 0 for a in state))
            rhs = [F(0)]*len(sets)
            for B,q in nu.items():
                rhs[pos[B | {x}]] += q
                for state,mass in laws[starts[B]].items():
                    C = frozenset(i for i,a in enumerate(state) if a == -1)
                    rhs[pos[C]] -= q*mass
            expected = [(1-p[x])/p[x]*sigma*(a-b) for a,b in zip(pi,mu)]
            assert rhs == expected
            vector_entries += len(sets)
            bx = sigma*sum(a-b for B,a,b in zip(sets,pi,mu) if x in B)
            assert bx > 0 and bx >= K[x][n]*prod(p)
            diagonal_checks += 1
            if all(h[i] > 0 for i in U):
                Kh = [[K[i][j]*h[j]/h[i] for j in U]+[K[i][n]/h[i]] for i in U]
                ph = [p[i] for i in U]
                initial = (1,)*len(U)
                transformed,_ = rn.stopped([initial],lambda state:rn.arw_step(Kh,ph,state),
                                           lambda state:all(a <= 0 for a in state))
                for B in sets:
                    if not B <= R:continue
                    terminal = tuple(-int(i in B) for i in U)
                    value = prod(h.values())/prod(h[i] for i in B)*transformed[initial].get(terminal,F(0))
                    assert value == nu.get(B,F(0))
                    doob_entries += 1
    assert F(1,4)*(-F(2,15)) == -F(1,30)
    return dict(singleton_vector_entries=vector_entries,
                singleton_diagonal_checks=diagonal_checks,
                no_visit_doob_entries=doob_entries,
                negative_individual_summand="-1/30")

def full_block_exit_audit():
    """Full occupation with mixed initial modes and weighted labelled exits."""
    from collections import defaultdict
    from math import prod
    checked = 0
    for n in (2,3):
        p = [F(1,3),F(2,5),F(3,7)][:n]
        K = [[F(0) if i == j else F(1,5) for j in range(n)]
             +[F(1,5),1-F(n,5)] for i in range(n)]
        V = frozenset(range(n))
        sets = [frozenset(C) for C in subsets(range(n))]
        d = {R:determinant([[F(i == j)-K[i][j] for j in sorted(R)]
                            for i in sorted(R)]) for R in sets}
        odds = [(1-q)/q for q in p]
        Z = d[V]*sum(prod(odds[i] for i in C)/d[V-C] for C in sets)
        theta = {C:prod(odds[i] for i in C)*d[V]/(d[V-C]*Z) for C in sets}
        assert sum(theta.values()) == 1
        extras = [(0,)*n]+[tuple(int(i == z) for i in range(n)) for z in range(n)]
        extras.append(tuple(int(i in (0,1)) for i in range(n)))
        for weights in ((F(1,2),F(2,3)),(F(0),F(1)),(F(0),F(0))):
            h = {}
            for C in sets:
                R = sorted(V-C)
                G = rn.inverse([[F(i == j)-K[i][j] for j in R] for i in R])
                h[C] = [F(0)]*n
                for i,row in zip(R,G):
                    h[C][i] = sum(q*sum(K[j][n+e]*weights[e] for e in range(2))
                                  for j,q in zip(R,row))
            fail = ("fail",)
            def step(state):
                z = next(i for i,a in enumerate(state) if a > 0)
                out = defaultdict(F)
                if state[z] == 1:
                    end = list(state);end[z] = -1
                    out[tuple(end)] += p[z]
                jump = (1-p[z]) if state[z] == 1 else F(1)
                for y,q in enumerate(K[z]):
                    if not q:continue
                    end = list(state);end[z] -= 1
                    if y < n:
                        end[y] = 2 if end[y] == -1 else end[y]+1
                        out[tuple(end)] += jump*q
                    else:
                        out[tuple(end)] += jump*q*weights[y-n]
                        out[fail] += jump*q*(1-weights[y-n])
                return {st:q for st,q in out.items() if q}
            starts = {}
            for S in sets:
                for eta in extras:
                    state = tuple(1+eta[i] if eta[i] or i not in S else -1
                                  for i in range(n))
                    starts[S,eta] = state
            laws,_ = rn.stopped(starts.values(),step,
                                lambda state:state == fail or all(a <= 0 for a in state))
            for (S,eta),state in starts.items():
                value = sum(theta[C]*prod(h[C][i]**eta[i] for i in range(n))
                            for C in sets if C <= S)
                assert value == laws[state].get((-1,)*n,F(0))
                checked += 1
    return dict(full_block_labelled_exit_entries=checked)

def dag_wake_field_audit():
    """Positive backward field and coherent random-graph dual on DAGs."""
    from collections import defaultdict
    from math import prod
    from random import Random
    rng = Random(2026101616)
    checked = graph_draws = 0
    for n in (3,4):
        V = frozenset(range(n))
        sets = [frozenset(C) for C in subsets(range(n))]
        K = []
        for i in range(n):
            row = [rng.randrange(1,4) if j > i else 0 for j in range(n)]
            row += [rng.randrange(1,3),rng.randrange(1,3)]
            K.append([F(a,sum(row)) for a in row])
        p = [rng.choice([F(1,3),F(2,5),F(4,5)]) for _ in range(n)]
        def dual(T, weights):
            result, potentials = {}, {}
            for C in sets:
                h = [F(0)]*n
                mass = F(1)
                for i in reversed(range(n)):
                    ell = sum(K[i][j]*h[j] for j in range(i+1,n))
                    ell += sum(K[i][n+e]*weights[e] for e in range(2))
                    alpha = p[i] if i in T else p[i]+(1-p[i])*ell
                    beta = 1-p[i] if i in T else (1-p[i])*(1-ell)
                    assert alpha+beta == 1 and min(alpha,beta) >= 0
                    mass *= beta if i in C else alpha
                    h[i] = F(0) if i in C else ell
                result[C], potentials[C] = mass, h
            assert sum(result.values()) == 1
            return result,potentials

        extras = [(0,)*n]+[tuple(int(i == z) for i in range(n)) for z in range(n)]
        extras.append(tuple(int(i in (0,1)) for i in range(n)))
        for weights in ((F(1),F(1)),(F(1,2),F(2,3)),(F(0),F(1))):
            fail = ("fail",)
            def step(state):
                z = next(i for i,a in enumerate(state) if a > 0)
                out = defaultdict(F)
                if state[z] == 1:
                    end = list(state);end[z] = -1
                    out[tuple(end)] += p[z]
                jump = (1-p[z]) if state[z] == 1 else F(1)
                for y,q in enumerate(K[z]):
                    if not q:continue
                    end = list(state);end[z] -= 1
                    if y < n:
                        end[y] = 2 if end[y] == -1 else end[y]+1
                        out[tuple(end)] += jump*q
                    else:
                        out[tuple(end)] += jump*q*weights[y-n]
                        out[fail] += jump*q*(1-weights[y-n])
                return {st:q for st,q in out.items() if q}
            starts = {(S,eta):tuple(1+eta[i] if eta[i] or i not in S else -1
                                   for i in range(n)) for S in sets for eta in extras}
            laws,_ = rn.stopped(starts.values(),step,
                                lambda state:state == fail or all(a <= 0 for a in state))
            for T in sets:
                theta,h = dual(T,weights)
                for (S,eta),state in starts.items():
                    value = sum(theta[C]*prod(h[C][i]**eta[i] for i in range(n))
                                for C in sets if C <= S)
                    raw = sum(q for end,q in laws[state].items()
                              if end != fail and all(end[i] == -1 for i in T))
                    assert value == raw
                    checked += 1

        def walks(i):
            result = [(frozenset(),sum(K[i][n:]))]
            for j in range(i+1,n):
                for visited,mass in walks(j):
                    if K[i][j]*mass:
                        result.append((visited | {j},K[i][j]*mass))
            return result
        empirical = {T:defaultdict(F) for T in sets}
        for A in sets:
            coinweight = prod(1-p[i] if i in A else p[i] for i in range(n))
            for paths in product(*(walks(i) for i in range(n))):
                mass = coinweight*prod(q for _,q in paths)
                if not mass:continue
                H = {}
                for T in sets:
                    marked = set()
                    for i in reversed(range(n)):
                        if i in A and (i in T or paths[i][0] & marked):
                            marked.add(i)
                    H[T] = frozenset(marked)
                    empirical[T][H[T]] += mass
                for T in sets:
                    for U in sets:
                        assert H[T | U] == H[T] | H[U]
                graph_draws += 1
        for T in sets:
            theta,_ = dual(T,(F(1),F(1)))
            assert all(empirical[T].get(C,F(0)) == theta[C] for C in sets)
    return dict(dag_weighted_field_entries=checked, coherent_graph_draws=graph_draws)


def reset_contact_obstruction():
    """Root-loop averaging does not remove the lower unread outside path."""
    from collections import defaultdict
    K = [[F(0),F(1,2),F(1,2),F(0)],
         [F(1),F(0),F(0),F(0)],
         [F(1,2),F(0),F(0),F(1,2)]]
    p = [F(1,2)]*3
    roots, full, sleeping = (1,1,0),(1,1,1),(1,1,-1)
    laws,_ = rn.stopped([roots,full,sleeping],lambda state:rn.arw_step(K,p,state),
                        lambda state:all(q<=0 for q in state))
    target = lambda law:sum(q for end,q in law.items() if end[0] == -1)
    base,upper,mixed = [target(laws[start]) for start in (roots,full,sleeping)]
    def first_record(bit,arrow=None):
        def step(state):
            counts,used = state
            z = next(i for i,a in enumerate(counts) if a>0)
            if z != 2 or used:
                return {(end,used):q for end,q in rn.arw_step(K,p,counts).items()}
            out=defaultdict(F)
            if counts[z] == 1 and bit:
                end=list(counts);end[z]=-1
                out[tuple(end),False]+=1
            else:
                destinations=[(arrow,F(1))] if arrow is not None else enumerate(K[z])
                for y,q in destinations:
                    if not q:continue
                    end=list(counts);end[z]-=1
                    if y<3:end[y]=2 if end[y]==-1 else end[y]+1
                    out[tuple(end),True]+=q
            return out
        start=(roots,False)
        laws,_=rn.stopped([start],step,lambda state:all(q<=0 for q in state[0]))
        return sum(q for (end,_),q in laws[start].items() if end[0] == -1)
    low_contact=first_record(0,0)
    low_sleep=first_record(1)
    low_sink=first_record(0,3)
    assert low_contact == F(10382,12675)
    assert upper == F(70459,93795)
    contact=upper-low_contact
    assert contact == -F(10613,156325)
    actual=upper-base
    assert actual == F(4414,93795)
    assert base == F(1,2)*low_sleep+F(1,4)*(low_contact+low_sink)
    no_contact=F(1,2)*(mixed-low_sleep)+F(1,4)*(base-low_sink)
    assert actual == no_contact+F(1,4)*contact
    assert no_contact > 0
    return dict(reset_contact_gain=str(contact), actual_gain=str(actual),
                no_contact_credit=str(no_contact),
                contact_probability="1/4", common_loop_sector_weights=["4/5","1/5"])


def two_site_wake_field(a, b, p, target):
    """Closed positive coefficients for a two-site block without self-loops."""
    from math import prod
    V = frozenset((0,1))
    sets = [frozenset(C) for C in subsets(range(2))]
    target = frozenset(target)
    c, d = a*b, 1-a*b
    assert 0 <= a <= 1 and 0 <= b <= 1 and d > 0
    if not target:
        return {C:F(not C) for C in sets}
    p0,p1 = p
    q0,q1 = 1-p0,1-p1
    ds = 1-c*q0*q1
    df = 1-c*(1-p0*p1)
    if target == V:
        return {frozenset():p0*p1/df,
                frozenset((0,)):q0*p1*d/df,
                frozenset((1,)):p0*q1*d/df,
                V:q0*q1*d/df}
    if target == frozenset((1,)):
        swapped = two_site_wake_field(b,a,[p1,p0],(0,))
        return {frozenset(1-i for i in C):q for C,q in swapped.items()}
    t = 1-b*q1
    ans = {
        frozenset():p0/ds+p0*p1*c*q0*t/(ds*df),
        frozenset((0,)):q0*d*t*(1-c*q1)/(ds*df),
        frozenset((1,)):p0*q0*q1*c*d*t/(ds*df),
        V:q0*q1*(b*(1-a)*df+c*q0*d*t)/(ds*df),
    }
    assert sum(ans.values()) == 1 and all(q >= 0 for q in ans.values())
    return ans


def two_site_field_audit():
    from random import Random
    rng = Random(2026101617)
    checked = 0
    for k in range(24):
        a = rng.choice([F(0),F(1,4),F(1,2),F(3,4),F(1)])
        b = rng.choice([F(0),F(1,3),F(2,3),F(1)])
        if a*b == 1:b=F(2,3)
        p = [rng.choice([F(1,5),F(1,2),F(4,5)]) for _ in range(2)]
        K = [[F(0),a,1-a],[b,F(0),1-b]]
        sets = [frozenset(C) for C in subsets(range(2))]
        starts = {S:tuple(-1 if i in S else 1 for i in range(2)) for S in sets}
        laws,_ = rn.stopped(starts.values(),lambda state:rn.arw_step(K,p,state),
                            lambda state:all(q <= 0 for q in state))
        for T in sets:
            theta = two_site_wake_field(a,b,p,T)
            moment = {S:sum(q for end,q in laws[starts[S]].items()
                            if all(end[i] == -1 for i in T)) for S in sets}
            for C in sets:
                assert theta[C] == sum((-1)**len(C-S)*moment[S] for S in sets if S <= C)
                checked += 1
    return dict(two_site_field_coefficients=checked)


def scc_pair_field_audit():
    """Boundary-target augmentation and reverse composition of two SCCs."""
    from collections import defaultdict
    from math import prod
    K = [[0,F(3,10),F(1,5),0,F(1,5),F(3,10)],
         [F(2,5),0,0,F(1,5),F(1,10),F(3,10)],
         [0,0,0,F(1,2),F(1,5),F(3,10)],
         [0,0,F(1,4),0,F(1,4),F(1,2)]]
    K = [[F(q) for q in row] for row in K]
    p = [F(1,3),F(2,5),F(3,7),F(4,5)]
    n = 4
    V = frozenset(range(n))
    sets = [frozenset(C) for C in subsets(range(n))]
    localsets = [frozenset(C) for C in subsets(range(2))]
    extras = [(0,)*n]+[tuple(int(i == z) for i in range(n)) for z in range(n)]
    extras.append((1,1,0,0))
    checked = 0
    for marks in ((F(1),F(1)),(F(1,2),F(2,3)),(F(0),F(0))):
        def local(block, T, boundary):
            a,b = K[block[0]][block[1]],K[block[1]][block[0]]
            pp = [p[i] for i in block]
            G = rn.inverse([[F(i == j)-K[block[i]][block[j]] for j in range(2)]
                            for i in range(2)])
            rhs = [sum(K[i][y]*boundary[y] for y in boundary) for i in block]
            h = [sum(q*w for q,w in zip(row,rhs)) for row in G]
            target = frozenset(j for j,i in enumerate(block) if i in T)
            if not any(h):
                theta = two_site_wake_field(a,b,pp,(0,1))
            else:
                assert min(h)>0
                ah,bh = a*h[1]/h[0],b*h[0]/h[1]
                theta = {C:F(0) for C in localsets}
                for D in localsets:
                    if D & target:continue
                    weight = prod(1-h[i] if i in D else h[i] for i in range(2) if i not in target)
                    sub = two_site_wake_field(ah,bh,pp,target | D)
                    for C,q in sub.items():theta[C] += weight*q
                assert sum(theta.values()) == 1
            potentials = {}
            for C in localsets:
                U = sorted(set(range(2))-C)
                Gc = rn.inverse([[F(i == j)-K[block[i]][block[j]] for j in U] for i in U])
                hc = [F(0)]*2
                for i,row in zip(U,Gc):hc[i]=sum(q*rhs[j] for j,q in zip(U,row))
                potentials[C] = hc
            return theta,potentials
        alltheta = {}
        potentials = {}
        for T in sets:
            theta = {}
            low,lh = local((2,3),T,{4:marks[0],5:marks[1]})
            for C2,w2 in low.items():
                boundary = {2:lh[C2][0],3:lh[C2][1],4:marks[0],5:marks[1]}
                high,hh = local((0,1),T,boundary)
                for C1,w1 in high.items():
                    C = C1 | frozenset(i+2 for i in C2)
                    theta[C] = w1*w2
                    value = hh[C1]+lh[C2]
                    if C in potentials:assert potentials[C] == value
                    potentials[C] = value
            assert sum(theta.values()) == 1 and all(q>=0 for q in theta.values())
            alltheta[T] = theta
        fail = ("fail",)
        def step(state):
            z = next(i for i,a in enumerate(state) if a>0)
            out = defaultdict(F)
            if state[z] == 1:
                end=list(state);end[z]=-1;out[tuple(end)]+=p[z]
            jump=1-p[z] if state[z]==1 else F(1)
            for y,q in enumerate(K[z]):
                if not q:continue
                end=list(state);end[z]-=1
                if y<n:
                    end[y]=2 if end[y]==-1 else end[y]+1
                    out[tuple(end)]+=jump*q
                else:
                    out[tuple(end)]+=jump*q*marks[y-n]
                    out[fail]+=jump*q*(1-marks[y-n])
            return {s:q for s,q in out.items() if q}
        starts={(S,eta):tuple(1+eta[i] if eta[i] or i not in S else -1 for i in range(n))
                for S in sets for eta in extras}
        laws,_=rn.stopped(starts.values(),step,lambda state:state==fail or all(a<=0 for a in state))
        for T,theta in alltheta.items():
            for (S,eta),start in starts.items():
                predicted=sum(theta[C]*prod(potentials[C][i]**eta[i] for i in range(n))
                              for C in sets if C<=S)
                raw=sum(q for end,q in laws[start].items() if end!=fail and all(end[i]==-1 for i in T))
                assert predicted == raw,(T,S,eta,predicted,raw)
                checked+=1
    return dict(two_scc_weighted_field_entries=checked)


def coherent_field_certificates():
    """Two exact finite union-coherent laws; no general cyclic construction."""
    from collections import defaultdict
    half = [[F(0),F(1,2),F(1,2)],[F(1,2),F(0),F(1,2)]]
    small = [[120,0,81,19],[81,9,20,16],[0,0,9,0],[19,0,16,0]]
    smalllaw = {(i,j):F(q,390) for i,row in enumerate(small)
                for j,q in enumerate(row) if q}
    data = """
000:7245328757904 001:200478698464 002:6882141565 003:6220743192
004:1064741081376 005:18396786000 006:58660687523 007:66887592728
020:1816188378858 024:108879591960 026:166682762296 030:455135220540
037:112257786012 044:22893292240 060:61229663640 064:9454853556
070:77459130144 077:13552472320 100:7912067238936 104:526137443208
105:484665351696 110:229955463864 111:110382149624 116:3676789552
120:665687897100 122:3888107040 123:15593248037 124:48181620631
130:1213513297974 134:12198877242 135:75689392956 220:213841409844
222:76308957455 224:8002700613 300:637334594496 303:65082324715
306:66666345206 320:107375898489 404:70739938080 440:37964311288
444:55940421592 500:448022710080 504:315655591560 550:85540684320
560:21644391339 570:31530541466 606:25831296621 660:30863085507
700:135714839040 706:38947183938
"""
    largelaw = {tuple(map(int,key)):F(int(value),25285973743827)
                for key,value in (word.split(":") for word in data.split())}
    case = rn.CASES[3]
    largeK = [[F(q) for q in row] for row in case["K"]]
    largep = list(map(F,case["p"]))
    checked = 0
    for K,p,certificate in [(half,[F(1,2)]*2,smalllaw),(largeK,largep,largelaw)]:
        n=len(p)
        masks=range(1<<n)
        starts={S:tuple(-1 if S & (1<<i) else 1 for i in range(n)) for S in masks}
        laws,_=rn.stopped(starts.values(),lambda state:rn.arw_step(K,p,state),
                          lambda state:all(q<=0 for q in state))
        assert sum(certificate.values()) == 1 and min(certificate.values()) > 0
        for T in masks:
            moment={S:sum(q for end,q in laws[start].items()
                          if all(end[i]==-1 for i in range(n) if T & (1<<i)))
                    for S,start in starts.items()}
            empirical=defaultdict(F)
            for relation,q in certificate.items():
                union=0
                for i in range(n):
                    if T & (1<<i):union |= relation[i]
                empirical[union]+=q
            for C in masks:
                theta=sum((-1)**(C.bit_count()-S.bit_count())*moment[S]
                          for S in masks if S & C == S)
                assert theta == empirical[C]
                checked += int(T != 0)
    return dict(finite_coherent_field_constraints=checked, cyclic_cases=2)


def creator_mixed_field_obstruction():
    """Full-root conditioning does not preserve arbitrary mixed initial flags."""
    from collections import defaultdict
    K=[[F(0),F(1,2),F(1,2)],[F(1,2),F(0),F(1,2)]]
    p=[F(1,2)]*2
    sets=[frozenset(C) for C in subsets(range(2))]
    starts={S:tuple(-1 if i in S else 1 for i in range(2)) for S in sets}
    raw,_=rn.stopped(starts.values(),lambda s:rn.arw_step(K,p,s),
                    lambda s:all(a<=0 for a in s))
    nu={C:F(4 if len(C)==2 else 3,13) for C in sets}
    out={S:defaultdict(F) for S in sets}
    for C in sets:
        if len(C)==2:
            for S in sets:out[S][frozenset((0,1))]+=nu[C]
            continue
        for arrows in [(2,2),(1,2),(2,0)]:
            def step(state):
                counts,used=state
                z=next(i for i,a in enumerate(counts) if a>0)
                if used & (1<<z):
                    return {(end,used):q for end,q in rn.arw_step(K,p,counts).items()}
                end=list(counts)
                if counts[z]==1 and z in C:
                    end[z]=-1
                    return {(tuple(end),used):F(1)}
                end[z]-=1
                y=arrows[z]
                if y<2:end[y]=2 if end[y]==-1 else end[y]+1
                return {(tuple(end),used | (1<<z)):F(1)}
            init={S:(counts,0) for S,counts in starts.items()}
            laws,_=rn.stopped(init.values(),step,lambda s:all(a<=0 for a in s[0]))
            for S,start in init.items():
                for (end,_),q in laws[start].items():
                    out[S][frozenset(i for i,a in enumerate(end) if a==-1)]+=nu[C]*F(1,3)*q
    physical={S:{frozenset(i for i,a in enumerate(end) if a==-1):q
                 for end,q in raw[start].items()} for S,start in starts.items()}
    assert all(out[frozenset()][B] == physical[frozenset()].get(B,0) for B in sets)
    S=frozenset((0,))
    assert [out[S][B] for B in sets] == [F(2,39),F(76,195),F(4,195),F(7,13)]
    assert [physical[S].get(B,0) for B in sets] == [F(5,78),F(68,195),F(19,390),F(7,13)]
    target=lambda law:sum(q for B,q in law.items() if 1 in B)
    false_coefficient=target(out[S])-target(out[frozenset()])
    assert false_coefficient == -F(1,195)
    return dict(unchanged_creator_mixed_field_coefficient=str(false_coefficient),
                correct_coefficient=str(target(physical[S])-target(physical[frozenset()])))


def target_generator_obstruction():
    """Two-site endpoint dual is not a homogeneous subset Markov evolution."""
    theta=[[F(1),F(0),F(0),F(0)],
           [F(22,39),F(21,65),F(3,130),F(7,78)],
           [F(22,39),F(3,130),F(21,65),F(7,78)],
           [F(4,13),F(3,13),F(3,13),F(3,13)]]
    inverse=rn.inverse(theta)
    stopped=[[F(i==j)-inverse[i][j] for j in range(4)] for i in range(4)]
    expected=[[0,0,0,0],[F(8,3),-F(11,3),-F(4,3),F(7,3)],
              [F(8,3),-F(4,3),-F(11,3),F(7,3)],[-4,6,6,-8]]
    assert stopped == expected
    # Exact spectral coefficients of log(theta)[V,empty].
    eigenvalues=[F(1),F(1,2),F(3,10),F(1,13)]
    weights={}
    for lam in eigenvalues:
        projector=[[F(i==j) for j in range(4)] for i in range(4)]
        for other in eigenvalues:
            if other == lam:continue
            factor=[[(theta[i][j]-other*F(i==j))/(lam-other) for j in range(4)]
                    for i in range(4)]
            projector=rn.multiply(projector,factor)
        assert rn.multiply(projector,theta) == [[lam*q for q in row] for row in projector]
        weights[lam]=projector[3][0]
    assert weights == {F(1):F(1),F(1,2):-F(16,11),F(3,10):F(0),F(1,13):F(5,11)}
    assert 2**16 < 13**5
    return dict(current_set_stopping_rate_entry=str(stopped[3][0]),
                fixed_time_log_entry="(16 log(2)-5 log(13))/11 < 0")


def marked_site_deletion_audit():
    """Exact signed transfer of marked field coefficients to a smaller graph."""
    from math import prod
    checked=0
    witnesses=[]
    for n in (2,3):
        K=[[F(1,2) if j == (i+1)%n or j == n else F(0)
            for j in range(n+1)] for i in range(n)]
        p=[F(1,2)]*n
        x=0
        R=tuple(range(1,n))
        target=frozenset((n-1,))
        subsetsR=[frozenset(C) for C in subsets(R)]
        G=rn.inverse([[F(i==j)-K[i][j] for j in R] for i in R])
        h={i:sum(q*K[j][n] for j,q in zip(R,row)) for i,row in zip(R,G)}
        Kh=[[K[i][j]*h[j]/h[i] for j in R]+[K[i][n]/h[i]] for i in R]
        ph=[p[i] for i in R]
        starts={B:tuple(1 if i==x else -int(i in B) for i in range(n)) for B in subsetsR}
        added,_=rn.stopped(starts.values(),lambda state:rn.arw_step(K,p,state),
                           lambda state:all(q<=0 for q in state))
        tilted={}
        for B,start in starts.items():
            value=sum(q for end,q in added[start].items() if all(end[i]==-1 for i in target))
            tilted[B]=(F(target <= B | {x})-value)/prod(h[i] for i in B)
        kappa={D:sum((-1)**len(D-B)*tilted[B] for B in subsetsR if B <= D)
               for D in subsetsR}

        def fields(kernel,sleep,labels):
            sets=[frozenset(C) for C in subsets(labels)]
            initial={S:tuple(-1 if i in S else 1 for i in labels) for S in sets}
            laws,_=rn.stopped(initial.values(),lambda state:rn.arw_step(kernel,sleep,state),
                              lambda state:all(q<=0 for q in state))
            coeff={}
            for T in sets:
                moment={S:sum(q for end,q in laws[start].items()
                              if all(end[j]==-1 for j,i in enumerate(labels) if i in T))
                        for S,start in initial.items()}
                coeff[T]={C:sum((-1)**len(C-S)*moment[S] for S in sets if S <= C)
                          for C in sets}
            return coeff
        full=fields(K,p,tuple(range(n)))
        smaller=fields(Kh,ph,R)
        for C in subsetsR:
            rhs=prod(h.values())*sum(kappa[D]*smaller[D][C] for D in subsetsR)
            assert full[target][C | {x}] == rhs
            checked+=1
        if n==2:
            assert kappa[frozenset()] == -F(2,15)
            assert prod(h.values())*kappa[frozenset((1,))] == F(7,39)
            assert full[target][frozenset((0,))] == F(3,130)
            assert full[target][frozenset((0,1))] == F(7,78)
        else:
            assert kappa[frozenset((1,))] == -F(54646,725571)
            assert full[target][frozenset((0,1))] == F(12211,2042348)
        witnesses.append(str(kappa[frozenset()] if n==2 else kappa[frozenset((1,))]))
    return dict(marked_site_deletion_entries=checked, negative_target_coefficients=witnesses)


if __name__ == "__main__":
    print(marked_site_deletion_audit(), flush=True)
    print(target_generator_obstruction(), flush=True)
    print(coherent_field_certificates(), flush=True)
    print(creator_mixed_field_obstruction(), flush=True)
    print(reset_contact_obstruction(), flush=True)
    print(two_site_field_audit(), flush=True)
    print(scc_pair_field_audit(), flush=True)
    print(full_wake_archived_audit(), flush=True)
    print(full_wake_sandpile_obstruction(), flush=True)
    print(singleton_field_audit(), flush=True)
    print(full_block_exit_audit(), flush=True)
    print(dag_wake_field_audit(), flush=True)
    if "--wake-screen" in __import__("sys").argv:
        print(full_wake_screen(), flush=True)

    print(surplus_reset_normalization_audit(), flush=True)

    print(periodic_loop_audit(), flush=True)
    print(common_retained_loop_obstruction(), flush=True)
    print(terminal_pair_flow_audit(), flush=True)
    print(marker_occupation_obstruction(), flush=True)
    print(loop_thinning_coefficient_audit(), flush=True)
    print(creator_exit_forest_obstruction(), flush=True)
    from verify_backward_kernels import cases
    print(toppling_and_partition_obstructions(), flush=True)
    print(full_inverse_audit(), flush=True)
    print(causal_charge_obstruction(), flush=True)
    print(inverse_diagonal_audit(), flush=True)
    for name, K, p in cases():
        print(name, audit(K, p), flush=True)
