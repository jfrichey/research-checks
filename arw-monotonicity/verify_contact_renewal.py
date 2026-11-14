"""Positive three-particle contact renewal, and its collision-kernel identity.

Exact Fraction calculations; standard library only. This verifies identities,
not the open add-one sign. Written by Codex during the ARW research task.
No files are written when this script runs.
"""

from fractions import Fraction as F
from itertools import combinations, combinations_with_replacement, product
import verify_response_renewal as rn


def contact_model(K, p):
    n = len(p)
    supports = [tuple(z for z, bit in enumerate(bits) if bit)
                for bits in product((0, 1), repeat=n) if sum(bits) <= 3]
    pos = {B: i for i, B in enumerate(supports)}
    size = len(supports)

    def walk(taboo):
        outside = [z for z in range(n) if z not in taboo]
        G = rn.inverse([[F(i == j) - (1-p[i])*K[i][j]
                         for j in outside] for i in outside])
        result = {}
        for a in range(n):
            if a in taboo:
                result[a] = ({}, F(0), {z: F(z == a) for z in taboo})
                continue
            row = G[outside.index(a)]
            sleep = {z: row[j]*p[z] for j, z in enumerate(outside)}
            sink = sum(row[j]*(1-p[z])*K[z][n] for j, z in enumerate(outside))
            hit = {y: sum(row[j]*(1-p[z])*K[z][y]
                          for j, z in enumerate(outside)) for y in taboo}
            assert sum(sleep.values()) + sink + sum(hit.values()) == 1
            result[a] = sleep, sink, hit
        return result

    walks = {J: walk(set(J)) for k in range(3) for J in combinations(range(n), k)}
    S = {a: walks[()][a][0] for a in range(n)}
    e = {a: walks[()][a][1] for a in range(n)}
    h = {(a, w): S[a][w]/S[w][w] for a in range(n) for w in range(n)}
    single = {}
    for a in range(n):
        row = [F(0)] * size
        row[pos[()]] = e[a]
        for w, q in S[a].items():
            row[pos[(w,)]] = q
        single[a] = row

    def nocollision(a, b):
        row = [e[a]*q for q in single[b]]
        for w, q in S[a].items():
            sl, ex, _ = walks[(w,)][b]
            row[pos[(w,)]] += q*ex
            for u, weight in sl.items():
                row[pos[tuple(sorted((u, w)))]] += q*weight
        return row

    M = [[S[t][u]*h[t, u] for u in range(n)] for t in range(n)]
    invM = rn.inverse([[F(t == u)-M[t][u] for u in range(n)] for t in range(n)])
    two_same = rn.multiply(invM, [nocollision(t, t) for t in range(n)])

    for v in range(n):
        d = [sum(q for C, q in zip(supports, two_same[t]) if v in C) - S[t][v]
             for t in range(n)]
        assert all(q >= 0 for q in d)
        for t in range(n):
            assert d[t] == sum(q*z for q, z in zip(M[t], d)) + S[t][v]*(1-h[t, v])

    def two(a, b):
        row = nocollision(a, b)
        for u in range(n):
            weight = S[a][u]*h[b, u]
            row = [q + weight*z for q, z in zip(row, two_same[u])]
        return row

    mixed = {}
    for t in range(n):
        for w in range(n):
            sl, ex, hit = walks[(w,)][t]
            row = [hit[w]*q for q in two_same[w]]
            row[pos[(w,)]] += ex
            for u, q in sl.items():
                row[pos[tuple(sorted((u, w)))]] += q
            mixed[t, w] = row

    H = {(w, t): {u: S[t][u]+S[t][w]*S[w][u]/(1-S[w][w])
                   for u in range(n) if u != w}
         for w in range(n) for t in range(n)}
    E = {(w, t): e[t]+S[t][w]*e[w]/(1-S[w][w])
         for w in range(n) for t in range(n)}
    states = [(t, w) for t in range(n) for w in range(n) if t != w]
    index = {state: i for i, state in enumerate(states)}
    T = [[F(0)] * len(states) for _ in states]
    reward = [[F(0)] * size for _ in states]
    cross = {}
    identity_checks = 0
    for i, (t, w) in enumerate(states):
        sl, ex, hit = walks[(w,)][t]
        for u, q in H[w, t].items():
            cross[t, w, u] = hit[w]*q
        reward[i] = [hit[w]*E[w, t]*q + ex*z
                     for q, z in zip(two_same[w], mixed[t, w])]
        for u, q in sl.items():
            sl2, ex2, hit2 = walks[tuple(sorted((u, w)))][t]
            T[i][index[u, w]] += q*hit2[u]
            cross[t, w, u] += q*hit2[w]
            reward[i][pos[tuple(sorted((u, w)))]] += q*ex2
            for z, weight in sl2.items():
                reward[i][pos[tuple(sorted((u, w, z)))]] += q*weight
        for u in range(n):
            if u == w:
                continue
            T[i][index[w, u]] += cross[t, w, u]
            sl_u = walks[(w,)][u][0]
            killed_M = sl[u]*sl[u]/sl_u[u]
            traced_M = M[t][u] + M[t][w]*M[w][u]/(1-M[w][w])
            assert cross[t, w, u]*h[w, u] == traced_M-killed_M
            identity_checks += 1
        assert all(q >= 0 for q in T[i]+reward[i])
        assert sum(T[i]) + sum(reward[i]) == 1
        assert sum(T[i]) < 1
    invT = rn.inverse([[F(i == j)-T[i][j] for j in range(len(states))]
                       for i in range(len(states))])
    B = rn.multiply(invT, reward)

    def triple_pair(t, x):
        row = [E[t, x]*q for q in two_same[t]]
        for w, q in H[t, x].items():
            row = [a+q*b for a, b in zip(row, B[index[t, w]])]
        return row

    def two_with_sleeper(a, b, w):
        sl, ex, hit = walks[(w,)][a]
        row = [hit[w]*q + ex*z
               for q, z in zip(triple_pair(w, b), mixed[b, w])]
        for u, q in sl.items():
            sl2, ex2, hit2 = walks[tuple(sorted((u, w)))][b]
            row = [r+q*(hit2[w]*x+hit2[u]*y)
                   for r, x, y in zip(row, B[index[w, u]], B[index[u, w]])]
            row[pos[tuple(sorted((u, w)))]] += q*ex2
            for z, weight in sl2.items():
                row[pos[tuple(sorted((u, w, z)))]] += q*weight
        return row

    triples = {}
    for a, b, c in combinations_with_replacement(range(n), 3):
        row = [e[c]*q for q in two(a, b)]
        for w, q in S[c].items():
            row = [r+q*z for r, z in zip(row, two_with_sleeper(a, b, w))]
        triples[a, b, c] = row
    return supports, states, M, T, reward, two_same, B, triples, identity_checks


def creator_flow(K, p, data=None):
    """Audit the localized debt identity; its terminal sign is not assumed."""
    if data is None:
        data = contact_model(K, p)
    supports, states, C, T, _, two_same, B, _, _ = data
    n = len(p)
    G = rn.inverse([[F(i == j)-(1-p[i])*K[i][j] for j in range(n)]
                    for i in range(n)])
    S = [[G[t][v]*p[v] for v in range(n)] for t in range(n)]
    h = {(t, w): S[t][w]/S[w][w] for t in range(n) for w in range(n)}
    U = [[sum(q for support, q in zip(supports, row) if v in support)
          for v in range(n)] for row in two_same]
    D = [[U[t][v]-S[t][v] for v in range(n)] for t in range(n)]
    assert all(q >= 0 for row in D for q in row)
    assert all(sum(row) < 1 for row in D)
    invC = rn.inverse([[F(t == u)-C[t][u] for u in range(n)] for t in range(n)])
    assert D == rn.multiply(invC, [[S[t][u]-C[t][u] for u in range(n)]
                                  for t in range(n)])
    A = [[F(u == w)*h[t, w]*(1-h[t, w]) for u in range(n)] for t, w in states]
    invT = rn.inverse([[F(i == j)-T[i][j] for j in range(len(states))]
                       for i in range(len(states))])
    J = rn.multiply(invT, A)
    L = [[q+F(u == w)*h[t, w] for u, q in enumerate(row)]
         for (t, w), row in zip(states, J)]
    defect = rn.multiply(L, [[F(u == v)-D[u][v] for v in range(n)]
                            for u in range(n)])
    beta = [[sum(q for support, q in zip(supports, row) if v in support)-U[t][v]
             for v in range(n)] for (t, w), row in zip(states, B)]
    for i, (t, w) in enumerate(states):
        for v in range(n):
            assert beta[i][v] == F(v == w)-defect[i][v]
    return D, L, beta


def audit(K, p):
    data = contact_model(K, p)
    supports, states, M, T, reward, two_same, B, triples, identity_checks = data
    n = len(p)
    predicted = {}
    for t in range(n):
        initial = tuple(2*int(z == t) for z in range(n))
        predicted[initial] = two_same[t]
    for (t, w), row in zip(states, B):
        initial = [0]*n
        initial[t], initial[w] = 2, -1
        predicted[tuple(initial)] = row
    for roots, row in triples.items():
        initial = tuple(roots.count(z) for z in range(n))
        predicted[initial] = row
    laws, _ = rn.stopped(predicted, lambda s: rn.arw_step(K, p, s),
                         lambda s: all(a <= 0 for a in s))
    checked = 0
    for initial, row in predicted.items():
        assert sum(row) == 1 and all(q >= 0 for q in row)
        for C, q in zip(supports, row):
            final = tuple(-int(z in C) for z in range(n))
            assert q == laws[initial].get(final, F(0)), (initial, C)
            checked += 1
    forcing_bad = []
    sleeper_bad = []
    for v in range(n):
        baseline = [sum(q for C, q in zip(supports, two_same[t]) if v in C)
                    for t, w in states]
        for i, (t, w) in enumerate(states):
            forcing = sum(q for C, q in zip(supports, reward[i]) if v in C)
            forcing += sum(q*z for q, z in zip(T[i], baseline)) - baseline[i]
            if forcing < 0:
                forcing_bad.append((t, w, v, str(forcing)))
            gain = sum(q for C, q in zip(supports, B[i]) if v in C) - baseline[i]
            if gain < 0:
                sleeper_bad.append((t, w, v, str(gain)))
    creator_flow(K, p, data)
    return dict(n=n, contact_states=len(states), raw_law_entries=checked,
                creator_flow_entries=len(states)*n,
                trace_killing_identities=identity_checks,
                negative_forcing=forcing_bad[:1], negative_sleeper_gains=sleeper_bad[:1])


if __name__ == "__main__":
    from verify_backward_kernels import cases
    for name, K, p in cases():
        print(name, audit(K, p), flush=True)
