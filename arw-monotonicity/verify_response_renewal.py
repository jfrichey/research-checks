"""Exact audits for response-renewal.md; standard library only.

The selected, byte-preserved source files in evidence/base supply the original raw-ARW solver and signed
renewal construction. This script builds the new positive operator separately,
then compares both constructions on every formal response coordinate.
AI-written research code; finite checks do not replace the general proofs.
"""

from collections import defaultdict
from fractions import Fraction as F
from itertools import product
import json
from pathlib import Path
import importlib.util

BASE = Path(__file__).parent / "evidence/base"
spec = importlib.util.spec_from_file_location("archived_avoidance_audit", BASE / "verify_avoidance_exact.py")
old = importlib.util.module_from_spec(spec)
spec.loader.exec_module(old)
CASES = json.loads((BASE / "activation_blocks_exact_verification.json").read_text())["block_inverse_audits"]

stopped, arw_step = old.stopped, old.arw_step
inverse, multiply = old.inverse, old.multiply


def positive_operator(K, p, bound):
    """Return occupied marked states, T, killed-walk Z, and no-return W."""
    n = len(K)
    subsets = [tuple(i for i in range(n) if bits[i])
               for bits in product((0, 1), repeat=n) if 0 < sum(bits) <= bound]
    states = [(B, x) for B in subsets for x in B]
    index = {state: i for i, state in enumerate(states)}
    Z = {}
    for B in subsets:
        vacant = [z for z in range(n) if z not in B]
        matrix = [[F(z == w) - (1-p[z])*K[z][w] for w in vacant]
                  for z in vacant]
        rhs = [[(1-p[z])*K[z][t] for t in B] for z in vacant]
        hit = multiply(inverse(matrix), rhs)
        Z[B] = {z: {t: F(z == t) for t in B} for z in B}
        for i, z in enumerate(vacant):
            Z[B][z] = dict(zip(B, hit[i]))

    W = {}
    for x in range(n):
        starts = {}
        for B in subsets:
            if x not in B:
                continue
            starts[B] = defaultdict(F)
            for z, weight in enumerate(K[x]):
                if not weight:
                    continue
                state = [-int(t in B and t != x) for t in range(n)]
                if z < n:
                    state[z] = 2 if state[z] == -1 else 1
                starts[B][tuple(state)] += weight
        initial = {s for law in starts.values() for s in law}
        laws, _ = stopped(initial, lambda s: arw_step(K, p, s),
                          lambda s: s[x] > 0 or all(a <= 0 for a in s))
        for B, first_step in starts.items():
            W[B, x] = defaultdict(F)
            for state, weight in first_step.items():
                for final, probability in laws[state].items():
                    if final[x] > 0:
                        continue
                    C = tuple(z for z in range(n) if final[z] == -1)
                    W[B, x][C] += weight * probability

    T = [[F(0)] * len(states) for _ in states]
    for row, (B, x) in enumerate(states):
        for z in range(n):
            for t, probability in Z[B][z].items():
                T[row][index[B, t]] += K[x][z] * probability
        for C, weight in W[B, x].items():
            if C:
                for t, probability in Z[C][x].items():
                    T[row][index[C, t]] += weight * probability
    assert all(weight >= 0 for row in T for weight in row)
    assert all(not T[i][j] or len(C) <= len(B)
               for i, (B, x) in enumerate(states)
               for j, (C, y) in enumerate(states))
    return states, T, Z, W


def audit(case, bound):
    K = [[F(q) for q in row] for row in case["K"]]
    p = list(map(F, case["p"]))
    n = len(K)
    full_subsets = [tuple(i for i in range(n) if bits[i])
                    for bits in case["subsets"]]
    subsets = [B for B in full_subsets if len(B) <= bound]
    si = {B: i for i, B in enumerate(subsets)}
    selected = [full_subsets.index(B) for B in subsets]
    E = [[F(case["inverse_matrix"][i][j]) for j in selected] for i in selected]
    P = [[F(case["stabilization_matrix"][i][j]) for j in selected] for i in selected]
    assert multiply(E, P) == [[F(i == j) for j in selected] for i in selected]
    occupied, T, Z, W = positive_operator(K, p, bound)

    # Independent raw absorption for all stable-background additions.
    starts = {}
    for B in subsets:
        for x in range(n):
            state = [-int(z in B) for z in range(n)]
            state[x] = 2 if state[x] == -1 else 1
            starts[B, x] = tuple(state)
    laws, raw_count = stopped(starts.values(), lambda s: arw_step(K, p, s),
                             lambda s: all(a <= 0 for a in s))
    plus = {(B, x): {tuple(z for z in range(n) if c[z] == -1): q
                     for c, q in laws[s].items()} for (B, x), s in starts.items()}
    response = {(B, x): {C: plus[B, x].get(C, F(0))-F(C == B)
                         for C in full_subsets} for B in subsets for x in range(n)}
    for B in subsets:
        for x in range(n):
            for C in full_subsets:
                jump = sum(K[x][z]*response[B, z][C] for z in range(n))
                if x not in B:
                    create = F(C == tuple(sorted(B+(x,))))-F(C == B)
                    rhs = p[x]*create+(1-p[x])*jump
                else:
                    rhs = jump
                    for D, weight in W[B, x].items():
                        create = F(C == tuple(sorted(D+(x,))))-F(C == D)
                        rhs += weight*(response[D, x][C]-create)
                assert response[B, x][C] == rhs, (B, x, C)

    # Killed creation matrices and E N_x identity, entrywise.
    killed = {}
    for x in range(n):
        active = {B: tuple(int(z in B) for z in range(n)) for B in subsets}
        laws, _ = stopped(active.values(), lambda s: arw_step(K, p, s),
                          lambda s: s[x] > 0 or all(a <= 0 for a in s))
        killed[x] = [[F(0)] * len(subsets) for _ in subsets]
        for B, start in active.items():
            for final, weight in laws[start].items():
                if final[x] <= 0:
                    C = tuple(z for z in range(n) if final[z] == -1)
                    killed[x][si[B]][si[C]] += weight
        expected = [[F(B == C) if x not in B else
                     -(1-p[x])/p[x]*W[B, x].get(C, F(0))
                     for C in subsets] for B in subsets]
        assert multiply(E, killed[x]) == expected

    # Rebuild the original signed renewal and compare the new factorization
    # on every basis vector, without assuming D is a realizable gradient.
    old_edges, old_M, H, _ = old.exact_renewal(K, p, bound)
    old_edges = [(tuple(z for z in range(n) if bits[z]), x) for bits, x in old_edges]
    edges = [(B, x) for B in subsets for x in range(n) if x not in B]
    ei = {edge: i for i, edge in enumerate(edges)}
    M = [[F(0)] * len(edges) for _ in edges]
    for i, edge in enumerate(old_edges):
        for j, other in enumerate(old_edges):
            M[ei[edge]][ei[other]] = old_M[i][j]
    for col, basis in enumerate(edges):
        d = {}
        for B in subsets:
            vacant = [z for z in range(n) if z not in B]
            direct = {z: F((B, z) == basis) for z in vacant}
            if B:
                source_hit = multiply(
                    inverse([[F(z == t)-K[z][t] for t in B] for z in B]),
                    [[K[z][t] for t in vacant] for z in B])
                direct.update({z: sum(source_hit[i][j]*direct[t]
                                      for j, t in enumerate(vacant))
                               for i, z in enumerate(B)})
            for x in range(n):
                d[B, x] = direct[x]
        g = {x: [sum(E[i][j]*d[B, x] for j, B in enumerate(subsets))
                 for i, A in enumerate(subsets)] for x in range(n)}
        for row, (A, x) in enumerate(edges):
            if not A:
                assert M[row][col] == 0
                continue
            bits = tuple(int(z in A) for z in range(n))
            value = sum(H[bits][x][y] * sum(
                (P[si[A]][j]-killed[y][si[A]][j])*g[y][j]
                for j in range(len(subsets))) for y in range(n) if y not in A)
            assert value == M[row][col], (A, x, basis)

    contractions = []
    for k in range(1, bound+1):
        indices = [i for i, (B, x) in enumerate(occupied) if len(B) == k]
        h = [sum(q for C, q in plus[occupied[i]].items() if len(C) < k)
             for i in indices]
        b = []
        for i in indices:
            B, x = occupied[i]
            b.append(sum(weight*sum(q for D, q in plus[C, x].items() if len(D) < k)
                         for C, weight in W[B, x].items() if len(C) == k-1))
        assert all(0 < value <= 1 for value in h)
        assert all(value > 0 for value in b)
        assert all(h[r]-sum(T[i][j]*h[c] for c, j in enumerate(indices)) == b[r]
                   for r, i in enumerate(indices))
        contractions.append({"k": k, "minimum_defect": str(min(b)),
                             "contraction": str(max(1-bv/hv for bv, hv in zip(b, h)))})
    return {"n": n, "common_bound": bound, "raw_states": raw_count,
            "response_identities": len(subsets)*n*len(full_subsets),
            "killed_inverse_entries": n*len(subsets)**2,
            "renewal_factorization_entries": len(edges)**2,
            "maximum_T_row_sum": str(max(map(sum, T))),
            "strict_contractions": contractions, "verified": True}


def creator_capacity_obstruction():
    """CX: an unrestricted lower-only repair has insufficient total image mass."""
    K = [[F(0), F(1), F(0), F(0)],
         [F(1, 2), F(0), F(0), F(1, 2)],
         [F(0), F(1), F(0), F(0)]]
    p = [F(1, 2), F(99, 100), F(1, 2)]
    single, pair = (1, 0, 0), (1, 0, 1)
    laws, count = stopped([single, pair], lambda s: arw_step(K, p, s),
                          lambda s: all(a <= 0 for a in s))
    occupation = {start: sum(q for end, q in law.items() if end[0] == -1)
                  for start, law in laws.items()}
    green = inverse([[F(i == j)-(1-p[i])*K[i][j] for j in range(3)]
                     for i in range(3)])
    other_hit = green[2][0]/green[0][0]
    branch = (1-occupation[single])*(1-other_hit)
    capacity = 1-occupation[pair]
    assert occupation[single] == F(200, 399)
    assert other_hit == F(1, 400)
    assert branch == F(199, 400)
    assert capacity == F(99301, 238602)
    assert branch-capacity == F(3880699, 47720400) > 0
    return {"check": "unrestricted_lower_creator_repair_capacity",
            "K": [[str(q) for q in row] for row in K], "p": list(map(str, p)),
            "target": 0, "roots": [0, 2], "raw_states": count,
            "branch_mass": str(branch), "all_empty_target_images": str(capacity),
            "excess": str(branch-capacity), "verified": True,
            "scope": "No claim against the upper-conditioned repair or main conjecture"}


def last_arrival_deletion_audit():
    """Check the distributional identity and a failed conditional-background order."""
    K = [[F(0), F(1, 2), F(0), F(0), F(1, 2)],
         [F(0), F(0), F(1, 3), F(1, 3), F(1, 3)],
         [F(1, 3), F(0), F(0), F(1, 3), F(1, 3)],
         [F(1, 2), F(0), F(0), F(0), F(1, 2)]]
    p, v = [F(2, 3), F(1, 2), F(9, 10), F(9, 10)], 1
    lower, upper = (0, 1, 0, 1), (0, 1, 1, 1)
    starts = [lower, upper, (0, 1, 2, 0), (1, 0, 1, 0)]
    laws, count = stopped(starts, lambda s: arw_step(K, p, s),
                          lambda s: all(a <= 0 for a in s))
    killed, _ = stopped(starts, lambda s: arw_step(K, p, s),
                        lambda s: s[v] > 0 or all(a <= 0 for a in s))
    _, _, _, W = positive_operator(K, p, 3)
    distributions = {}
    entries = 0
    for start in starts:
        mu = {tuple(z for z, a in enumerate(final) if a == -1): weight
              for final, weight in laws[start].items()}
        no_visit = {tuple(z for z, a in enumerate(final) if a == -1): weight
                    for final, weight in killed[start].items() if final[v] <= 0}
        for bits in product((0, 1), repeat=4):
            C = tuple(z for z, a in enumerate(bits) if a)
            if v in C:
                continue
            value = ((1-p[v])/p[v])*sum(weight*W[B, v].get(C, F(0))
                       for B, weight in mu.items() if v in B)
            assert mu.get(C, F(0))-no_visit.get(C, F(0)) == value
            entries += 1
        distributions[start] = mu
    q = [sum(w for B, w in distributions[start].items() if v in B)
         for start in (lower, upper)]
    joint = [sum(w for B, w in distributions[start].items() if v in B and 3 in B)
             for start in (lower, upper)]
    conditional_difference = joint[1]/q[1]-joint[0]/q[0]
    assert q[1] > q[0] and joint[1] > joint[0]
    assert conditional_difference == F(-567216996528315050522131684863826450,
                                       4778182373573871099409781457610716184503)
    return {"check": "last_arrival_deletion_and_conditional_backgrounds",
            "K": [[str(q) for q in row] for row in K], "p": list(map(str, p)),
            "target": v, "lower": lower, "upper": upper,
            "raw_states": count, "deletion_identity_entries": entries,
            "target_occupations": list(map(str, q)), "joint_with_site_3": list(map(str, joint)),
            "conditional_difference": str(conditional_difference), "verified": True,
            "scope": "Conditional occupied-background domination fails; target increment positive"}


def visited_set_obstruction():
    """A monotone function of one completed visited set cannot represent q_v."""
    case = CASES[1]
    subsets = [tuple(i for i, bit in enumerate(bits) if bit) for bits in case['subsets']]
    P = [[F(q) for q in row] for row in case['stabilization_matrix']]
    assert case['p'] == ['1/2', '1/2']
    index = {A: i for i, A in enumerate(subsets)}
    R = [[F(0)]*4 for _ in range(4)]
    R[index[()]][index[()]] = 1
    R[index[(0, 1)]][index[(0, 1)]] = 1
    for site in (0, 1):
        R[index[(site,)]][index[(site,)]] = F(3, 4)
        R[index[(site,)]][index[(0, 1)]] = F(1, 4)
    occupation = [[sum(P[i][j] for j, B in enumerate(subsets) if 0 in B)]
                  for i in range(4)]
    phi = multiply(inverse(R), occupation)
    expected = {(): F(0), (0,): F(34, 65), (1,): F(-2, 195), (0, 1): F(22, 39)}
    assert all(phi[index[A]][0] == value for A, value in expected.items())
    return {'check': 'visited_set_representation_obstruction', 'target': 0,
            'phi': {str(A): str(value) for A, value in expected.items()},
            'verified': True, 'scope': 'No nonnegative function of the visited set alone'}


if __name__ == "__main__":
    for case_index in (1, 2, 3, 4):
        print(json.dumps(audit(CASES[case_index], min(len(CASES[case_index]["p"])-1, 2)),
                         sort_keys=True), flush=True)
    print(json.dumps(creator_capacity_obstruction(), sort_keys=True), flush=True)
    print(json.dumps(last_arrival_deletion_audit(), sort_keys=True), flush=True)
    print(json.dumps(visited_set_obstruction(), sort_keys=True), flush=True)
