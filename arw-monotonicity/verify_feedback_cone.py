#!/usr/bin/env python3
"""Exact finite evidence for the complete-feedback cone question in RN28.

A single open-domain pass can leave C^0. Two fixed examples preserve it
after the whole excursion through an added portal vertex; a cyclic
domain, a later-entry DAG, and a source-entry DAG give exact counterexamples.
Standard library only.
Binary laws use the proved complete-block inverse and surplus-IDLA algebra;
archived matrices and raw absorption audit it. The cyclic counterexample
also checks both full four-site matrices with reversed toppling priority
and certifies that its reference law maximizes the bad-column violation.
The source-entry example has all sleep probabilities 1/2 and an independent
112-state complete-return renewal matching every feedback entry.
AI-written research verifier; no parameter search or generated files.
"""

from collections import defaultdict
from fractions import Fraction as F
from functools import lru_cache
from math import comb
import json
import sys

sys.dont_write_bytecode = True
import verify_response_renewal as rn
from verify_creation_cone import simplex
from verify_creator_potentials import idla_addition_matrices
from verify_backward_queue_cone import (binary_states, occupied, open_absorption,
                                       deletion_laws, upper_sets, ConeLP, bernoulli)


def binary_stabilization(K, p):
    """E(B,.) = delta_empty product_{i in B}(D_i+r_i(D_i-J_i))."""
    n = len(p)
    supports = tuple(frozenset(i for i, x in enumerate(b) if x)
                     for b in binary_states(n))
    m = len(supports)
    D = idla_addition_matrices(K, supports)
    sleep_insertions = []
    for x in range(n):
        r = (1-p[x])/p[x]
        sleep_insertions.append([
            [(1+r)*D[x][i][j]-r*(sum(K[x][y]*D[y][i][j] for y in range(n))
                                 +K[x][n]*F(i == j)) for j in range(m)]
            for i in range(m)])
    E = []
    for B in supports:
        row = [F(not C) for C in supports]
        for x in sorted(B):
            row = rn.multiply([row], sleep_insertions[x])[0]
        E.append(row)
    P = rn.inverse(E)
    assert all(sum(row) == 1 and min(row) >= 0 for row in P)
    assert rn.multiply(E, P) == [[F(i == j) for j in range(m)] for i in range(m)]
    return supports, P, E, D


def archived_audit():
    case = next(c for c in rn.CASES if len(c['K']) == 4)
    K = [[F(x) for x in row] for row in case['K']]
    p = tuple(map(F, case['p']))
    supports, P, E, _ = binary_stabilization(K, p)
    assert [tuple(int(i in B) for i in range(4)) for B in supports] == [tuple(b) for b in case['subsets']]
    assert P == [[F(x) for x in row] for row in case['stabilization_matrix']]
    assert E == [[F(x) for x in row] for row in case['inverse_matrix']]
    return dict(archived_stabilization_entries=256, archived_inverse_entries=256)


def raw_matrix_audit(K, p, supports, P, E, D, A):
    """Keep the portal bit and independently topple the highest active site."""
    n, m = len(p), len(supports)
    def step(state):
        x = max(i for i, count in enumerate(state) if count > 0)
        out = defaultdict(F)
        end = list(state)
        if end[x] == 1:
            end[x] = -1
        out[tuple(end)] += p[x]
        for y, weight in enumerate(K[x]):
            if not weight:
                continue
            end = list(state)
            end[x] -= 1
            if y < n:
                end[y] = 2 if end[y] == -1 else end[y]+1
            out[tuple(end)] += (1-p[x])*weight
        assert sum(out.values()) == 1
        return dict(out)

    active = {B: tuple(int(i in B) for i in range(n)) for B in supports}
    # Physical portal is zero; this includes backgrounds already containing it.
    addition = {B: (1+int(0 in B),)+tuple(-int(i in B) for i in range(1, n))
                for B in supports}
    terminal = {B: tuple(-int(i in B) for i in range(n)) for B in supports}
    raw, states = rn.stopped(list(active.values())+list(addition.values()), step,
                             lambda s: all(x <= 0 for x in s))
    assert all(sum(row.values()) == 1 and min(row.values()) >= 0 for row in raw.values())
    raw_P = [[raw[active[B]].get(terminal[C], F(0)) for C in supports] for B in supports]
    raw_A = [[raw[addition[B]].get(terminal[C], F(0)) for C in supports] for B in supports]
    assert raw_P == P and raw_A == A
    raw_E = rn.inverse(raw_P)
    assert raw_E == E
    assert rn.multiply(rn.multiply(raw_E, D), raw_P) == raw_A
    return dict(raw_binary_stabilization_entries=m*m, raw_binary_inverse_entries=m*m,
                raw_active_addition_entries=m*m, raw_addition_conjugacy_entries=m*m,
                raw_matrix_transient_states=states, raw_toppling_priority='highest active vertex')


def order_interval_obstruction(bits, pi0, feedback):
    """Deleting one initial sleeper need not preserve the reference budget."""
    assert bits == binary_states(2)
    pos = {b: i for i, b in enumerate(bits)}
    nu = list(pi0)
    nu[pos[(0,1)]] += nu[pos[(1,1)]]
    nu[pos[(1,1)]] = F(0)
    assert nu == [F(31,150), F(80,150), F(39,150), F(0)]
    for event in upper_sets(2):
        assert sum(nu[k] for k in event) <= sum(pi0[k] for k in event)
    target = [k for k, b in enumerate(bits) if b[1]]
    benchmark = sum(pi0[k] for k in target)
    value = sum(nu[i]*sum(feedback[i][k] for k in target) for i in range(len(bits)))
    assert benchmark == F(8,15)
    assert value-benchmark == F(294469,31635375) > 0
    return dict(initial_law=nu, dominated_by_reference=True, target=1,
                reference_probability=benchmark, feedback_probability=value,
                excess=value-benchmark,
                scope='ordinary stochastic domination is insufficient even for complete feedback')


def feedback_audit(name, K, p, portal_p, portal_routes, raw_check=False,
                   expected_failure=None, full_raw_check=False):
    """K has open vertices followed by portal and sink columns."""
    n = len(p)
    assert len(portal_routes) == n+1 and sum(portal_routes) == 1
    assert all(len(row) == n+2 and sum(row) == 1 for row in K)
    bits = binary_states(n)
    pos = {b: i for i, b in enumerate(bits)}
    W = deletion_laws(K, p, bits)
    columns, labels = [], []
    for i in range(n):
        r = (1-p[i])/p[i]
        for c in bits:
            if not c[i]:
                labels.append((i, c))
                columns.append([F(b == c) if not b[i]
                                else -r*W[b, i].get((c, 0), F(0)) for b in bits])
    full_active = ((1,)*n, 0)
    open_law, open_states = open_absorption(K, p, [full_active])
    z = sum(w for (c, H), w in open_law[full_active].items() if H == 0)
    assert z > 0
    pi0 = [open_law[full_active].get((b, 0), F(0))/z for b in bits]
    assert sum(pi0) == 1
    assert all(sum(x*a for x, a in zip(pi0, col)) == 0 for col in columns)

    # The new physical portal is vertex zero.
    full_K = [[F(0)]+list(portal_routes)]
    for row in K:
        full_K.append([row[n]]+list(row[:n])+[row[n+1]])
    full_p = (portal_p,)+tuple(p)
    supports, P, E, D = binary_stabilization(full_K, full_p)
    A = rn.multiply(rn.multiply(E, D[0]), P)
    assert all(sum(row) == 1 and min(row) >= 0 for row in A)
    full_raw = raw_matrix_audit(full_K, full_p, supports, P, E, D[0], A) if full_raw_check else None
    support_pos = {B: i for i, B in enumerate(supports)}
    feedback = []
    for b in bits:
        B = frozenset(i+1 for i, x in enumerate(b) if x)
        row = [F(0)]*len(bits)
        for C, w in zip(supports, A[support_pos[B]]):
            row[pos[tuple(int(i+1 in C) for i in range(n))]] += w
        assert min(row) >= 0 and sum(row) == 1
        feedback.append(row)
    raw_entries = 0
    if raw_check:
        starts = {b: (1,)+tuple(-x for x in b) for b in bits}
        raw, _ = rn.stopped(list(starts.values()), lambda s: rn.arw_step(full_K, full_p, s),
                            lambda s: all(x <= 0 for x in s))
        for b, row in zip(bits, feedback):
            actual = [F(0)]*len(bits)
            for end, w in raw[starts[b]].items():
                actual[pos[occupied(end)[1:]]] += w
            assert row == actual
            raw_entries += len(bits)

    # Independent first-contact vector identity for this continuation.
    G = rn.inverse([[F(i == j)-K[i][j] for j in range(n)] for i in range(n)])
    h = [sum(G[i][j]*K[j][n+1] for j in range(n)) for i in range(n)]
    sigma = portal_routes[n]+sum(portal_routes[i]*h[i] for i in range(n))
    r_portal = (1-portal_p)/portal_p
    V = frozenset(range(n+1))
    upper, lower = P[support_pos[V]], P[support_pos[V-{0}]]
    delta = [F(0)]*len(bits)
    for B, a, b in zip(supports, upper, lower):
        delta[pos[tuple(int(i+1 in B) for i in range(n))]] += a-b
    terminal = [sum(pi0[k]*feedback[k][l] for k in range(len(bits))) for l in range(len(bits))]
    assert terminal == [pi0[l]-r_portal*sigma/z*delta[l] for l in range(len(bits))]
    terminal_residuals = [sum(w*a for w, a in zip(terminal, col)) for col in columns]
    failure = None
    if expected_failure is None:
        assert min(terminal_residuals) >= 0
    else:
        label, expected_value = expected_failure
        assert [(lab, d) for lab, d in zip(labels, terminal_residuals) if d < 0] == [expected_failure]
        assert terminal_residuals[labels.index(label)] == expected_value < 0
        failure = dict(label=label, residual=expected_value, initial_reference=pi0)

    lp = ConeLP(columns)
    if expected_failure is None:
        for label, col in zip(labels, columns):
            residual = [sum(w*a for w, a in zip(row, col)) for row in feedback]
            lp.nonpositive([-x for x in residual], (name, 'complete_feedback', label))
    else:
        # Maximize the violation on C^0 with mass <= 1. Certify both LP sides.
        col = columns[labels.index(expected_failure[0])]
        objective = [-sum(w*a for w, a in zip(row, col)) for row in feedback]
        value, solution, pivots, dual = simplex(lp.A, lp.rhs, objective)
        assert min(solution) >= 0 and sum(solution) <= 1
        assert all(sum(w*a for w, a in zip(solution, column)) >= 0 for column in columns)
        assert sum(w*c for w, c in zip(solution, objective)) == value
        assert min(dual) >= 0 and sum(y*b for y, b in zip(dual, lp.rhs)) == value
        assert all(sum(row[k]*y for row, y in zip(lp.A, dual)) >= objective[k]
                   for k in range(len(bits)))
        assert solution == pi0 and sum(solution) == 1 and value == -expected_failure[1] > 0
        lp.checks += 1
        lp.max_pivots = max(lp.max_pivots, pivots)
        failure.update(maximum_bad_column_violation=value, lp_maximizer_is_reference=True,
                       lp_primal_dual_verified=True)
    reference_upper_events = 0
    if n <= 3:
        for event in upper_sets(n):
            assert sum(pi0[k]-terminal[k] for k in event) >= 0
            reference_upper_events += 1
    interval_failure = (order_interval_obstruction(bits, pi0, feedback)
                        if name == 'two_site_both_portals' else None)
    return dict(name=name, open_K=K, open_p=p, portal_p=portal_p,
                portal_routes_then_sink=portal_routes, open_raw_transient_states=open_states,
                feedback_matrix_entries=len(bits)**2, raw_feedback_entries=raw_entries,
                full_raw_matrix_audit=full_raw,
                first_contact_entries=len(bits), checked_LP_duals=lp.checks,
                max_pivots=lp.max_pivots, zero_reference_terminal_residuals=terminal_residuals.count(0),
                minimum_reference_terminal_residual=min(terminal_residuals),
                reference_upper_event_budgets=reference_upper_events, cone_counterexample=failure,
                order_interval_obstruction=interval_failure,
                scope=('fixed counterexample to cone preservation; reference occupation budgets still hold'
                       if failure else 'fixed complete-feedback cone preservation; no general theorem'))


def source_entry_queue_law(b, incoming, epsilon):
    """Exact topological pass on the fixed DAG d,0,1,2,3; retain portal count.

    Routes are d -> 0/sink, 0 -> 1/2, 1,2 -> 3, and 3 -> portal/sink.
    At a visited vertex, the final bit is a fresh Ber(1/2), and the
    number of departures is incoming + initial sleeper - final bit.
    Unvisited vertices retain their initial bit. No binary inverse is used.
    """
    assert len(b) == len(incoming) == 5
    states = {(tuple(incoming), (), 0): F(1)}
    for i in range(5):
        following = defaultdict(F)
        for (counts, zs, H), mass in states.items():
            if counts[i] == 0:
                following[counts, zs+(b[i],), H] += mass
                continue
            for z in (0, 1):
                departures = counts[i]+b[i]-z
                for k in range(departures+1):
                    if i == 0:
                        weight = (F(comb(departures, k))*epsilon**k
                                  *(1-epsilon)**(departures-k))
                    elif i in (1, 4):
                        weight = F(comb(departures, k), 2**departures)
                    else:
                        if k != departures:
                            continue
                        weight = F(1)
                    if not weight:
                        continue
                    remaining = list(counts)
                    remaining[i] = 0
                    exits = H
                    if i == 0:
                        remaining[1] += k
                    elif i == 1:
                        remaining[2] += k
                        remaining[3] += departures-k
                    elif i in (2, 3):
                        remaining[4] += departures
                    else:
                        exits += k
                    following[tuple(remaining), zs+(z,), exits] += mass*weight/2
        states = dict(following)
    law = defaultdict(F)
    for (counts, zs, H), mass in states.items():
        assert not any(counts)
        assert sum(zs)+H <= sum(b)+sum(incoming)
        law[zs, H] += mass
    assert min(law.values()) >= 0 and sum(law.values()) == 1
    return dict(law)


def source_entry_return_renewal(bits, epsilon):
    """Stabilize j, then the DAG, retaining every return to j.

    A state (b,M) has stable DAG bits b and M active particles at j.
    Its last particle leaves bit z ~ Ber(1/2); M-z particles enter d.
    If H>0 particles return, the next active count is H+z, including
    the awakened resident at j. If H=0, the whole system is stable.
    """
    @lru_cache(None)
    def pass_law(b, incoming):
        return source_entry_queue_law(b, (incoming, 0, 0, 0, 0), epsilon)

    def step(state):
        b, M = state
        assert M >= 1 and sum(b)+M <= 6
        out = defaultdict(F)
        for z in (0, 1):
            incoming = M-z
            if incoming == 0:
                out[b, 0] += F(1, 2)
                continue
            for (c, H), weight in pass_law(b, incoming).items():
                out[c, H+z if H else 0] += weight/2
        assert sum(out.values()) == 1
        return dict(out)

    starts = [(b, 1) for b in bits]
    laws, transient = rn.stopped(starts, step, lambda state: state[1] == 0)
    assert transient == 112
    rows = []
    for b in bits:
        assert all(M == 0 for c, M in laws[b, 1])
        row = [laws[b, 1].get((c, 0), F(0)) for c in bits]
        assert min(row) >= 0 and sum(row) == 1
        rows.append(row)
    return rows, transient


def source_entry_dag_audit():
    """One fixed source-entry cone obstruction; no LP or parameter search.

    The complete residual also has expansion -epsilon/1296+O(epsilon^2).
    At epsilon=0 the initial source bit and its last bit are independent
    Ber(1/2), so its first visited batch has mean one departure. For H>=2,
    the portal emits at least H-1>=1 particles back to d, resetting its
    final bit to Ber(1/2). The initial portal jump supplies the other 1/2.
    H=0,1 have zero open-diamond residual separately. Histories with no
    diamond entrance have zero residual for every initial source law;
    its epsilon-dependence therefore contributes no linear term.
    """
    epsilon = F(1, 4860)
    vertex_order = ('d', '0', '1', '2', '3')
    n = len(vertex_order)
    p = (F(1, 2),)*n
    portal_p = F(1, 2)
    portal_routes = (F(1), F(0), F(0), F(0), F(0), F(0))
    K = [[F(0)]*(n+2) for _ in range(n)]
    K[0][1], K[0][n+1] = epsilon, 1-epsilon
    K[1][2] = K[1][3] = F(1, 2)
    K[2][4] = K[3][4] = F(1)
    K[4][n] = K[4][n+1] = F(1, 2)
    assert all(sum(row) == 1 for row in K)
    assert all(not K[i][k] for i in range(n) for k in range(i+1))
    sources = [k for k in range(n) if all(not K[i][k] for i in range(n))]
    assert sources == [0]
    assert portal_routes[0] == 1 and not any(portal_routes[1:])

    bits = binary_states(n)
    pos = {b: i for i, b in enumerate(bits)}
    h = (1-epsilon/2,)+(F(1, 2),)*4
    assert all(h[i] == K[i][n+1]+sum(K[i][k]*h[k] for k in range(n))
               for i in range(n))
    beta = tuple(p[i]/(p[i]+(1-p[i])*h[i]) for i in range(n))
    assert beta == (F(9720, 19439),)+(F(2, 3),)*4
    pi0 = [bernoulli(b, beta) for b in bits]
    # Direct topological stabilization independently verifies the reference.
    root_law = source_entry_queue_law((0,)*n, (1,)*n, epsilon)
    zero_exit_mass = sum(w for (c, H), w in root_law.items() if H == 0)
    assert zero_exit_mass > 0
    assert [root_law.get((b, 0), F(0))/zero_exit_mass for b in bits] == pi0

    target = (0, 1, 0, 0, 1)
    pair = (0, 1, 1, 0, 0)
    triple = (0, 1, 1, 0, 1)
    phi = [F(b == target)-F(1, 2)*F(b == pair)-F(1, 4)*F(b == triple)
           for b in bits]
    # Deleting site labelled 1 (index 2) forces a jump to site 3 (index 4).
    # No return to 1 is possible in the DAG, and zero portal count gives W^0.
    for b, coefficient in zip(bits, phi):
        if b[2]:
            deleted = b[:2]+(0,)+b[3:]
            forced = source_entry_queue_law(deleted, (0, 0, 0, 0, 1), epsilon)
            assert coefficient == -forced.get((target, 0), F(0))
        else:
            assert coefficient == F(b == target)
    assert sum(w*a for w, a in zip(pi0, phi)) == 0

    full_K = [[F(0)]+list(portal_routes)]
    full_K += [[row[n]]+row[:n]+[row[n+1]] for row in K]
    full_p = (portal_p,)+p
    supports, P, E, D = binary_stabilization(full_K, full_p)
    A = rn.multiply(rn.multiply(E, D[0]), P)
    assert all(min(row) >= 0 and sum(row) == 1 for row in A)
    support_pos = {B: i for i, B in enumerate(supports)}
    feedback = []
    for b in bits:
        B = frozenset(i+1 for i, x in enumerate(b) if x)
        row = [F(0)]*len(bits)
        for C, w in zip(supports, A[support_pos[B]]):
            row[pos[tuple(int(i+1 in C) for i in range(n))]] += w
        assert min(row) >= 0 and sum(row) == 1
        feedback.append(row)
    renewal, transient = source_entry_return_renewal(bits, epsilon)
    assert feedback == renewal  # All 1,024 entries, not just the bad payoff.
    terminal = [sum(pi0[k]*feedback[k][l] for k in range(len(bits)))
                for l in range(len(bits))]
    expected = -F(
        815854470961087421045964392289833963226600657811975473605569834502612804828894789226329312052772014047225336192224724158763238259124960,
        5142631108291643266970096601611875478561455703618858290681854746260042104086137688579870306010848084895731462984290765898413531497344645974347)
    residual = sum(w*a for w, a in zip(terminal, phi))
    assert residual == expected < 0

    V = frozenset(range(n+1))
    delta = [F(0)]*len(bits)
    for B, upper, lower in zip(supports, P[support_pos[V]], P[support_pos[V-{0}]]):
        delta[pos[tuple(int(i+1 in B) for i in range(n))]] += upper-lower
    assert terminal == [pi0[k]-h[0]/zero_exit_mass*delta[k] for k in range(len(bits))]

    # Exact coefficients behind the analytic first-order calculation.
    by_count = [F(0)]*6
    for diamond_bits in binary_states(4):
        weight = bernoulli(diamond_bits, (F(2, 3),)*4)
        law = source_entry_queue_law((0,)+diamond_bits, (0, 1, 0, 0, 0), epsilon)
        for (c, H), probability in law.items():
            by_count[H] += weight*probability*phi[pos[c]]
    assert by_count == [F(0), F(0), -F(1, 432), -F(1, 1296), F(0), F(0)]
    mean_departures = sum(F(1, 4)*(1+b-z) for b in (0, 1) for z in (0, 1))
    assert mean_departures == 1
    first_order_weights = []
    for H in range(6):
        untouched = F(1) if H == 0 else F(1, 2) if H == 1 else F(0)
        weight = sum(F(1, 8)*(1+b-z)
                     *(untouched*(1-z)+(1-untouched)/2)
                     for b in (0, 1) for z in (0, 1))
        first_order_weights.append(weight)
    assert first_order_weights == [F(3, 8), F(5, 16)]+[F(1, 4)]*4
    linear = sum(w*a for w, a in zip(first_order_weights, by_count))
    assert linear == -F(1, 1296)

    # The weaker density question survives this same fixed cone obstruction.
    density_checks = 0
    for b in bits:
        for i in range(n):
            if not b[i]:
                c = b[:i]+(1,)+b[i+1:]
                assert terminal[pos[c]]/pi0[pos[c]] <= terminal[pos[b]]/pi0[pos[b]]
                density_checks += 1
    assert density_checks == 80
    return dict(name='source_entry_DAG_obstruction', epsilon=epsilon,
                open_vertex_order=vertex_order, open_K=K, open_p=p,
                portal_p=portal_p, portal_routes_then_sink=portal_routes,
                full_vertex_order=('j',)+vertex_order, full_K=full_K, full_p=full_p,
                dag_topological_order=vertex_order, unique_source='d',
                portal_jumps_only_to_source=True, h=h, reference_beta=beta,
                support_bit_order=bits, reference_law=pi0,
                full_root_zero_exit_mass=zero_exit_mass,
                cone_column=dict(deleted_vertex='1', target=('0', '3'), coefficients=phi),
                reference_cone_residual=F(0), exact_terminal_cone_residual=residual,
                complete_return_transient_states=transient,
                independently_matched_feedback_entries=len(bits)**2,
                first_contact_entries=len(bits),
                open_diamond_residuals_by_exit_count=by_count,
                source_mean_departures_at_epsilon_zero=mean_departures,
                first_order_weights_by_exit_count=first_order_weights,
                residual_linear_coefficient=linear,
                reference_density_cover_checks=density_checks,
                scope='source-entry DAG cone preservation fails; the fixed reference density decreases')


def main():
    print(json.dumps(archived_audit()), flush=True)
    K = [[F(0), F(1,2), F(1,3), F(1,6)],
         [F(1,3), F(0), F(1,6), F(1,2)]]
    print(json.dumps(feedback_audit('two_site_both_portals', K, (F(1,3), F(2,5)),
                                   F(1,2), (F(1), F(0), F(0)), True), default=str), flush=True)
    K = [[F(0)]*6 for _ in range(4)]
    K[0][1] = K[0][2] = F(1,2)
    K[1][3] = K[2][3] = F(1)
    K[3][4] = K[3][5] = F(1,2)
    print(json.dumps(feedback_audit('merging_diamond', K, (F(1,2),)*4, F(1,2),
                                   (F(1), F(0), F(0), F(0), F(0))), default=str), flush=True)
    K = [[F(0), F(0), F(1,2), F(1,4), F(1,4)],
         [F(1,2), F(0), F(0), F(1,4), F(1,4)],
         [F(15,16), F(3,64), F(0), F(1,128), F(1,128)]]
    failure = ((1, (1,0,1)), -F(429949477180782443029835986988779712426332577408224,
                               539009793644576310168618312663274330682808973998650717))
    print(json.dumps(feedback_audit('cyclic_three_site_obstruction', K, (F(1,2),)*3,
                                   F(1,2), (F(0), F(0), F(1), F(0)), True, failure,
                                   full_raw_check=True), default=str), flush=True)
    K = [[F(0), F(1), F(0), F(0)],
         [F(0), F(0), F(1,2), F(1,2)]]
    print(json.dumps(feedback_audit('later_entry_DAG_obstruction', K, (F(1,2),)*2,
                                   F(1,2), (F(0), F(1), F(0)), True,
                                   ((0, (0,1)), -F(2,105))), default=str), flush=True)
    print(json.dumps(source_entry_dag_audit(), default=str), flush=True)
    print('PASS: 36 feedback preservation certificates, two optimal counterexamples, '
          'and a source-entry DAG obstruction; raw and complete-return checks agree.', flush=True)


if __name__ == '__main__':
    main()
