#!/usr/bin/env python3
"""Bounded exact audits of RN Section 28; standard library only.

Two fixed paths and one outward fork audit weighted deletion identities,
Doob transforms, static cone bounds, and entrance return-count closure.
Every LP includes an independently checked rational primal/dual certificate.
Four-vertex path and fork feedback graphs are checked against every upper set.
Two fixed cyclic pairs audit static domination and portal-dependent closure.
A fixed directed three-cycle audits weighted intertwining and dominators.
Expected later-source, merge, and two-portal obstructions delimit the scope.
No parameter search, numerical approximation, or generated files are used.
AI-written research audit; the finite checks do not replace the proofs.
"""

from collections import defaultdict
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
from math import comb, prod
import json
import sys

sys.dont_write_bytecode = True
import verify_response_renewal as rn
from verify_creation_cone import simplex


def binary_states(n):
    return tuple(product((0, 1), repeat=n))


def bernoulli(bits, parameters):
    return prod((p if b else 1-p for b, p in zip(bits, parameters)), start=F(1))


def occupied(config):
    return tuple(int(x == -1) for x in config)


@lru_cache(None)
def upper_sets(n):
    """All upper sets, including the empty and full events."""
    states = binary_states(n)
    index = {b: k for k, b in enumerate(states)}
    covers = []
    for k, b in enumerate(states):
        for i in range(n):
            if not b[i]:
                c = list(b)
                c[i] = 1
                covers.append((k, index[tuple(c)]))
    return tuple(tuple(k for k in range(len(states)) if mask >> k & 1)
                 for mask in range(1 << len(states))
                 if all(not (mask >> k & 1) or (mask >> l & 1)
                        for k, l in covers))


def open_step(K, p, state):
    """One raw ARW instruction; the last two columns are portal and sink."""
    config, H = state
    n = len(p)
    i = next(i for i, x in enumerate(config) if x > 0)
    out = defaultdict(F)
    end = list(config)
    if end[i] == 1:
        end[i] = -1
    out[tuple(end), H] += p[i]
    for dest, weight in enumerate(K[i]):
        if not weight:
            continue
        end = list(config)
        end[i] -= 1
        if dest < n:
            end[dest] = 2 if end[dest] == -1 else end[dest]+1
        out[tuple(end), H+int(dest == n)] += (1-p[i])*weight
    return dict(out)


def open_absorption(K, p, starts, forbidden=None):
    def stop(state):
        config, _ = state
        return (all(x <= 0 for x in config)
                or (forbidden is not None and config[forbidden] > 0))

    raw, count = rn.stopped(list(dict.fromkeys(starts)),
                            lambda s: open_step(K, p, s), stop)
    laws = {}
    for start, row in raw.items():
        laws[start] = {(occupied(end), H): weight
                       for (end, H), weight in row.items()
                       if forbidden is None or end[forbidden] <= 0}
    return laws, count


def batch_law(b, source, M, p, routes):
    """Independent topological queue calculation, retaining exact exits."""
    states = {(0, 0, ()): F(1)}
    for i, (nxt, portal, sink) in enumerate(routes):
        out = defaultdict(F)
        for (incoming, H, zs), mass in states.items():
            incoming += M if i == source else 0
            if incoming == 0:
                out[0, H, zs+(b[i],)] += mass
                continue
            for z in (0, 1):
                L = incoming+b[i]-z
                for k in range(L+1):
                    for h in range(L-k+1):
                        weight = (comb(L, k)*comb(L-k, h)*nxt**k
                                  * portal**h*sink**(L-k-h))
                        if weight:
                            out[k, H+h, zs+(z,)] += mass*(p[i] if z else 1-p[i])*weight
        states = out
    assert all(k == 0 for k, H, z in states)
    law = {(z, H): mass for (k, H, z), mass in states.items() if mass}
    assert sum(law.values()) == 1
    return law


def deletion_laws(K, p, bits):
    """Force the marked first jump, then kill any return to that vertex."""
    n = len(p)
    W = {}
    for i in range(n):
        first = {}
        for b in bits:
            if not b[i]:
                continue
            row = defaultdict(F)
            for dest, weight in enumerate(K[i]):
                if not weight:
                    continue
                config = [-int(b[k] and k != i) for k in range(n)]
                if dest < n:
                    config[dest] = 2 if config[dest] == -1 else 1
                row[tuple(config), int(dest == n)] += weight
            first[b] = row
        laws, _ = open_absorption(K, p, [s for row in first.values() for s in row], i)
        for b, row in first.items():
            answer = defaultdict(F)
            for start, weight in row.items():
                for output, probability in laws[start].items():
                    answer[output] += weight*probability
            W[b, i] = dict(answer)
    return W


class ConeLP:
    def __init__(self, columns):
        self.columns = columns
        self.A = [[-a for a in col] for col in columns]+[[F(1)]*len(columns[0])]
        self.rhs = [F(0)]*len(columns)+[F(1)]
        self.checks = 0
        self.max_pivots = 0

    def nonpositive(self, objective, label):
        value, solution, pivots, dual = simplex(self.A, self.rhs, objective)
        assert all(x >= 0 for x in solution) and sum(solution) <= 1, label
        assert all(sum(a*x for a, x in zip(col, solution)) >= 0
                   for col in self.columns), label
        assert sum(c*x for c, x in zip(objective, solution)) == value, label
        assert all(x >= 0 for x in dual), label
        assert sum(y*b for y, b in zip(dual, self.rhs)) == value, label
        assert all(sum(row[k]*y for row, y in zip(self.A, dual)) >= objective[k]
                   for k in range(len(objective))), label
        # The zero measure is feasible, so the desired upper bound is exactly zero.
        assert value == 0, (label, value, solution)
        self.checks += 1
        self.max_pivots = max(self.max_pivots, pivots)


def path_audit(name, p, routes):
    n = len(p)
    bits = binary_states(n)
    index = {b: k for k, b in enumerate(bits)}
    sets = {b: tuple(i for i, x in enumerate(b) if x) for b in bits}
    K = [[F(0)]*(n+2) for _ in range(n)]
    h = [F(0)]*n
    for i in reversed(range(n)):
        nxt, portal, sink = routes[i]
        assert nxt+portal+sink == 1 and (i < n-1 or nxt == 0)
        if i < n-1:
            K[i][i+1] = nxt
        K[i][n], K[i][n+1] = portal, sink
        h[i] = sink+(nxt*h[i+1] if i < n-1 else 0)
        assert h[i] > 0
    beta = tuple(p[i]/(p[i]+(1-p[i])*h[i]) for i in range(n))
    mu_beta = [bernoulli(b, beta) for b in bits]
    hB = {b: prod((h[i] for i in sets[b]), start=F(1)) for b in bits}
    Khat = [[K[i][k]*h[k]/h[i] for k in range(n)]+[K[i][n+1]/h[i]]
            for i in range(n)]
    assert all(sum(row) == 1 for row in Khat)
    counts = defaultdict(int)

    inputs = [(0, M) for M in range(4)]+[(y, 1) for y in range(1, n)]
    starts = {}
    for source, M in inputs:
        for b in bits:
            config = [-x for x in b]
            if M:
                config[source] = M+b[source]
            starts[source, M, b] = (tuple(config), 0)
    raw, transient = open_absorption(K, p, starts.values())
    kernels = {}
    for (source, M, b), start in starts.items():
        row = raw[start]
        assert row == batch_law(b, source, M, p, routes), (name, source, M, b)
        kernels[source, M, b] = row
        counts['raw_vs_batch_entries'] += len(row)

    tilted, _ = rn.stopped(list(dict.fromkeys(s[0] for s in starts.values())),
                           lambda s: rn.arw_step(Khat, p, s),
                           lambda s: all(x <= 0 for x in s))
    for (source, M, b), start in starts.items():
        hat = {occupied(end): w for end, w in tilted[start[0]].items()}
        for c in bits:
            assert kernels[source, M, b].get((c, 0), F(0)) == (
                h[source]**M*hB[b]/hB[c]*hat.get(c, F(0)))
            counts['zero_addition_Doob_entries'] += 1

    W = deletion_laws(K, p, bits)
    _, _, _, What = rn.positive_operator(Khat, p, n)
    labels, columns = [], []
    for i in range(n):
        r = (1-p[i])/p[i]
        for c in bits:
            if c[i]:
                continue
            labels.append((i, c))
            column = []
            for b in bits:
                w0 = W.get((b, i), {}).get((c, 0), F(0))
                what = What.get((sets[b], i), {}).get(sets[c], F(0))
                assert w0 == hB[b]/hB[c]*what
                ell = F(b == c) if not b[i] else -r*w0
                ellhat = F(b == c) if not b[i] else -r*what
                assert ell == hB[b]/hB[c]*ellhat
                column.append(ell)
                counts['weighted_deletion_L0_entries'] += 1
            columns.append(column)
    assert all(sum(x*a for x, a in zip(mu_beta, col)) == 0 for col in columns)

    # Weighted last arrival (137), on all binary all-active creation inputs.
    active_starts = {eta: (eta, 0) for eta in bits}
    creation, _ = open_absorption(K, p, active_starts.values())
    for i in range(n):
        killed, _ = open_absorption(K, p, active_starts.values(), i)
        r = (1-p[i])/p[i]
        for eta, start in active_starts.items():
            row = creation[start]
            for H in range(sum(eta)+1):
                for c in bits:
                    if c[i]:
                        continue
                    correction = sum(row.get((b, H-k), F(0))
                                     *W[b, i].get((c, k), F(0))
                                     for b in bits if b[i] for k in range(H+1))
                    assert row.get((c, H), F(0)) == killed[start].get((c, H), F(0))+r*correction
                    counts['weighted_last_arrival_entries'] += 1
    normal = prod((p[i]+(1-p[i])*h[i] for i in range(n)), start=F(1))
    for b, w in zip(bits, mu_beta):
        assert creation[(tuple([1]*n), 0)].get((b, 0), F(0)) == normal*w
        counts['no_contact_product_entries'] += 1

    lp = ConeLP(columns)
    for i in range(n):
        for prefix in product((0, 1), repeat=i):
            objective = [(F(b[i])-beta[i]) if b[:i] == prefix else F(0) for b in bits]
            lp.nonpositive(objective, (name, 'prefix', i, prefix))
            counts['prefix_bound_duals'] += 1
    for event in upper_sets(n):
        if not event or len(event) == len(bits):
            continue
        benchmark = sum(mu_beta[k] for k in event)
        lp.nonpositive([F(k in event)-benchmark for k in range(len(bits))],
                       (name, 'static_upper_set', event))
        counts['static_upper_set_duals'] += 1
    for y in range(n):
        objective = [h[y]-sum(w for (c, H), w in kernels[y, 1, b].items() if H == 0)
                     for b in bits]
        lp.nonpositive(objective, (name, 'survival', y))
        counts['survival_duals'] += 1
    for M in range(4):
        for H in range(M+n+1):
            for (i, c), col in zip(labels, columns):
                residual = [sum(w*col[index[z]] for (z, hh), w in kernels[0, M, b].items()
                                if hh == H) for b in bits]
                if M and i == 0:
                    assert all(x >= 0 for x in residual), (name, 'first_constraint', M, H, c)
                    counts['point_mass_first_constraint_entries'] += len(residual)
                lp.nonpositive([-x for x in residual], (name, 'entrance_closure', M, H, i, c))
                counts['entrance_closure_duals'] += 1

    report = dict(name=name, p=p, routes_next_portal_sink=routes, h=h, beta=beta,
                  incoming_batches=list(range(4)), raw_transient_states=transient,
                  checks=dict(counts), checked_LP_duals=lp.checks, max_pivots=lp.max_pivots)
    data = dict(name=name, p=p, routes=routes, K=K, h=h, beta=beta, bits=bits, index=index,
                mu_beta=mu_beta, normal=normal, labels=labels, columns=columns, kernels=kernels)
    return report, data


def scope_audit(data):
    """A reachable entrance density failure, and a later-source cone failure."""
    bits, mu = data['bits'], data['mu_beta']
    kernels, index, columns = data['kernels'], data['index'], data['columns']

    def advance(law, source, M, H):
        return [sum(law[k]*kernels[source, M, b].get((c, H), F(0))
                    for k, b in enumerate(bits)) for c in bits]

    reachable = advance(advance(mu, 0, 1, 3), 0, 3, 1)
    probability = sum(reachable)
    normalized = [x/probability for x in reachable]
    assert all(sum(x*a for x, a in zip(normalized, col)) >= 0 for col in columns)
    increase = normalized[index[(1, 0)]]/mu[index[(1, 0)]]-normalized[index[(0, 0)]]/mu[index[(0, 0)]]
    assert increase == F(27047483, 39998901)

    later = advance(mu, 1, 1, 2)
    mass = sum(later)
    conditional = [x/mass for x in later]
    column = columns[data['labels'].index((0, (0, 1))) ]
    defect = sum(x*a for x, a in zip(conditional, column))
    assert mass == F(99, 39800)
    assert conditional == [F(1,199), F(0), F(198,199), F(0)]
    assert defect == -F(891, 99500)
    return dict(entrance_history=[3, 1], entrance_history_in_cone=True,
                entrance_history_density_increase=increase,
                later_source=1, later_source_batch=1, later_source_returns=2,
                later_source_return_probability=mass, later_source_conditional_law=conditional,
                later_source_residual_label=[0, [0, 1]], later_source_residual=defect,
                scope='later-source obstruction; not an entrance-path counterexample')


def fork_batch_law(b, M, p, K):
    """Independent root split and leaf queues for the fixed outward fork."""
    assert M >= 1 and K[0][0] == K[0][3] == K[0][4] == 0
    assert K[0][1]+K[0][2] == 1
    assert all(K[i][k] == 0 for i in (1, 2) for k in range(3))
    answer = defaultdict(F)
    for z0 in (0, 1):
        L = M+b[0]-z0
        for k in range(L+1):
            split = comb(L, k)*K[0][1]**k*K[0][2]**(L-k)
            states = {((z0,), 0): (p[0] if z0 else 1-p[0])*split}
            for i, incoming in ((1, k), (2, L-k)):
                out = defaultdict(F)
                for (zs, H), mass in states.items():
                    if incoming == 0:
                        out[zs+(b[i],), H] += mass
                        continue
                    for z in (0, 1):
                        emitted = incoming+b[i]-z
                        for h in range(emitted+1):
                            weight = comb(emitted, h)*K[i][3]**h*K[i][4]**(emitted-h)
                            out[zs+(z,), H+h] += mass*(p[i] if z else 1-p[i])*weight
                states = out
            for output, mass in states.items():
                answer[output] += mass
    law = {output: mass for output, mass in answer.items() if mass}
    assert sum(law.values()) == 1
    return law


def fork_audit():
    """One fixed heterogeneous fork; full-cone closure for root batches 1,2."""
    name = 'heterogeneous_outward_fork'
    n = 3
    p = (F(1,3), F(2,5), F(3,7))
    # Columns: root, leaf 1, leaf 2, portal, sink.
    K = [[F(0), F(1,3), F(2,3), F(0), F(0)],
         [F(0), F(0), F(0), F(1,4), F(3,4)],
         [F(0), F(0), F(0), F(3,5), F(2,5)]]
    assert all(sum(row) == 1 for row in K)
    bits = binary_states(n)
    index = {b: k for k, b in enumerate(bits)}
    sets = {b: tuple(i for i, x in enumerate(b) if x) for b in bits}
    h = [F(0)]*n
    for i in reversed(range(n)):
        h[i] = K[i][n+1]+sum(K[i][k]*h[k] for k in range(i+1, n))
    assert h == [F(31,60), F(3,4), F(2,5)]
    beta = tuple(p[i]/(p[i]+(1-p[i])*h[i]) for i in range(n))
    mu_beta = [bernoulli(b, beta) for b in bits]
    hB = {b: prod((h[i] for i in sets[b]), start=F(1)) for b in bits}
    Khat = [[K[i][k]*h[k]/h[i] for k in range(n)]+[K[i][n+1]/h[i]]
            for i in range(n)]
    assert all(sum(row) == 1 for row in Khat)
    W = deletion_laws(K, p, bits)
    _, _, _, What = rn.positive_operator(Khat, p, n)
    labels, columns = [], []
    counts = defaultdict(int)
    for i in range(n):
        r = (1-p[i])/p[i]
        for c in bits:
            if c[i]:
                continue
            labels.append((i, c))
            column = []
            for b in bits:
                w0 = W.get((b, i), {}).get((c, 0), F(0))
                what = What.get((sets[b], i), {}).get(sets[c], F(0))
                assert w0 == hB[b]/hB[c]*what
                ell = F(b == c) if not b[i] else -r*w0
                ellhat = F(b == c) if not b[i] else -r*what
                assert ell == hB[b]/hB[c]*ellhat
                column.append(ell)
                counts['weighted_deletion_L0_entries'] += 1
            columns.append(column)
    assert all(sum(x*a for x, a in zip(mu_beta, col)) == 0 for col in columns)

    starts = {(M, b): ((M+b[0], -b[1], -b[2]), 0)
              for M in (1, 2) for b in bits}
    raw, transient = open_absorption(K, p, starts.values())
    kernels = {}
    for (M, b), start in starts.items():
        row = raw[start]
        assert row == fork_batch_law(b, M, p, K), (name, M, b)
        assert all(H <= M+sum(b) for c, H in row)
        kernels[M, b] = row
        counts['raw_vs_batch_entries'] += len(row)
    normal = prod((p[i]+(1-p[i])*h[i] for i in range(n)), start=F(1))
    active = ((1,)*n, 0)
    creation, _ = open_absorption(K, p, [active])
    for b, w in zip(bits, mu_beta):
        assert creation[active].get((b, 0), F(0)) == normal*w
        counts['no_contact_product_entries'] += 1

    lp = ConeLP(columns)
    for i in range(n):
        for prefix in product((0, 1), repeat=i):
            objective = [(F(b[i])-beta[i]) if b[:i] == prefix else F(0) for b in bits]
            lp.nonpositive(objective, (name, 'prefix', i, prefix))
            counts['prefix_bound_duals'] += 1
    for event in upper_sets(n):
        if not event or len(event) == len(bits):
            continue
        benchmark = sum(mu_beta[k] for k in event)
        lp.nonpositive([F(k in event)-benchmark for k in range(len(bits))],
                       (name, 'static_upper_set', event))
        counts['static_upper_set_duals'] += 1
    for M in (1, 2):
        for H in range(M+n+1):
            for (i, c), col in zip(labels, columns):
                residual = [sum(w*col[index[z]] for (z, hh), w in kernels[M, b].items()
                                if hh == H) for b in bits]
                if i == 0:
                    assert all(x >= 0 for x in residual), (name, 'root_constraint', M, H, c)
                    counts['point_mass_root_constraint_entries'] += len(residual)
                lp.nonpositive([-x for x in residual], (name, 'entrance_closure', M, H, i, c))
                counts['entrance_closure_duals'] += 1
    report = dict(name=name, p=p, K=K, boundary_columns=['portal', 'sink'], h=h, beta=beta,
                  incoming_batches=[1, 2], raw_transient_states=transient, checks=dict(counts),
                  checked_LP_duals=lp.checks, max_pivots=lp.max_pivots,
                  scope='one fixed outward fork; no general-forest inference from this finite audit')
    data = dict(name=name, p=p, K=K, h=h, beta=beta, bits=bits, index=index,
                mu_beta=mu_beta, normal=normal, labels=labels, columns=columns)
    return report, data


def diamond_scope_audit():
    """Fixed diamond: product-beta entrance closure fails even after mixing H."""
    n = 4
    p = (F(1,2),)*n
    K = [[F(0)]*(n+2) for _ in range(n)]
    K[0][1] = K[0][2] = F(1,2)
    K[1][3] = K[2][3] = F(1)
    K[3][n] = K[3][n+1] = F(1,2)
    bits = binary_states(n)
    mu = {b: bernoulli(b, (F(2,3),)*n) for b in bits}
    W = deletion_laws(K, p, bits)
    checked = 0
    for i in range(n):
        for c in bits:
            if c[i]:
                continue
            residual = mu[c]-sum(mu[b]*W[b, i].get((c, 0), F(0))
                                 for b in bits if b[i])
            assert residual == 0
            checked += 1
    starts = {b: ((1+b[0],)+tuple(-x for x in b[1:]), 0) for b in bits}
    raw, transient = open_absorption(K, p, starts.values())
    assert all(sum(row.values()) == 1 for row in raw.values())
    # Four initial sleepers plus one entrant make these counts exhaustive.
    return_counts = tuple(range(n+2))
    assert all(H in return_counts for row in raw.values() for c, H in row)
    slices = {H: {c: sum(mu[b]*raw[starts[b]].get((c, H), F(0)) for b in bits)
                  for c in bits} for H in return_counts}
    nu = slices[3]
    C, B1, B2 = (1, 0, 0, 1), (1, 1, 0, 0), (1, 0, 1, 0)
    assert nu[C] == 0 and nu[B1] == nu[B2] == F(1,648)
    assert W[B1, 1].get((C, 0), F(0)) == F(1,2)
    contributions = {b: nu[b]*W[b, 1].get((C, 0), F(0))
                     for b in bits if b[1] and nu[b]*W[b, 1].get((C, 0), F(0))}
    assert contributions == {B1: F(1,1296)}
    defect = nu[C]-sum(contributions.values())
    assert defect == -F(1,1296) and sum(nu.values()) > 0
    residuals = tuple(slices[H][C]-sum(slices[H][b]*W[b, 1].get((C, 0), F(0))
                                      for b in bits if b[1]) for H in return_counts)
    assert residuals == (F(0), F(0), -F(1,432), -F(1,1296), F(0), F(0))
    unconditional = {c: sum(slices[H][c] for H in return_counts) for c in bits}
    assert sum(unconditional.values()) == 1
    unconditional_defect = unconditional[C]-sum(
        unconditional[b]*W[b, 1].get((C, 0), F(0)) for b in bits if b[1])
    assert unconditional_defect == sum(residuals) == -F(1,324)
    return dict(name='merging_diamond_scope_obstruction', p=p, K=K,
                initial_beta=[F(2,3)]*n, initial_zero_residual_columns=checked,
                source=0, incoming_batch=1, portal_returns=3, raw_transient_states=transient,
                return_probability=sum(nu.values()), output_C=list(C), output_C_mass=nu[C],
                predecessor_B1=list(B1), predecessor_B2=list(B2), predecessor_masses=[nu[B1], nu[B2]],
                deletion_vertex=1, deletion_weight=F(1,2), unnormalized_residual=defect,
                return_counts=return_counts, return_count_residuals=residuals,
                all_return_counts_residual=unconditional_defect,
                scope='entrance/product-law obstruction with a merge; not an outward-fork counterexample')


def cyclic_two_site_audit():
    """Fixed cyclic pairs: static pi_0 bounds and two distinct portal scopes."""
    n = 2
    p = (F(1,3), F(2,5))
    bits = binary_states(n)
    index = {b: k for k, b in enumerate(bits)}
    sets = {b: tuple(i for i, x in enumerate(b) if x) for b in bits}
    support_order = ((0, 0), (1, 0), (0, 1), (1, 1))
    cases = (
        ('two_site_cyclic_both_portals',
         [[F(0), F(1,2), F(1,3), F(1,6)],
          [F(1,3), F(0), F(1,6), F(1,2)]],
         (F(1,2), F(2,3)), (1,)),
        ('two_site_cyclic_partner_only_portal',
         [[F(0), F(1,2), F(0), F(1,2)],
          [F(1,3), F(0), F(1,6), F(1,2)]],
         (F(9,10), F(4,5)), (1, 2)),
    )
    reports = []
    for name, K, expected_h, batches in cases:
        assert all(sum(row) == 1 for row in K)
        d = 1-K[0][1]*K[1][0]
        h = tuple((K[i][n+1]+K[i][1-i]*K[1-i][n+1])/d for i in range(n))
        assert h == expected_h
        assert all(h[i] == K[i][n+1]+sum(K[i][k]*h[k] for k in range(n))
                   for i in range(n))
        hB = {b: prod((h[i] for i in sets[b]), start=F(1)) for b in bits}
        Khat = [[K[i][k]*h[k]/h[i] for k in range(n)]+[K[i][n+1]/h[i]]
                for i in range(n)]
        assert all(sum(row) == 1 for row in Khat)
        W = deletion_laws(K, p, bits)
        _, _, _, What = rn.positive_operator(Khat, p, n)
        labels, columns = [], []
        counts = defaultdict(int)
        for i in range(n):
            r = (1-p[i])/p[i]
            for c in bits:
                if c[i]:
                    continue
                labels.append((i, c))
                column = []
                for b in bits:
                    w0 = W.get((b, i), {}).get((c, 0), F(0))
                    what = What.get((sets[b], i), {}).get(sets[c], F(0))
                    assert w0 == hB[b]/hB[c]*what
                    ell = F(b == c) if not b[i] else -r*w0
                    ellhat = F(b == c) if not b[i] else -r*what
                    assert ell == hB[b]/hB[c]*ellhat
                    column.append(ell)
                    counts['weighted_deletion_L0_entries'] += 1
                columns.append(column)
        assert len(columns) == 4

        starts = {(M, b): ((M+b[0], -b[1]), 0) for M in batches for b in bits}
        active = ((1, 1), 0)
        raw, transient = open_absorption(K, p, list(starts.values())+[active])
        assert all(sum(row.values()) == 1 and all(w >= 0 for w in row.values())
                   for row in raw.values())
        kernels = {(M, b): raw[start] for (M, b), start in starts.items()}
        assert all(0 <= H and sum(c)+H <= M+sum(b)
                   for (M, b), row in kernels.items() for c, H in row)
        normal = sum(w for (c, H), w in raw[active].items() if H == 0)
        assert normal > 0
        pi0 = [raw[active].get((b, 0), F(0))/normal for b in bits]
        assert sum(pi0) == 1 and all(w > 0 for w in pi0)
        initial_residuals = [sum(w*a for w, a in zip(pi0, col)) for col in columns]
        assert initial_residuals == [F(0)]*len(columns)
        lp = ConeLP(columns)
        for event in upper_sets(n):
            if not event or len(event) == len(bits):
                continue
            benchmark = sum(pi0[k] for k in event)
            lp.nonpositive([F(k in event)-benchmark for k in range(len(bits))],
                           (name, 'static_pi0_upper_set', event))
            counts['static_upper_set_duals'] += 1
        assert counts['static_upper_set_duals'] == 4
        report = dict(name=name, p=p, K=K, boundary_columns=['portal', 'sink'], h=h,
                      bits_order=bits, pi0=pi0, no_portal_full_root_mass=normal,
                      initial_cone_residuals=initial_residuals,
                      source=0, incoming_batches=batches, raw_transient_states=transient)
        if K[0][n]:
            assert [pi0[index[b]] for b in support_order] == [F(x,150) for x in (31,39,38,42)]
            nu = [sum(pi0[k]*kernels[1, b].get((c, 2), F(0))
                      for k, b in enumerate(bits)) for c in bits]
            mass = sum(nu)
            assert mass == F(428,4725)
            conditional = [w/mass for w in nu]
            assert [conditional[index[b]] for b in support_order] == [F(x,428) for x in (260,21,147,0)]
            column = columns[labels.index((1, (1, 0)))]
            defect = sum(w*a for w, a in zip(nu, column))
            normalized_defect = sum(w*a for w, a in zip(conditional, column))
            assert defect == -F(1,1350) and normalized_defect == -F(7,856)
            report.update(portal_returns=2, return_probability=mass,
                          conditional_law=conditional, residual_label=[1, [1, 0]],
                          unnormalized_residual=defect, normalized_residual=normalized_defect,
                          scope='static pi0 domination holds; both-portal count closure fails')
        else:
            for M in batches:
                for H in range(M+n+1):
                    for (i, c), col in zip(labels, columns):
                        residual = [sum(w*col[index[z]] for (z, hh), w in kernels[M, b].items()
                                        if hh == H) for b in bits]
                        lp.nonpositive([-x for x in residual], (name, 'entrance_closure', M, H, i, c))
                        counts['entrance_closure_duals'] += 1
            assert counts['entrance_closure_duals'] == 36
            report.update(return_counts_by_batch={M: list(range(M+n+1)) for M in batches},
                          cone_columns=len(columns),
                          scope='fixed partner-only portal; source-0 batches 1,2, every return count')
        report.update(checks=dict(counts), checked_LP_duals=lp.checks, max_pivots=lp.max_pivots)
        reports.append(report)
    return reports


def cyclic_dominator_audit():
    """One fixed cycle; exact count polynomials, with no symbolic dependency."""
    n = 3
    p = (F(1,3), F(2,5), F(3,7))
    K = [[F(0), F(1,2), F(0), F(0), F(1,2)],
         [F(0), F(0), F(2,3), F(0), F(1,3)],
         [F(1,3), F(0), F(0), F(1,3), F(1,3)]]
    assert all(sum(row) == 1 for row in K)
    bits = binary_states(n)
    counts = defaultdict(int)

    # Sparse entries (degree, row, column) retain exact portal counts.
    def law_poly(rows, cols, laws):
        index = {c: k for k, c in enumerate(cols)}
        return {(H, k, index[c]): w for k, b in enumerate(rows)
                for (c, H), w in laws[b].items() if w}

    def add(left, right, scale=F(1)):
        answer = dict(left)
        for key, w in right.items():
            answer[key] = answer.get(key, F(0))+scale*w
        return {key: w for key, w in answer.items() if w}

    def mul(left, right):
        by_row = defaultdict(list)
        for (H, row, col), w in right.items():
            by_row[row].append((H, col, w))
        answer = defaultdict(F)
        for (H, row, middle), w in left.items():
            for k, col, v in by_row[middle]:
                answer[H+k, row, col] += w*v
        return {key: w for key, w in answer.items() if w}

    def coefficient(poly, H):
        return {key: w for key, w in poly.items() if key[0] == H}

    def check(left, right, rows, cols, degree, label):
        assert all(0 <= H <= degree and 0 <= row < rows and 0 <= col < cols
                   for poly in (left, right) for H, row, col in poly), label
        assert left == right, (label, add(left, right, -F(1)))
        counts[label] += rows*cols*(degree+1)

    def matrix_at_one(poly, rows, cols):
        answer = [[F(0)]*cols for _ in range(rows)]
        for (H, row, col), w in poly.items():
            answer[row][col] += w
        return answer

    starts = {}
    for y, M in [(y, 1) for y in range(n)]+[(0, 2)]:
        for b in bits:
            config = [-x for x in b]
            config[y] = M+b[y]
            starts[y, M, b] = (tuple(config), 0)
    active_starts = {b: (b, 0) for b in bits}
    raw, transient = open_absorption(K, p, list(starts.values())+list(active_starts.values()))
    assert all(sum(row.values()) == 1 and all(w >= 0 for w in row.values())
               for row in raw.values())
    assert all(0 <= H and sum(c)+H <= M+sum(b)
               for (y, M, b), start in starts.items() for c, H in raw[start])
    A = {y: law_poly(bits, bits, {b: raw[starts[y, 1, b]] for b in bits}) for y in range(n)}
    batch = {1: A[0], 2: law_poly(bits, bits, {b: raw[starts[0, 2, b]] for b in bits})}
    check(mul(A[0], A[0]), batch[2], len(bits), len(bits), 2*(n+1),
          'sequential_vs_batch_coefficients')

    # P(1) invertibility checks the finite-matrix input to the general argument.
    P = law_poly(bits, bits, {b: raw[active_starts[b]] for b in bits})
    P1 = matrix_at_one(P, len(bits), len(bits))
    inverse = rn.inverse(P1)
    identity = [[F(i == j) for j in range(len(bits))] for i in range(len(bits))]
    assert rn.multiply(P1, inverse) == rn.multiply(inverse, P1) == identity
    counts['binary_creation_inverse_entries'] = 2*len(bits)**2

    W = deletion_laws(K, p, bits)
    labels, columns = [], []
    for i in range(n):
        fresh = tuple(b for b in bits if not b[i])
        fresh_index = {b: k for k, b in enumerate(fresh)}
        r = (1-p[i])/p[i]
        # E keeps rows avoiding i; D removes i only on rows containing i.
        E, D = {}, {}
        for row, b in enumerate(bits):
            c = tuple(x if k != i else 0 for k, x in enumerate(b))
            (D if b[i] else E)[0, row, fresh_index[c]] = F(1)
        Wpoly = law_poly(bits, fresh, {b: W.get((b, i), {}) for b in bits})
        L = add(E, Wpoly, -r)
        L0 = coefficient(L, 0)
        for col, c in enumerate(fresh):
            labels.append((i, c))
            columns.append([L0.get((0, row, col), F(0)) for row in range(len(bits))])

        killed_starts = {start for (y, M, b), start in starts.items() if M == 1 and not b[i]}
        killed, _ = open_absorption(K, p, killed_starts | set(active_starts.values()), i)
        Ahat = {y: law_poly(fresh, fresh, {b: killed[starts[y, 1, b]] for b in fresh})
                for y in range(n)}
        assert Ahat[i] == {}
        no_visit = law_poly(bits, fresh, {b: killed[active_starts[b]] for b in bits})
        check(mul(P, L), no_visit, len(bits), len(fresh), 2*n,
              'weighted_last_arrival_coefficients')
        for y in range(n):
            check(mul(A[y], L), mul(L, Ahat[y]), len(bits), len(fresh), 2*n+1,
                  'weighted_intertwining_coefficients')

        J = {(H, k, k): w for H, w in ((0, K[i][n+1]), (1, K[i][n]))
             for k in range(len(fresh)) if w}
        for y in range(n):
            J = add(J, Ahat[y], K[i][y])
        check(mul(D, J), Wpoly, len(bits), len(fresh), n,
              'forced_jump_factorization_coefficients')

        # A graph check and an independent raw check of the dominator condition.
        reached = set() if i == 0 else {0}
        todo = list(reached)
        while todo:
            y = todo.pop()
            assert K[y][n] == 0, ('entry_to_portal_bypasses_mark', i, y)
            for dest in range(n):
                if dest != i and K[y][dest] and dest not in reached:
                    reached.add(dest)
                    todo.append(dest)
        counts['entry_portal_dominators'] += 1
        B = Ahat[0]
        assert all(H == 0 for H, row, col in B)
        counts['killed_entry_zero_positive_count_coefficients'] += n*len(fresh)**2
        bracket = add(mul(A[0], D), mul(D, B), -F(1))
        assert all(w >= 0 for w in bracket.values()), ('dominator_bracket', i)
        counts['pointwise_nonnegative_bracket_coefficients'] += (n+2)*len(bits)*len(fresh)
        delta_J = add(J, coefficient(J, 0), -F(1))
        assert all(w >= 0 and H > 0 for (H, row, col), w in delta_J.items())
        right = add(mul(L0, B), mul(bracket, delta_J), r)
        check(mul(A[0], L0), right, len(bits), len(fresh), 2*n+1,
              'dominator_identity_coefficients')

    lp = ConeLP(columns)
    for M in (1, 2):
        for H in range(M+n+1):
            for (i, c), col in zip(labels, columns):
                residual = [sum(batch[M].get((H, row, k), F(0))*col[k]
                                for k in range(len(bits))) for row in range(len(bits))]
                lp.nonpositive([-w for w in residual], ('cyclic_dominator_closure', M, H, i, c))
                counts['entrance_closure_duals'] += 1
    assert len(columns) == 12 and counts['entrance_closure_duals'] == 132
    return dict(name='directed_three_cycle_dominator', p=p, K=K,
                boundary_columns=['portal', 'sink'], source=0, portal_predecessor=2,
                incoming_batches=[1, 2], return_counts_by_batch={M: list(range(M+n+1)) for M in (1, 2)},
                raw_transient_states=transient, cone_columns=len(columns), checks=dict(counts),
                checked_LP_duals=lp.checks, max_pivots=lp.max_pivots,
                scope='one fixed cycle; all source/mark intertwining identities and entry-0 closure only')


def source_sharing_audit():
    """A three-site cone law cannot be peeled as a convex two-root mixture."""
    n = 3
    p = (F(1,2),)*n
    K = [[F(0) if i == j else F(1,4) for j in range(n)]
         +[F(0), F(1,2)] for i in range(n)]
    bits = binary_states(n)
    W = deletion_laws(K, p, bits)
    by_size = {0: F(1161,11571), 1: F(1518,11571),
               2: F(1952,11571), 3: F(0)}
    mu = {b: by_size[sum(b)] for b in bits}
    assert sum(mu.values()) == 1
    checked = 0
    for i in range(n):
        for c in bits:
            if not c[i]:
                residual = mu[c]-sum(mu[b]*W[b,i].get((c,0), F(0))
                                     for b in bits if b[i])
                assert residual == (F(16,133) if sum(c) == 2 else F(0))
                checked += 1
    roots = [b for b in bits if sum(b) == 2]
    raw, count = open_absorption(K, p, [(b,0) for b in roots])
    pair_matrix = [[raw[(b,0)].get((c,0), F(0)) for c in roots] for b in roots]
    inverse = rn.inverse(pair_matrix)
    coefficients = rn.multiply([[mu[c] for c in roots]], inverse)[0]
    assert coefficients == [F(61,133)]*3 and sum(coefficients) == F(183,133) > 1
    assert F(133,16) < 3*F(61,16)
    return dict(name='three_site_shared_source_mass', p=p, K=K,
                mass_per_support_by_size=by_size, checked_deletion_constraints=checked,
                raw_transient_states=count, pair_source_residual=F(16,133),
                forced_two_root_mixture_coefficients=coefficients,
                required_mixture_mass=sum(coefficients),
                scope='obstruction to convex root-law peeling; not to static domination')


def physical_audit(domain, entrance_probability=F(4,5)):
    """Close a three-vertex open domain through its root; portal numbered 0."""
    n = len(domain['p'])+1
    assert n == 4
    assert 0 <= entrance_probability <= 1
    p = (F(4,9),)+tuple(domain['p'])
    K = [[F(0)]*(n+1) for _ in range(n)]
    K[0][1], K[0][n] = entrance_probability, 1-entrance_probability
    for i, row in enumerate(domain['K'], start=1):
        for dest in range(n-1):
            K[i][dest+1] = row[dest]
        K[i][0], K[i][n] = row[n-1], row[n]
    assert all(sum(row) == 1 for row in K)
    full, missing, mixed = (1,)*n, (0,)+(1,)*(n-1), (-1,)+(1,)*(n-1)
    continuation = {b: (1,)+tuple(-x for x in b) for b in domain['bits']}
    starts = [full, missing, mixed]+list(continuation.values())
    raw, transient = rn.stopped(starts, lambda s: rn.arw_step(K, p, s),
                                lambda s: all(x <= 0 for x in s))
    bits = binary_states(n)
    def row(start):
        law = {occupied(end): w for end, w in raw[start].items()}
        assert sum(law.values()) == 1 and all(w >= 0 for w in law.values())
        return [law.get(b, F(0)) for b in bits]
    upper, lower, asleep = row(full), row(missing), row(mixed)
    continuation_rows = {b: row(start) for b, start in continuation.items()}
    U = [sum(w*continuation_rows[b][k] for b, w in zip(domain['bits'], domain['mu_beta']))
         for k in range(len(bits))]
    benchmark = [F(0)]*len(bits)
    index = {b: k for k, b in enumerate(bits)}
    for b, w in zip(domain['bits'], domain['mu_beta']):
        benchmark[index[(1,)+b]] = w
    r = (1-p[0])/p[0]
    sigma = K[0][n]+K[0][1]*domain['h'][0]
    assert sigma > 0
    for k in range(len(bits)):
        assert asleep[k]-upper[k] == r*sigma*(upper[k]-lower[k])
        assert domain['normal']*(benchmark[k]-U[k]) == r*sigma*(upper[k]-lower[k])
    events = upper_sets(n)
    assert len(events) == 168
    gains = []
    for event in events:
        gain = sum(upper[k]-lower[k] for k in event)
        budget = sum(benchmark[k]-U[k] for k in event)
        assert gain >= 0 and budget >= 0, ('physical_upper_set', event, gain, budget)
        if event and len(event) < len(bits):
            gains.append(gain)
    root_law = [sum(U[k] for k, b in enumerate(bits) if b[1:] == c) for c in domain['bits']]
    assert all(sum(x*a for x, a in zip(root_law, col)) >= 0 for col in domain['columns'])
    return dict(name=domain['name']+'_feedback', p=p, K=K, portal=0, raw_initial_states=len(starts),
                raw_transient_states=transient, first_contact_vector_entries=len(bits),
                generic_payoff_singleton_entries=len(bits), upper_set_comparisons=len(events),
                continuation_upper_set_budgets=len(events), final_root_cone_columns=len(domain['columns']),
                minimum_nontrivial_upper_set_gain=min(gains), sigma_portal=sigma,
                scope='full versus full-minus-designated-portal; all occupied-set upper events')


def main():
    cases = (
        ('sharp_two_buffer', (F(99,100),)*2,
         ((F(9,10), F(1,20), F(1,20)), (F(0), F(1,2), F(1,2)))),
        ('heterogeneous_three_buffer', (F(1,3), F(2,5), F(3,7)),
         ((F(1,2), F(1,3), F(1,6)), (F(2,5), F(1,5), F(2,5)),
          (F(0), F(3,5), F(2,5)))),
    )
    data = []
    for name, p, routes in cases:
        report, case_data = path_audit(name, p, routes)
        data.append(case_data)
        print(json.dumps(report, default=str), flush=True)
    print(json.dumps(scope_audit(data[0]), default=str), flush=True)
    print(json.dumps(physical_audit(data[1]), default=str), flush=True)
    fork_report, fork_data = fork_audit()
    print(json.dumps(fork_report, default=str), flush=True)
    print(json.dumps(physical_audit(fork_data, entrance_probability=F(1)), default=str), flush=True)
    print(json.dumps(diamond_scope_audit(), default=str), flush=True)
    for report in cyclic_two_site_audit():
        print(json.dumps(report, default=str), flush=True)
    print(json.dumps(cyclic_dominator_audit(), default=str), flush=True)
    print(json.dumps(source_sharing_audit(), default=str), flush=True)
    print('PASS: fixed path/fork/cyclic audits; later-source, merge, and two-portal obstructions verified.', flush=True)


if __name__ == '__main__':
    main()
