#!/usr/bin/env python3
"""Exact all-rate checks of four-site graphs and selected five-site motifs.

Enumerate simple sink graphs with four nonsink vertices and connected
internal graph, up to relabeling that fixes the sink. Disconnected internal
graphs factor into the already proved smaller cases. For each law, verify
all 174 disjoint increasing-event inequalities as univariate polynomials
in t=1/lambda. Negative coefficients alone are not a counterexample.

The optional five-site cases test whether coefficientwise positivity of
all NA-gap numerators survives branching, a cycle below a cut vertex, or
two conditional-correlation motifs sharing a vertex.
This stronger algebraic property could bridge the proved small- and
large-rate regimes. These are selected diagnostics, not a five-site census.

--endpoint-recursion checks two exact obstructions to closing endpoint
leaf recursion: rate-only conditional deletion and preservation under
active insertion of externally tilted trace laws. It checks only the
specified two-/three-site systems, including symbolic raw absorption;
it does not enumerate graphs or refute stationary endpoint-path NA.

--endpoint-transfer checks a sink-end transfer retaining the special
stationary traced input. It verifies the two-to-three-site transfer at
all common rates, raw absorption, and the covariance/Rayleigh identities.
The rank-pair obstruction concerns a termwise proof, not endpoint NA.
"""

import argparse
import itertools
import json
from pathlib import Path

from algebra import stationary, trace_kernel
from exact import (
    DomainMatrix, graph_kernel, increasing_events, sp, stable_mask,
    stationary_law, transitions,
)


def graph_types():
    n = 4
    edges = tuple(itertools.combinations(range(n+1), 2))
    edge_index = {edge: j for j, edge in enumerate(edges)}
    permutations = tuple(itertools.permutations(range(n)))
    seen = set()
    for bits in range(1 << len(edges)):
        chosen = [edge for j, edge in enumerate(edges) if bits >> j & 1]
        if not any(n in edge for edge in chosen):
            continue
        reachable = {0}
        for _ in range(n):
            reachable |= {b for a, b in chosen if a in reachable and b < n}
            reachable |= {a for a, b in chosen if b in reachable and b < n}
        if len(reachable) != n:
            continue
        labels = []
        for perm in permutations:
            pp = perm + (n,)
            labels.append(sum(1 << edge_index[tuple(sorted((pp[a], pp[b])))]
                              for a, b in chosen))
        canonical = min(labels)
        if canonical in seen:
            continue
        seen.add(canonical)
        representative = [edge for j, edge in enumerate(edges) if canonical >> j & 1]
        neighbors = tuple(tuple(b if a == v else a for a, b in representative if v in (a, b))
                          for v in range(n))
        yield canonical, neighbors


def polynomial_certificate(mu, variable):
    n = len(mu).bit_length()-1
    variables = tuple(variable) if isinstance(variable, (tuple, list)) else (variable,)
    fractions = [tuple(sp.Poly(v, *variables) for v in sp.fraction(sp.cancel(value)))
                 for value in mu]
    denominator = sp.Poly(1, *variables)
    for _, den in fractions:
        denominator = sp.lcm(denominator, den)
    if denominator.LC() < 0:
        denominator = -denominator
    denominator_positive = all(c >= 0 for c in denominator.coeffs())
    weights = [num*denominator.exquo(den) for num, den in fractions]
    zero = sp.Poly(0, *variables)
    assert sum(weights, zero) == denominator
    comparisons, failures, distinct = 0, [], set()
    for am in range(1, 1 << n):
        aa = tuple(i for i in range(n) if am >> i & 1)
        for bm in range(am+1, 1 << n):
            if am & bm:
                continue
            bb = tuple(i for i in range(n) if bm >> i & 1)
            for ea in increasing_events(aa):
                ia = [sum(((state >> i) & 1) << j for j, i in enumerate(aa)) in ea
                      for state in range(1 << n)]
                pa = sum((w for w, keep in zip(weights, ia) if keep), zero)
                for eb in increasing_events(bb):
                    ib = [sum(((state >> i) & 1) << j for j, i in enumerate(bb)) in eb
                          for state in range(1 << n)]
                    pb = sum((w for w, keep in zip(weights, ib) if keep), zero)
                    joint = sum((w for w, x, y in zip(weights, ia, ib) if x and y), zero)
                    gap = pa*pb-denominator*joint
                    comparisons += 1
                    distinct.add(tuple(gap.terms()))
                    if any(c < 0 for c in gap.coeffs()):
                        failures.append({"coordinates": [aa, bb], "events": [ea, eb],
                                         "gap_numerator": str(sp.factor(gap.as_expr()))})
    return {"comparisons": comparisons, "distinct_gap_polynomials": len(distinct),
            "denominator_positive_coefficients": denominator_positive,
            "all_rates_NA_certified": denominator_positive and not failures,
            "failed_coefficient_tests": failures,
            "denominator": str(sp.factor(denominator.as_expr())),
            "numerators": [str(sp.factor(w.as_expr())) for w in weights]}


def rayleigh_polynomial_certificate(case, variable):
    """All positive fields, by coefficient positivity of each Rayleigh gap.

    Multiply univariate parameter polynomials while keeping the remaining
    field exponents as small integer tuples. This avoids repeatedly
    expanding a large symbolic expression in all field variables.
    """
    weights = [sp.Poly(sp.sympify(value), variable) for value in case["numerators"]]
    n = len(weights).bit_length()-1
    zero = sp.Poly(0, variable)
    pairs = []
    for i, j in itertools.combinations(range(n), 2):
        remaining = [k for k in range(n) if k not in (i, j)]
        parts = [[], [], [], []]
        for mask, weight in enumerate(weights):
            which = ((mask >> i) & 1)+2*((mask >> j) & 1)
            parts[which].append((tuple((mask >> k) & 1 for k in remaining), weight))
        coefficients = {}
        for left, right, sign in ((1, 2, 1), (0, 3, -1)):
            for powers, ww in parts[left]:
                for other_powers, vv in parts[right]:
                    exponent = tuple(a+b for a, b in zip(powers, other_powers))
                    coefficients[exponent] = coefficients.get(exponent, zero)+sign*ww*vv
        failures = [{"field_exponents": list(powers), "coefficient_polynomial": str(value.as_expr())}
                    for powers, value in coefficients.items() if any(c < 0 for c in value.coeffs())]
        pairs.append({"pair": [i, j], "remaining_fields": remaining,
                      "field_coefficients_checked": len(coefficients),
                      "failed_coefficient_tests": failures})
    return {"all_rates_Rayleigh_certified": case["denominator_positive_coefficients"]
            and all(not item["failed_coefficient_tests"] for item in pairs), "pairs": pairs}


def _endpoint_symbolic_absorption(starts, kernel, sleep, reverse=False):
    """Raw absorption over a rational-function field, for transient starts.

    Keep this independent of the IDLA change of basis. Check every transient
    equation and row normalization before returning the requested start laws
    as columns. The endpoint certificate uses only fixed two-/three-site
    chains, with the specified legal toppling orders.
    """
    states = list(starts)
    index = {state: i for i, state in enumerate(states)}
    assert len(index) == len(states) and states
    assert all(any(value > 0 for value in state) for state in states)
    rows = []
    for state in states:
        row = transitions(state, kernel, sleep, reverse=reverse)
        rows.append(row)
        for nxt in row:
            if any(value > 0 for value in nxt) and nxt not in index:
                index[nxt] = len(states)
                states.append(nxt)
    left = sp.eye(len(states))
    right = sp.zeros(len(states), 1 << len(kernel))
    for i, row in enumerate(rows):
        for nxt, probability in row.items():
            if nxt in index:
                left[i, index[nxt]] -= probability
            else:
                right[i, stable_mask(nxt)] += probability
    aa = DomainMatrix.from_Matrix(left).to_field()
    bb = DomainMatrix.from_Matrix(right).convert_to(aa.domain)
    solution = aa.inv().matmul(bb)
    assert (aa.matmul(solution)-bb).is_zero_matrix
    matrix = solution.to_Matrix()
    assert all(sp.cancel(sum(matrix.row(i))-1) == 0 for i in range(matrix.rows))
    return matrix[:len(starts), :].T, len(states)


def endpoint_recursion_certificate():
    """Exact conditional-deletion and tilted-insertion obstructions.

    Let mu_m^(c) have rate lambda on the endpoint path except for terminal
    rate c*lambda. Deleting the terminal site leaves a holding probability
    1/2 at its predecessor; erasing holds gives terminal rate 2*lambda.
    The report's active-insertion identity therefore reads
    mu_m^(c) = I_m^(c) mu_(m-1)^(2). In I_m^(c), the OLD terminal's rate
    changes from 2*lambda in the erased trace representation to lambda in
    the restored enlarged path. This is not appending a leaf with all
    existing rates frozen. Intermediate growth uses c=2; the final
    common-rate path uses c=1.

    Rate rigidity holds at every size: conditional on an empty far leaf,
    the remaining generating polynomial has a zero with coordinates
    z_i=(k-i)/(k+1) for i<k, z_k=-x_k/(2*(k+1)), and z_i=0 for i>k,
    for each k through the new terminal. Here x_i is the original inverse
    sleep rate. These are the report's harmonic evaluations, not stability.
    If the conditional law were stationary on the smaller simple path
    with any rates, its corresponding roots force inverse rates x_i at
    interior sites and x_k/2 at the new terminal. Indeed, after fixing
    the other coordinates the polynomial is affine in z_k with strictly
    positive slope, so its root is unique. The three-site calculation
    below separates the conditional law from this forced trace law at
    every positive common rate.

    For tilted insertion, x and s are the restored old and new terminal
    inverse rates. The erased singleton trace has inverse rate x/2.
    Its occupation probability 2/(x+2), tilted by y>0, becomes
    2*y/(x+2*y). The exact insertion kernel has positive output covariance
    above the recorded field threshold. This is a preservation obstruction,
    not a counterexample to stationary endpoint-path Rayleigh or NA.
    """
    x, s, r, y = sp.symbols("x s r y", positive=True)

    def assert_positive(expression, *variables):
        numerator, denominator = sp.fraction(sp.cancel(expression))
        for polynomial in (numerator, denominator):
            coefficients = sp.Poly(polynomial, *variables).coeffs()
            assert all(value >= 0 for value in coefficients)
            assert any(value > 0 for value in coefficients)

    def generating(weights, z0, z1):
        return weights[0]+weights[1]*z0+weights[2]*z1+weights[3]*z0*z1

    pair_kernel = graph_kernel(((1, 2), (0,)))
    triple_neighbors = ((1, 3), (0, 2), (1,))
    triple_kernel = graph_kernel(triple_neighbors)
    denominator = x**4+10*x**3+32*x**2+44*x+16
    weights = sp.Matrix([
        x**2*(x**2+4*x+2), 2*x*(x**2+4*x+2),
        2*x*(x+4)*(2*x+1), 4*(x**2+8*x+4),
    ])
    conditional = weights/denominator
    assert sp.expand(sum(weights)-denominator) == 0
    for value in conditional:
        assert_positive(value, x)
    mu = sp.Matrix(stationary(triple_kernel, (x,)*3))
    assert (mu[:4, 0]/sum(mu[:4, 0])-conditional).applyfunc(sp.cancel) == sp.zeros(4, 1)
    traced = sp.Matrix(stationary(trace_kernel(triple_kernel, [2]), (x, x)))
    erased = sp.Matrix(stationary(pair_kernel, (x, x/2)))
    assert (traced-erased).applyfunc(sp.cancel) == sp.zeros(4, 1)
    harmonic_zeros = ((-x/2, sp.S.Zero), (sp.Rational(1, 2), -x/4))
    for z0, z1 in harmonic_zeros:
        assert sp.cancel(generating(conditional, z0, z1)) == 0
        assert sp.cancel(generating(traced, z0, z1)) == 0
    full_gap = 4*x**3/(denominator*(x**2+3*x+4))
    assert sp.cancel(traced[3]-4/(x**2+3*x+4)) == 0
    assert sp.cancel(conditional[3]-traced[3]-full_gap) == 0
    assert_positive(full_gap, x)
    raw, triple_states = _endpoint_symbolic_absorption(
        [(1, 1, 1)], triple_kernel, (1/(1+x),)*3, reverse=True)
    assert (raw-mu).applyfunc(sp.cancel) == sp.zeros(8, 1)
    raw_one, triple_raw_counts = stationary_law(triple_kernel, (sp.Rational(1, 2),)*3)
    assert sp.Matrix(raw_one) == mu.subs(x, 1)
    assert conditional.subs(x, 1) == sp.Matrix([7, 14, 30, 52])/103
    assert traced.subs(x, 1) == sp.Matrix([3, 6, 13, 22])/44
    deletion = {
        "neighbors": triple_neighbors,
        "parameter": "x=1/lambda>0; all three sleep rates equal lambda",
        "conditioning": "X_2=0",
        "conditional_numerators": [str(value) for value in weights],
        "conditional_denominator": str(denominator),
        "harmonic_zeros": [[str(value) for value in point] for point in harmonic_zeros],
        "forced_inverse_rates_on_smaller_path": ["x", "x/2"],
        "trace_law": [str(value) for value in traced],
        "conditional_minus_trace_full_occupation": str(full_gap),
        "all_positive_rates_strict_gap_certified": True,
        "rate_only_conditional_stationary_closure": False,
        "symbolic_raw_check": {
            "toppling_order": "rightmost", "transient_states": triple_states,
            "full_absorption_residual_and_normalization": True,
            "agrees_with_IDLA_stationary_law": True,
        },
        "rate_one_check": {
            "conditional_law": [str(value) for value in conditional.subs(x, 1)],
            "trace_law": [str(value) for value in traced.subs(x, 1)],
            "full_occupation_gap": str(full_gap.subs(x, 1)),
            "raw_transient_states": list(triple_raw_counts),
            "both_toppling_orders_and_every_addition_equation": True,
        },
    }
    print("endpoint conditional deletion: all-rate strict obstruction and raw checks passed",
          flush=True)

    cc = s*x+s+x+2
    dd = s*x+2*s+2*x+2
    insertion = sp.Matrix([
        [s*x/dd, s*x*(2*s*x+2*s+x+2)/(2*cc*dd)],
        [2*s/dd, s*(2*s*x+2*s+x+2)/(cc*dd)],
        [2*(x+1)/dd, s*x*(3*x+4)/(2*cc*dd)],
        [0, (x+2)/cc],
    ])
    for value in insertion:
        if value != 0:
            assert_positive(value, x, s)
    assert all(sp.cancel(sum(insertion[:, j])-1) == 0 for j in range(2))
    symbolic_orders = []
    for reverse in (False, True):
        raw, count = _endpoint_symbolic_absorption(
            [(0, 1), (-1, 1)], pair_kernel, (1/(1+x), 1/(1+s)), reverse=reverse)
        assert (raw-insertion).applyfunc(sp.cancel) == sp.zeros(4, 2)
        symbolic_orders.append({
            "toppling_order": "rightmost" if reverse else "leftmost",
            "transient_states": count,
            "full_absorption_residual_and_normalization": True,
            "agrees_with_insertion_kernel": True,
        })
    input_occupation = 2/(x+2)
    stationary_input = sp.Matrix([1-input_occupation, input_occupation])
    pair_mu = sp.Matrix(stationary(pair_kernel, (x, s)))
    assert (insertion*stationary_input-pair_mu).applyfunc(sp.cancel) == sp.zeros(4, 1)
    output = (insertion*sp.Matrix([1-r, r])).applyfunc(sp.factor)
    gap = s*(x+1)*(2*cc-r*(x+2))*(4*cc-r*(s+2)*(x+2)**2)/(2*cc**2*dd**2)
    assert sp.cancel(output[1]*output[2]-output[0]*output[3]-gap) == 0
    tilted_occupation = 2*y/(x+2*y)
    threshold = 2*cc/(s*x+2*x+4)
    assert sp.cancel(
        (4*cc-r*(s+2)*(x+2)**2).subs(r, tilted_occupation)
        -2*x*(s*x+2*x+4)*(threshold-y)/(x+2*y)) == 0
    # For 0<=r<=1 the first gap factor is positive; the second has exactly
    # the sign of threshold-y. The threshold is strictly greater than one.
    assert_positive(2*cc-(x+2), x, s)
    assert_positive(threshold-1, x, s)
    rational_cases = []
    for tip_inverse, expected_output, expected_covariance in (
        (sp.S.One, sp.Matrix([19, 38, 34, 84])/175, sp.Rational(304, 30625)),
        (sp.Rational(1, 2), sp.Matrix([7, 14, 23, 66])/110, sp.Rational(7, 605)),
    ):
        substitution = {x: 1, s: tip_inverse}
        stationary_output = (insertion*stationary_input).subs(substitution)
        raw_one, counts = stationary_law(
            pair_kernel, (sp.Rational(1, 2), 1/(1+tip_inverse)))
        assert sp.Matrix(raw_one) == stationary_output
        tilted_output = output.subs({**substitution, r: sp.Rational(4, 5)})
        assert tilted_output == expected_output
        covariance = sp.factor(tilted_output[0]*tilted_output[3]
                               -tilted_output[1]*tilted_output[2])
        assert covariance == expected_covariance and covariance > 0
        rational_cases.append({
            "restored_sleep_rates": ["1", str(1/tip_inverse)],
            "erased_trace_sleep_rate": "2",
            "trace_input_occupation": "2/3", "external_field": "2",
            "tilted_input_occupation": "4/5",
            "stationary_output": [str(value) for value in stationary_output],
            "tilted_output": [str(value) for value in tilted_output],
            "positive_covariance": str(covariance),
            "raw_transient_states": list(counts),
            "both_toppling_orders_and_every_addition_equation": True,
        })
    tilted = {
        "neighbors": ((1, 2), (0,)),
        "parameters": "x,s>0 are restored inverse rates; y>0 is the input external field",
        "erased_trace_inverse_rate": "x/2",
        "trace_input_occupation": str(input_occupation),
        "tilted_input_occupation": str(tilted_occupation),
        "insertion_input_states": [[0, 1], [-1, 1]],
        "insertion_matrix": [[str(value) for value in insertion.row(i)] for i in range(4)],
        "rayleigh_gap_for_Bernoulli_r_input": str(gap),
        "field_threshold": str(threshold),
        "threshold_strictly_greater_than_one_certified": True,
        "output_covariance_sign": "positive iff y exceeds field_threshold; zero at equality",
        "all_positive_inverse_rates_sign_factorization_certified": True,
        "symbolic_raw_checks": symbolic_orders,
        "symbolic_stationary_input_agrees_with_IDLA": True,
        "rate_one_checks": rational_cases,
    }
    print("endpoint tilted insertion: all-rate threshold and both raw orders passed",
          flush=True)
    return {
        "scope": "specific leaf-recursion obstructions; no stationary endpoint-path NA counterexample",
        "two_site_output_mask_order": [[], [0], [1], [0, 1]],
        "recursion": "mu_m^(c) = I_m^(c) mu_(m-1)^(2); c=2 for growth, c=1 for the common-rate target",
        "restored_old_terminal_rate": (
            "2*lambda in the erased trace representation changes to lambda in the enlarged path; "
            "insertion does not freeze existing rates"
        ),
        "conditional_leaf_deletion": deletion,
        "tilted_leaf_insertion": tilted,
    }


def _endpoint_inverse(matrix):
    field = DomainMatrix.from_Matrix(matrix).to_field()
    inverse = field.inv()
    assert field.matmul(inverse) == DomainMatrix.from_Matrix(sp.eye(matrix.rows)).convert_to(field.domain)
    return inverse.to_Matrix()


def endpoint_sink_transfer(addition, rate, killing, right):
    """Exact transfer with a held active pool at a new sink-end vertex.

    ``addition`` is the old endpoint addition kernel (columns are inputs).
    Old sink losses return to the new vertex, whose jump probabilities are
    killing, right, and 1-killing-right for a self-loop; its sleep rate is
    rate. Old physical rates stay fixed. Total retained count is h and
    the active pool is h-|S|. The domain I_h consists of |S|<h.

    R_h=(killing+right)I-right*T[I_h,I_h]+rate*P_(h-1).
    F_h=(D_h+killing*F_(h-1)*K_h)*R_h^-1, where D_h records
    sleeping at the new root, its last particle dying, or a no-loss old
    addition consuming that particle. Output masks use the new root as
    bit m. F_(m+1)*nu_m gives the enlarged stationary law, because after
    first stabilizing the saturated old path the pool is m+1-|S|.
    An arbitrary input with one new active particle is a different sampler.

    This sink-end construction differs from far-leaf insertion: in that
    erased trace representation the OLD terminal rate is 2*lambda, and
    restoring the enlarged path changes it to lambda, not a frozen append.
    """
    size = addition.rows
    m = size.bit_length()-1
    assert size == 1 << m and addition.cols == size
    ranks = [mask.bit_count() for mask in range(size)]
    assert all(sp.cancel(sum(addition[:, mask])-1) == 0 for mask in range(size))
    assert all(addition[y, s] == 0 for y in range(size) for s in range(size)
               if ranks[y] > ranks[s]+1)
    maps, domains, inverses = {}, {}, {}
    for h in range(1, m+3):
        indices = [s for s in range(size) if ranks[s] < h]
        domains[h] = indices
        boundary = sp.diag(*(int(ranks[s] == h-1) for s in indices))
        resolvent = ((killing+right)*sp.eye(len(indices))
                     - right*addition.extract(indices, indices)+rate*boundary)
        inverses[h] = _endpoint_inverse(resolvent)
        direct = sp.zeros(2*size, len(indices))
        for col, s in enumerate(indices):
            if ranks[s] == h-1:
                direct[s, col] += killing
                direct[size+s, col] += rate
                for y in range(size):
                    if ranks[y] == h:
                        direct[y, col] += right*addition[y, s]
            if h > 1 and s in domains[h-1]:
                direct[:, col] += killing*maps[h-1][:, domains[h-1].index(s)]
        maps[h] = (direct*inverses[h]).applyfunc(sp.factor)
        assert all(sp.cancel(sum(maps[h][:, col])-1) == 0 for col in range(len(indices)))
    new_addition = sp.zeros(2*size)
    for s in range(size):
        for occupied in (0, 1):
            h = ranks[s]+occupied+1
            new_addition[:, s+occupied*size] = maps[h][:, domains[h].index(s)]
    return maps, new_addition, inverses, domains


def endpoint_hold_one(addition, killing, right):
    """Stop surplus routing at the first stage with exactly one held active.

    Applied to nu_m after the first old stabilization, the result is the
    stationary law of the trace over the new root. Its internal walk matrix
    is Q+[right/(killing+right)]*k*e_0^T, with old rates unchanged.
    """
    size = addition.rows
    m = size.bit_length()-1
    ranks = [s.bit_count() for s in range(size)]
    maps = {}
    for h in range(1, m+2):
        indices = [s for s in range(size) if ranks[s] < h]
        inner = [s for s in indices if ranks[s] < h-1]
        current = sp.zeros(size, len(indices))
        for col, s in enumerate(indices):
            if ranks[s] == h-1:
                current[s, col] = 1
        if inner:
            resolvent = ((killing+right)*sp.eye(len(inner))
                         - right*addition.extract(inner, inner))
            direct = killing*maps[h-1]
            for col, s in enumerate(inner):
                for y in range(size):
                    if ranks[y] == h-1:
                        direct[y, col] += right*addition[y, s]
            stopped = (direct*_endpoint_inverse(resolvent)).applyfunc(sp.factor)
            for col, s in enumerate(inner):
                current[:, indices.index(s)] = stopped[:, col]
        assert all(sp.cancel(sum(current[:, col])-1) == 0 for col in range(len(indices)))
        maps[h] = current
    return maps[m+1]


def endpoint_transfer_certificate():
    """One exact diagnostic, not a sweep or an all-size sign certificate.

    Common inverse rate t>0: recover the 2->3 transfer, the held-one trace,
    every new-root addition column, a=(kappa*b+rho*Tplus*b)/lambda,
    and Delta_(root,j)P=(rho/lambda)*(B*d_j C-C*d_j B).
    Check symbolic raw saturated absorption independently of IDLA; at rate
    one also check raw additions in both orders and every stationarity
    equation. The negative extreme rank pair uses the actual generated b.
    """
    from algebra import basis_matrix, idla_operators
    from exact import absorption, add, stable

    def equal(left, right):
        assert (left-right).applyfunc(sp.cancel) == sp.zeros(*left.shape)

    def addition_data(kernel, inverse_rates, port):
        operators, jumps = idla_operators(kernel)
        basis = basis_matrix(operators, jumps, inverse_rates)
        inverse = _endpoint_inverse(basis)
        return ((inverse*operators[port]*basis).applyfunc(sp.factor), inverse[:, -1])

    def generating(vector, fields):
        return sum(value*sp.prod(z for j, z in enumerate(fields) if mask >> j & 1)
                   for mask, value in enumerate(vector))

    t = sp.Symbol("t", positive=True)
    rate, killing, right = 1/t, sp.Rational(1, 2), sp.Rational(1, 2)
    old_neighbors = ((1, 2), (0,))
    # Computational labeling: new root 2, old vertices 0--1, true sink 3.
    new_neighbors = ((1, 2), (0,), (0, 3))
    old_kernel, new_kernel = map(graph_kernel, (old_neighbors, new_neighbors))
    addition, nu = addition_data(old_kernel, (t, t), 0)
    maps, new_addition, inverses, domains = endpoint_sink_transfer(addition, rate, killing, right)
    mu = (maps[3]*nu).applyfunc(sp.factor)
    independent_addition, independent_mu = addition_data(new_kernel, (t,)*3, 2)
    equal(mu, independent_mu)
    equal(new_addition, independent_addition)
    equal(new_addition*mu, mu)
    held = (endpoint_hold_one(addition, killing, right)*nu).applyfunc(sp.factor)
    equal(held, sp.Matrix(stationary(trace_kernel(new_kernel, [2]), (t, t))))
    print("endpoint transfer: all-rate law, addition kernel, and held-one trace passed", flush=True)

    occupation = {3: inverses[3]*nu}
    for h in (2, 1):
        restriction = occupation[h+1].extract([domains[h+1].index(s) for s in domains[h]], [0])
        occupation[h] = killing*inverses[h]*restriction
    levels = []
    for r in range(3):
        level = sp.zeros(4, 1)
        for col, s in enumerate(domains[r+1]):
            if s.bit_count() == r:
                level[s] = sp.factor(rate*occupation[r+1][col])
        levels.append(level)
    b = sum(levels, sp.zeros(4, 1))
    a = mu[:4, 0]
    equal(b, mu[4:, 0])
    plus = sp.Matrix(4, 4, lambda y, s: addition[y, s]
                     if int(y).bit_count() == int(s).bit_count()+1 else 0)
    c = (plus*b).applyfunc(sp.factor)
    equal(a, (killing*b+right*c)/rate)
    beta, gamma = sum(b), sum(c)
    assert sp.cancel((rate+killing)*beta+right*gamma-rate) == 0
    g = sp.Matrix([[0, 0, 1, 1]])
    gap = sp.factor(beta*(g*c)[0]-gamma*(g*b)[0])
    covariance = sp.factor((g*b)[0]-beta*(g*(a+b))[0])
    assert sp.cancel(covariance+right*gap/rate) == 0
    z0, z1, zo = sp.symbols("z0 z1 zo")
    bb, cc = generating(b, (z0, z1)), generating(c, (z0, z1))
    polynomial = generating(mu, (z0, z1, zo))
    assert sp.cancel(polynomial-(zo+killing/rate)*bb-right*cc/rate) == 0
    for z in (z0, z1):
        delta = sp.diff(polynomial, zo)*sp.diff(polynomial, z)-polynomial*sp.diff(polynomial, zo, z)
        assert sp.cancel(delta-right*(bb*sp.diff(cc, z)-cc*sp.diff(bb, z))/rate) == 0
    print("endpoint transfer: all-rate boundary flux and covariance/Rayleigh identities passed", flush=True)

    raw_symbolic, symbolic_states = _endpoint_symbolic_absorption(
        [(1, 1, 1)], new_kernel, (1/(1+t),)*3, reverse=True)
    equal(raw_symbolic, mu)
    print("endpoint transfer: all-rate raw saturated absorption passed", flush=True)
    raw_checks = []
    starts = [add(stable(s, 3), 2) for s in range(8)]
    for reverse in (False, True):
        raw, count = absorption(starts, new_kernel, (sp.Rational(1, 2),)*3, reverse=reverse)
        equal(raw.T, new_addition.subs(t, 1))
        raw_checks.append({"order": "rightmost" if reverse else "leftmost", "transient_states": count})
    raw_mu, raw_counts = stationary_law(new_kernel, (sp.Rational(1, 2),)*3)
    equal(sp.Matrix(raw_mu), mu.subs(t, 1))
    levels_one = [v.subs(t, 1) for v in levels]
    plus_one = plus.subs(t, 1)
    terms = sp.Matrix(3, 3, lambda r, s:
                      sum(levels_one[r])*(g*plus_one*levels_one[s])[0]
                      - sum(plus_one*levels_one[s])*(g*levels_one[r])[0])
    assert sum(terms) == gap.subs(t, 1)
    assert terms[2, 0] == -sp.Rational(7792, 1156805) and terms[0, 2] == 0
    assert covariance.subs(t, 1) == -sp.Rational(81494093, 17595004050)
    assert terms[2, 0] == -levels_one[2][3]*levels_one[0][0]*plus_one[1, 0]
    print("endpoint transfer: rate-one raw checks and negative extreme rank pair passed", flush=True)
    return {
        "scope": "all-size transfer formulas; symbolic 2-to-3 verification only; all-size endpoint sign remains open",
        "parameter": "t=1/lambda>0; common physical sleep rate; kappa=rho=1/2",
        "old_neighbors": old_neighbors, "new_neighbors": new_neighbors,
        "new_root_bit": 2, "far_endpoint_bit": 1,
        "input": "nu_m with root pool m+1-|S|; stopped held-one input is the stationary trace law",
        "R_h": "(kappa+rho)I-rho*T[I_h,I_h]+lambda*P_(h-1); I_h={S: |S|<h}",
        "F_h": "(D_h+kappa*F_(h-1)*K_h)*R_h^(-1); D_h=[kappa*E_h*P_(h-1)+rho*P_h*T*E_h; lambda*E_h*P_(h-1)]",
        "law_recursion": "nu_(m+1)=F_(m+1)*nu_m",
        "addition_recursion": "T_(m+1)e_(beta,S)=F_(|S|+beta+1)e_S; relabel new root to 0 before next step",
        "old_terminal_rate_caution": "far-leaf erased trace uses 2*lambda at the old terminal; restoring that path uses lambda",
        "symbolic_law": [str(v) for v in mu],
        "symbolic_held_one_trace": [str(v) for v in held],
        "symbolic_b": [str(v) for v in b],
        "symbolic_no_loss_image": [str(v) for v in c],
        "boundary_flux": "a=(kappa*b+rho*Tplus*b)/lambda; Tplus(Y,S)=T(Y,S)*1_(|Y|=|S|+1)",
        "covariance_identity": "-Cov(X_root,g)=(rho/lambda)*(mass(b)*g^T*Tplus*b-mass(Tplus*b)*g^T*b)",
        "rayleigh_identity": "Delta_(root,j)P=(rho/lambda)*(B*partial_j(C)-C*partial_j(B))",
        "symbolic_far_endpoint_covariance": str(covariance),
        "symbolic_checks": {"law_and_all_root_addition_columns": True, "held_one_stationary_trace": True,
                            "boundary_flux_and_covariance": True, "both_root_Rayleigh_identities": True,
                            "raw_saturated_order": "rightmost", "raw_transient_states": symbolic_states},
        "rate_one": {"held_one_trace": [str(v) for v in held.subs(t, 1)],
                     "b": [str(v) for v in b.subs(t, 1)], "no_loss_image": [str(v) for v in c.subs(t, 1)],
                     "rank_pair_terms": [[str(v) for v in terms.row(r)] for r in range(3)],
                     "far_endpoint_gap": str(gap.subs(t, 1)), "far_endpoint_covariance": str(covariance.subs(t, 1)),
                     "raw_root_addition_checks": raw_checks, "raw_stationary_transient_states": raw_counts,
                     "all_addition_stationarity_checked": True},
        "all_size_extreme_rank_obstruction": "m>=2: G_(m,0)=-b^(m)(V)*b^(0)(empty)*sum_(j != far)Tplus({j},empty)<0; G_(0,m)=0",
        "remaining_gap": "on generated b, prove mass(b)*g^T*Tplus*b >= mass(Tplus*b)*g^T*b for the far endpoint at all m,lambda",
    }


def _endpoint_ordinary_cone_obstruction():
    """The normalized entrance cone fails even on the common-rate path.

    This rate-one two-to-three-site model was independently raw-validated
    by endpoint_transfer_certificate. Here only its proposed cone is tested.
    """
    from algebra import basis_matrix, idla_operators

    Q = sp.Rational
    kernel = (((1,Q(1,2)),(2,Q(1,2))),((0,sp.S.One),))
    A, jumps = idla_operators(kernel)
    M = basis_matrix(A,jumps,(sp.S.One,sp.S.One))
    D = sp.Matrix(4,4,lambda i,j: M[i,j]
                  if int(i).bit_count() == int(j).bit_count() else 0)
    U = D.inv()
    W = sp.eye(4)-M*U
    L = sp.Matrix(4,4,lambda i,j: A[0][i,j]
                  if int(i).bit_count() == int(j).bit_count()+1 else 0)
    J = sp.diag(*A[0].diagonal())
    P = sp.zeros(4)
    P[3,3] = 1
    R = A[1]-P
    F, Z = (sp.eye(4)+L)/2, sp.eye(4)-J/2
    H = sp.eye(4)+Z-F*Z.inv()*W-W*F
    q = sp.ones(1,4)*U
    p = sp.Matrix([[0,0,1,1]])*U
    normalizer = q*(sp.eye(4)+F)
    top = 1/H[3,3]
    eta = 1-Q(3,2)*q[3]*top
    v2 = eta/(normalizer[2]+eta*normalizer[0])
    candidate = sp.Matrix([eta*v2,0,v2,top])
    assert eta == Q(19,37)
    assert all(v >= 0 for v in candidate)
    assert all(v >= 0 for v in eta*candidate-R*candidate)
    assert (normalizer*candidate)[0] == 1
    assert H*candidate != P[:,3]
    gap = ((q*candidate)[0]*(p*L*candidate)[0]
           -(q*L*candidate)[0]*(p*candidate)[0])
    assert gap == -Q(146170705,13661635689)
    return {
        "scope": "ordinary three-site endpoint path, all sleep rates one",
        "eta": str(eta), "artificial_entrance": list(map(str,candidate)),
        "positive_and_R_subinvariant": True,
        "actual_top_mass_and_terminal_normalization": True,
        "endpoint_bracket": str(gap), "is_actual_generated_entrance": False,
        "interpretation": "normalized entrance cone does not suffice even in the main model",
    }


def endpoint_entrance_certificate():
    """One nonunit normalization check for the analytic path entrance proof.

    Two old sites have unequal sleep rates. The inserted root has a third
    rate, unequal sink/return probabilities, and holding. Compare every
    all-active input with raw ARW absorption in both toppling orders.
    The all-size theorem is analytic; this finite case checks its formulas.
    """
    from algebra import basis_matrix, idla_operators
    from exact import absorption, add, stable

    Q = sp.Rational
    kappa, rho, root_rate = Q(1, 3), Q(1, 4), Q(5, 4)
    x_o = 1/root_rate
    outside_rates = (Q(3, 2), Q(2, 3))
    old_kernel = (((1, Q(1, 2)), (2, Q(1, 2))), ((0, sp.S.One),))
    new_kernel = (
        ((3, kappa), (1, rho), (0, 1-kappa-rho)),
        ((0, Q(1, 2)), (2, Q(1, 2))),
        ((1, sp.S.One),),
    )

    def data(kernel, rates):
        operators, jumps = idla_operators(kernel)
        matrix = basis_matrix(operators, jumps, tuple(1/r for r in rates))
        inverse = _endpoint_inverse(matrix)
        assert matrix*inverse == sp.eye(matrix.rows)
        return matrix, inverse, operators

    M, N, A = data(old_kernel, outside_rates)
    rates = (root_rate,)+outside_rates
    Mp, Np, Ap = data(new_kernel, rates)
    D = sp.Matrix(4, 4, lambda i,j: M[i,j]
                  if int(i).bit_count() == int(j).bit_count() else 0)
    U = D.inv()
    W = sp.eye(4)-M*U
    assert all(v >= 0 for v in W)
    assert all(W[i,j] == 0 for i in range(4) for j in range(4)
               if i.bit_count() != j.bit_count()-1)
    assert N == U+N*W and D*N == sp.eye(4)+W+W**2
    L = sp.Matrix(4, 4, lambda i,j: A[0][i,j]
                  if int(i).bit_count() == int(j).bit_count()+1 else 0)
    J = sp.diag(*A[0].diagonal())
    assert A[0] == L+J
    F = kappa*sp.eye(4)+rho*L
    Z = (kappa+rho)*sp.eye(4)-rho*J
    T_idla = F*Z.inv()
    C = sp.eye(4)+x_o*Z
    vacant, occupied = [0,2,4,6], [1,3,5,7]
    order = vacant+occupied
    expected = sp.Matrix.vstack(sp.Matrix.hstack(D, -x_o*F*D),
                                sp.Matrix.hstack(-W*D, (C-T_idla*W)*D))
    assert Mp.extract(order,order) == expected
    assert Ap[0].extract(occupied,occupied) == T_idla
    H = C-T_idla*W-x_o*W*F
    q = sp.ones(1,4)*U
    assert all(v > 0 for v in q)
    assert sp.ones(1,4)*H == q*(sp.eye(4)+x_o*F)
    assert all(H[i,j] <= 0 for i in range(4) for j in range(4) if i != j)
    assert all(v > 0 for v in sp.ones(1,4)*H)
    entrance = D*Np.extract(occupied,order)
    assert entrance == sp.Matrix.hstack(H.inv()*W,H.inv())
    assert all(v >= 0 for v in entrance)

    # The same point also checks the far-end first-loss commutator. Its
    # all-size sign follows from the legal first-loss factorization.
    P = sp.zeros(4)
    P[3,3] = 1
    R = A[1]-P
    assert L*R == R*L and J*R-R*J == P*L
    C_f = P+W*R-R*W
    assert all(v >= 0 for v in C_f)
    T_f = N*A[1]*M
    assert D*T_f*U == R+(sp.eye(4)-W).inv()*C_f
    assert all(v >= 0 for v in D*T_f*U)
    B = (T_idla*R-R*T_idla)*W
    b = B[3,3]
    assert B == b*P and 0 <= b <= 1-q[3]
    commutator = H*R-R*H
    assert commutator == (1+x_o*kappa)*P-B-T_idla*C_f-x_o*C_f*F
    assert all(commutator[i,j] <= 0 for i in range(4) for j in range(3))
    assert commutator[:,3] == ((1+x_o*kappa)*(1-q[3])-b)*P[:,3]
    assert all(v >= 0 for v in commutator[:,3])
    c = 1+x_o*(kappa+rho)
    root_return = (T_idla*W+x_o*kappa*W+x_o*rho*J+x_o*rho*W*L)/c
    assert sp.eye(4)-root_return == H/c
    assert all(v >= 0 for v in root_return)

    # Decide whether the new entrance cone, together with its actual top
    # mass and terminal normalization, already implies the desired sign.
    # This constructed vector satisfies those constraints but is not d.
    d = H.inv()[:,3]
    eta = 1-(1+x_o*kappa)*q[3]*d[3]
    assert all(v >= 0 for v in eta*d-R*d)
    normalizer = q*(sp.eye(4)+x_o*F)
    v2 = eta/(normalizer[2]+eta*normalizer[0])
    relaxed = sp.Matrix([eta*v2,0,v2,d[3]])
    assert 0 < eta < 1 and all(v >= 0 for v in relaxed)
    assert all(v >= 0 for v in eta*relaxed-R*relaxed)
    assert (normalizer*relaxed)[0] == 1
    assert H*relaxed != P[:,3]
    far_row = sp.Matrix([[0,0,1,1]])*U
    relaxed_gap = ((q*relaxed)[0]*(far_row*L*relaxed)[0]
                   -(q*L*relaxed)[0]*(far_row*relaxed)[0])
    assert relaxed_gap < 0

    starts = [tuple(int(bool(mask & (1 << j))) for j in range(3)) for mask in range(8)]
    probabilities = tuple(rate/(1+rate) for rate in rates)
    counts = []
    for reverse in (False,True):
        raw, count = absorption(starts,new_kernel,probabilities,reverse=reverse)
        assert raw.T == Np
        counts.append(count)

    far_starts = [add(stable(mask,2),1) for mask in range(4)]
    far_counts = []
    for reverse in (False,True):
        raw, count = absorption(far_starts,old_kernel,
                                tuple(r/(1+r) for r in outside_rates),reverse=reverse)
        assert raw.T == T_f
        no_loss = sp.Matrix(4,4,lambda i,j: raw[j,i]
                            if int(i).bit_count() == int(j).bit_count()+1 else 0)
        assert no_loss == U*R*D
        far_counts.append(count)

    def rows(matrix):
        return [[str(matrix[i,j]) for j in range(matrix.cols)] for i in range(matrix.rows)]

    print("endpoint entrance: nonunit root, unequal outside rates, holds, and both raw orders passed",
          flush=True)
    return {
        "status": "exact normalization check passed; all-size proof is analytic in report",
        "scope": "one two-to-three-site path; all eight all-active inputs",
        "root_sleep_rate": str(root_rate), "root_inverse_sleep_rate": str(x_o),
        "outside_sleep_rates": list(map(str,outside_rates)),
        "root_sink_probability": str(kappa), "root_return_probability": str(rho),
        "root_hold_probability": str(1-kappa-rho),
        "rank_bidiagonal_identity": True, "first_loss_kernel_nonnegative": True,
        "root_block_identity": True, "strict_column_M_matrix_identity": True,
        "all_active_entrance_identity_and_nonnegativity": True,
        "all_columns_match_raw_ARW_both_orders": True, "raw_transient_counts": counts,
        "W": rows(W), "U": rows(U), "H": rows(H), "entrance": rows(entrance),
        "far_end_first_loss_commutator": {
            "identity_and_nonnegativity": True,
            "H_commutator_boundary_signs": True,
            "root_renewal_identity": True,
            "far_addition_matches_raw_both_orders": True,
            "raw_transient_counts": far_counts,
            "R": rows(R), "C_f": rows(C_f), "H_R_minus_R_H": rows(commutator),
        },
        "normalized_entrance_cone_obstruction": {
            "eta": str(eta), "artificial_entrance": list(map(str,relaxed)),
            "positive_and_R_subinvariant": True,
            "actual_top_mass_and_terminal_normalization": True,
            "endpoint_bracket": str(relaxed_gap),
            "is_actual_generated_entrance": False,
            "interpretation": "these necessary constraints alone do not imply endpoint NC",
        },
        "ordinary_normalized_entrance_cone_obstruction": _endpoint_ordinary_cone_obstruction(),
        "NA_or_covariance_theorem": False,
    }


def endpoint_return_reversibility_certificate():
    """One rate-one four-path discriminator for rank-level reversibility.

    The nonzero cycle products refute that proposed mechanism, not NA.
    Uses the established IDLA/insertion identities; no new raw absorption.
    """
    from algebra import basis_matrix, idla_operators

    kernel = graph_kernel(((1, 4), (0, 2), (1, 3), (2,)))
    operators, jumps = idla_operators(kernel)
    matrix = basis_matrix(operators, jumps, (sp.S.One,)*4)
    inverse = matrix.inv()
    assert matrix*inverse == sp.eye(16)
    assert min(inverse) >= 0
    assert all(sum(inverse[:, j]) == 1 for j in range(16))
    additions = [inverse*a*matrix for a in operators]
    for addition in additions:
        assert min(addition) >= 0
        assert all(sum(addition[:, j]) == 1 for j in range(16))
        assert addition*inverse[:, 15] == inverse[:, 15]
    occupied_insertion = additions[3][8:16, 0:8]
    returning = 2*sp.eye(8)-occupied_insertion.inv()
    assert min(returning) >= 0
    assert all(sum(returning[:, j]) < 1 for j in range(8))
    assert all(returning[i, j] == 0 for i in range(8) for j in range(8)
               if i.bit_count() > j.bit_count())
    saved = json.loads(Path(__file__).with_name("endpoint_path_certificates.json").read_text())["cases"]["4"]
    parameter = sp.Symbol("t")
    denominator = sp.sympify(saved["denominator"]).subs(parameter, 1)
    assert [sp.sympify(v).subs(parameter, 1)/denominator for v in saved["numerators"]] == list(inverse[:, 15])
    expected = {1: sp.Rational(-225, 316485232), 2: sp.Rational(-874237, 15320200968)}
    blocks = []
    for rank in (1, 2):
        ids = [i for i in range(8) if i.bit_count() == rank]
        a, b, c = ids
        assert all(returning[i, j] > 0 for i in ids for j in ids if i != j)
        difference = (returning[b, a]*returning[c, b]*returning[a, c]
                      -returning[c, a]*returning[b, c]*returning[a, b])
        assert difference == expected[rank]
        blocks.append({"rank": rank, "outside_masks": ids,
                       "block": [list(map(str, row)) for row in returning.extract(ids, ids).tolist()],
                       "cycle_difference": str(difference),
                       "reversible_with_positive_measure": False})
    return {"graph": "sink--0--1--2--3", "common_sleep_rate": "1",
            "inserted_vertex": 3, "matrix_convention": "row output, column input",
            "E": [list(map(str, row)) for row in returning.tolist()],
            "rank_blocks": blocks,
            "checks": {"basis_inverse_residual": True,
                       "N_and_all_addition_kernels_stochastic": True,
                       "all_four_addition_stationarity_equations": True,
                       "E_nonnegative_strictly_substochastic_rank_nonincreasing": True,
                       "matches_saved_stationary_certificate_at_t_1": True,
                       "new_raw_absorption_run": False}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", help="Optional JSON certificate file")
    parser.add_argument("--case", type=int, help="Check only this canonical graph bitmask")
    parser.add_argument("--five-site", choices=["rooted_fork", "cycle_below_cut", "shared_conditional_base"],
                        help="Test one selected five-site graph instead of the classification")
    parser.add_argument("--sitewise-three", action="store_true",
                        help="Test separate sleep rates on the three-site Rayleigh counterexample")
    parser.add_argument("--endpoint-path", type=int, choices=range(2, 6),
                        help="Certify NA and Rayleigh on a path with a sink at one endpoint")
    parser.add_argument("--endpoint-recursion", action="store_true",
                        help="Check exact endpoint-deletion and tilted-insertion obstructions only; "
                             "--output merges a dedicated key, preserving existing cases")
    parser.add_argument("--endpoint-transfer", action="store_true",
                        help="Check the sink-end transfer and its exact identities only; "
                             "--output merges a dedicated key, preserving existing cases")
    parser.add_argument("--endpoint-entrance", action="store_true",
                        help="Check the path first-loss and root-entrance formulas at one nonunit point; "
                             "--output preserves earlier certificate keys")
    parser.add_argument("--endpoint-return-reversibility", action="store_true",
                        help="Refute rank-level reversibility of the four-path return kernel at rate one")
    args = parser.parse_args()
    if args.endpoint_recursion or args.endpoint_transfer or args.endpoint_entrance or args.endpoint_return_reversibility:
        selector = ("--endpoint-return-reversibility" if args.endpoint_return_reversibility else
                    "--endpoint-entrance" if args.endpoint_entrance else
                    "--endpoint-transfer" if args.endpoint_transfer else "--endpoint-recursion")
        if (args.endpoint_path is not None or args.five_site is not None
                or args.case is not None or args.sitewise_three
                or sum((args.endpoint_recursion,args.endpoint_transfer,args.endpoint_entrance,
                        args.endpoint_return_reversibility)) != 1):
            parser.error(selector+" cannot be combined with other case selectors")
        if args.endpoint_return_reversibility:
            key, certificate = "endpoint_return_reversibility", endpoint_return_reversibility_certificate()
            print("rank-level return reversibility: FALSE; exact cycle products checked", flush=True)
        elif args.endpoint_entrance:
            key, certificate = "endpoint_entrance", endpoint_entrance_certificate()
        elif args.endpoint_transfer:
            key, certificate = "endpoint_transfer", endpoint_transfer_certificate()
        else:
            key, certificate = "endpoint_recursion_obstructions", endpoint_recursion_certificate()
        if args.output:
            path = Path(args.output)
            result = json.loads(path.read_text()) if path.exists() else {}
            result[key] = certificate
            path.write_text(json.dumps(result, indent=2)+"\n")
        return
    t = sp.Symbol("t")
    if args.endpoint_path is not None:
        assert args.five_site is None and args.case is None and not args.sitewise_three
        n = args.endpoint_path
        neighbors = tuple(tuple(([i-1] if i else [n])+([i+1] if i+1 < n else []))
                          for i in range(n))
        mu = stationary(graph_kernel(neighbors), (t,)*n)
        print("endpoint path", n, "exact stationary law solved", flush=True)
        case = polynomial_certificate(mu, t)
        case["neighbors"] = neighbors
        case["Rayleigh"] = rayleigh_polynomial_certificate(case, t)
        print("all-rate NA", case["all_rates_NA_certified"],
              "all-rate Rayleigh", case["Rayleigh"]["all_rates_Rayleigh_certified"], flush=True)
        if args.output:
            path = Path(args.output)
            result = json.loads(path.read_text()) if path.exists() else {"cases": {}}
            result["parameter"] = "t=1/lambda"
            result["scope"] = "specified endpoint path lengths; no all-length claim"
            result["cases"][str(n)] = case
            path.write_text(json.dumps(result, indent=2)+"\n")
        return
    if args.sitewise_three:
        assert args.five_site is None and args.case is None
        variables = sp.symbols("t0:3")
        neighbors = ((1, 3), (0, 2, 3), (1,))
        print("sitewise coefficient test on the three-site conditional-correlation example", flush=True)
        case = polynomial_certificate(stationary(graph_kernel(neighbors), variables), variables)
        assert case["comparisons"] == 15
        case["neighbors"] = neighbors
        print("all positive sleep rates: NA", case["all_rates_NA_certified"],
              "failed coefficient tests", len(case["failed_coefficient_tests"]), flush=True)
        if args.output:
            result = {"parameters": "t_i=1/lambda_i for i=0,1,2",
                      "scope": "one three-site graph with independent positive sleep rates",
                      "case": case}
            Path(args.output).write_text(json.dumps(result, indent=2) + "\n")
        return
    if args.five_site:
        from vacancy import adjacency
        edges, sinks = {
            "rooted_fork": ([(0, 1), (0, 2), (1, 3), (2, 4)], {0}),
            "cycle_below_cut": ([(0, 1), (1, 2), (1, 3), (2, 3), (3, 4)], {0}),
            # Each arm 0--1--2 and 0--3--4 has the original three-site
            # Rayleigh counterexample's sink attachments. This asks whether
            # sharing site 0 can convert conditional correlations into an
            # unconditional NA obstruction.
            "shared_conditional_base": ([(0, 1), (1, 2), (0, 3), (3, 4)], {0, 1, 3}),
        }[args.five_site]
        neighbors = adjacency(5, edges, sinks)
        print("targeted five-site all-rate coefficient test", args.five_site, flush=True)
        mu = stationary(graph_kernel(neighbors), (t,)*5)
        print("exact stationary law solved", flush=True)
        case = polynomial_certificate(mu, t)
        assert case["comparisons"] == 2280
        case["neighbors"] = neighbors
        print("all-rate NA", case["all_rates_NA_certified"],
              "failed coefficient tests", len(case["failed_coefficient_tests"]), flush=True)
        if args.output:
            output = Path(args.output)
            result = json.loads(output.read_text()) if output.exists() else {
                "parameter": "t=1/lambda", "scope": "selected motifs, not an exhaustive five-site classification",
                "cases": {}}
            assert result.get("scope") == "selected motifs, not an exhaustive five-site classification"
            result["cases"][args.five_site] = case
            output.write_text(json.dumps(result, indent=2) + "\n")
        return
    types = sorted(graph_types())
    result = {"parameter": "t=1/lambda", "graph_type_count": len(types), "cases": {}}
    print("connected internal four-site graph/sink types", len(types), flush=True)
    for index, (label, neighbors) in enumerate(types, 1):
        if args.case is not None and label != args.case:
            continue
        mu = stationary(graph_kernel(neighbors), (t,)*4)
        case = polynomial_certificate(mu, t)
        assert case["comparisons"] == 174
        case["neighbors"] = neighbors
        result["cases"][str(label)] = case
        print(index, "/", len(types), "graph", label,
              "all-rate NA", case["all_rates_NA_certified"],
              "failed coefficient tests", len(case["failed_coefficient_tests"]), flush=True)
        if args.output:
            # Preserve completed cases if this foreground research run is
            # interrupted. The coverage flag is recomputed from actual cases.
            result["complete_classification"] = len(result["cases"]) == len(types)
            Path(args.output).write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()
