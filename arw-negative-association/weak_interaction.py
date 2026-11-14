#!/usr/bin/env python3
"""Exact weak-interaction diagnostics in squarefree directed-edge variables.

Question: which leading round-trip coefficients survive when a proposed
trace-moment identity fails? The unrestricted identity is false on complete3
and diamond4. Paths, a star, and two directed cycles sharing one vertex
test its restricted use on chains of cycles. Finite checks are not a
general proof. The quotient sets z_e**2=0, not z_e**2=z_e.
"""

import argparse
import json
from pathlib import Path

from exact import sp


def nilpotent_operations(variables):
    """Canonical polynomial operations with every variable's square zero."""
    def cut(expression):
        polynomial = sp.Poly(sp.expand(expression), *variables)
        return sum(coefficient*sp.prod(variables[i] for i, exponent in enumerate(exponents) if exponent)
                   for exponents, coefficient in polynomial.terms() if max(exponents) <= 1)

    def cm(matrix):
        return matrix.applyfunc(cut)

    def inverse(matrix):
        constant = matrix.subs(dict.fromkeys(variables, 0))
        constant_inverse = constant.inv()
        step = cm(-constant_inverse*(matrix-constant))
        power = constant_inverse
        result = power
        for _ in variables:
            power = cm(step*power)
            result += power
        result = cm(result)
        assert cm(matrix*result-sp.eye(matrix.rows)) == sp.zeros(matrix.rows)
        return result

    return cut, cm, inverse


def nilpotent_stationary(transition, sleep_rates, variables):
    """Return all edge-monomial coefficient vectors of the exact law.

    Nonzero constant internal transitions are allowed. Matrix inversion
    uses a finite Neumann expansion around the constant matrix, so no
    rational functions of nilpotent variables enter the representation.
    """
    n, terms = len(sleep_rates), 1 << len(variables)
    size = 1 << n
    x = [1/sp.Rational(rate) for rate in sleep_rates]
    assert all(value > 0 for value in x)
    _, cm, inverse = nilpotent_operations(variables)
    assert transition.shape == (n, n+1)
    assert all(sp.expand(sum(transition[i, :])-1) == 0 for i in range(n))
    operators = [sp.zeros(size) for _ in range(n)]
    for mask in range(size):
        inside = [i for i in range(n) if mask >> i & 1]
        outside = [i for i in range(n+1) if i == n or not mask >> i & 1]
        if inside:
            kernel = transition.extract(inside, inside)
            green = inverse(sp.eye(len(inside))-kernel)
            hits = cm(green*transition.extract(inside, outside))
            assert cm((sp.eye(len(inside))-kernel)*hits
                      -transition.extract(inside, outside)) == sp.zeros(len(inside), len(outside))
        for vertex in range(n):
            if vertex not in inside:
                operators[vertex][mask | (1 << vertex), mask] = 1
            else:
                for j, target in enumerate(outside):
                    operators[vertex][mask if target == n else mask | (1 << target), mask] += hits[inside.index(vertex), j]
    jumps = [cm(transition[i, n]*sp.eye(size)
                +sum((transition[i, j]*operators[j] for j in range(n)), sp.zeros(size)))
             for i in range(n)]
    for i in range(n):
        assert cm(operators[i]*operators[i]-operators[i]*jumps[i]) == sp.zeros(size)
        for j in range(i):
            assert cm(operators[i]*operators[j]-operators[j]*operators[i]) == sp.zeros(size)
    basis_operators = [(1+x[i])*operators[i]-x[i]*jumps[i] for i in range(n)]
    columns = [sp.eye(size)[:, 0]]
    for mask in range(1, size):
        vertex = (mask & -mask).bit_length()-1
        columns.append(cm(basis_operators[vertex]*columns[mask ^ (1 << vertex)]))
    matrix = sp.Matrix.hstack(*columns)
    coefficients = [sp.zeros(size) for _ in range(terms)]
    for i in range(size):
        for j in range(size):
            for exponents, coefficient in sp.Poly(matrix[i, j], *variables).terms():
                coefficients[sum(exponent << k for k, exponent in enumerate(exponents))][i, j] = coefficient
    inverse_zero = coefficients[0].inv()
    law = [inverse_zero*sp.eye(size)[:, -1]]
    for mask in range(1, terms):
        rhs, subset = sp.zeros(size, 1), mask
        while subset:
            rhs += coefficients[subset]*law[mask ^ subset]
            subset = (subset-1) & mask
        law.append(-inverse_zero*rhs)
        assert sum(law[-1]) == 0
    assert sum(law[0]) == 1
    return law


def squarefree_moments(n, arcs, sleep_rates):
    """Compare every moment coefficient in Q[z_e]/(z_e**2 : e in arcs)."""
    arcs = tuple(arcs)
    assert len(set(arcs)) == len(arcs) and all(i != j for i, j in arcs)
    m, size = len(arcs), 1 << n
    variables = sp.symbols(f"z:{m}")
    terms = 1 << m
    x = [1/sp.Rational(rate) for rate in sleep_rates]
    transition = sp.zeros(n, n+1)
    for edge, (i, j) in enumerate(arcs):
        transition[i, j] = variables[edge]
    for i in range(n):
        transition[i, n] = 1-sum(transition[i, :n])
    law = nilpotent_stationary(transition, sleep_rates, variables)

    def vector(expression):
        result = [sp.S.Zero]*terms
        for exponents, coefficient in sp.Poly(expression, *variables).terms():
            if max(exponents) <= 1:
                result[sum(exponent << k for k, exponent in enumerate(exponents))] = coefficient
        return result

    def multiply(left, right):
        result = []
        for mask in range(terms):
            value, subset = sp.S.Zero, mask
            while True:
                value += left[subset]*right[mask ^ subset]
                if subset == 0:
                    break
                subset = (subset-1) & mask
            result.append(value)
        return result

    def inverse(values):
        result = [1/values[0]]
        for mask in range(1, terms):
            value, subset = sp.S.Zero, mask
            while subset:
                value += values[subset]*result[mask ^ subset]
                subset = (subset-1) & mask
            result.append(-value/values[0])
        return result

    laplacian = sp.eye(n)-transition[:, :n]
    determinants = {}
    for mask in range(size):
        vertices = [i for i in range(n) if mask >> i & 1]
        determinants[mask] = vector(sp.expand(laplacian.extract(vertices, vertices).det()))
    failures, moments = [], {}
    for occupied in range(1, size):
        denominator, subset = [sp.S.Zero]*terms, occupied
        while True:
            term = multiply(determinants[size-1], inverse(determinants[(size-1) ^ subset]))
            weight = sp.prod(x[i] for i in range(n) if subset >> i & 1)
            denominator = [a+weight*b for a, b in zip(denominator, term)]
            if subset == 0:
                break
            subset = (subset-1) & occupied
        traced = inverse(denominator)
        moments[occupied] = [sum(law[mask][state] for state in range(size) if state & occupied == occupied)
                             for mask in range(terms)]
        for mask, (actual, comparison) in enumerate(zip(moments[occupied], traced)):
            if actual != comparison:
                failures.append({"occupied_mask": occupied, "edge_mask": mask,
                                 "stationary_minus_trace_coefficient": str(actual-comparison)})
    pairs = []
    for i in range(n):
        for j in range(i+1, n):
            product = multiply(moments[1 << i], moments[1 << j])
            covariance = [a-b for a, b in zip(moments[(1 << i) | (1 << j)], product)]
            nonzero = [mask for mask in range(terms) if covariance[mask]]
            order = min((mask.bit_count() for mask in nonzero), default=None)
            pairs.append({"vertices": [i, j], "first_squarefree_order": order,
                          "leading_coefficients": [{"edge_mask": mask, "coefficient": str(covariance[mask])}
                                                   for mask in nonzero if mask.bit_count() == order]})
    return {"vertices": n, "directed_arcs_in_variable_order": arcs, "sleep_rates": sleep_rates,
            "quotient": "z_e squared equals zero for every directed edge",
            "moment_coefficients_compared": (size-1)*terms,
            "stationary_moments_equal_trace_full_moments_in_quotient": not failures,
            "disagreement_count": len(failures), "first_disagreements": failures[:12],
            "pair_covariance_leading_squarefree_terms": pairs}


def port_deletion_certificate():
    """Test cut-vertex deletion with a nontrivial, unscaled outside branch.

    Vertex 1 has the singleton branch {0} and the two-site branch {2,3}.
    Only the four directed arcs incident to 1 are nilpotent. Fixed loops
    and both transitions inside the second branch remain fully present.
    This tests the general branch argument, not just weak killing of
    every edge simultaneously.
    """
    variables = sp.symbols("a b c d")
    a, b, c, d = variables
    kernel = sp.Matrix([[sp.Rational(1, 7), a, 0, 0, 0],
                        [b, sp.Rational(1, 5), d, 0, 0],
                        [0, c, sp.Rational(1, 11), sp.Rational(1, 3), 0],
                        [0, 0, sp.Rational(1, 4), 0, 0]])
    for i in range(4):
        kernel[i, 4] = 1-sum(kernel[i, :4])
    rates = (2, 3, 5, 7)
    law = nilpotent_stationary(kernel, rates, variables)
    outside, destinations = [0, 2, 3], [0, 2, 3, 4]
    traced = (kernel.extract(outside, destinations)
              +kernel.extract(outside, [1])*kernel.extract([1], destinations)/(1-kernel[1, 1]))
    trace_law = nilpotent_stationary(traced, (rates[0], rates[2], rates[3]), variables)
    compared = 0
    for coefficient, traced_coefficient in zip(law, trace_law):
        for mask in range(8):
            original = sum(((mask >> k) & 1) << v for k, v in enumerate(outside))
            assert coefficient[original]+coefficient[original | 2] == traced_coefficient[mask]
            compared += 1
    return {"deleted_vertex": 1, "outside_branches": [[0], [2, 3]],
            "sleep_rates": rates, "variables": list(map(str, variables)),
            "transition_matrix_including_sink": [[str(value) for value in row] for row in kernel.tolist()],
            "outside_law_coefficients_compared": compared,
            "outside_marginal_equals_stationary_trace_modulo_port_squares": True}


def symbolic_kernel_certificate():
    """Check rational symbolic-kernel code against independent raw absorption.

    The three-site path distinguishes the distance-four leading term from
    the already proved second-order formula for adjacent vertices.
    """
    from algebra import stationary
    from exact import stationary_law

    parameter = sp.Symbol("epsilon")
    kernel = (((1, parameter), (3, 1-parameter)),
              ((0, parameter/2), (2, parameter/2), (3, 1-parameter)),
              ((1, parameter), (3, 1-parameter)))
    rates = (2, 3, 5)
    law = stationary(kernel, tuple(1/sp.Rational(rate) for rate in rates))
    fixed = tuple(tuple((j, probability.subs(parameter, sp.Rational(1, 2)))
                        for j, probability in row) for row in kernel)
    raw, counts = stationary_law(fixed, tuple(sp.Rational(rate, rate+1) for rate in rates))
    assert tuple(value.subs(parameter, sp.Rational(1, 2)) for value in law) == raw
    covariance = (sum(law[state] for state in range(8) if state & 5 == 5)
                  -sum(law[state] for state in range(8) if state & 1)
                  *sum(law[state] for state in range(8) if state & 4))
    leading = sp.series(covariance, parameter, 0, 5).removeO()
    assert leading == -sp.Rational(5, 648)*parameter**4
    return {"case": "path3_uniform_killing", "sleep_rates": rates,
            "internal_kernel": [["0", "epsilon", "0"],
                                ["epsilon/2", "0", "epsilon/2"],
                                ["0", "epsilon", "0"]],
            "raw_comparison_epsilon": "1/2", "raw_transient_state_counts": counts,
            "raw_both_toppling_orders_and_every_addition_checked": True,
            "symbolic_basis_residual_and_normalization_checked": True,
            "endpoint_covariance_through_order_four": str(leading)}


def nonlinear_event_certificate():
    """Check the influence factors for separated, genuinely multi-site events."""
    n = 5
    variables = sp.symbols("z:8")
    arcs = tuple(edge for i in range(4) for edge in ((i, i+1), (i+1, i)))
    transition = sp.zeros(n, n+1)
    for variable, (i, j) in zip(variables, arcs):
        transition[i, j] = variable
    for i in range(n):
        transition[i, n] = 1-sum(transition[i, :n])
    rates = (1, 2, 3, 5, 7)
    law = nilpotent_stationary(transition, rates, variables)
    left = [sum(value[state] for state in range(32) if state & 3) for value in law]
    right = [sum(value[state] for state in range(32) if state & 24 == 24) for value in law]
    joint = [sum(value[state] for state in range(32) if state & 3 and state & 24 == 24)
             for value in law]
    nonzero, compared = {}, 0
    for mask in range(256):
        if mask.bit_count() > 4:
            continue
        product, subset = 0, mask
        while True:
            product += left[subset]*right[mask ^ subset]
            if subset == 0:
                break
            subset = (subset-1) & mask
        covariance = joint[mask]-product
        compared += 1
        if covariance:
            nonzero[mask] = covariance
    parameters = tuple(sp.Rational(rate, rate+1) for rate in rates)
    expected = (-parameters[1]*(1-parameters[1])*parameters[3]*(1-parameters[3])
                *(1-parameters[0])*parameters[4])
    assert expected == -sp.Rational(35, 2592)
    assert nonzero == {60: expected}
    return {"vertices": n, "directed_arcs_in_variable_order": arcs,
            "sleep_rates": rates, "left_event": "X0 or X1", "right_event": "X3 and X4",
            "distance_between_essential_supports": 2,
            "coefficients_checked_through_degree_four": compared,
            "leading_edge_mask": 60, "leading_coefficient": str(expected),
            "product_law_influence_formula_matches": True}


def marked_exit_operator(kernel, sleep_rates, portal, initial_states=None, reverse=False):
    """Raw addition law weighted by exits at portal / its killing probability.

    Differentiate the finite absorption equations with respect to a mark
    on absorbing jumps from this site. This construction is independent
    of the IDLA basis and the self-loop perturbation identities.
    """
    from exact import DomainMatrix, add, stable, stable_mask, transitions

    n = len(kernel)
    probabilities = tuple(sp.Rational(rate, rate+1) for rate in sleep_rates)
    killing = dict(kernel[portal]).get(n, 0)
    assert killing > 0
    starts = ([add(stable(mask, n), portal) for mask in range(1 << n)]
              if initial_states is None else list(initial_states))
    assert all(any(value > 0 for value in state) for state in starts)
    states, index = [], {}

    def discover(state):
        if any(value > 0 for value in state) and state not in index:
            index[state] = len(states)
            states.append(state)

    for state in starts:
        discover(state)
    rows = []
    for state in states:
        row = transitions(state, kernel, probabilities, reverse=reverse)
        rows.append(row)
        for target in row:
            discover(target)
    size = len(states)
    matrix, rhs = sp.eye(size), sp.zeros(size, 1 << n)
    for i, row in enumerate(rows):
        for target, probability in row.items():
            if target in index:
                matrix[i, index[target]] -= probability
            else:
                rhs[i, stable_mask(target)] += probability
    inverse = DomainMatrix.from_Matrix(matrix).to_field().inv().to_Matrix()
    absorption = inverse*rhs
    assert matrix*absorption == rhs
    marked_rhs = sp.zeros(size, 1 << n)
    for i, state in enumerate(states):
        order = range(n-1, -1, -1) if reverse else range(n)
        active = next(j for j in order if state[j] > 0)
        if active != portal:
            continue
        jump = 1-probabilities[portal] if state[portal] == 1 else 1
        target = list(state)
        target[portal] -= 1
        target = tuple(target)
        if target in index:
            marked_rhs[i, :] = jump*absorption[index[target], :]
        else:
            marked_rhs[i, stable_mask(target)] = jump
    marked = inverse*marked_rhs
    assert matrix*marked == marked_rhs and all(value >= 0 for value in marked)
    return sp.Matrix.vstack(*(marked[index[state], :] for state in starts)).T, size


def bridge_response_certificate():
    """Check the rank-one cross response for two interacting two-site blocks."""
    responses, laws, greens, raw_records = [], [], [], []
    variable = sp.Symbol("z")
    blocks = (([[0, sp.Rational(1, 2)], [sp.Rational(1, 3), 0]], (1, 2), 1),
              ([[0, sp.Rational(1, 4)], [sp.Rational(1, 2), 0]], (3, 5), 0))
    for internal, rates, portal in blocks:
        transition = sp.zeros(2, 3)
        transition[:, :2] = sp.Matrix(internal)
        for i in range(2):
            transition[i, 2] = 1-sum(transition[i, :2])
        kernel = tuple(tuple((j, value) for j, value in enumerate(transition.row(i)) if value)
                       for i in range(2))
        green = (sp.eye(2)-sp.Matrix(internal)).inv()[portal, portal]
        transition[portal, portal] += variable
        transition[portal, 2] -= variable
        base, derivative = nilpotent_stationary(transition, rates, (variable,))
        marked, count = marked_exit_operator(kernel, rates, portal)
        assert sum(marked*base) == green
        assert derivative == green*base-marked*base
        laws.append(base)
        responses.append(derivative)
        greens.append(green)
        raw_records.append({"kernel": [[[j, str(value)] for j, value in row] for row in kernel],
                            "sleep_rates": rates, "portal": portal,
                            "Green_diagonal": str(green),
                            "stationary_self_loop_response": list(map(str, derivative)),
                            "raw_marked_exit_operator": [[str(value) for value in row] for row in marked.tolist()],
                            "raw_transient_states": count,
                            "marked_exit_response_identity_checked": True})
    a, b = sp.symbols("a b")
    transition = sp.Matrix([[0, sp.Rational(1, 2), 0, 0, 0],
                            [sp.Rational(1, 3), 0, a, 0, 0],
                            [0, b, 0, sp.Rational(1, 4), 0],
                            [0, 0, sp.Rational(1, 2), 0, 0]])
    for i in range(4):
        transition[i, 4] = 1-sum(transition[i, :4])
    law = nilpotent_stationary(transition, (1, 2, 3, 5), (a, b))
    assert law[1] == law[2] == sp.zeros(16, 1)
    delta = sp.Matrix(4, 4, lambda i, j: law[3][i+4*j])
    expected = (greens[1]*responses[0]*laws[1].T
                +greens[0]*laws[0]*responses[1].T-responses[0]*responses[1].T)
    assert delta == expected
    connected = (delta-(delta*sp.ones(4, 1))*laws[1].T
                 -laws[0]*(delta.T*sp.ones(4, 1)).T)
    assert connected == -responses[0]*responses[1].T and connected.rank() == 1
    return {"blocks": raw_records, "bridge_arcs": [[1, 2], [2, 1]],
            "mixed_law_response_formula_checked": True,
            "connected_response_is_negative_outer_product": True,
            "connected_response_rank": 1}


def marked_source_certificate():
    """Reject two proposed positive mixtures of all-active input laws.

    A cyclic block distinguishes this claim from the positive three-site
    saturated coefficients. F marks saturated-stabilization departures;
    D nu marks departures during an addition to a stationary input. Neither
    need have nonnegative coordinates in the all-active stabilization basis.
    """
    from algebra import basis_matrix, idla_operators
    from bk import upsets
    from exact import graph_kernel, stationary_law

    neighbors = ((1, 4), (0, 2, 3, 4), (1, 3), (1, 2))
    n, portal, rates = 4, 1, (1,)*4
    kernel = graph_kernel(neighbors)
    operators, jumps = idla_operators(kernel)
    matrix = basis_matrix(operators, jumps, rates)
    stabilization = matrix.inv()
    addition = stabilization*operators[portal]*matrix
    law = stabilization[:, -1]
    raw_law, raw_counts = stationary_law(kernel, (sp.Rational(1, 2),)*n)
    assert law == sp.Matrix(raw_law)
    assert all(value >= 0 for value in addition)
    assert addition*law == law
    saturated, count = marked_exit_operator(kernel, rates, portal, [(1,)*n])
    reverse, reverse_count = marked_exit_operator(kernel, rates, portal, [(1,)*n], reverse=True)
    assert reverse == saturated
    marked, marked_count = marked_exit_operator(kernel, rates, portal)
    variable = sp.Symbol("z")
    transition = sp.zeros(n, n+1)
    for i, row in enumerate(kernel):
        for j, value in row:
            transition[i, j] = value
    green = (sp.eye(n)-transition[:, :n]).inv()[portal, portal]
    transition[portal, portal] += variable
    transition[portal, n] -= variable
    base, response = nilpotent_stationary(transition, rates, (variable,))
    assert base == law
    assert response == green*law-marked*law
    assert response == (addition-sp.eye(1 << n))*saturated
    saturated_weights = matrix*saturated
    stationary_weights = matrix*marked*law
    assert saturated_weights[2] == -sp.Rational(3202506376937, 822513236688100)
    assert saturated_weights[-1] == 0
    responses = [sum(response[state] for state in range(1 << n) if event >> state & 1)
                 for event in upsets(n)[1:-1]]
    assert all(value > 0 for value in responses)
    return {"neighbors": neighbors, "sleep_rates": rates, "portal": portal,
            "raw_stationary_law_transient_states": raw_counts,
            "raw_saturated_marked_transient_states": [count, reverse_count],
            "raw_addition_marked_transient_states": marked_count,
            "both_toppling_orders_agree_for_saturated_marked_law": True,
            "both_response_identities_checked": True,
            "saturated_marked_law": list(map(str, saturated)),
            "saturated_all_active_basis_coefficients": list(map(str, saturated_weights)),
            "stationary_marked_all_active_basis_coefficients": list(map(str, stationary_weights)),
            "negative_saturated_coefficient_masks": [i for i, value in enumerate(saturated_weights) if value < 0],
            "negative_stationary_coefficient_masks": [i for i, value in enumerate(stationary_weights) if value < 0],
            "increasing_event_responses_checked": len(responses),
            "all_event_responses_strictly_positive": True}


def component_cycle_certificate():
    """Check cycle and chain coefficients with an interacting middle block.

    The two-site block stays unscaled. A directed three-block cycle tests
    the full marked-law formula; two cycles meeting in that block test
    the squared Green factor and propagation of a remote loop response.
    """
    a, b, c, d, z = sp.symbols("a b c d z")
    rates = (2, 3, 7, 5)
    transition = sp.Matrix([
        [z, a, 0, 0, 1-a-z],
        [b, 0, sp.Rational(1, 2), c, sp.Rational(1, 2)-b-c],
        [0, sp.Rational(1, 3), 0, 0, sp.Rational(2, 3)],
        [0, d, 0, 0, 1-d]])
    law = nilpotent_stationary(transition, rates, (a, b, c, d, z))
    first = [sum(vector[state] for state in range(16) if state & 1) for vector in law]
    last = [sum(vector[state] for state in range(16) if state & 8) for vector in law]
    covariance = []
    for mask in range(16):
        product, subset = 0, mask
        while True:
            product += first[subset]*last[mask ^ subset]
            if subset == 0:
                break
            subset = (subset-1) & mask
        covariance.append(sum(law[mask][state] for state in range(16)
                              if state & 1 and state & 8)-product)
    green = sp.Rational(6, 5)
    first_response, last_response = sp.Rational(2, 9), sp.Rational(5, 36)
    assert covariance[-1] == -first_response*last_response*green**2 == -sp.Rational(2, 45)
    assert all(value == 0 for value in covariance[:-1])
    assert first[15] == first_response*green**2 == sp.Rational(8, 25)
    assert last[15] == last_response*green**2 == sp.Rational(1, 5)
    assert all(last[mask | 16] == 0 for mask in range(15))
    assert last[31] == last_response*green**2 == sp.Rational(1, 5)

    transition[0, 0] = 0
    transition[0, 4] = 1-a
    transition[1, 0] = 0
    transition[1, 4] = sp.Rational(1, 2)-c
    transition[3, 1] = 0
    transition[3, 0] = d
    cycle_law = nilpotent_stationary(transition, rates, (a, c, d))
    block = sp.Matrix([[z, sp.Rational(1, 2), sp.Rational(1, 2)-z],
                       [sp.Rational(1, 3), 0, sp.Rational(2, 3)]])
    nu, response = nilpotent_stationary(block, rates[1:3], (z,))
    marked = green*nu-response
    initial = sp.Matrix([sp.Rational(1, 3), sp.Rational(2, 3)])
    final = sp.Matrix([sp.Rational(1, 6), sp.Rational(5, 6)])
    initial_marked = initial-sp.Matrix([-first_response, first_response])
    final_marked = final-sp.Matrix([-last_response, last_response])
    expected = sp.Matrix([
        green*initial[state & 1]*nu[(state >> 1) & 3]*final[(state >> 3) & 1]
        -initial_marked[state & 1]*marked[(state >> 1) & 3]*final_marked[(state >> 3) & 1]
        for state in range(16)])
    assert cycle_law[-1] == expected
    assert all(vector == sp.zeros(16, 1) for vector in cycle_law[1:-1])
    return {"blocks": [[0], [1, 2], [3]], "portals": [0, 1, 3], "sleep_rates": rates,
            "unscaled_internal_probabilities": {"1->2": "1/2", "2->1": "1/3"},
            "middle_Green_diagonal": str(green),
            "chain_weak_arcs_in_variable_order": [[0, 1], [1, 0], [1, 3], [3, 1]],
            "chain_covariance_coefficient": str(covariance[-1]),
            "chain_lower_squarefree_covariance_coefficients_zero": True,
            "chain_marginal_coefficients": [str(first[15]), str(last[15])],
            "remote_self_loop_response_coefficient": str(last[31]),
            "remote_lower_squarefree_response_coefficients_zero": True,
            "cycle_weak_arcs_in_variable_order": [[0, 1], [1, 3], [3, 0]],
            "cycle_mixed_law_coefficient": list(map(str, expected)),
            "full_cycle_mixed_law_formula_checked": True}


def self_loop_response_order_certificate():
    """Test the open stochastic sign of the self-loop response at specified points.

    The conditional-correlation examples test whether that obstruction
    also breaks the proposed response order. These are finite checks,
    not a proof at arbitrary rates or for arbitrary graphs.
    """
    from bk import minimal_witnesses, upsets
    from exact import graph_kernel

    variable = sp.Symbol("z")
    examples = (("conditional3", ((1, 3), (0, 2, 3), (1,)), (64,)*3, (0, 1)),
                ("endpoint4", ((1, 4), (0, 2), (1, 3), (2,)), (1,)*4, (0,)),
                ("conditional_tree5", ((1, 5), (0, 2, 3), (1,), (1, 4), (3,)), (64,)*5, (0,)))
    records = []
    for name, neighbors, rates, portals in examples:
        n = len(neighbors)
        base = sp.zeros(n, n+1)
        for i, row in enumerate(graph_kernel(neighbors)):
            for j, probability in row:
                base[i, j] = probability
        for portal in portals:
            transition = base.copy()
            transition[portal, portal] += variable
            transition[portal, n] -= variable
            derivative = nilpotent_stationary(transition, rates, (variable,))[1]
            values = [(sum(derivative[state] for state in range(1 << n) if bits >> state & 1), bits)
                      for bits in upsets(n) if bits not in (0, (1 << (1 << n))-1)]
            minimum, event = min(values)
            assert minimum > 0
            records.append({"case": name, "neighbors": neighbors, "sleep_rates": rates,
                            "portal": portal, "increasing_events_checked": len(values),
                            "minimum_event_derivative": str(minimum),
                            "minimum_event_witnesses": minimal_witnesses(event, n),
                            "all_nonconstant_increasing_event_derivatives_positive": True})
    return {"scope": "specified graph/portal/rate choices only", "cases": records}


def symbolic_response_order_certificate():
    """All common rates on the three-site stationary Rayleigh counterexample."""
    from algebra import stationary
    from exact import graph_kernel, na_gaps

    variable, inverse_rate = sp.symbols("z t")
    neighbors = ((1, 3), (0, 2, 3), (1,))
    base = graph_kernel(neighbors)
    records, curvature_obstruction = [], None
    from bk import minimal_witnesses, upsets
    for portal in (0, 1):
        kernel = []
        for i, row in enumerate(base):
            probabilities = dict(row)
            if i == portal:
                probabilities[i] = probabilities.get(i, 0)+variable
                probabilities[3] -= variable
            kernel.append(tuple(probabilities.items()))
        law = stationary(tuple(kernel), (inverse_rate,)*3)
        original = tuple(value.subs(variable, 0) for value in law)
        for record in na_gaps(original):
            numerator, denominator = sp.fraction(sp.factor(record[-1]))
            if sp.Poly(denominator, inverse_rate).LC() < 0:
                numerator, denominator = -numerator, -denominator
            assert numerator != 0
            assert all(value >= 0 for value in sp.Poly(numerator, inverse_rate).all_coeffs())
            assert all(value >= 0 for value in sp.Poly(denominator, inverse_rate).all_coeffs())
        response = [sp.factor(sp.diff(value, variable).subs(variable, 0)) for value in law]
        if portal == 0:
            event = sum(law[state] for state in range(8) if state & 1 or state & 6 == 6)
            second = sp.factor(sp.diff(event, variable, 2).subs(variable, 0)/2)
            high_rate = sp.limit(second/inverse_rate**2, inverse_rate, 0)
            fixed_rate = second.subs(inverse_rate, sp.Rational(1, 64))
            assert high_rate == -sp.Rational(149, 1500) and fixed_rate < 0
            curvature_obstruction = {"portal": 0, "event": "X0 or (X1 and X2)",
                                     "self_loop_coefficient_order": 2,
                                     "leading_high_rate_order": 2,
                                     "leading_high_rate_coefficient": str(high_rate),
                                     "coefficient_at_lambda_64": str(fixed_rate)}
        denominator = sp.lcm([sp.denom(value) for value in response])
        weights = [sp.Poly(sp.cancel(value*denominator), inverse_rate) for value in response]
        denominator = sp.Poly(denominator, inverse_rate)
        if denominator.LC() < 0:
            denominator = -denominator
            weights = [-value for value in weights]
        assert all(value >= 0 for value in denominator.all_coeffs())
        checked = 0
        for event in upsets(3):
            if event in (0, 255):
                continue
            value = sum((weights[state] for state in range(8) if event >> state & 1),
                        sp.Poly(0, inverse_rate))
            assert not value.is_zero
            assert all(coefficient >= 0 for coefficient in value.all_coeffs()), minimal_witnesses(event, 3)
            checked += 1
        records.append({"portal": portal, "increasing_events": checked,
                        "denominator_degree": denominator.degree(),
                        "denominator_coefficients_descending": list(map(str, denominator.all_coeffs())),
                        "response_numerator_coefficients_descending": [list(map(str, value.all_coeffs())) for value in weights],
                        "all_common_rates_strict_response_order_certified": True,
                        "base_law_strict_NA_at_all_common_rates_checked": True})
    return {"neighbors": neighbors, "variable": "t=1/lambda", "cases": records,
            "higher_response_coefficient_counterexample": curvature_obstruction}


def two_site_response_certificate():
    """Strict response order for every irreducible two-site kernel and rates.

    Remove fixed self-loops first. At the perturbed vertex 0, write a for
    the other-site jump probability (0<a<1), b for the reverse probability
    (0<b<=1), and p,r for sleep probabilities. The loop derivative is
    p(1-p) partial_p + a partial_a. Positive odds turn its four increasing
    event derivatives into polynomials with strictly positive coefficients.
    """
    p, r, a, b = sp.symbols("p r a b")
    q, s = 1-p, 1-r
    normalizer, harmonic = 1-a*b+a*b*p*r, 1-a*b*q*s
    both = p*r/normalizer
    first = p*s*(1-a*b+a*b*r*(1-b)*q)/(normalizer*harmonic)
    second = r*q*(1-a*b+a*b*p*(1-a)*s)/(normalizer*harmonic)
    x, y, aa, bb = sp.symbols("x y A B", positive=True)
    records = []
    for name, probability in (("X0", first+both), ("X1", second+both),
                              ("AND", both), ("OR", first+second+both)):
        derivative = sp.factor(p*q*sp.diff(probability, p)+a*sp.diff(probability, a))
        cases = []
        for boundary in (False, True):
            substituted = derivative.subs(b, 1) if boundary else derivative
            changes = {p: x/(1+x), r: y/(1+y), a: aa/(1+aa)}
            if not boundary:
                changes[b] = bb/(1+bb)
            numerator, denominator = sp.fraction(sp.factor(substituted.subs(changes)))
            variables = (x, y, aa) if boundary else (x, y, aa, bb)
            numerator, denominator = (sp.Poly(value, *variables) for value in (numerator, denominator))
            assert all(value > 0 for value in numerator.coeffs())
            assert all(value > 0 for value in denominator.coeffs())
            cases.append({"reverse_probability": "b=1" if boundary else "b=B/(1+B)",
                          "positive_numerator_terms": len(numerator.terms()),
                          "positive_denominator_terms": len(denominator.terms())})
        records.append({"event": name, "cases": cases})
    return {"scope": "0<a<1, 0<b<=1, all positive site rates; fixed loops can be erased",
            "substitution": "p=x/(1+x), r=y/(1+y), a=A/(1+A)",
            "all_four_increasing_event_derivatives_strictly_positive": True,
            "cases": records}


def reinjection_flow_certificate():
    """Validate the all-rate rank flow under uniform sink reinjection.

    The directed, looped, unequal-rate example tests the scope beyond
    reversible graphs and the factor 1/(1-z). The endpoint path tests
    the sole-killed-port self-loop specialization. Exact differentiation
    of the IDLA algebra is compared with independent raw addition laws,
    both raw toppling orders, and every addition-stationarity equation.
    """
    from exact import absorption, add, graph_kernel, stable

    rational = sp.Rational
    directed = (
        ((0, rational(1, 6)), (1, rational(1, 3)),
         (2, rational(1, 6)), (3, rational(1, 3))),
        ((0, rational(1, 4)), (2, rational(1, 4)), (3, rational(1, 2))),
        ((0, rational(1, 5)), (1, rational(2, 5)),
         (2, rational(1, 10)), (3, rational(3, 10))))
    examples = (
        ("directed_looped3", directed, (rational(1, 2), 2, 3),
         (rational(1, 2), rational(1, 3), rational(1, 6)),
         (sp.S.Zero, rational(2, 5))),
        ("endpoint_path3", graph_kernel(((1, 3), (0, 2), (1,))),
         (1, 1, 1), (sp.S.One, sp.S.Zero, sp.S.Zero),
         (sp.S.Zero, rational(1, 3))))
    variable = sp.Symbol("z")
    records = []
    for name, original, rates, drive, levels in examples:
        n, size = len(original), 1 << len(original)
        base = sp.zeros(n, n+1)
        for i, row in enumerate(original):
            for j, value in row:
                base[i, j] = value
        direction = sp.zeros(n, n+1)
        for i in range(n):
            for j in range(n):
                direction[i, j] = base[i, n]*drive[j]
            direction[i, n] = -base[i, n]
        sleep = tuple(rational(rate, rate+1) for rate in rates)
        original_law = original_additions = original_addition = None
        for level in levels:
            fixed = base+level*direction
            law, derivative = nilpotent_stationary(
                fixed+variable*direction, rates, (variable,))
            kernel = tuple(tuple((j, fixed[i, j]) for j in range(n+1)
                                 if fixed[i, j]) for i in range(n))
            starts = [(1,)*n]+[
                add(stable(mask, n), vertex)
                for vertex in range(n) for mask in range(size)]
            raw, raw_count = absorption(starts, kernel, sleep)
            reverse, reverse_count = absorption(
                [(1,)*n], kernel, sleep, reverse=True)
            assert raw[0, :].T == law == reverse.T
            additions = []
            for vertex in range(n):
                operator = raw[1+vertex*size:1+(vertex+1)*size, :].T
                assert operator*law == law
                additions.append(operator)
            addition = sum((drive[i]*additions[i] for i in range(n)),
                           sp.zeros(size))
            if level == 0:
                original_law = law
                original_additions = additions
                original_addition = addition
            assert original_law is not None
            # Finite-parameter feedback transform, independently using the
            # raw additions at both kernels, not just its derivative at zero.
            feedback_basis_step = (sp.eye(size)-level*original_addition)/(1-level)
            feedback_step = (1-level)*sp.eye(size)+level*addition
            basis_powers, inverse_powers = [sp.eye(size)], [sp.eye(size)]
            for _ in range(n):
                basis_powers.append(basis_powers[-1]*feedback_basis_step)
                inverse_powers.append(inverse_powers[-1]*feedback_step)
            feedback_basis = sp.Matrix.hstack(*[
                basis_powers[n-mask.bit_count()][:, mask] for mask in range(size)])
            feedback_inverse = sp.Matrix.hstack(*[
                inverse_powers[n-mask.bit_count()][:, mask] for mask in range(size)])
            assert feedback_basis*law == original_law
            assert feedback_basis*feedback_inverse == sp.eye(size)
            assert all(value >= 0 for value in feedback_inverse)
            resolvent = (sp.eye(size)-level*original_addition).inv()
            for before, after in zip(original_additions, additions):
                assert (feedback_basis*after
                        == (1-level)*before*resolvent*feedback_basis)
            marked_rank = sp.Matrix([
                mask.bit_count()*law[mask] for mask in range(size)])
            assert derivative == (marked_rank-addition*marked_rank)/(1-level)
            assert sum(derivative) == 0
            tails, crossings = [], []
            for threshold in range(1, n+1):
                tail = sum(derivative[mask] for mask in range(size)
                           if mask.bit_count() >= threshold)
                center = rational(2*threshold-1, 2)
                crossing = sp.S.Zero
                for before in range(size):
                    for after in range(size):
                        term = ((before.bit_count()-center)
                                *(int(before.bit_count() >= threshold)
                                  -int(after.bit_count() >= threshold)))
                        assert term >= 0
                        crossing += law[before]*addition[after, before]*term
                assert tail == crossing/(1-level) and tail > 0
                tails.append(str(tail))
                crossings.append(str(crossing))
            mean_derivative = sum(mask.bit_count()*derivative[mask]
                                  for mask in range(size))
            squared_increment = sum(
                law[before]*addition[after, before]
                *(after.bit_count()-before.bit_count())**2
                for before in range(size) for after in range(size))
            assert mean_derivative == squared_increment/(2*(1-level)) > 0
            records.append({
                "case": name, "original_kernel": [
                    [[j, str(value)] for j, value in row] for row in original],
                "sleep_rates": list(map(str, rates)),
                "reinjection_distribution": list(map(str, drive)),
                "reinjection_probability": str(level),
                "stationary_law": list(map(str, law)),
                "stationary_derivative": list(map(str, derivative)),
                "rank_tail_derivatives": tails,
                "nonnegative_crossing_expectations": crossings,
                "mean_rank_derivative": str(mean_derivative),
                "stationary_squared_rank_increment": str(squared_increment),
                "raw_transient_state_counts": [raw_count, reverse_count],
                "both_toppling_orders_and_every_addition_checked": True,
                "full_vector_flow_and_rank_crossing_identities_checked": True,
                "finite_feedback_basis_law_and_stochastic_inverse_checked": True,
                "finite_feedback_intertwining_checked_for_every_addition": True})
    return {
        "scope": "exact checks supporting report (50)-(52), (56)-(57)",
        "cases": records}


def distributed_bridge_certificate():
    """Check total-loss routing between two interacting components.

    Both killed vertices in each two-site block can send a particle to
    a nonuniform entrance distribution in the other block. Raw loss
    conservation gives the responses; nilpotent IDLA algebra independently
    checks the mixed law. A self-reinjection variable in the second block
    tests the nonconstant port normalization g(s)=1/(1-s).
    """
    from exact import absorption, add, stable

    rational = sp.Rational
    specs = (
        (sp.Matrix([[0, rational(1, 2)], [rational(1, 3), 0]]),
         (1, 2), (rational(2, 3), rational(1, 3))),
        (sp.Matrix([[0, rational(1, 4)], [rational(1, 2), 0]]),
         (3, 5), (rational(1, 4), rational(3, 4))))
    laws, responses, records = [], [], []
    for internal, rates, drive in specs:
        transition = internal.row_join(sp.ones(2, 1)-internal*sp.ones(2, 1))
        kernel = tuple(tuple((j, value) for j, value in enumerate(transition.row(i))
                             if value) for i in range(2))
        sleep = tuple(rational(rate, rate+1) for rate in rates)
        starts = [(1, 1)]+[add(stable(mask, 2), v)
                           for v in range(2) for mask in range(4)]
        raw, count = absorption(starts, kernel, sleep)
        reverse, reverse_count = absorption([(1, 1)], kernel, sleep, reverse=True)
        nu = raw[0, :].T
        assert reverse.T == nu
        additions = [raw[1+4*v:1+4*(v+1), :].T for v in range(2)]
        assert all(operator*nu == nu for operator in additions)
        addition = sum((drive[v]*additions[v] for v in range(2)), sp.zeros(4))
        loss = sp.Matrix(4, 4, lambda after, before:
                         (1+int(before).bit_count()-int(after).bit_count())*addition[after, before])
        assert all(value >= 0 for value in loss) and sum(loss*nu) == 1
        response = nu-loss*nu
        marked_rank = sp.Matrix([mask.bit_count()*nu[mask] for mask in range(4)])
        assert response == marked_rank-addition*marked_rank
        tails = [sum(response[mask] for mask in range(4) if mask.bit_count() >= r)
                 for r in (1, 2)]
        assert min(tails) > 0
        mean_response = sum(mask.bit_count()*response[mask] for mask in range(4))
        loss_variance = sum(
            nu[before]*addition[after, before]*(before.bit_count()-after.bit_count())**2
            for before in range(4) for after in range(4))
        assert mean_response == loss_variance/2
        laws.append(nu)
        responses.append(response)
        records.append({
            "internal_kernel": [[str(value) for value in row] for row in internal.tolist()],
            "sleep_rates": rates, "entrance_distribution": list(map(str, drive)),
            "stationary_law": list(map(str, nu)),
            "total_loss_response": list(map(str, response)),
            "rank_tail_responses": list(map(str, tails)),
            "mean_rank_response": str(mean_response),
            "stationary_avalanche_loss_variance": str(loss_variance),
            "raw_transient_state_counts": [count, reverse_count],
            "both_toppling_orders_every_addition_and_loss_conservation_checked": True})
    a, b, s = sp.symbols("a b s")
    transition = sp.zeros(4, 5)
    for block, (internal, _, drive) in enumerate(specs):
        offset = 2*block
        other_offset = 2*(1-block)
        coupling = a if block == 0 else b
        other_drive = specs[1-block][2]
        for i in range(2):
            killing = 1-sum(internal[i, :])
            for j in range(2):
                transition[offset+i, offset+j] = internal[i, j]
                transition[offset+i, other_offset+j] = coupling*killing*other_drive[j]
                if block == 1:
                    transition[offset+i, offset+j] += s*killing*drive[j]
            transition[offset+i, 4] = killing*(1-coupling-(s if block == 1 else 0))
    law = nilpotent_stationary(transition, (1, 2, 3, 5), (a, b, s))
    product = sp.Matrix([laws[0][mask & 3]*laws[1][mask >> 2] for mask in range(16)])
    assert law[0] == product
    for mask in (1, 2, 5, 6):
        assert law[mask] == sp.zeros(16, 1)
    delta = sp.Matrix(4, 4, lambda i, j: law[3][i+4*j])
    expected = (responses[0]*laws[1].T+laws[0]*responses[1].T
                -responses[0]*responses[1].T)
    assert delta == expected
    connected = delta-(delta*sp.ones(4, 1))*laws[1].T-laws[0]*(delta.T*sp.ones(4, 1)).T
    assert connected == -responses[0]*responses[1].T
    count = sp.Matrix([mask.bit_count() for mask in range(4)])
    covariance = (count.T*connected*count)[0]
    assert covariance == -sp.prod(sp.Rational(record["stationary_avalanche_loss_variance"])
                                  for record in records)/4 < 0
    remote = sum((mask & 3).bit_count()*law[7][mask] for mask in range(16))
    assert remote == (count.T*responses[0])[0] > 0
    return {
        "scope": "exact support for distributed-port cycle law and block-count covariance",
        "blocks": records,
        "mixed_covariance_coefficient": str(covariance),
        "mixed_covariance_equals_minus_quarter_loss_variance_product": True,
        "full_mixed_law_and_negative_outer_product_checked": True,
        "remote_self_reinjection_marginal_coefficient": str(remote),
        "remote_coefficient_checks_port_g_prime_equals_one_not_zero": True}


def feedback_audit_certificate():
    """Check the finite feedback cocycle, generator, and marking convention.

    A reducible kernel with a zero killing row and a non-reaching drive
    tests the limiting scope beyond the harmonic eigenbasis proof. A
    direct-killing example tests harmonic scaling with unequal rates.
    The one-site witness distinguishes the sink-marked terminal measure
    from the original-attempt-marked measure despite their equal masses.
    These are targeted exact checks, not a general proof.
    """
    from exact import absorption, add, stable

    rational = sp.Rational

    def raw(internal, rates, drive, level):
        n, size = internal.rows, 1 << internal.rows
        killing = sp.ones(n, 1)-internal*sp.ones(n, 1)
        inside = internal+level*killing*drive.T
        transition = inside.row_join((1-level)*killing)
        kernel = tuple(tuple((j, transition[i, j])
                             for j in range(n+1) if transition[i, j])
                       for i in range(n))
        starts = [(1,)*n]+[
            add(stable(mask, n), v) for v in range(n) for mask in range(size)]
        sleep = tuple(rational(rate, rate+1) for rate in rates)
        result, count = absorption(starts, kernel, sleep)
        reverse, reverse_count = absorption([(1,)*n], kernel, sleep, reverse=True)
        law = result[0, :].T
        assert law == reverse.T
        individual = [result[1+v*size:1+(v+1)*size, :].T for v in range(n)]
        assert all(operator*law == law for operator in individual)
        addition = sum((drive[v]*individual[v] for v in range(n)), sp.zeros(size))
        return law, individual, addition, kernel, (count, reverse_count)

    def columns(matrix, n):
        return sp.Matrix.hstack(*[
            (matrix**(n-mask.bit_count()))[:, mask] for mask in range(1 << n)])

    def harmonic(internal, rates, subset):
        n = internal.rows
        killing = sp.ones(n, 1)-internal*sp.ones(n, 1)
        active = sp.zeros(n, 1)
        if subset:
            values = ((sp.eye(len(subset))-internal.extract(subset, subset)).inv()
                      *killing.extract(subset, [0]))
            for i, v in enumerate(subset):
                active[v] = values[i]
        first_step = killing+internal*active
        sleeping = sp.Matrix([
            active[v] if v in subset else -first_step[v]/rates[v]
            for v in range(n)])
        row = sp.Matrix([[
            sp.prod(sleeping[v] for v in range(n) if mask & (1 << v))
            for mask in range(1 << n)]])
        return active, first_step, sleeping, row

    examples = (
        ("reducible_zero_killing_nonreaching_drive",
         sp.Matrix([[0, 1], [0, 0]]), (2, 3), sp.Matrix([0, 1])),
        ("irreducible_direct_killing_unequal_rates",
         sp.Matrix([[rational(1, 5), rational(2, 5)],
                    [rational(1, 3), rational(1, 6)]]),
         (2, 3), sp.Matrix([rational(2, 3), rational(1, 3)])))
    levels = (sp.S.Zero, rational(1, 4), rational(2, 3))
    records = []
    for name, internal, rates, drive in examples:
        n, size = internal.rows, 1 << internal.rows
        identity = sp.eye(size)
        killing = sp.ones(n, 1)-internal*sp.ones(n, 1)
        data = {level: raw(internal, rates, drive, level) for level in levels}
        original_law, original_individual, original_addition, _, _ = data[0]
        bases, inverses, level_records = {}, {}, []
        original_rows = [
            harmonic(internal, rates, [v for v in range(n) if mask & (1 << v)])[3]
            for mask in range(size)]
        harmonic_rank = sp.Matrix.vstack(*original_rows).rank()
        if name == "reducible_zero_killing_nonreaching_drive":
            assert harmonic_rank < size
        else:
            assert harmonic_rank == size
        for level in levels:
            law, individual, addition, _, counts = data[level]
            bases[level] = columns(
                (identity-level*original_addition)/(1-level), n)
            inverses[level] = columns((1-level)*identity+level*addition, n)
            assert bases[level]*law == original_law
            assert bases[level]*inverses[level] == inverses[level]*bases[level] == identity
            assert all(value >= 0 for value in inverses[level])
            assert sp.ones(1, size)*inverses[level] == sp.ones(1, size)
            for v in range(n):
                assert (bases[level]*individual[v]
                        == (1-level)*original_individual[v]
                        *(identity-level*original_addition).inv()*bases[level])
            current = internal+level*killing*drive.T
            for mask in range(size):
                subset = [v for v in range(n) if mask & (1 << v)]
                active, first_step, sleeping, row = harmonic(internal, rates, subset)
                active_z, step_z, sleep_z, row_z = harmonic(current, rates, subset)
                scale = (1-level)/(1-level*(drive.T*active)[0])
                assert active_z == scale*active
                assert step_z == scale*first_step and sleep_z == scale*sleeping
                assert row*bases[level] == scale**(-n)*row_z
                for v in range(n):
                    assert row*original_individual[v] == active[v]*row
                    assert row_z*individual[v] == active_z[v]*row_z
            level_records.append({
                "reinjection_probability": str(level),
                "stationary_law": list(map(str, law)),
                "raw_transient_state_counts": list(counts),
                "harmonic_subsets_checked": size})
        a, b = levels[1:]
        delta = (b-a)/(1-a)
        intermediate_basis = columns((identity-delta*data[a][2])/(1-delta), n)
        intermediate_inverse = columns((1-delta)*identity+delta*data[b][2], n)
        assert bases[a]*intermediate_basis == bases[b]
        assert intermediate_inverse*inverses[a] == inverses[b]
        assert intermediate_basis*intermediate_inverse == identity
        variable = sp.Symbol("z")
        basis = columns((identity-variable*original_addition)/(1-variable), n)
        inverse = basis.inv().applyfunc(sp.cancel)
        inverse_derivative = inverse.diff(variable)
        count_operator = sp.diag(*[mask.bit_count() for mask in range(size)])
        for level in levels:
            generator = ((data[level][2]-identity)*(n*identity-count_operator)
                         /(1-level))
            assert inverse.subs(variable, level) == inverses[level]
            assert inverse_derivative.subs(variable, level) == generator*inverses[level]
            assert all(generator[i, j] >= 0 for i in range(size)
                       for j in range(size) if i != j)
            assert sp.ones(1, size)*generator == sp.zeros(1, size)
        records.append({
            "case": name,
            "original_internal_kernel": [
                list(map(str, row)) for row in internal.tolist()],
            "original_killing": list(map(str, killing)),
            "sleep_rates": rates,
            "reinjection_distribution": list(map(str, drive)),
            "original_harmonic_row_rank": harmonic_rank,
            "stable_state_count": size,
            "levels": level_records,
            "cocycle_parameters": {
                "a": str(a), "b": str(b), "relative_reinjection": str(delta)},
            "both_raw_toppling_orders_and_every_addition_checked": True,
            "active_first_step_sleeping_and_harmonic_row_scaling_checked": True,
            "harmonic_rows_are_eigenrows_of_every_raw_addition": True,
            "finite_law_transform_and_both_inverse_identities_checked": True,
            "individual_addition_intertwining_checked": True,
            "inverse_nonnegativity_and_column_sums_checked": True,
            "basis_and_stochastic_inverse_cocycles_checked": True,
            "generator_matches_exact_inverse_derivative_with_raw_addition": True,
            "generator_offdiagonal_signs_and_column_sums_checked": True})

    # Original K=0 means every jump attempt is an original exit attempt.
    # After feedback, both self-routed jumps and sink jumps receive this mark.
    level = rational(1, 2)
    law, _, addition, kernel, counts = raw(
        sp.zeros(1), (1,), sp.Matrix([1]), level)
    assert law == sp.Matrix([rational(1, 3), rational(2, 3)])
    assert addition == sp.Matrix.hstack(law, law)
    sink_operator, marked_count = marked_exit_operator(kernel, (1,), 0)
    conservation_operator = sp.Matrix(
        2, 2, lambda after, before:
        (1+int(before)-int(after))*addition[after, before]/(1-level))
    assert sink_operator == conservation_operator
    sink_measure = sink_operator*law

    # The two transient states have one and two active particles.
    # t marks each jump, including a jump that leaves the active count unchanged.
    mark = sp.Symbol("t")
    left = sp.Matrix([[1-mark/4, 0], [-mark/2, 1-mark/2]])
    rhs = sp.Matrix([[mark/4, rational(1, 2)], [0, 0]])
    generating_law = (left.inv()*rhs).applyfunc(sp.cancel)
    assert (left*generating_law-rhs).applyfunc(sp.cancel) == sp.zeros(2)
    assert generating_law.subs(mark, 1).T == addition
    attempt_operator = generating_law.diff(mark).subs(mark, 1).T
    assert all(value >= 0 for value in attempt_operator)
    attempt_measure = attempt_operator*law
    assert sink_measure == sp.Matrix([rational(10, 9), rational(8, 9)])
    assert attempt_measure == sp.Matrix([rational(8, 9), rational(10, 9)])
    assert sum(sink_measure) == sum(attempt_measure) == 2
    assert sink_measure != attempt_measure
    count_operator = sp.diag(0, 1)
    response = (count_operator*law-addition*count_operator*law)/(1-level)
    assert response == 2*law-sink_measure
    assert response != 2*law-attempt_measure
    witness = {
        "original_internal_kernel": [["0"]],
        "sleep_rates": [1],
        "reinjection_distribution": ["1"],
        "reinjection_probability": str(level),
        "stable_state_order": ["empty", "occupied"],
        "stationary_law": list(map(str, law)),
        "raw_transient_state_counts": list(counts),
        "sink_marking_transient_states": marked_count,
        "normalized_sink_marked_operator": [
            list(map(str, row)) for row in sink_operator.tolist()],
        "original_attempt_marked_operator": [
            list(map(str, row)) for row in attempt_operator.tolist()],
        "normalized_sink_marked_terminal_measure": list(map(str, sink_measure)),
        "original_attempt_marked_terminal_measure": list(map(str, attempt_measure)),
        "common_marked_mass": str(sum(sink_measure)),
        "actual_self_reinjection_response": list(map(str, response)),
        "incorrect_response_from_attempt_marking": list(map(str, 2*law-attempt_measure)),
        "attempt_generating_law_by_initial_active_count": [
            list(map(str, row)) for row in generating_law.tolist()],
        "raw_sink_marking_equals_particle_conservation": True,
        "attempt_marking_verified_by_finite_absorption_equations": True,
        "equal_masses_but_distinct_terminal_measures_checked": True}
    return {
        "scope": "targeted exact checks of report (56)-(58) and original-exit normalization",
        "cases": records,
        "one_site_marking_witness": witness}


def multiblock_feedback_certificate():
    """Check the multiblock sampler, inverse, and sharp killed boundary.

    One interacting two-site block and one singleton test the full operator
    identities with nonuniform entrance and killing, and diagonal/cross routing.
    A killed two-singleton boundary checks that zero residual killing makes
    the sampler and its corresponding original operator singular.
    These targeted exact checks supplement the finite proofs.
    """
    from exact import absorption, add, stable

    R = sp.Rational


    def raw(internal, rates):
        n, size = internal.rows, 1 << internal.rows
        kill = sp.ones(n, 1) - internal * sp.ones(n, 1)
        transition = internal.row_join(kill)
        kernel = tuple(tuple((j, transition[i, j])
                             for j in range(n + 1) if transition[i, j])
                       for i in range(n))
        starts = [(1,) * n] + [
            add(stable(mask, n), v) for v in range(n) for mask in range(size)]
        sleep = tuple(R(rate, rate + 1) for rate in rates)
        laws, count = absorption(starts, kernel, sleep)
        reverse, reverse_count = absorption([(1,) * n], kernel, sleep, reverse=True)
        nu = laws[0, :].T
        assert nu == reverse.T
        additions = [laws[1 + v * size:1 + (v + 1) * size, :].T
                     for v in range(n)]
        assert all(T * nu == nu for T in additions)
        assert all(T * U == U * T for T in additions for U in additions)
        return nu, additions, (count, reverse_count)


    def harmonic(internal, rates, subset):
        n = internal.rows
        kill = sp.ones(n, 1) - internal * sp.ones(n, 1)
        a = sp.zeros(n, 1)
        if subset:
            values = (sp.eye(len(subset)) - internal.extract(subset, subset)).inv() * kill.extract(subset, [0])
            for j, v in enumerate(subset):
                a[v] = values[j]
        first = kill + internal * a
        b = sp.Matrix([a[v] if v in subset else -first[v] / rates[v]
                       for v in range(n)])
        row = sp.Matrix([[
            sp.prod(b[v] for v in range(n) if mask & (1 << v))
            for mask in range(1 << n)]])
        return a, first, b, row


    blocks = ((0, 1), (2,))
    internals = (sp.Matrix([[0, R(1, 2)], [R(1, 3), 0]]), sp.zeros(1))
    rates_by_block = ((1, 2), (3,))
    rates = (1, 2, 3)
    entrances = (sp.Matrix([R(2, 3), R(1, 3)]), sp.Matrix([1]))
    Z = sp.Matrix([[R(1, 5), R(1, 4)], [R(1, 3), R(1, 6)]])
    q = sp.ones(2, 1) - Z * sp.ones(2, 1)
    base = sp.diag(*internals)
    coupled = base.copy()
    for i, vertices in enumerate(blocks):
        kill = sp.ones(len(vertices), 1) - internals[i] * sp.ones(len(vertices), 1)
        for j, targets in enumerate(blocks):
            for u, x in enumerate(vertices):
                for v, y in enumerate(targets):
                    coupled[x, y] += kill[u] * Z[i, j] * entrances[j][v]

    isolated = [raw(K, local_rates)
                for K, local_rates in zip(internals, rates_by_block)]
    pi = sp.Matrix([
        isolated[0][0][mask & 3] * isolated[1][0][mask >> 2]
        for mask in range(8)])
    nu0, original_additions, base_counts = raw(base, rates)
    nuZ, coupled_additions, coupled_counts = raw(coupled, rates)
    assert nu0 == pi
    I = sp.eye(8)
    A, T = [], []
    for vertices, rho in zip(blocks, entrances):
        A.append(sum((rho[j] * original_additions[v]
                      for j, v in enumerate(vertices)), sp.zeros(8)))
        T.append(sum((rho[j] * coupled_additions[v]
                      for j, v in enumerate(vertices)), sp.zeros(8)))
    completion = [q[i] * I + sum((Z[i, j] * T[j] for j in range(2)), sp.zeros(8))
         for i in range(2)]
    assert completion[0] * completion[1] == completion[1] * completion[0]
    assert all(sp.ones(1, 8) * step == sp.ones(1, 8) for step in completion)
    assert all(value >= 0 for step in completion for value in step)
    vacancies = [(2 - (mask & 3).bit_count(), 1 - (mask >> 2).bit_count())
                 for mask in range(8)]
    C = sp.Matrix.hstack(*[
        (completion[0] ** vacancies[mask][0] * completion[1] ** vacancies[mask][1])[:, mask]
        for mask in range(8)])
    assert C * pi == nuZ
    assert all(value >= 0 for value in C)
    assert sp.ones(1, 8) * C == sp.ones(1, 8)

    # A 2-by-2 matrix whose entries are 8-by-8 commuting original operators.
    operator_matrix = sp.BlockMatrix([
        [(I if i == j else sp.zeros(8)) - Z[i, j] * A[j]
         for j in range(2)] for i in range(2)]).as_explicit()
    operator_rhs = sp.Matrix.vstack(*(q[i] * I for i in range(2)))
    operator_solution = operator_matrix.inv() * operator_rhs
    assert operator_matrix * operator_solution == operator_rhs
    Rs = [operator_solution[i * 8:(i + 1) * 8, :] for i in range(2)]
    assert Rs[0] * Rs[1] == Rs[1] * Rs[0]
    for i in range(2):
        assert Rs[i] == q[i] * I + sum(
            (Z[i, j] * A[j] * Rs[j] for j in range(2)), sp.zeros(8))
        assert all(value >= 0 for value in Rs[i])
        assert sp.ones(1, 8) * Rs[i] == sp.ones(1, 8)
        assert Rs[i].det() > 0
    Rinv = [operator.inv() for operator in Rs]
    M = sp.Matrix.hstack(*[
        (Rinv[0] ** vacancies[mask][0] * Rinv[1] ** vacancies[mask][1])[:, mask]
        for mask in range(8)])
    assert M * C == C * M == I
    assert M * nuZ == pi
    for i, vertices in enumerate(blocks):
        assert M * completion[i] == Rs[i] * M
        for v in vertices:
            assert M * coupled_additions[v] == Rs[i] * original_additions[v] * M

    harmonic_records = []
    for mask in range(8):
        subsets = [[v for v in range(len(vertices)) if mask & (1 << vertices[v])]
                   for vertices in blocks]
        local_harmonic = [harmonic(K, local_rates, subset)
                          for K, local_rates, subset in zip(internals, rates_by_block, subsets)]
        h = sp.Matrix([(rho.T * local[0])[0]
                       for rho, local in zip(entrances, local_harmonic)])
        r = (sp.eye(2) - Z * sp.diag(*h)).inv() * q
        assert all(q[i] <= r[i] <= 1 for i in range(2))
        a0, first0, b0, row0 = harmonic(
            base, rates, [v for v in range(3) if mask & (1 << v)])
        az, firstz, bz, rowz = harmonic(
            coupled, rates, [v for v in range(3) if mask & (1 << v)])
        for i, vertices in enumerate(blocks):
            for v in vertices:
                assert az[v] == r[i] * a0[v]
                assert firstz[v] == r[i] * first0[v]
                assert bz[v] == r[i] * b0[v]
            assert row0 * Rs[i] == r[i] * row0
            assert rowz * completion[i] == r[i] * rowz
        gamma = r[0] ** 2 * r[1]
        assert rowz * C == gamma * row0
        assert row0 * M == rowz / gamma
        harmonic_records.append({
            "global_subset_mask": mask,
            "base_entrance_eigenvalues": list(map(str, h)),
            "block_scales": list(map(str, r))})

    record = {
        "block_vertices": blocks,
        "original_internal_kernels": [
            [list(map(str, row)) for row in K.tolist()] for K in internals],
        "sleep_rates": rates,
        "entrances": [list(map(str, rho)) for rho in entrances],
        "routing_matrix": [list(map(str, row)) for row in Z.tolist()],
        "residual_killing_fractions": list(map(str, q)),
        "coupled_internal_kernel": [list(map(str, row)) for row in coupled.tolist()],
        "base_raw_transient_counts": base_counts,
        "coupled_raw_transient_counts": coupled_counts,
        "product_base_law": list(map(str, pi)),
        "coupled_stationary_law": list(map(str, nuZ)),
        "sampler_determinant": str(C.det()),
        "original_operator_determinants": [str(operator.det()) for operator in Rs],
        "harmonic_scales": harmonic_records,
        "both_raw_toppling_orders_and_all_addition_equations_checked": True,
        "sampler_stochasticity_and_full_law_identity_checked": True,
        "both_inverse_identities_and_law_inverse_checked": True,
        "individual_addition_and_completion_intertwinings_checked": True,
        "all_eight_harmonic_scalings_and_rows_checked": True}
    # The killed boundary routing has q_1=0. Its sampler still works,
    # but the all-empty harmonic row is a nonzero left null vector.
    boundary_base = sp.zeros(2)
    boundary_kernel = sp.Matrix([[0, 1], [0, 0]])
    boundary_pi, boundary_base_additions, boundary_base_counts = raw(
        boundary_base, (1, 1))
    boundary_law, boundary_additions, boundary_counts = raw(
        boundary_kernel, (1, 1))
    assert boundary_pi == boundary_law == sp.ones(4, 1)/4
    boundary_completion = boundary_additions[1]
    boundary_sampler = sp.Matrix.hstack(*[
        (boundary_completion**(1-(mask & 1)))[:, mask]
        for mask in range(4)])
    expected_boundary = sp.Matrix([
        [R(1, 2), 0, R(1, 2), 0],
        [0, 1, 0, 0],
        [R(1, 2), 0, R(1, 2), 0],
        [0, 0, 0, 1]])
    assert boundary_sampler == expected_boundary
    assert boundary_sampler*boundary_pi == boundary_law
    assert sp.ones(1, 4)*boundary_sampler == sp.ones(1, 4)
    assert all(value >= 0 for value in boundary_sampler)
    assert boundary_sampler.det() == 0 and boundary_sampler.rank() == 3
    boundary_null = sp.Matrix([[1, 0, -1, 0]])
    assert boundary_null*boundary_sampler == sp.zeros(1, 4)
    boundary_original_R = boundary_base_additions[1]
    assert boundary_original_R.det() == 0
    assert sp.Matrix([[1, -1, -1, 1]])*boundary_original_R == sp.zeros(1, 4)
    record['boundary_singularity'] = {
        'original_internal_kernels': [[['0']], [['0']]],
        'sleep_rates': [1, 1],
        'entrances': [['1'], ['1']],
        'routing_matrix': [['0', '1'], ['0', '0']],
        'residual_killing_fractions': ['0', '1'],
        'stable_state_order': ['00', '10', '01', '11'],
        'base_raw_transient_counts': boundary_base_counts,
        'coupled_raw_transient_counts': boundary_counts,
        'product_base_law': list(map(str, boundary_pi)),
        'coupled_stationary_law': list(map(str, boundary_law)),
        'sampler_matrix': [list(map(str, row)) for row in boundary_sampler.tolist()],
        'sampler_rank': boundary_sampler.rank(),
        'sampler_determinant': str(boundary_sampler.det()),
        'coupled_empty_harmonic_left_null_row': list(map(str, boundary_null)),
        'first_original_operator_rank': boundary_original_R.rank(),
        'first_original_operator_determinant': str(boundary_original_R.det()),
        'raw_sampler_identity_and_stochasticity_checked': True,
        'nonzero_left_null_rows_and_singular_operators_checked': True}
    return record


def two_block_full_occupation_certificate():
    """Test the round-trip formula without identifying traces with marginals.

    Symbolic routing verifies every Green minor on an interacting block plus
    a singleton. Two routings with equal product independently check equal
    full occupation and unequal actual marginals by raw stabilization.
    """
    from exact import absorption, add, stable

    R = sp.Rational
    a, b = sp.symbols("a b", nonnegative=True)
    base = sp.Matrix([[0, R(1, 2), 0], [R(1, 3), 0, 0], [0, 0, 0]])
    coupled = sp.Matrix([
        [0, R(1, 2), a/2], [R(1, 3), 0, 2*a/3], [2*b/3, b/3, 0]])
    rates = (1, 2, 3)
    entrances = (sp.Matrix([R(2, 3), R(1, 3)]), sp.Matrix([1]))
    blocks = ((0, 1), (2,))
    base_green = (sp.eye(3)-base).inv()
    green = (sp.eye(3)-coupled).inv().applyfunc(sp.factor)
    minor_records, inverse_full = [], sp.S.Zero
    for mask in range(8):
        selected = [v for v in range(3) if mask & (1 << v)]
        h = []
        block_minors = []
        for vertices, rho in zip(blocks, entrances):
            remaining = [v for v in vertices if v not in selected]
            within = [v for v in vertices if v in selected]
            block_minors.append(base_green.extract(within, within).det())
            kill = (sp.ones(3, 1)-base*sp.ones(3, 1)).extract(remaining, [0])
            rho_remaining = rho.extract(
                [vertices.index(v) for v in remaining], [0])
            h.append((rho_remaining.T*(sp.eye(len(remaining))
                        -base.extract(remaining, remaining)).inv()*kill)[0]
                     if remaining else sp.S.Zero)
        actual = sp.factor(green.extract(selected, selected).det())
        predicted = sp.prod(block_minors)*(1-a*b*h[0]*h[1])/(1-a*b)
        assert sp.cancel(actual-predicted) == 0
        inverse_full += sp.prod(R(1, rates[v]) for v in selected)/predicted
        minor_records.append(dict(mask=mask, no_hit_probabilities=list(map(str, h)),
                                  green_minor=str(actual)))
    full = sp.factor(1/inverse_full)
    assert sp.cancel(a*sp.diff(full, a)-b*sp.diff(full, b)) == 0
    records = []
    for av, bv in ((R(1, 2), R(1, 3)), (R(2, 3), R(1, 4))):
        kernel = coupled.subs({a: av, b: bv})
        transition = kernel.row_join(sp.ones(3, 1)-kernel*sp.ones(3, 1))
        raw_kernel = tuple(tuple((j, transition[i, j]) for j in range(4)
                                 if transition[i, j]) for i in range(3))
        starts = [(1, 1, 1)]+[add(stable(mask, 3), v)
                               for v in range(3) for mask in range(8)]
        sleep = tuple(R(rate, rate+1) for rate in rates)
        laws, count = absorption(starts, raw_kernel, sleep)
        reverse, reverse_count = absorption(starts, raw_kernel, sleep, reverse=True)
        assert laws == reverse
        law = laws[0, :].T
        assert sum(law) == 1
        assert all(laws[1+8*v:1+8*(v+1), :].T*law == law for v in range(3))
        assert law[-1] == full.subs({a: av, b: bv})
        means = [sum(law[mask]*sum(bool(mask & (1 << v)) for v in vertices)
                     for mask in range(8)) for vertices in blocks]
        records.append(dict(routing=[str(av), str(bv)], product=str(av*bv),
                            stationary_law=list(map(str, law)),
                            full_occupation=str(law[-1]),
                            actual_block_mean_counts=list(map(str, means)),
                            raw_transient_counts=[count, reverse_count]))
    assert records[0]["full_occupation"] == records[1]["full_occupation"]
    assert records[0]["actual_block_mean_counts"] != records[1]["actual_block_mean_counts"]
    return dict(scope="exact formula and actual-marginal distinction, not an NA proof",
                base_kernel=[list(map(str, row)) for row in base.tolist()],
                rates=rates, entrances=[list(map(str, rho)) for rho in entrances],
                symbolic_green_minors=minor_records,
                symbolic_full_occupation=str(full), raw_cases=records,
                all_green_minor_identities_checked=True,
                both_raw_orders_and_all_addition_equations_checked=True)


def two_block_covariance_certificate():
    """Check finite two-block responses and distinguish two sign obstructions.

    Raw two-singleton ARW and squarefree derivatives verify the full mixed
    covariance identity and its positive-curvature witness. A separate
    abstract rank chain shows why linear rank response alone cannot imply
    the loss-PGF sign; it is not claimed to be an ARW addition kernel.
    """
    from exact import absorption, add, stable

    R = sp.Rational
    a, b, u = sp.symbols('a b u')
    t = a*b
    den = (4-3*t)*(4-t)
    law = sp.Matrix([
        (4-9*t+3*t*t+t*(a+b))/den,
        (4-3*t-t*b)/den,
        (4-3*t-t*a)/den,
        1/(4-3*t)])
    harmonic_rows = sp.Matrix([
        [1, 1, 1, 1],
        [1, -(1-a), -(1-b), (1-a)*(1-b)],
        [1, 1-a, -(1-t), -(1-a)*(1-t)],
        [1, -(1-t), 1-b, -(1-b)*(1-t)]])
    assert (harmonic_rows*law-sp.Matrix([1, 0, 0, 0])).applyfunc(sp.cancel) == sp.zeros(4, 1)
    f, g = sp.Matrix([[0, 1, 0, 1]]), sp.Matrix([[0, 0, 1, 1]])
    fg = sp.Matrix([[0, 0, 0, 1]])
    covariance = sp.factor((fg*law)[0]-(f*law)[0]*(g*law)[0])
    expected_covariance = -4*t*(4-2*(a+b)-3*t+t*(a+b)+t*t)/den**2
    assert sp.cancel(covariance-expected_covariance) == 0
    level = R(3, 4)
    point = {a: level, b: level}
    nu = law.subs(point)
    nu_a = law.diff(a).subs(point)
    nu_b = law.diff(b).subs(point)
    nu_ab = law.diff(a, b).subs(point)
    curvature = sp.factor(sp.diff(covariance, a, b).subs(point))
    assert covariance.subs(point) == -R(576, 34225)
    assert curvature == R(7851490304, 141733425625) > 0

    # Independent raw stabilization and squarefree IDLA derivatives at this point.
    kernel = (((1, level), (2, 1-level)), ((0, level), (2, 1-level)))
    starts = [(1, 1)]+[add(stable(mask, 2), v) for v in range(2) for mask in range(4)]
    raw, raw_count = absorption(starts, kernel, (R(1, 2), R(1, 2)))
    reverse, reverse_count = absorption([(1, 1)], kernel, (R(1, 2), R(1, 2)), reverse=True)
    assert raw[0, :].T == reverse.T == nu
    T1, T2 = [raw[1+4*v:1+4*(v+1), :].T for v in range(2)]
    assert T1*nu == T2*nu == nu and T1*T2 == T2*T1
    eps, eta = sp.symbols('eps eta')
    transition = sp.Matrix([[0, level+eps, 1-level-eps],
                            [level+eta, 0, 1-level-eta]])
    coefficients = nilpotent_stationary(transition, (1, 1), (eps, eta))
    assert coefficients == [nu, nu_a, nu_b, nu_ab]

    I = sp.eye(4)
    q1 = q2 = 1-level
    Q1, Q2 = q1*I+level*T2, q2*I+level*T1
    Delta = Q1*Q2-level**2*T1*T2
    assert Delta == q1*q2*I+q1*level*T1+q2*level*T2
    W = Delta.inv()
    V = [sp.diag(1, 0, 1, 0), sp.diag(1, 1, 0, 0)]
    F1 = W*(Q2*V[0]+level*T1*V[1])*nu
    F2 = W*(level*T2*V[0]+Q1*V[1])*nu
    for v, measure in enumerate((F1, F2)):
        marked, _ = marked_exit_operator(kernel, (1, 1), v, initial_states=[(1, 1)])
        assert marked == measure
    assert Q1*F1-level*T1*F2 == V[0]*nu
    assert Q2*F2-level*T2*F1 == V[1]*nu
    assert (T2-I)*F1 == nu_a and (T1-I)*F2 == nu_b

    B_a = [(T2-I)*W*Q2, (T2-I)*W*(level*T1)]
    B_b = [(T1-I)*W*(level*T2), (T1-I)*W*Q1]
    E = [T2*Q2*(T1*Q1-I)*W**2,
         T1*Q1*(T2*Q2-I)*W**2]
    G_a = sum((B_a[i]*V[i] for i in range(2)), sp.zeros(4))
    G_b = sum((B_b[i]*V[i] for i in range(2)), sp.zeros(4))
    H_ab = G_a*G_b+G_b*G_a
    H_ab -= sum((B_a[i]*B_b[j]*V[i]*V[j]
                 for i in range(2) for j in range(2)), sp.zeros(4))
    H_ab += sum((E[i]*V[i] for i in range(2)), sp.zeros(4))
    assert G_a*nu == nu_a and G_b*nu == nu_b
    assert H_ab*nu == nu_ab
    mf, mg = (f*nu)[0], (g*nu)[0]
    centered_product = sp.Matrix([[(f[j]-mf)*(g[j]-mg) for j in range(4)]])
    operator_curvature = ((centered_product*H_ab*nu)[0]
                          -(f*nu_a)[0]*(g*nu_b)[0]
                          -(f*nu_b)[0]*(g*nu_a)[0])
    assert operator_curvature == curvature

    # The one-way boundary identity is symbolic in the forward routing strength.
    refresh_law = sp.Matrix([R(1, 2), R(1, 2)])
    refresh = sp.Matrix.hstack(refresh_law, refresh_law)
    psi = sp.expand(sum(refresh_law[x]*refresh[y, x]*(y-R(1, 2))*u**(1+x-y)
                        for x in range(2) for y in range(2)))
    assert psi == (1-u**2)/8
    boundary_derivative = sp.factor(sp.diff(covariance, b).subs(b, 0))
    assert sp.cancel(boundary_derivative+a*(2-a)/32) == 0
    assert sp.cancel(boundary_derivative+R(1, 4)*psi.subs(u, 1-a)) == 0

    # Stationarity, reversible lazy dynamics, skip-free upward counts, and strict
    # linear rank response do not by themselves imply the required loss-PGF order.
    abstract_law = sp.Matrix([R(1, 10), R(1, 10), R(4, 5)])
    abstract_T = sp.Matrix([[R(3, 4), R(1, 4), 0],
                           [R(1, 4), R(37, 50), R(1, 800)],
                           [0, R(1, 100), R(799, 800)]])
    assert abstract_T*abstract_law == abstract_law
    assert sp.ones(1, 3)*abstract_T == sp.ones(1, 3)
    assert abstract_T*sp.diag(*abstract_law) == sp.diag(*abstract_law)*abstract_T.T
    assert all(abstract_T[i, i] >= R(1, 2) for i in range(3))
    assert all(abstract_T[y, x] == 0 for x in range(3) for y in range(3) if y > x+1)
    N = sp.diag(0, 1, 2)
    rank_response = (sp.eye(3)-abstract_T)*N*abstract_law
    assert [sum(rank_response[r:]) for r in (1, 2)] == [R(1, 40), R(1, 1000)]
    mean = (sp.Matrix([[0, 1, 2]])*abstract_law)[0]
    abstract_psi = sp.factor(sum(
        abstract_law[x]*abstract_T[y, x]*(y-mean)*u**(1+x-y)
        for x in range(3) for y in range(3) if abstract_T[y, x]))
    assert sp.expand(abstract_psi-(1-u)*(108*u-43)/2500) == 0
    assert abstract_psi.subs(u, R(1, 4)) == -R(3, 625)

    return {
        'actual_arw_case': {
            'base_kernels': [[[0]], [[0]]], 'sleep_rates': [1, 1],
            'a': str(level), 'b': str(level),
            'stationary_law': list(map(str, nu)),
            'covariance': str(covariance.subs(point)),
            'mixed_covariance_derivative': str(curvature),
            'raw_transient_counts': [raw_count, reverse_count],
            'raw_stationarity_harmonic_and_squarefree_derivatives_checked': True,
            'raw_marked_loss_balance_and_first_flows_checked': True,
            'full_mixed_operator_identity_and_centered_covariance_checked': True,
            'one_way_covariance_derivative': str(boundary_derivative)},
        'abstract_rank_chain_only': {
            'stationary_law': list(map(str, abstract_law)),
            'column_transition_matrix': [list(map(str, row)) for row in abstract_T.tolist()],
            'strict_rank_tail_responses': ['1/40', '1/1000'],
            'terminal_count_loss_pgf_covariance': str(abstract_psi),
            'value_at_u_one_quarter': str(abstract_psi.subs(u, R(1, 4))),
            'not_claimed_to_be_an_ARW_addition_kernel': True}}


def loss_universal_bound_certificate():
    """Check sharpness and the restricted-entrance endpoint repair.

    The general bound is proved by stationarity and the exact nonnegative
    slack in report (67u). Reuse the saved abstract obstruction; the one
    product ARW case discriminates the invalid L>=3 strictness witness
    from the valid A=2,B=1,L=2 witness. No parameter sweep is performed.
    """
    u, epsilon = sp.symbols("u epsilon")
    R = sp.Rational
    T = sp.Matrix([[R(3, 4), R(1, 4), 0],
                   [R(1, 4), R(37, 50), R(1, 800)],
                   [0, R(1, 100), R(799, 800)]])
    law = sp.Matrix([R(1, 10), R(1, 10), R(4, 5)])
    assert T*law == law
    mean = sum(j*law[j] for j in range(3))
    variance = sum(law[a]*T[b, a]*(a-b)**2
                   for a in range(3) for b in range(3))
    psi = sp.factor(sum(law[a]*T[b, a]*(b-mean)*u**(1+a-b)
                        for a in range(3) for b in range(3) if T[b, a]))
    bound = variance*(1-u)*(1-mean*(1-u))/2
    assert sp.expand(psi-bound-(1-u)**2/1000) == 0
    sharp_mean = 2-3*epsilon
    sharp_flows = {(0, 1): epsilon, (1, 0): epsilon, (2, 2): 1-2*epsilon}
    sharp_psi = sum(weight*(b-sharp_mean)*u**(1+a-b)
                    for (a, b), weight in sharp_flows.items())
    assert sp.expand(sharp_psi-epsilon*(1-u)*(1-sharp_mean*(1-u))) == 0

    # K=0 on two sites, both rates 3, entrance only at site 0.
    p = R(3, 4)
    joint = []
    for old0 in (0, 1):
        for old1 in (0, 1):
            initial = (p if old0 else 1-p)*(p if old1 else 1-p)
            for new0 in (0, 1):
                weight = initial*(p if new0 else 1-p)
                joint.append((old0+old1, new0+old1, 1+old0-new0, weight))
    assert sum(weight for a, b, loss, weight in joint) == 1
    arw_mean = sum(b*weight for a, b, loss, weight in joint)
    arw_variance = sum((loss-1)**2*weight for a, b, loss, weight in joint)
    arw_psi = sp.factor(sum((b-arw_mean)*u**loss*weight
                            for a, b, loss, weight in joint))
    eta = sum(weight for a, b, loss, weight in joint if (a, b, loss) == (2, 1, 2))
    arw_bound = arw_variance*(1-u)*(1-arw_mean*(1-u))/2
    assert arw_mean == R(3, 2) and eta == R(9, 64)
    assert max(loss for a, b, loss, weight in joint) == 2
    assert sp.expand(arw_psi-arw_bound-eta*(1-u)**2) == 0
    assert arw_psi.subs(u, R(1, 3)) == R(1, 16)
    return {
        "scope": "universal analytic bound; existing abstract obstruction and one restricted-entrance strictness discriminator",
        "bound": "Psi(u)>=Var(L)*(1-u)*(1-m*(1-u))/2; m=E terminal count",
        "exact_slack": "(1-u)*(E[B*J]+m*((1-u)*Var(L)/2-E[J])); J=d-sum_{r=1}^d u^r on Delta=-d<0, zero otherwise",
        "canonical_ARW_strict_interval": "max(0,1-1/m)<=u<1; no conditional loss order inferred",
        "abstract_obstruction": {"mean": str(mean), "variance_loss": str(variance),
                                 "Psi": str(psi), "bound_residual": "(1-u)^2/1000"},
        "sharp_abstract_flow": {"flow": "epsilon on 0->1 and1->0; self-loop mass 1-2epsilon at2",
                                "mean": str(sharp_mean), "Psi": str(sp.factor(sharp_psi)),
                                "bound_is_equality": True, "not_an_ARW_counterexample": True},
        "restricted_entrance_ARW": {"kernel": [[0, 0], [0, 0]], "rates": [3, 3],
                                   "entrance": [1, 0], "mean": str(arw_mean),
                                   "variance_loss": str(arw_variance), "maximum_loss": 2,
                                   "A2_B1_L2_probability": str(eta), "Psi": str(arw_psi),
                                   "exact_slack": "9*(1-u)^2/64", "Psi_at_u_one_third": "1/16"},
        "one_way_singleton_consequence": "strict count NC for fixed 0<a<=min(1,1/m) and sufficiently small positive reverse routing; threshold depends on fixed data"}


def loss_pgf_certificate():
    """Two fixed all-u loss-PGF checks, conditional orders, and source obstructions."""
    from exact import absorption, add, graph_kernel, stable

    R = sp.Rational
    u, s = sp.symbols("u s")
    cases = (
        ("rayleigh_obstruction3_lambda64",
         ((1, 3), (0, 2, 3), (1,)), (64, 64, 64)),
        ("endpoint3_lambda1",
         ((1, 3), (0, 2), (1,)), (1, 1, 1)))
    records = []
    for name, neighbors, rates in cases:
        n, size = len(neighbors), 1 << len(neighbors)
        kernel = graph_kernel(neighbors)
        weights = sp.Matrix([len(row) for row in neighbors])
        killing = sp.Matrix([dict(row).get(n, 0) for row in kernel])
        sink_flux = sp.diag(*weights)*killing
        drive = sink_flux/sum(sink_flux)
        internal = sp.zeros(n)
        for i, row in enumerate(kernel):
            for j, probability in row:
                if j < n:
                    internal[i, j] = probability
        assert sp.diag(*weights)*internal == internal.T*sp.diag(*weights)
        starts = [(1,)*n]+[add(stable(mask, n), v)
                           for v in range(n) for mask in range(size)]
        sleep = tuple(R(rate, rate+1) for rate in rates)
        raw, raw_count = absorption(starts, kernel, sleep)
        reverse, reverse_count = absorption(starts, kernel, sleep, reverse=True)
        assert raw == reverse
        nu = raw[0, :].T
        additions = [raw[1+v*size:1+(v+1)*size, :].T for v in range(n)]
        assert all(T*nu == nu for T in additions)
        T = sum((drive[v]*additions[v] for v in range(n)), sp.zeros(size))
        ranks = [mask.bit_count() for mask in range(size)]
        mean = sum(nu[mask]*ranks[mask] for mask in range(size))
        joint = sp.zeros(n+2, n+1)
        loss_operator = sp.zeros(size)
        for before in range(size):
            for after in range(size):
                loss = 1+ranks[before]-ranks[after]
                if T[after, before]:
                    assert loss >= 0
                    joint[loss, ranks[after]] += nu[before]*T[after, before]
                    loss_operator[after, before] = T[after, before]*u**loss
        assert sum(joint) == 1
        assert sum(ell*sum(joint.row(ell)) for ell in range(n+2)) == 1
        psi = sp.Poly(sum((j-mean)*joint[ell, j]*u**ell
                         for ell in range(n+2) for j in range(n+1)), u)
        assert psi.eval(1) == 0
        loss_variance = sum((ell-1)**2*sum(joint.row(ell)) for ell in range(n+2))
        assert psi.diff().eval(1) == -loss_variance/2
        quotient = sp.Poly(sp.cancel(psi.as_expr()/(1-u)), u)
        assert quotient.eval(1) == loss_variance/2 > 0
        degree = quotient.degree()
        odds = sp.Poly(sp.cancel((1+s)**degree
                                *quotient.as_expr().subs(u, s/(1+s))), s)
        ordinary_coefficients = [quotient.nth(i) for i in range(degree+1)]
        odds_coefficients = [odds.nth(i) for i in range(degree+1)]
        assert all(x >= 0 for x in ordinary_coefficients)
        assert all(x > 0 for x in odds_coefficients)
        content, primitive = quotient.primitive()
        odds_content, odds_primitive = odds.primitive()
        assert quotient.eval(0) == psi.eval(0)
        # The power coefficients are exactly Cov(N(Y),1_{L<=k}).
        for threshold, coefficient in enumerate(ordinary_coefficients):
            direct = sum((j-mean)*joint[ell, j]
                         for ell in range(threshold+1) for j in range(n+1))
            assert coefficient == direct
        rank_law = [sum(joint[:, j]) for j in range(n+1)]
        assert rank_law == [sum(nu[mask] for mask in range(size) if ranks[mask] == j)
                            for j in range(n+1)]
        conditional_pgfs = [
            sp.Poly(sum(joint[ell, j]*u**ell for ell in range(n+2))/rank_law[j], u)
            for j in range(n+1)]
        regression = []
        for j in range(1, n+1):
            gap = sp.Poly(sp.cancel((conditional_pgfs[j].as_expr()
                                    -conditional_pgfs[j-1].as_expr())/(1-u)), u)
            gap_content, gap_primitive = gap.primitive()
            coefficients = [gap.nth(k) for k in range(gap.degree()+1)]
            assert all(x > 0 for x in coefficients)
            for threshold, coefficient in enumerate(coefficients):
                cdf_difference = (sum(joint[ell, j] for ell in range(threshold+1))/rank_law[j]
                                  -sum(joint[ell, j-1] for ell in range(threshold+1))/rank_law[j-1])
                assert coefficient == cdf_difference
            regression.append({
                "upper_rank": j,
                "content": str(gap_content),
                "primitive_coefficients_ascending": [
                    str(gap_primitive.nth(k)) for k in range(gap.degree()+1)],
                "all_nontrivial_loss_cdf_differences_positive": all(x > 0 for x in coefficients)})
        # A positive-weight first-loss source: start full, add at killed vertex 0,
        # kill the surplus immediately, then stabilize one active at 0 and sleepers
        # at every other site. Its mean exceeds the stationary mean, so a pointwise
        # positive-source argument cannot prove the covariance sign.
        assert drive[0] > 0 and killing[0] > 0
        post_loss = (1,)+(-1,)*(n-1)
        post_law, post_count = absorption([post_loss], kernel, sleep)
        post_reverse, post_reverse_count = absorption([post_loss], kernel, sleep, reverse=True)
        assert post_law == post_reverse
        post_mean = sum(post_law[0, mask]*ranks[mask] for mask in range(size))
        assert post_mean > mean
        # Independently check the loss-conjugation/Perron identity.
        tilted_nu = sp.Matrix([u**(-ranks[mask])*nu[mask] for mask in range(size)])
        assert (loss_operator*tilted_nu-u*tilted_nu).applyfunc(sp.cancel) == sp.zeros(size, 1)
        record = {
            "case": name, "neighbors_including_sink": neighbors, "sleep_rates": rates,
            "reversible_vertex_weights": list(map(str, weights)),
            "original_killing": list(map(str, killing)),
            "entrance_proportional_sink_conductance": list(map(str, drive)),
            "raw_transient_counts": [raw_count, reverse_count],
            "stationary_law": list(map(str, nu)), "mean_count": str(mean),
            "stationary_loss_variance": str(loss_variance),
            "loss_terminal_rank_joint_rows": [
                list(map(str, row)) for row in joint.tolist()],
            "psi_polynomial": str(sp.factor(psi.as_expr())),
            "quotient_content": str(content),
            "quotient_primitive_coefficients_ascending": [
                str(primitive.nth(i)) for i in range(degree+1)],
            "odds_content": str(odds_content),
            "odds_primitive_coefficients_ascending": [
                str(odds_primitive.nth(i)) for i in range(degree+1)],
            "bernstein_coefficients": [
                str(odds.nth(i)/sp.binomial(degree, i)) for i in range(degree+1)],
            "all_power_coefficients_nonnegative": all(x >= 0 for x in ordinary_coefficients),
            "all_odds_coefficients_positive": all(x > 0 for x in odds_coefficients),
            "conditional_loss_pgfs_by_terminal_rank": [
                str(sp.factor(p.as_expr())) for p in conditional_pgfs],
            "adjacent_conditional_loss_cdf_comparison": regression,
            "negative_first_loss_source_bracket": {
                "post_loss_state": post_loss,
                "positive_initial_path_weight": str(nu[-1]*drive[0]*killing[0]),
                "post_loss_expected_terminal_count": str(post_mean),
                "stationary_mean_minus_post_loss_mean": str(mean-post_mean),
                "both_raw_toppling_orders_checked": True,
                "raw_transient_counts": [post_count, post_reverse_count]},
            "both_raw_orders_all_addition_equations_and_loss_conjugation_checked": True,
            "coefficient_loss_threshold_covariance_identity_checked": True}
        records.append(record)
    return {"scope": "two prescribed exact all-u checks, not a general sign theorem",
                      "cases": records}


def loss_entrance_certificate():
    """One new entrance and its exact mixture segment on the fixed rate-64 block."""
    from exact import absorption, add, graph_kernel, stable

    R = sp.Rational
    neighbors = ((1, 3), (0, 2, 3), (1,))
    rates = (64, 64, 64)
    n, size = 3, 8
    kernel = graph_kernel(neighbors)
    sleep = tuple(R(rate, rate+1) for rate in rates)
    # Evaluate this one kernel once in each order; reuse its addition matrices
    # for the old sink-conductance entrance and the single new leaf entrance.
    starts = [(1,)*n]+[add(stable(mask, n), v)
                       for v in range(n) for mask in range(size)]
    raw, count = absorption(starts, kernel, sleep)
    reverse, reverse_count = absorption(starts, kernel, sleep, reverse=True)
    assert raw == reverse
    nu = raw[0, :].T
    additions = [raw[1+v*size:1+(v+1)*size, :].T for v in range(n)]
    assert all(T*nu == nu for T in additions)
    ranks = [mask.bit_count() for mask in range(size)]
    mean = sum(ranks[mask]*nu[mask] for mask in range(size))
    rank_law = [sum(nu[mask] for mask in range(size) if ranks[mask] == j)
                for j in range(n+1)]
    u, s, t = sp.symbols('u s t')

    def joint_law(T):
        joint = sp.zeros(n+2, n+1)
        for x in range(size):
            for y in range(size):
                if T[y, x]:
                    loss = 1+ranks[x]-ranks[y]
                    assert loss >= 0
                    joint[loss, ranks[y]] += nu[x]*T[y, x]
        assert sum(joint) == 1
        assert [sum(joint[:, j]) for j in range(n+1)] == rank_law
        assert sum(ell*sum(joint.row(ell)) for ell in range(n+2)) == 1
        return joint

    baseline = joint_law((additions[0]+additions[1])/2)
    joint = joint_law(additions[2])
    mixture = (1-t)*baseline+t*joint
    baseline_gaps, gaps = [], []
    for j in range(1, n+1):
        def gap_polynomial(table):
            return sp.Poly(sum(
                (sum(table[ell, j] for ell in range(k+1))/rank_law[j]
                 -sum(table[ell, j-1] for ell in range(k+1))/rank_law[j-1])*u**k
                for k in range(n+1)), u)

        gap0, gap1, gap_t = map(gap_polynomial, (baseline, joint, mixture))
        assert sp.expand(gap_t.as_expr()-(1-t)*gap0.as_expr()-t*gap1.as_expr()) == 0
        assert gap0.degree() == gap1.degree() == n+1-j
        assert all(gap0.nth(k) > 0 and gap1.nth(k) > 0
                   for k in range(gap1.degree()+1))
        for table, gap in ((baseline, gap0), (joint, gap1)):
            pgf_difference = sum((table[ell, j]/rank_law[j]
                                  -table[ell, j-1]/rank_law[j-1])*u**ell
                                 for ell in range(n+2))
            assert sp.expand(pgf_difference-(1-u)*gap.as_expr()) == 0
        for gap, target in ((gap0, baseline_gaps), (gap1, gaps)):
            content, primitive = gap.primitive()
            target.append({
                'upper_rank': j, 'content': str(content),
                'primitive_coefficients_ascending': [
                    str(primitive.nth(k)) for k in range(gap.degree()+1)],
                'nonzero_conditional_cdf_gaps_positive': True,
                'negative_cdf_thresholds': []})

    psi = sp.Poly(sum((j-mean)*joint[ell, j]*u**ell
                     for ell in range(n+2) for j in range(n+1)), u)
    quotient = sp.Poly(sp.cancel(psi.as_expr()/(1-u)), u)
    loss_variance = sum((ell-1)**2*sum(joint.row(ell)) for ell in range(n+2))
    assert psi.eval(1) == 0 and psi.diff().eval(1) == -loss_variance/2
    content, primitive = quotient.primitive()
    odds = sp.Poly(sp.cancel((1+s)**quotient.degree()
                            *quotient.as_expr().subs(u, s/(1+s))), s)
    odds_content, odds_primitive = odds.primitive()
    assert all(quotient.nth(k) > 0 and odds.nth(k) > 0
               for k in range(quotient.degree()+1))
    return {
        'scope': 'one new entrance and its mixture with the original entrance; not arbitrary entrances',
        'case': 'rayleigh_obstruction3_lambda64',
        'neighbors_including_sink': neighbors, 'sleep_rates': rates,
        'entrance': ['0', '0', '1'],
        'entrance_vertex_has_original_killing': False,
        'raw_transient_counts': [count, reverse_count],
        'both_raw_orders_and_site_addition_stationarity_checked': True,
        'stationary_law': list(map(str, nu)),
        'joint_loss_terminal_rank_rows': [list(map(str, row)) for row in joint.tolist()],
        'conditional_loss_cdf_comparisons': gaps,
        'quotient_content': str(content),
        'quotient_primitive_coefficients_ascending': [
            str(primitive.nth(k)) for k in range(quotient.degree()+1)],
        'odds_content': str(odds_content),
        'odds_primitive_coefficients_ascending': [
            str(odds_primitive.nth(k)) for k in range(quotient.degree()+1)],
        'all_power_coefficients_positive': True,
        'all_odds_coefficients_positive': True,
        'mixture_segment': {
            'entrance': ['(1-t)/2', '(1-t)/2', 't'], 'parameter_interval': ['0', '1'],
            'baseline_entrance': ['1/2', '1/2', '0'],
            'baseline_conditional_loss_cdf_comparisons': baseline_gaps,
            'conditional_cdf_gaps_affine_in_entrance_checked': True,
            'all_nine_nontrivial_gaps_strict_on_entire_segment': True,
            'unequal_weights_on_vertices_zero_and_one_not_claimed': True}}


def directed_cycle_loss_certificate():
    """Exact cycle feedback construction and all-rate two-site loss order.

    Serializing local inputs gives iid refresh bits eta_i(k) and cumulative
    forward maps f_i(k)=B_i(k+x_i-eta_i(k)), which are nondecreasing. Cutting
    the last edge gives F=f_(n-1) o ... o f_0. Releasing buffered returns
    iterates r -> F(1+r), hence chooses the least fixed point r_*. Since
    m-F(m) increases by at most one per step, M=1+r_* is its first passage
    to one. At that stopping index Y_i=eta_i(k_i(M)), and conservation gives
    L=1+|X|-|Y|. This identifies the joint law for every cycle length, but
    the stopping index biases the terminal bits. A single particle from
    empty sleeps at its source with probability p_0/(1-c), where
    c=product_i((1-p_i)*t_i); freezing one bit per site instead gives p_0.

    Only the two-site conditional loss order is proved here: harmonic
    interpolation gives its exact stationary addition kernel, and every
    nontrivial conditional CDF gap has positive numerator and denominator
    coefficients after a=A/(1+A), b=B/(1+B). Relabeling and linearity give
    arbitrary entrance; continuity handles one zero residual killing row.
    """
    from exact import absorption, add, stable

    x,y,a,b,A,B,u=sp.symbols('x y a b A B u')
    weights=[(-(1-a)/x,-(1-b)/y),
             (1-a,-(1-a*b)/y),
             (-(1-a*b)/x,1-b),(1,1)]
    H=sp.Matrix([[1,c,d,c*d] for c,d in weights])
    Hi=H.inv().applyfunc(sp.factor)
    nu=Hi[:,3]
    T=(Hi*sp.diag(0,1-a,0,1)*H).applyfunc(sp.factor)
    assert (T*nu-nu).applyfunc(sp.cancel)==sp.zeros(4,1)
    assert (sp.ones(1,4)*T-sp.ones(1,4)).applyfunc(sp.cancel)==sp.zeros(1,4)
    ranks=[0,1,1,2]
    joint=sp.zeros(4,3)
    for j in range(4):
     for i in range(4):
      loss=1+ranks[i]-ranks[j]
      if loss<0:
       assert T[j,i]==0
      else:
       joint[loss,ranks[j]]+=T[j,i]*nu[i]
    joint=joint.applyfunc(sp.factor)
    ranklaw=[nu[0],nu[1]+nu[2],nu[3]]
    checks=[]
    for j in (1,2):
     for k in range(3):
      gap=sp.factor(sum(joint[ell,j] for ell in range(k+1))/ranklaw[j]-sum(joint[ell,j-1] for ell in range(k+1))/ranklaw[j-1])
      if gap==0:
       continue
      # Map each forward probability into (0,1) by a=A/(1+A), b=B/(1+B).
      numerator,denominator=sp.fraction(sp.factor(gap.subs({a:A/(1+A),b:B/(1+B)})))
      np=sp.Poly(numerator,x,y,A,B)
      dp=sp.Poly(denominator,x,y,A,B)
      if all(c<0 for c in dp.coeffs()):
       np,dp=-np,-dp
      assert all(c>0 for c in np.coeffs())
      assert all(c>0 for c in dp.coeffs())
      # numerator/denominator were obtained by exact factor/substitution above;
      # Poly only changes their representation, not their rational function.
      checks.append({'upper_rank':j,'loss_threshold':k,'factored_gap':str(gap),
                     'odds_numerator':str(np.as_expr()),'odds_denominator':str(dp.as_expr()),
                     'numerator_terms':len(np.terms()),'denominator_terms':len(dp.terms()),
                     'numerator_all_positive':all(c>0 for c in np.coeffs()),
                     'denominator_all_positive':all(c>0 for c in dp.coeffs()),
                     'min_numerator_coefficient':str(min(np.coeffs()))})
    # One asymmetric raw check detects an orientation/rate error in the harmonic
    # representation. No other graphs or numerical parameters are evaluated.
    point={x:sp.Rational(2),y:sp.Rational(3),a:sp.Rational(1,3),b:sp.Rational(2,3)}
    kernel=(((1,sp.Rational(1,3)),(2,sp.Rational(2,3))),
            ((0,sp.Rational(2,3)),(2,sp.Rational(1,3))))
    sleep=(sp.Rational(2,3),sp.Rational(3,4))
    starts=[(1,1)]+[add(stable(mask,2),0) for mask in range(4)]
    raw,raw_count=absorption(starts,kernel,sleep)
    reverse,reverse_count=absorption(starts,kernel,sleep,reverse=True)
    assert raw==reverse
    assert raw[0,:].T==nu.subs(point)
    assert raw[1:,:].T==T.subs(point)
    assert H*nu==sp.Matrix([0,0,0,1]) or (H*nu-sp.Matrix([0,0,0,1])).applyfunc(sp.cancel)==sp.zeros(4,1)
    # Check the positive determinant by explicit row operations, avoiding a
    # large generic rational determinant expansion.
    d=(1+x)*(1+y)-a*b
    e=x*y+(1-a*b)*(1+x+y)
    scaled=H*sp.diag(1,x,y,x*y)
    reduced=sp.Matrix([list((scaled[i,:]-scaled[0,:])[0,1:]) for i in range(1,4)])
    reduced[0,:]=reduced[0,:]/(1-a)
    reduced[1,:]=reduced[1,:]/(1-b)
    reduced[2,:]=reduced[2,:]-reduced[0,:]-reduced[1,:]
    assert (reduced[2,:]-sp.Matrix([[0,0,e]])).applyfunc(sp.cancel)==sp.zeros(1,3)
    assert sp.cancel(reduced[:2,:2].det()-d)==0
    harmonic_determinant=sp.factor((1-a)*(1-b)*d*e/(x*x*y*y))
    # The same raw addition matrix checks the terminal-coin obstruction.
    p0=x/(1+x)
    p1=y/(1+y)
    no_sleep_round=(1-p0)*a*(1-p1)*b
    assert sp.cancel(T[1,0]-p0/(1-no_sleep_round))==0
    terminal_coin_gap=sp.factor(T[1,0]-p0)
    assert sp.cancel(terminal_coin_gap-x*a*b/((1+x)*d))==0
    result={'rates':['x','y'],'kernel':[['0','a'],['b','0']],
            'entrance':['1','0'],'stationary_law':list(map(str,nu)),
            'addition_matrix':[list(map(str,row)) for row in T.tolist()],
            'joint_loss_rank_rows':[list(map(str,row)) for row in joint.tolist()],
            'conditional_loss_gaps':checks,
            'harmonic_basis_determinant':str(harmonic_determinant),
            'raw_validation':{'a':'1/3','b':'2/3','rates':['2','3'],
                              'transient_counts':[raw_count,reverse_count],
                              'stationary_law_and_entrance_addition_both_orders_checked':True}}
    return {
        'scope': 'all-length exact cycle flow; all-parameter conditional loss order only for two vertices',
        'cycle_flow': {
            'hypotheses': 'after self-loop erasure, forward probabilities t_i in (0,1], product t_i<1, sleep probabilities p_i in (0,1)',
            'local_processor': 'eta_i(0)=x_i; eta_i(k), k>=1, iid Bernoulli(p_i); independent B_i counts iid Bernoulli(t_i) destinations',
            'forward_map': 'f_i(k)=B_i(k+x_i-eta_i(k)); F=f_(n-1) o ... o f_0',
            'least_fixed_point': 'r_*=min{r>=0:F(1+r)=r}; monotone iteration r_0=0, r_(h+1)=F(1+r_h)',
            'first_passage': 'M=1+r_*=min{m>=1:m-F(m)>=1}; M-F(M)=1',
            'terminal_configuration': 'k_0=M; k_(i+1)=f_i(k_i); Y_i=eta_i(k_i)',
            'total_loss': 'L=1+|X|-|Y|',
            'status': 'exact construction for every cycle length; not a conditional-order proof',
            'remaining_inequality': 'pi_j sum_(i<=j+k) W_(i,j+1) >= pi_(j+1) sum_(i<=j+k-1) W_(i,j), W_(i,j)=Pr(|X|=i,|Y|=j)'},
        'one_terminal_coin_obstruction': {
            'initial_state': 'empty', 'source': '0',
            'all_length_round_survival': 'c=product_i((1-p_i)*t_i)',
            'actual_probability_of_sleep_at_source': 'p_0/(1-c)',
            'frozen_one_coin_probability': 'p_0',
            'strict_difference': 'p_0*c/(1-c)>0',
            'two_site_symbolic_difference': str(terminal_coin_gap),
            'same_raw_addition_matrix_checks_identity': True,
            'does_not_rule_out_other_weighted_terminal_bit_constructions': True},
        'two_site_all_parameter': result,
        'two_site_theorem_scope': {
            'rates': 'arbitrary positive finite rates; self-loops erased rowwise',
            'routing': '0<a,b<=1, ab<1',
            'entrance': 'arbitrary; second pure entrance by relabeling, mixtures by linearity',
            'conditional_loss_order': 'L conditional on terminal rank decreases stochastically',
            'strictness': 'all nontrivial CDF gaps strict when 0<a,b<1; Cov(N(Y),u^L)>0 for every 0<=u<1 throughout the killed domain',
            'boundary_extension': 'continuity at a=1 or b=1 with ab<1',
            'longer_cycle_order_not_claimed': True}}


def directed_cycle_record_survival_certificate(two_site):
    """Check marked reversal and finite last exits at the existing raw point.

    The argument is the freshly computed two_site_all_parameter certificate;
    no saved certificate or scratch module is imported. Its raw validation
    already checks the closed-cycle stationary law and source-0 addition in
    both toppling orders. The independent open-path construction below uses
    the same rates and routing, with complementary sleep probabilities.

    For terminal y, m=n-|y| and D=c+|s|-m. Surviving first steps lie on
    D=1. Retain 0<=c<=m+1 and D>=1; absorb success at c>=m+2. Although
    c>=m+1 guarantees survival, D=1 can still recur at (empty,m+1), so
    that earlier cutoff would lose the last boundary visit. In row convention
    h=b(I-Q)^-1 g; boundary occupation times the probability of success
    without another boundary visit gives the joint survival/loss masses.

    Retaining those visits gives W(z,w) at c=1+|z|-|w|, where z=Cy.
    Each boundary escape factor equals beta(w)=h(Cw). In row convention
    pi_p(x) P_p(x,y)=pi_p(y) W(Cy,Cx); W beta=beta and pi_q W=pi_q.
    Its Doob transform W(z,w) beta(w)/beta(z) is the complemented
    stationary reversal of P_p. The closed harmonic addition matrix T
    below uses columns, so P_p(x,y)=T[y,x].

    This is one exact consistency check of the proved representation, not
    a proof of the open all-strength normalized loss-order inequality.
    """
    from collections import defaultdict

    def product_law(probabilities):
        return [sp.prod(probabilities[i] if mask >> i & 1 else 1-probabilities[i]
                        for i in range(len(probabilities)))
                for mask in range(1 << len(probabilities))]

    def marked_path(probabilities, forward_probabilities):
        """Column input; outputs (mask, physical losses, buffered returns)."""
        n = len(probabilities)
        columns = []
        for old in range(1 << n):
            law = defaultdict(lambda: sp.S.Zero)

            def visit(i, incoming, new, loss, weight):
                if i == n:
                    law[new, loss, incoming] += weight
                    return
                oldbit = old >> i & 1
                if incoming == 0:
                    visit(i+1, 0, new | oldbit << i, loss, weight)
                    return
                for bit in (0, 1):
                    emit = incoming + oldbit - bit
                    bitweight = probabilities[i] if bit else 1-probabilities[i]
                    for forward in range(emit+1):
                        route = (sp.binomial(emit, forward)
                                 *forward_probabilities[i]**forward
                                 *(1-forward_probabilities[i])**(emit-forward))
                        if route:
                            visit(i+1, forward, new | bit << i, loss+emit-forward,
                                  weight*bitweight*route)

            visit(0, 1, 0, 0, sp.S.One)
            assert sum(law.values()) == 1
            for new, loss, returns in law:
                assert 1+old.bit_count()-new.bit_count() == loss+returns
            columns.append(dict(law))
        return columns

    def finite_survival(terminal, n, complementary_path):
        m = n-terminal.bit_count()
        states = [(s, c) for c in range(m+2) for s in range(1 << n)
                  if c+s.bit_count() >= m+1]
        index = {state: i for i, state in enumerate(states)}
        size = len(states)
        Q, g, initial = sp.zeros(size), sp.zeros(size, 1), sp.zeros(1, size)
        for i, (old, c) in enumerate(states):
            failure = sp.S.Zero
            for (new, loss, _returns), weight in complementary_path[old].items():
                if c+loss >= m+2:
                    g[i] += weight
                elif c+loss+new.bit_count() >= m+1:
                    Q[i, index[new, c+loss]] += weight
                else:
                    failure += weight
            assert sum(Q[i, :])+g[i]+failure == 1
        # The unadvanced initial phase has D=0 and is outside E_m.
        for (new, loss, _returns), weight in complementary_path[((1 << n)-1)^terminal].items():
            assert loss < m+2
            if loss+new.bit_count() >= m+1:
                assert loss+new.bit_count() == m+1
                initial[index[new, loss]] += weight
        boundary = [i for i, (s, c) in enumerate(states) if c+s.bit_count() == m+1]
        interior = [i for i in range(size) if i not in boundary]
        G = (sp.eye(size)-Q).inv()
        assert all(entry >= 0 for entry in G)
        occupation = initial*G
        if interior:
            continuation = ((sp.eye(len(interior))-Q.extract(interior, interior)).inv()
                            *g.extract(interior, [0]))
        last = [sp.S.Zero]*(n+2)
        boundary_records = []
        for i in boundary:
            escape = g[i]
            if interior:
                escape += (Q.extract([i], interior)*continuation)[0]
            assert 0 <= escape <= 1
            last[states[i][1]] += occupation[i]*escape
            state, loss = states[i]
            boundary_records.append((state, loss, occupation[i], escape))
        h = (occupation*g)[0]
        assert sum(last) == h
        assert 0 < h <= 1
        return h, last, len(states), len(boundary), boundary_records

    raw = two_site['raw_validation']
    assert raw['stationary_law_and_entrance_addition_both_orders_checked']
    rates = tuple(sp.Rational(value) for value in raw['rates'])
    t = (sp.Rational(raw['a']), sp.Rational(raw['b']))
    assert rates == (2, 3) and t == (sp.Rational(1, 3), sp.Rational(2, 3))
    n, full = 2, 3
    p = tuple(rate/(1+rate) for rate in rates)
    q = tuple(1-value for value in p)
    R = sp.prod(t)
    pi = product_law(p)
    forward, complementary = marked_path(p, t), marked_path(q, t)
    mark_checks = 0
    for old in range(1 << n):
        for new in range(1 << n):
            marks = {(loss, returns) for (state, loss, returns) in forward[old]
                     if state == new}
            marks |= {(loss, returns) for (state, loss, returns) in complementary[full^new]
                      if state == full^old}
            for loss, returns in sorted(marks):
                assert (pi[old]*forward[old].get((new, loss, returns), 0)
                        == pi[new]*complementary[full^new].get((full^old, loss, returns), 0))
                mark_checks += 1
    x, y, a, b = sp.symbols('x y a b')
    point = {x: rates[0], y: rates[1], a: t[0], b: t[1]}
    nu = sp.Matrix([sp.sympify(value).subs(point) for value in two_site['stationary_law']])
    T = sp.Matrix([[sp.sympify(value).subs(point) for value in row]
                   for row in two_site['addition_matrix']])
    assert T*nu == nu
    weights = [(-(1-t[0])/rates[0], -(1-t[1])/rates[1]),
               (1-t[0], -(1-R)/rates[1]), (-(1-R)/rates[0], 1-t[1]), (1, 1)]
    H = sp.Matrix([[1, c, d, c*d] for c, d in weights])
    assert H*nu == sp.Matrix([0, 0, 0, 1])
    records = []
    weighted_h = sp.S.Zero
    W, beta = sp.zeros(1 << n), sp.zeros(1 << n, 1)
    boundary_escapes = []
    for terminal in range(1 << n):
        h, joint, dimension, boundary_dimension, visits = finite_survival(terminal, n, complementary)
        start = full^terminal
        beta[start] = h
        for state, loss, occupation, escape in visits:
            assert loss == 1+start.bit_count()-state.bit_count()
            W[start, state] = occupation
            boundary_escapes.append((state, escape))
        assert h == (1-R)*nu[terminal]/pi[terminal]
        expected_joint = [sp.S.Zero]*(n+2)
        for old in range(1 << n):
            loss = 1+old.bit_count()-terminal.bit_count()
            if loss < 0:
                assert T[terminal, old] == 0
            else:
                expected_joint[loss] += (1-R)*nu[old]*T[terminal, old]/pi[terminal]
        assert joint == expected_joint
        weighted_h += pi[terminal]*h
        records.append({
            'terminal_mask': terminal, 'initial_holes': n-terminal.bit_count(),
            'phase_dimension': dimension, 'boundary_dimension': boundary_dimension,
            'h': str(h), 'survival_and_loss': list(map(str, joint)),
            'conditional_loss': [str(value/h) for value in joint],
            'all_harmonic_and_stationary_joint_entries_match': True})
    assert weighted_h == 1-R
    # W is constructed from the same finite visits, independently of T.
    for state, escape in boundary_escapes:
        assert escape == beta[state]
    closed_reversal_checks = 0
    for old in range(1 << n):
        for terminal in range(1 << n):
            assert pi[old]*T[terminal, old] == pi[terminal]*W[full^terminal, full^old]
            if 1+old.bit_count()-terminal.bit_count() < 0:
                assert W[full^terminal, full^old] == 0
            closed_reversal_checks += 1
    pi_q = sp.Matrix([product_law(q)])
    assert pi_q == sp.Matrix([[pi[full^state] for state in range(1 << n)]])
    assert W*beta == beta
    assert pi_q*W == pi_q
    assert (pi_q*beta)[0] == 1-R
    doob = sp.Matrix(1 << n, 1 << n, lambda z, w: W[z, w]*beta[w]/beta[z])
    complemented_reversal = sp.Matrix(
        1 << n, 1 << n, lambda z, w: nu[full^w]*T[full^z, full^w]/nu[full^z])
    assert doob == complemented_reversal
    assert all(entry >= 0 for entry in doob)
    assert doob*sp.ones(1 << n, 1) == sp.ones(1 << n, 1)
    rho = sp.Matrix([[nu[full^state] for state in range(1 << n)]])
    assert rho == pi_q*sp.diag(*beta)/(1-R)
    assert rho*doob == rho
    return {
        'scope': 'one previously raw-certified asymmetric two-site point',
        'forward_probabilities': list(map(str, t)), 'rates': list(map(str, rates)),
        'entrance': '0', 'marked_DAG_reversal_entries_checked': mark_checks,
        'record_density': str(weighted_h), 'records': records,
        'survival_formula': 'h(y)=b_y (I-Q_m)^-1 g_m, m=n-|y|',
        'last_exit_formula': 'J_l(y)=sum_(s:(s,l) in B_m) (b_y (I-Q_m)^-1)(s,l) e_m(s,l)',
        'success_cutoff': 'physical loss >= m+2; retain the possible last visit at (empty,m+1)',
        'same_raw_stationary_and_addition_certificate_reused': True,
        'sign_inequality_tested': False,
        'normalized_order_status': 'OPEN for general cycles; these positive formulas do not sign normalized comparisons',
        'closed_reversal': {
            'matrix_convention': 'row initial, column final; binary-mask order',
            'W_construction': 'W(z,w)=(b_z G_|z|)(w,1+|z|-|w|), or zero if the loss is negative',
            'W': [list(map(str, row)) for row in W.tolist()],
            'beta': list(map(str, beta)), 'pi_q': list(map(str, pi_q)),
            'boundary_escape_entries_checked': len(boundary_escapes),
            'boundary_escape_equals_beta': True,
            'marked_closed_reversal_entries_checked': closed_reversal_checks,
            'W_beta_equals_beta': True, 'pi_q_W_equals_pi_q': True,
            'pi_q_beta': str((pi_q*beta)[0]),
            'doob_kernel': [list(map(str, row)) for row in doob.tolist()],
            'doob_equals_complemented_stationary_reversal': True,
            'doob_row_stochastic': True,
            'invariant_complemented_stationary_law': list(map(str, rho))}}


def three_cycle_loss_certificate():
    """All common rates and strengths on the homogeneous directed three-cycle.

    Rotation is transitive on every stable rank layer, so uniform entrance
    has an exact four-state rank quotient. Averaged harmonic functions give
    its four distinct eigenvalues and determine its stationary law/addition
    matrix. The stationary count-pair law for any pure entrance is the same
    by rotation. All nine nontrivial conditional loss CDF gaps have positive
    polynomial numerators/denominators after t=z/(1+z), for x,z>0. This does
    not assert a rank quotient at larger cycle lengths or unequal rates.
    """
    from itertools import combinations
    from exact import absorption, add, stable

    x,t,z=sp.symbols('x t z')
    weights=[[-(1-t)/x]*3,
             [1-t,-(1-t)/x,-(1-t*t)/x],
             [1-t*t,1-t,-(1-t**3)/x],[1,1,1]]
    H=sp.Matrix([[sum(sp.prod(b[i] for i in S) for S in combinations(range(3),j))/sp.binomial(3,j)
                  for j in range(4)] for b in weights])
    Hi=H.inv().applyfunc(sp.factor)
    nu=Hi[:,3]
    eigen=sp.diag(0,(1-t)/3,(2-t-t*t)/3,1)
    T=(Hi*eigen*H).applyfunc(sp.factor)
    assert (T*nu-nu).applyfunc(sp.cancel)==sp.zeros(4,1)
    assert (sp.ones(1,4)*T-sp.ones(1,4)).applyfunc(sp.cancel)==sp.zeros(1,4)
    joint=sp.zeros(5,4)
    for before in range(4):
     for after in range(4):
      loss=1+before-after
      if loss<0:
       assert T[after,before]==0
      else:
       joint[loss,after]=sp.factor(nu[before]*T[after,before])
    checks=[]
    for j in range(1,4):
     for k in range(4):
      gap=sp.factor(sum(joint[ell,j] for ell in range(k+1))/nu[j]
                    -sum(joint[ell,j-1] for ell in range(k+1))/nu[j-1])
      if gap==0:
       continue
      numerator,denominator=sp.fraction(sp.factor(gap.subs(t,z/(1+z))))
      np,dp=sp.Poly(numerator,x,z),sp.Poly(denominator,x,z)
      if all(c<0 for c in dp.coeffs()):
       np,dp=-np,-dp
      checks.append({'upper_rank':j,'loss_threshold':k,'gap':str(gap),
                     'positive_odds_numerator':all(c>0 for c in np.coeffs()),
                     'positive_odds_denominator':all(c>0 for c in dp.coeffs()),
                     'numerator_terms':len(np.terms()),'denominator_terms':len(dp.terms()),
                     'minimum_numerator_coefficient':str(min(np.coeffs())),
                     'odds_numerator':str(np.as_expr()),'odds_denominator':str(dp.as_expr())})
    assert len(checks)==9 and all(c['positive_odds_numerator'] and c['positive_odds_denominator'] for c in checks)
    # One independent raw quotient check, not a parameter search.
    level=sp.Rational(2,3)
    kernel=tuple((((i+1)%3,level),(3,1-level)) for i in range(3))
    sleep=(sp.Rational(2,3),)*3
    starts=[(1,1,1)]+[add(stable(mask,3),v) for v in range(3) for mask in range(8)]
    raw,raw_count=absorption(starts,kernel,sleep)
    reverse,reverse_count=absorption(starts,kernel,sleep,reverse=True)
    assert raw==reverse
    raw_nu=raw[0,:].T
    additions=[raw[1+v*8:1+(v+1)*8,:].T for v in range(3)]
    assert all(A*raw_nu==raw_nu for A in additions)
    assert sp.Matrix([sum(raw_nu[m] for m in range(8) if m.bit_count()==j) for j in range(4)])==nu.subs({x:2,t:level})
    uniform=sum(additions,sp.zeros(8))/3
    rank_T=sp.zeros(4)
    for before in range(4):
     columns=[]
     for mask in range(8):
      if mask.bit_count()==before:
       columns.append(sp.Matrix([sum(uniform[out,mask] for out in range(8) if out.bit_count()==after)
                                 for after in range(4)]))
     assert all(c==columns[0] for c in columns)
     rank_T[:,before]=columns[0]
    assert rank_T==T.subs({x:2,t:level})
    raw_joint=[]
    for A in additions:
     J=sp.zeros(5,4)
     for before in range(8):
      for after in range(8):
       if A[after,before]:
        loss=1+before.bit_count()-after.bit_count()
        assert loss>=0
        J[loss,after.bit_count()]+=raw_nu[before]*A[after,before]
     raw_joint.append(J)
    assert all(J==joint.subs({x:2,t:level}) for J in raw_joint)
    result={'scope':'homogeneous directed 3-cycle; common rate x>0, forward probability 0<t<1; any entrance by rotation',
            'harmonic_rank_basis':[list(map(str,row)) for row in H.tolist()],
            'rank_stationary_law':list(map(str,nu)),
            'rank_addition_matrix':[list(map(str,row)) for row in T.tolist()],
            'conditional_loss_gaps':checks,
            'all_nine_gaps_positive_odds_polynomial':len(checks)==9 and all(c['positive_odds_numerator'] and c['positive_odds_denominator'] for c in checks),
            'raw_check':{'common_rate':'2','forward_probability':'2/3',
                         'transient_counts':[raw_count,reverse_count],
                         'quotient_lumpability_stationarity_all_pure_sources_both_orders_checked':True}}
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", choices=("path3", "complete3", "path4", "diamond4", "star4", "bowtie5"))
    parser.add_argument("--validate", action="store_true")
    parser.add_argument("--port-deletion", action="store_true")
    parser.add_argument("--bridge-response", action="store_true")
    parser.add_argument("--marked-source", action="store_true")
    parser.add_argument("--component-cycles", action="store_true")
    parser.add_argument("--response-order-symbolic", action="store_true")
    parser.add_argument("--reinjection-flow", action="store_true")
    parser.add_argument("--distributed-bridge", action="store_true")
    parser.add_argument("--feedback-audit", action="store_true")
    parser.add_argument("--multiblock-feedback", action="store_true")
    parser.add_argument("--two-block-full-occupation", action="store_true")
    parser.add_argument("--two-block-covariance", action="store_true")
    parser.add_argument("--loss-pgf", action="store_true")
    parser.add_argument("--loss-universal-bound", action="store_true")
    parser.add_argument("--loss-entrance", action="store_true")
    parser.add_argument("--directed-cycle-loss", action="store_true")
    parser.add_argument("--three-cycle-loss", action="store_true")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    examples = {
        "path3": (3, ((0, 1), (1, 0), (1, 2), (2, 1)), (1, 2, 3)),
        "complete3": (3, tuple((i, j) for i in range(3) for j in range(3) if i != j), (2, 3, 5)),
        "path4": (4, ((0, 1), (1, 0), (1, 2), (2, 1), (2, 3), (3, 2)), (1, 2, 3, 5)),
        "diamond4": (4, ((0, 1), (1, 0), (1, 3), (3, 1), (0, 2), (2, 0), (2, 3), (3, 2)), (1, 2, 3, 5)),
        "star4": (4, ((0, 1), (1, 0), (0, 2), (2, 0), (0, 3), (3, 0)), (1, 2, 3, 5)),
        "bowtie5": (5, ((0, 1), (1, 2), (2, 0), (2, 3), (3, 4), (4, 2)), (1, 2, 3, 5, 7))}
    result = {}
    if args.case or not (args.validate or args.port_deletion or args.bridge_response or args.marked_source or args.component_cycles or args.response_order_symbolic or args.reinjection_flow or args.distributed_bridge or args.feedback_audit or args.multiblock_feedback or args.two_block_full_occupation or args.two_block_covariance or args.loss_pgf or args.loss_universal_bound or args.loss_entrance or args.directed_cycle_loss or args.three_cycle_loss):
        name = args.case or "path3"
        result[name] = squarefree_moments(*examples[name])
    if args.validate:
        result["validation"] = symbolic_kernel_certificate()
        result["nonlinear_events"] = nonlinear_event_certificate()
    if args.port_deletion:
        result["port_deletion"] = port_deletion_certificate()
    if args.bridge_response:
        result["bridge_response"] = bridge_response_certificate()
        result["self_loop_response_order"] = self_loop_response_order_certificate()
    if args.marked_source:
        result["marked_source_obstruction"] = marked_source_certificate()
    if args.component_cycles:
        result["component_cycles"] = component_cycle_certificate()
    if args.response_order_symbolic:
        result["symbolic_response_order"] = symbolic_response_order_certificate()
        result["two_site_response_order"] = two_site_response_certificate()
    if args.reinjection_flow:
        result["reinjection_flow"] = reinjection_flow_certificate()
    if args.distributed_bridge:
        result["distributed_bridge"] = distributed_bridge_certificate()
    if args.feedback_audit:
        result["feedback_audit"] = feedback_audit_certificate()
    if args.multiblock_feedback:
        result["multiblock_feedback"] = multiblock_feedback_certificate()
    if args.two_block_full_occupation:
        result["two_block_full_occupation"] = two_block_full_occupation_certificate()
    if args.two_block_covariance:
        result["two_block_covariance"] = two_block_covariance_certificate()
    if args.loss_universal_bound:
        result["loss_universal_bound"] = loss_universal_bound_certificate()
    if args.loss_pgf:
        result["loss_pgf"] = loss_pgf_certificate()
    if args.loss_entrance:
        result["loss_entrance"] = loss_entrance_certificate()
    if args.directed_cycle_loss:
        result["directed_cycle_loss"] = directed_cycle_loss_certificate()
        result["directed_cycle_loss"]["record_survival"] = directed_cycle_record_survival_certificate(
            result["directed_cycle_loss"]["two_site_all_parameter"])
    if args.three_cycle_loss:
        result["three_cycle_loss"] = three_cycle_loss_certificate()
    if args.output:
        saved = json.loads(args.output.read_text()) if args.output.exists() else {}
        saved.update(result)
        args.output.write_text(json.dumps(saved, indent=2)+"\n")
        print(json.dumps({"saved": str(args.output), "updated_certificates": list(result)}), flush=True)
    else:
        print(json.dumps(result, indent=2), flush=True)


if __name__ == "__main__":
    main()
