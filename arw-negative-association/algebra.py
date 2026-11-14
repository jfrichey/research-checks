#!/usr/bin/env python3
"""Small exact linear systems for stationary ARW and high-rate series.

Purpose: resolve the zero leading vacancy gaps at articulation vertices.
The IDLA multiplication matrices implement a quotient algebra, not an ARW
simulation. Checks against raw stabilization and its symbolic certificates
are therefore required before using the resulting higher coefficients.
"""

from exact import DomainMatrix, graph_kernel, sp


def idla_operators(kernel):
    """Column-stochastic operators: launch a walk to the first empty site."""
    n, size = len(kernel), 1 << len(kernel)
    operators = [sp.zeros(size) for _ in range(n)]
    transition = sp.zeros(n, n+1)
    for i, row in enumerate(kernel):
        for j, prob in row:
            transition[i, j] += prob
    symbolic_kernel = bool(transition.free_symbols)

    def canonical(matrix):
        return matrix.applyfunc(sp.cancel) if symbolic_kernel else matrix

    for mask in range(size):
        inside = [i for i in range(n) if mask >> i & 1]
        outside = [i for i in range(n+1) if i == n or not mask >> i & 1]
        if inside:
            laplacian = sp.eye(len(inside))-transition.extract(inside, inside)
            exits = transition.extract(inside, outside)
            harmonic = canonical(laplacian.inv()*exits)
            assert canonical(laplacian*harmonic-exits) == sp.zeros(*exits.shape)
            assert all(sp.cancel(sum(harmonic[i, :])-1) == 0 for i in range(len(inside)))
        for v in range(n):
            if v not in inside:
                operators[v][mask | (1 << v), mask] = 1
            else:
                row = inside.index(v)
                for j, w in enumerate(outside):
                    target = mask if w == n else mask | (1 << w)
                    operators[v][target, mask] += harmonic[row, j]
    jumps = []
    for v in range(n):
        jump = transition[v, n]*sp.eye(size)
        for w in range(n):
            jump += transition[v, w]*operators[w]
        jumps.append(canonical(jump))
        assert canonical(operators[v]*operators[v]-operators[v]*jump) == sp.zeros(size)
    for v in range(n):
        for w in range(v):
            assert canonical(operators[v]*operators[w]-operators[w]*operators[v]) == sp.zeros(size)
    return operators, jumps


def basis_matrix(operators, jumps, inverse_rates):
    n, size = len(operators), operators[0].rows
    empty = sp.zeros(size, 1)
    empty[0] = 1
    columns = [empty]
    for mask in range(1, size):
        v = (mask & -mask).bit_length()-1
        prev = columns[mask ^ (1 << v)]
        columns.append((operators[v]*prev
                        + inverse_rates[v]*(operators[v]-jumps[v])*prev).applyfunc(sp.expand))
    return sp.Matrix.hstack(*columns)


def stationary(kernel, inverse_rates):
    operators, jumps = idla_operators(kernel)
    matrix = basis_matrix(operators, jumps, inverse_rates)
    target = sp.zeros(matrix.rows, 1)
    target[-1] = 1
    aa = DomainMatrix.from_Matrix(matrix).to_field()
    bb = DomainMatrix.from_Matrix(target).convert_to(aa.domain)
    solution = aa.lu_solve(bb)
    assert (aa.matmul(solution)-bb).is_zero_matrix
    law = tuple(sp.factor(v) for v in solution.to_Matrix())
    assert sp.cancel(sum(law)-1) == 0
    return law


def stationary_orbit_polynomials(kernel, orbits, parameter):
    """Common-rate law using symmetry and fraction-free particle-count blocks.

    Return integer polynomials (weights, denominator). Orbits partition
    the stable masks and are intended to come from graph automorphisms.
    The final full 2^n-state residual verifies the reduction independently
    of that intended symmetry. No symbolic rational-function factorization
    is needed. This is useful for the specified six-site tree certificate.
    """
    n = len(kernel)
    orbits = sorted((tuple(sorted(orbit)) for orbit in orbits),
                    key=lambda orbit: (orbit[0].bit_count(), orbit[0]))
    assert sorted(mask for orbit in orbits for mask in orbit) == list(range(1 << n))
    assert all(len({mask.bit_count() for mask in orbit}) == 1 for orbit in orbits)
    additions, jumps = idla_operators(kernel)
    matrix = basis_matrix(additions, jumps, (parameter,)*n)
    zero, one = (sp.Poly(value, parameter, domain=sp.QQ) for value in (0, 1))
    reduced = [[sp.Poly(sum(matrix[row, source[0]] for row in target), parameter, domain=sp.QQ)
                for source in orbits] for target in orbits]
    denominator, numerators = one, [zero]*len(orbits)
    for rank in range(n, -1, -1):
        current = [i for i, orbit in enumerate(orbits) if orbit[0].bit_count() == rank]
        higher = [i for i, orbit in enumerate(orbits) if orbit[0].bit_count() > rank]
        size = len(current)
        augmented = [[reduced[i][j] for j in current]
                     +[(denominator if rank == n else zero)
                       -sum((reduced[i][j]*numerators[j] for j in higher), zero)]
                     for i in current]
        previous = one
        for k in range(size-1):
            pivot = augmented[k][k]
            assert not pivot.is_zero
            for i in range(k+1, size):
                for j in range(k+1, size+1):
                    augmented[i][j] = (pivot*augmented[i][j]
                                       -augmented[i][k]*augmented[k][j]).exquo(previous)
                augmented[i][k] = zero
            previous = pivot
        determinant = augmented[-1][-2]
        solution = [zero]*size
        for i in range(size-1, -1, -1):
            solution[i] = (determinant*augmented[i][-1]
                           -sum((augmented[i][j]*solution[j] for j in range(i+1, size)), zero)
                          ).exquo(augmented[i][i])
        for j in higher:
            numerators[j] *= determinant
        for i, j in enumerate(current):
            numerators[j] = solution[i]
        denominator *= determinant
        common = denominator
        for value in numerators:
            if not value.is_zero:
                common = sp.gcd(common, value)
        numerators = [value.exquo(common) for value in numerators]
        denominator = denominator.exquo(common)
    weights = [zero]*(1 << n)
    for orbit, value in zip(orbits, numerators):
        for mask in orbit:
            weights[mask] = value.exquo_ground(len(orbit))
    scale = sp.ilcm(*(value.q for polynomial in weights+[denominator]
                       for value in polynomial.all_coeffs()))
    if denominator.LC() < 0:
        scale = -scale
    weights = [sp.Poly(value*scale, parameter, domain=sp.ZZ) for value in weights]
    denominator = sp.Poly(denominator*scale, parameter, domain=sp.ZZ)
    content = sp.igcd(*(int(value) for polynomial in weights+[denominator]
                        for value in polynomial.all_coeffs()))
    weights = [value.exquo_ground(content) for value in weights]
    denominator = denominator.exquo_ground(content)
    assert sum(weights, sp.Poly(0, parameter)) == denominator
    for row in range(1 << n):
        residual = sum((sp.Poly(matrix[row, col], parameter)*weights[col]
                        for col in range(1 << n) if matrix[row, col] != 0),
                       sp.Poly(0, parameter))
        target = denominator if row == (1 << n)-1 else sp.Poly(0, parameter)
        assert (residual-target).is_zero
    return weights, denominator


def full_occupation_normalizer(kernel, inverse_rates):
    """Reciprocal full-occupation probability from principal Green minors.

    D(x)=sum_A product_{v in A} x_v / det(G_A), with G=(I-K)^(-1).
    The constant-term proof in report.md applies to any finite killed
    kernel, including nonreversible walks and self-loops.
    """
    n = len(kernel)
    transition = sp.zeros(n)
    for i, row in enumerate(kernel):
        for j, probability in row:
            if j < n:
                transition[i, j] += probability
    green = (sp.eye(n)-transition).inv()
    value = sp.S.One
    for mask in range(1, 1 << n):
        sites = [i for i in range(n) if mask >> i & 1]
        value += sp.prod(inverse_rates[i] for i in sites)/green.extract(sites, sites).det()
    return sp.factor(value)


def trace_kernel(kernel, removed):
    """Trace the killed walk over the removed vertices; keep remaining order."""
    n = len(kernel)
    removed = sorted(set(removed))
    assert all(0 <= v < n for v in removed)
    remaining = [v for v in range(n) if v not in removed]
    assert remaining
    transition = sp.zeros(n, n+1)
    for i, row in enumerate(kernel):
        for j, probability in row:
            transition[i, j] += probability
    destinations = remaining+[n]
    traced = transition.extract(remaining, destinations)
    if removed:
        traced += (transition.extract(remaining, removed)
                   *(sp.eye(len(removed))-transition.extract(removed, removed)).inv()
                   *transition.extract(removed, destinations))
    assert all(value >= 0 for value in traced)
    assert all(sum(traced[i, :]) == 1 for i in range(len(remaining)))
    return tuple(tuple((j, value) for j, value in enumerate(traced.row(i)) if value)
                 for i in range(len(remaining)))


def trace_comparison_certificate():
    """Finite diagnostics for the proposed conditional trace sandwich.

    If nu is the stationary trace law after removing v, the candidate is
    mu(outside | X_v=1) <=st nu <=st mu(outside | X_v=0).
    Also test the stronger upper bound mu(outside) <=st nu. These are not
    consequences of the proved full-occupation identity by themselves.
    """
    from bk import upsets

    records = []
    examples = (("conditional_base3", ((1, 3), (0, 2, 3), (1,))),
                ("endpoint_path4", ((1, 4), (0, 2), (1, 3), (2,))),
                ("conditional_tree5", ((1, 5), (0, 2, 3), (1,), (1, 4), (3,))))
    for name, neighbors in examples:
        n = len(neighbors)
        kernel = graph_kernel(neighbors)
        for rate in (1, 64):
            law = stationary(kernel, (sp.Rational(1, rate),)*n)
            failures = [[], [], []]
            comparisons = 0
            for v in range(n):
                remaining = [i for i in range(n) if i != v]
                nu = stationary(trace_kernel(kernel, [v]), (sp.Rational(1, rate),)*(n-1))
                parts = [[sp.S.Zero]*(1 << (n-1)) for _ in range(2)]
                for mask, mass in enumerate(law):
                    target = sum(((mask >> i) & 1) << j for j, i in enumerate(remaining))
                    parts[(mask >> v) & 1][target] += mass
                totals = tuple(map(sum, parts))
                differences = (tuple(nu[s]-parts[0][s]-parts[1][s] for s in range(len(nu))),
                               tuple(nu[s]-parts[1][s]/totals[1] for s in range(len(nu))),
                               tuple(parts[0][s]/totals[0]-nu[s] for s in range(len(nu))))
                for event in upsets(n-1):
                    comparisons += 1
                    for which, values in enumerate(differences):
                        gap = sum(values[s] for s in range(len(nu)) if event >> s & 1)
                        if gap < 0:
                            failures[which].append({"removed_vertex": v, "event_bits": event,
                                                    "gap": str(gap)})
            records.append({"case": name, "neighbors": neighbors, "lambda": rate,
                            "comparisons_for_each_inequality_including_constants": comparisons,
                            "unconditional_upper_bound_failures": failures[0],
                            "occupied_below_trace_failures": failures[1],
                            "empty_above_trace_failures": failures[2]})
    # Symbolic single-rate endpoints, independent of the event comparisons.
    x = sp.Symbol("x")
    kernel = graph_kernel(examples[0][1])
    endpoints = []
    for v in range(3):
        rates = [sp.Rational(1, 2), sp.Rational(1, 3), sp.Rational(1, 5)]
        rates[v] = x
        law = stationary(kernel, tuple(rates))
        remaining = [i for i in range(3) if i != v]
        nu = stationary(trace_kernel(kernel, [v]), tuple(rates[i] for i in remaining))
        for mask, mass in enumerate(law):
            target = sum(((mask >> i) & 1) << j for j, i in enumerate(remaining))
            assert mass.subs(x, 0) == (nu[target] if mask >> v & 1 else 0)
            numerator, denominator = map(lambda value: sp.Poly(value, x), sp.fraction(mass))
            assert numerator.degree() <= denominator.degree()
            limit = (numerator.LC()/denominator.LC()
                     if numerator.degree() == denominator.degree() else sp.S.Zero)
            assert limit == (0 if mask >> v & 1 else nu[target])
        endpoints.append({"removed_vertex": v, "both_rate_endpoints_verified": True})
    return {"scope": "specified finite cases; stochastic inequalities remain conjectural",
            "comparisons": records, "symbolic_endpoint_checks": endpoints}


def vertex_excursion_operators(kernel, inverse_rates, vertex):
    """Exact one-active-vertex insertion and its return-excursion resolvent.

    Columns start with vertex active and a stable outside configuration.
    B records outputs with vertex asleep; A records outputs with it empty.
    B=(I+x R)^(-1), A=x J B, R=I-E, and J=H R, as proved in report.md.
    H is one addition from the removed vertex into the traced ARW system.
    """
    rates = tuple(map(sp.Rational, inverse_rates))
    assert all(value > 0 for value in rates)
    n = len(kernel)
    outside = [i for i in range(n) if i != vertex]
    size = 1 << (n-1)
    empty = [sum(((mask >> j) & 1) << i for j, i in enumerate(outside))
             for mask in range(size)]
    occupied = [mask | (1 << vertex) for mask in empty]
    additions, jumps = idla_operators(kernel)
    mm = basis_matrix(additions, jumps, rates)
    nn = DomainMatrix.from_Matrix(mm).to_field().inv().to_Matrix()
    insertion = nn*additions[vertex]*mm
    aa = insertion.extract(empty, empty)
    bb = insertion.extract(occupied, empty)
    bi = DomainMatrix.from_Matrix(bb).to_field().inv().to_Matrix()
    rr = (bi-sp.eye(size))/rates[vertex]
    jj = aa*bi/rates[vertex]
    ee = sp.eye(size)-rr
    assert all(value >= 0 for value in ee) and all(value >= 0 for value in jj)
    assert all(sum(ee[:, col])+sum(jj[:, col]) == 1 for col in range(size))
    assert all(sum(jj[:, col]) > 0 for col in range(size))
    assert all(ee[row, col] == 0 for row in range(size) for col in range(size)
               if row.bit_count() > col.bit_count())
    traced = trace_kernel(kernel, [vertex])
    ta, tj = idla_operators(traced)
    tm = basis_matrix(ta, tj, tuple(rates[i] for i in outside))
    tn = DomainMatrix.from_Matrix(tm).to_field().inv().to_Matrix()
    nu = tn[:, -1]
    row = dict(kernel[vertex])
    hh = row.get(n, 0)*sp.eye(size)
    for j, site in enumerate(outside):
        hh += row.get(site, 0)*tn*ta[j]*tm
    hh /= 1-row.get(vertex, 0)
    assert hh*rr == jj and hh*nu == nu
    mu = nn[:, -1]
    a, b = mu.extract(empty, [0]), mu.extract(occupied, [0])
    assert aa*nu == a and bb*nu == b
    assert a == rates[vertex]*jj*b and a+hh*b == nu
    for threshold in range(n+1):
        original_tail = sum(mu[mask] for mask in range(1 << n) if mask.bit_count() >= threshold)
        lower = sum(nu[mask] for mask in range(size) if mask.bit_count() >= threshold)
        upper = sum(nu[mask] for mask in range(size) if mask.bit_count() >= threshold-1)
        assert lower <= original_tail <= upper
    return {"R": rr, "J": jj, "H": hh, "A": aa, "B": bb, "nu": nu,
            "outside_vertices": outside, "input_masks": empty,
            "empty_output_masks": empty, "occupied_output_masks": occupied}


def vertex_reduction_certificate():
    """Validate deletion identities, the generator signs, and exact insertion."""
    from exact import absorption, add, stable

    directed = (((0, sp.Rational(1, 4)), (1, sp.Rational(1, 2)), (3, sp.Rational(1, 4))),
                ((1, sp.Rational(1, 3)), (2, sp.Rational(1, 2)), (3, sp.Rational(1, 6))),
                ((0, sp.Rational(1, 3)), (1, sp.Rational(1, 3)), (3, sp.Rational(1, 3))))
    examples = (("endpoint_path3", graph_kernel(((1, 3), (0, 2), (1,))), (1, 1, 1)),
                ("conditional_base3", graph_kernel(((1, 3), (0, 2, 3), (1,))), (2, 3, 5)),
                ("directed_looped3", directed, (2, 3, 5)),
                ("conditional_tree5", graph_kernel(((1, 5), (0, 2, 3), (1,), (1, 4), (3,))),
                 (64,)*5))
    records = []
    for name, kernel, sleep_rates in examples:
        n = len(kernel)
        inverse = tuple(1/sp.Rational(value) for value in sleep_rates)
        for vertex in (range(n) if n == 3 else (0, 2)):
            result = vertex_excursion_operators(kernel, inverse, vertex)
            entry = {"case": name, "sleep_rates": sleep_rates, "vertex": vertex,
                     "kernel": [[[j, str(probability)] for j, probability in row] for row in kernel],
                     "insertion_reproduces_stationary_law": True,
                     "R_has_nonpositive_off_diagonal_and_positive_column_sums": True,
                     "return_excursions_do_not_increase_count": True,
                     "rank_interlacing_with_trace": True,
                     "trace_and_last_departure_identities": True}
            if n == 3:
                starts = [add(stable(mask, n), vertex) for mask in result["input_masks"]]
                probabilities = tuple(sp.Rational(rate)/(1+sp.Rational(rate)) for rate in sleep_rates)
                for reverse in (False, True):
                    raw, count = absorption(starts, kernel, probabilities, reverse=reverse)
                    assert raw[:, result["empty_output_masks"]].T == result["A"]
                    assert raw[:, result["occupied_output_masks"]].T == result["B"]
                entry["both_raw_toppling_orders_checked"] = True
                entry["raw_transient_states_last_order"] = count
                for key in ("R", "J", "H"):
                    entry[key] = [[str(value) for value in row] for row in result[key].tolist()]
            if name == "directed_looped3" and vertex == 0:
                changed = list(inverse)
                changed[vertex] *= 2
                other = vertex_excursion_operators(kernel, changed, vertex)
                assert all(other[key] == result[key] for key in ("R", "J", "H", "nu"))
                entry["excursion_operators_independent_of_inserted_sleep_rate"] = True
            records.append(entry)
    # Arbitrary NA input cannot replace the traced stationary input.
    kernel = graph_kernel(((1, 2), (0, 2)))
    raw, states = absorption([(1, -1)], kernel, (sp.Rational(1, 2),)*2)
    law = list(raw.row(0))
    covariance = law[3]-(law[1]+law[3])*(law[2]+law[3])
    assert covariance == sp.Rational(1333, 76050)
    return {"cases": records, "generic_insertion_counterexample": {
                "neighbors": ((1, 2), (0, 2)), "lambda": 1,
                "initial": "vertex 0 active, vertex 1 asleep",
                "law": list(map(str, law)), "covariance": str(covariance),
                "raw_transient_states": states}}


def excursion_tail_certificate():
    """Finite tests of the two proposed excursion-tail inequalities.

    Set w_k=E^k nu, s_k=sum(w_k). The candidate is
    w_k <=st H w_k <=st s_k nu, for equal-mass measures. If true for
    every k, geometric summation proves all three trace comparisons at
    every inserted-vertex rate. A bounded k check is not that proof.
    """
    from bk import upsets

    examples = (
        ("conditional_base3", ((1, 3), (0, 2, 3), (1,)), 0, 1),
        ("conditional_base3", ((1, 3), (0, 2, 3), (1,)), 1, 64),
        ("endpoint_path4", ((1, 4), (0, 2), (1, 3), (2,)), 0, 1),
        ("conditional_tree5", ((1, 5), (0, 2, 3), (1,), (1, 4), (3,)), 2, 64),
        ("conditional_tree5", ((1, 5), (0, 2, 3), (1,), (1, 4), (3,)), 3, 1))
    records = []
    for name, neighbors, vertex, rate in examples:
        n, last = len(neighbors), 12
        operators = vertex_excursion_operators(
            graph_kernel(neighbors), (sp.Rational(1, rate),)*n, vertex)
        nu, hh = operators["nu"], operators["H"]
        ee = sp.eye(nu.rows)-operators["R"]
        events = [(bits, tuple(s for s in range(nu.rows) if bits >> s & 1))
                  for bits in upsets(n-1)]
        probabilities = [sum(nu[s] for s in event) for _, event in events]
        vector = nu
        failures = [[], []]
        for k in range(last+1):
            terminal = hh*vector
            mass = sum(vector)
            assert sum(terminal) == mass
            for (bits, event), probability in zip(events, probabilities):
                current = sum(vector[s] for s in event)
                final = sum(terminal[s] for s in event)
                for which, gap in enumerate((final-current, mass*probability-final)):
                    if gap < 0:
                        failures[which].append({"returns": k, "event_bits": bits,
                                                "gap": str(gap)})
            vector = ee*vector
        records.append({"case": name, "neighbors": neighbors, "vertex": vertex,
                        "outside_sleep_rate": rate, "maximum_returns": last,
                        "comparisons_per_inequality_including_constants": (last+1)*len(events),
                        "continuation_increase_failures": failures[0],
                        "terminal_tail_domination_failures": failures[1]})

    # A generic stochastic-order invariant cannot prove these orbit claims.
    operators = vertex_excursion_operators(graph_kernel(examples[0][1]), (1,)*3, 0)
    nu = operators["nu"]
    ee = sp.eye(4)-operators["R"]
    pp = ee+operators["J"]
    event = (2, 3)  # Original vertex 2 occupied.
    gap = sum(pp[s, 3]-pp[s, 2] for s in event)
    assert gap == -sp.Rational(50443, 926640)
    probability = sum(nu[s] for s in event)
    test = sp.Matrix([int(s in event)-probability for s in range(4)])
    values = ee.T*test
    # Extreme downward transports of nu characterize its stochastic lower cone.
    transport = [max((s for s in range(4) if s & eta == s), key=lambda s: values[s])
                 for eta in range(4)]
    below = sp.zeros(4, 1)
    for eta, target in enumerate(transport):
        below[target] += nu[eta]
    excess = (test.T*ee*below)[0]
    assert excess == sp.Rational(21563, 5475600)
    return {"scope": "specified graphs, outside rates, vertices, and 0..12 returns only",
            "cases": records, "generic_order_obstructions": {
                "case": "conditional_base3", "vertex": 0, "outside_sleep_rate": 1,
                "outside_vertex_order": [1, 2], "event_masks": event,
                "P_equals_E_plus_J_monotonicity_gap_from_mask_2_to_3": str(gap),
                "downward_transport_target_masks": transport,
                "law_stochastically_below_nu": list(map(str, below)),
                "E_image_excess_over_its_mass_times_nu": str(excess)}}


def excursion_tail_symbolic_certificate():
    """All return counts on 0--1--2, sink adjacent to 0 and 1, deleting 0.

    The outside sleep rates are both 1/x, x>0; the inserted rate is arbitrary.
    For each increasing event, its two gap rows f satisfy
      f(I-zE)^(-1)nu = z*(alpha*D_{0,1}+beta*z*D_{0}+gamma*z**2)/D,
    where D_S=det(I-zE[S,S]) and D=D_{0,1,2,3}. The fixed decomposition
    has nonnegative rational coefficients for every x>0. For any E>=0,
    D_S/D=exp(sum_{r>=1} z**r*(tr(E**r)-tr(E[S,S]**r))/r) has
    nonnegative Taylor coefficients: the trace difference counts closed
    walks visiting outside S. This proves both inequalities for every k,
    without a return-count cutoff or a search over determinant subsets.
    """
    from exact import absorption, add, stable, stationary_law

    x, z = sp.symbols("x z", positive=True)
    neighbors = ((1, 3), (0, 2, 3), (1,))
    kernel = graph_kernel(neighbors)
    identity = sp.eye(4)
    empty, occupied = [0, 2, 4, 6], [1, 3, 5, 7]

    def inverse(matrix):
        return DomainMatrix.from_Matrix(matrix).to_field().inv().to_Matrix()

    def nonnegative_coefficients(value, strict=False):
        numerator, denominator = (sp.Poly(p, x, domain=sp.QQ)
                                  for p in sp.fraction(sp.cancel(value)))
        if denominator.LC() < 0:
            numerator, denominator = -numerator, -denominator
        return (all(c >= 0 for c in numerator.all_coeffs())
                and all(c >= 0 for c in denominator.all_coeffs())
                and (not strict or not numerator.is_zero))

    additions, jumps = idla_operators(kernel)
    matrix = basis_matrix(additions, jumps, (1, x, x))
    insertion = (inverse(matrix)*additions[0]*matrix).applyfunc(sp.cancel)
    aa, bb = insertion.extract(empty, empty), insertion.extract(occupied, empty)
    bi = inverse(bb)
    rr = (bi-identity).applyfunc(sp.factor)
    ee = (identity-rr).applyfunc(sp.factor)
    jj = (aa*bi).applyfunc(sp.factor)
    traced = trace_kernel(kernel, [0])
    ta, tj = idla_operators(traced)
    tm = basis_matrix(ta, tj, (x, x))
    tn = inverse(tm)
    nu = tn[:, -1].applyfunc(sp.factor)
    # The first jump from 0 goes to the sink or outside vertex 1 equally.
    hh = ((identity+tn*ta[0]*tm)/2).applyfunc(sp.factor)
    assert (hh*rr-jj).applyfunc(sp.cancel) == sp.zeros(4)
    assert (hh*nu-nu).applyfunc(sp.cancel) == sp.zeros(4, 1)
    assert sp.cancel(sum(nu)-1) == 0
    for addition in ta:
        assert (tn*addition*tm*nu-nu).applyfunc(sp.cancel) == sp.zeros(4, 1)
    for operator in (ee, jj, hh):
        assert all(nonnegative_coefficients(value) for value in operator)
    assert all(nonnegative_coefficients(value, strict=True) for value in nu)
    for col in range(4):
        assert sp.cancel(sum(hh[:, col])-1) == 0
        assert sp.cancel(sum(ee[:, col])+sum(jj[:, col])-1) == 0
        assert nonnegative_coefficients(sum(jj[:, col]), strict=True)
    assert all(ee[row, col] == 0 for row in range(4) for col in range(4)
               if row.bit_count() > col.bit_count())

    # Certify arbitrary inserted inverse rate y without inverting a matrix
    # over Q(x,y): M_y * [y J; I] = a_0 M_y[:,empty] * (I+y R).
    y = sp.Symbol("y", positive=True)
    my = basis_matrix(additions, jumps, (y, x, x))
    residual = (my[:, empty]*(y*jj)+my[:, occupied]
                -additions[0]*my[:, empty]*(identity+y*rr))
    assert residual.applyfunc(sp.cancel) == sp.zeros(8, 4)
    print("symbolic excursion operators: positivity, trace stationarity, and arbitrary inserted rate checked",
          flush=True)

    # Rank blocks give D and a degree-three vector t=adj(I-zE)nu.
    block = ee.extract([1, 2], [1, 2])
    qb = (sp.eye(2)-z*block).det()
    dd = (1-z*ee[0, 0])*(1-z*ee[3, 3])*qb
    middle = ((sp.eye(2)-z*block).adjugate()
              *((1-z*ee[3, 3])*nu.extract([1, 2], [0])
                +z*ee.extract([1, 2], [3])*nu[3]))
    numerator = sp.Matrix([
        (1-z*ee[3, 3])*qb*nu[0]+z*(ee[0, 1]*middle[0]+ee[0, 2]*middle[1])
        +z*ee[0, 3]*qb*nu[3],
        (1-z*ee[0, 0])*middle[0],
        (1-z*ee[0, 0])*middle[1],
        (1-z*ee[0, 0])*qb*nu[3]])
    tc = [numerator.applyfunc(lambda value: sp.factor(sp.expand(value).coeff(z, k)))
          for k in range(4)]
    dc = [sp.factor(sp.expand(dd).coeff(z, k)) for k in range(5)]
    determinant = (identity-z*ee).det(method="berkowitz")
    for k in range(5):
        assert sp.cancel(sp.expand(determinant).coeff(z, k)-dc[k]) == 0
        residual = ((tc[k] if k < 4 else sp.zeros(4, 1))
                    -(ee*tc[k-1] if k else sp.zeros(4, 1))-dc[k]*nu)
        assert residual.applyfunc(sp.cancel) == sp.zeros(4, 1)

    subset, singleton = (0, 1), 0
    sub = ee.extract(subset, subset)
    ds = (sp.eye(2)-z*sub).det()
    records = []
    for event in ((3,), (1, 3), (2, 3), (1, 2, 3)):
        indicator = sp.Matrix([[int(i in event) for i in range(4)]])
        probability = sum(nu[i] for i in event)
        for kind in ("left", "right"):
            ff = (indicator*(hh-identity) if kind == "left"
                  else probability*sp.ones(1, 4)-indicator*hh)
            cc = [sp.factor((ff*column)[0]) for column in tc]
            assert cc[0] == 0
            alpha = cc[1]
            beta = sp.factor(cc[2]+alpha*sp.trace(sub))
            gamma = sp.factor(cc[3]-alpha*sub.det()+beta*ee[singleton, singleton])
            coefficients = (alpha, beta, gamma)
            assert all(nonnegative_coefficients(value) for value in coefficients)
            assert all((value == 0) == (index == 2 and event == (3,) and kind == "left")
                       for index, value in enumerate(coefficients))
            decomposition = z*(alpha*ds+beta*z*(1-z*ee[singleton, singleton])+gamma*z**2)
            for k in range(4):
                assert sp.cancel(cc[k]-sp.expand(decomposition).coeff(z, k)) == 0
            records.append({"event": event, "kind": kind,
                            "coefficients": list(map(str, coefficients))})
    assert len(records) == 8
    print("excursion-tail certificate: eight fixed determinant decompositions and resolvent identity checked",
          flush=True)

    # Independent raw insertion, including recovery of E,J,H from A,B.
    expected = {key: value.subs(x, 1) for key, value in
                (("A", aa), ("B", bb), ("R", rr), ("E", ee), ("J", jj), ("H", hh))}
    starts = [add(stable(mask, 3), 0) for mask in empty]
    counts = []
    for reverse in (False, True):
        raw, count = absorption(starts, kernel, (sp.Rational(1, 2),)*3, reverse=reverse)
        raw_a, raw_b = raw[:, empty].T, raw[:, occupied].T
        raw_r = inverse(raw_b)-identity
        raw_j = raw_a*inverse(raw_b)
        raw_h = raw_j*inverse(raw_r)
        recovered = {"A": raw_a, "B": raw_b, "R": raw_r, "E": identity-raw_r,
                     "J": raw_j, "H": raw_h}
        assert recovered == expected
        counts.append(count)
    raw_nu, trace_counts = stationary_law(traced, (sp.Rational(1, 2),)*2)
    assert sp.Matrix(raw_nu) == nu.subs(x, 1)
    # Check all binary all-active inputs as well, for the mixture obstruction.
    binary_starts = [tuple((mask >> i) & 1 for i in range(2)) for mask in range(4)]
    for reverse in (False, True):
        raw_trace, _ = absorption(binary_starts, traced, (sp.Rational(1, 2),)*2,
                                  reverse=reverse)
        assert raw_trace.T == tn.subs(x, 1)
    mixture = inverse(raw_trace.T)*recovered["E"]*sp.Matrix(raw_nu)
    conjugate_entry = (recovered["H"]*recovered["E"]*inverse(recovered["H"]))[0, 2]
    assert mixture == sp.Matrix([sp.Rational(19, 810), -sp.Rational(311, 21060),
                                sp.Rational(1763, 21060), sp.Rational(7, 52)])
    assert conjugate_entry == -sp.Rational(8341, 1853280)
    print("raw x=1 insertion agrees in both orders:", counts,
          "; trace stationarity/all additions:", trace_counts, flush=True)
    return {
        "case": "conditional_base3", "neighbors": neighbors, "removed_vertex": 0,
        "parameter": "x=1/common outside sleep rate > 0; inserted sleep rate arbitrary > 0",
        "outside_vertices": [1, 2], "states": [[], [1], [2], [1, 2]],
        "scope": "E^k nu <=st H E^k nu <=st (1^T E^k nu) nu for every integer k>=0",
        "operators": {"E": [[str(v) for v in row] for row in ee.tolist()],
                      "H": [[str(v) for v in row] for row in hh.tolist()],
                      "nu": list(map(str, nu))},
        "decomposition": {
            "subset": subset, "singleton": singleton, "coefficient_order": ["alpha", "beta", "gamma"],
            "D_S": "det(I-z E[S,S]); D=det(I-z E)",
            "identity": "f(I-zE)^(-1)nu = z*(alpha*D_{0,1}+beta*z*D_{0}+gamma*z**2)/D",
            "left_row": "1_U^T(H-I)", "right_row": "nu(U)*1^T-1_U^T H",
            "positivity": "D_S/D has nonnegative Taylor coefficients for every nonnegative E",
            "records": records},
        "symbolic_checks": {
            "E_J_H_nonnegative": True, "E_plus_J_column_mass_one": True,
            "J_column_mass_strictly_positive": True, "H_column_stochastic": True,
            "H_nu_equals_nu": True, "trace_addition_equations": 2,
            "arbitrary_inserted_rate_identity": True,
            "determinant_and_resolvent_coefficients": 5, "gap_decompositions": 8,
            "triple_numerator_and_denominator_coefficients_nonnegative": True,
            "only_zero_coefficient": "gamma, left comparison, event {3}"},
        "raw_check": {
            "x": 1, "inserted_sleep_rate": 1,
            "insertion_transient_states_forward_reverse": counts,
            "trace_transient_states_forward_all_additions_reverse_saturated": trace_counts,
            "both_insertion_orders_and_recovered_E_J_H_checked": True,
            "trace_all_binary_inputs_both_orders_checked": True,
            "trace_law": list(map(str, raw_nu)),
            "trace_binary_input_mixture_coefficients_for_E_nu": list(map(str, mixture)),
            "H_E_H_inverse_entry_0_2": str(conjugate_entry)}}


def determinantal_obstruction():
    """Rule out even a complex L-ensemble on the four-site endpoint path.

    Pair and triple masses already give incompatible principal-minor
    equations. This is not a refutation of a determinant representation
    for the one-variable rank polynomial alone.
    """
    from itertools import combinations
    from exact import stationary_law

    neighbors = ((1, 4), (0, 2), (1, 3), (2,))
    kernel = graph_kernel(neighbors)
    law = stationary(kernel, (1,)*4)
    raw, states = stationary_law(kernel, (sp.Rational(1, 2),)*4)
    assert raw == law
    weights = [value/law[0] for value in law]
    pairs = {(i, j): weights[1 << i]*weights[1 << j]-weights[(1 << i) | (1 << j)]
             for i, j in combinations(range(4), 2)}
    assert all(value != 0 for value in pairs.values())
    triples = {tri: (sp.prod(weights[1 << i] for i in tri)
                     -sum(weights[1 << i]*pairs[tuple(j for j in tri if j != i)] for i in tri)
                     -weights[sum(1 << i for i in tri)])
               for tri in combinations(range(4), 3)}
    u, v, w, t = sp.symbols("u v w t")

    def primitive(expression, variable):
        return sp.Poly(expression, variable).clear_denoms()[1].primitive()[1].as_expr()

    quadratics = {}
    for (i, j), variable in (((1, 2), u), ((2, 3), v), ((1, 3), w)):
        quadratics[i, j] = primitive(pairs[0, j]*variable**2+triples[0, i, j]*variable
                                    +pairs[0, i]*pairs[i, j], variable)
    # Diagonal similarity permits L[0,i]=-1, L[i,0]=-pairs[0,i].
    # Then u=L[1,2], v=L[2,3], w=L[1,3], all nonzero, and t=u*v/w.
    polynomial = sp.resultant(quadratics[1, 2], u*v-t*w, u)
    polynomial = sp.resultant(quadratics[2, 3], polynomial, v)
    polynomial = sp.resultant(quadratics[1, 3], polynomial, w)
    polynomial = sp.Poly(primitive(polynomial, t), t)
    last = sp.Poly(primitive(pairs[1, 3]*t**2+triples[1, 2, 3]*t
                            +pairs[1, 2]*pairs[2, 3], t), t)
    assert polynomial.degree() == 8 and last.degree() == 2
    assert polynomial.gcd(last).degree() == 0
    prime = 101
    reduced = sp.Poly(polynomial, t, modulus=prime)
    reduced_last = sp.Poly(last, t, modulus=prime)
    assert reduced.degree() == 8 and reduced_last.degree() == 2
    ss, tt, gg = sp.gcdex(reduced, reduced_last)
    assert gg.as_expr() == 1 and ss*reduced+tt*reduced_last == gg
    return {"neighbors": neighbors, "lambda": 1, "law": list(map(str, law)),
            "raw_absorption_states": states, "both_orders_and_all_additions_checked": True,
            "all_pair_products_nonzero": True, "necessary_polynomial_degrees": [8, 2],
            "rational_polynomial_gcd_degree": 0,
            "modular_certificate": {"prime": prime,
                "polynomials_descending_coefficients": [list(map(int, reduced.all_coeffs())),
                                                         list(map(int, reduced_last.all_coeffs()))],
                "bezout_multipliers_descending_coefficients": [list(map(int, ss.all_coeffs())),
                                                                 list(map(int, tt.all_coeffs()))],
                "identity": "s(t)*F(t)+t_multiplier(t)*G(t)=1 modulo 101"},
            "conclusion": "No complex 4x4 matrix has principal minors mu(S)/mu(empty)."}


def nested_basis_entry(kernel, inverse_rates, target, source):
    """Closed entry M[target,source] when target is a subset of source.

    This uses the affine constant-term identity proved in report.md,
    independently of the IDLA multiplication-matrix construction.
    """
    assert target & ~source == 0
    n = len(kernel)
    transition, killing = sp.zeros(n), sp.zeros(n, 1)
    for i, row in enumerate(kernel):
        for j, probability in row:
            if j < n:
                transition[i, j] += probability
            else:
                killing[i] += probability
    laplacian = sp.eye(n)-transition
    inside = [i for i in range(n) if target >> i & 1]
    outside = [i for i in range(n) if source >> i & 1 and not target >> i & 1]
    result = sp.S.Zero
    for choice in range(1 << len(inside)):
        selected = [inside[j] for j in range(len(inside)) if choice >> j & 1]
        remaining = [i for i in inside if i not in selected]
        submatrix = laplacian.extract(remaining, remaining)
        harmonic = submatrix.inv()*killing.extract(remaining, [0]) if remaining else sp.zeros(0, 1)
        escape_product = sp.prod(killing[v]+sum(transition[v, w]*harmonic[j]
                                                for j, w in enumerate(remaining)) for v in outside)
        result += sp.prod(inverse_rates[v] for v in selected)*escape_product/submatrix.det()
    result *= (-1)**len(outside)*sp.prod(inverse_rates[v] for v in outside)*laplacian.extract(inside, inside).det()
    return sp.factor(result)


def homogenized_hessian_certificate(law):
    """Check quadratic derivatives of sum_B law[B] z_B s**(n-|B|).

    With full support this is the Lorentzian criterion. Real symmetric
    Hessians have only real eigenvalues, so Descartes' sign variations
    count their positive eigenvalues exactly, without numerical tolerances.
    A failure also rules out preservation from a Lorentzian input.
    """
    n = len(law).bit_length()-1
    checked, failures = 0, []
    for selected in range(1 << n):
        remaining = n-selected.bit_count()
        if remaining < 2:
            continue
        sites = [i for i in range(n) if not selected >> i & 1]
        hessian = sp.zeros(remaining+1)
        hessian[0, 0] = sp.factorial(remaining)*law[selected]
        for a, i in enumerate(sites, 1):
            hessian[0, a] = hessian[a, 0] = sp.factorial(remaining-1)*law[selected | (1 << i)]
            for b, j in enumerate(sites, 1):
                if i != j:
                    hessian[a, b] = sp.factorial(remaining-2)*law[selected | (1 << i) | (1 << j)]
        coefficients = hessian.charpoly().all_coeffs()
        signs = [sp.sign(value) for value in coefficients if value]
        positive = sum(a != b for a, b in zip(signs, signs[1:]))
        checked += 1
        if positive > 1:
            failures.append({"occupied_derivative_mask": selected,
                             "dummy_derivative_order": remaining-2,
                             "positive_eigenvalues": positive,
                             "characteristic_polynomial_coefficients": list(map(str, coefficients))})
    return {"quadratic_derivatives_checked": checked, "hessian_condition": not failures,
            "full_support": all(value > 0 for value in law), "failures": failures}


def high_rate_series(kernel, order, alpha=None):
    """Coefficients through t**order, with lambda_v=alpha_v/t."""
    n = len(kernel)
    alpha = (sp.S.One,)*n if alpha is None else tuple(map(sp.Rational, alpha))
    t = sp.Symbol("t")
    operators, jumps = idla_operators(kernel)
    matrix = basis_matrix(operators, jumps, tuple(t/a for a in alpha))
    terms = [matrix.applyfunc(lambda v: v.coeff(t, k)) for k in range(n+1)]
    assert terms[0] == sp.eye(matrix.rows)
    law = [sp.zeros(matrix.rows, 1)]
    law[0][-1] = 1
    for degree in range(1, order+1):
        value = sp.zeros(matrix.rows, 1)
        for k in range(1, min(n, degree)+1):
            value -= terms[k]*law[degree-k]
        assert sum(value) == 0
        law.append(value)
    return tuple(tuple(row) for row in law)


def cylinder_gap_series(law, left, right):
    """Vacancy-cylinder gap, given the coefficient vectors of a law."""
    assert not left & right
    def vacancy(mask):
        return [sum(vector[state] for state in range(len(vector)) if state & mask == 0)
                for vector in law]
    aa, bb, cc = vacancy(left), vacancy(right), vacancy(left | right)
    return [sum(aa[j]*bb[k-j] for j in range(k+1))-cc[k] for k in range(len(law))]


def wake_jump_operators(kernel, inverse_rates):
    """Rate-flow kernels: wake a sleeper, force its first jump, stabilize.

    Returns the all-active stabilization matrix and the stochastic kernels.
    Inputs here are rational so that nonnegativity can be checked exactly.
    """
    inverse_rates = tuple(map(sp.Rational, inverse_rates))
    operators, jumps = idla_operators(kernel)
    matrix = basis_matrix(operators, jumps, inverse_rates)
    stabilization = matrix.inv()
    n, size = len(kernel), matrix.rows
    result = []
    for v in range(n):
        remove, put = sp.zeros(size), sp.zeros(size)
        for mask in range(size):
            if mask >> v & 1:
                remove[mask ^ (1 << v), mask] = 1
            else:
                put[mask | (1 << v), mask] = 1
        jump = stabilization*jumps[v]*matrix
        wake = sp.eye(size)+(jump-put)*remove
        assert all(sum(wake[:, mask]) == 1 for mask in range(size))
        assert all(value >= 0 for value in wake)
        result.append(wake)
    return stabilization, tuple(result)


def diagnostics():
    """Check the rate identity and disprove generic NA preservation by N.

    These discriminate between an evolution of the stationary family and
    the stronger, false claim that stabilization preserves product inputs.
    """
    import json
    from pathlib import Path
    from exact import absorption, stationary_law

    # A different exact formula for the full-state mass, checked against
    # every saved four-site law and raw absorption on a directed looped
    # graph. The latter also discriminates between general killed kernels
    # and an identity restricted to reversible walks.
    saved = json.loads(Path(__file__).with_name("four_site_certificates.json").read_text())
    t = sp.Symbol("t")
    for entry in saved["cases"].values():
        normalizer = full_occupation_normalizer(graph_kernel(entry["neighbors"]), (t,)*4)
        full_mass = sp.sympify(entry["numerators"][-1])/sp.sympify(entry["denominator"])
        assert sp.cancel(full_mass*normalizer-1) == 0
    directed = (((0, sp.Rational(1, 4)), (1, sp.Rational(1, 2)), (3, sp.Rational(1, 4))),
                ((1, sp.Rational(1, 3)), (2, sp.Rational(1, 2)), (3, sp.Rational(1, 6))),
                ((0, sp.Rational(1, 3)), (1, sp.Rational(1, 3)), (3, sp.Rational(1, 3))))
    variables = sp.symbols("x0:3")
    normalizer = full_occupation_normalizer(directed, variables)
    additions, jumps = idla_operators(directed)
    directed_basis = basis_matrix(additions, jumps, variables)
    assert sp.expand(directed_basis[-1, -1]-normalizer) == 0
    nested_count, sign_count = 0, 0
    for target in range(8):
        for source in range(8):
            signed = (-1)**((source & ~target).bit_count())*directed_basis[target, source]
            assert all(coefficient >= 0 for coefficient in sp.Poly(signed, *variables).coeffs())
            sign_count += 1
            if target & ~source == 0:
                closed = nested_basis_entry(directed, variables, target, source)
                assert sp.expand(closed-directed_basis[target, source]) == 0
                nested_count += 1
    point = dict(zip(variables, (sp.Rational(1, 2), sp.Rational(1, 3), sp.Rational(1, 5))))
    direct_law, direct_count = stationary_law(directed, (sp.Rational(2, 3), sp.Rational(3, 4), sp.Rational(5, 6)))
    assert direct_law[-1]*normalizer.subs(point) == 1
    print("Green-minor full-occupation formula: 44 saved laws and directed looped raw absorption agree", flush=True)
    full_check = {"four_site_laws_checked": len(saved["cases"]),
                  "directed_kernel": [[[j, str(prob)] for j, prob in row] for row in directed],
                  "inverse_rate_variables": list(map(str, variables)),
                  "normalizer": str(normalizer), "sitewise_basis_identity_checked": True,
                  "nested_basis_entries_checked": nested_count,
                  "basis_coefficient_signs_checked": sign_count,
                  "raw_check_sleep_rates": [2, 3, 5], "raw_absorption_states": direct_count,
                  "raw_full_occupation_probability": str(direct_law[-1])}

    # The homogenized stationary polynomial being Lorentzian would imply
    # ULC of ranks, and hence NA on complete graphs. This is a distinct,
    # still conjectural property; it does not imply general-graph NA.
    hessian_cases = []
    five = json.loads(Path(__file__).with_name("five_site_certificates.json").read_text())
    for family, collection in (("four_site", saved["cases"]), ("five_site", five["cases"])):
        for label, entry in collection.items():
            weights = list(map(sp.sympify, entry["numerators"]))
            for rate in (1, 64):
                result = homogenized_hessian_certificate([value.subs(t, sp.Rational(1, rate)) for value in weights])
                result.update({"family": family, "case": label, "lambda": rate})
                hessian_cases.append(result)
    print("stationary homogenization Hessians:", sum(c["quadratic_derivatives_checked"] for c in hessian_cases),
          "checked; failing parameter cases", sum(not c["hessian_condition"] for c in hessian_cases), flush=True)

    case = json.loads(Path(__file__).with_name("sitewise_certificates.json").read_text())["case"]
    variables = sp.symbols("t0:3")
    point = dict(zip(variables, map(sp.Rational, (1, 2, 3))))
    law = sp.Matrix([sp.sympify(w)/sp.sympify(case["denominator"]) for w in case["numerators"]])
    matrix, wakes = wake_jump_operators(graph_kernel(case["neighbors"]), (1, 2, 3))
    evaluated = law.subs(point)
    assert evaluated == matrix[:, -1]
    for v in range(3):
        assert law.diff(variables[v]).subs(point) == (wakes[v]-sp.eye(8))*evaluated/(1+point[variables[v]])
    print("exact rate derivatives agree with all three wake-jump kernels", flush=True)

    kernel = graph_kernel(((1, 4), (0, 2), (1, 3), (2,)))
    operators, jumps = idla_operators(kernel)
    t = sp.Symbol("t")
    polynomial = basis_matrix(operators, jumps, (t,)*4)
    uniform = sp.ones(16, 1)/16
    derivative = -polynomial.applyfunc(lambda value: value.coeff(t, 1))*uniform
    marginals = [sum(derivative[mask] for mask in range(16) if mask >> i & 1) for i in (1, 2)]
    joint_derivative = sum(derivative[mask] for mask in range(16) if mask >> 1 & 1 and mask >> 2 & 1)
    slope = joint_derivative-sum(marginals)/2
    assert slope == sp.Rational(1, 192)

    stabilization = basis_matrix(operators, jumps, (sp.Rational(1, 64),)*4).inv()
    output = stabilization*uniform
    marginal = [sum(output[mask] for mask in range(16) if mask >> i & 1) for i in (1, 2)]
    joint = sum(output[mask] for mask in range(16) if mask >> 1 & 1 and mask >> 2 & 1)
    covariance = joint-sp.prod(marginal)
    assert covariance > 0
    starts = [tuple((mask >> i) & 1 for i in range(4)) for mask in range(16)]
    raw, count = absorption(starts, kernel, (sp.Rational(64, 65),)*4)
    assert raw.T == stabilization
    print("product active input loses NA; all 16 columns checked by raw absorption", flush=True)

    deterministic = stabilization[:, 3]
    tail_marginals = [sum(deterministic[mask] for mask in range(16) if mask >> i & 1) for i in (2, 3)]
    tail_joint = sum(deterministic[mask] for mask in range(16) if mask & 12 == 12)
    tail_covariance = tail_joint-sp.prod(tail_marginals)
    assert tail_covariance > 0
    initial = sp.zeros(16, 1)
    initial[3] = 1
    first_matrix = polynomial.applyfunc(lambda value: value.coeff(t, 1))
    second_matrix = polynomial.applyfunc(lambda value: value.coeff(t, 2))
    first = -first_matrix*initial
    second = -first_matrix*first-second_matrix*initial
    first_marginals = [sum(first[mask] for mask in range(16) if mask >> i & 1) for i in (2, 3)]
    second_joint = sum(second[mask] for mask in range(16) if mask & 12 == 12)
    tail_leading = second_joint-sp.prod(first_marginals)
    assert first_marginals == [sp.Rational(3, 4), 0] and tail_leading == sp.Rational(1, 16)
    deterministic_hessian = homogenized_hessian_certificate(tuple(deterministic))
    assert not deterministic_hessian["hessian_condition"]
    print("deterministic active input {0,1} also loses NA and the homogenized Lorentzian property", flush=True)
    return {"full_occupation_formula": full_check,
            "trace_comparison": trace_comparison_certificate(),
            "vertex_reduction": vertex_reduction_certificate(),
            "excursion_tails": excursion_tail_certificate(),
            "determinantal_obstruction": determinantal_obstruction(),
            "stationary_homogenization": {"scope": "specified finite graphs at lambda=1,64 only",
                                          "cases": hessian_cases},
            "deterministic_input_counterexample": {
                "graph": [[1, 4], [0, 2], [1, 3], [2]], "initial_active_sites": [0, 1],
                "pair": [2, 3], "parameter": "t=1/lambda", "leading_covariance": str(tail_leading),
                "first_marginal_coefficients": list(map(str, first_marginals)),
                "second_joint_coefficient": str(second_joint), "finite_witness_lambda": 64,
                "finite_witness_covariance": str(tail_covariance), "raw_column_checked": 3,
                "homogenized_hessian_check": deterministic_hessian},
            "rate_flow": {"graph": case["neighbors"], "inverse_rates": [1, 2, 3],
                           "exact_derivative_checks": 3, "stochastic_kernels_checked": True},
            "product_input_counterexample": {
                "graph": [[1, 4], [0, 2], [1, 3], [2]], "input": "independent Bernoulli(1/2) active particles",
                "pair": [1, 2], "parameter": "t=1/lambda", "leading_covariance": str(slope),
                "marginal_derivatives": list(map(str, marginals)), "joint_derivative": str(joint_derivative),
                "finite_witness_lambda": 64, "finite_witness_covariance": str(covariance),
                "raw_absorption_states": count, "raw_columns_checked": 16}}


def full_event_sector_certificate():
    """Test the sufficient one-vacancy source sign on four specified cases.

    The exact top-sector identity holds generally, but nonnegativity of
    its source vector is open. Reuse three saved, independently validated
    three-site operators; the saved branching tree tests loss of symmetry
    at the known conditional-correlation obstruction. No graph/rate sweep.
    """
    import json
    from pathlib import Path

    saved = json.loads(Path(__file__).with_name("algebra_certificates.json").read_text())
    records = []

    def check(label, vertex, rates, ee, hh, nu, outside, resolvent):
        size, m = nu.rows, len(outside)
        full = size-1
        holes = [full ^ (1 << i) for i in range(m)]
        assert hh*nu == nu and sum(nu) == 1
        assert all(ee[i, j] == 0 for i in range(size) for j in range(size)
                   if i.bit_count() > j.bit_count())
        assert all(hh[i, j] == 0 for i in range(size) for j in range(size)
                   if i.bit_count() > j.bit_count()+1)
        dd, bb = ee.extract(holes, holes), ee.extract(holes, [full])
        edge, vv, ww = ee[full, full], nu.extract(holes, [0]), nu[full]
        hh_top = hh.extract([full], holes)
        source = (dd-edge*sp.eye(m))*vv+bb*ww
        assert (hh_top*vv)[0] == (1-hh[full, full])*ww
        for k in (1, 2, 3):
            geometric = sum((edge**ell*dd**(k-1-ell) for ell in range(k)), sp.zeros(m))
            assert ((hh-sp.eye(size))*ee**k*nu)[full] == (hh_top*geometric*source)[0]
        assert all(value > 0 for value in source)
        record = {"case": label, "removed_vertex": vertex,
                  "sleep_rates": list(map(str, rates)), "outside_vertices": outside,
                  "source": list(map(str, source)), "source_strictly_positive": True,
                  "return_count_identities": [1, 2, 3]}
        if resolvent:
            z = sp.Symbol("z")
            reduced = sp.factor(z*(hh_top*(sp.eye(m)-z*dd).inv()*source)[0]/(1-z*edge))
            direct = ((hh-sp.eye(size))*(sp.eye(size)-z*ee).inv()*nu)[full]
            assert sp.cancel(direct-reduced) == 0
            record["full_left_resolvent"] = str(reduced)
            record["full_resolvent_identity_checked"] = True
        records.append(record)
        print(label, "one-vacancy source strictly positive; top-sector identities passed", flush=True)

    for label, vertex in (("conditional_base3", 0), ("endpoint_path3", 1),
                          ("directed_looped3", 0)):
        case = next(c for c in saved["vertex_reduction"]["cases"]
                    if c["case"] == label and c["vertex"] == vertex)
        kernel = tuple(tuple((j, sp.Rational(value)) for j, value in row)
                       for row in case["kernel"])
        rates = tuple(map(sp.Rational, case["sleep_rates"]))
        outside = [i for i in range(len(kernel)) if i != vertex]
        nu = sp.Matrix(stationary(trace_kernel(kernel, [vertex]),
                                  tuple(1/rates[i] for i in outside)))
        check(label, vertex, rates, sp.eye(nu.rows)-sp.Matrix(case["R"]),
              sp.Matrix(case["H"]), nu, outside, True)

    neighbors = ((1, 5), (0, 2, 3), (1,), (1, 4), (3,))
    operators = vertex_excursion_operators(graph_kernel(neighbors), (sp.Rational(1, 64),)*5, 2)
    check("conditional_tree5", 2, (64,)*5,
          sp.eye(16)-operators["R"], operators["H"], operators["nu"],
          operators["outside_vertices"], False)
    records[-1]["neighbors"] = [list(row) for row in neighbors]
    return {"scope": "general identity; source positivity proved only in the four specified exact cases",
            "cases": records}


def _field_idla_operators(transition):
    """IDLA harmonic construction in a rational-function field.

    Keeping arithmetic in DomainMatrix avoids repeated expansion of the
    weighted three-site audit. Harmonic, stochastic, and algebra identities
    are all checked, just as in idla_operators.
    """
    domain = transition.domain
    n, size = transition.shape[0], 1 << transition.shape[0]
    identity = DomainMatrix.eye(size, domain)
    arrays = [[[domain.zero for _ in range(size)] for _ in range(size)] for _ in range(n)]
    for mask in range(size):
        inside = [i for i in range(n) if mask >> i & 1]
        outside = [i for i in range(n+1) if i == n or not mask >> i & 1]
        if inside:
            left = DomainMatrix.eye(len(inside), domain)-transition.extract(inside, inside)
            exits = transition.extract(inside, outside)
            harmonic = left.inv().matmul(exits)
            assert (left.matmul(harmonic)-exits).is_zero_matrix
            ones = DomainMatrix([[domain.one] for _ in outside], (len(outside), 1), domain).to_sparse()
            expected = DomainMatrix([[domain.one] for _ in inside], (len(inside), 1), domain).to_sparse()
            assert (harmonic.matmul(ones)-expected).is_zero_matrix
        for i in range(n):
            if i not in inside:
                arrays[i][mask | (1 << i)][mask] = domain.one
            else:
                for j, site in enumerate(outside):
                    target = mask if site == n else mask | (1 << site)
                    arrays[i][target][mask] += harmonic[inside.index(i), j].element
    additions = [DomainMatrix(array, (size, size), domain).to_sparse() for array in arrays]
    jumps = []
    for i in range(n):
        jump = identity.scalarmul(transition[i, n].element)
        for j in range(n):
            jump += additions[j].scalarmul(transition[i, j].element)
        jumps.append(jump)
        assert (additions[i].matmul(additions[i])-additions[i].matmul(jump)).is_zero_matrix
        for j in range(i):
            assert (additions[i].matmul(additions[j])-additions[j].matmul(additions[i])).is_zero_matrix
    return additions, jumps


def _field_trace(transition, removed):
    n, domain = transition.shape[0], transition.domain
    keep = [i for i in range(n) if i not in removed]
    targets = keep+[n]
    harmonic = (DomainMatrix.eye(len(removed), domain)
                -transition.extract(removed, removed)).inv().matmul(transition.extract(removed, targets))
    return transition.extract(keep, targets)+transition.extract(keep, removed).matmul(harmonic)


def _field_basis_taylor(additions, jumps, variable_sites):
    """Coefficients through t^2, inverse rates t or one as specified."""
    n, domain = len(additions), additions[0].domain
    size = 1 << n
    zero = DomainMatrix.zeros((size, 1), domain)
    initial = DomainMatrix([[domain.one]]+[[domain.zero]]*(size-1), (size, 1), domain).to_sparse()
    columns = [[initial, zero, zero]]
    for mask in range(1, size):
        i = (mask & -mask).bit_length()-1
        previous = columns[mask ^ (1 << i)]
        difference = additions[i]-jumps[i]
        if i in variable_sites:
            column = [additions[i].matmul(previous[k])
                      +(difference.matmul(previous[k-1]) if k else zero) for k in range(3)]
        else:
            column = [(additions[i]+difference).matmul(previous[k]) for k in range(3)]
        columns.append(column)
    return [columns[0][k].hstack(*(columns[mask][k] for mask in range(1, size))) for k in range(3)]


def full_event_dominated_coefficients(case):
    """The two reversible three-site traces needed for dominated sites.

    Conductances: root--i=u, root--j=v, root--sink=1; w is i--j
    for case 'same', and j--sink for case 'separate'. Only the source
    at i is asserted to have nonnegative mixed second coefficient.
    """
    u, v, w = sp.symbols("u v w", positive=True)
    domain = sp.QQ.frac_field(u, v, w)
    if case == "same":
        transition = sp.Matrix([[0, u/(1+u+v), v/(1+u+v), 1/(1+u+v)],
                                [u/(u+w), 0, w/(u+w), 0],
                                [v/(v+w), w/(v+w), 0, 0]])
    else:
        assert case == "separate"
        transition = sp.Matrix([[0, u/(1+u+v), v/(1+u+v), 1/(1+u+v)],
                                [1, 0, 0, 0], [v/(v+w), 0, 0, w/(v+w)]])
    transition = DomainMatrix.from_Matrix(transition).convert_to(domain)
    additions, jumps = _field_idla_operators(transition)
    matrix = _field_basis_taylor(additions, jumps, {1, 2})
    empty, occupied = [0, 2, 4, 6], [1, 3, 5, 7]
    pp = [a.extract(empty, empty) for a in matrix]
    qq = [a.extract(empty, occupied) for a in matrix]
    rr = [a.extract(occupied, empty) for a in matrix]
    ss = [a.extract(occupied, occupied) for a in matrix]
    identity = DomainMatrix.eye(4, domain)
    assert (pp[0]-identity).is_zero_matrix and rr[0].is_zero_matrix
    root = additions[0].extract(occupied, occupied)
    f1, f2 = pp[1]+root.matmul(rr[1]), pp[2]+root.matmul(rr[2])
    z1 = ss[1]-rr[1].matmul(qq[0])
    z2 = ss[2]-rr[1].matmul(qq[1])-rr[2].matmul(qq[0])+rr[1].matmul(pp[1]).matmul(qq[0])
    # B^-1=(P+A_root,oo R)^-1(S-R P^-1 Q); E=2I-B^-1.
    ee = [identity.scalarmul(domain(2))-ss[0], f1.matmul(ss[0])-z1,
          -z2+f1.matmul(z1)-(f1.matmul(f1)-f2).matmul(ss[0])]
    ta, tj = _field_idla_operators(_field_trace(transition, [0]))
    tm = _field_basis_taylor(ta, tj, {0, 1})
    assert (tm[0]-identity).is_zero_matrix
    rows = list(range(4))
    nu = [identity.extract(rows, [3]), -tm[1].extract(rows, [3]),
          (tm[1].matmul(tm[1])-tm[2]).extract(rows, [3])]
    holes, sources = [2, 1], []
    for order in (1, 2):
        value = DomainMatrix.zeros((2, 1), domain)
        for j in range(order+1):
            value += (ee[j].extract(holes, holes)
                      -DomainMatrix.eye(2, domain).scalarmul(ee[j][3, 3].element)).matmul(
                          nu[order-j].extract(holes, [0]))
            value += ee[j].extract(holes, [3]).scalarmul(nu[order-j][3, 0].element)
        sources.append(value)
    assert sources[0][0, 0].element == 0
    trace = _field_trace(transition, [2]).to_Matrix()
    c, d, a, b = trace[0, 0], trace[1, 1], trace[0, 1], trace[1, 0]
    assert sp.cancel(1-d-b) == 0
    self_term = sp.factor(a*b*(1-d-a*b/(1-c)))
    mixed = sp.factor(sources[1].to_Matrix()[0]-self_term)
    if case == "same":
        positive = (u**4*w*(v+w)**2
                    +u**3*(v+w)*(v*(v**2+3*v*w+4*w**2)+2*w*(v+w))
                    +u**2*w*(v**2*(4*v**2+10*v*w+7*w**2)
                             +2*v*(2*v+3*w)*(v+w)+(v+w)**2)
                    +u*v*w**2*(v**2*(5*v+6*w)+8*v*(v+w)+3*v+2*w)
                    +2*v**2*w**3*(v+1)**2)
        expected = u*positive/((u+w)**2*(v+w)**2*(u+v+1)**3*(u*v+u*w+u+v*w+w))
    else:
        positive = sp.S.One
        expected = u*v**4/((v+1)*(v+w)**2*(u+v+1)**3)
    assert sp.cancel(mixed-expected) == 0
    numerator, denominator = sp.fraction(mixed)
    assert all(value >= 0 for part in (numerator, denominator)
               for value in sp.Poly(part, u, v, w).coeffs())
    check_point = {u: 1, v: 1, w: 0 if case == "same" else 1}
    assert mixed.subs(check_point) == (sp.Rational(1, 54) if case == "same" else sp.Rational(1, 216))
    print(case, "dominated-site mixed coefficient: symbolic positive factorization verified", flush=True)
    return {"case": case, "parameters": "u,v,w>0; boundary conductances by continuity",
            "dominated_site": 1, "linear_coefficient": "0", "self_coefficient": str(self_term),
            "mixed_coefficient": str(expected), "positive_polynomial": str(positive),
            "symbolic_identity_verified": True, "all_numerator_denominator_coefficients_nonnegative": True,
            "harmonic_stochastic_and_IDLA_relations_checked_in_rational_field": True,
            "purposeful_rational_check": str(mixed.subs(check_point))}


def _raw_root_excursion_absorption(kernel, sleep, starts, reverse=False):
    """Raw E/J absorption, forcing one original jump from root vertex zero.

    Stabilize the outside before every further root jump. Stop with a
    stable outside when the root pool is one (E) or zero (J). A root
    self-loop at the first jump is an E return, not an erased launch.
    Sleep probability one is allowed: lone particles sleep immediately,
    while surplus particles still jump after ineffective sleeps are skipped.
    Return E above J, as columns, plus the transient-state count. This
    construction uses no IDLA change of basis or insertion identity.
    """
    n = len(kernel)
    assert all(sp.cancel(sum(probability for _, probability in row)-1) == 0
               for row in kernel)

    def jump(state, site):
        row = {}
        for destination, probability in kernel[site]:
            after = list(state)
            after[site] -= 1
            if destination < n:
                after[destination] = 2 if after[destination] == -1 else after[destination]+1
            after = tuple(after)
            row[after] = row.get(after, 0)+probability
        return row

    def terminal(state):
        if any(value > 0 for value in state[1:]) or state[0] > 1:
            return None
        assert state[0] in (0, 1)
        mask = sum(1 << (i-1) for i in range(1, n) if state[i] == -1)
        return mask+(0 if state[0] else 1 << (n-1))

    states, index = [], {}

    def discover(state):
        if terminal(state) is None and state not in index:
            index[state] = len(states)
            states.append(state)

    assert all(state[0] == 1 and all(value <= 0 for value in state[1:]) for state in starts)
    launched = [jump(state, 0) for state in starts]
    for row in launched:
        for state in row:
            discover(state)
    rows = []
    for state in states:
        order = range(n-1, 0, -1) if reverse else range(1, n)
        site = next((i for i in order if state[i] > 0), 0)
        if site == 0 or state[site] > 1:
            row = jump(state, site)
        else:
            row = {after: (1-sleep[site])*value for after, value in jump(state, site).items()}
            after = list(state)
            after[site] = -1
            after = tuple(after)
            row[after] = row.get(after, 0)+sleep[site]
            row = {after: value for after, value in row.items() if value != 0}
        assert sp.cancel(sum(row.values())-1) == 0
        rows.append(row)
        for after in row:
            discover(after)
    left = sp.eye(len(states))
    right = sp.zeros(len(states), 1 << n)
    for i, row in enumerate(rows):
        for after, probability in row.items():
            out = terminal(after)
            if out is None:
                left[i, index[after]] -= probability
            else:
                right[i, out] += probability
    field = DomainMatrix.from_Matrix(left).to_field()
    source = DomainMatrix.from_Matrix(right).convert_to(field.domain)
    solution = field.inv().matmul(source)
    assert (field.matmul(solution)-source).is_zero_matrix
    solved = solution.to_Matrix()
    assert all(sp.cancel(sum(solved.row(i))-1) == 0 for i in range(solved.rows))
    result = sp.zeros(1 << n, len(starts))
    for col, row in enumerate(launched):
        for after, probability in row.items():
            out = terminal(after)
            if out is None:
                result[:, col] += probability*solved[index[after], :].T
            else:
                result[out, col] += probability
    result = result.applyfunc(sp.factor)
    assert all(sp.cancel(sum(result[:, col])-1) == 0 for col in range(result.cols))
    return result, len(states)


def full_event_raw_highrate_audit():
    """One symbolic two-site family and one permanent-site trace audit.

    The latter fixes a reversible weighted three-site walk with loops at
    every site, keeps site 2 occupied at infinite rate, and leaves the
    inverse rate x at site 1 symbolic. Both legal outside priorities must
    give the full E/J trace, including the root's affine loop contribution.
    These checks audit the tracing step, not a finite-size inference.
    """
    a, b, x = sp.symbols("a b x", positive=True)

    def equal(left, right):
        assert (left-right).applyfunc(sp.cancel) == sp.zeros(*left.shape)

    g, denominator = a*b, 1+(1-a*b)*x
    expected = sp.Matrix([
        [g*x/(1+x), g*x*(2-b-g+(2-b)*(1-g)*x)/((1+x)*denominator)],
        [0, g/denominator]])
    raw, pair_states = _raw_root_excursion_absorption(
        (((1, a), (2, 1-a)), ((0, b), (2, 1-b))),
        (0, 1/(1+x)), [(1, 0), (1, -1)])
    equal(raw[:2, :], expected)
    law = sp.Matrix([x*(1-g), 1])/denominator
    source = (raw[:2, :]*law)[0]-raw[1, 1]*law[0]
    formula = g*x*((1-b)+(1-g)*x)/((1+x)*denominator)
    assert sp.cancel(source-formula) == 0

    kernel = (
        ((0, sp.Rational(1, 4)), (1, sp.Rational(1, 4)),
         (2, sp.Rational(1, 4)), (3, sp.Rational(1, 4))),
        ((0, sp.Rational(1, 3)), (1, sp.Rational(1, 6)),
         (2, sp.Rational(1, 3)), (3, sp.Rational(1, 6))),
        ((0, sp.Rational(1, 4)), (1, sp.Rational(1, 4)),
         (2, sp.Rational(1, 4)), (3, sp.Rational(1, 4))))
    traced = trace_kernel(kernel, [2])
    assert traced == (
        ((0, sp.Rational(1, 3)), (1, sp.Rational(1, 3)), (2, sp.Rational(1, 3))),
        ((0, sp.Rational(4, 9)), (1, sp.Rational(5, 18)), (2, sp.Rational(5, 18))))
    pair, traced_states = _raw_root_excursion_absorption(
        traced, (0, 1/(1+x)), [(1, 0), (1, -1)])
    orders = []
    for reverse in (False, True):
        full, count = _raw_root_excursion_absorption(
            kernel, (0, 1/(1+x), sp.S.One), [(1, 0, -1), (1, -1, -1)], reverse)
        assert full.extract([0, 1, 4, 5], [0, 1]) == sp.zeros(4, 2)
        equal(full.extract([2, 3, 6, 7], [0, 1]), pair)
        orders.append({"order": "rightmost" if reverse else "leftmost", "transient_states": count})
    c, d = sp.Rational(1, 3), sp.Rational(5, 18)
    aa, bb = sp.Rational(1, 3), sp.Rational(4, 9)
    erased = raw.subs({a: aa/(1-c), b: bb/(1-d), x: (1-d)*x}, simultaneous=True)
    equal(pair[:2, :], c*sp.eye(2)+(1-c)*erased[:2, :])
    equal(pair[2:, :], (1-c)*erased[2:, :])
    killing = 1-d-bb
    trace_killing = 1-d-aa*bb/(1-c)
    trace_law = sp.Matrix([x*trace_killing, 1])/(1+x*trace_killing)
    loop_source = sp.factor((pair[:2, :]*trace_law)[0]-pair[1, 1]*trace_law[0])
    assert sp.diff(loop_source, x).subs(x, 0) == aa*bb*killing/(1-d) == sp.Rational(20, 351)
    assert loop_source.subs(x, 1) == sp.Rational(184, 3627)
    print("Raw two-site source and permanent-site E/J trace passed, including both outside orders and root loops.",
          flush=True)
    return {
        "scope": "symbolic two-site family and one fixed looped permanent-site trace; no graph/rate sweep",
        "parameter": "x=1/lambda_1>0; site 2 occupied at infinite rate",
        "raw_transient_equations_and_normalization_checked": True,
        "loopless_symbolic_source_checked": True,
        "loopless_raw_transient_states": pair_states,
        "kernel": [[[j, str(probability)] for j, probability in row] for row in kernel],
        "traced_kernel": [[[j, str(probability)] for j, probability in row] for row in traced],
        "traced_pair_raw_transient_states": traced_states,
        "permanent_site_orders": orders,
        "permanent_site_occupation_and_full_E_J_trace_checked": True,
        "root_loop_E_affine_and_J_scaling_checked": True,
        "looped_source": str(loop_source),
        "looped_linear_coefficient": "20/351",
    }


def full_event_highrate_certificate():
    """Symbolic two-site identity behind the general first-order source.

    This checks every a,b,x symbolically and both weighted three-site
    mixed coefficient families, plus independent raw two-site and
    permanent-site tracing equations. The general-graph extension uses the
    tracing/analyticity argument in the report, not finite graph inference.
    """
    a, b, x = sp.symbols("a b x", positive=True)
    kernel = (((1, a), (2, 1-a)), ((0, b), (2, 1-b)))
    additions, jumps = idla_operators(kernel)
    matrix = basis_matrix(additions, jumps, (1, x))
    inverse = DomainMatrix.from_Matrix(matrix).to_field().inv().to_Matrix()
    insertion = (inverse*additions[0]*matrix).applyfunc(sp.cancel)
    sleeping = insertion.extract([1, 3], [0, 2])
    excursion = (2*sp.eye(2)-sleeping.inv()).applyfunc(sp.factor)
    g = a*b
    expected = sp.Matrix([
        [g*x/(1+x), g*x*(2-b-g+(2-b)*(1-g)*x)/((1+x)*(1+(1-g)*x))],
        [0, g/(1+(1-g)*x)]])
    assert (excursion-expected).applyfunc(sp.cancel) == sp.zeros(2)
    law = sp.Matrix([x*(1-g), 1])/(1+x*(1-g))
    # trace_kernel deliberately asserts numerical probability inequalities;
    # its symbolic two-site trace here is the explicit singleton loop ab.
    trace_law = sp.Matrix(stationary((((0, g), (1, 1-g)),), (x,)))
    assert (law-trace_law).applyfunc(sp.cancel) == sp.zeros(2, 1)
    source = sp.factor((excursion[0, 0]-excursion[1, 1])*law[0]
                       +excursion[0, 1]*law[1])
    formula = a*b*x*((1-b)+(1-a*b)*x)/((1+x)*(1+(1-a*b)*x))
    assert sp.cancel(source-formula) == 0
    assert sp.cancel(sp.diff(source, x).subs(x, 0)-a*b*(1-b)) == 0
    assert sp.cancel(source.subs(b, 1)-a*(1-a)*x**2/((1+x)*(1+(1-a)*x))) == 0

    c, d = sp.symbols("c d", nonnegative=True)
    loop_source = (1-c)*formula.subs({a:a/(1-c), b:b/(1-d), x:(1-d)*x}, simultaneous=True)
    linear = sp.factor(sp.diff(loop_source, x).subs(x, 0))
    assert sp.cancel(linear-a*b*(1-d-b)/(1-d)) == 0
    print("Two-site excursion/source and generic loop-erasure coefficient verified.", flush=True)
    raw_audit = full_event_raw_highrate_audit()
    dominated = [full_event_dominated_coefficients(case) for case in ("same", "separate")]
    return {"scope": "symbolic identities for the general undirected high-outside-rate full-event source proof",
            "assumptions_without_loops": "0<a,b<=1, ab<1, x>0",
            "excursion": [[str(sp.factor(value)) for value in expected.row(i)] for i in range(2)],
            "trace_law": list(map(str, law)), "source": str(sp.factor(formula)),
            "source_identity_checked": True, "strictly_positive_for_ab_positive": True,
            "loop_erasure_linear_coefficient": str(linear),
            "raw_two_site_and_permanent_site_audit": raw_audit,
            "dominated_two_site_source": str(sp.factor(source.subs(b, 1))),
            "dominated_site_mixed_coefficients": dominated,
            "general_dominated_site_quadratic_status": "self term strictly positive; all mixed terms nonnegative by the two traced families"}


def source_reversible_components_check(transition, inverse_rates, block, forcing, scaling):
    """Audit reversible trace components and their mixing on the saved triangle.

    The symbolic decomposition uses all seven retained subsets of the existing
    graph. The single x1=x2=1 boundary specialization tests whether the mixing
    correction has a favorable sign; it does not. No graph or rate sweep.
    """
    n = 3
    pi = list(map(sp.Integer, (6, 7, 8)))
    laplacian = sp.eye(n)-transition[:, :n]

    def determinant(sites):
        return laplacian.extract(sites, sites).det()

    def trace(sites):
        removed = [i for i in range(n) if i not in sites]
        result = transition.extract(sites, sites+[n])
        if removed:
            result += (transition.extract(sites, removed)
                       *laplacian.extract(removed, removed).inv()
                       *transition.extract(removed, sites+[n]))
        return result

    minors = [determinant([j for j in range(n) if j != i]) for i in range(n)]
    killed = scaling*block*scaling.inv()
    source = scaling*forcing
    generator = sp.zeros(n+1)
    generator[0, 0] = -sum(source)
    generator[1:, 0] = source
    generator[0, 1:] = sp.ones(1, n)*killed
    generator[1:, 1:] = -killed
    component_sum = sp.zeros(n+1)
    for mask in range(1, 1 << n):
        sites = [i for i in range(n) if mask >> i & 1]
        traced = trace(sites)
        scale = sp.prod(inverse_rates[i] for i in sites) / determinant(
            [i for i in range(n) if i not in sites])
        component = sp.zeros(n+1)
        weights = sp.ones(n+1, 1)
        for a, i in enumerate(sites):
            escape = 1-traced[a, a]
            weights[i+1] = inverse_rates[i]*pi[i]*escape/minors[i]
            component[0, i+1] = sp.cancel(scale*minors[i]*traced[a, -1]
                                          /(inverse_rates[i]*escape))
            component[i+1, 0] = scale*pi[i]*traced[a, -1]
            for b, j in enumerate(sites):
                if i != j:
                    component[j+1, i+1] = sp.cancel(scale*minors[i]*traced[a, b]
                                                  /(inverse_rates[i]*escape))
        for j in range(n+1):
            component[j, j] = -sum(component[:, j])
        symmetric = component*sp.diag(*weights)
        assert (symmetric-symmetric.T).applyfunc(sp.cancel) == sp.zeros(n+1)
        component_sum += component
    assert (component_sum-generator).applyfunc(sp.cancel) == sp.zeros(n+1)

    boundary = dict(zip(inverse_rates, (0, 1, 1)))
    killed0 = killed.subs(boundary)
    odds0 = killed0.inv()*source.subs(boundary)
    derivative = killed0.inv()*(source.diff(inverse_rates[0]).subs(boundary)
                               -killed.diff(inverse_rates[0]).subs(boundary)*odds0)
    assert odds0[0] == 0
    green = killed0[1:, 1:].inv()
    baseline = sp.Matrix([1, odds0[1], odds0[2]])
    records = []
    separate_force = sp.zeros(3, 1)
    for mask in range(4):
        selected = [j+1 for j in range(2) if mask >> j & 1]
        retained = [0]+selected
        traced = trace(retained)
        escapes = [1-traced[j, j] for j in range(len(retained))]
        # a_S includes the root in V\S; gamma cancels its Schur factor.
        scale = 1/determinant([j for j in range(n) if j not in selected])
        gamma = scale*minors[0]*escapes[0]
        exits = sp.Matrix([traced[0, -1]/escapes[0]]
                          +[traced[0, retained.index(j)]/escapes[0]
                            if j in selected else 0 for j in (1, 2)])
        weights = {0: sp.S.One}
        weights.update({j: pi[j]*escapes[retained.index(j)]/minors[j]
                        for j in selected})
        ratio = sp.Matrix([baseline[j]/weights[j] if j in weights else 0
                           for j in range(3)])
        incoming = pi[0]*escapes[0]/minors[0]*(exits.T*ratio)[0]
        conductance = sp.zeros(3)
        if selected:
            smaller = trace(selected)
            for a, i in enumerate(selected):
                value = pi[i]*smaller[a, -1]
                conductance[i, i] += value
                conductance[0, 0] += value
                conductance[i, 0] -= value
                conductance[0, i] -= value
                for b in range(a+1, len(selected)):
                    j = selected[b]
                    value = pi[i]*smaller[a, b]
                    conductance[i, i] += value
                    conductance[j, j] += value
                    conductance[i, j] -= value
                    conductance[j, i] -= value
        separate_force -= scale*escapes[0]*conductance*ratio
        records.append((gamma, exits, incoming))
    total = sum(row[0] for row in records)
    mean_exits = sum((g*h for g, h, k in records), sp.zeros(3, 1))/total
    mean_incoming = sum(g*k for g, h, k in records)/total
    mixing = -sum((g*(h-mean_exits)*(k-mean_incoming) for g, h, k in records),
                  sp.zeros(3, 1))
    assert sum(separate_force) == sum(mixing) == 0
    rows = []
    for a, target in enumerate((1, 2)):
        row = sp.Matrix([[0, green[a, 0], green[a, 1]]])
        direct, correction = (row*separate_force)[0], (row*mixing)[0]
        assert direct+correction == derivative[target] > 0
        assert direct > 0 and correction < 0
        rows.append(dict(target=target, dirichlet=str(direct), mixing=str(correction),
                         actual_scaled_odds_derivative=str(derivative[target]),
                         green_exit_means=[str((row*h)[0]) for g, h, k in records]))
    probabilities = [g/total for g, h, k in records]
    gap = probabilities[0]*probabilities[3]-probabilities[1]*probabilities[2]
    assert probabilities == [sp.Rational(k, 149) for k in (56, 31, 31, 31)]
    assert gap == sp.Rational(775, 22201)
    return dict(scope="same weighted triangle; symbolic component identity and one root-boundary point",
                all_seven_component_reversibility_and_sum_identities=True,
                root=0, outside_inverse_rates=[1, 1],
                subset_order=[[], [1], [2], [1, 2]],
                subset_probabilities=list(map(str, probabilities)), FKG_square_gap=str(gap),
                incoming_ratios=[str(k) for g, h, k in records], rows=rows,
                all_derivative_identities=True,
                consequence="Mixing correction can be negative; FKG of retained-set weights alone does not prove the source sign.")


def one_vacancy_block_certificate():
    """Audit the common scaling and a proposed marginal-bound supersolution.

    The fixed weighted triangle has an asymmetric reversible kernel and
    killing only at vertex zero. All inverse rates remain independent.
    The general proof is the trace/forest argument in the report; this
    finite calculation checks its coefficient formula against the IDLA basis.
    """
    inverse_rates = sp.symbols("x0:3", positive=True)
    kernel = (
        ((1, sp.Rational(1, 3)), (2, sp.Rational(1, 2)), (3, sp.Rational(1, 6))),
        ((0, sp.Rational(2, 7)), (2, sp.Rational(5, 7))),
        ((0, sp.Rational(3, 8)), (1, sp.Rational(5, 8))))
    transition = sp.zeros(3, 4)
    for i, row in enumerate(kernel):
        for j, probability in row:
            transition[i, j] = probability
    laplacian = sp.eye(3)-transition[:, :3]
    determinant = lambda indices: laplacian.extract(indices, indices).det()
    minors = [determinant([j for j in range(3) if j != i]) for i in range(3)]
    scaling = sp.diag(*(sp.Rational(weight)/minor
                        for weight, minor in zip((6, 7, 8), minors)))
    additions, jumps = idla_operators(kernel)
    matrix = basis_matrix(additions, jumps, inverse_rates)
    holes = [6, 5, 3]
    block = matrix.extract(holes, holes)
    records = []
    for mask in range(8):
        selected = [i for i in range(3) if mask >> i & 1]
        remaining = [i for i in range(3) if i not in selected]
        monomial = sp.prod(inverse_rates[i] for i in selected)
        coefficient = block.applyfunc(
            lambda value: sp.Poly(value, *inverse_rates).coeff_monomial(monomial))
        expected = sp.zeros(3)
        for j in remaining:
            retained = selected+[j]
            removed = [i for i in range(3) if i not in retained]
            targets = retained+[3]
            trace = transition.extract(retained, targets)
            if removed:
                trace += (transition.extract(retained, removed)
                          *laplacian.extract(removed, removed).inv()
                          *transition.extract(removed, targets))
            expected[j, j] = minors[j]/determinant([i for i in remaining if i != j])
            for i in selected:
                expected[i, j] = -minors[i]*trace[retained.index(i), retained.index(j)]/determinant(remaining)
        assert coefficient == expected
        scaled = scaling*coefficient*scaling.inv()
        assert all(scaled[i, j] <= 0 for i in range(3) for j in range(3) if i != j)
        assert all(sum(scaled[:, j]) >= 0 for j in range(3))
        records.append({"selected_sites": selected,
                        "coefficient": [[str(value) for value in coefficient.row(i)] for i in range(3)],
                        "scaled_column_sums": [str(sum(scaled[:, j])) for j in range(3)]})

    def nonnegative_coefficients(value):
        return all(c >= 0 for c in sp.Poly(sp.expand(value), *inverse_rates).coeffs())

    assert nonnegative_coefficients(block.det())
    scaled_adjugate = scaling*block.adjugate()*scaling.inv()
    assert all(nonnegative_coefficients(value) for value in scaled_adjugate)
    forcing = -matrix.extract(holes, [7])
    assert all(nonnegative_coefficients(value) for value in forcing)

    # If this vector were a supersolution, inverse positivity would give
    # every single-site full-occupation marginal trace bound. It is not.
    full_denominator = matrix[7, 7]
    all_one = dict.fromkeys(inverse_rates, sp.S.One)
    proposed = sp.Matrix([
        (inverse_rates[i]*sp.diff(full_denominator, inverse_rates[i])
         /full_denominator.subs(inverse_rates[i], 0)).subs(all_one)
        for i in range(3)])
    residual = block.subs(all_one)*proposed-forcing.subs(all_one)
    assert residual[0] == -sp.Rational(105705905740999, 1777113292980480)
    boundary = dict(all_one)
    boundary[inverse_rates[0]] = 0
    block0 = block.subs(boundary)
    odds0 = block0.inv()*forcing.subs(boundary)
    all_zero = dict.fromkeys(inverse_rates, sp.S.Zero)
    marked_block = block.diff(inverse_rates[0]).subs(all_zero)
    marked_forcing = forcing.diff(inverse_rates[0]).subs(all_zero)
    marked_derivative = block0.inv()*(marked_forcing-marked_block*odds0)
    full_derivative = block0.inv()*(forcing.diff(inverse_rates[0]).subs(boundary)
                                    -block.diff(inverse_rates[0]).subs(boundary)*odds0)
    assert list(marked_derivative[1:]) == [
        -sp.Rational(973451, 88692399), -sp.Rational(2694595, 384333729)]
    assert all(value > 0 for value in full_derivative[1:])
    # Audit r_i=nu_full*d/dx_root[mu(V\i)/mu(V)] at x_root=0
    # against three previously raw-validated excursion kernels.
    import json
    from pathlib import Path
    saved = json.loads(Path(__file__).with_name("algebra_certificates.json").read_text())
    derivative_checks = []
    root_inverse = sp.symbols("root_inverse", nonnegative=True)
    for label, root in (("conditional_base3", 0), ("endpoint_path3", 1),
                        ("directed_looped3", 0)):
        case = next(c for c in saved["vertex_reduction"]["cases"]
                    if c["case"] == label and c["vertex"] == root)
        case_kernel = tuple(tuple((j, sp.Rational(value)) for j, value in row)
                            for row in case["kernel"])
        rates = tuple(map(sp.Rational, case["sleep_rates"]))
        outside = [i for i in range(3) if i != root]
        inverse = [root_inverse if i == root else 1/rates[i] for i in range(3)]
        aa, jj = idla_operators(case_kernel)
        mm = basis_matrix(aa, jj, inverse)
        hh, bb = mm.extract(holes, holes), -mm.extract(holes, [7])
        hh0 = hh.subs(root_inverse, 0)
        odds0 = hh0.inv()*bb.subs(root_inverse, 0)
        odds_derivative = hh0.inv()*(bb.diff(root_inverse)-hh.diff(root_inverse)*odds0)
        nu = sp.Matrix(stationary(trace_kernel(case_kernel, [root]),
                                  [1/rates[i] for i in outside]))
        ee = sp.eye(4)-sp.Matrix(case["R"])
        outside_holes = [2, 1]
        source = ((ee.extract(outside_holes, outside_holes)-ee[3, 3]*sp.eye(2))
                  *nu.extract(outside_holes, [0])+ee.extract(outside_holes, [3])*nu[3])
        assert (source-nu[3]*odds_derivative.extract(outside, [0])).applyfunc(sp.cancel) == sp.zeros(2, 1)
        assert sp.cancel(ee[3, 3]-1+mm[7, 7].diff(root_inverse)
                         /mm[7, 7].subs(root_inverse, 0)) == 0
        derivative_checks.append({"case": label, "root": root,
                                  "source": list(map(str, source)),
                                  "odds_derivative_and_diagonal_identities_checked": True})

    # The stronger condition that the total derivative forcing is positive
    # even before the trace Green inverse is false. This pendant geometry
    # was chosen from the mixed high-rate coefficient, not from a sweep.
    ray = sp.symbols("t", positive=True)
    forcing_kernel = (((1, sp.Rational(1, 3)), (2, sp.Rational(1, 3)), (3, sp.Rational(1, 3))),
                      ((0, sp.S.One),), ((0, sp.Rational(1, 2)), (3, sp.Rational(1, 2))))
    fa, fj = idla_operators(forcing_kernel)
    fm = basis_matrix(fa, fj, (root_inverse, ray, 10*ray))
    fh, fb = fm.extract(holes, holes), -fm.extract(holes, [7])
    fh0, fh1 = fh.subs(root_inverse, 0), fh.diff(root_inverse)
    fb0, fb1 = fb.subs(root_inverse, 0), fb.diff(root_inverse)
    fw0 = (fh0.inv()*fb0).applyfunc(sp.cancel)
    fraw = (fb1-fh1*fw0).applyfunc(sp.cancel)
    theta = sp.cancel(fraw[0]/fh0[0, 0])
    total = (fraw[1:, :]-fh0[1:, 0]*theta).applyfunc(sp.factor)
    trace_block = fh0[1:, 1:]
    ray_source = (trace_block.inv()*total/fm[7, 7].subs(root_inverse, 0)).applyfunc(sp.factor)
    expected_first = (ray**2*(112500*ray**4+188500*ray**3+12550*ray**2-1329*ray-51)
                      /(1620*(ray+1)*(10*ray+1)*(5*ray**2+9*ray+1)))
    assert sp.cancel(total[0]-expected_first) == 0
    assert total[0].subs(ray, sp.Rational(1, 100)) == -sp.Rational(502763, 157016296800)
    assert sp.limit(total[0]/ray**2, ray, 0) == -sp.Rational(17, 540)
    assert sp.limit(ray_source[0]/ray**2, ray, 0) == sp.Rational(133, 540)
    assert sp.limit(ray_source[1]/ray, ray, 0) == sp.Rational(5, 6)
    for value in ray_source:
        assert all(c > 0 for part in sp.fraction(value) for c in sp.Poly(part, ray).coeffs())
    print("One-vacancy coefficient trace formulas, scaling and forest signs passed; proposed supersolution fails.",
          flush=True)
    return {"scope": "one weighted triangle, independent positive inverse rates; general proof is in report",
            "edge_conductances": {"0-1": 2, "0-2": 3, "1-2": 5, "0-sink": 1},
            "scaling": [str(scaling[i, i]) for i in range(3)],
            "subset_coefficients": records,
            "determinant_adjugate_and_forcing_coefficients_nonnegative": True,
            "source_derivative_checks": derivative_checks,
            "reversible_component_decomposition": source_reversible_components_check(
                transition, inverse_rates, block, forcing, scaling),
            "total_root_forcing_obstruction": {
                "edges": ["0-1", "0-2", "0-sink", "2-sink"],
                "root": 0, "outside_inverse_rates": ["t", "10t"],
                "trace_block": [[str(value) for value in trace_block.row(i)] for i in range(2)],
                "total_forcing": list(map(str, total)), "source": list(map(str, ray_source)),
                "first_forcing_at_t_1_100": "-502763/157016296800",
                "leading_forcing_and_sources": ["-17/540*t^2", "133/540*t^2", "5/6*t"],
                "source_positive_at_every_t": True,
                "status": "total eliminated forcing can be negative; its final trace Green average is positive"},
            "separate_selected_subset_obstruction": {
                "selected_subset": [0], "root_inverse_rate": 0,
                "outside_inverse_rates": [1, 1],
                "marked_derivative": list(map(str, marked_derivative)),
                "full_derivative": list(map(str, full_derivative)),
                "status": "individual selected-set derivatives can be negative; total source sign is not refuted"},
            "failed_supersolution": {
                "definition": "g_i=x_i*(partial_i D)/(D at x_i=0); candidate H*g>=b",
                "inverse_rates": [1, 1, 1],
                "residual": list(map(str, residual)),
                "status": "sufficient supersolution fails; does not refute the marginal trace bound"}}


def full_event_triangle_certificate():
    """All-rate source signs on a weighted triangle killed only at the root.

    Uses the proved one-vacancy trace formulas, Cramer's derivative and
    coefficient positivity in u,v,w,x,y. Optional python-flint accelerates
    exact polynomial arithmetic; no numerical sweep or approximation enters.
    """
    import hashlib
    from itertools import permutations
    from flint import fmpq_mpoly_ctx, fmpq
    u, v, w = sp.symbols('u v w', positive=True)
    wv = (u, v, w)
    K = sp.Matrix([[0, u / (1 + u + v), v / (1 + u + v), 1 / (1 + u + v)], [u / (u + w), 0, w / (u + w), 0], [v / (v + w), w / (v + w), 0, 0]])
    L = sp.eye(3) - K[:, :3]
    dets = {mask: sp.factor(L.extract([i for i in range(3) if mask >> i & 1], [i for i in range(3) if mask >> i & 1]).det()) for mask in range(8)}
    minor = [dets[7 ^ 1 << i] for i in range(3)]
    traces = {}
    for mask in range(1, 8):
        keep = [i for i in range(3) if mask >> i & 1]
        removed = [i for i in range(3) if i not in keep]
        q = K.extract(keep, keep + [3])
        if removed:
            q += K.extract(keep, removed) * L.extract(removed, removed).inv() * K.extract(removed, keep + [3])
        traces[mask] = (keep, q.applyfunc(sp.factor))
    hcoef = []
    bcoef = []
    dcoef = []
    for mask in range(8):
        hu = sp.zeros(3)
        bu = sp.zeros(3, 1)
        selected = [i for i in range(3) if mask >> i & 1]
        remain = 7 ^ mask
        for j in range(3):
            if j not in selected:
                keep, q = traces[mask | 1 << j]
                hu[j, j] = minor[j] / dets[remain ^ 1 << j]
                for i in selected:
                    hu[i, j] = -minor[i] * q[keep.index(i), keep.index(j)] / dets[remain]
            else:
                keep, q = traces[mask]
                bu[j] = minor[j] * q[keep.index(j), -1] / dets[remain]
        hcoef.append(hu.applyfunc(sp.factor))
        bcoef.append(bu.applyfunc(sp.factor))
        dcoef.append(sp.factor(dets[7] / dets[remain]))
    all_coeff = [z for a in hcoef + bcoef for z in a] + dcoef
    common = sp.lcm([sp.fraction(z)[1] for z in all_coeff])
    common = sp.Poly(common, *wv)
    assert all((c > 0 for c in common.coeffs()))
    print('Positive graph denominator cleared for the triangle source calculation.', flush=True)
    ctx = fmpq_mpoly_ctx.get(('u', 'v', 'w', 'x', 'y'))
    one = ctx.constant(1)
    zero = ctx.constant(0)
    gens = ctx.gens()

    def poly(z):
        p = sp.Poly(sp.cancel(z * common.as_expr()), *wv)
        return ctx.from_dict({m + (0, 0): fmpq(str(c)) for m, c in p.terms()})
    H = [[[zero for j in range(3)] for i in range(3)] for k in range(2)]
    b = [[zero for i in range(3)] for k in range(2)]
    F = zero
    for mask in range(8):
        mono = gens[3] ** int(bool(mask & 2)) * gens[4] ** int(bool(mask & 4))
        k = mask & 1
        for i in range(3):
            b[k][i] += poly(bcoef[mask][i]) * mono
            for j in range(3):
                H[k][i][j] += poly(hcoef[mask][i, j]) * mono
        if not k:
            F += poly(dcoef[mask]) * mono
    lpoly = ctx.from_dict({m + (0, 0): fmpq(str(c)) for m, c in common.terms()})
    perms = [(p, (-1) ** sum((p[i] > p[j] for i in range(3) for j in range(i + 1, 3)))) for p in permutations(range(3))]

    def det(a):
        out = zero
        for p, sign in perms:
            out += sign * a[0][p[0]] * a[1][p[1]] * a[2][p[2]]
        return out

    def replace(a, col, values):
        return [[values[i] if j == col else a[i][j] for j in range(3)] for i in range(3)]

    def deriv(a0, a1):
        return sum((det(replace(a0, j, [a1[i][j] for i in range(3)])) for j in range(3)), zero)
    D0 = det(H[0])
    D1 = deriv(*H)
    rows = []
    for i in (1, 2):
        a0 = replace(H[0], i, b[0])
        a1 = replace(H[1], i, b[1])
        n = (deriv(a0, a1) * D0 - det(a0) * D1) * lpoly
        d = D0 ** 2 * F
        g = n.gcd(d)
        n = n // g
        d = d // g
        if list(d.to_dict().values())[0] < 0:
            n = -n
            d = -d
        assert all((c > 0 for p in [n, d] for c in p.to_dict().values()))
        row = {'hole': i, 'numerator_terms': len(n), 'denominator_terms': len(d), 'all_coefficients_positive': True, 'numerator_sha256': hashlib.sha256(str(n).encode()).hexdigest(), 'denominator_sha256': hashlib.sha256(str(d).encode()).hexdigest()}
        print('Triangle source', i, ':', len(n), 'positive numerator coefficients;', len(d), 'positive denominator coefficients.', flush=True)
        rows.append(row)

        def evaluate(p):
            from fractions import Fraction
            return sum((Fraction(str(c)) * 2 ** int(m[0]) * 3 ** int(m[1]) * 5 ** int(m[2]) for m, c in p.to_dict().items()), Fraction(0))
        row['at_weights_2_3_5_rates_1'] = str(evaluate(n) / evaluate(d))
    from exact import absorption
    fixed_kernel = (((1, sp.Rational(1, 3)), (2, sp.Rational(1, 2)), (3, sp.Rational(1, 6))), ((0, sp.Rational(2, 7)), (2, sp.Rational(5, 7))), ((0, sp.Rational(3, 8)), (1, sp.Rational(5, 8))))
    starts = [(1, -int(bool(mask & 1)), -int(bool(mask & 2))) for mask in range(4)]
    raw0, count0 = _raw_root_excursion_absorption(fixed_kernel, (0, sp.Rational(1, 2), sp.Rational(1, 2)), starts)
    raw1, count1 = _raw_root_excursion_absorption(fixed_kernel, (0, sp.Rational(1, 2), sp.Rational(1, 2)), starts, True)
    assert raw0 == raw1
    traced = trace_kernel(fixed_kernel, [0])
    trace_raw, trace_count = absorption([(1, 1)], traced, (sp.Rational(1, 2),) * 2)
    nu = trace_raw[0, :].T
    assert nu == sp.Matrix(stationary(traced, (1, 1)))
    excursion = raw0[:4, :]
    holes = [2, 1]
    source = (excursion.extract(holes, holes) - excursion[3, 3] * sp.eye(2)) * nu.extract(holes, [0]) + excursion.extract(holes, [3]) * nu[3]
    for i, record in enumerate(rows):
        assert source[i] == sp.Rational(record['at_weights_2_3_5_rates_1'])
        assert record['numerator_terms'] == 4053 and record['denominator_terms'] == 10242
    print('Triangle source matches raw excursions in both outside orders.', flush=True)
    return {'scope': 'all positive u,v,w,x,y; weighted triangle killed only at root 0', 'conductances': {'0-1': 'u', '0-2': 'v', '1-2': 'w', '0-sink': '1'}, 'outside_inverse_rates': ['x', 'y'], 'root_sleep_rate': 'arbitrary positive finite; does not enter the source', 'method': 'audited trace coefficient formulas; fraction-free Cramer derivative; every polynomial coefficient checked', 'positive_graph_denominator': str(sp.factor(common.as_expr())), 'sources': rows, 'raw_check': {'weights': [2, 3, 5], 'outside_inverse_rates': [1, 1], 'both_orders_equal': True, 'transient_states': [count0, count1], 'trace_raw_states': trace_count, 'sources': list(map(str, source))}, 'consequence': 'strict left all-return-count full-event comparison and strict marginal trace bound; not full NA'}


def full_event_tree_certificate():
    """Four independent outside rates on the saved branching tree at root one.

This targeted optional certificate uses the audited H_U,b_U trace formulas,
integer FLINT polynomials, fraction-free Cramer derivatives, and every
coefficient sign. One full-IDLA rational derivative checks all four sources.
No full excursion inverse, graph sweep, or default diagnostic is run.
"""
    from fractions import Fraction
    import hashlib
    import json
    from flint import fmpz_mpoly_ctx
    neighbors = ((1, 5), (0, 2, 3), (1,), (1, 4), (3,))
    n = 5
    root = 1
    outside = [0, 2, 3, 4]
    full = (1 << n) - 1
    names = tuple((f'x{i}' for i in outside))
    K = sp.zeros(n, n + 1)
    for i, row in enumerate(neighbors):
        for j in row:
            K[i, j] = sp.Rational(1, len(row))
    L = sp.eye(n) - K[:, :n]

    def sites(mask):
        return [i for i in range(n) if mask >> i & 1]
    dets = {mask: L.extract(sites(mask), sites(mask)).det() for mask in range(1 << n)}
    assert all((v > 0 for v in dets.values()))
    minor = [dets[full ^ 1 << i] for i in range(n)]
    traces = {}
    for mask in range(1, 1 << n):
        keep = sites(mask)
        removed = sites(full ^ mask)
        q = K.extract(keep, keep + [n])
        if removed:
            q += K.extract(keep, removed) * L.extract(removed, removed).inv() * K.extract(removed, keep + [n])
        assert all((sum(q.row(j)) == 1 for j in range(len(keep))))
        traces[mask] = (keep, q)
    hcoef = []
    bcoef = []
    dcoef = []
    for mask in range(1 << n):
        hu = sp.zeros(n)
        bu = sp.zeros(n, 1)
        U = sites(mask)
        remain = full ^ mask
        for j in range(n):
            if j not in U:
                keep, q = traces[mask | 1 << j]
                hu[j, j] = minor[j] / dets[remain ^ 1 << j]
                for i in U:
                    hu[i, j] = -minor[i] * q[keep.index(i), keep.index(j)] / dets[remain]
            else:
                keep, q = traces[mask]
                bu[j] = minor[j] * q[keep.index(j), -1] / dets[remain]
        hcoef.append(hu)
        bcoef.append(bu)
        dcoef.append(dets[full] / dets[remain])
    common = int(sp.ilcm(*[int(v.q) for mat in hcoef + bcoef for v in mat], *[int(v.q) for v in dcoef]))
    assert common > 0
    print(f'trace coefficients constructed; positive graph denominator {common}', flush=True)
    ctx = fmpz_mpoly_ctx.get(names)
    gens = ctx.gens()
    zero = ctx.constant(0)
    one = ctx.constant(1)
    H = [[[zero for j in range(n)] for i in range(n)] for order in range(2)]
    b = [[zero for i in range(n)] for order in range(2)]
    normalizer = zero
    for mask in range(1 << n):
        mon = one
        for k, i in enumerate(outside):
            if mask >> i & 1:
                mon *= gens[k]
        order = mask >> root & 1
        for i in range(n):
            b[order][i] += int(common * bcoef[mask][i]) * mon
            for j in range(n):
                H[order][i][j] += int(common * hcoef[mask][i, j]) * mon
        if not order:
            normalizer += int(common * dcoef[mask]) * mon
    assert all((c > 0 for c in normalizer.coeffs()))

    def determinant(matrix):
        size = len(matrix)
        if not size:
            return one
        aa = [list(row) for row in matrix]
        previous = one
        sign = 1
        for k in range(size - 1):
            pivotrow = next((i for i in range(k, size) if aa[i][k]), None)
            if pivotrow is None:
                return zero
            if pivotrow != k:
                aa[k], aa[pivotrow] = (aa[pivotrow], aa[k])
                sign = -sign
            pivot = aa[k][k]
            for i in range(k + 1, size):
                for j in range(k + 1, size):
                    numerator = pivot * aa[i][j] - aa[i][k] * aa[k][j]
                    aa[i][j] = numerator // previous
                aa[i][k] = zero
            previous = pivot
        return sign * aa[-1][-1]

    def replace(matrix, col, values):
        return [[values[i] if j == col else matrix[i][j] for j in range(n)] for i in range(n)]

    def determinant_derivative(a0, a1):
        return sum((determinant(replace(a0, j, [a1[i][j] for i in range(n)])) for j in range(n) if any((a1[i][j] for i in range(n)))), zero)

    def evaluate(poly, point):
        return sum((int(c) * sp.prod((Fraction(point[k]) ** int(e) for k, e in enumerate(exponents))) for exponents, c in poly.to_dict().items()), Fraction(0))

    def digest(poly):
        terms = sorted(((tuple(map(int, m)), str(c)) for m, c in poly.to_dict().items()))
        return hashlib.sha256(json.dumps(terms, separators=(',', ':')).encode()).hexdigest()

    def polynomial_record(poly):
        terms = poly.to_dict()
        negatives = sorted(((tuple(map(int, m)), int(c)) for m, c in terms.items() if c < 0), key=lambda z: (sum(z[0]), z[0]))
        return {'terms': len(terms), 'total_degree': int(poly.total_degree()), 'individual_degrees': list(map(int, poly.degrees())), 'minimum_total_degree': int(min(map(sum, terms))), 'minimum_exponents': [int(min((m[j] for m in terms))) for j in range(len(names))], 'negative_terms': len(negatives), 'all_nonzero_coefficients_positive': not negatives, 'first_negative': {'exponents': list(negatives[0][0]), 'coefficient': str(negatives[0][1])} if negatives else None, 'canonical_terms_sha256': digest(poly)}
    D0 = determinant(H[0])
    D1 = determinant_derivative(H[0], H[1])
    assert D0 and all((c > 0 for c in D0.coeffs()))
    assert all((c > 0 for c in D1.coeffs()))
    print(f'determinant terms: D0={len(D0)}, D1={len(D1)}', flush=True)
    positive_denominator = D0 ** 2 * normalizer
    point = (1, 2, 3, 5)
    records = []
    evaluations = {}
    reduced = []
    for i in outside:
        a0 = replace(H[0], i, b[0])
        a1 = replace(H[1], i, b[1])
        N0 = determinant(a0)
        N1 = determinant_derivative(a0, a1)
        numerator = common * (N1 * D0 - N0 * D1)
        denominator = positive_denominator
        print(f'hole {i}: before gcd numerator={len(numerator)}, denominator={len(denominator)}', flush=True)
        g = numerator.gcd(denominator)
        nn = numerator // g
        dd = denominator // g
        if evaluate(dd, (1, 1, 1, 1)) < 0:
            nn = -nn
            dd = -dd
        assert nn * denominator == dd * numerator
        nr = polynomial_record(nn)
        dr = polynomial_record(dd)
        assert dr['all_nonzero_coefficients_positive']
        assert nr['minimum_exponents'][outside.index(i)] >= 1
        value = evaluate(nn, point) / evaluate(dd, point)
        record = {'hole': i, 'numerator': nr, 'denominator': dr, 'at_outside_inverse_rates_1_2_3_5': str(value)}
        records.append(record)
        evaluations[i] = value
        reduced.append((nn, dd))
        print(f'hole {i}: reduced numerator={len(nn)}, denominator={len(dd)}, negative numerator terms={nr['negative_terms']}', flush=True)
    kernel = tuple((tuple(((j, K[i, j]) for j in range(n + 1) if K[i, j])) for i in range(n)))
    assert all((record['numerator']['all_nonzero_coefficients_positive'] for record in records))
    assert [(record['numerator']['terms'], record['denominator']['terms']) for record in records] == [(3584, 5832), (3571, 5832), (3127, 5184), (3077, 5184)]
    A, J = idla_operators(kernel)
    inverses = {i: sp.Rational(point[k]) for k, i in enumerate(outside)}
    empty = sp.eye(1 << n)[:, 0]
    nil = sp.zeros(1 << n, 1)
    columns = [(empty, nil)]
    for mask in range(1, 1 << n):
        i = (mask & -mask).bit_length() - 1
        prev = columns[mask ^ 1 << i]
        base = A[i] if i == root else A[i] + inverses[i] * (A[i] - J[i])
        columns.append((base * prev[0], base * prev[1] + ((A[i] - J[i]) * prev[0] if i == root else nil)))
    M0 = sp.Matrix.hstack(*(p[0] for p in columns))
    M1 = sp.Matrix.hstack(*(p[1] for p in columns))
    DM0 = DomainMatrix.from_Matrix(M0).to_field()
    inv = DM0.inv()
    efull = sp.eye(1 << n)[:, -1]
    mu0 = inv.matmul(DomainMatrix.from_Matrix(efull).convert_to(DM0.domain)).to_Matrix()
    mu1 = -inv.matmul(DomainMatrix.from_Matrix(M1 * mu0).convert_to(DM0.domain)).to_Matrix()
    assert M0 * mu0 == efull and M0 * mu1 + M1 * mu0 == sp.zeros(1 << n, 1)
    assert sum(mu0) == 1 and sum(mu1) == 0
    assert all((mu0[mask] == 0 for mask in range(1 << n) if not mask >> root & 1))
    assert sp.Rational(str(evaluate(normalizer, point) / common)) == M0[-1, -1]
    checked = []
    for i in outside:
        value = sp.factor(mu1[full ^ 1 << i] - mu0[full ^ 1 << i] * mu1[full] / mu0[full])
        assert value == sp.Rational(str(evaluations[i]))
        checked.append(str(value))
    print('independent full-IDLA rational source check passed', flush=True)
    result = {'scope': 'ONE saved five-site branching tree; four independent positive outside inverse rates; root=1', 'neighbors': neighbors, 'root': root, 'outside_vertices': outside, 'variables': names, 'positive_graph_denominator': common, 'recipe': {'Hcoef': 'H_U[j,j]=d_j/det L_(V minus (U union j)); H_U[i,j]=-d_i Q^(U union j)[i,j]/det L_(V minus U) for i in U,j notin U', 'bcoef': 'b_U[i]=d_i*k_i^U/det L_(V minus U), i in U', 'D0': 'det(common*H_at_root0)', 'D1': 'sum of determinants replacing one H0 column by its H1 column', 'Ni0_Ni1': 'same Cramer determinants with column i replaced by beta0,beta1', 'source': 'common*(Ni1*D0-Ni0*D1)/(D0^2*(common*full_normalizer_at_root0))', 'arithmetic': 'integer FLINT fmpz_mpoly; exact Bareiss divisions; polynomial gcd; denominator positive coefficientwise', 'hash_format': 'sha256 of sorted (exponent_tuple,decimal_coefficient) JSON, separators comma/colon'}, 'determinant_at_root0': polynomial_record(D0), 'normalizer_at_root0_scaled': str(normalizer), 'sources': records, 'all_four_source_numerators_coefficientwise_nonnegative': all((r['numerator']['negative_terms'] == 0 for r in records)), 'rational_check': {'outside_inverse_rates': point, 'root_inverse_rate': 0, 'full_basis_dimension': 1 << n, 'method': 'independent full IDLA basis and stationary derivative; no E inverse', 'normalization_and_derivative_residuals': True, 'source': checked}}
    result['root_sleep_rate'] = 'arbitrary positive finite; does not enter the source'
    result['consequence'] = 'strict left all-return-count full-event comparison and strict marginal trace bound on the specified tree/root; not full NA'
    return result


def full_event_equal_row_certificate():
    """One exact normalization check for the equal-row source theorem.

    Fix n=3, K_ij=kappa=1/4 including loops, root 0, and outside inverse
    rates (2/3,4/9). Only the root inverse t remains symbolic. Polynomial
    monomial integrals construct M u=aD independently of the full eight-state
    IDLA basis. Compare the boundary forcing with S4 and the resulting
    source with vertex_excursion_operators at this same outside point.

    The loopless complete graph is the same physical model: its inverse
    rates are 3/4 times the holding inverses. Thus lambda_hold=3/4*lambda
    and source_loopless=4/3*source_hold. All comparisons retain this scaling.
    This finite check discriminates normalizations; the all-size theorem
    uses the analytic invariant-cone and covariance proof in the report.
    """
    t, z = sp.symbols("t z", nonnegative=True)
    n, root, full = 3, 0, 7
    outside = [1, 2]
    holes, outside_holes = [6, 5, 3], [2, 1]
    r = sp.Rational(1, 4)
    kappa, d = 1-n*r, 1-(n-1)*r
    a = kappa/r
    inverse = (t, sp.Rational(2, 3), sp.Rational(4, 9))
    kernel = tuple(tuple((j, r) for j in range(n+1)) for _ in range(n))

    def clean(matrix):
        return matrix.applyfunc(sp.cancel)

    def moment(polynomial, shift=0, logarithm=False):
        """Integral z^(a-1+shift)*polynomial, optionally times -log(z)."""
        power = 2 if logarithm else 1
        return sp.expand(sum(coefficient/(a+shift+degree[0])**power
                             for degree, coefficient in sp.Poly(
                                 sp.expand(polynomial), z).terms()))

    def exact_solve(matrix, target):
        field_matrix = DomainMatrix.from_Matrix(matrix).to_field()
        field_target = DomainMatrix.from_Matrix(target).convert_to(field_matrix.domain)
        solution = field_matrix.lu_solve(field_target)
        assert (field_matrix.matmul(solution)-field_target).is_zero_matrix
        return solution.to_Matrix().applyfunc(sp.factor)

    def rows(matrix):
        return [list(map(str, row)) for row in matrix.tolist()]

    products = [sp.prod(1+inverse[j]*z for j in range(n) if j != i)
                for i in range(n)]
    A = sp.Matrix([moment((1-z)*product) for product in products])
    D = sp.Matrix([moment(product, 1, logarithm=True) for product in products])
    I = sp.Matrix([moment(product, 1) for product in products])
    B = sp.zeros(n)
    for i in range(n):
        for j in range(n):
            if i != j:
                product = sp.prod(1+inverse[k]*z for k in range(n) if k not in (i, j))
                B[i, j] = moment((1-z)*product, 1)
        assert sp.expand(I[i]-a*A[i]-sum(inverse[j]*B[i, j] for j in range(n))) == 0
    M = sp.Matrix(n, n, lambda i, j: I[i] if i == j else -inverse[j]*B[i, j])
    X = sp.diag(*inverse)
    assert clean(X*M-(X*M).T) == sp.zeros(n)
    u = exact_solve(M, a*D)

    # Independent full-state IDLA construction, with the same symbolic t.
    additions, jumps = idla_operators(kernel)
    basis = basis_matrix(additions, jumps, inverse)
    target = sp.zeros(1 << n, 1)
    target[full] = 1
    law = exact_solve(basis, target)
    assert sp.cancel(sum(law)-1) == 0
    assert clean(basis*law-target) == sp.zeros(1 << n, 1)
    block = basis.extract(holes, holes)
    injection = -basis.extract(holes, [full])
    # W=(r/d)I, C=H, and the original odds obey Cw=(d/r)s.
    assert clean(block*X-(d/r)*X*M) == sp.zeros(n)
    assert clean(injection-(d/r)*a*X*D) == sp.zeros(n, 1)
    odds = clean(law.extract(holes, [0])/law[full])
    assert clean(odds-X*u) == sp.zeros(n, 1)
    normalizer = full_occupation_normalizer(kernel, inverse)
    assert sp.cancel(law[full]*normalizer-1) == 0

    M0, u0 = M.subs(t, 0), u.subs(t, 0)
    assert all(a/(a+1) < value < 1 for value in u0)
    assert u0[root] > u0[2] > u0[1]
    differential_forcing = clean((a*D.diff(t)-M.diff(t)*u0).subs(t, 0))
    outside_matrix = M0.extract(outside, outside)
    assert M0.extract(outside, [root]) == sp.zeros(2, 1)
    assert all(value > 0 for value in outside_matrix*sp.ones(2, 1))
    assert all(value >= 0 for value in outside_matrix.inv())
    u_prime = exact_solve(outside_matrix, differential_forcing.extract(outside, [0]))
    assert clean(u.diff(t).subs(t, 0).extract(outside, [0])-u_prime) == sp.zeros(2, 1)

    # In Z_i Cov(z,h), the self term cancels. All remaining integrands
    # are polynomials or polynomials times -log(z), so no quadrature or
    # rational-log integration is used. Check the fixed-point centering too.
    theta = [sp.cancel(value.subs(t, 0)*z/(1+value.subs(t, 0)*z)) for value in inverse]
    forcing_records = []
    for local, i in enumerate(outside):
        product0 = products[i].subs(t, 0)
        polynomial = sp.cancel(product0*(-a*u0[i]
                               +sum(theta[j]*(u0[j]-u0[i]) for j in range(n))))
        centered = (moment((1-z)*polynomial)
                    +a*moment(product0, 1, logarithm=True))
        assert centered == 0
        covariance_term = (moment((1-z)*polynomial, 1)
                           +a*moment(product0, 2, logarithm=True))
        direct_root_term = B[i, root].subs(t, 0)*(u0[root]-u0[i])
        forcing = direct_root_term+covariance_term
        assert forcing == differential_forcing[i]
        assert direct_root_term > 0 and covariance_term > 0 and forcing > 0
        assert u_prime[local] > 0
        forcing_records.append({
            "site": i, "root_term": str(direct_root_term),
            "Z_times_covariance": str(covariance_term),
            "forcing": str(forcing), "u_derivative": str(u_prime[local]),
            "fixed_point_centering_checked": True,
            "S4_equals_differentiation": True, "both_terms_strictly_positive": True})
    odds_prime = clean(odds.diff(t).subs(t, 0))
    assert odds_prime[root] == u0[root]
    for local, i in enumerate(outside):
        assert odds_prime[i] == inverse[i]*u_prime[local]

    traced = trace_kernel(kernel, [root])
    expected_trace = tuple(tuple((j, sp.Rational(1, 3)) for j in range(n))
                           for _ in outside)
    assert traced == expected_trace
    assert (1-(n-1)*sp.Rational(1, 3))/sp.Rational(1, 3) == a
    nu = sp.Matrix(stationary(traced, [inverse[i] for i in outside]))
    occupied = [1+sum(((mask >> j) & 1) << i for j, i in enumerate(outside))
                for mask in range(4)]
    empty = [mask-1 for mask in occupied]
    assert law.subs(t, 0).extract(occupied, [0]) == nu
    assert law.subs(t, 0).extract(empty, [0]) == sp.zeros(4, 1)
    source_from_odds = nu[3]*odds_prime.extract(outside, [0])

    def excursion_source(operators):
        E = sp.eye(4)-operators["R"]
        trace_law = operators["nu"]
        source = ((E.extract(outside_holes, outside_holes)-E[3, 3]*sp.eye(2))
                  *trace_law.extract(outside_holes, [0])
                  +E.extract(outside_holes, [3])*trace_law[3])
        return E, source

    # This API takes rational rates. Root t=1 extracts E once; the next
    # residual checks its occupied-root resolvent for the entire symbolic t.
    excursions = vertex_excursion_operators(kernel, (sp.S.One, *inverse[1:]), root)
    E, source = excursion_source(excursions)
    assert excursions["nu"] == nu and source == source_from_odds
    assert all(value > 0 for value in source)
    occupied_law = law.extract(occupied, [0])
    assert clean((sp.eye(4)+t*excursions["R"])*occupied_law-nu) == sp.zeros(4, 1)
    assert clean(law.extract(empty, [0])-t*excursions["J"]*occupied_law) == sp.zeros(4, 1)
    assert sp.cancel(E[3, 3]-1+normalizer.diff(t)/normalizer.subs(t, 0)) == 0

    # Same physical point after erasing the inserted self-loops. The basis
    # identity checks all t before the independent rational excursion check.
    simple_kernel = graph_kernel(((1, 2, 3), (0, 2, 3), (0, 1, 3)))
    scale = 1-r
    simple_inverse = tuple(scale*value for value in inverse)
    simple_additions, simple_jumps = idla_operators(simple_kernel)
    assert simple_additions == additions
    for i in range(n):
        assert jumps[i] == r*additions[i]+scale*simple_jumps[i]
    simple_basis = basis_matrix(simple_additions, simple_jumps, simple_inverse)
    assert simple_basis == basis
    simple_excursions = vertex_excursion_operators(
        simple_kernel, (scale, *simple_inverse[1:]), root)
    simple_E, simple_source = excursion_source(simple_excursions)
    assert simple_excursions["nu"] == nu
    assert simple_excursions["R"] == excursions["R"]/scale
    assert simple_excursions["J"] == excursions["J"]/scale
    assert simple_excursions["H"] == excursions["H"]
    assert simple_source == source/scale
    assert simple_E[3, 3] == 1-(1-E[3, 3])/scale
    original_t = sp.Symbol("original_root_inverse", nonnegative=True)
    simple_odds_prime = clean(odds.subs(t, original_t/scale).diff(original_t).subs(original_t, 0))
    assert simple_source == nu[3]*simple_odds_prime.extract(outside, [0])

    print("Equal-row integral/IDLA, positive S4 forcing, excursion source, and holding-scale checks passed.",
          flush=True)
    return {
        "scope": "one n=3 equal-row point with symbolic root inverse; no graph or outside-rate sweep",
        "root": root, "outside_vertices": outside,
        "holding_kernel_probability": str(r), "killing": str(kappa), "a": str(a),
        "holding_inverse_rates": list(map(str, inverse)),
        "method": "exact polynomial monomial and log-monomial integrals; independent full eight-state IDLA basis and excursion operators",
        "integral_matrix_M": rows(M), "integral_D": list(map(str, D)),
        "normalized_odds_u": list(map(str, u)),
        "full_IDLA_stationary_law": list(map(str, law)),
        "full_IDLA_residual_and_normalization_checked": True,
        "original_odds_and_scaled_block_checked": True,
        "scaling": {"W": str(r/d), "d_over_r": str(d/r),
                    "identity": "Cw=(d/r)s; C*diag(x)=(d/r)*diag(x)*M"},
        "boundary": {
            "u": list(map(str, u0)), "outside_matrix": rows(outside_matrix),
            "forcing": forcing_records, "odds_derivative": list(map(str, odds_prime)),
            "trace_equal_row_probability": "1/3", "trace_a": str(a),
            "trace_stationary_law": list(map(str, nu)), "nu_full": str(nu[3]),
            "source": list(map(str, source)), "source_strictly_positive": True,
            "source_matches_independent_excursion_operators": True,
            "symbolic_root_resolvent_and_full_normalizer_checked": True,
            "excursion_full_diagonal": str(E[3, 3])},
        "holding_conversion": {
            "sleep_rate_factor_hold_over_loopless": str(scale),
            "loopless_inverse_rates": list(map(str, simple_inverse)),
            "full_symbolic_IDLA_basis_equal": True,
            "trace_stationary_laws_equal": True,
            "source_factor_loopless_over_hold": str(1/scale),
            "loopless_source": list(map(str, simple_source)),
            "loopless_source_matches_odds_derivative": True,
            "return_and_exit_operators_have_required_scaling": True},
        "analytic_inequality": {
            "statement": "v(z)-z*v'(z)=integral_0^1 (1-s)*(z/(s+(1-s)*z))^2 ds < 1/2 for 0<z<1",
            "status": "proved analytic identity and strict bound in the theorem; not a finite test"},
        "theorem_scope": "equal-row source positivity for a>=1 uses the analytic proof; this single certificate checks its normalizations, not full NA"}


def equal_row_vacancies_certificate():
    """Reproduce the saved pair-sector check and one saved top-ULC point.

    The asymmetric four-site point exposes positive disjoint-pair entries
    and checks 66 integral coefficients against the full IDLA basis. The
    homogeneous check reads the existing complete-graph n=4, lambda=1 rank
    law; it computes only scalar integrals at holding inverse rate 5/4.
    No homogeneous rank law is recomputed and no parameter sweep is used.

    The all-sector identity and all-n homogeneous top-ULC theorem are proved
    analytically in the report. These two fixed checks verify normalization;
    they do not establish full ULC, NA, or unequal-rate conditional NC.
    """
    import json
    from itertools import combinations
    from pathlib import Path

    n, a = 4, sp.S.One
    x = tuple(map(sp.Integer, (1, 2, 3, 4)))
    r = 1/(a+n)
    kernel = tuple(tuple((j, r) for j in range(n))+((n, a*r),) for _ in range(n))
    additions, jumps = idla_operators(kernel)
    M = basis_matrix(additions, jumps, x)
    full = (1 << n)-1
    pairs = list(combinations(range(n), 2))
    z = sp.Symbol("z")

    def mask(sites):
        return sum(1 << i for i in sites)

    def moment(polynomial, shift=0, log_power=0):
        """Integral z^(a-1+shift)*polynomial*(-log(z))^log_power."""
        return sp.factor(sum(coefficient*sp.factorial(log_power)
                             /(a+shift+degree[0])**(log_power+1)
                             for degree, coefficient in sp.Poly(
                                 sp.cancel(polynomial), z).terms()))

    def integral_entry(I, J):
        I, J = set(I), set(J)
        s, m = len(I), len(I)-len(J)
        if m < 0:
            return sp.S.Zero
        delta = len(J-I)
        polynomial = ((1-z)**delta
                      *sp.prod(1+x[k]*z for k in range(n) if k not in I | J))
        return sp.factor((-1)**len(I-J)*(a+s)*a**m
                         *sp.prod(x[i] for i in I-J)*moment(polynomial, s, m)
                         /sp.factorial(m))

    columns = [()]+[(i,) for i in range(n)]+pairs
    for I in pairs:
        for J in columns:
            assert integral_entry(I, J) == M[full ^ mask(I), full ^ mask(J)]
    N = M.inv()
    law = N[:, full]
    assert min(law) > 0 and sum(law) == 1
    assert M*law == sp.eye(1 << n)[:, full]
    w = sp.Matrix([law[full ^ (1 << i)]/law[full] for i in range(n)])
    u = sp.Matrix([w[i]/x[i] for i in range(n)])
    q = sp.Matrix([law[full ^ mask(I)]/law[full] for I in pairs])
    X = sp.diag(*(sp.prod(x[i] for i in I) for I in pairs))
    row_indices = [full ^ mask(I) for I in pairs]
    H = M.extract(row_indices, row_indices)
    N2 = N.extract(row_indices, row_indices)
    assert H*N2 == sp.eye(len(pairs))
    assert min(N2) >= 0 and all(sum(N2[:, j]) <= 1 for j in range(len(pairs)))
    L2 = X.inv()*H*X/(a+2)
    inverse_L2 = (a+2)*X.inv()*N2*X
    assert L2*inverse_L2 == sp.eye(len(pairs)) and min(inverse_L2) >= 0
    forcing = -X.inv()*(M.extract(row_indices, [full ^ (1 << i) for i in range(n)])*w
                       +M.extract(row_indices, [full]))/(a+2)
    v = X.inv()*q
    assert L2*v == forcing
    product = sp.Matrix([u[i]*u[j] for i, j in pairs])
    residual = L2*product-forcing
    gaps = product-v
    disjoint_row, disjoint_col = pairs.index((0, 1)), pairs.index((2, 3))
    assert M[full ^ mask((0, 1)), full ^ mask((2, 3))] == sp.Rational(1, 5)
    assert L2[disjoint_row, disjoint_col] == sp.Rational(2, 5) > 0
    assert min(gaps) > 0 and min(residual) > 0
    # Preserve all values from the existing single-point scratch certificate.
    pair_record = dict(
        scope='One exact example, chosen to expose the disjoint-pair sign at the smallest possible size.',
        n=n, a=str(a), kernel_entry=str(r), inverse_rates=list(map(str, x)),
        integral_entries_checked=len(pairs)*len(columns),
        independent_IDLA_basis_agrees=True,
        normalized_two_vacancy_system_verified=True,
        pair_block_inverse_nonnegative=True,
        pair_order=[list(I) for I in pairs],
        one_hole_u=list(map(str, u)), two_hole_v=list(map(str, v)),
        normalized_conditional_NC_gaps=list(map(str, gaps)),
        product_supersolution_residual=list(map(str, residual)),
        normalized_forcing=list(map(str, forcing)),
        disjoint_pair_entry=dict(target_holes=[0, 1], source_holes=[2, 3],
                                 raw_M='1/5', normalized_L2='2/5'),
        all_conditional_NC_gaps_positive=bool(min(gaps) > 0),
        all_product_residuals_nonnegative=bool(min(residual) >= 0))

    # Reuse the already saved rank law. Its entries are rank probabilities,
    # so divide by n and binomial(n,2) before comparing per-configuration odds.
    certificate_path = Path(__file__).with_name("complete_graph_certificates.json")
    saved = json.loads(certificate_path.read_text())
    matches = [case for case in saved if case.get("n") == n
               and sp.Rational(case.get("lambda", "0")) == 1]
    assert len(matches) == 1
    saved_case = matches[0]
    assert saved_case["independent_refresh_check"]
    ranks = tuple(map(sp.Rational, saved_case["rank"]))
    assert len(ranks) == n+1 and sum(ranks) == 1 and min(ranks) > 0
    original_sleep = sp.Rational(saved_case["lambda"])
    holding_sleep = sp.Rational(n, n+1)*original_sleep
    common_inverse = 1/holding_sleep
    assert common_inverse == sp.Rational(5, 4)
    m = n-2
    P = (1+common_inverse*z)**m
    Y = (1-z)**2*P
    transfer = common_inverse*(1-z)/(1+common_inverse*z)
    assert sp.cancel(Y.diff(z, 2)/2-P*(1-2*m*transfer+sp.binomial(m, 2)*transfer**2)) == 0
    assert sp.cancel(-Y.diff(z)/(1-z)-P*(2-m*transfer)) == 0
    eta_mass = moment(Y)
    eta1_mass = moment((1-z)*P)
    one_product = (1+common_inverse*z)*P
    one_mass = moment((1-z)*one_product)
    homogeneous_u = moment(one_product, 1, 1)/one_mass
    L = moment(P*(1-2*m*transfer+sp.binomial(m, 2)*transfer**2), 2)
    assert L == a*(a+1)*eta_mass/2
    scalar_forcing = (a*homogeneous_u*moment(P*(2-m*transfer), 2, 1)
                      -a**2*moment(P, 2, 2)/2)
    homogeneous_b = scalar_forcing/L
    saved_u = ranks[n-1]/(n*common_inverse*ranks[n])
    saved_b = ranks[n-2]/(sp.binomial(n, 2)*common_inverse**2*ranks[n])
    assert homogeneous_u == saved_u and homogeneous_b == saved_b

    # All voltage moments reduce to polynomial log-moments after multiplying
    # by the eta density; the factors (1-z)^2 cancel the voltage denominators.
    mean_voltage = moment((1-z)*P, 1, 1)/eta_mass
    mean_voltage_squared = moment(P, 2, 2)/eta_mass
    mean_z_voltage_prime = (moment(P, 1, 1)-moment((1-z)*P, 1))/eta_mass
    assert homogeneous_b == (2*homogeneous_u*(a*mean_voltage+mean_z_voltage_prime)
                              -a*mean_voltage_squared)/(a+1)
    variance_voltage = mean_voltage_squared-mean_voltage**2
    centered_square = mean_voltage_squared-2*homogeneous_u*mean_voltage+homogeneous_u**2
    assert centered_square == variance_voltage+(mean_voltage-homogeneous_u)**2
    bracket = homogeneous_u-2*mean_z_voltage_prime
    ulc_gap = homogeneous_u**2-homogeneous_b
    assert (a+1)*ulc_gap == a*centered_square+homogeneous_u*bracket
    assert variance_voltage > 0 and centered_square > 0 and bracket > 0 and ulc_gap > 0

    # The finitely many components are the monomials of this fixed P, not
    # additional graph/rate cases. Check both beta identities componentwise.
    beta_records = []
    beta2_total = beta3_total = beta2_mean = beta3_mean = sp.S.Zero
    beta3_relative_times_value = sp.S.Zero
    for (degree,), coefficient in sorted(sp.Poly(P, z).terms()):
        shape = a+degree
        mass2 = 1/(shape*(shape+1))
        mass3 = 2/(shape*(shape+1)*(shape+2))
        mean2 = 1/(shape+1)**2/mass2
        twice_mean3 = 2*(1/(shape+1)**2-1/(shape+1)+1/(shape+2))/mass3
        assert mean2 == twice_mean3 == shape/(shape+1)
        relative = mass2/mass3
        assert relative == (shape+2)/2
        beta2_total += coefficient*mass2
        beta3_total += coefficient*mass3
        beta2_mean += coefficient*mass2*mean2
        beta3_mean += coefficient*mass3*twice_mean3
        beta3_relative_times_value += coefficient*mass3*relative*twice_mean3
        beta_records.append({"degree": degree, "coefficient": str(coefficient),
                             "shape_A": str(shape), "beta_A_2": str(mass2),
                             "beta_A_3": str(mass3), "common_mean": str(mean2),
                             "relative_weight": str(relative)})
    assert beta2_total == eta1_mass and beta3_total == eta_mass
    mean_eta1_voltage = moment(P, 1, 1)/eta1_mass
    assert beta2_mean/beta2_total == mean_eta1_voltage
    assert beta3_mean/beta3_total == 2*mean_z_voltage_prime
    beta_covariance = (beta3_relative_times_value/beta3_total
                       -(beta2_total/beta3_total)*(beta3_mean/beta3_total))
    beta_gain = mean_eta1_voltage-2*mean_z_voltage_prime
    assert beta_gain == beta_covariance/(beta2_total/beta3_total) > 0
    mean_eta1_z = moment((1-z)*P, 1)/eta1_mass
    mean_eta1_z_voltage = moment(P, 2, 1)/eta1_mass
    tilt_covariance = mean_eta1_z_voltage-mean_eta1_z*mean_eta1_voltage
    tilt_gain = homogeneous_u-mean_eta1_voltage
    assert tilt_gain == common_inverse*tilt_covariance/(1+common_inverse*mean_eta1_z) > 0
    assert bracket == beta_gain+tilt_gain
    rank_gap = ranks[n-1]**2/n**2-ranks[n-2]*ranks[n]/sp.binomial(n, 2)
    assert rank_gap == ranks[n]**2*common_inverse**2*ulc_gap > 0

    print("Equal-row vacancies: all 66 entries, pair inverse/signs, saved homogeneous law, and top-ULC integral identities passed.",
          flush=True)
    return {
        "scope": "the existing asymmetric four-site point and the existing homogeneous n=4 lambda=1 rank certificate only",
        "two_vacancy": pair_record,
        "same_rank_block_stabilization_and_substochasticity_checked": True,
        "homogeneous_top_ulc": {
            "source_file": certificate_path.name,
            "saved_case": {"n": n, "lambda": str(original_sleep)},
            "saved_rank_law": list(map(str, ranks)),
            "new_rank_law_computed": False,
            "a": str(a), "holding_sleep_rate": str(holding_sleep),
            "holding_inverse_rate": str(common_inverse),
            "u": str(homogeneous_u), "b": str(homogeneous_b),
            "scalar_L": str(L), "scalar_forcing": str(scalar_forcing),
            "saved_top_rank_normalizations_match": True,
            "eta_mass": str(eta_mass), "eta1_mass": str(eta1_mass),
            "mean_voltage": str(mean_voltage),
            "mean_voltage_squared": str(mean_voltage_squared),
            "mean_z_voltage_prime": str(mean_z_voltage_prime),
            "variance_voltage": str(variance_voltage),
            "centered_square": str(centered_square),
            "positive_T2_terms": [str(a*centered_square), str(homogeneous_u*bracket)],
            "beta_components": beta_records,
            "beta_mixture_covariance": str(beta_covariance),
            "beta_mixture_gain": str(beta_gain),
            "tilt_covariance": str(tilt_covariance), "tilt_gain": str(tilt_gain),
            "bracket": str(bracket), "u_squared_minus_b": str(ulc_gap),
            "normalized_rank_gap": str(rank_gap),
            "integration_by_parts_beta_and_variance_identities_checked": True},
        "analytic_theorem": {
            "status": "proved analytically in the report; separate from these finite checks",
            "sector_formula_scope": "all n, a>0, and positive site inverse rates",
            "top_ULC_scope": "n>=2, a>0, common positive inverse rate x; strict top normalized-rank inequality",
            "limitations": "full ULC/NA and unequal-rate conditional NC are not inferred"}}


def third_sector_lc_certificate():
    """Strict third-sector LC for arbitrary two-monomial support.

    eta_s is proportional to z^(a-1)(1-z)^s(1+x*z)^(3-s)*P(z).
    Certify b_2^2>b_1*b_3 for P=z^k, a>0,k>=0,x>=0, and for
    P=c0+c2*z^2, a>0,x>=0,c0,c2>=0,c0+c2>0. Retain both certificates
    and add arbitrary support {i,i+k}, i,k>=0. No assertion for mixtures
    with three or more distinct powers or all sizes of the physical
    binomial family is inferred. No parameter sweep.
    """
    from hashlib import sha256
    from math import comb, factorial, prod
    import json
    try:
        from flint import fmpz_mpoly_ctx
    except ImportError as error:
        raise SystemExit("--third-sector-lc requires optional python-flint") from error

    names = ("Z1", "N1", "Z2", "M2", "H2", "Z3", "H3", "J3", "K3")

    def linear_forms(A, a, x, support, reciprocal):
        """Nine linear forms, optionally multiplied by a common denominator.

        The measure has base exponent A-1; a remains the operator parameter.
        At support index k, the adjoint of D+a is -D-(A+k-a).
        reciprocal(d) is either 1/d or the exact polynomial quotient L/d.
        """
        def Z(s, k):
            return factorial(s)*reciprocal(prod(A+k+j for j in range(s+1)))

        def N1(k):
            return reciprocal((A+k+1)**2)

        def M2(k):
            h = A+k-a
            return (h+2)*reciprocal((A+k+2)**2)-h*reciprocal((A+k+1)**2)

        def H2(k):
            return 2*reciprocal((A+k+2)**3)

        def H3(k):
            return 6*reciprocal((A+k+3)**4)

        def J3(k):
            h = A+k-a
            return 2*((h+3)*reciprocal((A+k+3)**3)-h*reciprocal((A+k+2)**3))

        def K3(k):
            h = A+k-a
            return ((h+2)*(h+3)*reciprocal((A+k+3)**2)
                    -2*h*(h+2)*reciprocal((A+k+2)**2)
                    +h*(h-1)*reciprocal((A+k+1)**2))

        def T(r, kernel):
            return sum(weight*comb(r, ell)*x**ell*kernel(k+ell)
                       for k, weight in support for ell in range(r+1))

        return {"Z1": T(2, lambda k: Z(1, k)), "N1": T(2, N1),
                "Z2": T(1, lambda k: Z(2, k)), "M2": T(1, M2), "H2": T(1, H2),
                "Z3": T(0, lambda k: Z(3, k)), "H3": T(0, H3),
                "J3": T(0, J3), "K3": T(0, K3)}

    def values(forms, a):
        Z1, N1, Z2, M2, H2, Z3, H3, J3, K3 = (forms[name] for name in names)
        u = sp.cancel(N1/Z1)
        b = sp.cancel((2*u*M2-a*H2)/((a+1)*Z2))
        c = sp.cancel((a**2*H3-3*a*u*J3+3*b*K3)/((a+1)*(a+2)*Z3))
        return u, b, c

    def quintic(forms, a):
        Z1, N1, Z2, M2, H2, Z3, H3, J3, K3 = (forms[name] for name in names)
        B = 2*N1*M2-a*Z1*H2
        numerator = ((a+2)*Z3*B**2-a**2*(a+1)*N1*Z1*Z2**2*H3
                     +3*a*(a+1)*N1**2*Z2**2*J3-3*N1*B*Z2*K3)
        denominator = (a+1)**2*(a+2)*Z1**2*Z2**2*Z3
        return numerator, denominator

    # Independent integral engine: differentiate first, cancel the whole
    # integrand, then integrate polynomial powers of z and ell=-log(z).
    z, ell = sp.symbols("z ell")
    v = z*ell/(1-z)

    def D(expression):
        return sp.cancel(z*sp.diff(expression, z)-sp.diff(expression, ell))

    def integral(expression, A):
        polynomial = sp.Poly(sp.cancel(expression), z, ell)
        return sp.factor(sum(coefficient*factorial(log_power)/(A+power)**(log_power+1)
                             for (power, log_power), coefficient in polynomial.terms()))

    def direct_forms(A, a, x, P):
        expressions = {
            "Z1": (1-z)*(1+x*z)**2,
            "N1": (1-z)*(1+x*z)**2*v,
            "Z2": (1-z)**2*(1+x*z),
            "M2": (1-z)**2*(1+x*z)*(a*v+D(v)),
            "H2": (1-z)**2*(1+x*z)*v**2,
            "Z3": (1-z)**3,
            "H3": (1-z)**3*v**3,
            "J3": (1-z)**3*(a*v**2+D(v**2)),
            "K3": (1-z)**3*(a*(a+1)*v+(2*a+1)*D(v)+D(D(v))),
        }
        return {name: integral(P*expr, A) for name, expr in expressions.items()}

    # The cleared sign problem is homogeneous quintic in the P coefficients.
    a = sp.Symbol("a", positive=True)
    formal = dict(zip(names, sp.symbols(" ".join(names))))
    u, b, c = values(formal, a)
    formal_num, formal_den = quintic(formal, a)
    assert sp.cancel(b**2-u*c-formal_num/formal_den) == 0
    assert sp.Poly(formal_num, *formal.values()).homogeneous_order() == 5

    # One component. Keep A independent until the positivity shift A=a+k.
    A, x = sp.symbols("A x", positive=True)
    k = sp.Symbol("k", nonnegative=True)
    component = linear_forms(A, a, x, [(0, sp.S.One)], lambda d: 1/d)
    direct_component = direct_forms(A, a, x, sp.S.One)
    assert all(sp.cancel(component[name]-direct_component[name]) == 0 for name in names)
    u, b, c = values(component, a)
    assert all(sp.cancel(left-right) == 0
               for left, right in zip((u, b, c), values(direct_component, a)))
    gap = sp.factor(b**2-u*c)
    numerator, denominator = sp.fraction(gap)
    raw = sp.Poly(numerator, A, a, x)
    num = sp.Poly(sp.expand(numerator.subs(A, a+k)), a, k, x)
    den = sp.Poly(sp.expand(denominator.subs(A, a+k)), a, k, x)

    def component_digest(polynomial):
        data = ";".join(",".join(map(str, powers))+":"+str(coefficient)
                        for powers, coefficient in polynomial.terms())
        return sha256(data.encode("ascii")).hexdigest()

    assert all(coefficient > 0 for polynomial in (num, den) for coefficient in polynomial.coeffs())
    assert len(num.terms()) == 1074 and len(den.terms()) == 1404
    numerator_hash, denominator_hash = component_digest(num), component_digest(den)
    assert numerator_hash == "ecc717f937c3eb6cbc2f3bb345a1fe1a8dfb874ff618db0073280fecac80ec9f"
    assert denominator_hash == "86a9d29ba52a35a3013f70459ec02abaca4ca5e6b9fa27c4a41b98afd611ec94"
    base = a**4*(a+1)**4*(a+2)**3*(a+3)**3*(a**3+13*a**2+45*a+47)
    assert sp.expand(num.as_expr().subs({k: 0, x: 0})-base) == 0

    # Retain all eleven independently audited component moments.
    Z3 = 6/(A*(A+1)*(A+2)*(A+3))
    B2 = lambda y: 1/(y*(y+1))
    B3 = lambda y: 2/(y*(y+1)*(y+2))
    F2 = lambda y: a*(1/(y+1)**2-1/(y+2)**2)+1/((y+1)**2*(y+2))
    Mv = A*(3*A**2+12*A+11)/(3*(A+1)*(A+2)*(A+3))
    Mw = A*(3*A+5)/(6*(A+1)*(A+2))
    Mkappa = 2*A/(3*(A+1))
    expected_moments = {
        "Z1": sum(comb(2, j)*x**j*B2(A+j) for j in range(3)),
        "N1": sum(comb(2, j)*x**j/(A+j+1)**2 for j in range(3)),
        "Z2": B3(A)+x*B3(A+1),
        "R2": sum(x**j*2/(A+j+2)**3 for j in range(2)),
        "F2": F2(A), "Mv": Mv, "Mw": Mw, "Mkappa": Mkappa,
        "Mv3": A*(A+1)*(A+2)/(A+3)**3,
        "MF": A*(A+1)*(A+2)*(A+3)/3*((A-a+3)/(A+3)**3-(A-a)/(A+2)**3),
        "MG": a*(a+1)*Mv+2*a*Mw+Mkappa,
    }
    derived_moments = {name: direct_component[name] for name in ("Z1", "N1", "Z2")}
    derived_moments.update({
        "R2": direct_component["H2"], "F2": integral((1-z)**2*(a*v+D(v)), A),
        "Mv": integral((1-z)**3*v, A)/Z3,
        "Mw": integral((1-z)**3*D(v), A)/Z3,
        "Mkappa": integral((1-z)**3*(D(D(v))+D(v)), A)/Z3,
        "Mv3": direct_component["H3"]/Z3,
        "MF": direct_component["J3"]/Z3, "MG": direct_component["K3"]/Z3,
    })
    assert all(sp.cancel(derived_moments[name]-expr) == 0 for name, expr in expected_moments.items())
    component_record = {
        "scope": "P=z^k, a>0,k>=0,x>=0; no mixture closure inferred",
        "variables": ["a", "k", "x"], "shift": "A=a+k",
        "raw_A_a_x_numerator_terms": len(raw.terms()),
        "raw_A_a_x_negative_terms": sum(1 for coefficient in raw.coeffs() if coefficient < 0),
        "shifted_a_k_x_numerator_terms": len(num.terms()),
        "shifted_a_k_x_negative_terms": 0,
        "denominator_terms": len(den.terms()), "denominator_negative_terms": 0,
        "denominator_factored": str(sp.factor(denominator)),
        "numerator_hash": numerator_hash, "denominator_hash": denominator_hash,
        "hash_format": "SymPy Poly.terms descending lex order; semicolon-separated exponent_csv:integer_coefficient",
        "numerator_at_k0_x0": str(base),
        "gap_at_x0": str(sp.factor(gap.subs(x, 0))),
        "gap_limit_x_infinity": str(sp.factor(sp.Poly(numerator, x).LC()/sp.Poly(denominator, x).LC())),
        "direct_integral_audit": {
            "moment_identities": list(expected_moments), "all_eleven_moments_match": True,
            "all_nine_linear_forms_match": True, "all_three_sector_solutions_match": True},
    }

    # Two support points. All variables remain formal in integer arithmetic.
    ctx = fmpz_mpoly_ctx.get(["a", "x", "w"])
    aa, xx, ww = ctx.gens()
    L = aa*(aa+1)**2*(aa+2)**3*(aa+3)**4*(aa+4)**3*(aa+5)**4

    def cleared_reciprocal(d):
        quotient = L/d
        assert quotient*d == L
        return quotient

    cleared = linear_forms(aa, aa, xx, [(0, ctx.constant(1)), (2, ww)], cleared_reciprocal)
    cleared_num, cleared_den = quintic(cleared, aa)
    data = cleared_num.to_dict()
    assert all(coefficient > 0 for coefficient in data.values())
    rows = sorted([*map(int, powers), int(coefficient)] for powers, coefficient in data.items())
    digest = sha256(json.dumps(rows, separators=(",", ":")).encode("ascii")).hexdigest()
    assert len(rows) == 2890 and list(map(int, cleared_num.degrees())) == [72, 6, 5]
    assert digest == "8634f8671e0d5afaf41bd5077741c76e186341ab012ed2a535a93d8d5ba5bbdb"
    endpoint_terms = {(4, 0, 0): 2545769580488294400000000000000000000,
                      (6, 0, 5): 1671768834048000000000000000000}
    assert all(data[powers] == coefficient for powers, coefficient in endpoint_terms.items())

    def evaluate(polynomial, point):
        return sum((sp.Integer(coefficient)*prod(value**power for value, power in zip(point, powers))
                    for powers, coefficient in polynomial.to_dict().items()), sp.S.Zero)

    direct = direct_forms(sp.S.One, sp.S.One, sp.S.One, 1+z**2)
    expected = dict(zip(names, map(sp.Rational, [
        "17/15", "141/200", "7/15", "141/400", "1567/12000", "4/15",
        "97/3456", "11669/108000", "187/450"])))
    assert direct == expected
    point = [sp.S.One]*3
    assert all(evaluate(cleared[name], point)/evaluate(L, point) == direct[name] for name in names)
    u, b, c = values(direct, sp.S.One)
    direct_gap = sp.cancel(b**2-u*c)
    assert (u, b, c, direct_gap) == tuple(map(sp.Rational, [
        "423/680", "125651/380800", "15280457/102816000", "7146358579/435025920000"]))
    assert evaluate(cleared_num, point)/evaluate(cleared_den, point) == direct_gap

    w = sp.Symbol("w", nonnegative=True)
    initial = values(direct_forms(sp.S.One, sp.S.One, sp.S.Zero, 1+w*z**2), sp.S.One)
    initial_gap = sp.cancel(initial[1]**2-initial[0]*initial[2])
    expected_initial = (
        3*(w+4)/(4*(w+6)),
        (1971*w**2+20156*w+48000)/(3600*(w+6)*(w+10)),
        (84403*w**3+1584093*w**2+9361350*w+17212500)/(216000*(w+6)*(w+10)*(w+15)),
        (43353*w**5+6634746*w**4+108069938*w**3+671800320*w**2
         +1817167500*w+1788750000)/(6480000*(w+6)**2*(w+10)**2*(w+15)),
    )
    assert all(sp.cancel(left-right) == 0 for left, right in zip((*initial, initial_gap), expected_initial))
    point = [sp.S.One, sp.S.Zero, w]
    assert sp.cancel(evaluate(cleared_num, point)/evaluate(cleared_den, point)-initial_gap) == 0
    initial_num, initial_den = sp.fraction(initial_gap)
    assert sp.degree(initial_num, w) == sp.degree(initial_den, w) == 5
    assert sp.gcd(initial_num, initial_den) == 1
    two_support_record = {
        "scope": "P=c0+c2*z^2; a>0,x>=0,c0,c2>=0,c0+c2>0",
        "variables": ["a", "x", "w"], "substitution": "P=1+w*z^2; homogenize to c0+c2*z^2",
        "clearing_factor_L": "a*(a+1)^2*(a+2)^3*(a+3)^4*(a+4)^3*(a+5)^4",
        "polynomial": "L^5 times the quintic numerator; no factor or content removed",
        "degrees": [72, 6, 5], "positive_terms": len(rows), "negative_terms": 0,
        "minimum_coefficient": min(row[-1] for row in rows), "sha256": digest,
        "hash_format": "Sorted rows [a_power,x_power,w_power,integer_coefficient]; compact ASCII JSON without newline",
        "terms_by_w_degree": [sum(row[2] == degree for row in rows) for degree in range(6)],
        "strict_endpoint_terms": [{"powers": list(powers), "coefficient": coefficient}
                                  for powers, coefficient in endpoint_terms.items()],
        "strictness": "After homogenization the numerator is >= C0*a^4*c0^5+C2*a^6*c2^5>0",
        "direct_integral_audit": {
            "a": 1, "x": 1, "P": "1+z^2", "all_nine_forms_match": True,
            "forms": {name: str(value) for name, value in direct.items()},
            "u": str(u), "b": str(b), "c": str(c), "gap": str(direct_gap),
            "cleared_quintic_matches": True},
        "initial_discriminator": {
            "a": 1, "x": 0, "P": "1+w*z^2", "symbolic_identity_in_w": True,
            "u_b_c_gap": list(map(str, expected_initial)),
            "numerator_degree": 5, "denominator_degree": 5, "gcd": 1},
    }

    # Arbitrary support {i,i+k}. Clear each component separately before
    # mixing: this uses degree-10 D(A), not a product of both denominators.
    # For physical coefficients ci>0,cj>=0, w=cj*Di/(ci*Dj) is nonnegative,
    # and the mixed cleared forms equal (Di/ci) times the physical forms.
    generic_ctx = fmpz_mpoly_ctx.get(["a", "i", "k", "x", "w"])
    ga, gi, gk, gx, gw = generic_ctx.gens()

    def cleared_component(index):
        origin = ga+index
        scale = origin*(origin+1)**2*(origin+2)**3*(origin+3)**4

        def reciprocal(d):
            quotient = scale/d
            assert quotient*d == scale
            return quotient

        return scale, linear_forms(origin, ga, gx, [(0, generic_ctx.constant(1))], reciprocal)

    Di, lower = cleared_component(gi)
    Dj, upper = cleared_component(gi+gk)
    mixed = {name: lower[name]+gw*upper[name] for name in names}
    generic_num, generic_den = quintic(mixed, ga)
    generic_data = generic_num.to_dict()
    assert all(coefficient > 0 for coefficient in generic_data.values())
    generic_rows = sorted([*map(int, powers), int(coefficient)]
                          for powers, coefficient in generic_data.items())
    generic_hash = sha256(json.dumps(generic_rows, separators=(",", ":")).encode("ascii")).hexdigest()
    generic_degrees = list(map(int, generic_num.degrees()))
    assert len(generic_rows) == 296701 and generic_degrees == [37, 34, 34, 6, 5]
    assert generic_hash == "3515e1db73a3cebd6680d7535d84a33653da2eb7823bcc0ae198beaeadab0fe9"
    by_weight = [sum(powers[-1] == degree for powers in generic_data) for degree in range(6)]
    assert by_weight == [5044, 36965, 55850, 64350, 67632, 66860]

    # The w^0 and w^5 coefficients are exactly the two pure quintics.
    for degree, forms in ((0, lower), (5, upper)):
        endpoint = generic_ctx.from_dict({(*powers[:-1], 0): coefficient
                                          for powers, coefficient in generic_data.items()
                                          if powers[-1] == degree})
        assert endpoint == quintic(forms, ga)[0]

    # Symbolic rescaling identity, independent of the coefficient arithmetic.
    ci, cj, di, dj, f0, f1 = sp.symbols("ci cj di dj f0 f1", positive=True)
    weight_change = cj*di/(ci*dj)
    assert sp.cancel((ci/di)*(di*f0+weight_change*dj*f1)-(ci*f0+cj*f1)) == 0

    # Reuse the retained direct-integral point a=x=ci=cj=1, i=0,k=2.
    # The transformed weight is 4/1125; the component scales are positive.
    rational_point = [sp.S.One, sp.S.Zero, sp.Integer(2), sp.S.One, sp.Rational(4, 1125)]
    di_value, dj_value = evaluate(Di, rational_point), evaluate(Dj, rational_point)
    assert (di_value, dj_value) == (27648, 7776000)
    assert rational_point[-1] == di_value/dj_value
    assert all(evaluate(mixed[name], rational_point)/di_value == direct[name] for name in names)

    def known_point(polynomial):
        # Integer evaluation avoids hundreds of thousands of rational additions.
        # a=x=1, i=0,k=2, w=4/1125; the weight degree is at most five.
        powers_two = [2**j for j in range(int(polynomial.degrees()[2])+1)]
        weight_scales = [4**j*1125**(5-j) for j in range(6)]
        total = sum(int(coefficient)*powers_two[powers[2]]*weight_scales[powers[4]]
                    for powers, coefficient in polynomial.to_dict().items() if powers[1] == 0)
        return sp.Rational(total, 1125**5)

    generic_gap = known_point(generic_num)/known_point(generic_den)
    assert generic_gap == direct_gap
    generic_record = {
        "scope": "P=ci*z^i+cj*z^(i+k), a>0,x>=0,i,k>=0,ci,cj>=0,ci+cj>0; strict b^2>u*c",
        "variables": ["a", "i", "k", "x", "w"],
        "component_clearing": "D(A)=A*(A+1)^2*(A+2)^3*(A+3)^4; Mi=D(a+i)*fi, Mj=D(a+i+k)*fj",
        "weight_change": "w=cj*D(a+i)/(ci*D(a+i+k)) for ci>0",
        "rescaling_identity": "ci*fi+cj*fj=(ci/Di)*(Mi+w*Mj); the quintic and its denominator scale by (ci/Di)^5",
        "polynomial": "Phi(Mi+w*Mj), with the general homogeneous quintic Phi",
        "degrees": generic_degrees, "terms": len(generic_rows),
        "positive_terms": len(generic_rows), "negative_count": 0,
        "minimum_coefficient": min(row[-1] for row in generic_rows),
        "terms_by_weight_degree": by_weight,
        "sha256": generic_hash,
        "hash_format": "Sorted rows [a_power,i_power,k_power,x_power,w_power,integer_coefficient]; compact ASCII JSON without newline",
        "exact_component_divisions_verified": True,
        "pure_weight_coefficients_match_component_quintics": True,
        "symbolic_weight_rescaling_verified": True,
        "strictness": "For ci>0, nonnegative mixed coefficients give Phi(Mi+w*Mj)>=Phi(Mi)>0 by the monomial certificate; ci=0 is the pure upper monomial",
        "direct_integral_audit": {
            "a": 1, "x": 1, "i": 0, "k": 2, "ci": 1, "cj": 1,
            "Di": str(di_value), "Dj": str(dj_value), "w": str(rational_point[-1]),
            "all_nine_rescaled_forms_match_retained_integrals": True,
            "gap": str(generic_gap), "quintic_ratio_matches": True},
        "limitations": "Does not cover three or more distinct support indices or all-size next ULC"}
    print("Third-sector LC: 1074 monomial, 2890 support-{0,2}, and 296701 arbitrary-two-support positive numerator coefficients; all audits passed.", flush=True)
    return {
        "scope": "Base certificate for monomials and arbitrary two-power mixtures; the all-support theorem is certified separately by third_sector_lc_all",
        "definitions": {
            "v": "-z*log(z)/(1-z)", "D": "z*partial_z",
            "eta_s": "normalized z^(a-1)*(1-z)^s*(1+x*z)^(3-s)*P(z), s=1,2,3",
            "u": "E_eta1[v]",
            "b": "(2*u*E_eta2[(D+a)v]-a*E_eta2[v^2])/(a+1)",
            "c": "(a^2*E_eta3[v^3]-3*a*u*E_eta3[(D+a)v^2]+3*b*E_eta3[(D+a)(D+a+1)v])/((a+1)*(a+2))"},
        "general_quintic": {
            "identity_verified": True, "homogeneous_degree_in_P_coefficients": 5,
            "B": "2*N1*M2-a*Z1*H2",
            "numerator": "(a+2)*Z3*B^2-a^2*(a+1)*N1*Z1*Z2^2*H3+3*a*(a+1)*N1^2*Z2^2*J3-3*N1*B*Z2*K3",
            "denominator": "(a+1)^2*(a+2)*Z1^2*Z2^2*Z3"},
        "monomial": component_record, "two_support": two_support_record,
        "arbitrary_two_support": generic_record}


def third_sector_lc_all_certificate():
    """Six three-support quintic coefficients closing the next LC inequality.

    Build in independent A_j, then compose A_j=a+i+successive gaps.
    Keep full, unreduced polynomials and stream their coefficients/hashes.
    The one/two-support certificate is a separate input rerun by this flag.
    The audited compact-polytope reduction extends three supports to every
    finite nonnegative P. This proves the next ULC inequality, not all ranks.
    """
    from collections import Counter
    from hashlib import sha256
    from itertools import permutations, product
    from math import comb, factorial, prod
    import json
    try:
        from flint import fmpz_mpoly_ctx
    except ImportError as error:
        raise SystemExit("--third-sector-lc-all requires optional python-flint") from error

    names = ("Z1", "N1", "Z2", "M2", "H2", "Z3", "H3", "J3", "K3")
    # Counts/hashes identify the original full polynomials. Rational point
    # coefficients come from the independent differentiated-integral audit.
    specifications = (
        ((3, 1, 1), 340009, [37, 34, 16, 8, 6],
         "3802c78b6e2ef443b429a561b5b5bed79baaf75d751d037d768973dd91dd1334",
         "156371144851/25604812800000", 4819528811309696318177280),
        ((1, 3, 1), 445471, [37, 34, 30, 8, 6],
         "54f3a6a333092b7390ecce4893c9bc59db0185d0043ea08fb61b8470eccce16e",
         "105514794817/171460800000000", 329348457925448564736000000),
        ((1, 1, 3), 666730, [37, 34, 30, 23, 6],
         "7a068e98bcb9116a300e128543055dd6e2c19c48bbf6532d2186b89e34b4a087",
         "261863696477629/3842322220800000000", 189241099455663237144792268800),
        ((2, 2, 1), 416123, [37, 34, 23, 8, 6],
         "9b8f7a2b2145cc6dd69e08a892b86168935594c2155e8d40071c476c2512ccc5",
         "12059370474181/3840721920000000", 64528284416470567354368000),
        ((2, 1, 2), 559671, [37, 34, 23, 16, 6],
         "0917bfa67b84d8116216f582979c1a08cd7b7548e5f6585ce01fa6346493c0e6",
         "1465850298759937/1317367618560000000", 1647154724245456853202370560),
        ((1, 2, 2), 614743, [37, 34, 30, 16, 6],
         "cdad14d644561d6347b1a03aab0e7503ba32392eb4f2fafa80a85bee558e386e",
         "1033985458127261/3293419046400000000", 12102859014066631534510080000),
    )
    patterns = {row[0] for row in specifications}
    assert patterns == {powers for powers in product(range(1, 6), repeat=3) if sum(powers) == 5}

    def seven_products(a):
        return (
            (4*(a+2), ("Z3", "N1", "N1", "M2", "M2")),
            (-4*a*(a+2), ("Z3", "N1", "M2", "Z1", "H2")),
            (a*a*(a+2), ("Z3", "Z1", "Z1", "H2", "H2")),
            (-a*a*(a+1), ("N1", "Z1", "Z2", "Z2", "H3")),
            (3*a*(a+1), ("N1", "N1", "Z2", "Z2", "J3")),
            (-6, ("N1", "N1", "M2", "Z2", "K3")),
            (3*a, ("N1", "Z1", "H2", "Z2", "K3")),
        )

    def normalized_gap(forms, a):
        Z1, N, Z2, M2, H2, Z3, H3, J3, K3 = (forms[name] for name in names)
        u = N/Z1
        b = (2*u*M2-a*H2)/((a+1)*Z2)
        c = (a*a*H3-3*a*u*J3+3*b*K3)/((a+1)*(a+2)*Z3)
        return b*b-u*c, (a+1)**2*(a+2)*Z1**2*Z2**2*Z3

    # Popper's independent audit: derive the 27 moments from v and its
    # derivatives, not from the typed rational kernels or the seven products.
    sa = sp.Symbol("a", positive=True)
    formal = dict(zip(names, sp.symbols(" ".join(names))))
    formal_gap, formal_denominator = normalized_gap(formal, sa)
    seven = sum(prefactor*sp.prod(formal[name] for name in factors)
                for prefactor, factors in seven_products(sa))
    assert sp.cancel(formal_gap*formal_denominator-seven) == 0
    z, ell = sp.symbols("z ell")
    voltage = z*ell/(1-z)

    def D(expression):
        return sp.cancel(z*sp.diff(expression, z)-sp.diff(expression, ell))

    def integral(expression):
        polynomial = sp.Poly(sp.cancel(expression), z, ell)
        assert all(power >= 0 and logs >= 0 for (power, logs), _ in polynomial.terms())
        return sp.cancel(sum(coefficient*factorial(logs)/sp.Integer(power+1)**(logs+1)
                             for (power, logs), coefficient in polynomial.terms()))

    kernels = (
        (1-z)*(1+z)**2,
        (1-z)*(1+z)**2*voltage,
        (1-z)**2*(1+z),
        (1-z)**2*(1+z)*(voltage+D(voltage)),
        (1-z)**2*(1+z)*voltage**2,
        (1-z)**3,
        (1-z)**3*voltage**3,
        (1-z)**3*(voltage**2+D(voltage**2)),
        (1-z)**3*(2*voltage+3*D(voltage)+D(D(voltage))),
    )
    indices = (0, 1, 3)
    direct = [dict(zip(names, (integral(z**index*kernel) for kernel in kernels)))
              for index in indices]
    assert all(value > 0 for forms in direct for value in forms.values())
    weights = sp.symbols("p0 p1 p2")
    mixture = {name: sum(weights[j]*direct[j][name] for j in range(3)) for name in names}
    direct_gap, direct_denominator = normalized_gap(mixture, sp.S.One)
    direct_quintic = sp.Poly(sp.cancel(direct_gap*direct_denominator), *weights)
    assert len(direct_quintic.terms()) == 21
    assert all(sum(powers) == 5 for powers, _ in direct_quintic.terms())

    # Independent A_j keep intermediate polynomial products small. The
    # positive component denominator clears all nine moments at once.
    ctx = fmpz_mpoly_ctx.get(["a", "A0", "A1", "A2", "x"])
    a, A0, A1, A2, x = ctx.gens()
    zero = ctx.constant(0)

    def component(A):
        index = A-a
        denominator = A*(A+1)**2*(A+2)**3*(A+3)**4

        def over(d):
            quotient = denominator/d
            assert quotient*d == denominator
            return quotient

        def Z(s, j):
            return factorial(s)*over(prod(a+j+h for h in range(s+1)))

        def N(j):
            return over((a+j+1)**2)

        def M(j):
            return (j+2)*over((a+j+2)**2)-j*over((a+j+1)**2)

        def H(j):
            return 2*over((a+j+2)**3)

        def V(j):
            return 6*over((a+j+3)**4)

        def J(j):
            return 2*((j+3)*over((a+j+3)**3)-j*over((a+j+2)**3))

        def K(j):
            return ((j+2)*(j+3)*over((a+j+3)**2)-2*j*(j+2)*over((a+j+2)**2)
                    +j*(j-1)*over((a+j+1)**2))

        def T(r, function):
            return sum((comb(r, h)*x**h*function(index+h) for h in range(r+1)), zero)

        return dict(zip(names, (T(2, lambda j: Z(1, j)), T(2, N),
                               T(1, lambda j: Z(2, j)), T(1, M), T(1, H),
                               Z(3, index), V(index), J(index), K(index))))

    components = [component(A) for A in (A0, A1, A2)]

    def at_integer_point(polynomial, point):
        powers = [[value**j for j in range(int(degree)+1)]
                  for value, degree in zip(point, polynomial.degrees())]
        return sum(int(polynomial.coefficient(index))*prod(
            powers[j][int(exponent)] for j, exponent in enumerate(polynomial.monomial(index)))
                   for index in range(len(polynomial)))

    component_point = [1, 1, 2, 4, 1]
    scales = [(index+1)*(index+2)**2*(index+3)**3*(index+4)**4 for index in indices]
    assert scales == [27648, 720000, 51861600]
    for j in range(3):
        for name in names:
            assert sp.Rational(at_integer_point(components[j][name], component_point), scales[j]) == direct[j][name]

    target = fmpz_mpoly_ctx.get(["a", "i", "d", "e", "x"])
    aa, ii, dd, ee, xx = target.gens()
    compose_with = (aa, aa+ii, aa+ii+dd, aa+ii+dd+ee, xx)
    point = [1, 0, 1, 2, 1]
    records = []
    for pattern, expected_count, expected_degrees, expected_hash, expected_direct, expected_value in specifications:
        label = "".join(map(str, pattern))
        print(f"Third-sector pattern {label}: constructing full coefficient in independent A_j.", flush=True)
        labels = tuple(j for j, count in enumerate(pattern) for _ in range(count))
        assignments = sorted(set(permutations(labels)))
        assignment_count = factorial(5)//prod(factorial(count) for count in pattern)
        assert len(assignments) == assignment_count
        total = zero
        group_counts = []
        for prefactor, factors in seven_products(a):
            grouped = Counter(tuple(sorted(zip(assignment, factors))) for assignment in assignments)
            assert sum(grouped.values()) == assignment_count
            subtotal = zero
            for key, multiplicity in grouped.items():
                value = ctx.constant(multiplicity)
                for j, name in key:
                    value = value*components[j][name]
                subtotal = subtotal+value
            total = total+prefactor*subtotal
            group_counts.append(len(grouped))
        assert group_counts == ([8, 20, 8, 13, 8, 13, 20] if max(pattern) == 3
                                else [11, 30, 11, 18, 11, 18, 30])
        # No factor is removed: counts and hashes refer to the full polynomial.
        answer = total.compose(*compose_with, ctx=target)
        del total
        degrees = list(map(int, answer.degrees()))
        assert len(answer) == expected_count and degrees == expected_degrees
        point_powers = [[value**j for j in range(degree+1)] for value, degree in zip(point, degrees)]
        digest = sha256()
        digest.update(b"[")
        previous, minimum, point_value, checked = None, None, 0, 0
        for offset, index in enumerate(range(len(answer)-1, -1, -1)):
            powers = tuple(map(int, answer.monomial(index)))
            coefficient = int(answer.coefficient(index))
            if previous is not None:
                assert previous < powers
            previous = powers
            assert coefficient > 0, (pattern, powers, coefficient)
            minimum = coefficient if minimum is None else min(minimum, coefficient)
            if offset:
                digest.update(b",")
            digest.update(json.dumps([*powers, coefficient], separators=(",", ":")).encode("ascii"))
            point_value += coefficient*prod(point_powers[j][power] for j, power in enumerate(powers))
            checked += 1
        digest.update(b"]")
        assert checked == expected_count and digest.hexdigest() == expected_hash
        direct_coefficient = direct_quintic.coeff_monomial(pattern)
        denominator_scale = prod(scales[j]**pattern[j] for j in range(3))
        assert direct_coefficient == sp.Rational(expected_direct) > 0
        assert point_value == expected_value == direct_coefficient*denominator_scale
        records.append({
            "pattern": list(pattern), "terms": checked, "positive_terms": checked,
            "negative_count": 0, "degrees": degrees, "minimum_coefficient": minimum,
            "sha256": digest.hexdigest(), "removed_factor": "1",
            "streamed_ascending_monomial_order_verified": True,
            "assignments_per_product": assignment_count, "group_counts_by_product": group_counts,
            "direct_original_weight_coefficient": str(direct_coefficient),
            "positive_component_denominator_scale": str(denominator_scale),
            "cleared_coefficient_at_audit_point": str(point_value),
            "full_polynomial_evaluation_matches_direct_integrals": True})
        print(f"Third-sector pattern {label}: {checked} positive coefficients; hash and direct integral coefficient match.", flush=True)
        del answer
    total_terms = sum(record["terms"] for record in records)
    assert total_terms == 3042747
    return {
        "status": "proved by audited support reduction and exact one/two/three-support coefficient certificates",
        "scope": "a>0,x>=0; every nonzero finite nonnegative combination P of nonnegative real powers; strict b^2>u*c",
        "physical_consequence": "P=(1+x*z)^(n-3), n>=3: strict next homogeneous ULC inequality; not all-rank ULC or NA",
        "definitions_reference": "third_sector_lc.definitions",
        "prior_certificate_dependency": "third_sector_lc is rerun unchanged by --third-sector-lc-all",
        "variables": ["a", "i", "d", "e", "x"], "indices": ["i", "i+d", "i+d+e"],
        "component_denominator": "D(A)=A*(A+1)^2*(A+2)^3*(A+3)^4",
        "weight_rescaling": "q_j=p_j/D(a+k_j); sum p_j*f_j=sum q_j*M_j, M_j=D(a+k_j)*f_j",
        "construction": "Seven quintic products; distinct multiset assignments with grouped multiplicities; independent A_j then A_j=a+i+successive gaps",
        "hash_format": "Compact ASCII JSON of ascending [a_power,i_power,d_power,e_power,x_power,integer_coefficient] rows; streamed via monomial(index)/coefficient(index)",
        "full_unreduced_polynomials": True, "all_six_multiplicity_patterns_verified": True,
        "positive_coefficient_total": total_terms, "negative_coefficient_total": 0, "cases": records,
        "independent_integral_audit": {
            "status": "received and reproduced", "a": 1, "x": 1, "indices": list(indices),
            "parameter_point": point, "component_denominators": list(map(str, scales)),
            "direct_components": [{"index": index, "moments": {name: str(direct[j][name]) for name in names}}
                                  for j, index in enumerate(indices)],
            "all_27_components_match": True, "all_six_weight_coefficients_match": True,
            "seven_products_equal_normalized_gap_symbolically": True,
            "full_composed_polynomials_match_the_independent_point_values": True},
        "support_reduction": {
            "status": "independently audited analytic input; no numerical LP test used",
            "constraints": ["Z3=1", "N1-u*Z1=0", "2*u*M2-a*H2-(a+1)*b*Z2=0"],
            "argument": "On the original finite support, positive Z3 coefficients make the feasible set compact. These three linear equalities fix u,b. Maximize c=(a^2*H3-3*a*u*J3+3*b*K3)/((a+1)*(a+2)); a maximizer uses at most three support points. Since u>0, its strict gap bounds the original gap below.",
            "strictness": "Nonzero mixtures have a positive pure fifth-power contribution; one/two-support boundary cases are covered by the prior certificate.",
            "limitations": "Finite supports only; no infinite-series closure or lower-sector ULC assertion"}}


def fourth_sector_lc_certificate():
    """Two exact fourth-sector certificates and their normalization audits.

    The monomial problem has eta_r proportional to
    z^(a+k-1)(1-z)^r(1+x*z)^(4-r), a>0,k,x>=0. The separate two-power
    result fixes a=1,x=0 and P=1+C*z^k, k,C>=0. Both prove b3^2>b2*b4
    on their stated domains, with no arbitrary-mixture or all-rank claim.
    No earlier certificate flags are invoked. The physical check reads the
    already saved n=4 rank law, with loopless lambda=1 and holding x=5/4.
    """
    from hashlib import sha256
    from math import comb, factorial
    from pathlib import Path
    import json
    try:
        from flint import fmpz_mpoly_ctx
    except ImportError as error:
        raise SystemExit("--fourth-sector-lc requires optional python-flint") from error

    sa, sk, z, ell = sp.symbols("a k z ell")
    voltage = z*ell/(1-z)  # ell=-log(z)

    def derivative(expression):
        return sp.cancel(z*sp.diff(expression, z)-sp.diff(expression, ell))

    def integral(polynomial, exponent):
        return sp.cancel(sum(coefficient*factorial(logs)/(exponent+power)**(logs+1)
                             for (power, logs), coefficient in polynomial.terms()))

    direct_polynomials, direct_moments, checked = {}, {}, []
    for r in range(1, 5):
        for t in range(r+1):
            m = r-t
            expression = voltage**m
            for j in range(t):
                expression = sp.cancel(derivative(expression)+(sa+j)*expression)
            polynomial = sp.Poly(sp.cancel((1-z)**r*expression), z, ell)
            direct = integral(polynomial, sa+sk)
            adjoint = factorial(m)*sum(
                (-1)**(t+j)*comb(t, j)*sp.ff(sk, t-j)*sp.rf(sk+m+1, j)
                /(sa+sk+m+j)**(m+1) for j in range(t+1))
            assert sp.cancel(direct-adjoint) == 0, (r, t)
            if t == r:
                assert sp.cancel(direct-sp.rf(sa, r)*factorial(r)/sp.rf(sa+sk, r+1)) == 0
            direct_polynomials[r, t], direct_moments[r, t] = polynomial, direct
            checked.append([r, t])

    def direct_rows(P, inverse_rate):
        # Independent original, unscaled equations at a=1. Keep the tilt
        # inside the differential integral and solve one four-row system.
        rows, normalizers = {}, {}
        matrix, rhs = sp.zeros(4), sp.zeros(4, 1)
        for r in range(1, 5):
            weight = sp.sympify(P)*(1+inverse_rate*z)**(4-r)
            normalizers[r] = integral(sp.Poly(sp.expand(weight*(1-z)**r), z, ell), sp.S.One)
            for t in range(r+1):
                polynomial = sp.Poly(sp.expand(weight*direct_polynomials[r, t].as_expr().subs(sa, 1)), z, ell)
                rows[r, t] = integral(polynomial, sp.S.One)
            assert rows[r, r] == factorial(r)*normalizers[r]
            rhs[r-1] = -rows[r, 0]
            for t in range(1, r+1):
                matrix[r-1, t-1] = (-1)**t*comb(r, t)*rows[r, t]
        return [sp.S.One, *(matrix.inv()*rhs)], rows

    # All-parameter monomial numerator, using the original integer basis.
    ctx = fmpz_mpoly_ctx.get(["a", "k", "x"])
    a, k, x = ctx.gens()
    zero, one = ctx.constant(0), ctx.constant(1)
    A = a+k

    def product(values):
        answer = one
        for value in values:
            answer *= value
        return answer

    def rising(value, length):
        return product(value+j for j in range(length))

    def falling(value, length):
        return product(value-j for j in range(length))

    common = product((A+j)**(j+1) for j in range(5))

    def over(denominator):
        quotient = common/denominator
        assert quotient*denominator == common
        return quotient

    Z, J, H = {}, {}, {}
    for r in range(1, 5):
        Z[r] = sum((comb(4-r, j)*x**j*factorial(r)
                    *over(product(A+j+h for h in range(r+1)))
                    for j in range(5-r)), zero)
        H[r] = rising(a+1, r-1)*Z[r]
        for t in range(r):
            m = r-t
            J[r, t] = sum((
                comb(4-r, j)*x**j*factorial(m)*(-1)**(t+h)*comb(t, h)
                *falling(k+j, t-h)*rising(k+j+m+1, h)
                *over((A+j+m+h)**(m+1))
                for j in range(5-r) for h in range(t+1)), zero)
    N, Q = {0: one}, {0: one}
    for r in range(1, 5):
        N[r] = sum(((-1)**(r+t+1)*comb(r, t)*a**(r-t-1)
                    *N[t]*J[r, t]*product(H[h] for h in range(t+1, r))
                    for t in range(r)), zero)
        Q[r] = Q[r-1]*H[r]
    Gamma = H[4]*N[3]**2-H[3]*N[2]*N[4]
    monomial_data = Gamma.to_dict()
    assert len(monomial_data) == 52302
    assert all(coefficient > 0 for coefficient in monomial_data.values())
    monomial_rows = sorted([*map(int, powers), int(coefficient)]
                          for powers, coefficient in monomial_data.items())
    monomial_hash = sha256(json.dumps(monomial_rows, separators=(",", ":")).encode()).hexdigest()
    assert monomial_hash == "3a7653b6cbc5c73f9c91d7c7c0a224ecf1b2d166a1e8cf5f4a6251da7895dc60"
    assert list(map(int, Gamma.degrees())) == [89, 80, 12]
    base = sum(int(coefficient)*sa**int(powers[0])
               for powers, coefficient in monomial_data.items()
               if powers[1] == 0 and powers[2] == 0)
    expected_base = (3456*sa**6*(sa+1)**10*(sa+2)**15*(sa+3)**20*(sa+4)**30
                     *(sa**8+38*sa**7+615*sa**6+5452*sa**5+28975*sa**4
                       +94954*sa**3+188241*sa**2+207188*sa+97228))
    assert sp.expand(base-expected_base) == 0
    base_terms = sorted((int(powers[0]), int(coefficient))
                        for powers, coefficient in monomial_data.items()
                        if powers[1] == 0 and powers[2] == 0)
    assert list(base_terms[0]) == [6, 44262899067425937016127246828084788199424]

    point_b, unused = direct_rows(1, sp.S.One)
    expected_point_b = list(map(sp.Rational, ["1", "973/1560", "25193/74880", "649/4160", "779/12480"]))
    assert point_b == expected_point_b
    for r in range(1, 5):
        assert point_b[r] == sp.Rational(int(N[r](1, 0, 1)), int(Q[r](1, 0, 1)))
    point_gap = sp.factor(point_b[3]**2-point_b[2]*point_b[4])
    assert point_gap == sp.Rational(3119507, 934502400)
    assert point_gap == sp.Rational(int(Gamma(1, 0, 1)), int((Q[3]**2*H[4])(1, 0, 1)))

    # Two-power boundary: a=1,x=0,P=1+C*z^k. Positive rescaling
    # w=C*D4(1)/D4(k+1) makes each row kernel an integer polynomial.
    bk, bw = sp.symbols("k w", nonnegative=True)
    boundary_common = sp.Poly((bk+1)*sp.prod((bk+j+1)**(j+1) for j in range(1, 5)), bk, domain=sp.ZZ)
    D0 = boundary_common.eval(0)
    assert D0 == 86400000
    pzero, pone = sp.Poly(0, bk, bw, domain=sp.ZZ), sp.Poly(1, bk, bw, domain=sp.ZZ)
    L = {}
    for r in range(1, 5):
        for t in range(r+1):
            m = r-t
            cleared = sp.Poly(0, bk, domain=sp.ZZ)
            for j in range(t+1):
                numerator = factorial(m)*(-1)**(t+j)*comb(t, j)*sp.ff(bk, t-j)*sp.rf(bk+m+1, j)
                denominator = sp.Poly((bk+1+m+j)**(m+1), bk, domain=sp.ZZ)
                cleared += sp.Poly(numerator, bk, domain=sp.ZZ)*boundary_common.exquo(denominator)
            assert sp.cancel(cleared.as_expr()/boundary_common.as_expr()
                             -direct_moments[r, t].subs({sa: 1, sk: bk})) == 0
            assert cleared.eval(0) == D0*factorial(r)/sp.Integer(r+1)**(r-t+1)
            terms = {(degree[0], 1): coefficient for degree, coefficient in cleared.terms()}
            terms[0, 0] = cleared.eval(0)
            L[r, t] = sp.Poly.from_dict(terms, (bk, bw), domain=sp.ZZ)
    diagonals = {r: L[r, r] for r in range(1, 5)}
    numerators, products = {0: pone}, {0: pone}
    for r in range(1, 5):
        numerator = pzero
        for t in range(r):
            after_t = pone
            for j in range(t+1, r):
                after_t *= diagonals[j]
            numerator += (-1)**(r+t+1)*comb(r, t)*L[r, t]*numerators[t]*after_t
        numerators[r] = numerator
        products[r] = products[r-1]*diagonals[r]
    boundary_gap_numerator = numerators[3]**2*diagonals[4]-numerators[2]*numerators[4]*diagonals[3]
    boundary_gap_denominator = diagonals[1]**2*diagonals[2]**2*diagonals[3]**2*diagonals[4]
    content, primitive = boundary_gap_numerator.primitive()
    boundary_rows = sorted([*map(int, powers), int(coefficient)]
                           for powers, coefficient in primitive.terms())
    assert content == 165888 and len(boundary_rows) == 348
    assert all(row[-1] > 0 for row in boundary_rows)
    assert list(map(int, primitive.degree_list())) == [80, 7]
    boundary_hash = sha256(json.dumps(boundary_rows, separators=(",", ":")).encode()).hexdigest()
    assert boundary_hash == "fd88490e254c248c15ca8d2642981b1d8dfde7f024cfb070181561380c07566d"
    constant_term = primitive.coeff_monomial((0, 0))
    assert constant_term == 195186952277065728000000000000000000000000000000
    weight_counts = [sum(row[1] == j for row in boundary_rows) for j in range(8)]
    assert weight_counts == [1, 14, 27, 39, 51, 62, 73, 81]
    base_b = list(map(sp.Rational, ["1/2", "2/9", "17/192", "319/10000"]))
    for r, value in enumerate(base_b, 1):
        for variable in (bk, bw):
            assert sp.cancel(numerators[r].as_expr().subs(variable, 0)
                             /products[r].as_expr().subs(variable, 0)-value) == 0
    base_gap = sp.factor(base_b[2]**2-base_b[1]*base_b[3])
    assert base_gap == sp.Rational(17297, 23040000)

    # Parent's single interior check uses the original C, not w.
    boundary_b, boundary_moments = direct_rows(1+3*z**2, sp.S.Zero)
    rescale = 3*D0/boundary_common.eval(2)
    assert rescale == sp.Rational(100, 50421)
    boundary_point = {bk: 2, bw: rescale}
    for r in range(1, 5):
        assert boundary_b[r] == sp.cancel(numerators[r].as_expr().subs(boundary_point)
                                         /products[r].as_expr().subs(boundary_point))
        for t in range(r+1):
            assert sp.cancel(L[r, t].as_expr().subs(boundary_point)/D0) == boundary_moments[r, t]
    expected_boundary_b = list(map(sp.Rational, ["1", "7/12", "14023/46800", "572521/4212000",
                                                "133032315503/2427122880000"]))
    assert boundary_b == expected_boundary_b
    boundary_gap = sp.factor(boundary_b[3]**2-boundary_b[2]*boundary_b[4])
    assert boundary_gap == sp.Rational(699467433935021, 340768052352000000)
    assert boundary_gap == sp.cancel(boundary_gap_numerator.as_expr().subs(boundary_point)
                                     /boundary_gap_denominator.as_expr().subs(boundary_point)) > 0

    # Read the saved law; do not recompute ranks or invoke an older flag.
    holding_x = sp.Rational(5, 4)
    physical_b, unused = direct_rows(1, holding_x)
    rank_source = Path(__file__).with_name("complete_graph_certificates.json")
    matches = [case for case in json.loads(rank_source.read_text())
               if case["n"] == 4 and case["lambda"] == "1"]
    assert len(matches) == 1
    saved = matches[0]
    rank = list(map(sp.Rational, saved["rank"]))
    expected_physical_b = [rank[4-r]/(sp.binomial(4, r)*holding_x**r*rank[4]) for r in range(5)]
    assert physical_b == expected_physical_b

    def at_rational_point(polynomial, values):
        return sum((int(coefficient)*sp.prod(value**int(power) for value, power in zip(values, powers))
                    for powers, coefficient in polynomial.to_dict().items()), sp.S.Zero)

    for r in range(1, 5):
        assert physical_b[r] == (at_rational_point(N[r], (sp.S.One, sp.S.Zero, holding_x))
                                 /at_rational_point(Q[r], (sp.S.One, sp.S.Zero, holding_x)))
    physical_gap = sp.factor(physical_b[3]**2-physical_b[2]*physical_b[4])
    assert physical_gap == sp.Rational(2707657252163, 649428141796875)
    rank_gap = sp.factor((rank[1]/4)**2-(rank[2]/6)*rank[0])
    assert rank_gap == sp.Rational(saved["ULC_gaps"][0]) == holding_x**6*rank[4]**2*physical_gap

    print("Fourth-sector certificates: 52302 monomial and 348 boundary coefficients positive; direct moments and saved-law normalization passed.", flush=True)
    return {
        "status": "proved on the two stated domains by exact coefficient certificates",
        "target": "b3^2>b2*b4, with b0=1 and b_r defined by (35ac)",
        "limitations": "No arbitrary-mixture, general two-support a,x, or all-rank ULC claim",
        "kernel_audit": {
            "formula": "J_(r,t)(k)=(r-t)! sum_(ell=0)^t (-1)^(t+ell) binom(t,ell) falling(k,t-ell) rising(k+r-t+1,ell)/(a+k+r-t+ell)^(r-t+1)",
            "checked": checked, "all_fourteen_general_a_k_differential_moments_match": True,
            "diagonal_normalization": "J_(r,r)(k)=(a)_r*B(a+k,r+1)",
            "common_denominator": "D4(A)=A(A+1)^2(A+2)^3(A+3)^4(A+4)^5",
            "analytic_status": "Adjoint identity and endpoint terms independently audited for all r,t; finite exact checks retained here"},
        "monomial": {
            "scope": "a>0,k>=0,x>=0; eta_r proportional to z^(a+k-1)(1-z)^r(1+x*z)^(4-r)",
            "variables": ["a", "k", "x"],
            "gap_numerator": "Gamma=H4*N3^2-H3*N2*N4; b_r=N_r/Q_r, Q_r=product_(j<=r)H_j",
            "gap_denominator": "Q3^2*H4>0, H_r=(a+1)_(r-1)*D4*Z_r",
            "numerator_terms": len(monomial_data), "degrees": list(map(int, Gamma.degrees())),
            "negative_count": 0, "positive_coefficient_certificate": True,
            "numerator_sha256": monomial_hash,
            "numerator_at_k0_x0": str(sp.factor(base)), "strict_base_monomial": list(base_terms[0]),
            "strictness": "Gamma(a,k,x)>=Gamma(a,0,0)>0 for a>0,k,x>=0",
            "point_check": {"a": 1, "k": 0, "x": 1, "b": list(map(str, point_b)), "gap": str(point_gap)}},
        "two_power_boundary": {
            "scope": "a=1,x=0,P(z)=1+C*z^k, real k>=0, finite C>=0",
            "rescaling": "w=C*D4(1)/D4(k+1)", "D4_at_one": str(D0),
            "positive_denominator": "d1^2*d2^2*d3^2*d4, d_r=M_(r,r)(0)+w*M_(r,r)(k)>0",
            "gap_formula": "(n3^2*d4-n2*n4*d3)/(d1^2*d2^2*d3^2*d4)",
            "gap_numerator_content": str(content), "primitive_degrees_k_w": [80, 7],
            "primitive_terms": len(boundary_rows), "negative_count": 0,
            "coefficient_sha256": boundary_hash, "constant_term_k_w": str(constant_term),
            "terms_by_weight_degree": weight_counts,
            "base_cases_C_zero_or_k_zero": {"b1_to_b4": list(map(str, base_b)), "gap": str(base_gap),
                                            "identities_in_the_other_parameter": True},
            "all_fourteen_kernels_match_direct_moments": True,
            "unscaled_direct_integral_check": {
                "a": 1, "x": 0, "k": 2, "C": 3, "w": str(rescale),
                "b": list(map(str, boundary_b)), "gap": str(boundary_gap),
                "all_fourteen_moments_and_four_values_match": True}},
        "saved_physical_normalization": {
            "source": rank_source.name, "n": 4, "original_sleep_rate": "1",
            "holding_sleep_rate": "4/5", "holding_inverse_rate": str(holding_x),
            "a": 1, "P": "1", "b": list(map(str, physical_b)),
            "all_sector_values_match_saved_rank_law_and_monomial_polynomials": True,
            "physical_gap": str(physical_gap), "saved_rank_gap": str(rank_gap),
            "gap_conversion": "saved lowest normalized rank gap=x^6*r4^2*(b3^2-b2*b4)",
            "new_rank_law_computed": False}}


def fourth_sector_lc_two_support_certificate():
    """All six mixed coefficients for arbitrary fourth-sector two-power laws.

    P=c0*z^i+c1*z^(i+d), a>0, i,d,x>=0 and nonzero nonnegative weights.
    Build in independent A0,A1, extract one binary-weight coefficient at a
    time, then compose A0=a+i,A1=a+i+d. Stream all 7,586,990 coefficients.
    Pure endpoints use the separately saved monomial certificate. This
    optional calculation takes several minutes and does not rerun it.
    """
    import gc
    from fractions import Fraction
    from hashlib import sha256
    from math import comb, factorial, prod
    import json
    import resource
    try:
        from flint import fmpz_mpoly_ctx
    except ImportError as error:
        raise SystemExit("--fourth-sector-lc-two-support requires optional python-flint") from error

    # Full unreduced coefficients; never change these hashes by removing
    # positive factors. The last value is the original-weight integral
    # coefficient at the already audited a=1,i=0,d=1,x=1 point.
    specifications = (
        (1, 125220, [9, 72, 13, 12], 633642, [89, 80, 13, 12], 3456,
         "845ee2f59dd9ca01226c65d88ea656ac5f0e2358100318d7d9169173a765436b",
         "695951365939/2531250000000"),
        (2, 204421, [9, 61, 26, 12], 1047819, [89, 80, 26, 12], 3456,
         "e2bf3ab79b08ee2aeb045edf594198768ee0607fba04d4b6834b4fb7381be845",
         "343373377541411/1093500000000000"),
        (3, 242676, [9, 50, 38, 12], 1305905, [89, 80, 38, 12], 3456,
         "b0e57f83204725b728a96fcb2ae61eee2bbf7db86680d5e256e6b7c48efea06c",
         "167359636865849/874800000000000"),
        (4, 242676, [9, 38, 50, 12], 1466738, [89, 80, 50, 12], 6912,
         "b8cbaf9598bf9abd50c3b0f0b062e3e2138b566272bb46fc48ff298b196a276d",
         "2175492075657707/32805000000000000"),
        (5, 204421, [9, 26, 61, 12], 1548627, [89, 80, 61, 12], 3456,
         "fe1c471b51c061e9ff1fcbcaab7a4d91649177a908669a7bd5886b863d6852d2",
         "1509140610961987/118098000000000000"),
        (6, 125220, [9, 13, 72, 12], 1584259, [89, 80, 72, 12], 3456,
         "4ae27bdc4be0f1156c454da029a9a12cf7c70d6598bdb693bb9de86833c8e40d",
         "88833574108289/73811250000000000"),
    )
    assert [row[0] for row in specifications] == list(range(1, 7))
    z, ell, C = sp.symbols("z ell C")
    voltage = z*ell/(1-z)

    def derivative(expression):
        return sp.cancel(z*sp.diff(expression, z)-sp.diff(expression, ell))

    def integral(expression):
        polynomial = sp.Poly(sp.cancel(expression), z, ell)
        assert all(power >= 0 and logs >= 0 for (power, logs), _ in polynomial.terms())
        return sp.cancel(sum(coefficient*factorial(logs)/sp.Integer(power+1)**(logs+1)
                             for (power, logs), coefficient in polynomial.terms()))

    def direct_audit(a_value, x_value, indices):
        # Original, unrescaled moments from differentiation, followed by
        # a direct triangular linear solve for every original weight C.
        forms = [{}, {}]
        for r in range(1, 5):
            for j, index in enumerate(indices):
                weight = z**(a_value+index-1)*(1-z)**r*(1+x_value*z)**(4-r)
                forms[j][r, "Z"] = integral(weight)
            for t in range(r):
                kernel = voltage**(r-t)
                for h in range(t):
                    kernel = sp.cancel(derivative(kernel)+(a_value+h)*kernel)
                for j, index in enumerate(indices):
                    weight = z**(a_value+index-1)*(1-z)**r*(1+x_value*z)**(4-r)
                    forms[j][r, t] = integral(weight*kernel)
        mixture = {key: forms[0][key]+C*forms[1][key] for key in forms[0]}
        H = {r: sp.rf(a_value+1, r-1)*mixture[r, "Z"] for r in range(1, 5)}
        matrix, rhs = sp.zeros(4), sp.zeros(4, 1)
        for r in range(1, 5):
            rhs[r-1] = -a_value**r*mixture[r, 0]
            for t in range(1, r):
                matrix[r-1, t-1] = (-1)**t*comb(r, t)*a_value**(r-t)*mixture[r, t]
            matrix[r-1, r-1] = (-1)**r*sp.rf(a_value, r)*mixture[r, "Z"]
        b = [sp.S.One, *[sp.cancel(value) for value in matrix.inv()*rhs]]
        N, Q = {0: sp.S.One}, {0: sp.S.One}
        for r in range(1, 5):
            N[r] = sp.expand(sum((-1)**(r+t+1)*comb(r, t)*a_value**(r-t-1)
                                *N[t]*mixture[r, t]*sp.prod(H[h] for h in range(t+1, r))
                                for t in range(r)))
            Q[r] = Q[r-1]*H[r]
            assert sp.cancel(b[r]-N[r]/Q[r]) == 0
        Gamma = sp.Poly(sp.expand(H[4]*N[3]**2-H[3]*N[2]*N[4]), C)
        assert Gamma.degree() == 7 and len(Gamma.terms()) == 8
        assert all(value > 0 for value in Gamma.all_coeffs())
        assert sp.cancel((b[3]**2-b[2]*b[4])*Q[3]**2*H[4]-Gamma.as_expr()) == 0
        scales = [sp.prod((a_value+index+j)**(j+1) for j in range(5)) for index in indices]
        w = sp.Symbol("w")
        rescaled = sp.Poly(sp.expand(scales[0]**7*Gamma.as_expr().subs(C, w*scales[1]/scales[0])), w)
        coefficients = [Gamma.coeff_monomial(C**j) for j in range(8)]
        assert all(rescaled.coeff_monomial(w**j) == coefficients[j]*scales[0]**(7-j)*scales[1]**j
                   for j in range(8))
        return {"forms": forms, "b": b, "coefficients": coefficients, "scales": scales}

    unit = direct_audit(sp.S.One, sp.S.One, (0, 1))
    nonunit = direct_audit(sp.Integer(2), sp.Rational(1, 2), (1, 3))
    expected_nonunit = [
        "128803898074669271/595710034935750000000000",
        "44899783100541024876271/84739132931174104320000000000",
        "109536099746451751957509835531/205000910387096393170944000000000000",
        "234030647931497758147323431657/820003641548385572683776000000000000",
        "1017992333327130617067631102829/11808052438296752246646374400000000000",
        "1089610561448249325164670681773/75571535605099214378536796160000000000",
        "28529933501622834170623880017/24182891393631748601131774771200000000",
        "1359067605269311362149881/42991806922011997513123155148800000"]
    assert list(map(str, nonunit["coefficients"])) == expected_nonunit
    assert nonunit["scales"] == [130691232000, 14932726824960]
    nonunit_b = [sp.factor(value.subs(C, 2)) for value in nonunit["b"]]
    assert list(map(str, nonunit_b)) == ["1", "3672923/4504500", "1566683187287/2440538100000",
                                       "4913691023858387/10113589886400000", "293346822237130217/832203967795200000"]
    nonunit_gap = sp.factor(nonunit_b[3]**2-nonunit_b[2]*nonunit_b[4])
    assert nonunit_gap == sp.Rational(35975619260950347550813343034857,
                                     3682249214050525136578560000000000)

    ctx = fmpz_mpoly_ctx.get(["a", "A0", "A1", "x"])
    a, A0, A1, x = ctx.gens()
    zero, one = ctx.constant(0), ctx.constant(1)

    def product(values):
        value = one
        for item in values:
            value *= item
        return value

    def rising(value, length):
        return product(value+j for j in range(length))

    def falling(value, length):
        return product(value-j for j in range(length))

    def component(A):
        index = A-a
        denominator = product((A+j)**(j+1) for j in range(5))

        def over(q):
            quotient = denominator/q
            assert quotient*q == denominator
            return quotient

        Z, J, H = {}, {}, {}
        for r in range(1, 5):
            Z[r] = sum((comb(4-r, j)*x**j*factorial(r)
                        *over(product(A+j+h for h in range(r+1)))
                        for j in range(5-r)), zero)
            H[r] = rising(a+1, r-1)*Z[r]
            for t in range(r):
                m = r-t
                J[r, t] = sum((comb(4-r, j)*x**j*factorial(m)*(-1)**(t+h)*comb(t, h)
                               *falling(index+j, t-h)*rising(index+j+m+1, h)
                               *over((A+j+m+h)**(m+1))
                               for j in range(5-r) for h in range(t+1)), zero)
        return denominator, Z, J, H

    components = [component(A) for A in (A0, A1)]

    def at_rational_point(polynomial, point):
        powers = [[value**j for j in range(int(degree)+1)]
                  for value, degree in zip(point, polynomial.degrees())]
        return sum((Fraction(int(polynomial.coefficient(index)))*prod(
            powers[j][int(power)] for j, power in enumerate(polynomial.monomial(index)))
                    for index in range(len(polynomial))), Fraction(0))

    for audit, point in [(unit, tuple(map(Fraction, [1, 1, 2, 1]))),
                         (nonunit, (Fraction(2), Fraction(3), Fraction(5), Fraction(1, 2)))]:
        for j in range(2):
            assert at_rational_point(components[j][0], point) == audit["scales"][j]
            for r in range(1, 5):
                assert at_rational_point(components[j][1][r], point)/int(audit["scales"][j]) == audit["forms"][j][r, "Z"]
                for t in range(r):
                    assert at_rational_point(components[j][2][r, t], point)/int(audit["scales"][j]) == audit["forms"][j][r, t]
    print("Fourth-sector two-support: both direct-integral audits and all 28 component entries per point passed.", flush=True)

    target = fmpz_mpoly_ctx.get(["a", "i", "d", "x"])
    aa, ii, dd, xx = target.gens()

    def selected_coefficient(specification):
        weight, expected_independent, expected_independent_degrees, expected_count, expected_degrees, expected_minimum, expected_hash, expected_direct = specification
        cap = min(weight, 4)
        H = {r: [components[j][3][r] for j in range(2)] for r in range(1, 5)}
        J = {(r, t): [components[j][2][r, t] for j in range(2)]
             for r in range(1, 5) for t in range(r)}

        def add_series(left, right):
            length = min(cap+1, max(len(left), len(right)))
            return [(left[j] if j < len(left) else zero)+(right[j] if j < len(right) else zero)
                    for j in range(length)]

        def mul_series(left, right):
            result = [zero for _ in range(min(cap+1, len(left)+len(right)-1))]
            for i, first in enumerate(left):
                for j, second in enumerate(right):
                    if i+j < len(result):
                        result[i+j] = result[i+j]+first*second
            return result

        N = {0: [one]}
        for r in range(1, 5):
            total = [zero]
            for t in range(r):
                term = mul_series(N[t], J[r, t])
                for h in range(t+1, r):
                    term = mul_series(term, H[h])
                factor = (-1)**(r+t+1)*comb(r, t)*a**(r-t-1)
                total = add_series(total, [factor*value for value in term])
            N[r] = total

        def coefficient_of_product(factors):
            result = zero

            def visit(position, remaining, accumulated):
                nonlocal result
                if position == len(factors):
                    if remaining == 0:
                        result = result+accumulated
                    return
                for j, value in enumerate(factors[position]):
                    if j <= remaining:
                        visit(position+1, remaining-j, accumulated*value)

            visit(0, weight, one)
            return result

        total = coefficient_of_product([H[4], N[3], N[3]])-coefficient_of_product([H[3], N[2], N[4]])
        independent_count = len(total)
        independent_degrees = list(map(int, total.degrees()))
        assert independent_count == expected_independent and independent_degrees == expected_independent_degrees
        del N, H, J, term
        gc.collect()
        answer = total.compose(aa, aa+ii, aa+ii+dd, xx, ctx=target)
        del total
        gc.collect()
        degrees = list(map(int, answer.degrees()))
        assert len(answer) == expected_count and degrees == expected_degrees
        print(f"Fourth-sector two-support w^{weight}: streaming {expected_count} coefficients.", flush=True)

        direct = unit["coefficients"][weight]
        assert direct == sp.Rational(expected_direct)
        scale = unit["scales"][0]**(7-weight)*unit["scales"][1]**weight
        expected_point = direct*scale
        assert expected_point.q == 1 and answer(1, 0, 1, 1) == int(expected_point)
        nonunit_scale = nonunit["scales"][0]**(7-weight)*nonunit["scales"][1]**weight
        nonunit_expected = nonunit["coefficients"][weight]*nonunit_scale
        # At (a,i,d,x)=(2,1,2,1/2), multiplying by 2^12 makes each
        # monomial integral. This checks the full composed polynomial at
        # the parent's nonunit-a point during the same coefficient pass.
        digest = sha256()
        digest.update(b"[")
        previous, minimum, checked_count = None, None, 0
        unit_value, nonunit_scaled_value = 0, 0
        for offset, index in enumerate(range(len(answer)-1, -1, -1)):
            powers = tuple(map(int, answer.monomial(index)))
            coefficient = int(answer.coefficient(index))
            if previous is not None:
                assert previous < powers
            previous = powers
            assert coefficient > 0, (weight, powers, coefficient)
            minimum = coefficient if minimum is None else min(minimum, coefficient)
            if offset:
                digest.update(b",")
            digest.update(json.dumps([*powers, coefficient], separators=(",", ":")).encode("ascii"))
            if powers[1] == 0:
                unit_value += coefficient
            nonunit_scaled_value += coefficient << (powers[0]+powers[2]+12-powers[3])
            checked_count += 1
            if checked_count % 500000 == 0:
                print(f"Fourth-sector two-support w^{weight}: {checked_count} positive coefficients checked.", flush=True)
        digest.update(b"]")
        assert checked_count == expected_count and minimum == expected_minimum
        assert digest.hexdigest() == expected_hash
        assert unit_value == expected_point
        assert sp.Rational(nonunit_scaled_value, 2**12) == nonunit_expected
        record = {
            "weight_degree": weight, "pattern": [7-weight, weight],
            "independent_A_terms": independent_count, "independent_A_degrees": independent_degrees,
            "terms": checked_count, "degrees": degrees,
            "streamed_positive_terms_before_stop": checked_count,
            "first_negative": None, "all_coefficients_positive": True,
            "minimum_positive_coefficient": minimum, "sha256": digest.hexdigest(),
            "hash_format": "Compact ASCII JSON of ascending [a,i,d,x,coefficient] rows, streamed via FLINT indices",
            "direct_original_weight_coefficient_at_1_0_1_1": str(direct),
            "component_denominator_scale_at_audit_point": str(scale),
            "cleared_coefficient_at_audit_point": str(expected_point),
            "direct_integral_check_passed": True,
            "direct_original_weight_coefficient_at_2_1_2_half": str(nonunit["coefficients"][weight]),
            "cleared_coefficient_at_nonunit_point": str(nonunit_expected),
            "full_polynomial_nonunit_integral_check_passed": True,
            "peak_rss_mib": round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024, 2)}
        del answer
        gc.collect()
        print(f"Fourth-sector two-support w^{weight}: count, hash and both integral coefficients passed.", flush=True)
        return record

    records = {}
    for specification in specifications:
        weight = specification[0]
        print(f"Fourth-sector two-support w^{weight}: constructing the unreduced coefficient.", flush=True)
        records[str(weight)] = selected_coefficient(specification)
    total_terms = sum(record["terms"] for record in records.values())
    assert total_terms == 7586990
    return {
        "formula_id": "Gamma=H4*N3^2-H3*N2*N4;degree7;D4",
        "scope": "P=c0*z^i+c1*z^(i+d), a>0, i,d,x>=0, c0,c1>=0 not both zero",
        "status": "strict b3^2>b2*b4 by all six mixed coefficients and the saved monomial endpoints",
        "base_certificate_scope": "Unchanged parent fields describe the base certificates; this subkey adds arbitrary two supports",
        "variables": ["a", "i", "d", "x"],
        "positive_weight_change": "w=c1*D4(a+i)/(c0*D4(a+i+d)), for c0>0",
        "moment_rescaling": "component_j = D4(A_j) times original linear moment vector",
        "gap_numerator": "Gamma=H4*N3^2-H3*N2*N4; homogeneous degree seven in the original moments",
        "gap_denominator": "Q3^2*H4>0, Qr=product_(h<=r)Hh; Hr=(a+1)_(r-1)*Zr",
        "pure_endpoints": "strict monomial theorem; endpoint coefficients nonnegative after k=i or k=i+d",
        "monomial_dependency": "fourth_sector_lc.monomial; validated saved hash, not rerun by this flag",
        "monomial_dependency_hash": "3a7653b6cbc5c73f9c91d7c7c0a224ecf1b2d166a1e8cf5f4a6251da7895dc60",
        "construction": "Independent A0,A1; truncated nonnegative-power series for N_r; extract each mixed weight coefficient before substitution",
        "series_truncation": "For weight degree j, cap=min(j,4) retains every N_r coefficient that can contribute; N_r has degree r<=4",
        "full_unreduced_polynomials": True, "mixed": records,
        "positive_coefficient_total": total_terms, "all_six_mixed_complete_and_positive": True,
        "negative_coefficient_found": False,
        "nonunit_audit": {
            "a": "2", "x": "1/2", "indices": [1, 3], "symbolic_original_weight": "C",
            "all_28_independent_component_moments_match_at_nonunit_a": True,
            "triangular_matrix_solution_equals_recursive_numerators": True,
            "degree_seven_gap_identity_and_all_weight_scalings_verified": True,
            "direct_original_weight_coefficients_ascending": list(map(str, nonunit["coefficients"])),
            "component_scales": list(map(str, nonunit["scales"])),
            "point_C": 2, "point_b": list(map(str, nonunit_b)), "point_gap": str(nonunit_gap),
            "all_six_full_polynomials_match_the_nonunit_integral_coefficients": True},
        "unit_audit": {"a": 1, "x": 1, "indices": [0, 1],
                       "all_28_component_moments_match": True,
                       "direct_original_weight_coefficients_ascending": list(map(str, unit["coefficients"])),
                       "component_scales": list(map(str, unit["scales"]))},
        "limitations": "Two supports only; no arbitrary mixture or all-rank ULC assertion"}


def fourth_sector_chebyshev_obstruction_certificate():
    """A fixed, exact obstruction to the five-function T-system shortcut.

    At a=1,x=0,P=z, solve the first three original differentiated-integral
    equations for b1,b2,b3. The functions (Z4,f1,f2,f3,g4) then have ordered
    evaluation determinants of opposite signs on k>=0. This rules out a
    Chebyshev-system argument using this family; it does not disprove a
    three-support reduction by other means or any sector LC inequality.
    Only small univariate rational functions are computed by this flag.
    """
    from math import comb, factorial

    k, z, ell = sp.symbols("k z ell")
    A = k+1
    voltage = z*ell/(1-z)  # ell=-log(z), so z*d/dz = z*partial_z-partial_ell

    def integral(expression):
        polynomial = sp.Poly(sp.cancel(expression), z, ell)
        return sp.cancel(sum(coefficient*factorial(logs)/(A+power)**(logs+1)
                             for (power, logs), coefficient in polynomial.terms()))

    # Differentiate before integrating. This does not use the adjoint
    # rational-kernel formula to obtain the fixed b values.
    direct, kernels, normalizers = {}, {}, {}
    for r in range(1, 5):
        normalizers[r] = factorial(r)/sp.rf(A, r+1)
        assert sp.cancel(integral((1-z)**r)-normalizers[r]) == 0
        for t in range(r+1):
            m = r-t
            expression = voltage**m
            for j in range(t):
                expression = sp.cancel(z*sp.diff(expression, z)-sp.diff(expression, ell)
                                       +(1+j)*expression)
            direct[r, t] = integral((1-z)**r*expression)
            kernels[r, t] = factorial(m)*sum(
                (-1)**(t+j)*comb(t, j)*sp.ff(k, t-j)*sp.rf(k+m+1, j)
                /(A+m+j)**(m+1) for j in range(t+1))
            assert sp.cancel(direct[r, t]-kernels[r, t]) == 0, (r, t)
        assert sp.cancel(direct[r, r]-factorial(r)*normalizers[r]) == 0

    matrix, rhs = sp.zeros(3), sp.zeros(3, 1)
    for r in range(1, 4):
        rhs[r-1] = -direct[r, 0].subs(k, 1)  # P=z
        for t in range(1, r+1):
            matrix[r-1, t-1] = (-1)**t*comb(r, t)*direct[r, t].subs(k, 1)
    fixed_b = matrix.inv()*rhs
    assert matrix*fixed_b == rhs
    assert list(fixed_b) == [sp.Rational(2, 3), sp.Rational(61, 144),
                             sp.Rational(6203, 24000)]
    print("fourth-sector Chebyshev obstruction: 14 differentiated moments and the fixed b values checked",
          flush=True)
    u, v, w = fixed_b
    functions = [
        normalizers[4],
        kernels[1, 0]-u*normalizers[1],
        kernels[2, 0]-2*u*kernels[2, 1]+2*v*normalizers[2],
        kernels[3, 0]-3*u*kernels[3, 1]+3*v*kernels[3, 2]-6*w*normalizers[3],
        -kernels[4, 0]+4*u*kernels[4, 1]-6*v*kernels[4, 2]+4*w*kernels[4, 3]]
    functions = list(map(sp.cancel, functions))
    assert all(function.subs(k, 1) == 0 for function in functions[1:4])

    clearing = sp.prod((A+j)**(j+1) for j in range(5))
    polynomials, scales = [], []
    for function in functions:
        denominator, polynomial = sp.Poly(sp.cancel(clearing*function), k).clear_denoms(convert=True)
        content, primitive = polynomial.primitive()
        assert denominator > 0 and content > 0
        polynomials.append(primitive.as_expr())
        scales.append(sp.Rational(denominator, content))
    assert scales == [sp.Rational(1, 24), 3, 36, 6000, 1500]
    scale_product = sp.prod(scales)

    # If p_j=c_j*D*f_j, the Leibniz matrix has diagonal D. Hence
    # W(p)=(prod c_j)*D^5*W(f). Verify every entry of that identity
    # before using the polynomial determinant to determine its sign.
    raw_columns = []
    for function in functions:
        column = [function]
        for _ in range(4):
            column.append(sp.cancel(sp.diff(column[-1], k)))
        raw_columns.append(column)
    polynomial_derivatives = sp.Matrix([[sp.diff(polynomial, k, r) for polynomial in polynomials]
                                        for r in range(5)])
    leibniz = sp.Matrix(5, 5, lambda r, t: comb(r, t)*sp.diff(clearing, k, r-t)
                        if t <= r else 0)
    # Cancel after each derivative and work in QQ(k) for the products;
    # repeated expression expansion would obscure this small identity.
    field = sp.QQ.frac_field(k)
    for j, column in enumerate(raw_columns):
        raw = list(map(field.from_sympy, column))
        for r in range(5):
            expected = sum((field.from_sympy(leibniz[r, t])*raw[t]
                            for t in range(r+1)), field.zero)*field.from_sympy(scales[j])
            assert field.from_sympy(polynomial_derivatives[r, j]) == expected
    assert leibniz.det() == clearing**5
    wronskian = sp.Poly(polynomial_derivatives.det(method="domain-ge"), k)
    prefactor = 1152*(k+4)**3*(k+5)**12
    quotient, remainder = sp.div(wronskian, sp.Poly(prefactor, k), auto=False)
    assert remainder.is_zero and quotient.domain == sp.ZZ
    assert quotient.degree() == 30
    assert quotient.TC() == -65335420077170702484936960000
    assert quotient.LC() == 801221751
    assert wronskian.eval(0) < 0 and wronskian.LC() > 0
    print("fourth-sector Chebyshev obstruction: positive Wronskian scaling and R30 checked", flush=True)

    evaluations = []
    for nodes, expected_sign in ((tuple(range(5)), -1), (tuple(range(100, 105)), 1)):
        raw = sp.Matrix([[function.subs(k, node) for function in functions]
                         for node in nodes]).det()
        cleared = sp.Matrix([[polynomial.subs(k, node) for polynomial in polynomials]
                             for node in nodes]).det()
        determinant_scale = scale_product*sp.prod(clearing.subs(k, node) for node in nodes)
        assert determinant_scale > 0 and cleared == raw*determinant_scale
        assert sp.sign(raw) == expected_sign
        evaluations.append({"nodes": list(nodes), "determinant": str(raw),
                            "sign": expected_sign, "positive_clearing_identity_checked": True})

    print("fourth-sector Chebyshev obstruction: direct b solve, Wronskian scaling, and both determinant signs checked",
          flush=True)
    return {
        "status": "exact counterexample to this five-function Chebyshev shortcut",
        "parameters": {"a": 1, "x": 0, "P": "z"},
        "function_order": ["Z4", "f1", "f2", "f3", "g4"],
        "function_formulas": ["Z4", "J10-b1*Z1", "J20-2*b1*J21+2*b2*Z2",
                              "J30-3*b1*J31+3*b2*J32-6*b3*Z3",
                              "-J40+4*b1*J41-6*b2*J42+4*b3*J43"],
        "kernel_definition": "Jrt=integral_0^1 z^k*(1-z)^r*prod_(j=0)^(t-1)(z*d/dz+1+j)*(-z*log(z)/(1-z))^(r-t) dz; Zr=integral_0^1 z^k*(1-z)^r dz",
        "direct_integral_check": {
            "all_14_kernel_identities": True, "all_four_normalizers_and_diagonals": True,
            "matrix": [[str(entry) for entry in row] for row in matrix.tolist()],
            "rhs": list(map(str, rhs)), "fixed_b": list(map(str, fixed_b)),
            "constraints_vanish_at_k1": True},
        "positive_clearing": {"D4": str(clearing), "column_scales": list(map(str, scales)),
                              "column_scale_product": str(scale_product),
                              "primitive_columns": list(map(str, polynomials)),
                              "wronskian_identity": "W(c_j*D4*f_j)=(prod c_j)*D4^5*W(f_j)",
                              "all_25_leibniz_entries_checked": True},
        "wronskian": {"positive_prefactor": str(prefactor), "remaining_degree": 30,
                      "R30_coefficients_ascending": list(map(str, reversed(quotient.all_coeffs()))),
                      "constant": str(quotient.TC()), "leading": str(quotient.LC()),
                      "sign_at_zero": -1, "sign_for_sufficiently_large_k": 1},
        "ordered_evaluation_determinants": evaluations,
        "conclusion": "Opposite signs on ordered nodes rule out a T-system on k>=0 for this family",
        "limitations": "Does not disprove support<=3 by another argument, sector LC, or ARW negative association"}


def fourth_sector_four_support_4111_x0_certificate():
    """Reproduce only the ordered (4,1,1,1) coefficient at x=0.

    Domain: a>0 and i,d,e,f>=0, with support indices
    (i,i+d,i+d+e,i+d+e+f); multiplicity four is on the lowest support.
    This is not a certificate for other orientations, positive x, all
    four-support mixtures, or fourth-sector LC for arbitrary mixtures.

    The passed computation is preserved: exact factor tracking, full
    five-variable composition, all-coefficient scan/hash, and the direct
    differentiated-integral check at a=2, indices 0,1,2,3. No composition
    truncation, derivative-based replacement, or intermediate sign test is
    introduced. The removed factors are strictly positive on the domain.
    Strictness at zero gaps follows from 210 times the monomial numerator.

    This optional run previously took about 5 minutes and 1.4 GiB peak RSS.
    It uses optional python-flint and writes no files except through main.
    A temporary soft address-space cap of 1900 MiB retains any lower limit;
    the hard limit is unchanged and the original limits are restored finally.
    Mathematical body retained from the passed script with SHA256
    9191113d32218588645498029ed22f25b672de9b233ad0e882b06504948fb714.
    """
    import gc, hashlib, itertools, json, math, resource, time
    try:
        from flint import fmpz_mpoly_ctx
    except ImportError as error:
        raise SystemExit("--fourth-sector-four-support-4111-x0 requires optional python-flint") from error
    original_limits = resource.getrlimit(resource.RLIMIT_AS)
    soft_cap = 1900*1024**2
    for limit in original_limits:
        if limit != resource.RLIM_INFINITY:
            soft_cap = min(soft_cap, limit)
    resource.setrlimit(resource.RLIMIT_AS, (soft_cap, original_limits[1]))
    try:
        started = time.monotonic()
        ctx=fmpz_mpoly_ctx.get(['a','A0','A1','A2','A3'])
        a,*As=ctx.gens();zero,one=ctx.constant(0),ctx.constant(1)
        descriptors=[('a' if q==0 else f'a+{q}',a+q,0) for q in range(4)]
        for j,A in enumerate(As):descriptors.extend((f'A{j}+{q}',A+q,j+1) for q in range(5))
        empty=(0,)*len(descriptors);divisions=0;maximum=0
        def report(phase,**kw):
            d={'status':'in progress','phase':phase,'seconds':round(time.monotonic()-started,2),
               'peak_rss_mib':round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024,2),**kw}
            print(json.dumps(d),flush=True)
        def product(values):
            out=one
            for v in values:out=out*v
            return out
        def rising(v,n):return product(v+j for j in range(n))
        def falling(v,n):return product(v-j for j in range(n))
        class Factored:
            __slots__=('poly','valuation')
            def __init__(self,poly,valuation=empty,normalize=True):
                nonlocal divisions,maximum
                if not poly:self.poly,self.valuation=zero,empty;return
                maximum=max(maximum,len(poly));powers=list(valuation)
                if normalize:
                    for j,(_,factor,index) in enumerate(descriptors):
                        while poly.degrees()[index]>0:
                            quotient,remainder=divmod(poly,factor)
                            if remainder:break
                            assert quotient*factor==poly
                            poly=quotient;powers[j]+=1;divisions+=1
                self.poly,self.valuation=poly,tuple(powers)
            def __add__(self,other):
                if not self.poly:return other
                if not other.poly:return self
                common=tuple(min(p,q) for p,q in zip(self.valuation,other.valuation))
                def lift(obj):
                    out=obj.poly
                    for (_,factor,_),p,q in zip(descriptors,obj.valuation,common):
                        if p>q:out=out*factor**(p-q)
                    return out
                return Factored(lift(self)+lift(other),common)
            def __mul__(self,other):
                if not self.poly or not other.poly:return Fzero
                return Factored(self.poly*other.poly,tuple(p+q for p,q in zip(self.valuation,other.valuation)),False)
        Fzero,Fone=Factored(zero),Factored(one)
        def component(A):
            k=A-a;den=product((A+j)**(j+1) for j in range(5))
            def over(d):
                q=den/d;assert q*d==den;return q
            H,J={},{}
            for r in range(1,5):
                H[r]=Factored(rising(a+1,r-1)*math.factorial(r)*over(product(A+j for j in range(r+1))))
                for t in range(r):
                    m=r-t
                    J[r,t]=Factored(sum((math.factorial(m)*(-1)**(t+ell)*math.comb(t,ell)
                       *falling(k,t-ell)*rising(k+m+1,ell)*over((A+m+ell)**(m+1)) for ell in range(t+1)),zero))
            return H,J
        components=[component(A) for A in As]
        keys=((0,0,0),(1,0,0),(0,1,0),(0,0,1))
        H={r:{key:components[j][0][r] for j,key in enumerate(keys)} for r in range(1,5)}
        J={(r,t):{key:components[j][1][r,t] for j,key in enumerate(keys)} for r in range(1,5) for t in range(r)}
        def add_series(left,right):return {key:left.get(key,Fzero)+right.get(key,Fzero) for key in set(left)|set(right)}
        def mul_series(left,right):
            out={}
            for p,v in left.items():
                for q,w in right.items():
                    key=tuple(i+j for i,j in zip(p,q))
                    if max(key)<=1:out[key]=out.get(key,Fzero)+v*w
            return out
        N={0:{(0,0,0):Fone}}
        for r in range(1,5):
            total={}
            for t in range(r):
                term=mul_series(N[t],J[r,t])
                for h in range(t+1,r):term=mul_series(term,H[h])
                scalar=Factored((-1)**(r+t+1)*math.comb(r,t)*a**(r-t-1))
                total=add_series(total,{key:scalar*v for key,v in term.items()})
            N[r]=total;report('sector ready',r=r,coefficient_terms={str(k):len(v.poly) for k,v in total.items()})
        def selected(factors):
            out=Fzero
            def visit(pos,remaining,acc):
                nonlocal out
                if pos==len(factors):
                    if remaining==(0,0,0):out=out+acc
                    return
                for key,v in factors[pos].items():
                    if all(k<=r for k,r in zip(key,remaining)):visit(pos+1,tuple(r-k for k,r in zip(key,remaining)),acc*v)
            visit(0,(1,1,1),Fone);return out
        answer=selected([H[4],N[3],N[3]])+Factored(ctx.constant(-1))*selected([H[3],N[2],N[4]])
        valuation=answer.valuation;reduced=answer.poly
        factor_record=[{'factor':name,'power':p} for (name,_,_),p in zip(descriptors,valuation) if p]
        report('independent coefficient factored',terms=len(reduced),factors=factor_record)
        del answer,N,H,J,components,total,term;gc.collect()
        target=fmpz_mpoly_ctx.get(['a','i','d','e','f'])
        aa,ii,dd,ee,ff=target.gens()
        polynomial=reduced.compose(aa,aa+ii,aa+ii+dd,aa+ii+dd+ee,aa+ii+dd+ee+ff,ctx=target)
        del reduced;gc.collect()
        report('ordered coefficient composed',terms=len(polynomial),degrees=list(map(int,polynomial.degrees())))
        # Direct differentiated integrals at a=2, indices0,1,2,3, x=0.
        z,L=sp.symbols('z L');ps=sp.symbols('p0:4');v=z*L/(1-z)
        def Dz(expr):return sp.cancel(z*sp.diff(expr,z)-sp.diff(expr,L))
        def integral(expr):return sum(c*sp.factorial(q)/sp.Integer(p+1)**(q+1) for (p,q),c in sp.Poly(sp.cancel(expr),z,L).terms())
        ZH,JJ={},{}
        for r in range(1,5):
            ZH[r]=sp.Poly(sp.rf(3,r-1)*sum(ps[j]*integral(z**(j+1)*(1-z)**r) for j in range(4)),*ps)
            for t in range(r):
                kernel=v**(r-t)
                for h in range(t):kernel=Dz(kernel)+(2+h)*kernel
                JJ[r,t]=sp.Poly(sum(ps[j]*integral(z**(j+1)*(1-z)**r*kernel) for j in range(4)),*ps)
        NN={0:sp.Poly(1,*ps)}
        for r in range(1,5):
            NN[r]=sum(((-1)**(r+t+1)*math.comb(r,t)*2**(r-t-1)*NN[t]*JJ[r,t]
                   *math.prod(ZH[h] for h in range(t+1,r)) for t in range(r)),sp.Poly(0,*ps))
        gamma=ZH[4]*NN[3]**2-ZH[3]*NN[2]*NN[4]
        direct=gamma.coeff_monomial((4,1,1,1))
        scales=[math.prod((2+k+j)**(j+1) for j in range(5)) for k in range(4)]
        scale=scales[0]**4*math.prod(scales[1:])
        point=(2,0,1,1,1);independent=(2,2,3,4,5)
        factor_value=math.prod(int(factor(*independent))**p for (_,factor,_),p in zip(descriptors,valuation))
        assert int(polynomial(*point))*factor_value==direct*scale
        digest=hashlib.sha256();digest.update(b'[');last=None;negative=None;checked=0
        for offset,index in enumerate(range(len(polynomial)-1,-1,-1)):
            powers=tuple(map(int,polynomial.monomial(index)));co=int(polynomial.coefficient(index))
            assert last is None or last<powers;last=powers
            if co<0:negative={'powers':powers,'coefficient':co};break
            assert co>0
            if offset:digest.update(b',')
            digest.update(json.dumps([*powers,co],separators=(',',':')).encode());checked+=1
        digest.update(b']')
        out={'status':'positive' if negative is None else 'coefficient obstruction','scope':'Only fourth-sector pattern4111 at x=0',
             'domain':'a>0,i,d,e,f>=0','indices':['i','i+d','i+d+e','i+d+e+f'],'pattern':[4,1,1,1],
             'factors':factor_record,'exact_linear_divisions':divisions,'maximum_intermediate_terms':maximum,
             'terms':len(polynomial),'positive_terms_checked':checked,'negative':negative,
             'degrees':list(map(int,polynomial.degrees())),'sha256':digest.hexdigest() if negative is None else None,
             'direct_integral_audit':{'a':2,'indices':[0,1,2,3],'x':0,'coefficient':str(direct),'scale':str(scale),'removed_factor':str(factor_value),'matches':True},
             'peak_rss_mib':round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024,2),'seconds':round(time.monotonic()-started,2)}
        # Lock the completed certificate, without changing its mathematical body.
        assert out["status"] == "positive" and out["negative"] is None
        assert out["terms"] == out["positive_terms_checked"] == 4955029
        assert out["degrees"] == [63, 56, 38, 26, 13]
        assert out["exact_linear_divisions"] == 598
        assert out["maximum_intermediate_terms"] == 556741
        assert out["sha256"] == "06866d17e911643996c5ba835a3bdf12993ff2a68cb9a24b61cd76689b618cde"
        assert out["factors"] == [
            {"factor": "a+1", "power": 1}, {"factor": "a+2", "power": 1},
            {"factor": "A0+2", "power": 2}, {"factor": "A0+3", "power": 7},
            {"factor": "A0+4", "power": 15}]
        assert out["direct_integral_audit"]["coefficient"] == "2740516640336675465263/164701910459036160000000000"
        return out
    finally:
        resource.setrlimit(resource.RLIMIT_AS, original_limits)


def fourth_sector_lc_three_support_certificate():
    """All fifteen positive three-power fourth-sector weight coefficients.

    Domain: a>0 and i,d,e,x>=0, with powers (i,i+d,i+d+e).
    The preserved exact engine checks all thirteen x-slices separately,
    including zero slices, with positive-factor removal, ordered composition,
    streamed signs/hashes, and an independent formal-x integral audit.
    Together with the saved monomial and two-support results this proves
    strict b3^2>b2*b4 for at most three powers. Four supports are excluded.

    This is a long optional computation. Its mathematical body is retained
    from the completed standalone run; only local scoping and checkpoint
    output change. Original 511 metadata are embedded as an audited input,
    so reproduction requires no temporary files. Writes occur only in main.
    A temporary soft address-space cap of 1900 MiB retains any lower limit;
    the hard limit is unchanged, and the prior limits are restored finally.
    """
    import gc
    import hashlib
    import itertools
    import json
    import math
    import resource
    import time
    from types import SimpleNamespace
    import sympy as sp
    try:
        from flint import fmpz_mpoly_ctx
    except ImportError as error:
        raise SystemExit("--fourth-sector-lc-three-support requires optional python-flint") from error

    original_limits = resource.getrlimit(resource.RLIMIT_AS)
    soft_cap = 1900*1024**2
    for limit in original_limits:
        if limit != resource.RLIM_INFINITY:
            soft_cap = min(soft_cap, limit)
    resource.setrlimit(resource.RLIMIT_AS, (soft_cap, original_limits[1]))
    try:
        started = time.monotonic()
        args = SimpleNamespace(memory_mib=1900)
        old_511_path = "Original standalone 511 certificate, preserved hash and common factor"
        old_511 = {'status': 'positive reduced coefficient certificate', 'sha256': 'c6c75acbe7ec7baa47ab7711a746a64ea6fab078f87f990f08fde5b8340b7199', 'reduced_terms_after_composition': 3514838, 'removed_factor_list': [{'factor': 'a+1', 'power': 1}, {'factor': 'a+2', 'power': 1}, {'factor': 'A0+2', 'power': 4}, {'factor': 'A0+3', 'power': 9}, {'factor': 'A0+4', 'power': 16}]}
        assert old_511["status"] == "positive reduced coefficient certificate"
        assert old_511["sha256"] == 'c6c75acbe7ec7baa47ab7711a746a64ea6fab078f87f990f08fde5b8340b7199'

        representatives = ((5, 1, 1), (4, 2, 1), (3, 3, 1), (3, 2, 2))
        patterns = sorted({p for rep in representatives for p in itertools.permutations(rep)})
        assert set(patterns) == {(p, q, 7-p-q) for p in range(1, 6)
                                 for q in range(1, 6) if 7-p-q >= 1}
        assert len(patterns) == 15
        state = {
            'status': 'in progress',
            'scope': 'All 15 three-positive-weight degree-seven coefficients of the fourth-sector gap; no four-support theorem',
            'domain': 'a>0, i,d,e,x>=0',
            'indices': ['i', 'i+d', 'i+d+e'],
            'gap_numerator': 'Gamma=H4*N3^2-H3*N2*N4',
            'component_rescaling': 'M_j=D4(A_j)*F(A_j-a), D4(A)=A(A+1)^2(A+2)^3(A+3)^4(A+4)^5',
            'x_degree_coverage_required': list(range(13)),
            'series_x_degree_bounds': {'N1': 3, 'N2': 5, 'N3': 6, 'N4': 6, 'Gamma': 12},
            'x_handling': 'x remains an explicit nonnegative series exponent; no x factor is removed and no x=0 normalization is used',
            'factor_handling': 'Each slice stores all removed listed linear factors; common factor is their coordinatewise minimum over nonzero slices',
            'hash_format': 'Compact ASCII JSON of ascending [a,i,d,e,integer_coefficient] rows for the fully reduced x-slice',
            'original_511_preserved': {
                'path': str(old_511_path),
                'full_common-factor-reduced_hash': old_511['sha256'],
                'full_common-factor-reduced_terms': old_511['reduced_terms_after_composition'],
                'removed_factor_list': old_511['removed_factor_list'],
                'new_slice_hashes_may_use_additional_positive_factors': True,
            },
            'address_space_limit_mib': args.memory_mib,
            'patterns': {},
            'representative_runs': [],
        }


        def rss():
            return round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024, 2)


        def save(phase=None):
            if phase is not None:
                state['current_phase'] = phase
            state['seconds'] = round(time.monotonic()-started, 2)
            state['peak_rss_mib'] = rss()


        def report(phase, **details):
            save(phase)
            print(json.dumps({'phase': phase, 'seconds': state['seconds'],
                              'peak_rss_mib': rss(), **details}), flush=True)


        # Independent audit with x formal, a=1, and exponents 0,1,2.
        # This is a small rational polynomial at one fixed index/rate point.
        z, ell, sx, p0, p1, p2 = sp.symbols('z ell x p0 p1 p2')
        voltage = z*ell/(1-z)
        weights = (p0, p1, p2)
        poly_gens = (p0, p1, p2, sx)


        def Dz(expression):
            return sp.cancel(z*sp.diff(expression, z)-sp.diff(expression, ell))


        def integral(expression):
            poly = sp.Poly(sp.cancel(expression), z, ell)
            assert all(j >= 0 and h >= 0 for (j, h), value in poly.terms())
            return sum(value*math.factorial(h)/sp.Integer(j+1)**(h+1)
                       for (j, h), value in poly.terms())


        direct_Z, direct_J = {}, {}
        for r in range(1, 5):
            direct_Z[r] = sp.Poly(sum(
                weights[j]*sum(math.comb(4-r, h)*sx**h
                               *integral(z**(j+h)*(1-z)**r)
                               for h in range(5-r)) for j in range(3)), *poly_gens)
            for t in range(r):
                kernel = voltage**(r-t)
                for h in range(t):
                    kernel = sp.cancel(Dz(kernel)+(1+h)*kernel)
                direct_J[r, t] = sp.Poly(sum(
                    weights[j]*sum(math.comb(4-r, h)*sx**h
                                   *integral(z**(j+h)*(1-z)**r*kernel)
                                   for h in range(5-r)) for j in range(3)), *poly_gens)
        direct_H = {r: math.factorial(r)*direct_Z[r] for r in range(1, 5)}
        direct_N = {0: sp.Poly(1, *poly_gens)}
        for r in range(1, 5):
            value = sp.Poly(0, *poly_gens)
            for t in range(r):
                term = (-1)**(r+t+1)*math.comb(r, t)*direct_N[t]*direct_J[r, t]
                for h in range(t+1, r):
                    term *= direct_H[h]
                value += term
            direct_N[r] = value
        direct_gamma = direct_H[4]*direct_N[3]**2-direct_H[3]*direct_N[2]*direct_N[4]
        assert direct_gamma.degree(sx) <= 12
        assert all(sum(powers[:3]) == 7 for powers, value in direct_gamma.terms())
        direct_coefficients = {p: [direct_gamma.coeff_monomial((*p, h)) for h in range(13)]
                               for p in patterns}
        scales = [math.prod((A+j)**(j+1) for j in range(5)) for A in (1, 2, 3)]
        state['direct_audit'] = {
            'method': 'Differential kernels integrated term by term at a=1, indices 0,1,2, with x kept formal',
            'point': {'a': 1, 'i': 0, 'd': 1, 'e': 1},
            'component_denominators': list(map(str, scales)),
            'original_weight_and_x_coefficients': {
                ''.join(map(str, p)): list(map(str, values)) for p, values in direct_coefficients.items()},
        }
        del direct_Z, direct_J, direct_H, direct_N, direct_gamma
        gc.collect()
        report('formal-x direct-integral audit table ready')

        # Native symbolic arithmetic has only four variables. x is never one of
        # the polynomial generators; it is a separate series coordinate.
        ctx = fmpz_mpoly_ctx.get(['a', 'A0', 'A1', 'A2'])
        a, A0, A1, A2 = ctx.gens()
        zero, one = ctx.constant(0), ctx.constant(1)
        descriptors = [('a' if q == 0 else f'a+{q}', a+q, 0)
                       for q in range(4)]
        for j, variable in enumerate((A0, A1, A2)):
            descriptors.extend((f'A{j}' if q == 0 else f'A{j}+{q}', variable+q, j+1)
                               for q in range(5))
        empty = (0,)*len(descriptors)
        division_count = 0
        max_tracked_polynomial_terms = 0


        def product(values):
            out = one
            for value in values:
                out = out*value
            return out


        def rising(value, length):
            return product(value+j for j in range(length))


        def falling(value, length):
            return product(value-j for j in range(length))


        class Factored:
            __slots__ = ('poly', 'valuation')

            def __init__(self, poly, valuation=empty, normalize=True):
                nonlocal division_count, max_tracked_polynomial_terms
                if not poly:
                    self.poly, self.valuation = zero, empty
                    return
                max_tracked_polynomial_terms = max(max_tracked_polynomial_terms, len(poly))
                powers = list(valuation)
                if normalize:
                    for j, (_, factor, index) in enumerate(descriptors):
                        while poly.degrees()[index] > 0:
                            quotient, remainder = divmod(poly, factor)
                            if remainder:
                                break
                            assert quotient*factor == poly
                            poly = quotient
                            powers[j] += 1
                            division_count += 1
                self.poly, self.valuation = poly, tuple(powers)

            def __add__(self, other):
                if not self.poly:
                    return other
                if not other.poly:
                    return self
                common = tuple(min(p, q) for p, q in zip(self.valuation, other.valuation))

                def lift(value):
                    poly = value.poly
                    for (_, factor, _), exponent, base in zip(descriptors, value.valuation, common):
                        if exponent > base:
                            poly = poly*factor**(exponent-base)
                    return poly

                return Factored(lift(self)+lift(other), common)

            def __mul__(self, other):
                if not self.poly or not other.poly:
                    return Fzero
                return Factored(self.poly*other.poly,
                                tuple(p+q for p, q in zip(self.valuation, other.valuation)), False)


        Fzero, Fone = Factored(zero), Factored(one)


        def component(A):
            index = A-a
            denominator = product((A+j)**(j+1) for j in range(5))

            def over(value):
                quotient = denominator/value
                assert quotient*value == denominator
                return quotient

            H, J = {}, {}
            for r in range(1, 5):
                for h in range(5-r):
                    Z = (math.comb(4-r, h)*math.factorial(r)
                         *over(product(A+h+j for j in range(r+1))))
                    H[r, h] = Factored(rising(a+1, r-1)*Z)
                    for t in range(r):
                        m = r-t
                        J[r, t, h] = Factored(sum((
                            math.comb(4-r, h)*math.factorial(m)
                            *(-1)**(t+ell)*math.comb(t, ell)
                            *falling(index+h, t-ell)*rising(index+h+m+1, ell)
                            *over((A+h+m+ell)**(m+1)) for ell in range(t+1)), zero))
            return H, J


        components = [component(A) for A in (A0, A1, A2)]
        weight_keys = ((0, 0), (1, 0), (0, 1))
        H = {r: {(p, q, h): components[j][0][r, h]
                 for j, (p, q) in enumerate(weight_keys) for h in range(5-r)}
             for r in range(1, 5)}
        J = {(r, t): {(p, q, h): components[j][1][r, t, h]
                      for j, (p, q) in enumerate(weight_keys) for h in range(5-r)}
             for r in range(1, 5) for t in range(r)}
        del components
        report('per-x factored component vectors ready', exact_linear_divisions=division_count)

        target = fmpz_mpoly_ctx.get(['a', 'i', 'd', 'e'])
        aa, ii, dd, ee = target.gens()
        ordered_A = (aa+ii, aa+ii+dd, aa+ii+dd+ee)


        def pattern_assignment(rep, pattern):
            available = list(range(3))
            out = []
            for multiplicity in rep:
                chosen = next(j for j in available if pattern[j] == multiplicity)
                out.append(chosen)
                available.remove(chosen)
            return tuple(out)


        def relabel_valuation(valuation, assignment):
            out = list(valuation[:4])+[0]*15
            for j in range(3):
                for q in range(5):
                    out[4+5*assignment[j]+q] = valuation[4+5*j+q]
            return tuple(out)


        def factor_list(valuation):
            return [{'factor': name, 'power': exponent}
                    for (name, _, _), exponent in zip(descriptors, valuation) if exponent]


        def add_series(left, right):
            return {key: left.get(key, Fzero)+right.get(key, Fzero)
                    for key in set(left)|set(right)}


        def mul_series(left, right, caps):
            out = {}
            for (i, j, h), first in left.items():
                for (p, q, t), second in right.items():
                    if i+p <= caps[0] and j+q <= caps[1]:
                        key = (i+p, j+q, h+t)
                        out[key] = out.get(key, Fzero)+first*second
            return out


        def selected_product(factors, wanted):
            out = Fzero
            first, second, third = factors
            for key1, value1 in first.items():
                for key2, value2 in second.items():
                    remaining = tuple(wanted[j]-key1[j]-key2[j] for j in range(3))
                    if min(remaining) < 0:
                        continue
                    value3 = third.get(remaining)
                    if value3 is not None:
                        out = out+value1*value2*value3
            return out


        def check_slice(pattern, assignment, x_degree, selected):
            label = ''.join(map(str, pattern))
            local_valuation = relabel_valuation(selected.valuation, assignment)
            original = direct_coefficients[pattern][x_degree]
            scale = math.prod(scales[j]**pattern[j] for j in range(3))
            expected_raw = original*scale
            assert expected_raw.q == 1
            record = {
                'x_degree': x_degree,
                'zero_slice': not bool(selected.poly),
                'removed_factor_list': factor_list(local_valuation),
                'factor_valuation': list(local_valuation),
                'independent_reduced_terms': len(selected.poly),
                'direct_original_weight_and_x_coefficient': str(original),
                'positive_component_denominator_scale': str(scale),
                'raw_cleared_value_at_audit_point': str(expected_raw),
            }
            if not selected.poly:
                assert expected_raw == 0
                record.update({'terms': 0, 'degrees': [0]*4, 'positive_terms_checked': 0,
                               'first_negative': None, 'sha256': hashlib.sha256(b'[]').hexdigest(),
                               'direct_integral_reconstruction_passed': True})
                return record

            save(f'composing pattern {label}, x-degree {x_degree}')
            answer = selected.poly.compose(aa, *(ordered_A[j] for j in assignment), ctx=target)
            record['terms'] = len(answer)
            record['degrees'] = list(map(int, answer.degrees()))
            audit_values = (1, 1, 2, 3)
            factor_value = math.prod(int(factor(*audit_values))**exponent
                                    for (_, factor, _), exponent in zip(descriptors, local_valuation))
            assert answer(1, 0, 1, 1)*factor_value == int(expected_raw)
            record['removed_factor_at_audit_point'] = str(factor_value)
            record['direct_integral_reconstruction_passed'] = True

            digest = hashlib.sha256()
            digest.update(b'[')
            previous, minimum, negative, checked = None, None, None, 0
            for offset, index in enumerate(range(len(answer)-1, -1, -1)):
                powers = tuple(map(int, answer.monomial(index)))
                value = int(answer.coefficient(index))
                if previous is not None:
                    assert previous < powers
                previous = powers
                if value < 0:
                    negative = {'powers': list(powers), 'coefficient': value}
                    break
                assert value > 0
                minimum = value if minimum is None else min(minimum, value)
                if offset:
                    digest.update(b',')
                digest.update(json.dumps([*powers, value], separators=(',', ':')).encode('ascii'))
                checked += 1
            digest.update(b']')
            record.update({'positive_terms_checked': checked, 'first_negative': negative,
                           'minimum_positive_coefficient': minimum,
                           'sha256': digest.hexdigest() if negative is None else None})
            del answer
            gc.collect()
            return record


        def finish_pattern(pattern):
            label = ''.join(map(str, pattern))
            record = state['patterns'][label]
            assert set(record['slices']) == set(map(str, range(13)))
            nonzero = [s for s in record['slices'].values() if not s['zero_slice']]
            assert nonzero
            common = tuple(min(s['factor_valuation'][j] for s in nonzero)
                           for j in range(len(descriptors)))
            record['common_removed_factor_list'] = factor_list(common)
            for s in nonzero:
                additional = tuple(v-c for v, c in zip(s['factor_valuation'], common))
                assert min(additional) >= 0
                s['additional_slice_factor_beyond_common'] = factor_list(additional)
            record['zero_x_degrees'] = [h for h in range(13) if record['slices'][str(h)]['zero_slice']]
            record['total_reduced_slice_coefficients'] = sum(s['terms'] for s in record['slices'].values())
            record['all_x_degrees_0_through_12_checked'] = True
            record['all_coefficients_nonnegative'] = all(s['first_negative'] is None for s in record['slices'].values())
            # Strictness of the full selected weight coefficient follows from x=0:
            # at i=d=e=0 all components coincide, giving 7!/prod(m_j!) times
            # the strictly positive monomial Gamma. Positive factors remain positive.
            assert not record['slices']['0']['zero_slice']
            record['strictness'] = 'x^0 slice: coincident components give multinomial(7;pattern) times the strict monomial Gamma; nonnegative index coefficients extend strictness'
            if pattern == (5, 1, 1):
                assert record['common_removed_factor_list'] == old_511['removed_factor_list']
                record['common_factor_matches_original_511'] = True


        try:
            for rep in representatives:
                rep_label = ''.join(map(str, rep))
                caps = rep[1:]
                orientations = sorted(set(itertools.permutations(rep)))
                for pattern in orientations:
                    label = ''.join(map(str, pattern))
                    state['patterns'][label] = {
                        'pattern': list(pattern), 'canonical_representative': list(rep),
                        'canonical_to_ordered_component': list(pattern_assignment(rep, pattern)),
                        'slices': {},
                    }
                report('building representative sector series', representative=list(rep), caps=list(caps))
                N = {0: {(0, 0, 0): Fone}}
                for r in range(1, 5):
                    total = {}
                    for t in range(r):
                        term = mul_series(N[t], J[r, t], caps)
                        for h in range(t+1, r):
                            term = mul_series(term, H[h], caps)
                        scalar = Factored((-1)**(r+t+1)*math.comb(r, t)*a**(r-t-1))
                        total = add_series(total, {key: scalar*value for key, value in term.items()})
                    assert all(key[2] <= r*(7-r)//2 and key[0]+key[1] <= r for key in total)
                    N[r] = total
                    report('representative sector row ready', representative=list(rep), r=r,
                           series_entries=len(total), max_reduced_terms=max(len(v.poly) for v in total.values()))

                slice_independent = []
                for x_degree in range(13):
                    save(f'building representative {rep_label}, x-degree {x_degree}')
                    wanted = (rep[1], rep[2], x_degree)
                    selected = (selected_product([H[4], N[3], N[3]], wanted)
                                +Factored(ctx.constant(-1))*selected_product([H[3], N[2], N[4]], wanted))
                    if selected.poly:
                        for _, factor, _ in descriptors:
                            assert divmod(selected.poly, factor)[1]
                    slice_independent.append({'x_degree': x_degree,
                                              'reduced_terms': len(selected.poly),
                                              'factor_valuation': list(selected.valuation)})
                    counts = {}
                    for pattern in orientations:
                        label = ''.join(map(str, pattern))
                        assignment = pattern_assignment(rep, pattern)
                        result = check_slice(pattern, assignment, x_degree, selected)
                        state['patterns'][label]['slices'][str(x_degree)] = result
                        save()
                        counts[label] = result['terms']
                        if result['first_negative'] is not None:
                            state['status'] = 'stopped at first negative reduced coefficient; no inequality counterexample inferred'
                            state['first_negative'] = {'pattern': list(pattern), 'x_degree': x_degree,
                                                       **result['first_negative']}
                            save('negative coefficient found')
                            print(json.dumps(state['first_negative']), flush=True)
                            raise SystemExit(2)
                    report('representative x-slice complete', representative=list(rep), x_degree=x_degree,
                           independent_reduced_terms=len(selected.poly), ordered_term_counts=counts)
                    del selected
                    gc.collect()

                for pattern in orientations:
                    finish_pattern(pattern)
                state['representative_runs'].append({'representative': list(rep),
                                                     'orientations': list(map(list, orientations)),
                                                     'independent_slices': slice_independent})
                report('representative complete', representative=list(rep),
                       patterns_completed=sum(bool(r.get('all_coefficients_nonnegative'))
                                              for r in state['patterns'].values()))
                del N, total, term
                gc.collect()

        except MemoryError:
            state['status'] = 'stopped at resource error'
            save('MemoryError under the 1900 MiB maximum address-space ceiling')
            raise

        assert set(state['patterns']) == {''.join(map(str, p)) for p in patterns}
        assert all(r['all_coefficients_nonnegative'] for r in state['patterns'].values())
        state['status'] = 'all 15 three-support coefficients proved positive'
        state['all_15_complete'] = True
        state['all_195_x_slices_covered'] = True
        state['total_reduced_slice_coefficients'] = sum(r['total_reduced_slice_coefficients']
                                                       for r in state['patterns'].values())
        state['exact_successful_linear_divisions'] = division_count
        state['maximum_tracked_independent_polynomial_terms'] = max_tracked_polynomial_terms
        state['four_support_status'] = 'not attempted'
        save('complete')
        print(json.dumps({'status': state['status'],
                          'total_reduced_slice_coefficients': state['total_reduced_slice_coefficients'],
                          'peak_rss_mib': state['peak_rss_mib'],
                          'seconds': state['seconds']}, indent=2), flush=True)

        # Lock every recorded slice, factor and hash to the completed run.
        assert hashlib.sha256(json.dumps(state["patterns"], sort_keys=True,
                                         separators=(",", ":")).encode()).hexdigest() == 'af2e5c92820b246af34e1e3cf61b0f54b110af9696eaf63441c8c504b6ac2110'
        state["canonical_integration"] = {'standalone_source_sha256': 'de31e7c2edb7557226c81b3b74a68280de2847d262973bae763c00e0299aa496', 'mathematical_body_ast_sha256': '71fa3e592b649cdb4e06b4180aa25be65573b59968715ff43c296241546d74d9', 'completed_patterns_sha256': 'af2e5c92820b246af34e1e3cf61b0f54b110af9696eaf63441c8c504b6ac2110', 'body_preservation': 'AST identical after removing standalone checkpoint writes and converting global counters to nonlocal', 'normalization_check': 'Separate one-point differentiated-integral audit; full coefficient run not repeated'}
        state.update({'mixture_consequence': 'Strict b3^2>b2*b4 for a>0, x>=0 and every nonzero nonnegative mixture of at most three nonnegative real powers', 'lower_support_dependencies': 'Saved fourth_sector_lc.monomial and fourth_sector_lc.arbitrary_two_support; neither is rerun by this flag', 'limitations': 'No four-support, arbitrary-finite-mixture or all-rank ULC claim'})
        state["address_space_limit_mib"] = soft_cap/1024**2
        return state
    finally:
        resource.setrlimit(resource.RLIMIT_AS, original_limits)


def main():
    import argparse
    import json
    from pathlib import Path

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--excursion-symbolic", action="store_true",
                        help="Only the all-rate, all-return-count excursion-tail certificate")
    parser.add_argument("--full-event-sector", action="store_true",
                        help="Only four exact checks of the general full-occupation reduction")
    parser.add_argument("--full-event-highrate", action="store_true",
                        help="Only the symbolic two-site and weighted three-site source checks")
    parser.add_argument("--one-vacancy-block", action="store_true",
                        help="Only the coefficient trace/scaling and forest-sign audit")
    parser.add_argument("--full-event-triangle", action="store_true",
                        help="All-rate weighted-triangle source certificate (requires python-flint)")
    parser.add_argument("--full-event-tree", action="store_true",
                        help="All-rate branching-tree source certificate (requires python-flint)")
    parser.add_argument("--full-event-equal-row", action="store_true",
                        help="Only the fixed three-site equal-row integral/source normalization check")
    parser.add_argument("--equal-row-vacancies", action="store_true",
                        help="Only the saved four-site vacancy check and saved homogeneous top-ULC normalization")
    parser.add_argument("--third-sector-lc", action="store_true",
                        help="Only monomial and arbitrary-two-support third-sector LC certificates (requires python-flint)")
    parser.add_argument("--third-sector-lc-all", action="store_true",
                        help="Prior third-sector certificate plus all six three-support coefficients (requires python-flint)")
    parser.add_argument("--fourth-sector-lc", action="store_true",
                        help="Only fourth-sector monomial and two-power boundary certificates (requires python-flint)")
    parser.add_argument("--fourth-sector-lc-two-support", action="store_true",
                        help="All six arbitrary-two-support fourth-sector coefficients; several minutes (requires python-flint)")
    parser.add_argument("--fourth-sector-lc-three-support", action="store_true",
                        help="All fifteen three-support fourth-sector coefficients in thirteen x-slices; long exact run (requires python-flint)")
    parser.add_argument("--fourth-sector-chebyshev-obstruction", action="store_true",
                        help="Only the fixed five-function Chebyshev obstruction and direct integral audit")
    parser.add_argument("--fourth-sector-four-support-4111-x0", action="store_true",
                        help="Only ordered 4111 at x=0; full composition, several minutes and about 1.4 GiB (requires python-flint)")
    parser.add_argument("--output", type=Path, help="Optional JSON certificate; merge existing keys")
    args = parser.parse_args()
    if args.full_event_tree:
        result = {"full_event_tree": full_event_tree_certificate()}
    elif args.full_event_triangle:
        result = {"full_event_triangle": full_event_triangle_certificate()}
    elif args.one_vacancy_block:
        result = {"one_vacancy_block": one_vacancy_block_certificate()}
    elif args.full_event_highrate:
        result = {"full_event_highrate": full_event_highrate_certificate()}
    elif args.full_event_sector:
        result = {"full_event_sector": full_event_sector_certificate()}
    elif args.excursion_symbolic:
        result = {"excursion_tails_symbolic": excursion_tail_symbolic_certificate()}
    elif args.full_event_equal_row:
        result = {"full_event_equal_row": full_event_equal_row_certificate()}
    elif args.equal_row_vacancies:
        result = {"equal_row_vacancies": equal_row_vacancies_certificate()}
    elif args.third_sector_lc_all:
        result = {"third_sector_lc": third_sector_lc_certificate(),
                  "third_sector_lc_all": third_sector_lc_all_certificate()}
    elif args.third_sector_lc:
        result = {"third_sector_lc": third_sector_lc_certificate()}
    elif args.fourth_sector_lc:
        result = {"fourth_sector_lc": fourth_sector_lc_certificate()}
    elif args.fourth_sector_lc_two_support:
        # Retain the audited base without invoking its flag. Prefer values
        # already in the selected output while retaining canonical fields
        # absent from a partial output. Nothing is written to the source.
        source = Path(__file__).with_name("algebra_certificates.json")
        base = json.loads(source.read_text()).get("fourth_sector_lc", {}) if source.exists() else {}
        if args.output and args.output.exists():
            existing = json.loads(args.output.read_text()).get("fourth_sector_lc", {})
            assert isinstance(existing, dict), "fourth_sector_lc must be an object"
            base = {**base, **existing}
        assert isinstance(base, dict), "fourth_sector_lc must be an object"
        monomial = base.get("monomial", {})
        if (monomial.get("numerator_sha256") != "3a7653b6cbc5c73f9c91d7c7c0a224ecf1b2d166a1e8cf5f4a6251da7895dc60"
                or not monomial.get("positive_coefficient_certificate")):
            parser.error("The two-support run requires the saved --fourth-sector-lc monomial certificate")
        result = {"fourth_sector_lc": {
            **base, "arbitrary_two_support": fourth_sector_lc_two_support_certificate()}}
    elif args.fourth_sector_lc_three_support:
        # Preserve every existing nested certificate; no earlier checker runs.
        source = Path(__file__).with_name("algebra_certificates.json")
        base = json.loads(source.read_text()).get("fourth_sector_lc", {}) if source.exists() else {}
        assert isinstance(base, dict), "fourth_sector_lc must be an object"
        if args.output and args.output.exists():
            existing = json.loads(args.output.read_text()).get("fourth_sector_lc", {})
            assert isinstance(existing, dict), "fourth_sector_lc must be an object"
            base = {**base, **existing}
        monomial = base.get("monomial", {})
        two_support = base.get("arbitrary_two_support", {})
        if (monomial.get("numerator_sha256") != "3a7653b6cbc5c73f9c91d7c7c0a224ecf1b2d166a1e8cf5f4a6251da7895dc60"
                or not monomial.get("positive_coefficient_certificate")
                or not two_support.get("all_six_mixed_complete_and_positive")):
            parser.error("The three-support run requires the saved monomial and arbitrary-two-support certificates")
        result = {"fourth_sector_lc": {
            **base, "arbitrary_three_support": fourth_sector_lc_three_support_certificate()}}
    elif args.fourth_sector_chebyshev_obstruction:
        # Preserve every saved fourth-sector child without invoking its
        # checker. The selected output can also be a partial certificate.
        source = Path(__file__).with_name("algebra_certificates.json")
        base = json.loads(source.read_text()).get("fourth_sector_lc", {}) if source.exists() else {}
        assert isinstance(base, dict), "fourth_sector_lc must be an object"
        if args.output and args.output.exists():
            existing = json.loads(args.output.read_text()).get("fourth_sector_lc", {})
            assert isinstance(existing, dict), "fourth_sector_lc must be an object"
            base = {**base, **existing}
        result = {"fourth_sector_lc": {
            **base, "chebyshev_obstruction": fourth_sector_chebyshev_obstruction_certificate()}}
    elif args.fourth_sector_four_support_4111_x0:
        # Retain every prior child, including the obstruction, without
        # rerunning any other mathematical certificate.
        source = Path(__file__).with_name("algebra_certificates.json")
        base = json.loads(source.read_text()).get("fourth_sector_lc", {}) if source.exists() else {}
        assert isinstance(base, dict), "fourth_sector_lc must be an object"
        if args.output and args.output.exists():
            existing = json.loads(args.output.read_text()).get("fourth_sector_lc", {})
            assert isinstance(existing, dict), "fourth_sector_lc must be an object"
            base = {**base, **existing}
        result = {"fourth_sector_lc": {
            **base, "four_support_4111_x0": fourth_sector_four_support_4111_x0_certificate()}}
    else:
        result = diagnostics()
    if args.output:
        # Both modes preserve certificates produced only by optional runs.
        saved = json.loads(args.output.read_text()) if args.output.exists() else {}
        if "fourth_sector_lc" in result:
            # A base-only rerun must also retain a completed optional
            # arbitrary_two_support certificate nested under this key.
            previous = saved.get("fourth_sector_lc", {})
            assert isinstance(previous, dict), "fourth_sector_lc must be an object"
            result["fourth_sector_lc"] = {**previous, **result["fourth_sector_lc"]}
        saved.update(result)
        args.output.write_text(json.dumps(saved, indent=2) + "\n")


if __name__ == "__main__":
    main()
