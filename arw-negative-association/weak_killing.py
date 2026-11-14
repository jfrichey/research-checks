#!/usr/bin/env python3
"""Exact leading vacancy coefficients for K_epsilon=(1-epsilon)P.

P is closed and irreducible; sleep rates are fixed and positive. This is
neither weak internal jumps nor a large-sleep expansion. The raw calculation
uses closed ARW Green operators followed by successive single-loss operators.
The independent IDLA calculation rescales each vacancy by epsilon.

Purposeful checks: --audit uses one rational three-site path; --path3-symbolic
proves all leading cylinder gaps for rates (r,h,r); --cycle4-raw reproduces the
previous four-cycle calculation in opposite legal priorities. No simulations.
--uniform-renewal audits the positive uniform-destination loss product at n=3.
"""

import argparse
from collections import defaultdict
import json
from pathlib import Path

from exact import DomainMatrix, absorption, sp


def particle_count(state):
    return sum(1 if x == -1 else x for x in state)


def stable_mask(state):
    return sum(1 << i for i, x in enumerate(state) if x == -1)


def active(state):
    return any(x > 0 for x in state)


def closed_rows(state, transition, rates, reverse=False):
    """Closed transition row and the coefficient of one killed jump."""
    order = range(len(state)-1, -1, -1) if reverse else range(len(state))
    i = next(i for i in order if state[i] > 0)
    zero = defaultdict(lambda: sp.S.Zero)
    jump = 1/(1+rates[i]) if state[i] == 1 else sp.S.One
    if state[i] == 1:
        sleep = list(state)
        sleep[i] = -1
        zero[tuple(sleep)] += rates[i]/(1+rates[i])
    for j, probability in enumerate(transition[i]):
        if probability == 0:
            continue
        nxt = list(state)
        nxt[i] -= 1
        nxt[j] = 2 if nxt[j] == -1 else nxt[j]+1
        zero[tuple(nxt)] += jump*probability
    lost = list(state)
    lost[i] -= 1
    assert sp.cancel(sum(zero.values())-1) == 0
    return dict(zero), (tuple(lost), jump, i)


def solve(left, right):
    aa = DomainMatrix.from_Matrix(left).to_field()
    bb = DomainMatrix.from_Matrix(right).convert_to(aa.domain)
    answer = aa.lu_solve(bb)
    assert (aa.matmul(answer)-bb).is_zero_matrix
    return answer.to_Matrix().applyfunc(sp.factor)


def leading_coefficients(transition, rates, reverse=False, marked=False):
    """Return c[A] for vacancy mask A, with optional first-loss post-states.

    A row with m particles and terminal occupied set B stores the coefficient
    of epsilon**(m-|B|). At the leading order, internal-jump factors 1-epsilon
    can be replaced by one. Stable inputs cannot undergo further losses.
    """
    transition = tuple(tuple(map(sp.sympify, row)) for row in transition)
    rates = tuple(map(sp.sympify, rates))
    n = len(transition)
    start = (1,)*n
    states, discovered, transitions = [start], {start}, {}
    for state in states:
        zero, loss = closed_rows(state, transition, rates, reverse)
        transitions[state] = zero, loss
        for nxt in list(zero)+[loss[0]]:
            if active(nxt) and nxt not in discovered:
                discovered.add(nxt)
                states.append(nxt)

    solutions, blocks = {}, {}

    def terminal(state):
        if active(state):
            return solutions[state]
        result = sp.zeros(1, 1 << n)
        result[0, stable_mask(state)] = 1
        return result

    for number in range(1, n+1):
        block = [s for s in states if particle_count(s) == number]
        index = {s:i for i, s in enumerate(block)}
        left, right = sp.eye(len(block)), sp.zeros(len(block), 1 << n)
        for state, i in index.items():
            zero, (lost, weight, _) = transitions[state]
            for nxt, probability in zero.items():
                if active(nxt):
                    left[i, index[nxt]] -= probability
                else:
                    right[i, stable_mask(nxt)] += probability
            right[i, :] += weight*terminal(lost)
        result = solve(left, right)
        for state, i in index.items():
            solutions[state] = result[i, :]
            assert all(result[i, a] == 0 for a in range(1 << n)
                       if a.bit_count() > number)
            assert sp.cancel(sum(result[i, a] for a in range(1 << n)
                                 if a.bit_count() == number)-1) == 0
        blocks[number] = (block, left)

    full = (1 << n)-1
    coefficients = tuple(solutions[start][0, full ^ a] for a in range(1 << n))
    assert coefficients[0] == 1
    out = {"coefficients": coefficients,
           "block_sizes": {m:len(b[0]) for m, b in blocks.items()}}
    if marked:
        block, left = blocks[n]
        initial = sp.zeros(len(block), 1)
        initial[block.index(start), 0] = 1
        visits = solve(left.T, initial)
        weights = defaultdict(lambda: sp.S.Zero)
        departures = [sp.S.Zero]*n
        for k, state in enumerate(block):
            lost, weight, vertex = transitions[state][1]
            weights[lost] += visits[k, 0]*weight
            departures[vertex] += visits[k, 0]*weight
        weights = {s:sp.factor(w) for s, w in weights.items() if w != 0}
        for a in range(1, 1 << n):
            assert sp.cancel(sum(w*terminal(s)[0, full ^ a]
                                 for s, w in weights.items())-coefficients[a]) == 0
        out["marked"] = [(s, w, terminal(s)) for s, w in weights.items()]
        out["departure_mass"] = tuple(map(sp.factor, departures))
    return out


def closed_potentials(transition):
    """Proper-set determinants, exit probabilities, and mean exit times."""
    p = sp.Matrix(transition)
    n, result = p.rows, {}
    for mask in range((1 << n)-1):
        inside = [i for i in range(n) if mask >> i & 1]
        outside = [i for i in range(n) if not mask >> i & 1]
        laplacian = sp.eye(len(inside))-p.extract(inside, inside)
        inverse = laplacian.inv()
        h = inverse*sp.ones(len(inside), 1)
        exits = inverse*p.extract(inside, outside)
        result[mask] = (inside, outside, laplacian.det(), h, exits)
    return p, result


def rescaled_idla(transition, inverse_rates):
    """Construct Mhat without symbolic epsilon; columns index occupied sets.

    Ahat has mean-exit-time diagonal and closed exit probabilities off the
    diagonal. On proper columns Jhat_i=I+sum_j P_ij Ahat_j. Full columns are
    set to zero: they are never used in <=n operations from the empty state.
    These truncated operators must not be assumed to commute globally.
    """
    p, potential = closed_potentials(transition)
    n, size = p.rows, 1 << p.rows
    additions = [sp.zeros(size) for _ in range(n)]
    for mask, (inside, outside, _, h, exits) in potential.items():
        for i in range(n):
            if i not in inside:
                additions[i][mask | (1 << i), mask] = 1
            else:
                row = inside.index(i)
                additions[i][mask, mask] = h[row]
                for k, j in enumerate(outside):
                    additions[i][mask | (1 << j), mask] = exits[row, k]
    operators = []
    for i in range(n):
        jump = sp.eye(size)
        for j in range(n):
            jump += p[i, j]*additions[j]
        jump[:, -1] = sp.zeros(size, 1)
        operators.append((1+inverse_rates[i])*additions[i]-inverse_rates[i]*jump)
    empty = sp.zeros(size, 1)
    empty[0] = 1
    columns = [empty]
    for mask in range(1, size):
        i = (mask & -mask).bit_length()-1
        columns.append(operators[i]*columns[mask ^ (1 << i)])
    matrix = sp.Matrix.hstack(*columns)
    return matrix, operators


def one_vacancy_source(transition, inverse_rates):
    """u_i=-Mhat[V minus {i},V], and the first full-normalizer derivative.

    This algebraic source is not the chronological first-loss departure law.
    Its subsequent closed (n-1)-particle stabilization gives c_single.

    The returned derivative is also the expected total number M of jumps
    in closed saturated stabilization, including self-loop jumps. With
    L=I-P and kappa=sum_i det(L[V minus {i}]),
      M=kappa*sum_{A nonempty} x_A/det(L[A complement]).
    Thus M is multiaffine in the inverse sleep rates. The expected number
    of jumps from vertex i is M*pi_i, where pi is stationary for P.
    """
    p, potential = closed_potentials(transition)
    n, full = p.rows, (1 << p.rows)-1
    x = tuple(map(sp.sympify, inverse_rates))
    source = []
    for i in range(n):
        b = full ^ (1 << i)
        value = sp.S.Zero
        a = b
        while True:
            remaining = b ^ a
            inside, _, determinant, h, _ = potential[remaining]
            return_time = 1+sum(p[i, j]*h[k] for k, j in enumerate(inside))
            value += sp.prod(x[j] for j in range(n) if a >> j & 1)*return_time/determinant
            if a == 0:
                break
            a = (a-1) & b
        source.append(sp.factor(x[i]*potential[b][2]*value))
    cofactor_sum = sum(potential[full ^ (1 << i)][2] for i in range(n))
    normalizer_derivative = cofactor_sum*sum(
        sp.prod(x[i] for i in range(n) if a >> i & 1)/potential[full ^ a][2]
        for a in range(1, full+1))
    assert sp.cancel(sum(source)-normalizer_derivative) == 0
    return sp.Matrix(source), sp.factor(normalizer_derivative)


def positive_coefficients(value, variables):
    numerator, denominator = sp.fraction(sp.factor(value))
    return (all(c > 0 for c in sp.Poly(numerator, *variables).coeffs())
            and all(c > 0 for c in sp.Poly(denominator, *variables).coeffs()))


def path3_record(coefficients, r, h, block_sizes):
    """Audit already-computed coefficients without solving again."""
    c = tuple(coefficients)
    assert c[0] == 1 and c[1] == c[4] and c[3] == c[6]
    types = {"leaf_center": (1, 2), "leaf_leaf": (1, 4),
             "leaf_pair": (1, 6), "center_pair": (2, 5)}
    gaps = {}
    for name, (a, b) in types.items():
        gap = sp.factor(c[a]*c[b]-c[a | b])
        assert positive_coefficients(gap, (r, h))
        gaps[name] = str(gap)
    checked = 0
    for a in range(1, 8):
        for b in range(a+1, 8):
            if not a & b:
                assert positive_coefficients(c[a]*c[b]-c[a | b], (r, h))
                checked += 1
    assert checked == 6
    # All three-site NA types: singleton vs singleton, AND, or OR.
    # The OR leading gap is the sum of the two singleton-pair gaps.
    for i in range(3):
        others = [j for j in range(3) if j != i]
        value = sum(c[1 << i]*c[1 << j]-c[(1 << i) | (1 << j)] for j in others)
        assert positive_coefficients(value, (r, h))
    return {"rates": ["r", "h", "r"], "parameters": "r,h>0",
            "method": "raw closed Green/loss hierarchy; exact residuals checked",
            "coefficients_by_vacancy_mask": list(map(str, c)),
            "distinct_cylinder_gaps": gaps, "positive_coefficients": True,
            "disjoint_cylinder_pairs_checked": checked,
            "NA_event_types": {"order_two": 6, "order_three": 3},
            "conclusion": "NA for sufficiently small killing, for every fixed r,h>0",
            "block_sizes": block_sizes}


def path3_symbolic():
    r, h = sp.symbols("r h", positive=True)
    p = ((0, 1, 0), (sp.Rational(1, 2), 0, sp.Rational(1, 2)), (0, 1, 0))
    result = leading_coefficients(p, (r, h, r))
    return path3_record(result["coefficients"], r, h, result["block_sizes"])


def rescaled_audit():
    p = ((0, 1, 0), (sp.Rational(1, 2), 0, sp.Rational(1, 2)), (0, 1, 0))
    matrix, operators = rescaled_idla(p, (1, 1, 1))
    target = sp.zeros(8, 1)
    target[7] = 1
    leading = solve(matrix, target)
    expected = tuple(sp.Rational(x) for x in
                     ("1", "33/5", "24/5", "169/10", "33/5", "27", "169/10", "219/5"))
    assert tuple(leading[7 ^ a] for a in range(8)) == expected
    # Reverse the factor order only in admissible products from empty.
    empty = sp.zeros(8, 1)
    empty[0] = 1
    for mask in range(8):
        alternate = empty
        for i in range(3):
            if mask >> i & 1:
                alternate = operators[i]*alternate
        assert alternate == matrix[:, mask]
    source, derivative = one_vacancy_source(p, (1, 1, 1))
    holes = [6, 5, 3]
    assert source == -matrix.extract(holes, [7]) == sp.Matrix([6, 6, 6])
    kernel = tuple(tuple((j, value) for j, value in enumerate(row) if value) for row in p)
    starts = [tuple((mask >> i) & 1 for i in range(3)) for mask in holes]
    raw, count = absorption(starts, kernel, (sp.Rational(1, 2),)*3)
    closed = raw.extract(range(3), holes).T
    assert matrix.extract(holes, holes)*closed == sp.eye(3)
    assert closed*source == sp.Matrix([expected[1], expected[2], expected[4]])
    assert derivative == sum(expected[1 << i] for i in range(3)) == 18

    # One rational marked-loss check, also auditing the failed local bound.
    marked = leading_coefficients(p, (1, 1, 1), marked=True)
    assert marked["coefficients"] == expected
    assert marked["departure_mass"] == (sp.Rational(9, 2), sp.S(9), sp.Rational(9, 2))
    state, weight, row = next(item for item in marked["marked"] if item[0] == (0, 2, 0))
    local_gap = sp.factor((expected[1]*row[0, 3]+expected[4]*row[0, 6])/2-row[0, 2])
    assert weight == 1 and row[0, 6] == row[0, 3] == sp.Rational(2, 5)
    assert row[0, 2] == sp.Rational(8, 3) and local_gap == -sp.Rational(2, 75)

    # Parent-supplied raw kernel: no duplicate raw solve for these rates.
    nonreversible = sp.Matrix([[130, 34, 22], [33, 105, 27], [69, 93, 183]])/232
    assert all(sum(nonreversible[:, j]) == 1 for j in range(3))
    cycle_gap = (nonreversible[0, 1]*nonreversible[1, 2]*nonreversible[2, 0]
                 -nonreversible[1, 0]*nonreversible[2, 1]*nonreversible[0, 2])
    assert cycle_gap == -sp.Rational(9, 26912)
    return {"case": "closed path3, rates (1,1,1)",
            "rescaled_coefficients_match_raw": True,
            "admissible_product_orders_agree": True,
            "one_vacancy_source": list(map(str, source)),
            "closed_source_factorization_checked": True,
            "closed_two_particle_transient_states": count,
            "normalizer_derivative": str(derivative),
            "normalizer_derivative_interpretation": "expected total closed saturated jumps, including self-loops",
            "first_loss_departure_mass": list(map(str, marked["departure_mass"])),
            "departure_mass_interpretation": "also the expected closed saturated jumps from each vertex",
            "failed_local_bound": {"state": state, "first_loss_weight": str(weight),
                                   "u_0": "2/5", "u_2": "2/5", "v_02": "8/3",
                                   "gap": str(local_gap)},
            "closed_nonreversibility": {
                "rates": [1, 2, 3], "occupied_mask_order": [3, 5, 6],
                "kernel": [[str(x) for x in row] for row in nonreversible.tolist()],
                "cycle_gap": str(cycle_gap),
                "validation": "parent raw kernel in both priorities; column sums and cycle arithmetic checked here"}}


def cycle4_raw():
    half = sp.Rational(1, 2)
    p = ((0, half, 0, half), (half, 0, half, 0),
         (0, half, 0, half), (half, 0, half, 0))
    forward = leading_coefficients(p, (1,)*4)
    reverse = leading_coefficients(p, (1,)*4, reverse=True)
    c = forward["coefficients"]
    assert c == reverse["coefficients"]
    for a in range(1, 16):
        expected = {1:sp.Rational(61, 6), 3:sp.Rational(1079, 6), 4:sp.S(440)}
        if a.bit_count() == 2:
            value = sp.S(57) if a in (5, 10) else sp.Rational(205, 4)
        else:
            value = expected[a.bit_count()]
        assert c[a] == value
    gaps = [c[a]*c[b]-c[a | b] for a in range(1, 16)
            for b in range(a+1, 16) if not a & b]
    assert len(gaps) == 25 and min(gaps) > 0
    return {"rates": [1]*4, "opposite_priorities_agree": True,
            "coefficient_classes": {"singleton": "61/6", "neighbor_pair": "205/4",
                                    "opposite_pair": "57", "triple": "1079/6", "full": "440"},
            "strict_disjoint_cylinder_gaps": len(gaps),
            "conclusion": "NA for sufficiently small uniform killing at fixed sleep rate one",
            "NA_reduction": "strict leading disjoint-cylinder gaps and minimum-size event witnesses; no further event enumeration",
            "block_sizes": forward["block_sizes"],
            "matches_previous_independent_IDLA_algebra": True}


def uniform_loss_kernel(n, h, inverse_rate):
    """Return the positive polynomial matrix B_h for uniform closed ARW.

    There are m=n-h+1 particles. Rows index input sleeper counts 0,...,m;
    columns index counts 0,...,m-1 after one loss. The stable last row is
    zero. After scaling out n*x per loss, the one-loss operator is B_h/h.
    Here x=inverse_rate is the fixed holding-model inverse sleep rate.

    The nonzero block is S*C, with S[s,v]=1[v>=s] and
    C[s,j]=(1+x)**(s-j)*binom(h+j-1,j)/binom(h+s,s) for j<=s.
    These factors also give total nonnegativity for x>=0.
    """
    if n < 1 or not 1 <= h <= n:
        raise ValueError("require n>=1 and 1<=h<=n")
    m = n-h+1
    t = 1+sp.sympify(inverse_rate)
    lower = sp.Matrix(m, m, lambda s, j:
                      t**(s-j)*sp.binomial(h+j-1, j)/sp.binomial(h+s, s)
                      if j <= s else 0)
    suffix = sp.Matrix(m, m, lambda s, v: 1 if v >= s else 0)
    core = (suffix*lower).applyfunc(sp.expand)
    return core.col_join(sp.zeros(1, m))


def uniform_loss_polynomials(n, inverse_rate):
    """Return d_0,...,d_n and B_1,...,B_n for the actual uniform family.

    The specified-hole coefficient is c_A=(n*x)**len(A)*d[len(A)].
    Exactly d_r=e_0.T*B_1*...*B_r*1 / falling_factorial(n,r).
    Every d_r has positive coefficients through degree r*(2*n-r-1)/2.
    This positivity does not by itself prove all-rank log-concavity.
    """
    if n < 1:
        raise ValueError("require n>=1")
    row = sp.zeros(1, n+1)
    row[0, 0] = 1
    coefficients, kernels = [sp.S.One], []
    for r in range(1, n+1):
        kernel = uniform_loss_kernel(n, r, inverse_rate)
        kernels.append(kernel)
        row = (row*kernel).applyfunc(sp.expand)
        coefficients.append(sp.expand(sum(row)/sp.ff(n, r)))
    return coefficients, kernels


def uniform_renewal_certificate():
    """Audit one symbolic n=3 product against the saved independent values."""
    x = sp.symbols("x", nonnegative=True)
    n = 3
    d, kernels = uniform_loss_polynomials(n, x)
    for h, kernel in enumerate(kernels, 1):
        m = n-h+1
        thinning = sp.Matrix(m, m, lambda s, j:
                             sp.binomial(h+j-1, j)/sp.binomial(h+s, s)
                             if j <= s else 0)
        cumulative = sp.Matrix(m, m, lambda s, k:
                              sp.binomial(h+k, k)/sp.binomial(h+s, s)
                              if k < s else 0)
        lower = sp.Matrix(m, m, lambda s, j:
                          (1+x)**(s-j)*sp.binomial(h+j-1, j)/sp.binomial(h+s, s)
                          if j <= s else 0)
        shift = sp.Matrix(m, m, lambda s, j: 1 if j == s+1 else 0)
        core = kernel[:m, :]
        assert thinning*sp.ones(m, 1) == sp.ones(m, 1)
        assert ((sp.eye(m)-x*cumulative)*lower-thinning).applyfunc(sp.expand) == sp.zeros(m)
        assert (((1+x)*sp.eye(m)-shift-x*thinning)*core-thinning).applyfunc(sp.expand) == sp.zeros(m)
        assert kernel[m, :] == sp.zeros(1, m)
        for s in range(m):
            for j in range(m):
                marked = sp.Rational(n, h)*thinning[s, j]
                preceding = sp.Rational(n, h)*thinning[s-1, j] if s else 0
                assert sp.cancel((h+s)*marked-s*preceding-n*int(s == j)) == 0

    expected = [sp.S.One, 1+x/2+x*x/9,
                sp.Rational(3, 4)+16*x/27+7*x*x/36+x**3/54,
                sp.Rational(49, 108)+23*x/54+17*x*x/108+x**3/54]
    assert all(sp.expand(u-v) == 0 for u, v in zip(d, expected))
    degrees, leading, coefficient_lower_checks = [], [], []
    for r in range(1, n+1):
        degree = r*(2*n-r-1)//2
        lc = 1/(sp.ff(n, r)*sp.prod(sp.binomial(n, h) for h in range(1, r+1)))
        polynomial = sp.Poly(d[r], x)
        assert polynomial.degree() == degree and polynomial.LC() == lc
        assert all(polynomial.nth(k) > 0 for k in range(degree+1))
        remainder = sp.Poly(sp.expand(d[r]-lc*(1+x)**degree), x)
        assert all(c >= 0 for c in remainder.all_coeffs())
        degrees.append(degree)
        leading.append(str(lc))
        coefficient_lower_checks.append(True)
    gap_degrees = []
    for s in range(1, n):
        gap = sp.Poly(sp.expand(d[s]**2-d[s-1]*d[s+1]), x)
        assert gap.degree() == 2*degrees[s-1]
        assert gap.LC() == sp.sympify(leading[s-1])**2
        gap_degrees.append(gap.degree())
    # Independent chronological-loss values already audited at n=3,x=1.
    raw_saved = [1, sp.Rational(29, 6), 14, sp.Rational(57, 2)]
    assert all(sp.cancel(n**r*d[r].subs(x, 1)-raw_saved[r]) == 0
               for r in range(n+1))
    return {
        "scope": "actual saturated uniform closed ARW; fixed holding inverse sleep x",
        "selected_n": n,
        "d_from_positive_product": [str(v) for v in d],
        "all_stage_thinning_normalizations_pass": True,
        "all_stage_marked_loss_recursions_pass": True,
        "all_stage_triangular_and_macro_chain_residuals_pass": True,
        "degrees": degrees,
        "leading_coefficients": leading,
        "coefficientwise_lower_bounds_pass": coefficient_lower_checks,
        "gap_degrees": gap_degrees,
        "saved_raw_x1_values_matched_not_rerun": [str(v) for v in raw_saved],
        "matrix_dimensions": [[k.rows, k.cols] for k in kernels],
        "limitations": "TN of the individual kernels does not prove all-x LC across loss levels",
    }


def uniform_rank_quotient_certificate():
    """Audit the positive quotient F_n/(1+y/n), report (35ao), at n=3.

    The general proof is the last-column identity and the rank-one
    dimension identity displayed in the report. Interlacing remains open
    beyond the small symbolic case checked here.
    """
    x, y, t = sp.symbols("x y t")
    n = 3
    row = sp.zeros(1, n+1)
    row[0, 0] = 1
    previous_total, previous_stable = sp.S.One, sp.S.Zero
    rank, quotient = sp.S.One, sp.S.One
    rows = []
    for h in range(1, n+1):
        kernel = uniform_loss_kernel(n, h, x)
        assert kernel[:, -1] == sp.Matrix([sp.Rational(h, n)]*(n-h+1)+[0])
        row = (row*kernel).applyfunc(sp.expand)
        total, stable_mass = sum(row), row[-1]
        assert sp.expand(stable_mass-sp.Rational(h, n)*(previous_total-previous_stable)) == 0
        rank += total*y**h/sp.factorial(h)
        quotient += (total-stable_mass)*y**h/sp.factorial(h)
        if h < n:
            addition = sp.ones(n-h+1, 1)*sp.Matrix(1, n-h, lambda _, j:
                sp.binomial(h+j-1, j)*(1+x)**(n-h-j)/sp.binomial(n, h))
            assert (kernel[:-1, :-1]-uniform_loss_kernel(n-1, h, x)
                    -addition).applyfunc(sp.expand) == sp.zeros(n-h+1, n-h)
        rows.append({"losses": h, "total": str(sp.expand(total)),
                     "stable_mass": str(stable_mass), "active_mass": str(sp.expand(total-stable_mass))})
        previous_total, previous_stable = total, stable_mass
    assert sp.expand(rank-(1+y/n)*quotient) == 0
    assert all(c > 0 for c in sp.Poly(quotient, x, y).coeffs())
    quotient_t = sp.factor(quotient.subs(x, t-1))
    first = sp.factor(quotient_t.subs(y, -2))
    second = sp.factor(quotient_t.subs(y, -2/(t+2)))
    assert first == t*(t+1)*(2*t+3)/9
    assert second == -t*(t+1)*(4*t+3)/(9*(t+2)**2)
    assert quotient_t.subs(y, 0) == 1
    rank2 = (1+y/2)*(1+(t+2)*y/2)
    assert sp.expand(quotient_t-rank2-y*t*(t+1)*(12+(2*t+9)*y)/36) == 0
    # Independent level updates do not preserve real roots. Here only the
    # second size-two kernel receives its actual size-three augmentation.
    first_kernel = uniform_loss_kernel(2, 1, 0)
    second_kernel = uniform_loss_kernel(2, 2, 0)/2
    second_kernel += sp.ones(2, 1)/6
    start = sp.Matrix([[1, 0, 0]])
    mixed = 1+(start*first_kernel*sp.ones(2, 1))[0]*y
    mixed += (start*first_kernel*second_kernel)[0]*y*y
    assert mixed == 1+2*y+sp.Rational(13, 12)*y*y
    assert sp.discriminant(mixed, y) == -sp.Rational(1, 3)
    partial = 1+4*y+sp.Rational(21, 4)*y*y+sp.Rational(125, 54)*y**3+sp.Rational(2137, 6912)*y**4
    residuals = [sp.factor(sum(sp.binomial(4-j, k)*k**j/sp.factorial(j)
                    *sp.diff(partial, y, j).subs(y, -k) for j in range(5-k)))
                 for k in range(1, 5)]
    assert residuals == [0, 3, 0, 0]
    partial_quotient = sp.cancel(partial/(1+y/4))
    assert sp.discriminant(partial_quotient, y) == -sp.Rational(9991, 6912)
    return {
        "scope": "General positive quotient and active-kernel identity; symbolic n=3 audit",
        "last_column": "B_h(s,n-h)=h/n for s<=n-h, zero in the stable input row",
        "stable_mass_recurrence": "z_h=(h/n)*(Z_(h-1)-z_(h-1)), z_0=0",
        "positive_quotient": "F_n(y)=(1+y/n)*G_n(y), G_n=sum_(h=0)^(n-1)(Z_h-z_h)*y^h/h!",
        "active_kernel_identity": "Btilde_h^(n)=B_h^(n-1)+1*[binom(h+j-1,j)*t^(n-h-j)/binom(n,h)]_j",
        "selected_n": n, "rows": rows, "rank_polynomial": str(sp.factor(rank)),
        "quotient_polynomial": str(sp.expand(quotient)),
        "last_column_recurrence_and_rankone_identity_pass": True,
        "G3_at_minus2": str(first), "G3_at_minus2_over_tplus2": str(second),
        "small_interlacing": "For every formal t>0, the roots of G_3 and F_2 alternate, and F_2 strictly interlaces F_3",
        "stagewise_obstruction": {"n": 3, "t": 1, "polynomial": str(mixed),
                                  "discriminant": "-1/3", "actual_family_counterexample": False},
        "partial_interpolation_obstruction": {"n": 4, "t": 1, "polynomial": str(partial),
            "row_residuals": list(map(str, residuals)), "quotient_discriminant": "-9991/6912",
            "actual_family_counterexample": False},
        "all_size_interlacing_or_LC_proved": False}


def uniform_parking_jensen_certificate():
    """One exact n=16 Jensen discriminator at x=0, and one prefix obstruction.

    This checks the new entire-function route without a size/rate sweep.
    The triangular parking recurrence is compared with the independent
    positive loss-kernel product at the SAME size. It is finite evidence,
    not a proof that all Jensen polynomials have negative real roots.
    """
    import hashlib
    n=16;d=[sp.S.One]
    for r in range(1,n+1):
        coefficient=-sum((-1)**(r-j)*d[r-j]
            /(sp.factorial(r-j)*sp.factorial(j)*sp.Integer(r)**j)
            for j in range(1,r+1))
        d.append((-1)**r*sp.factorial(r)*coefficient)
    assert all(value>0 for value in d)
    independent,_=uniform_loss_polynomials(n,0)
    assert d==independent
    y=sp.Symbol('y')
    polynomial=sp.Poly(sum(sp.binomial(n,r)*d[r]*y**r for r in range(n+1)),y)
    assert polynomial.eval(-n)==0
    _,integer=polynomial.clear_denoms();integer=integer.primitive()[1]
    gcd=sp.gcd(integer,integer.diff())
    negative=int(integer.count_roots(-sp.oo,0));positive=int(integer.count_roots(0,sp.oo))
    assert negative==n and positive==0 and gcd.degree()==0
    rows=[int(integer.nth(r)) for r in range(n+1)]
    prefix=sp.Poly(1,y)
    for j in (1,2,3):
        prefix=sp.Poly(sum(prefix.nth(k)*sum((sp.Rational(1,j*(j+1))*y)**ell/sp.factorial(ell)
                         for ell in range(j-k+1))*y**k for k in range(j+1)),y)
    assert prefix.as_expr()==1+sp.Rational(3,4)*y+sp.Rational(5,32)*y*y+sp.Rational(109,10368)*y**3
    discriminant=sp.discriminant(prefix.as_expr(),y)
    assert discriminant==sp.Rational(-149,1990656)<0
    return {'status':'finite n=16 Jensen check passed; all-size real-rootedness OPEN',
        'degree':n,'physical_inverse_rate':0,'negative_real_roots':negative,
        'positive_real_roots':positive,'nonreal_roots':n-negative-positive,
        'gcd_with_derivative_degree':int(gcd.degree()),'fixed_root_minus_n':True,
        'independent_positive_loss_product_matches_all_d_values':True,
        'primitive_integer_coefficients_ascending':rows,
        'coefficient_sha256':hashlib.sha256(json.dumps(rows,separators=(',',':')).encode()).hexdigest(),
        'method':'Exact rational Sturm root count; gcd certifies simplicity',
        'all_size_real_rootedness_proved':False,
        'prefix_approximation_obstruction':{'degree':3,'polynomial':str(prefix.as_expr()),
            'discriminant':str(discriminant),'actual_Jensen_or_D_counterexample':False},
        'root_size_or_rate_sweep':False}


def uniform_parking_spectrum_certificate():
    """Bounded checks for the physical x=0 compact parking operator.

    The analytic operator proof identifies the nearest zero as uniquely
    negative and simple. These exact checks fix the first three chain
    coefficients, the regularized determinant convention, and the stated
    quadratic/Gram obstructions. They do not prove reality of other zeros.
    No root count, size/rate grid, or earlier certificate helper is run.
    """
    import hashlib
    import time

    started = time.monotonic()
    r = sp.Symbol("r", positive=True)
    L, v, z = sp.symbols("L v z", real=True)
    # K^n 1=r^n W_n(v-log r), where V_n'=W_n and V_n(0)=0.
    # The exponential-polynomial integral is a finite derivative sum.
    polynomials = [L]
    for n in range(1, 4):
        previous = polynomials[-1]
        derivative = sum(sp.diff(previous, L, j)/sp.Integer(n)**(j+1)
                         for j in range(sp.degree(previous, L)+1))
        current = sp.integrate(derivative, L)
        current = sp.expand(current-current.subs(L, 0))
        assert sp.expand(n*sp.diff(current, L)-sp.diff(current, L, 2)-previous) == 0
        assert current.subs(L, 0) == 0
        polynomials.append(current)
    W = [sp.diff(poly, L) for poly in polynomials]
    functions = [r**n*W[n].subs(L, v-sp.log(r)) for n in range(4)]
    for n in range(1, 4):
        assert sp.simplify(sp.diff(functions[n], r)
                           -r**(n-1)*polynomials[n-1].subs(L, v-sp.log(r))) == 0
        assert sp.limit(functions[n], r, 0, dir="+") == 0
    assert W[:3] == [1, L+1, L**2/4+3*L/4+sp.Rational(3, 8)]
    assert sp.expand(W[3]-(L**3/36+11*L**2/72+49*L/216+sp.Rational(49, 648))) == 0

    def triangle_integral(poly):
        # integral_Omega r^j L^k = k!/(j+1)^(k+2), L=v-log r.
        terms = sp.Poly(sp.expand(poly.subs(v, L+sp.log(r))), r, L).terms()
        return sum(coefficient*sp.factorial(k)/sp.Integer(j+1)**(k+2)
                   for (j, k), coefficient in terms)

    area = triangle_integral(sp.S.One)
    chain = [triangle_integral(functions[n]) for n in range(3)]
    assert area == 1 and chain == [1, sp.Rational(3, 8), sp.Rational(49, 648)]
    assert chain == [W[n].subs(L, 0) for n in range(1, 4)]
    # Independent triangular normalization, restricted to these three degrees.
    gamma = [sp.S.One]
    for n in range(1, 4):
        gamma.append(-sum(gamma[ell]/(sp.factorial(n-ell)*sp.Integer(n)**(n-ell))
                          for ell in range(n)))
    assert chain == [(-1)**n*gamma[n] for n in range(1, 4)]
    a1, a2, a3 = chain
    norm_K_squared, norm_T_squared = a2, area**2-a2
    assert (norm_K_squared, norm_T_squared) == (sp.Rational(3, 8), sp.Rational(5, 8))

    # No ordinary trace of T is used: T is Hilbert--Schmidt, not trace class.
    p2, p3 = 1-2*a2, 1-3*a2+3*a3
    assert (p2, p3) == (sp.Rational(1, 4), sp.Rational(11, 108))
    determinant_series = sp.exp(z-p2*z**2/2+p3*z**3/3).series(z, 0, 4).removeO().expand()
    assert determinant_series == 1+a1*z+a2*z**2+a3*z**3
    # The product anomaly exponent is lambda*m+lambda*(1-m)=lambda.
    lam, scalar_trace = sp.symbols("lambda scalar_trace")
    assert sp.expand(lam*scalar_trace+lam*(1-scalar_trace)-lam) == 0

    # One fixed rectangle witness; both the pairing and direct-antiderivative
    # constructions are bounded exact checks of the same function.
    u, y = sp.symbols("u y", real=True)
    trigonometric = [sp.cos(2*sp.pi*u), sp.sin(2*sp.pi*u)]
    antiderivatives = [sp.sin(2*sp.pi*u)/(2*sp.pi),
                       (1-sp.cos(2*sp.pi*u))/(2*sp.pi)]
    for function, primitive in zip(trigonometric, antiderivatives):
        assert sp.simplify(sp.diff(primitive, u)-function) == 0
        assert primitive.subs(u, 0) == 0
    pairings = sp.Matrix(2, 2, lambda i, j:
                         sp.simplify(sp.integrate(trigonometric[i]*antiderivatives[j], (u, 0, 1))))
    assert pairings == sp.Matrix([[0, -1/(4*sp.pi)], [1/(4*sp.pi), 0]])
    quadratic_unit = sp.simplify(sum(value**2 for value in pairings))
    assert quadratic_unit == 1/(8*sp.pi**2)
    g = sp.cos(2*sp.pi*(u-y))
    tensor_primitive = (g-sp.cos(2*sp.pi*u)-sp.cos(2*sp.pi*y)+1)/(4*sp.pi**2)
    assert sp.simplify(sp.diff(tensor_primitive, u, y)-g) == 0
    assert sp.simplify(tensor_primitive.subs(u, 0)) == 0
    assert sp.simplify(tensor_primitive.subs(y, 0)) == 0
    assert sp.simplify(sp.integrate(g, (u, 0, 1))) == 0
    assert sp.simplify(sp.integrate(g, (y, 0, 1))) == 0
    assert sp.integrate(sp.cos(2*sp.pi*u)**2, (u, 0, 1))/(4*sp.pi**2) == quadratic_unit
    rectangle_area = sp.log(2)/4
    quadratic = sp.simplify(-rectangle_area**2*quadratic_unit)
    assert quadratic == -sp.log(2)**2/(128*sp.pi**2)

    # The original integral observable would force this two-vector Gram
    # matrix under a positive symmetrizing form on the cyclic space.
    beta = [sp.S.One, 1-a2, 1-2*a2+a3]
    assert beta == [1, sp.Rational(5, 8), sp.Rational(211, 648)]
    observable_gram = sp.Matrix([[beta[0], beta[1]], [beta[1], beta[2]]])
    gram_determinant = sp.factor(observable_gram.det())
    assert gram_determinant == -sp.Rational(337, 5184)
    f_to_three = 1-a1*z+a2*z**2-a3*z**3
    reciprocal_transfer = ((1/f_to_three-1)/z).series(z, 0, 3).removeO().expand()
    assert reciprocal_transfer == sum(beta[j]*z**j for j in range(3))
    actual_norm_squared = sp.simplify(triangle_integral((a2-functions[1])**2))
    assert actual_norm_squared == sp.Rational(359, 5184) > 0

    # One strict three-point poset pattern, represented by (log r,v).
    # This avoids numerical exponential comparisons and holds on open sets.
    points = [(sp.Integer(-4), sp.Integer(-3)),
              (sp.Integer(-2), sp.Integer(-1)),
              (sp.Integer(-5), -sp.Rational(1, 2))]
    assert all(log_r < height < 0 for log_r, height in points)
    sampled_K = sp.Matrix(3, 3, lambda i, j:
                          1 if points[j][0] < points[i][0] and points[j][1] < points[i][1] else 0)
    assert sampled_K == sp.Matrix([[0, 0, 0], [1, 0, 0], [0, 0, 0]])
    sampled_T = sp.ones(3)-sampled_K
    pair_minors = [sampled_T.extract(pair, pair).det()
                   for pair in ((0, 1), (0, 2), (1, 2))]
    assert pair_minors == [1, 0, 0] and sampled_T.det() == 0
    gram_u, gram_v, conjugate_u, conjugate_v = sp.symbols("u v ubar vbar")
    formal_gram = sp.Matrix([[1, 0, gram_u], [0, 1, gram_v],
                            [conjugate_u, conjugate_v, 1]])
    phase_identity = sp.expand(formal_gram.det())
    assert phase_identity == 1-gram_u*conjugate_u-gram_v*conjugate_v
    forced_determinant = phase_identity.subs(
        {gram_u*conjugate_u: 1, gram_v*conjugate_v: 1}, simultaneous=True)
    assert forced_determinant == -1

    polynomial_rows = [[str(sp.Poly(poly, L).nth(j)) for j in range(sp.degree(poly, L)+1)]
                       for poly in W]
    return {
        "status": "PASS: bounded exact normalization and obstruction checks",
        "physical_inverse_rate": 0,
        "sequence_convention": "a_n=d_n(0)/n!; D(z)=sum a_n*z^n; f(lambda)=D(-lambda)",
        "analytical_theorem": {
            "proved_scope": "D has exactly one zero of minimum modulus; it is negative and simple",
            "nearest_zero": "z_*=-1/r(T)",
            "proof_input": "The report's all-degree chain normalization, regularized determinant identity, and strict-square-kernel Perron argument; these are distinct from the finite checks here",
            "remaining_zero_reality": "OPEN",
            "all_zero_reality_proved": False,
            "all_size_Jensen_real_rootedness_proved": False,
            "positive_inverse_rate_extension_proved": False,
            "stationary_ULC_or_NA_proved_here": False},
        "operator": {
            "space": "L2(Omega,dr dv), Omega={0<r<exp(v),v<0}",
            "area": str(area),
            "K_kernel": "1_(s<r,w<v)", "P": "1 integral_Omega", "T": "P-K",
            "hilbert_schmidt_norm_squares": {"K": str(norm_K_squared), "T": str(norm_T_squared)},
            "ordinary_trace_or_determinant_of_T_used": False},
        "normalization": {
            "degrees_checked": [1, 2, 3],
            "V_polynomials": [str(poly) for poly in polynomials],
            "W_polynomials": [str(poly) for poly in W],
            "polynomial_recurrence": "n V_nprime-V_ndoubleprime=V_(n-1), V_n(0)=0",
            "polynomial_recurrence_checked": True,
            "original_coordinate_antiderivatives_and_zero_limits_checked": True,
            "triangle_moment_formula": "integral_Omega r^j*(v-log r)^k=k!/(j+1)^(k+2)",
            "chain_integrals_a1_a2_a3": list(map(str, chain)),
            "corner_values_and_chain_integrals_agree": True,
            "independent_triangular_parking_recurrence_checked": True,
            "W_coefficient_sha256": hashlib.sha256(json.dumps(polynomial_rows, separators=(",", ":")).encode()).hexdigest()},
        "regularized_determinant": {
            "identity": "D(z)=exp(z)*det_2(I+z*T); f(lambda)=exp(-lambda)*det_2(I-lambda*T)",
            "trace_T_squared": str(p2), "trace_T_cubed": str(p3),
            "series_through_degree_3": str(determinant_series),
            "log_series_convention": "z+sum_(k>=2)(-1)^(k+1)*tr(T^k)*z^k/k",
            "rank_one_product_anomaly_exponent_checked": True},
        "rectangle_quadratic_obstruction": {
            "rectangle": "(1/4,1/2) times (-log(2),0)",
            "rescaled_test_function": "cos(2*pi*(u-y))",
            "rectangle_area": str(rectangle_area),
            "volterra_cos_sin_pairings": [[str(value) for value in row] for row in pairings.tolist()],
            "direct_tensor_antiderivative_and_zero_mean_checked": True,
            "quadratic_value": str(quadratic),
            "scope": "Negative full-space quadratic form; does not automatically restrict to the cyclic space"},
        "cyclic_observable_metric_obstruction": {
            "beta_0_1_2": list(map(str, beta)),
            "required_Gram_matrix": [[str(value) for value in row] for row in observable_gram.tolist()],
            "Hankel_determinant": str(gram_determinant),
            "test_polynomial": "t-5/8",
            "hypothetical_symmetric_norm_squared": str(gram_determinant),
            "actual_L2_norm_squared": str(actual_norm_squared),
            "reciprocal_transfer_series": str(reciprocal_transfer),
            "transfer_identity": "L((I-lambda*T)^(-1)1)=(1/f(lambda)-1)/lambda",
            "scope": "No positive Hermitian form can both symmetrize T and retain [1,g]=c*integral(g), c>0; these beta moments are not the trace moments p_k"},
        "pointwise_Gram_obstruction": {
            "points_in_log_r_v_coordinates": [[str(value) for value in point] for point in points],
            "sample_poset": "A<C; B incomparable to both",
            "pair_principal_minors": list(map(str, pair_minors)),
            "sampled_T_determinant": "0",
            "Hermitian_phase_identity": "det=1-|u|^2-|v|^2",
            "forced_normalized_Gram_determinant": str(forced_determinant),
            "scope": "No positive Gram kernel can preserve these pointwise pair minors; a different integrated determinant representation is not excluded",
            "integrated_determinant_obstructed": False},
        "conditional_bridges": {
            "Hankel": {
                "status": "CONDITIONAL; positivity not proved",
                "required_forms": ["(p_(i+j+2))_(i,j>=0)", "(p_(i+j+3))_(i,j>=0)"],
                "p_k": "tr(T^k), k>=2",
                "consequence_if_proved": "Compactly supported Stieltjes representation continues -f'/f and excludes nonreal zeros",
                "positive_individual_traces_suffice": False},
            "Ferrers": {
                "status": "CONDITIONAL root implication; required real-rootedness not proved",
                "partitions": "lambda_(N,j)=floor(N*exp(-j/N)), positive rows only",
                "proved_limit_input": "M_(lambda_N)(z/N^2) converges locally uniformly to D(z)",
                "unproved_antecedent": "Real-rootedness of these ordinary Ferrers non-nesting rook polynomials",
                "consequence_if_proved": "All zeros of D are negative, hence the earlier all-size x=0 Jensen statement",
                "Ferrers_shapes_computed": 0}},
        "root_search": False, "root_counts_performed": 0,
        "size_or_rate_sweep": False, "new_size_or_rate_cases": 0,
        "all_real_spectrum_counterexample": False,
        "unrelated_certificate_helpers_called": False,
        "elapsed_seconds": time.monotonic()-started,
    }


def weighted_columns_path3_certificate():
    """One closed path at rate one: weighted no-loss bounds and transport.

    This tests averaging on the graph with the existing pointwise first-loss
    obstruction. It makes no universal kernel/rate assertion.
    """
    from exact import add, stable
    transition = ((0, 1, 0), (sp.Rational(1, 2), 0, sp.Rational(1, 2)), (0, 1, 0))
    kernel = tuple(tuple((j, p) for j, p in enumerate(row) if p) for row in transition)
    raw = leading_coefficients(transition, (1, 1, 1), marked=True)
    coefficients = raw["coefficients"]
    assert coefficients == leading_coefficients(transition, (1, 1, 1), reverse=True)["coefficients"]
    matrix, _ = rescaled_idla(transition, (1, 1, 1))
    inverse = solve(matrix, sp.eye(8))
    assert tuple(inverse[7 ^ holes, 7] for holes in range(8)) == coefficients
    total = sum(coefficients[1 << i] for i in range(3))
    assert total == 18
    starts = [add(stable(occupied, 3), vertex) for vertex in range(3) for occupied in range(7)]
    addition, count = absorption(starts, kernel, (sp.Rational(1, 2),)*3)
    opposite, opposite_count = absorption(starts, kernel, (sp.Rational(1, 2),)*3, reverse=True)
    assert addition == opposite
    pairs, singles = [6, 5, 3], [1, 2, 4]
    one_loss = inverse.extract(singles, pairs)
    source = -matrix.extract(pairs, [7])
    correction = matrix.extract(singles, [7])
    one_hole = inverse.extract(pairs, [7])
    rows, transport = [], []
    for vertex in range(3):
        for occupied in range(7):
            row = addition[7*vertex+occupied, :]
            assert sum(row) == 1 and min(row) >= 0
            assert all(row[output] == 0 for output in range(8)
                       if output.bit_count() != occupied.bit_count()+1)
        for holes in range(1, 8):
            rank = holes.bit_count()
            weighted = sum(coefficients[initial]*addition[7*vertex+(7 ^ initial), 7 ^ holes]
                           for initial in range(1, 8) if initial.bit_count() == rank+1)
            gap = total*coefficients[holes]-weighted
            stronger_gap = total*coefficients[holes]/(rank+1)-weighted
            assert gap > 0 and stronger_gap > 0
            rows.append({"entrance": vertex, "holes_mask": holes, "B": str(weighted),
                         "D": str(gap), "stronger_margin": str(stronger_gap),
                         "normalized_ratio": str(weighted/(total*coefficients[holes]))})
        restored = addition.extract([7*vertex+s for s in singles], pairs).T
        slack = one_hole*sp.ones(1, 3)-restored*one_loss
        assert min(slack) > 0
        assert slack*source+restored*inverse.extract(singles, singles)*correction == (
            total*one_hole-restored*inverse.extract(singles, [7]))
        transport.append({"entrance": vertex,
                          "closed_addition": [list(map(str, row)) for row in restored.tolist()],
                          "one_loss_transport_slack": [list(map(str, row)) for row in slack.tolist()]})
    assert min(sp.Rational(row["D"]) for row in rows) == sp.Rational(3509, 50)
    assert min(sp.Rational(row["stronger_margin"]) for row in rows) == sp.Rational(1349, 50)
    assert max(sp.Rational(row["normalized_ratio"]) for row in rows) == sp.Rational(11, 45)
    return {"scope": "one closed three-site path, uniform killing, common sleep rate one",
            "transition": [list(map(str, row)) for row in transition],
            "leading_coefficients_by_hole_mask": list(map(str, coefficients)),
            "M": str(total), "weighted_column_bounds": rows,
            "one_loss_all_active_columns": [list(map(str, row)) for row in one_loss.tolist()],
            "positive_source": list(map(str, source)),
            "positive_correction": list(map(str, correction)),
            "transport_checks": transport,
            "checks": {"raw_leading_recursion_both_priorities": True,
                       "independent_rescaled_IDLA_residual": True,
                       "all_21_closed_additions_both_priorities": True,
                       "addition_transient_counts": [count, opposite_count],
                       "closed_additions_normalized_nonnegative_conserve_particles": True,
                       "all_21_D_and_stronger_bounds_strict": True,
                       "all_27_transport_slacks_strict": True,
                       "positive_source_transport_identity": True},
            "general_kernel_transport_lemma_proved_by_computation": False}


def weighted_sources_path4_certificate():
    """One fixed source/cylinder discriminator; no raw absorption claim."""
    transition = ((0,1,0,0), (sp.Rational(1,2),0,sp.Rational(1,2),0),
                  (0,sp.Rational(1,2),0,sp.Rational(1,2)), (0,0,1,0))
    p, potential = closed_potentials(transition)
    matrix, operators = rescaled_idla(transition, (1,)*4)
    empty = sp.zeros(16,1)
    empty[0] = 1
    reversed_columns = [empty]
    for mask in range(1,16):
        vertex = mask.bit_length()-1
        reversed_columns.append(operators[vertex]*reversed_columns[mask ^ (1<<vertex)])
    assert matrix == sp.Matrix.hstack(*reversed_columns)
    inverse = solve(matrix, sp.eye(16))
    coefficients = [inverse[15 ^ holes,15] for holes in range(16)]
    assert coefficients[0] == 1 and min(coefficients) > 0
    one = [15 ^ (1<<i) for i in range(4)]
    pairs = [(i,j) for i in range(4) for j in range(i+1,4)]
    two = [15 ^ (1<<i) ^ (1<<j) for i,j in pairs]
    n1, n2 = inverse.extract(one,one), inverse.extract(two,two)
    for kernel in (n1,n2):
        assert min(kernel) >= 0 and sp.ones(1,kernel.rows)*kernel == sp.ones(1,kernel.cols)
    source = -matrix.extract(one,[15])
    explicit_source, total = one_vacancy_source(transition, (1,)*4)
    assert source == explicit_source and total == sum(source) == 38
    c1, c2 = inverse.extract(one,[15]), inverse.extract(two,[15])
    transport_source = -matrix.extract(two,one)
    correction = matrix.extract(two,[15])
    h = inverse.extract(two,one)
    paired_source = transport_source*c1
    assert h == n2*transport_source*n1 and c1 == n1*source
    assert c2 == h*source-n2*correction
    product_gaps = [source[i]*source[j]-paired_source[k] for k,(i,j) in enumerate(pairs)]
    assert min(product_gaps) == sp.Rational(8557,336) > 0
    rows = []
    for vertex in range(4):
        idla = sp.zeros(4,6)
        for col,(i,j) in enumerate(pairs):
            occupied = 15 ^ (1<<i) ^ (1<<j)
            inside, outside, _, _, exits = potential[occupied]
            for final, filled in ((i,j),(j,i)):
                idla[final,col] = int(vertex == filled) if vertex in outside else exits[inside.index(vertex),outside.index(filled)]
        assert sp.ones(1,4)*idla == sp.ones(1,6)
        bracket = total*source-idla*paired_source
        slack = n1*bracket
        gap = slack+n1*idla*correction
        assert min(bracket) > 0 and min(slack) > 0 and min(gap) > 0
        rows.append({'entrance':vertex, 'source_margin':list(map(str,bracket)),
                     'weighted_transport_slack':list(map(str,slack)), 'one_hole_D':list(map(str,gap))})
    assert min(sp.Rational(v) for row in rows for v in row['source_margin']) == sp.Rational(3631,21)
    assert min(sp.Rational(v) for row in rows for v in row['weighted_transport_slack']) == sp.Rational(34359515,150528)
    assert min(sp.Rational(v) for row in rows for v in row['one_hole_D']) == sp.Rational(12007189,50176)
    cylinders = [{'I':i,'J':j,'gap':str(coefficients[i]*coefficients[j]-coefficients[i|j])}
                 for i in range(1,16) for j in range(i+1,16) if not i&j]
    assert len(cylinders) == 25
    assert min(sp.Rational(row['gap']) for row in cylinders) == sp.Rational(7756751,235053) > 0
    k = matrix.extract(one,one)
    assert sp.ones(1,4)*k == sp.ones(1,4)
    assert all(k[i,j] < 0 for i in range(4) for j in range(4) if i!=j)
    generator = sp.eye(4)-k
    cycle = generator[0,1]*generator[1,2]*generator[2,0]-generator[0,2]*generator[2,1]*generator[1,0]
    inverse_cycle = n1[0,1]*n1[1,2]*n1[2,0]-n1[0,2]*n1[2,1]*n1[1,0]
    assert cycle == sp.Rational(1,8) and inverse_cycle == sp.Rational(1955,14450688)
    d = [potential[15 ^ (1<<i)][2] for i in range(4)]
    kappa, G = sum(d), sp.zeros(4,1)
    pi = [v/kappa for v in d]
    mixture_source = sp.zeros(4,1)
    mixture_mass = 0
    for target in range(1,16):
        inside, _, determinant, hitting, _ = potential[15 ^ target]
        weight = kappa/determinant
        mixture_mass += weight
        for vertex in range(4):
            if target >> vertex & 1:
                return_time = 1+sum(p[vertex,j]*hitting[t] for t,j in enumerate(inside))
                mixture_source[vertex] += weight*pi[vertex]*return_time
            else:
                G[vertex] += weight*hitting[inside.index(vertex)]
    assert mixture_source == source and mixture_mass == total
    scalar = sp.zeros(4,1)
    for target in range(4):
        inside, _, _, hitting, _ = potential[15 ^ (1<<target)]
        for t,vertex in enumerate(inside):
            scalar[vertex] += source[target]*hitting[t]
    assert scalar-G == sp.ones(4,1)*sp.Rational(93,2)
    return {'scope':'one closed path4, common rate one, physical uniform killing; exact leading data only',
            'leading_coefficients_by_hole_mask':list(map(str,coefficients)),
            'u':list(map(str,source)), 'M':str(total), 'G':list(map(str,G)),
            'source_potential_constant':'93/2', 'pair_order':pairs,
            'paired_source':list(map(str,paired_source)), 'product_source_gaps':list(map(str,product_gaps)),
            'entrance_checks':rows, 'disjoint_cylinder_gaps':cylinders,
            'one_hole_inverse_matrix':[list(map(str,row)) for row in k.tolist()],
            'one_hole_matrix':[list(map(str,row)) for row in n1.tolist()],
            'generator_cycle_difference':str(cycle), 'stabilization_cycle_difference':str(inverse_cycle),
            'checks':{'rescaled_IDLA_full_residual':True, 'both_admissible_factor_orders':True,
                      'separate_explicit_source_formula':True, 'source_mixture_and_potential':True,
                      'block_residuals_with_audited_intertwining':True, 'all_25_cylinder_gaps_positive':True,
                      'raw_absorption_performed':False, 'general_sign_theorem':False}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--audit", action="store_true", help="one rational rescaling/source audit")
    parser.add_argument("--path3-symbolic", action="store_true", help="two-rate path cylinder and NA proof")
    parser.add_argument("--cycle4-raw", action="store_true", help="raw four-cycle validation in two priorities")
    parser.add_argument("--uniform-renewal", action="store_true", help="positive uniform loss product, one symbolic n=3 audit")
    parser.add_argument("--uniform-rank-quotient", action="store_true", help="positive rank quotient and dimension identity, one symbolic n=3 audit")
    parser.add_argument("--uniform-parking-jensen", action="store_true", help="one exact n=16 x=0 Jensen root count and finite-prefix obstruction")
    parser.add_argument("--uniform-parking-spectrum", action="store_true", help="bounded compact-operator normalization, nearest-zero scope, and Gram obstructions at x=0")
    parser.add_argument("--weighted-columns-path3", action="store_true", help="one exact weighted no-loss/transport discriminator on the closed path at rate one")
    parser.add_argument("--weighted-sources-path4", action="store_true", help="one exact source/cylinder check and closed-stabilization reversibility obstruction")
    parser.add_argument("--output", type=Path,
                        default=Path(__file__).with_name("weak_killing_certificates.json"))
    args = parser.parse_args()
    if not (args.audit or args.path3_symbolic or args.cycle4_raw or args.uniform_renewal or args.uniform_rank_quotient or args.uniform_parking_jensen or args.uniform_parking_spectrum or args.weighted_columns_path3 or args.weighted_sources_path4):
        parser.error("choose a declared check, including --audit or --weighted-columns-path3")
    records = json.loads(args.output.read_text()) if args.output.exists() else {}
    for enabled, name, check in (
        (args.audit, "rescaled_audit", rescaled_audit),
        (args.path3_symbolic, "path3_symbolic", path3_symbolic),
        (args.cycle4_raw, "cycle4_raw", cycle4_raw),
        (args.uniform_renewal, "uniform_renewal", uniform_renewal_certificate),
        (args.uniform_rank_quotient, "uniform_rank_quotient", uniform_rank_quotient_certificate),
        (args.uniform_parking_jensen, "uniform_parking_jensen", uniform_parking_jensen_certificate),
        (args.uniform_parking_spectrum, "uniform_parking_spectrum", uniform_parking_spectrum_certificate),
        (args.weighted_columns_path3, "weighted_columns_path3", weighted_columns_path3_certificate),
        (args.weighted_sources_path4, "weighted_sources_path4", weighted_sources_path4_certificate),
    ):
        if enabled:
            records[name] = check()
            args.output.write_text(json.dumps(records, indent=2)+"\n")
            print(name+": exact checks passed", flush=True)


if __name__ == "__main__":
    main()
