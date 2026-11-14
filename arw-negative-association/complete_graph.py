#!/usr/bin/env python3
"""Two independent exact descriptions of complete-graph stationary counts.

Harmonic interpolation is derived in report.md. The binomial-refresh chain
is Proposition 5 of Junge--Kaufman--Meisel, arXiv:2604.04747v2 (check source
label prop:stop rather than relying on numbering). Its self-loop model is
converted to simple random walk by changing the sleep probability.

Purpose: select an all-size proof subgoal: ultra-log-concavity of rank
probabilities, not real-rootedness, which already fails at n=8, lambda=1.
"""

import argparse
import json

from exact import absorption, add, graph_kernel, sp, stable, stationary_law


def interpolate(n, rate):
    rate = sp.Rational(rate)
    matrix = sp.ones(n+1, n+1)
    for j in range(n):
        aa = sp.Rational(1, n+1-j)
        bb = -sp.Rational(n+1, n*(n+1-j))/rate
        for k in range(n+1):
            matrix[j, k] = sum(
                sp.binomial(j, ell)*sp.binomial(n-j, k-ell)*aa**ell*bb**(k-ell)
                for ell in range(k+1)) / sp.binomial(n, k)
    target = sp.zeros(n+1, 1)
    target[n] = 1
    rank = matrix.inv()*target
    assert matrix*rank == target
    assert sum(rank) == 1 and all(p > 0 for p in rank)
    return tuple(rank)


def rank_connection_identity(n, rate, rank):
    """Verify the all-size polynomial reformulation of harmonic interpolation.

    K_k(x)=[u^k](1+u)^(n-x)(1-beta*u)^x / binomial(n,k).
    The desired rank coefficients are the connection coefficients between
    product_{j=1}^n(j-x)/n! and (x+1)^(n-k) K_k(x). Their ULC remains open.
    """
    x = sp.Symbol("x")
    beta = sp.Rational(n+1, n)/sp.Rational(rate)
    polynomials = [sp.Poly(1, x)]
    if n:
        polynomials.append(sp.Poly((n-(1+beta)*x)/n, x))
    for k in range(1, n):
        polynomials.append(((n-(1+beta)*x+(beta-1)*k)*polynomials[k]
                            -beta*k*polynomials[k-1]).exquo_ground(n-k))
    value = sum((sp.Poly((x+1)**(n-k), x)*polynomials[k]*rank[k]
                 for k in range(n+1)), sp.Poly(0, x))
    target = sp.Poly(sp.prod(j-x for j in range(1, n+1))/sp.factorial(n), x)
    assert value == target
    # Positive lowest-order coefficients at x=-1 make this basis invertible.
    assert all(polynomials[k].eval(-1)
               == sum(sp.binomial(n+1, k-j)*beta**j for j in range(k+1))/sp.binomial(n, k)
               for k in range(n+1))
    assert rank[-1] == (n+1)*beta/((1+beta)**(n+1)-1)
    return {"identity_verified": True, "basis_positive_at_minus_one": True,
            "full_occupation_formula_verified": True}


def binomial_refresh(n, rate):
    """Exact boundary hitting law of a birth/death chain and a falling bound."""
    # Loop erasure changes sleep rate by (n+1)/n. Thus start the looped
    # model at rate n*rate/(n+1), with sink probability 1/(n+1).
    prob = n*sp.Rational(rate)/(n+1+n*sp.Rational(rate))
    initial = [sp.binomial(n, k)*prob**k*(1-prob)**(n-k) for k in range(n+1)]
    result = [sp.S.Zero] * (n+1)
    result[n] = initial[n]
    incoming = sp.Matrix(initial[:n])
    # Each coordinate refreshes at rate one; sink events have total rate one.
    for boundary in range(n, 0, -1):
        killed = sp.zeros(boundary)
        for y in range(boundary):
            birth, death = prob*(n-y), (1-prob)*y
            killed[y, y] = birth+death+1
            if y+1 < boundary:
                killed[y, y+1] = -birth
            if y > 0:
                killed[y, y-1] = -death
        visits = killed.T.inv()*incoming
        assert killed.T*visits == incoming
        assert all(v >= 0 for v in visits)
        result[boundary] += visits[boundary-1]*prob*(n-boundary+1)
        result[boundary-1] += visits[boundary-1]
        incoming = visits[:boundary-1, :]
    assert sum(result) == 1
    return tuple(result)


def boundary_profiles(n, prob):
    """Backward survival densities and their boundary diagonal.

    Concrete stronger subclaim: a_k=f_{k+1}(k) is log-concave.  The algebra
    in report.md shows that this implies ULC of the exit law.  These finite
    checks do not prove that subclaim for arbitrary n and prob.
    """
    prob = sp.Rational(prob)
    qq = 1-prob
    row = sp.ones(n+1, 1)
    diagonal = [None]*n + [sp.S.One]
    for boundary in range(n, 0, -1):
        killed = sp.zeros(boundary)
        for y in range(boundary):
            birth, death = prob*(n-y), qq*y
            killed[y, y] = 1+birth+death
            if y+1 < boundary:
                killed[y, y+1] = -birth
            if y:
                killed[y, y-1] = -death
        row = killed.inv()*row[:boundary, :]
        diagonal[boundary-1] = row[boundary-1]
        # Closed determinant formula, proved by the continuant recurrence.
        series = sum(sp.binomial(n-boundary+j, j)*prob**j
                     for j in range(boundary+1))
        assert killed.det() == sp.factorial(boundary)*series
    relative_flux = tuple(diagonal[k] + (qq*k*diagonal[k-1] if k else 0)
                          for k in range(n+1))
    return tuple(diagonal), relative_flux


def rank_basis(n, beta):
    """Symmetric IDLA change of basis, with beta=(n+1)/(n*lambda)."""
    addition = sp.zeros(n+1)
    for j in range(n+1):
        addition[j, j] = sp.Rational(1, n+1-j)
        if j < n:
            addition[j+1, j] = sp.Rational(n-j, n+1-j)
    powers = [sp.eye(n+1)]
    for _ in range(n):
        powers.append(addition*powers[-1])
    basis = sp.zeros(n+1)
    for k in range(n+1):
        basis[:, k] = sum((sp.binomial(k, ell)*(1+beta)**ell*(-beta)**(k-ell)
                           *powers[k-ell][:, ell] for ell in range(k+1)), sp.zeros(n+1, 1)).applyfunc(sp.expand)
    return basis


def active_rank_kernel(n, rate):
    """All-active stabilization on the symmetric (n+1)-dimensional space.

    Column k starts with one active at each of k sites. The reduction from
    the full IDLA algebra is proved in report.md; its top column is the
    stationary rank law, but the other columns allow a stronger ULC lemma.
    """
    beta = sp.Rational(n+1, n)/sp.Rational(rate)
    basis = rank_basis(n, beta)
    kernel = basis.inv()
    assert basis*kernel == sp.eye(n+1)
    assert all(value >= 0 for value in kernel)
    assert all(sum(kernel[:, k]) == 1 for k in range(n+1))
    return kernel


def lc_preservation_certificate(kernel):
    """Sufficient finite-matrix criterion for every ULC input rank law.

    Group the output LC-gap quadratic form by k+l. LC input products
    increase toward the middle, so nonnegative inward coefficient sums
    suffice by summation by parts (Wang--Yeh's LC-positivity criterion).
    """
    n = kernel.rows-1
    normalized = sp.Matrix(n+1, n+1, lambda j, k: kernel[j, k]*sp.binomial(n, k)/sp.binomial(n, j))
    failures, checked, zero = [], 0, 0
    for j in range(1, n):
        for total in range(2*n+1):
            cumulative = sp.S.Zero
            for k in range(total//2, max(0, total-n)-1, -1):
                ell = total-k
                if k == ell:
                    coefficient = normalized[j, k]**2-normalized[j-1, k]*normalized[j+1, k]
                else:
                    coefficient = (2*normalized[j, k]*normalized[j, ell]
                                   -normalized[j-1, k]*normalized[j+1, ell]
                                   -normalized[j-1, ell]*normalized[j+1, k])
                cumulative += coefficient
                checked += 1
                zero += int(cumulative == 0)
                if cumulative < 0:
                    failures.append({"row": j, "total": total, "cut": k, "value": str(cumulative)})
    return {"coefficient_sums_checked": checked, "zero_sums": zero,
            "ULC_preservation_certified": not failures, "failed_sums": failures,
            "rank_transition_matrix": [[str(kernel[j, k]) for k in range(n+1)] for j in range(n+1)]}


def symbolic_lc_certificate(n):
    """Certify the same all-input ULC criterion for every beta>0.

    Purpose: determine whether a coefficientwise argument can establish
    complete-graph NA uniformly in the rate despite failed real stability.
    """
    from exact import DomainMatrix

    beta = sp.Symbol("beta")
    kernel = DomainMatrix.from_Matrix(rank_basis(n, beta)).to_field().inv().to_Matrix()
    normalized = sp.Matrix(n+1, n+1, lambda j, k: sp.cancel(
        kernel[j, k]*sp.binomial(n, k)/sp.binomial(n, j)))
    parts = [tuple(sp.Poly(value, beta) for value in sp.fraction(entry)) for entry in normalized]
    denominator = sp.Poly(1, beta)
    for _, den in parts:
        denominator = sp.lcm(denominator, den)
    if denominator.LC() < 0:
        denominator = -denominator
    numerators = [[parts[j*(n+1)+k][0]*denominator.exquo(parts[j*(n+1)+k][1])
                   for k in range(n+1)] for j in range(n+1)]
    failures, checked = [], 0
    for j in range(1, n):
        for total in range(2*n+1):
            cumulative = sp.Poly(0, beta)
            for k in range(total//2, max(0, total-n)-1, -1):
                ell = total-k
                if k == ell:
                    term = numerators[j][k]**2-numerators[j-1][k]*numerators[j+1][k]
                else:
                    term = (2*numerators[j][k]*numerators[j][ell]
                            -numerators[j-1][k]*numerators[j+1][ell]
                            -numerators[j-1][ell]*numerators[j+1][k])
                cumulative += term
                checked += 1
                if any(value < 0 for value in cumulative.coeffs()):
                    failures.append({"row": j, "total": total, "cut": k,
                                     "polynomial": str(sp.factor(cumulative.as_expr()))})
    positive_denominator = all(value >= 0 for value in denominator.coeffs())
    return {"n": n, "parameter": "beta=(n+1)/(n*lambda)", "coefficient_sums_checked": checked,
            "positive_denominator_coefficients": positive_denominator,
            "all_rates_ULC_preservation_certified": positive_denominator and not failures,
            "failed_coefficient_tests": failures,
            "normalized_kernel_denominator": str(sp.factor(denominator.as_expr())),
            "normalized_kernel_numerators": [[str(sp.factor(value.as_expr())) for value in row]
                                              for row in numerators]}


def check(n, rate):
    rank = interpolate(n, rate)
    assert rank == binomial_refresh(n, rate)
    prob = n*sp.Rational(rate)/(n+1+n*sp.Rational(rate))
    diagonal, relative_flux = boundary_profiles(n, prob)
    assert all(rank[k] == sp.binomial(n, k)*prob**k*(1-prob)**(n-k)*relative_flux[k]
               for k in range(n+1))
    diagonal_gaps = [diagonal[k]**2-diagonal[k-1]*diagonal[k+1]
                     for k in range(1, n)]
    if n <= 4:
        adjacency = tuple(tuple(j for j in range(n+1) if j != i) for i in range(n))
        p = sp.Rational(rate)/(1+sp.Rational(rate))
        law, _ = stationary_law(graph_kernel(adjacency), (p,)*n)
        assert all(rank[k] == sum(v for mask, v in enumerate(law) if mask.bit_count() == k)
                   for k in range(n+1))
    normalized = [rank[k]/sp.binomial(n, k) for k in range(n+1)]
    gaps = [normalized[k]**2-normalized[k-1]*normalized[k+1] for k in range(1, n)]
    z = sp.Symbol("z")
    poly = sp.Poly(sum(v*z**k for k, v in enumerate(rank)), z)
    # A count of distinct real roots rules out real-rootedness only after
    # excluding repeated factors. Record that exact check explicitly.
    assert poly.gcd(poly.diff()).degree() == 0
    real_roots = int(poly.count_roots(-sp.oo, sp.oo))
    result = {"n": n, "lambda": str(rate), "rank": [str(v) for v in rank],
            "independent_refresh_check": True, "square_free": True, "real_roots": real_roots,
            "survival_diagonal": [str(v) for v in diagonal],
            "survival_diagonal_LC": all(v >= 0 for v in diagonal_gaps),
            "ULC": all(v >= 0 for v in gaps), "ULC_gaps": [str(v) for v in gaps]}
    result["rank_connection_identity"] = rank_connection_identity(n, rate, rank)
    if n == 3 and rate == 1:
        # A narrower dynamics-preservation route also fails: even an
        # exchangeable deterministic start can acquire positive covariance.
        adjacency = tuple(tuple(j for j in range(n+1) if j != i) for i in range(n))
        output, count = absorption([add(stable((1 << n)-1, n), 0)],
                                   graph_kernel(adjacency), (sp.Rational(1, 2),)*n)
        ranks = [sum(output[0, mask] for mask in range(1 << n) if mask.bit_count() == k)
                 for k in range(n+1)]
        marginal = sum(k*ranks[k]/n for k in range(n+1))
        pair = sum(k*(k-1)*ranks[k]/(n*(n-1)) for k in range(n+1))
        covariance = pair-marginal**2
        assert covariance == sp.Rational(7366296671, 853665123600)
        result["uniform_addition_counterexample"] = {
            "initial": "all three sites asleep", "lambda": "1",
            "rank": list(map(str, ranks)), "pair_covariance": str(covariance),
            "raw_absorption_states": count,
            "scope": "one uniform addition, not the stationary law"}
    kernel = active_rank_kernel(n, rate)
    assert tuple(kernel[:, -1]) == rank
    result["all_active_ULC_preservation"] = lc_preservation_certificate(kernel)
    if n == 3 and rate == 1:
        initial = [tuple(int(i < k) for i in range(n)) for k in range(n+1)]
        raw, _ = absorption(initial, graph_kernel(adjacency), (sp.Rational(1, 2),)*n)
        for k in range(n+1):
            assert tuple(sum(raw[k, mask] for mask in range(1 << n) if mask.bit_count() == j)
                         for j in range(n+1)) == tuple(kernel[:, k])
        result["all_active_ULC_preservation"]["independent_raw_rank_columns_checked"] = n+1
    return result


def renewal_diagonal(n, probability, boundary_rate=1):
    """Independent passage-count / backward-queue form of the diagonal.

    G_j(u)=beta_j/[beta_j+delta_j+alpha(1-u)-delta_j G_{j-1}(u)]
    is the PGF of the Poisson(alpha * passage_time(j,j+1)) count.
    Only coefficients through degree n are needed to calculate the queue
    B_j=(B_{j+1}+1-Z_j)_+, B_n=0, exactly.
    """
    probability, alpha = map(sp.Rational, (probability, boundary_rate))
    assert n >= 1 and 0 < probability <= 1 and alpha > 0
    coefficients = []
    previous = [sp.S.Zero]*(n+1)
    for j in range(n):
        birth, death = probability*(n-j), (1-probability)*j
        denominator = [birth+death+alpha-death*previous[0]]
        denominator += [-death*previous[r]-(alpha if r == 1 else 0)
                        for r in range(1, n+1)]
        assert denominator[0] > 0
        current = [birth/denominator[0]]
        for r in range(1, n+1):
            current.append(-sum(denominator[s]*current[r-s]
                                for s in range(1, r+1))/denominator[0])
        assert all(value >= 0 for value in current) and sum(current) <= 1
        coefficients.append(current)
        previous = current
    queue = [sp.S.One]
    diagonal = [sp.S.Zero]*n+[sp.S.One]
    queues = [None]*n+[queue]
    for k in range(n-1, -1, -1):
        count = coefficients[k]
        updated = [sum(queue[s]*(1-sum(count[:s+1]))
                       for s in range(len(queue)))]
        updated += [
            sum(queue[s]*count[s+1-r] for s in range(r-1, len(queue)))
            for r in range(1, len(queue)+1)]
        assert sum(updated) == 1 and all(value >= 0 for value in updated)
        diagonal[k] = updated[0]
        queues[k] = updated
        queue = updated
    return tuple(diagonal), coefficients, queues


def renewal_certificate():
    """Check the new representation, and rule out generic rank-flow ULC.

    Two boundary clock rates discriminate the Poisson parameter and queue
    indexing; these support the general derivation, not an all-size LC claim.
    """
    records = []
    n, probability = 4, sp.Rational(2, 5)
    for alpha in (sp.S.One, sp.Rational(2)):
        diagonal, coefficients, queues = renewal_diagonal(n, probability, alpha)
        row = sp.ones(n+1, 1)
        direct = [None]*n+[sp.S.One]
        for boundary in range(n, 0, -1):
            matrix = sp.zeros(boundary)
            for y in range(boundary):
                birth, death = probability*(n-y), (1-probability)*y
                matrix[y, y] = alpha+birth+death
                if y+1 < boundary:
                    matrix[y, y+1] = -birth
                if y:
                    matrix[y, y-1] = -death
            row = matrix.inv()*row[:boundary, :]*alpha
            direct[boundary-1] = row[boundary-1]
        assert diagonal == tuple(direct)
        if alpha == 1:
            assert diagonal == boundary_profiles(n, probability)[0]
        records.append({
            "n": n, "refresh_probability": str(probability),
            "boundary_rate": str(alpha), "diagonal": list(map(str, diagonal)),
            "passage_count_coefficients_through_degree_n": [
                list(map(str, row)) for row in coefficients],
            "backward_queue_distributions": [list(map(str, row)) for row in queues],
            "killed_resolvent_diagonal_agrees": True})
    # The actual stationary reinjection family is not an arbitrary ULC input.
    n, probability, density = 2, sp.Rational(1, 2), sp.Rational(3, 4)
    addition = sp.eye(n+1)
    for j in range(n+1):
        up, down = probability*(n-j)/n, (1-probability)*j/n
        addition[j, j] -= up+down
        if j < n:
            addition[j+1, j] += up
        if j:
            addition[j-1, j] += down
    rank = sp.Matrix([sp.binomial(n, j)*density**j*(1-density)**(n-j)
                      for j in range(n+1)])
    field = (sp.eye(n+1)-addition)*sp.diag(*range(n+1))*rank
    weight = [rank[j]/sp.binomial(n, j) for j in range(n+1)]
    velocity = [field[j]/sp.binomial(n, j) for j in range(n+1)]
    assert weight[1]**2 == weight[0]*weight[2]
    slope = 2*weight[1]*velocity[1]-velocity[0]*weight[2]-weight[0]*velocity[2]
    assert field == sp.Matrix([-sp.Rational(3, 32), -sp.Rational(3, 8),
                              sp.Rational(15, 32)])
    assert slope == -sp.Rational(3, 64)
    return {
        "scope": "finite checks of the proved renewal identity; LC remains open",
        "renewal_cases": records,
        "reinjection_ULC_cone_obstruction": {
            "n": n, "refresh_probability": str(probability),
            "input_binomial_density": str(density),
            "input_rank_law": list(map(str, rank)),
            "vector_field": list(map(str, field)),
            "normalized_LC_gap_derivative": str(slope),
            "input_is_not_the_stationary_law_of_this_addition_kernel": True}}


def reinjection_rank(n, probability, feedback):
    """Rank law and addition kernel for K(i,j)=feedback/n.

    The reference K=0 law is binomial. The finite reinjection transform
    supplies the current rank law and the intertwined addition kernel.
    This routine is used only at exact rational parameter points.
    """
    probability, feedback = map(sp.Rational, (probability, feedback))
    initial_addition = sp.eye(n+1)
    for j in range(n+1):
        up, down = probability*(n-j)/n, (1-probability)*j/n
        initial_addition[j, j] -= up+down
        if j < n:
            initial_addition[j+1, j] += up
        if j:
            initial_addition[j-1, j] += down
    base = (sp.eye(n+1)-feedback*initial_addition)/(1-feedback)
    matrix = sp.Matrix.hstack(*[(base**(n-j))[:, j] for j in range(n+1)])
    initial = sp.Matrix([sp.binomial(n, j)*probability**j*(1-probability)**(n-j)
                         for j in range(n+1)])
    inverse = matrix.inv()
    rank = inverse*initial
    addition = inverse*((1-feedback)*initial_addition
                        *(sp.eye(n+1)-feedback*initial_addition).inv())*matrix
    derivative = ((sp.eye(n+1)-addition)*sp.diag(*range(n+1))*rank
                  /(1-feedback))
    assert sum(rank) == 1 and all(value > 0 for value in rank)
    assert all(value >= 0 for value in addition)
    assert addition*rank == rank
    assert all(sum(addition[:, j]) == 1 for j in range(n+1))
    return rank, addition, derivative


def reverse_boundary(n, probability, tau):
    """Rank law, feedback derivative, and reverse mean without addition matrices.

    The exact rational inputs satisfy n>=1, 0<probability<1, tau>0,
    with feedback z=n*tau/(1+n*tau). The returned derivative is with
    respect to z at fixed probability. Two scalar tridiagonal solves
    per boundary require O(n^2) arithmetic and O(n) working storage.
    Their sensitivities are -tau*partial_tau of the survival profiles.
    No concavity claim follows from positivity of these recursions.
    """
    probability, tau = map(sp.Rational, (probability, tau))
    assert n >= 1 and 0 < probability < 1 and tau > 0
    qq = 1-probability
    birth = [probability*(n-j) for j in range(n+1)]
    death = [qq*j for j in range(n+1)]
    diagonal = [1+tau*(birth[j]+death[j]) for j in range(n+1)]
    pivots = [diagonal[0]]
    for j in range(1, n):
        pivots.append(diagonal[j]-tau**2*death[j]*birth[j-1]/pivots[-1])
    assert all(value > 0 for value in pivots)

    def solve(rhs):
        size = len(rhs)
        forward = [rhs[0]]
        for j in range(1, size):
            forward.append(rhs[j]+tau*death[j]*forward[-1]/pivots[j-1])
        answer = [sp.S.Zero]*size
        answer[-1] = forward[-1]/pivots[size-1]
        for j in range(size-2, -1, -1):
            answer[j] = (forward[j]+tau*birth[j]*answer[j+1])/pivots[j]
        return answer

    survival = [None]*n+[sp.S.One]
    sensitivity = [None]*n+[sp.S.Zero]
    previous = [sp.S.One]*(n+1)
    previous_sensitivity = [sp.S.Zero]*(n+1)
    for boundary in range(n, 0, -1):
        current = solve(previous[:boundary])
        current_sensitivity = solve([
            previous_sensitivity[j]+previous[j]-current[j]
            for j in range(boundary)])
        assert all(0 < current[j] <= previous[j] for j in range(boundary))
        assert all(value >= 0 for value in current_sensitivity)
        survival[boundary-1] = current[-1]
        sensitivity[boundary-1] = current_sensitivity[-1]
        previous, previous_sensitivity = current, current_sensitivity

    density = [survival[j]+(tau*death[j]*survival[j-1] if j else 0)
               for j in range(n+1)]
    scaled_derivative = [
        -sensitivity[j]+(tau*death[j]*(survival[j-1]-sensitivity[j-1]) if j else 0)
        for j in range(n+1)]
    initial = [sp.binomial(n, j)*probability**j*qq**(n-j) for j in range(n+1)]
    feedback = n*tau/(1+n*tau)
    rank = sp.Matrix([initial[j]*density[j] for j in range(n+1)])
    derivative = sp.Matrix([
        initial[j]*scaled_derivative[j]*(1+n*tau)**2/(n*tau)
        for j in range(n+1)])
    reverse_mean = sp.Matrix([
        j-scaled_derivative[j]/(feedback*density[j]) for j in range(n+1)])
    assert sum(rank) == 1 and sum(derivative) == 0
    assert all(value > 0 for value in rank)
    return rank, derivative, reverse_mean


def reverse_boundary_certificate():
    """One exact audit of scalar rank, sensitivity, and conditional moments.

    The interior point n=4,p=2/5,tau=3/2 exercises a nonzero q-p term
    in the moment hierarchy. It checks identities, not curvature signs.
    """
    n, probability, tau = 4, sp.Rational(2, 5), sp.Rational(3, 2)
    feedback = n*tau/(1+n*tau)
    rank, derivative, reverse_mean = reverse_boundary(n, probability, tau)
    reference, addition, reference_derivative = reinjection_rank(n, probability, feedback)
    assert rank == reference and derivative == reference_derivative
    count = sp.diag(*range(n+1))
    weights = sp.diag(*rank)
    reverse_weights = sp.diag(*[1/value for value in rank])
    assert reverse_mean == reverse_weights*addition*count*rank

    # Conjugating D diag(b)=I+tau U+tau D^2 V gives the exact
    # stationary reverse conditional identity for every test function.
    averaged = (1-feedback)*sp.eye(n+1)+feedback*addition
    backward = reverse_weights*averaged*weights
    assert backward*sp.ones(n+1, 1) == sp.ones(n+1, 1)
    lower = [probability*(n-j+1)*rank[j-1]/rank[j] if j else sp.S.Zero
             for j in range(n+1)]
    upper = [(1-probability)*(j+1)*rank[j+1]/rank[j] if j < n else sp.S.Zero
             for j in range(n+1)]
    for degree in range(3):
        moment = backward*sp.Matrix([j**degree for j in range(n+1)])
        next_moment = backward*sp.Matrix([j**(degree+1) for j in range(n+1)])
        weighted_two_step = backward**2*sp.Matrix([
            upper[j]*(j+1)**degree for j in range(n+1)])
        left = (1+tau*n*probability)*moment+tau*(1-2*probability)*next_moment
        right = sp.Matrix([j**degree+tau*lower[j]*(j-1)**degree
                           for j in range(n+1)])+tau*weighted_two_step
        assert left == right
    return {
        "scope": "one exact audit of all-size scalar identities; no curvature sign claim",
        "n": n, "p": str(probability), "tau": str(tau), "z": str(feedback),
        "rank_law": list(map(str, rank)),
        "derivative_wrt_feedback": list(map(str, derivative)),
        "reverse_addition_mean": list(map(str, reverse_mean)),
        "rank_and_derivative_match_reinjection_rank": True,
        "reverse_mean_matches_addition_kernel": True,
        "conditional_moment_degrees_checked": [0, 1, 2],
        "all_moment_residuals_zero": True}


def rank_top_structure(n, probability, tau):
    """Top three masses and curvature from positive resolvent integrals.

    Inputs are exact, n>=2, 0<p<1 and tau>0; here x=(1-p)/p.
    The logarithmic derivative is in tau at fixed p. These formulas
    do not by themselves prove that the top curvature is nonnegative.
    """
    probability, tau = map(sp.Rational, (probability, tau))
    assert n >= 2 and 0 < probability < 1 and tau > 0
    x = (1-probability)/probability
    t = sp.Symbol("t", positive=True)
    resolvents = [sp.factorial(a)*t**a*sum(
        sp.binomial(n-a, k)*x**k
        /sp.prod(1+t*(k+ell) for ell in range(a+1))
        for k in range(n-a+1)) for a in range(3)]
    first = -t*sp.diff(resolvents[0], t)
    second = (first+t*sp.diff(first, t)
              -2*t*first*sp.diff(resolvents[1], t)/resolvents[1])/2
    evaluate = lambda value: sp.cancel(value.subs(t, tau))
    top = [evaluate(1/resolvents[0]),
           evaluate(first/(resolvents[0]*resolvents[1])),
           evaluate(second/(resolvents[0]*resolvents[2]))]
    log_slope = evaluate(2*sp.diff(first, t)/first
                         +sp.diff(resolvents[2], t)/resolvents[2]
                         -sp.diff(second, t)/second
                         -2*sp.diff(resolvents[1], t)/resolvents[1])
    ww, vv, uu = resolvents[1]/t, first/(n*x*t), resolvents[2]/(2*t**2)
    pair = sp.diff(vv, t)*ww-2*vv*sp.diff(ww, t)
    return {"top_rank": top, "resolvents": list(map(evaluate, resolvents)),
            "F": evaluate(first), "J": evaluate(second),
            "U": evaluate(uu), "V": evaluate(vv), "W": evaluate(ww),
            "Q": evaluate(pair), "log_ratio_tau_derivative": log_slope,
            "reverse_mean_top_curvature": (1+n*tau)*log_slope/n}


def rank_top_structure_certificate():
    """Symbolic pair-kernel identity and one rational top-formula audit."""
    i, j, t, sigma, product = sp.symbols("i j t sigma product")
    den = lambda k: 1+t*k
    theta = lambda k: den(k)/den(k+1)
    gg = lambda k: k/den(k)+(k+1)/den(k+1)
    hh = lambda k: (k+1)/den(k+1)
    kernel = theta(i)*(gg(j)-hh(i))+theta(j)*(gg(i)-hh(j))
    numerator = (
        sigma+t*(2*sigma**2+3*sigma+2*product)
        +t**2*(sigma**3+3*sigma**2+6*sigma*product+3*sigma+12*product)
        +t**3*(4*sigma**2*product+2*sigma**2+18*sigma*product
                +sigma+4*product**2+10*product)
        +t**4*(3*sigma**2*product+5*sigma*product**2+9*sigma*product
                +12*product**2+4*product)
        +t**5*(sigma**2*product+3*sigma*product**2+sigma*product
                +2*product**3+2*product**2))
    denominator = den(i)*den(j)*den(i+1)**2*den(j+1)**2
    assert sp.cancel(kernel-numerator.subs({sigma: i+j, product: i*j})/denominator) == 0
    assert all(value >= 0 for value in sp.Poly(numerator, sigma, product, t).coeffs())

    # Reuse the existing scalar audit point; do not run its certificate
    # or the small-size log-curvature certificates again.
    n, probability, tau = 4, sp.Rational(2, 5), sp.Rational(3, 2)
    values = rank_top_structure(n, probability, tau)
    rank, derivative, mean = reverse_boundary(n, probability, tau)
    assert values["top_rank"] == [rank[n], rank[n-1], rank[n-2]]
    assert values["reverse_mean_top_curvature"] == mean[n-2]+mean[n]-2*mean[n-1]
    assert values["log_ratio_tau_derivative"] == n/(1+n*tau)**2*(
        2*derivative[n-1]/rank[n-1]-derivative[n-2]/rank[n-2]-derivative[n]/rank[n])
    x = (1-probability)/probability
    weights = [sp.binomial(n-1, k)*x**k/((1+tau*k)*(1+tau*(k+1)))
               for k in range(n)]
    pair_sum = sum(weights[a]*weights[b]*kernel.subs({i: a, j: b, t: tau})
                   for a in range(n) for b in range(n))
    assert values["Q"] == pair_sum
    assert 2*(n-1)*x*values["U"]*values["V"]**2/(values["W"]*values["Q"]) == (
        (rank[n-1]/n)**2/(rank[n]*rank[n-2]/sp.binomial(n, 2)))
    return {
        "scope": "all-index pair-kernel identity and one rational formula audit; all-rate curvature remains open",
        "pair_kernel": {"sigma": "i+j", "product": "i*j", "parameter": "t=tau",
                        "numerator": str(numerator), "denominator": str(denominator),
                        "symbolic_identity_verified": True,
                        "numerator_coefficients_nonnegative": True},
        "rational_audit": {
            "n": n, "p": str(probability), "tau": str(tau),
            "z": str(n*tau/(1+n*tau)), "rank_indices": [n, n-1, n-2],
            "top_rank_law": list(map(str, values["top_rank"])),
            "resolvents_Z0_Z1_Z2": list(map(str, values["resolvents"])),
            "log_ratio_tau_derivative": str(values["log_ratio_tau_derivative"]),
            "reverse_mean_top_curvature": str(values["reverse_mean_top_curvature"]),
            "positive_pair_sum": str(values["Q"]),
            "rank_and_curvature_match_scalar_boundary": True,
            "pair_sum_and_top_ratio_match": True}}


def rank_top_coefficients_certificate():
    """All-index mixed kernels proving top completed-row coefficient TP2.

    For D=product_(i=1)^n(1+i*tau), the reduced top rows are
    f=sum b_i*(D/d_(i+1))^2 and g=sum b_i*D^2/(d_i*d_(i+1)).
    Their coefficient minors are quadratic in arbitrary nonnegative b_i.
    Removing common linear factors reduces the mixed terms to these
    two symbolic tables. Toeplitz convolution by those factors, and by
    product_(b=1)^(n-2) Delta_b, preserves every mixed minor coefficient.
    This proves coefficient TP2 for the top pair at every size and rate;
    the separate n=5,p=1/2 check audits only the numerator formulas.
    """
    from pathlib import Path
    a, h, tau, u, v = sp.symbols("a h tau u v", nonnegative=True)
    dd = lambda k: 1+k*tau
    cores = [
        ("overlap", dd(a)*(u*dd(a+2)**2+v*dd(a+1)**2),
         dd(a+1)*dd(a+2)*(u*dd(a+2)+v*dd(a))),
        ("separated", dd(a)*dd(a+2+h)*(u*dd(a+3+h)**2+v*dd(a+1)**2),
         dd(a+1)*dd(a+3+h)*(u*dd(a+2+h)*dd(a+3+h)+v*dd(a)*dd(a+1)))]
    tables = {}
    for label, first, second in cores:
        first_poly, second_poly = sp.Poly(first, tau), sp.Poly(second, tau)
        degree = max(first_poly.degree(), second_poly.degree())
        entries = {}
        for lo in range(degree+1):
            for hi in range(lo+1, degree+1):
                mixed = sp.expand(first_poly.nth(lo)*second_poly.nth(hi)
                                  -first_poly.nth(hi)*second_poly.nth(lo)).coeff(u, 1).coeff(v, 1)
                assert all(value >= 0 for value in sp.Poly(mixed, a, h).coeffs())
                entries[f"{lo},{hi}"] = str(sp.factor(mixed))
        assert len(entries) == (6 if label == "overlap" else 10)
        tables[label] = {"F": str(first), "G": str(second), "mixed_minors": entries,
                         "coefficients_nonnegative": True}

    # Compare only the two formulas against the saved exact coefficient
    # rows. Do not rerun the 600-minor fixed-size certificate here.
    reference = Path(__file__).with_name("complete_graph_symbolic_certificates.json")
    saved = json.loads(reference.read_text())["rank_likelihood_order"]
    n, probability = 5, sp.Rational(1, 2)
    assert saved["n"] == n and saved["p"] == str(probability)
    generator = sp.zeros(n+1)
    for j in range(n+1):
        birth, death = probability*(n-j), (1-probability)*j
        generator[j, j] = -birth-death
        if j < n:
            generator[j, j+1] = birth
        if j:
            generator[j, j-1] = death
    matrix = sp.eye(n+1)-tau*generator.T
    dd_product = sp.prod(dd(j) for j in range(1, n+1))
    common = sp.prod(matrix[:b, :b].det() for b in range(1, n-1))
    weights = [sp.binomial(n-1, j)*probability**(n-1-j)*(1-probability)**j
               for j in range(n)]
    first = sum(weights[j]*sp.cancel(dd_product/dd(j+1))**2 for j in range(n))
    second = sum(weights[j]*sp.cancel(dd_product**2/(dd(j)*dd(j+1))) for j in range(n))
    for j, value in ((n-1, first), (n, second)):
        polynomial = sp.Poly(sp.expand(common*value), tau)
        assert [str(polynomial.nth(k)) for k in range(len(saved["coefficient_rows"][j]))] == (
            saved["coefficient_rows"][j])
    return {
        "scope": "finite symbolic identity proof of top completed-row coefficient TP2 for every n>=2,0<p<1; lower rows open",
        "variables": "a=i>=0; overlap j=i+1; separated j=i+2+h,h>=0; u,v formal mixture weights",
        "minor_orientation": "[u*v]([tau^lo]F*[tau^hi]G-[tau^hi]F*[tau^lo]G), lo<hi",
        "mixed_minor_tables": tables,
        "mixed_minors_checked": 16,
        "all_index_symbolic_coefficients_nonnegative": True,
        "proof_completion": "single-index pair reduces to d_i,d_(i+1); common linear factors preserve coefficient minors by Cauchy-Binet",
        "one_formula_audit": {
            "n": n, "p": str(probability), "rank_rows": [n-1, n],
            "reference_key": "rank_likelihood_order",
            "matches_saved_completed_coefficient_rows": True}}


def rank_next_coefficients_certificate():
    """Proved two-positive-support cases and actual-factor low orders only.

    With b_i=u,b_j=v, write the next numerator pair as F=2*C*f,
    G=2*C*g. Each coefficient minor is quartic in u,v. Adjacent and
    gap-two indices have positive reduced tables at every order.
    For separated indices, five u^3*v reduced coefficients are left
    uncontrolled; all other coefficients are certified nonnegative.
    The actual D^3 factor restores those terms through degree six.
    Neither full separated-support TP2 nor full-binomial TP2 is claimed.
    """
    from collections import Counter
    a, h, tau, size = sp.symbols("a h tau size", nonnegative=True)
    dd = lambda k: 1+k*tau
    parameter_poly = lambda value: sp.Poly(value, a, h, domain=sp.QQ)
    zero, one = parameter_poly(0), parameter_poly(1)

    # The complementary two-support case {0,i} has only a linear-in-weight
    # core after its common positive factor 2*i*v is removed.
    u, v = sp.symbols("u v", nonnegative=True)
    zero_support = {}
    for label, i in (("i_one", sp.Integer(1)), ("i_two", sp.Integer(2)),
                     ("i_at_least_three", a+3)):
        cubic = 1+(2*i+3)*tau+(i**2+3*i+3)*tau**2+(2*i+1)*tau**3
        if i == 1:
            first = u*(1+3*tau)+v*dd(1)
            second = u*dd(2)**2+v*dd(1)**2
        elif i == 2:
            first = u*cubic+v*dd(1)**3
            second = dd(2)*(u*dd(3)**2+v*dd(1)**2)
        else:
            first = dd(i-1)*(u*cubic+v*dd(1)**3)
            second = dd(i)*dd(1)*(u*dd(i+1)**2+v*dd(1)**2)
        fp, gp = sp.Poly(first, tau), sp.Poly(second, tau)
        degree = max(fp.degree(), gp.degree())
        table = {}
        for lo in range(degree+1):
            for hi in range(lo+1, degree+1):
                gap = sp.Poly(fp.nth(lo)*gp.nth(hi)-fp.nth(hi)*gp.nth(lo), a, u, v)
                assert all(value >= 0 for value in gap.coeffs())
                if (lo, hi) == (0, 1):
                    assert sp.expand(gap.as_expr()-(u+v)**2) == 0
                table[f"{lo},{hi}"] = str(gap.as_expr())
        assert len(table) == {"i_one": 3, "i_two": 6, "i_at_least_three": 10}[label]
        zero_support[label] = dict(i=str(i), core_F=str(sp.expand(first)),
                                   core_G=str(sp.expand(second)), minors=table)

    def cores(spacing, cancelled_power):
        i, j = a+1, a+1+spacing
        sigma, product = i+j, i*j
        # T_ij=d_i*d_j*d_(i+1)^2*d_(j+1)^2*K_ij, whose identity
        # is already preserved by --rank-top-structure.
        numerator = (
            sigma+tau*(2*sigma**2+3*sigma+2*product)
            +tau**2*(sigma**3+3*sigma**2+6*sigma*product+3*sigma+12*product)
            +tau**3*(4*sigma**2*product+2*sigma**2+18*sigma*product
                     +sigma+4*product**2+10*product)
            +tau**4*(3*sigma**2*product+5*sigma*product**2+9*sigma*product
                     +12*product**2+4*product)
            +tau**5*(sigma**2*product+3*sigma*product**2+sigma*product
                     +2*product**3+2*product**2))
        first = [i*dd(i-1)*dd(j-1)*dd(j)**2*dd(j+1)**3,
                 dd(i-1)*dd(j-1)*numerator,
                 j*dd(j-1)*dd(i-1)*dd(i)**2*dd(i+1)**3]
        second = [i*dd(i)*dd(j-1)*dd(j)**2*dd(j+1)**3,
                  i*dd(i)*dd(i+1)**2*dd(j-1)*dd(j)**2*dd(j+1)
                  +j*dd(j)*dd(j+1)**2*dd(i-1)*dd(i)**2*dd(i+1),
                  j*dd(j)*dd(i-1)*dd(i)**2*dd(i+1)**3]
        divisor = sp.Poly(dd(i+1)**cancelled_power, tau, a, h)
        reduce = lambda value: sp.Poly(
            sp.Poly(value, tau, a, h).exquo(divisor).as_expr(), tau)
        return list(map(reduce, first)), list(map(reduce, second))

    uncontrolled = {(3, 4), (3, 5), (4, 5), (4, 6), (5, 6)}
    cases = {}
    separated_first = separated_second = None
    for label, spacing, cancelled_power in (
            ("adjacent", sp.Integer(1), 2),
            ("gap_two", sp.Integer(2), 1),
            ("separated", h+3, 0)):
        first, second = cores(spacing, cancelled_power)
        degree = max(value.degree() for value in first+second)
        fc = [[parameter_poly(value.nth(k)) for k in range(degree+1)] for value in first]
        gc = [[parameter_poly(value.nth(k)) for k in range(degree+1)] for value in second]
        entries, checked = {}, 0
        for lo in range(degree+1):
            for hi in range(lo+1, degree+1):
                coefficients = []
                for v_power in range(5):
                    if label == "separated" and v_power == 1 and (lo, hi) in uncontrolled:
                        # These reduced signs are not a proved positivity claim.
                        # The low-order restoration below uses the full factor.
                        coefficients.append(None)
                        continue
                    gap = sum((fc[r][lo]*gc[v_power-r][hi]-fc[r][hi]*gc[v_power-r][lo]
                               for r in range(max(0, v_power-2), min(2, v_power)+1)), zero)
                    assert all(value >= 0 for value in gap.coeffs())
                    coefficients.append(str(gap.as_expr()))
                    checked += 1
                entries[f"{lo},{hi}"] = coefficients
        assert checked == {"adjacent": 75, "gap_two": 105, "separated": 135}[label]
        cases[label] = {
            "i": "a+1", "j": str(a+1+spacing), "core_degree": degree,
            "cancelled_d_i_plus_one_power": cancelled_power,
            "F_u2_uv_v2": [str(value.as_expr()) for value in first],
            "G_u2_uv_v2": [str(value.as_expr()) for value in second],
            "minor_coefficients_u4_u3v_u2v2_uv3_v4": entries,
            "nonnegative_coefficient_polynomials_checked": checked}
        if label == "separated":
            separated_first, separated_second = first[:2], second[:2]

    # For minimal n=j+1, compute the actual multiplier C through degree
    # six by Newton identities. Extra ambient factors preserve this prefix.
    i, j, n, limit = a+1, a+h+4, a+h+5, 6
    powers = [None]
    for r in range(1, limit+1):
        faulhaber = (sp.bernoulli(r+1, size+1)-sp.bernoulli(r+1, 1))/(r+1)
        powers.append(parameter_poly(
            3*faulhaber.subs(size, n)-(i-1)**r-2*i**r-3*(i+1)**r
            -(j-1)**r-2*j**r-3*(j+1)**r))
    multiplier = [one]
    for k in range(1, limit+1):
        multiplier.append(sum(((-1)**(r-1)*multiplier[k-r]*powers[r]
                               for r in range(1, k+1)), zero).mul_ground(sp.Rational(1, k)))

    def restore(value):
        return [sum((multiplier[r]*parameter_poly(value.nth(k-r))
                     for r in range(k+1)), zero) for k in range(limit+1)]

    first = list(map(restore, separated_first))
    second = list(map(restore, separated_second))
    restored = {}
    for lo in range(limit+1):
        for hi in range(lo+1, limit+1):
            gap = (first[0][lo]*second[1][hi]+first[1][lo]*second[0][hi]
                   -first[0][hi]*second[1][lo]-first[1][hi]*second[0][lo])
            assert all(value >= 0 for value in gap.coeffs())
            if (lo, hi) == (0, 1):
                assert gap == parameter_poly(2*i*(i+j))
            restored[f"{lo},{hi}"] = str(gap.as_expr())
    assert len(restored) == 21

    # One existing-witness audit of the power sums against actual factors.
    # This does not repeat the full restored-minor or fixed-size LR checks.
    removed = Counter({1: 2, 2: 3, 67: 1, 68: 2, 69: 3})
    direct = sp.Poly(1, tau)
    for k in range(1, 70):
        direct *= sp.Poly(dd(k), tau)**(3-removed[k])
    assert all(value.eval({a: 0, h: 64}) == direct.nth(k)
               for k, value in enumerate(multiplier))
    return {
        "scope": "proved two-positive-support coefficient cases and actual-factor restoration through degree six only; full next-row and full-binomial TP2 remain open",
        "variables": "a=i-1>=0; separated j=i+3+h,h>=0; tau formal; b_i=u,b_j=v, all other weights zero",
        "normalization": "F=2*C*f,G=2*C*g; stored restored minors use C*f,C*g, so actual F,G minors are four times these",
        "zero_and_positive_support": zero_support,
        "table_convention": "minor row lo,hi lists coefficients of u^4,u^3v,u^2v^2,uv^3,v^4; null denotes an uncontrolled reduced coefficient",
        "reduced_cases": cases,
        "reduced_nonnegative_polynomials_checked": 315,
        "uncontrolled_reduced_u3v_degree_pairs": [list(pair) for pair in sorted(uncontrolled)],
        "order_lemma": "all restored minors with lo<=2 are nonnegative for every hi>lo; the five uncontrolled inputs have first degree at least three",
        "low_degree_restoration": {
            "minimal_n": "a+h+5", "max_coefficient_degree": limit,
            "common_multiplier": "D^3/[d_(i-1)*d_i^2*d_(i+1)^3*d_(j-1)*d_j^2*d_(j+1)^3]",
            "multiplier_coefficients_0_to_6": [str(value.as_expr()) for value in multiplier],
            "restored_u3v_minors": restored, "nonnegative_polynomials_checked": 21,
            "extension": "extra ambient factors and the remaining killed continuants preserve this prefix by Cauchy-Binet"},
        "same_witness_multiplier_audit": {
            "i": 1, "j": 68, "n": 69, "degrees_checked": list(range(limit+1)),
            "power_sums_match_actual_factor_product": True},
        "remaining_gap": "separated two-support lo>=3,hi>=7 and general three/four-distinct-weight quartic terms"}


def rank_next_pointwise_certificate():
    """Strict next-row pointwise LR for any two positive support indices.

    Write f=(V'W-2VW')/2, g=U*V, and Delta=g'*f-g*f'. For weights
    b_i=u,b_j=v, the common denominator E makes E*f and E*g quadratic
    polynomial mixtures. Every quartic weight coefficient of E^2*Delta
    has positive coefficients after i=a+1,j=a+h+2, with a,h>=0.
    This covers every positive index pair and every tau>=0, not general
    three/four-index supports or coefficient TP2. No old check is rerun.
    """
    from hashlib import sha256
    from sympy.polys.rings import ring

    parameters, a, h = ring("a,h", sp.QQ)
    polynomials, tau = ring("tau", parameters)
    i, j = a+1, a+h+2
    dd = lambda k: 1+tau*k
    sigma, product = i+j, i*j
    numerator = (
        polynomials(sigma)+tau*(2*sigma**2+3*sigma+2*product)
        +tau**2*(sigma**3+3*sigma**2+6*sigma*product+3*sigma+12*product)
        +tau**3*(4*sigma**2*product+2*sigma**2+18*sigma*product
                 +sigma+4*product**2+10*product)
        +tau**4*(3*sigma**2*product+5*sigma*product**2+9*sigma*product
                 +12*product**2+4*product)
        +tau**5*(sigma**2*product+3*sigma*product**2+sigma*product
                 +2*product**3+2*product**2))
    first = [dd(i-1)*dd(j-1)*dd(j)**2*dd(j+1)**3*i,
             dd(i-1)*dd(j-1)*numerator,
             dd(j-1)*dd(i-1)*dd(i)**2*dd(i+1)**3*j]
    second = [dd(i)*dd(j-1)*dd(j)**2*dd(j+1)**3*i,
              dd(i)*dd(i+1)**2*dd(j-1)*dd(j)**2*dd(j+1)*i
              +dd(j)*dd(j+1)**2*dd(i-1)*dd(i)**2*dd(i+1)*j,
              dd(j)*dd(i-1)*dd(i)**2*dd(i+1)**3*j]

    # Quadratic identity g-f=tau*V*U1+tau^2*Q, with U1=sum b_i*u_i/d_i
    # and Q=sum_{i<j} b_i*b_j*(j-i)^2/(d_i*d_j*d_(i+1)^3*d_(j+1)^3).
    assert second[0]-first[0] == tau*dd(j-1)*dd(j)**2*dd(j+1)**3*i
    assert second[2]-first[2] == tau*dd(i-1)*dd(i)**2*dd(i+1)**3*j
    mixed_loss = (tau*dd(i+1)**2*dd(j-1)*dd(j)**2*dd(j+1)*i
                  +tau*dd(j+1)**2*dd(i-1)*dd(i)**2*dd(i+1)*j
                  +tau**2*dd(i-1)*dd(i)*dd(j-1)*dd(j)*(j-i)**2)
    assert second[1]-first[1] == mixed_loss

    expected_counts = [637, 647, 620, 556, 303]
    summaries, term_tables = [], []
    for power in range(5):
        wronskian = sum((first[z]*second[power-z].diff(tau)
                         -second[power-z]*first[z].diff(tau)
                         for z in range(max(0, power-2), min(2, power)+1)),
                        polynomials.zero)
        terms = []
        for (tau_power,), coefficient_poly in wronskian.items():
            for (a_power, h_power), coefficient in coefficient_poly.items():
                assert coefficient.denominator == 1 and coefficient > 0
                terms.append([int(tau_power), int(a_power), int(h_power), int(coefficient)])
        terms.sort()
        assert wronskian.degree() == 12 and len(terms) == expected_counts[power]
        if power == 0:
            assert wronskian == (dd(j-1)*dd(j)**2*dd(j+1)**3*i)**2
        elif power == 4:
            assert wronskian == (dd(i-1)*dd(i)**2*dd(i+1)**3*j)**2
        summaries.append({
            "v_power": power, "tau_degree": int(wronskian.degree()),
            "positive_integer_coefficients": len(terms), "negative_coefficients": 0,
            "minimum_coefficient": min(term[3] for term in terms),
            "coefficient_sha256": sha256(
                json.dumps(terms, separators=(",", ":")).encode("ascii")).hexdigest()})
        term_tables.append(terms)

    # One independent V,W,U audit, reusing the separated support and its
    # natural scale from the common-core obstruction. No size/rate sweep.
    t = sp.Symbol("tau")
    ii, jj = 1, 134
    d = lambda k: 1+k*t
    vv = d(ii+1)**-2+d(jj+1)**-2
    ww = 1/(d(ii)*d(ii+1))+1/(d(jj)*d(jj+1))
    uu = ii/(d(ii-1)*d(ii)*d(ii+1))+jj/(d(jj-1)*d(jj)*d(jj+1))
    ff, gg = (sp.diff(vv, t)*ww-2*vv*sp.diff(ww, t))/2, uu*vv
    at = sp.Rational(1, jj)
    direct = sp.factor((sp.diff(gg, t)*ff-gg*sp.diff(ff, t)).subs(t, at))
    denominator = d(ii-1)*d(ii)**2*d(ii+1)**3*d(jj-1)*d(jj)**2*d(jj+1)**3
    # Here a=0,h=j-i-1=132 and u=v=1.
    from_core = sum(sp.Integer(coefficient)*at**tau_power*sp.Integer(132)**h_power
                    for terms in term_tables
                    for tau_power, a_power, h_power, coefficient in terms if a_power == 0)
    assert direct > 0
    assert sp.factor(from_core/denominator.subs(t, at)**2-direct) == 0

    return {
        "scope": "strict pointwise next-row LR for arbitrary two positive support indices; general three/four-index supports and next-row coefficient TP2 remain open",
        "variables": "i=a+1,j=a+h+2; a,h>=0; tau>=0; b_i=u,b_j=v>=0,u+v>0; all other weights zero",
        "normalization": "f=(V'*W-2*V*W')/2,g=U*V; original F=2*D^3*f,G=2*D^3*g; G'*F-G*F'=4*D^6*Delta",
        "wronskian": "Delta=g'*f-g*f'; E^2*Delta=sum_{r=0}^4 u^(4-r)*v^r*Z_r; Z_r=sum_{a+b=r}(F_a*G_b'-G_b*F_a')",
        "common_denominator": "E=d_(i-1)*d_i^2*d_(i+1)^3*d_(j-1)*d_j^2*d_(j+1)^3; d_k=1+k*tau",
        "core_F_u2_uv_v2": [
            "i*d_(i-1)*d_(j-1)*d_j^2*d_(j+1)^3",
            "d_(i-1)*d_(j-1)*T_ij",
            "j*d_(j-1)*d_(i-1)*d_i^2*d_(i+1)^3"],
        "core_G_u2_uv_v2": [
            "i*d_i*d_(j-1)*d_j^2*d_(j+1)^3",
            "i*d_i*d_(i+1)^2*d_(j-1)*d_j^2*d_(j+1)+j*d_j*d_(j+1)^2*d_(i-1)*d_i^2*d_(i+1)",
            "j*d_j*d_(i-1)*d_i^2*d_(i+1)^3"],
        "T_ij_variables": "sigma=i+j,rho=i*j",
        "T_ij_coefficients_tau_0_to_5": [
            "sigma", "2*sigma^2+3*sigma+2*rho",
            "sigma^3+3*sigma^2+6*sigma*rho+3*sigma+12*rho",
            "4*sigma^2*rho+2*sigma^2+18*sigma*rho+sigma+4*rho^2+10*rho",
            "3*sigma^2*rho+5*sigma*rho^2+9*sigma*rho+12*rho^2+4*rho",
            "sigma^2*rho+3*sigma*rho^2+sigma*rho+2*rho^3+2*rho^2"],
        "loss_identity": "g-f=tau*V*U1+tau^2*Q; U1=sum b_i*u_i/d_i; Q=sum_{i<j} b_i*b_j*(j-i)^2/(d_i*d_j*d_(i+1)^3*d_(j+1)^3)",
        "loss_identity_verified": True,
        "hash_encoding": "SHA256 of ASCII compact JSON (separators comma,colon), lexicographically sorted [tau_power,a_power,h_power,integer_coefficient] rows; no trailing newline",
        "quartic_coefficients": summaries,
        "positive_integer_coefficients_checked": sum(expected_counts),
        "strictness": "Z_0=i^2*(d_(j-1)*d_j^2*d_(j+1)^3)^2; Z_4=j^2*(d_(i-1)*d_i^2*d_(i+1)^3)^2; both identities verified",
        "single_direct_audit": {
            "i": ii, "j": jj, "weights": [1, 1], "tau": str(at),
            "direct_wronskian": str(direct), "strictly_positive": True,
            "matches_core_divided_by_E_squared": True}}


def rank_next_pointwise_all_certificate():
    """Three/four-index quartics closing arbitrary-weight next-row LR.

    With h=g-f, certify the coefficients of Delta=g*h'-g'*h for all
    distinct index triples and quadruples, separately including zero.
    Positive index gaps remain symbolic. The one/two-index certificates
    are prior inputs and are not rerun. Only this helper requires flint.
    """
    from hashlib import sha256
    from math import prod
    try:
        from flint import fmpz_mpoly_ctx
    except ImportError as error:
        raise SystemExit("--rank-next-pointwise-all requires optional python-flint") from error

    ctx = fmpz_mpoly_ctx.get(["t", "a", "b", "c", "d"])
    t, a, b, c, gap_d = ctx.gens()
    zero, one = ctx.constant(0), ctx.constant(1)
    dd = lambda r: 1+t*r
    ee = lambda r: dd(r-1)*dd(r)**2*dd(r+1)**3

    def positive_pair(r, s):
        # B_rs=e_r*e_s*(g_rs,h_rs), with unordered quadratic coefficients.
        first = (r*dd(r)*dd(r+1)**2*dd(s-1)*dd(s)**2*dd(s+1)
                 +s*dd(s)*dd(s+1)**2*dd(r-1)*dd(r)**2*dd(r+1))
        second = (t*r*dd(r+1)**2*dd(s-1)*dd(s)**2*dd(s+1)
                  +t*s*dd(s+1)**2*dd(r-1)*dd(r)**2*dd(r+1)
                  +t**2*(s-r)**2*dd(r-1)*dd(r)*dd(s-1)*dd(s))
        return first, second

    def zero_pair(r):
        # e_0=d_1^3, u_0=0: never introduce d_(-1).
        return (r*dd(1)*dd(r)*dd(r+1)**2,
                t*r*dd(1)*dd(r+1)**2+t**2*r**2*dd(r-1)*dd(r))

    def scale(pair, factor):
        return pair[0]*factor, pair[1]*factor

    def polarized(first, second):
        G, H = first
        A, B = second
        return (G*B.derivative("t")+A*H.derivative("t")
                -G.derivative("t")*B-A.derivative("t")*H)

    def summarize(polynomial, variables, count, degree, tau_zero, constant):
        data = polynomial.to_dict()
        assert all(coefficient > 0 for coefficient in data.values())
        assert all(not any(powers[variables:]) for powers in data)
        terms = sorted([*map(int, powers[:variables]), int(coefficient)]
                       for powers, coefficient in data.items())
        degrees = list(map(int, polynomial.degrees()))[:variables]
        assert len(terms) == count and degrees[0] == degree
        assert ctx.from_dict({powers: value for powers, value in data.items()
                              if powers[0] == 0}) == tau_zero
        assert data[(0, 0, 0, 0, 0)] == constant and constant > 0
        return {
            "degrees": degrees, "positive_terms": count, "negative_terms": 0,
            "minimum_coefficient": min(row[-1] for row in terms),
            "constant_coefficient": constant, "tau_zero_identity_verified": True,
            "sha256": sha256(json.dumps(terms, separators=(",", ":")).encode("ascii")).hexdigest()
        }, terms

    def at_half(terms, degrees, parameters):
        # Integer homogenization evaluates t=1/2 without floating point.
        powers = [[base**exponent for exponent in range(degrees[position+1]+1)]
                  for position, base in enumerate(parameters)]
        scaled = sum(row[-1]*2**(degrees[0]-row[0])
                     *prod(powers[q][power] for q, power in enumerate(row[1:-1]))
                     for row in terms)
        return sp.Rational(scaled, 2**degrees[0])

    def direct_audit(support, numerators, monomials):
        """Differentiate V,W,U directly; do not use pair kernels or S."""
        at = sp.Rational(1, 2)
        weights = sp.symbols(f"b0:{len(support)}")
        V = Vp = Vpp = W = Wp = Wpp = U = Up = sp.S.Zero
        E = R = sp.S.One
        d = lambda r: 1+r*at
        for r, weight in zip(support, weights):
            v, vp, vpp = d(r+1)**-2, -2*(r+1)*d(r+1)**-3, 6*(r+1)**2*d(r+1)**-4
            w = 1/(d(r)*d(r+1))
            rate = r/d(r)+(r+1)/d(r+1)
            wp = -w*rate
            wpp = w*(rate**2+r**2/d(r)**2+(r+1)**2/d(r+1)**2)
            V += weight*v
            Vp += weight*vp
            Vpp += weight*vpp
            W += weight*w
            Wp += weight*wp
            Wpp += weight*wpp
            if r:
                u = r/(d(r-1)*d(r)*d(r+1))
                up = -u*((r-1)/d(r-1)+r/d(r)+(r+1)/d(r+1))
                U += weight*u
                Up += weight*up
                E *= d(r-1)*d(r)**2*d(r+1)**3
                R *= d(r-1)*d(r)*d(r+1)
        if 0 in support:
            E *= d(1)**3
            R *= d(1)
        f, fp = (Vp*W-2*V*Wp)/2, (Vpp*W-Vp*Wp-2*V*Wpp)/2
        g, gp = U*V, Up*V+U*Vp
        direct = sp.Poly(gp*f-g*fp, *weights)
        divisor = E**2 if len(support) == 3 else E*R
        values = []
        for numerator, monomial in zip(numerators, monomials):
            actual = direct.coeff_monomial(tuple(monomial))
            assert actual > 0 and actual == numerator/divisor
            values.append(str(actual))
        result = {"support": support, "tau": str(at)}
        if len(support) == 3:
            result.update(three_repeated_index_coefficients=values,
                          matches_direct_V_W_U_derivatives=True)
        else:
            result.update(direct_squarefree_quartic_coefficient=values[0],
                          matches_N_divided_by_E_R=True)
        return result

    three = {}
    repeated_patterns = [(0, 1, 2), (1, 0, 2), (2, 0, 1)]
    for label in ("positive", "zero"):
        if label == "positive":
            indices = [a+1, a+b+2, a+b+c+3]
            factors = [ee(r) for r in indices]
            pairs = {}
            for p, r in enumerate(indices):
                complement = factors[(p+1) % 3]*factors[(p+2) % 3]
                pairs[p, p] = scale((r*dd(r), t*r), complement)
                for q in range(p+1, 3):
                    pairs[p, q] = scale(positive_pair(r, indices[q]), factors[3-p-q])
            substitutions = ["i=a+1", "j=a+b+2", "k=a+b+c+3"]
            denominator = "E=e_i*e_j*e_k"
            counts, constants, degree = [22899, 22195, 20335], [34, 46, 58], 24
            support = [1, 2, 3]
        else:
            j, k = a+1, a+b+2
            indices = [zero, j, k]
            factors = [dd(1)**3, ee(j), ee(k)]
            pairs = {(0, 0): (zero, zero)}
            for p, r in ((1, j), (2, k)):
                pairs[p, p] = scale((r*dd(r), t*r), factors[0]*factors[3-p])
                pairs[0, p] = scale(zero_pair(r), factors[3-p])
            pairs[1, 2] = scale(positive_pair(j, k), factors[0])
            substitutions = ["i=0", "j=a+1", "k=a+b+2"]
            denominator = "E=d_1^3*e_j*e_k"
            counts, constants, degree = [1537, 1505, 1376], [4, 10, 16], 18
            support = [0, 1, 2]
        summaries, audit_values, monomials = [], [], []
        for repeated, s, u in repeated_patterns:
            block = (polarized(pairs[repeated, repeated], pairs[s, u])
                     +polarized(pairs[tuple(sorted((repeated, s)))],
                                pairs[tuple(sorted((repeated, u)))]))
            r_index, s_index, u_index = indices[repeated], indices[s], indices[u]
            tau_zero = 2*(r_index*(s_index+u_index)+(r_index+s_index)*(r_index+u_index))
            summary, terms = summarize(block, 4, counts[repeated], degree,
                                       tau_zero, constants[repeated])
            summary["repeated_position"] = repeated
            summary["tau_zero_identity"] = "2*(r*(s+u)+(r+s)*(r+u))"
            summaries.append(summary)
            audit_values.append(at_half(terms, summary["degrees"], [0, 0, 0]))
            monomial = [1, 1, 1]
            monomial[repeated] = 2
            monomials.append(monomial)
        three[label] = {"substitutions": substitutions, "denominator": denominator,
                        "patterns": summaries,
                        "one_direct_audit": direct_audit(support, audit_values, monomials)}

    four = {}
    partitions = [((0, 1), (2, 3)), ((0, 2), (1, 3)), ((0, 3), (1, 2))]
    for label in ("positive", "zero"):
        if label == "positive":
            indices = [a+1, a+b+2, a+b+c+3, a+b+c+gap_d+4]
            substitutions = ["i=a+1", "j=a+b+2", "k=a+b+c+3", "l=a+b+c+d+4"]
            count, degree, constant = 122110, 24, 140
            support, parameters = [2, 5, 8, 11], [1, 2, 2, 2]
        else:
            indices = [zero, a+1, a+b+2, a+b+c+3]
            substitutions = ["i=0", "j=a+1", "k=a+b+2", "l=a+b+c+3"]
            count, degree, constant = 9587, 19, 44
            support, parameters = [0, 3, 6, 9], [2, 2, 2, 0]
        # R contains each formal pole occurrence; overlaps of shifted
        # factors are retained. The positive factor E/R is cancelled.
        factors, groups = [], []
        for p, r in enumerate(indices):
            group = []
            roots = [(1, 3)] if label == "zero" and p == 0 else [(r-1, 1), (r, 2), (r+1, 3)]
            for root, multiplicity in roots:
                group.append((len(factors), root, multiplicity))
                factors.append(dd(root))
            groups.append(group)
        R = prod(factors, start=one)
        quotients = [prod((factor for q, factor in enumerate(factors) if q != p), start=one)
                     for p in range(len(factors))]
        logs = [sum((quotients[position]*root*multiplicity
                     for position, root, multiplicity in group), zero) for group in groups]
        bases = {(p, q): (zero_pair(indices[q]) if label == "zero" and p == 0
                         else positive_pair(indices[p], indices[q]))
                 for p in range(4) for q in range(p+1, 4)}
        ordinary, correction = zero, zero
        for left, right in partitions:
            first, second = bases[left], bases[right]
            ordinary += polarized(first, second)
            signed_log = logs[left[0]]+logs[left[1]]-logs[right[0]]-logs[right[1]]
            correction += signed_log*(first[0]*second[1]-second[0]*first[1])
        numerator = R*ordinary+correction
        tau_zero = 4*sum((indices[p]*indices[q] for p in range(4) for q in range(p+1, 4)), zero)
        summary, terms = summarize(numerator, 5, count, degree, tau_zero, constant)
        value = at_half(terms, summary["degrees"], parameters)
        audit = direct_audit(support, [value], [(1, 1, 1, 1)])
        audit["index_parameter_values"] = parameters
        summary.update(substitutions=substitutions, one_direct_audit=audit,
                       tau_zero_identity="4*(i*j+i*k+i*l+j*k+j*l+k*l)")
        four[label] = summary

    three_count = sum(pattern["positive_terms"] for case in three.values() for pattern in case["patterns"])
    four_count = sum(case["positive_terms"] for case in four.values())
    assert (three_count, four_count) == (69847, 131697)
    return {
        "scope": "all-index three/four-distinct-weight quartic coefficient proof; with prior one/two-index inputs, arbitrary finite nonnegative weights have strict next-row pointwise LR whenever U>0",
        "not_claimed": "coefficient TP2, all other adjacent rank ratios, or full ULC/NA",
        "optional_dependency": "python-flint, imported locally by this helper only",
        "prior_inputs_not_rerun": [
            {"certificate_key": "rank_next_pointwise", "function": "rank_next_pointwise_certificate", "scope": "one/two positive support indices"},
            {"certificate_key": "rank_next_coefficients", "function": "rank_next_coefficients_certificate", "scope": "zero_and_positive_support cases {0,i}"}],
        "normalization": "f=(V'*W-2*V*W')/2,g=U*V,h=g-f,Delta=g'*f-g*f'=g*h'-g'*h; original F=2*D^3*f,G=2*D^3*g",
        "quadratic_kernels": {
            "definitions": "d_r=1+r*t; e_r=d_(r-1)*d_r^2*d_(r+1)^3 for r>=1; e_0=d_1^3; B_rs=e_r*e_s*(g_rs,h_rs)",
            "positive_pair_G": "r*d_r*d_(r+1)^2*d_(s-1)*d_s^2*d_(s+1)+s*d_s*d_(s+1)^2*d_(r-1)*d_r^2*d_(r+1)",
            "positive_pair_H": "t*r*d_(r+1)^2*d_(s-1)*d_s^2*d_(s+1)+t*s*d_(s+1)^2*d_(r-1)*d_r^2*d_(r+1)+t^2*(s-r)^2*d_(r-1)*d_r*d_(s-1)*d_s",
            "zero_pair": "B_0r=(r*d_1*d_r*d_(r+1)^2,t*r*d_1*d_(r+1)^2+t^2*r^2*d_(r-1)*d_r)",
            "diagonal": "e_r*(g_rr,h_rr)=(r*d_r,r*t); (g_00,h_00)=(0,0)",
            "polarization": "S((G,H),(A,B))=G*B'+A*H'-G'*B-A'*H"},
        "hash_encoding": "SHA256 of ASCII compact JSON (separators comma,colon), sorted exponent rows followed by integer coefficient; three-index rows use t,a,b,c, four-index rows use t,a,b,c,d; no trailing newline",
        "three_index": {
            "variables": ["t", "a", "b", "c"], "nonnegative_index_parameters": True,
            "normalization": "E^2*[b_r^2*b_s*b_u]Delta=S(P_rr,P_su)+S(P_rs,P_ru); P_rs=E*(g_rs,h_rs)",
            "cases": three, "positive_integer_coefficients_checked": three_count},
        "four_index": {
            "variables": ["t", "a", "b", "c", "d"], "nonnegative_index_parameters": True,
            "normalization": "[b_i*b_j*b_k*b_l]Delta=N/(E*R); E=product e_r, R=product rho_r; rho_r=d_(r-1)*d_r*d_(r+1),rho_0=d_1",
            "reduced_numerator": "N=R*sum_pairings S(B_rs,B_uv)+sum_pairings R*(ell_r+ell_s-ell_u-ell_v)*det(B_rs,B_uv); ell_r=e_r'/e_r",
            "cancellation": "sum_pairings S(P_rs,P_uv)=(E/R)*N; E/R=product_(r>=1) d_r*d_(r+1)^2, with additional d_1^2 when index zero occurs",
            "cases": four, "positive_integer_coefficients_checked": four_count},
        "positive_integer_coefficients_checked": three_count+four_count,
        "proof_completion": "Delta is homogeneous quartic in the weights, so the prior one/two-index signs and these three/four-index signs exhaust every monomial; all coefficients are nonnegative for each t>=0",
        "strictness": "Delta>=sum_(r>=1) b_r^4*r^2/(d_(r-1)^2*d_r^4*d_(r+1)^6)>0 whenever U>0; (f/g)'=-Delta/g^2<0"}


def rank_likelihood_order_certificate():
    """Strict stationary LR order at fixed n=5,p=1/2, every tau>0.

    Put P_j=H*r_j/pi_j, H=product of the leading killed determinants.
    If a[j,k] is the coefficient of tau^k in P_j, then
    P'_j P_(j-1)-P_j P'_(j-1) is the sum over lo<hi of
    (hi-lo)*(a[j-1,lo]*a[j,hi]-a[j-1,hi]*a[j,lo])*tau^(hi+lo-1).
    We check every such minor exactly. Nonnegative coefficients of P_j
    alone, proved separately for every size, do not imply these minors.
    """
    n, probability = 5, sp.Rational(1, 2)
    tau = sp.Symbol("tau")
    generator = sp.zeros(n+1)
    for j in range(n+1):
        birth, death = probability*(n-j), (1-probability)*j
        generator[j, j] = -birth-death
        if j < n:
            generator[j, j+1] = birth
        if j:
            generator[j, j-1] = death
    matrix = sp.eye(n+1)-tau*generator.T
    transform = sp.Matrix.hstack(*[(matrix**(n-j))[:, j] for j in range(n+1)])
    initial = sp.Matrix([sp.binomial(n, j)*probability**j*(1-probability)**(n-j)
                         for j in range(n+1)])
    common = sp.Poly(sp.prod(matrix[:k, :k].det() for k in range(1, n+1)), tau)
    rank = transform.inv(method="DM")*initial
    polynomials = [sp.Poly(sp.cancel(common.as_expr()*rank[j]/initial[j]), tau)
                   for j in range(n+1)]
    numerators = sp.Matrix([initial[j]*polynomials[j].as_expr() for j in range(n+1)])
    assert all(sp.expand(value) == 0
               for value in transform*numerators-common.as_expr()*initial)
    assert sp.expand(sum(numerators)-common.as_expr()) == 0
    assert all(value > 0 for polynomial in polynomials for value in polynomial.all_coeffs())

    degree = common.degree()
    coefficients = [[polynomial.nth(k) for k in range(degree+1)]
                    for polynomial in polynomials]
    checked, zero, adjacent = 0, 0, []
    for j in range(n):
        row = []
        for lo in range(degree+1):
            for hi in range(lo+1, degree+1):
                gap = (coefficients[j][lo]*coefficients[j+1][hi]
                       -coefficients[j][hi]*coefficients[j+1][lo])
                assert gap >= 0
                checked += 1
                zero += int(gap == 0)
                if hi == lo+1:
                    row.append(str(gap))
        assert (coefficients[j][0]*coefficients[j+1][1]
                -coefficients[j][1]*coefficients[j+1][0]) == 1
        adjacent.append(row)
    return {
        "scope": "strict stationary LR order at fixed n=5,p=1/2, every tau>0; no all-n or all-p claim",
        "n": n, "p": str(probability), "parameter": "tau=z/[n(1-z)]",
        "criterion": "nonnegative adjacent-row coefficient minors imply increasing r_j/r_(j-1)",
        "general_status": {
            "polynomial_coefficient_nonnegativity": "proved analytically for every n and 0<p<1, separately from this check",
            "coefficient_cross_minors": "open for general n,p; coefficient nonnegativity alone does not suffice"},
        "common_denominator": str(common.as_expr()),
        "normalization": "P_j=common_denominator*r_j/pi_j",
        "coefficient_order": "ascending powers of tau, padded to common denominator degree",
        "coefficient_rows": [[str(value) for value in row] for row in coefficients],
        "degrees": [polynomial.degree() for polynomial in polynomials],
        "adjacent_column_minors": adjacent,
        "adjacent_row_all_column_pairs_checked": checked, "zero_minors": zero,
        "all_checked_minors_nonnegative": True,
        "constant_Wronskian_coefficient_each_row": 1,
        "defining_polynomial_system_residual_zero": True,
        "normalization_residual_zero": True}


def reinjection_logcurvature_symbolic(n):
    """All-positive-parameter ULC-ratio monotonicity at a specified size.

    x is the base sleep rate and y=z/(1-z) is the feedback odds. The
    coefficient test concerns the actual stationary trajectory, not a
    generic ULC-preserving vector field. An independent symmetric IDLA
    residual checks the symbolic law over the whole parameter domain.
    """
    x, y = sp.symbols("x y", positive=True)
    scaled = (1+x)*sp.eye(n+1)
    for j in range(n+1):
        if j < n:
            value = x*(n-j)*y/n
            scaled[j+1, j] -= value
            scaled[j, j] += value
        if j:
            value = sp.Integer(j)*y/n
            scaled[j-1, j] -= value
            scaled[j, j] += value
    powers = [sp.eye(n+1)]
    for _ in range(n):
        powers.append((scaled*powers[-1]).applyfunc(sp.expand))
    matrix = sp.Matrix.hstack(*[powers[n-j][:, j] for j in range(n+1)])
    target = sp.Matrix([sp.binomial(n, j)*x**j for j in range(n+1)])
    solution = matrix.inv(method="DM")*target
    rank = sp.Matrix([sp.factor(solution[j]/(1+x)**j) for j in range(n+1)])
    assert sp.cancel(sum(rank)-1) == 0
    initial = sp.Matrix([sp.binomial(n, j)*x**j/(1+x)**n for j in range(n+1)])
    assert (rank.subs(y, 0)-initial).applyfunc(sp.cancel) == sp.zeros(n+1, 1)

    # Independent IDLA algebra: after an original jump, absorption at the
    # sink has weight n, while each vacant vertex has weight y.
    jump = sp.zeros(n+1)
    for j in range(n+1):
        jump[j, j] = n/(n+(n-j)*y)
        if j < n:
            jump[j+1, j] = (n-j)*y/(n+(n-j)*y)
    jump_powers = [sp.eye(n+1)]
    for _ in range(n):
        jump_powers.append((jump*jump_powers[-1]).applyfunc(sp.cancel))
    idla_basis = sp.Matrix.hstack(*[
        sum((sp.binomial(k, ell)*(1+1/x)**ell*(-1/x)**(k-ell)
             *jump_powers[k-ell][:, ell] for ell in range(k+1)),
            sp.zeros(n+1, 1)).applyfunc(sp.factor)
        for k in range(n+1)])
    full = sp.zeros(n+1, 1)
    full[n] = 1
    assert (idla_basis*rank-full).applyfunc(sp.cancel) == sp.zeros(n+1, 1)

    def positive_coefficients(expression):
        numerator, denominator = sp.fraction(sp.factor(expression))
        polynomials = [sp.Poly(part, x, y) for part in (numerator, denominator)]
        assert all(all(coefficient >= 0 for coefficient in poly.coeffs())
                   and not poly.is_zero for poly in polynomials)
        return numerator, denominator, polynomials

    for value in rank:
        positive_coefficients(value)
    gaps = []
    for j in range(1, n):
        ratio = sp.factor(rank[j]**2/(rank[j-1]*rank[j+1]))
        slope = sp.factor(sp.diff(ratio, y)/ratio)
        numerator, denominator, polynomials = positive_coefficients(slope)
        gaps.append({
            "rank": j, "log_ratio_derivative_numerator": str(numerator),
            "log_ratio_derivative_denominator": str(denominator),
            "positive_coefficient_test": True,
            "polynomial_term_counts": [len(poly.terms()) for poly in polynomials]})
        print(n, j, "stationary log-ULC ratio increases for all x,y>0", flush=True)
    return {"n": n, "parameters": "x=lambda_0>0; y=z/(1-z)>0",
            "rank_law": list(map(str, rank)), "gaps": gaps,
            "binomial_initial_condition": True,
            "independent_symbolic_IDLA_residual": True}


def reinjection_logcurvature_certificate():
    """New invariant at n=3,4, plus a failed sufficient curvature condition."""
    symbolic = [reinjection_logcurvature_symbolic(n) for n in (3, 4)]
    n, probability, feedback = 3, sp.Rational(2, 5), sp.Rational(4, 9)
    rank, addition, derivative = reinjection_rank(n, probability, feedback)
    kernel = tuple(tuple([(j, feedback/n) for j in range(n)]
                         +[(n, 1-feedback)]) for _ in range(n))
    raw, raw_counts = stationary_law(kernel, (probability,)*n)
    raw_rank = [sum(raw[mask] for mask in range(1 << n) if mask.bit_count() == j)
                for j in range(n+1)]
    assert raw_rank == list(rank)
    starts = [add(stable((1 << j)-1, n), v) for j in range(n+1) for v in range(n)]
    forward, nf = absorption(starts, kernel, (probability,)*n)
    backward, nb = absorption(starts, kernel, (probability,)*n, reverse=True)
    assert forward == backward
    for j in range(n+1):
        for k in range(n+1):
            observed = sum(forward[j*n+v, mask] for v in range(n)
                           for mask in range(1 << n) if mask.bit_count() == k)/n
            assert observed == addition[k, j]
    for j in range(1, n):
        reverse_mean = [
            sum(a*rank[a]*addition[b, a] for a in range(n+1))/rank[b]
            for b in (j-1, j, j+1)]
        log_slope = 2*derivative[j]/rank[j]-derivative[j-1]/rank[j-1]-derivative[j+1]/rank[j+1]
        assert log_slope == (reverse_mean[0]+reverse_mean[2]-2*reverse_mean[1])/(1-feedback)

    # This lies at the simple-graph endpoint z=n/(n+1); hence the failed
    # sufficient criterion cannot be rescued by restricting to that path.
    n, probability, feedback = 4, sp.Rational(1, 65), sp.Rational(4, 5)
    rank, addition, derivative = reinjection_rank(n, probability, feedback)
    assert list(rank) == list(interpolate(n, sp.Rational(5, 256)))
    birth = [sum(derivative[j+1:])/rank[j] for j in range(n)]+[sp.S.Zero]
    incoming = [sp.S.Zero]+[sp.Rational(j, n-j+1)*birth[j-1] for j in range(1, n+1)]
    concavity = 2*incoming[3]-incoming[2]-incoming[4]
    assert concavity == -sp.Rational(
        255117176659230280777328753275519375,
        1304606218104828027066952281877315611264)
    assert all(2*derivative[j]/rank[j]-derivative[j-1]/rank[j-1]
               -derivative[j+1]/rank[j+1] > 0 for j in range(1, n))
    return {
        "scope": "all-feedback, all-sleep-rate ratio monotonicity at n=3,4 only",
        "symbolic_cases": symbolic,
        "raw_check": {"n": 3, "p": "2/5", "z": "4/9",
                      "stationary_transient_state_counts": list(raw_counts),
                      "addition_transient_state_counts": [nf, nb],
                      "both_orders_every_addition_equation": True,
                      "reverse_regression_identity": True},
        "pure_birth_curvature_obstruction": {
            "n": n, "p": str(probability), "z": str(feedback),
            "simple_graph_sleep_rate": "5/256", "rank": 3,
            "incoming_concavity_gap": str(concavity),
            "rank_law": list(map(str, rank)),
            "pure_birth_rates": list(map(str, birth)),
            "actual_ULC_ratio_derivatives_positive": True}}


def rank_dimension_transfer_certificate():
    """Check the retained n=2 -> 3, beta=1 dimension transfer and ULC obstruction.

    At fixed holding inverse beta, H_N=M_N^-1 Q_N M_N is ARW addition,
    E e_k=e_k and S e_k=e_(k+1). The all-size Pascal identity gives
    P_N=(S+beta H_N E)/(1+beta)=M_N^-1 S M_(N-1), a stochastic map of
    stationary ranks. Positivity uses the canonical ARW/IDLA intertwining.
    Holding beta fixed changes the loopless sleep rate (N+1)/(N*beta).
    This single exact case checks the identities and a nonstationary ULC
    input whose image is not ULC; it does not refute stationary ULC.
    """
    def q_matrix(size):
        result = sp.zeros(size+1)
        for k in range(size+1):
            result[k, k] = sp.Rational(1, size+1-k)
            if k < size:
                result[k+1, k] = sp.Rational(size-k, size+1-k)
        return result

    def formal_coefficient(size, degree, variable, beta):
        # Keep degree=size+1: the formal endpoint is not identically zero.
        return sp.expand(sum(
            sp.ff(size-variable, degree-j)/sp.factorial(degree-j)
            *sp.ff(variable, j)/sp.factorial(j)*(-beta)**j
            for j in range(degree+1)))

    n, N, beta = 2, 3, sp.S.One
    b = 1+beta
    X, z, h = sp.symbols("X z h")
    checks = {}

    def verify(name, condition):
        assert bool(condition), name
        checks[name] = True

    Qold, Q = q_matrix(n), q_matrix(N)
    Mold, M = rank_basis(n, beta), rank_basis(N, beta)
    E, S = sp.zeros(N+1, n+1), sp.zeros(N+1, n+1)
    for k in range(n+1):
        E[k, k], S[k+1, k] = 1, 1
    H = M.inv()*Q*M
    P = (S+beta*H*E)/b
    expected_P = sp.Matrix([
        [sp.Rational(1, 5), sp.Rational(1, 55), sp.Rational(4, 495)],
        [sp.Rational(4, 5), sp.Rational(14, 55), sp.Rational(34, 495)],
        [0, sp.Rational(8, 11), sp.Rational(149, 495)],
        [0, 0, sp.Rational(28, 45)]])
    verify("Q_shift_intertwining", Q*S == S*Qold)
    verify("Pascal_basis_intertwining", M*S+beta*Q*M*E == b*S*Mold)
    verify("two_transfer_forms", P == M.inv()*S*Mold)
    verify("retained_transfer_matrix", P == expected_P)
    verify("H_nonnegative", all(value >= 0 for value in H))
    verify("H_column_stochastic", all(sum(H[:, k]) == 1 for k in range(N+1)))
    verify("H_addition_rank_support", all(H[j, k] == 0
           for k in range(N+1) for j in range(k+2, N+1)))
    verify("P_nonnegative", all(value >= 0 for value in P))
    verify("P_column_stochastic", all(sum(P[:, k]) == 1 for k in range(n+1)))
    rold, rnew = Mold.inv()[:, n], M.inv()[:, N]
    verify("stationary_dimension_transition", P*rold == rnew)
    verify("stationary_vectors_normalized_positive",
           sum(rold) == sum(rnew) == 1 and all(v > 0 for v in list(rold)+list(rnew)))

    Cold = [formal_coefficient(n, k, X, beta) for k in range(n+2)]
    Cnew = [formal_coefficient(N, k, X, beta) for k in range(N+1)]
    Rold = sp.prod(j-X for j in range(1, n+1))/sp.factorial(n)
    endpoint = -b**(n+1)*X*Rold/(n+1)
    verify("formal_endpoint_coefficient", sp.expand(Cold[n+1]-endpoint) == 0)
    verify("formal_endpoint_not_identically_zero", Cold[n+1] != 0)
    verify("dimension_identity_including_endpoint", all(
        sp.expand(Cnew[k]-Cold[k]-(Cold[k-1] if k else 0)) == 0
        for k in range(N+1)))
    verify("three_term_recurrence_including_old_endpoint", all(
        sp.expand((k+1)*Cold[k+1]-(n-b*X+(beta-1)*k)*Cold[k]
                  +beta*(n-k+1)*(Cold[k-1] if k else 0)) == 0
        for k in range(n+1)))
    Aold = [sp.expand(z**k*Cold[k].subs(X, 1/z-1)) for k in range(n+1)]
    Anew = [sp.expand(z**k*Cnew[k].subs(X, 1/z-1)) for k in range(N+1)]
    verify("polynomial_lift_each_column", all(
        sp.expand(b*((N+1)*z-1)*Aold[k]-(k+1)*Anew[k+1]
                  -beta*(N-k)*z*Anew[k]) == 0 for k in range(n+1)))
    verify("multiplication_by_z_each_column", all(
        sp.expand(z*Anew[k]-sum(sp.binomial(N, k)/sp.binomial(N, j)
                  *H[j, k]*Anew[j] for j in range(N+1))) == 0
        for k in range(n+1)))
    verify("direct_natural_coefficient_map_each_column", all(
        sp.expand(((N+1)*z-1)/N*Aold[k]/sp.binomial(n, k)
                  -sum(P[j, k]*Anew[j]/sp.binomial(N, j)
                       for j in range(N+1))) == 0 for k in range(n+1)))
    for size, polys, rank in [(n, Aold, rold), (N, Anew, rnew)]:
        target = sp.prod((j+1)*z-1 for j in range(1, size+1))/sp.factorial(size)
        verify(f"stationary_polynomial_identity_n{size}", sp.expand(
            sum(rank[k]*polys[k]/sp.binomial(size, k) for k in range(size+1))
            -target) == 0)

    raw_input = sp.Matrix([1, 2*h, h**2])
    raw_output = P*raw_input
    raw_normalized = [raw_output[j]/sp.binomial(N, j) for j in range(N+1)]
    raw_gaps = [sp.factor(raw_normalized[j]**2-raw_normalized[j-1]*raw_normalized[j+1])
                for j in range(1, N)]
    retained_raw_gaps = [
        -(632*h**4-450*h**3-7299*h**2+14256*h-156816)/sp.Integer(2205225),
        -h**2*(9215*h**2+18288*h-152496)/sp.Integer(2205225)]
    verify("retained_boundary_family_gap_polynomials", all(
        sp.expand(a-c) == 0 for a, c in zip(raw_gaps, retained_raw_gaps)))
    witness = sp.Matrix([sp.Rational(1, 25), sp.Rational(8, 25), sp.Rational(16, 25)])
    output = P*witness
    expected_output = sp.Matrix([
        sp.Rational(47, 2475), sp.Rational(1948, 12375),
        sp.Rational(5264, 12375), sp.Rational(448, 1125)])
    input_normalized = [witness[k]/sp.binomial(n, k) for k in range(n+1)]
    input_gap = sp.factor(input_normalized[1]**2-input_normalized[0]*input_normalized[2])
    output_normalized = [output[k]/sp.binomial(N, k) for k in range(N+1)]
    output_gaps = [sp.factor(output_normalized[j]**2
                            -output_normalized[j-1]*output_normalized[j+1])
                   for j in range(1, N)]
    verify("witness_is_Binomial_2_four_fifths", witness == raw_input.subs(h, 4)/25)
    verify("witness_is_ULC", input_gap == 0 and all(v > 0 for v in witness) and sum(witness) == 1)
    verify("retained_witness_output", output == expected_output)
    verify("retained_negative_ULC_gap", output_gaps[1] == sp.Rational(-1089536, 1378265625))
    verify("witness_is_not_stationary", witness != rold)
    verify("gap_family_specialization", all(
        sp.cancel(raw_gaps[j].subs(h, 4)/625-output_gaps[j]) == 0 for j in range(2)))
    return {
        "scope": "One exact n=2 -> 3, holding beta=1 check of the all-size dimension identity",
        "parameters": {"old_size": n, "new_size": N, "holding_inverse_beta": str(beta),
                       "loopless_sleep_rate_old": str(sp.Rational(n+1, n)/beta),
                       "loopless_sleep_rate_new": str(sp.Rational(N+1, N)/beta)},
        "formulas": {"H_N": "M_N^(-1) Q_N M_N",
                     "P_N": "(S+beta H_N E)/(1+beta)=M_N^(-1) S M_(N-1)",
                     "stationary_transition": "r^(N)=P_N r^(N-1)",
                     "formal_endpoint": "C_(n,n+1)(X)=-(1+beta)^(n+1) X R_n(X)/(n+1)"},
        "transfer_matrix": [[str(P[j, k]) for k in range(n+1)] for j in range(N+1)],
        "stationary_rank_old": list(map(str, rold)),
        "stationary_rank_new": list(map(str, rnew)),
        "formal_endpoint_at_checked_point": str(sp.factor(endpoint)),
        "ULC_preservation_obstruction": {
            "input": list(map(str, witness)),
            "input_binomially_normalized": list(map(str, input_normalized)),
            "input_gap": str(input_gap), "output": list(map(str, output)),
            "output_normalized_gaps_j1_j2": list(map(str, output_gaps)),
            "raw_input_family": "(1,2*h,h^2)",
            "raw_output_normalized_gaps": list(map(str, raw_gaps)),
            "interpretation": "P_N does not preserve arbitrary ULC inputs. The witness is nonstationary and is not a counterexample to stationary ULC."},
        "checks": checks, "checks_passed": len(checks),
        "proof_status": {
            "all_size_positive_transition": "Audited Pascal identity and canonical ARW/IDLA intertwining; not inferred from this finite check",
            "generic_ULC_preservation": "False by the exact witness",
            "stationary_ULC": "Open"}}


def survival_top_moment_certificate():
    """Fixed analytic certificate for the truncated-exponential mixture lemma."""
    ell, dd, u, v, p, slack = sp.symbols("L D u v p C", positive=True)
    gg = v/u
    gap = gg-2*u+u*u+u*(gg*gg/(2*p)+slack)
    moment_gap = v+v*v/2-u*u*(2-u)
    assert sp.factor(u*gap-moment_gap-(1-p)*v*v/(2*p)-u*u*slack) == 0
    mean = 1/ell-1/dd
    second = 2/ell**2-(1+2/ell)/dd
    numerator = sp.factor((second+second**2/2-mean**2*(2-mean))*2*ell**4*dd**3)
    exponentials = {
        3: 2*ell+4,
        2: -2*ell**4+4*ell**3-10*ell**2-14*ell-12,
        1: ell**4+2*ell**3+24*ell**2+22*ell+12,
        0: -ell**4-6*ell**3-14*ell**2-10*ell-4}
    symbol = sp.Symbol("E")
    assert sp.expand(numerator.subs(dd, symbol-1)
                     -sum(poly*symbol**rate for rate, poly in exponentials.items())) == 0
    phi = sp.sqrt(1+4*u*u-2*u**3)-1
    curvature = 4-6*u-8*u**3+3*u**4
    assert sp.simplify(sp.diff(phi, u, 2)
                       -curvature/(1+4*u*u-2*u**3)**sp.Rational(3, 2)) == 0
    assert curvature.subs(u, sp.Rational(1, 2)) == sp.Rational(3, 16)
    assert sp.expand(sp.diff(curvature, u)+6+12*u*u*(2-u)) == 0
    initial = []
    for rate in [0]*5+[1]*5:
        initial.append(sum(poly.subs(ell, 0) for poly in exponentials.values()))
        exponentials = {a: sp.expand(sp.diff(poly, ell)+(a-rate)*poly)
                        for a, poly in exponentials.items()}
    assert initial == [0, 0, 0, 0, 0, 0, 0, 140, 700, 4088]
    assert exponentials[0] == exponentials[1] == 0
    assert exponentials[3] == 15552*ell+95904
    assert exponentials[2] == -64*ell**4-1792*ell**3-16640*ell**2-60928*ell-76864
    taylor = sum(ell**j/sp.factorial(j) for j in range(5))
    forcing = sp.Poly(sp.expand(exponentials[3]*taylor+exponentials[2]), ell)
    assert forcing.all_coeffs() == [648, 6524, 21968, 46864, 50528, 19040]
    assert moment_gap.subs({u: sp.Rational(1, 2), v: sp.Rational(1, 3)}) == sp.Rational(1, 72)
    return {"scope": "Every normalized mixture of exp(-L theta) on [0,1], finite L>=0",
            "strict_moment_inequality": "v+v^2/2>u^2(2-u)",
            "gap_decomposition": "uH=K(u,v)+(1-p)v^2/(2p)+u^2 C, C>=0",
            "differential_rates": [0]*5+[1]*5,
            "initial_values": list(map(int, initial)),
            "positive_forcing_coefficients_descending": list(map(int, forcing.all_coeffs())),
            "positive_remainder_factor": "15552L+95904",
            "curvature_numerator_minimum_on_0_to_half": "3/16",
            "uniform_component_gap": "1/72",
            "analytic_strictness": "Positive forcing, backward integration and convex-mixture closure"}


def survival_top_lc_certificate():
    """All-size top survival LC at boundary rate one, report (55g).

    This proves only a_(m-1)^2>a_(m-2)*a_m at infinite time; this function
    makes no assertion about earlier gaps or finite-time LC. Three preceding
    first-step equations bound the actual passage value and derivative;
    a small exact polynomial signs the resulting uniform lower bound.
    """
    import hashlib

    p, w, j = sp.symbols("p w j")
    t, u, z = sp.symbols("t u z")
    q, h = 1-p, w+2  # m=w+4; h=m-2.
    gg = 2*p/(1+2*p+h*q*(1-j))
    ff = p/(1+p+(h+1)*q*(1-gg))
    lower_gap = sp.factor(gg-2*ff+ff**2
                          +ff*gg**2/(2*p)*(1+h*q*j**2/(3*p)))
    preceding_lower = 4*p/(1+4*p+w*q)
    lower = 3*p/(1+3*p+(w+1)*q*(1-preceding_lower))
    upper = 3*p/(1+3*p)
    assert sp.cancel(upper-lower-(
        3*p*q*(w+1)*(1+w*q)
        /((1+3*p)*((1+3*p)*(1+4*p+w*q)+q*(w+1)*(1+w*q))))) == 0
    substituted = sp.factor(lower_gap.subs(j, (lower+u*upper)/(1+u))
                            .subs(p, t/(1+t)))
    numerator, denominator = sp.fraction(substituted)
    polynomial = sp.Poly(numerator, t, w, u)
    rows = [[*powers, int(value)] for powers, value in sorted(polynomial.terms())]
    assert len(rows) == 389 and all(row[-1] > 0 for row in rows)
    assert polynomial.total_degree() == 17
    # Coefficients of u^k, k=0,...,3, are positive polynomials in t,w.
    counts = [len(sp.Poly(numerator.expand().coeff(u, k), t, w).terms())
              for k in range(4)]
    assert counts == [95, 98, 98, 98]
    assert all(c >= 0 for c in sp.Poly(denominator, t, w, u).coeffs())
    assert denominator.subs({t: 0, w: 0, u: 0}) > 0
    assert numerator.expand().coeff(w, 0).coeff(u, 0).coeff(t, 1) > 0
    digest = hashlib.sha256(json.dumps(rows, separators=(",", ":")).encode()).hexdigest()

    def passages(m, probability):
        values = []
        for index in range(m):
            previous = values[-1] if values else 0
            values.append(sp.factor(probability*(m-index)/(z+probability*(m-index)
                          +(1-probability)*index*(1-previous))))
        return values

    small = {}
    for m in (2, 3):
        values = passages(m, p)
        fz, gz = values[-1], values[-2]
        gap = sp.factor((gz-2*fz+fz**2-fz*sp.diff(gz, z)).subs(z, 1))
        num, den = sp.fraction(sp.factor(gap.subs(p, t/(1+t))))
        assert all(c > 0 for c in sp.Poly(num, t).coeffs())
        assert all(c > 0 for c in sp.Poly(den, t).coeffs())
        assert gap.subs(p, 1) == sp.Rational(1, 36)
        small[str(m)] = {"gap": str(gap), "numerator_t": str(num), "denominator_t": str(den)}

    # These two rational points discriminate the m=4 endpoint of the
    # passage-value interval from an interior m>4 case. The independent
    # killed-resolvent calculation also checks the queue-gap normalization.
    direct = []
    for m, probability in ((4, sp.Rational(2, 5)), (6, sp.Rational(3, 7))):
        values = passages(m, probability)
        fvalue, gvalue, jvalue = [v.subs(z, 1) for v in values[-1:-4:-1]]
        gprime = sp.diff(values[-2], z).subs(z, 1)
        gap = sp.factor(gvalue-2*fvalue+fvalue**2-fvalue*gprime)
        diagonal = boundary_profiles(m, probability)[0]
        assert diagonal[-2] == 1-fvalue
        assert diagonal[-3] == 1-gvalue+fvalue*gprime
        assert gap == diagonal[-2]**2-diagonal[-3] > 0
        bounds = {p: probability, w: m-4}
        lo, hi = lower.subs(bounds), upper.subs(bounds)
        assert lo <= jvalue < hi
        uvalue = sp.factor((jvalue-lo)/(hi-jvalue))
        bound = sp.factor(lower_gap.subs({**bounds, j: jvalue}))
        assert gap >= bound > 0
        assert sp.factor(substituted.subs({t: probability/(1-probability), w: m-4,
                                          u: uvalue})-bound) == 0
        direct.append({"m": m, "p": str(probability), "actual_J": str(jvalue),
                       "interpolation_u": str(uvalue), "gap": str(gap),
                       "lower_bound": str(bound), "independent_resolvent_matches": True})
    moment_certificate = survival_top_moment_certificate()
    # Return-state renewal supplies an independent moment normalization at
    # the first already selected resolvent point, with no additional sweep.
    probability = sp.Rational(2, 5)
    def return_transform(m):
        return sum(sp.binomial(m, j)*probability**(m-j)*(1-probability)**j/(z+j)
                   for j in range(m+1))
    rm, rm1, rm2 = [return_transform(m) for m in (4, 3, 2)]
    fz = sp.cancel(probability/(1-probability)*(rm1/rm-1))
    gz = sp.cancel(probability/(1-probability)*(rm2-2*rm1+rm)/(rm1-rm))
    uu = fz.subs(z, 1)
    vv = uu*gz.subs(z, 1)
    mm = sp.factor(vv+vv**2/2-uu**2*(2-uu))
    cc = sp.factor(-sp.diff(gz, z).subs(z, 1)-gz.subs(z, 1)**2/(2*probability))
    assert mm > 0 and cc >= 0
    assert sp.factor((mm+(1-probability)*vv**2/(2*probability)+uu**2*cc)/uu
                     -sp.Rational(direct[0]["gap"])) == 0
    moment_certificate["return_moment_audit"] = {
        "m": 4, "p": "2/5", "K": str(mm), "C": str(cc),
        "same_gap_as_independent_resolvent": True}
    return {
        "status": "proved by exact nonnegative-coefficient certificate and passage bounds",
        "scope": "Every integer m>=2, 0<p<1; boundary rate one and infinite horizon only",
        "inequality": "a_(m-1)^2>a_(m-2)*a_m, a_m=1",
        "parameters": "m=w+4,p=t/(1+t),J=(lower+u*upper)/(1+u), t>0,w,u>=0",
        "lower_J": str(lower), "upper_J": str(upper), "lower_gap": str(lower_gap),
        "polynomial": {"variables": ["t", "w", "u"], "positive_terms": len(rows),
                       "total_degree": 17, "counts_by_u_degree": counts,
                       "sha256": digest, "ascending_coefficient_rows": rows,
                       "denominator": str(denominator), "denominator_positive": True},
        "small_dimensions": small, "direct_audits": direct,
        "independent_moment_proof": moment_certificate,
        "p1_endpoint_gap": "1/36", "all_diagonal_gaps_or_stationary_ULC_proved": False}


def survival_next_lc_certificate():
    """All-size second infinite-horizon survival gap, report (55i).

    The last four passage jets give a quadratic with a positive leading
    coefficient. Keeping two more first-step equations bounds its constant
    term by a positive polynomial. This does not assert stationary ULC.
    The optional python-flint dependency handles only exact composition.
    """
    import hashlib
    from flint import fmpz_mpoly_ctx

    p, w, k, t, u, v, z = sp.symbols("p w k t u v z")
    q = 1-p
    alpha, beta = q*(w+1), q*(w+2)  # Initially m=w+4.
    jj = 3*p/(1+3*p+alpha*(1-k))
    gg = 2*p/(1+2*p+beta*(1-jj))
    ff = p/(1+p+q*(w+3)*(1-gg))
    aa = 1-ff
    d0, d1 = gg**2/(2*p), beta*gg**2/(2*p)
    b0 = 1-gg-ff*d0
    c0 = b0**2-aa*(1-jj)
    c1 = -2*b0*ff*d1+aa*(gg+ff*d0)
    c2 = ff**2*d1**2+aa*ff*(gg/jj+d1)
    yfloor = jj**2/(3*p)*(1+alpha*k**2/(4*p))
    efloor = alpha*jj**2*k**3/(48*p**3)
    lower_constant = c0+c1*yfloor+c2*yfloor**2+aa*gg*ff*efloor
    lower_linear = c1+2*c2*yfloor
    lower = 4*p/(1+4*p+w*q)
    upper = 4*p/(1+4*p)
    interval = (lower+u*upper)/(1+u)

    def coefficient_data(poly):
        rows = [[*map(int, powers), int(value)]
                for powers, value in sorted(poly.terms())]
        assert all(row[-1] > 0 for row in rows)
        return {"variables": list(map(str, poly.gens)),
                "positive_terms": len(rows), "total_degree": int(poly.total_degree()),
                "degrees": list(map(int, poly.degree_list())),
                "sha256": hashlib.sha256(json.dumps(rows, separators=(",", ":")).encode()).hexdigest(),
                "ascending_coefficient_rows": rows}

    def positive_denominator(expression):
        scalar, factors = sp.factor_list(expression)
        assert scalar > 0
        data = []
        for factor, multiplicity in factors:
            poly = sp.Poly(factor, t, w, u)
            assert all(c > 0 for c in poly.coeffs())
            assert factor.subs({t: 0, w: 0, u: 0}) > 0
            data.append({"factor": str(factor), "multiplicity": int(multiplicity)})
        return {"positive_scalar": str(scalar), "positive_factors": data}

    expressions = {}
    denominators = {}
    for name, value in (("linear", lower_linear), ("constant", lower_constant)):
        rational = sp.factor(sp.factor(value).subs(k, interval).subs(p, t/(1+t)))
        numerator, denominator = sp.fraction(rational)
        expressions[name] = sp.Poly(numerator, t, w, u)
        denominators[name] = positive_denominator(denominator)
        print("Next survival gap: " + name + " numerator reconstructed", flush=True)

    linear = coefficient_data(expressions["linear"])
    assert linear["positive_terms"] == 1170 and linear["total_degree"] == 25
    old = expressions["constant"]
    assert len(old.terms()) == 1285 and sum(1 for c in old.coeffs() if c < 0) == 65
    degree = int(old.degree(u))
    assert degree == 5

    # For m=w+6 retain the birth-6p, birth-5p and birth-4p passages,
    # replacing only the transform preceding those by its lower bound zero.
    refined = sp.S.Zero
    for birth in range(6, 3, -1):
        refined = sp.cancel(birth*t/(1+(birth+1)*t+(w+6-birth)*(1-refined)))
    old_lower = 4*t/(1+5*t+w+2)
    upper_t = 4*t/(1+5*t)
    u0 = sp.factor((refined-old_lower)/(upper_t-refined))
    a, b = sp.fraction(u0)
    assert sp.factor(u0-5*t*(5*t+1)*(7*t+w+1)
                     /((5*t+w+3)*(7*t**2+2*t*w+9*t+w**2+3*w+2))) == 0
    assert all(c > 0 for c in sp.Poly(a, t, w).coeffs())
    assert all(c > 0 for c in sp.Poly(b, t, w).coeffs())
    assert b.subs({t: 0, w: 0}) > 0
    # K=(refined+v*upper)/(1+v) gives old_u=u0+(1+u0)*v.
    assert sp.cancel(((old_lower+u*upper_t)/(1+u)).subs(u, u0+(1+u0)*v)
                     -(refined+v*upper_t)/(1+v)) == 0

    ctx = fmpz_mpoly_ctx.get(["t", "w", "v"])
    tt, ww, vv = ctx.gens()
    def convert(expression):
        return ctx.from_dict({(i, j, 0): int(c)
                              for (i, j), c in sp.Poly(expression, t, w).terms()})
    an, bn = convert(a), convert(b)
    combined = an+(an+bn)*vv
    parts = [{} for _ in range(degree+1)]
    for (i, j, ell), value in old.terms():
        parts[ell][i, j, 0] = int(value)
    refined_numerator = ctx.constant(0)
    for ell in range(degree+1):
        refined_numerator += (ctx.from_dict(parts[ell]).compose(tt, ww+2, vv, ctx=ctx)
                              *combined**ell*bn**(degree-ell))
    refined_poly = sp.Poly.from_dict(
        {tuple(map(int, refined_numerator.monomial(i))): int(refined_numerator.coefficient(i))
         for i in range(len(refined_numerator))}, (t, w, v))
    refined_data = coefficient_data(refined_poly)
    assert refined_data["positive_terms"] == 3895 and refined_data["degrees"] == [36, 30, 5]
    assert refined_data["sha256"] == "07bf0de0d0a198ac1110368c95dde732bdf3a083d405a92be77da3a2d29483a4"
    assert refined_poly.coeff_monomial(t) > 0

    # Initial passages are exponential. Dimensions four and five attain the
    # corresponding exact continued-fraction endpoint; no limit is needed.
    base4 = sp.Poly(old.as_expr().subs({w: 0, u: 0}), t)
    a5, b5 = 5*t*(1+5*t), (5*t+2)*(t+1)
    base5 = sp.Poly(sum(old.as_expr().coeff(u, ell).subs(w, 1)
                       *a5**ell*b5**(degree-ell) for ell in range(degree+1)), t)
    base4_data, base5_data = coefficient_data(base4), coefficient_data(base5)
    assert base4_data["positive_terms"] == 21 and base5_data["positive_terms"] == 31
    assert base4.coeff_monomial(t) > 0 and base5.coeff_monomial(t) > 0
    trans = []
    for index in range(3):
        prior = trans[-1] if trans else 0
        trans.append(sp.factor(p*(3-index)/(z+p*(3-index)+q*index*(1-prior))))
    f3, g3, j3 = [value.subs(z, 1) for value in trans[::-1]]
    y3, z3 = [-sp.diff(trans[0], z).subs(z, 1), sp.diff(trans[0], z, 2).subs(z, 1)/2]
    d3 = -sp.diff(trans[1], z).subs(z, 1)
    gap3 = sp.factor((1-g3-f3*d3)**2-(1-f3)*(1-j3-y3*g3-y3*d3*f3-z3*g3*f3))
    num3, den3 = sp.fraction(sp.factor(gap3.subs(p, t/(1+t))))
    base3_data = coefficient_data(sp.Poly(num3, t))
    assert all(c > 0 for c in sp.Poly(den3, t).coeffs())
    assert gap3.subs(p, 1) == sp.Rational(53, 10368)

    def passage_jets(m, probability):
        jets = []
        previous, first, second = sp.S.Zero, sp.S.Zero, sp.S.Zero
        for index in range(m):
            birth, death = probability*(m-index), (1-probability)*index
            value = sp.cancel(birth/(1+birth+death*(1-previous)))
            slope = 1-death*first
            derivative = -value**2/birth*slope
            curvature = 2*value**3/birth**2*slope**2+value**2/birth*death*second
            previous, first, second = map(sp.cancel, (value, derivative, curvature))
            jets.append((previous, first, second))
        return jets

    direct = []
    for m, probability in ((3, sp.Rational(2, 5)), (4, sp.Rational(2, 5)),
                           (5, sp.Rational(2, 5)), (6, sp.Rational(3, 7)),
                           (7, sp.Rational(2, 5))):
        jets = passage_jets(m, probability)
        fv, gv, jv = [entry[0] for entry in jets[-1:-4:-1]]
        y, zeta, d = -jets[-3][1], jets[-3][2]/2, -jets[-2][1]
        av, bv, cv = 1-fv, 1-gv-fv*d, 1-jv-y*gv-y*d*fv-zeta*gv*fv
        diagonal = boundary_profiles(m, probability)[0]
        assert tuple(diagonal[-4:-1]) == (cv, bv, av)
        actual = sp.cancel(bv*bv-av*cv)
        assert actual > 0
        audit = {"m": m, "p": str(probability), "actual_gap": str(actual),
                 "three_survival_values_match_independent_resolvent": True}
        if m >= 4:
            kv = jets[-4][0]
            substitution = {p: probability, w: m-4, k: kv}
            yf, ef = map(lambda expr: sp.cancel(expr.subs(substitution)), (yfloor, efloor))
            e = sp.cancel(zeta-y*y/jv)
            assert y >= yf and e >= ef >= 0
            lc, ll = map(lambda expr: sp.cancel(expr.subs(substitution)),
                         (lower_constant, lower_linear))
            assert lc > 0 and ll > 0
            cc2 = sp.cancel(c2.subs(substitution))
            assert actual == sp.cancel(lc+ll*(y-yf)+cc2*(y-yf)**2+av*gv*fv*(e-ef))
            audit.update({"derivative_floor": str(yf), "variance_floor": str(ef),
                          "constant_lower_bound": str(lc), "shifted_linear_coefficient": str(ll)})
            if m >= 6:
                tv, wv = probability/(1-probability), m-6
                lo, hi = [sp.cancel(expr.subs({t: tv, w: wv})) for expr in (refined, upper_t)]
                assert lo <= kv < hi
                vv_actual = sp.cancel((kv-lo)/(hi-kv))
                assert vv_actual >= 0 and (vv_actual == 0) == (m == 6)
                audit["refined_interpolation_v"] = str(vv_actual)
        direct.append(audit)
    return {
        "status": "proved by passage bounds and exact positive-coefficient certificates",
        "scope": "Every integer m>=3, 0<p<1; boundary rate one and infinite horizon only",
        "inequality": "a_(m-2)^2>a_(m-3)*a_(m-1)",
        "quadratic": "gap=C0+C1*s+c2*s^2+A*G*F*(e-efloor), s=y-yfloor>=0",
        "floors": {"y": "J^2/(3p)*(1+q(m-3)K^2/(4p))", "e": "q(m-3)J^2*K^3/(48p^3)"},
        "linear_polynomial": linear, "original_positive_denominators": denominators,
        "refined_parameterization": "m=w+6,p=t/(1+t),old_w=w+2,old_u=u0+(1+u0)*v; t>0,w,v>=0",
        "refined_lower_K": str(refined), "u0": str(u0),
        "constant_polynomial": refined_data,
        "positive_homogenizing_denominator": str(b**degree),
        "strict_t_coefficient": str(refined_poly.coeff_monomial(t)),
        "small_dimensions": {"3": {"gap": str(gap3), "numerator": base3_data, "denominator_t": str(den3)},
                             "4": base4_data, "5": base5_data},
        "direct_audits": direct, "p1_endpoint_gap": "53/10368",
        "all_diagonal_gaps_or_stationary_ULC_proved": False}


def survival_third_reduction_certificate():
    """Third survival-gap reduction and two exact relaxation obstructions.

    Actual zeroth passage transforms, m>=4, 0<p<1, boundary rate one,
    infinite horizon. The strict R0 floor is the proved input from (55i);
    its large coefficient certificate is not rerun. Local identities,
    normalized passage moments and one asymptotic endpoint are checked.
    Neither H0 at the actual derivative nor the third gap is proved positive.
    """
    from itertools import product
    from math import comb

    def zero(value):
        assert sp.cancel(value) == 0

    F, G, J, K, x, v, w, y, z2, d = sp.symbols("F G J K x v w y z2 d")
    alpha, r, d0, d1, eps, eta, x0, shift = sp.symbols("alpha r d0 d1 epsilon eta x0 s")
    p, a, b = sp.symbols("p a b")
    # Normalized jets are Pr(Poisson(T)=j)=(-1)^j P^(j)(1)/j!.
    jets = ((K, x, v, w), (J, y, z2), (G, d), (F,))
    prefixes, failure = {}, sp.S.Zero
    for length in range(1, 5):
        words = [word for word in product(range(length), repeat=length)
                 if sum(word) == length-1
                 and all(sum(word[:j]) >= j for j in range(1, length))]
        prefixes[str(length)] = [list(word) for word in words]
        failure += sum(sp.prod(jets[j][count] for j, count in enumerate(word)) for word in words)
    assert [len(words) for words in prefixes.values()] == [1, 1, 2, 5]
    A = 1-F
    Braw = 1-G-F*d
    Craw = 1-J-y*G-y*d*F-z2*G*F
    Draw = 1-K-x*(1-Craw)-v*(J*(1-Braw)+y*G*F)-w*J*G*F
    zero(Draw-(1-failure))

    # Invert the local first-step denominator through order three. Here
    # lj=(-1)^j L^(j)(1)/j!, lambda=4p and c=q(m-4); omit cL at m=4.
    lam, c, l1, l2, l3 = sp.symbols("lambda c l1 l2 l3")
    denominator = (lam/K, 1+c*l1, -c*l2, c*l3)
    series = [K]
    for j in range(1, 4):
        series.append(sp.expand(-K/lam*sum(denominator[i]*series[j-i] for i in range(1, j+1))))
    first, slack2, slack3 = K*K*(1+c*l1)/lam, c*K*K*l2/lam, c*K*K*l3/lam
    zero(-series[1]-first)
    zero(series[2]-first**2/K-slack2)
    zero(-series[3]-first**3/K**2-2*first*slack2/K-slack3)
    al, rr = J*J/(3*p), a*J*J/(3*p)
    zero(J*J*(1+a*x)/(3*p)-(al+rr*x))
    zero(J**3*(1+a*x)**2/(3*p)**2+a*J*J*v/(3*p)-(al+rr*x)**2/J-rr*v)
    zero(G*G*(1+b*y)/(2*p)-(G*G/(2*p)+b*G*G*y/(2*p)))

    yy, vv, ww = alpha+r*x, x*x/K+eps, x**3/K**2+2*x*eps/K+eta
    zz, dd = yy**2/J+r*vv, d0+d1*yy
    substitution = {y: yy, z2: zz, d: dd, v: vv, w: ww}
    B, C, D = [sp.expand(value.subs(substitution, simultaneous=True)) for value in (Braw, Craw, Draw)]
    beta, Q = F*G*r, J*G*F
    B0, B1 = 1-G-F*(d0+d1*alpha), F*d1*r
    U0 = 1-J-(G+F*d0)*alpha-F*(d1+G/J)*alpha**2
    U1 = r*(G+F*d0+2*F*(d1+G/J)*alpha)
    U2 = r*r*F*(d1+G/J)+beta/K
    U = U0-U1*x-U2*x*x
    T0, T1 = J*(1-B0)+G*F*alpha, J*B1+beta
    D0, D1, D2, D3 = 1-K, 1-U0, U1+T0/K, U2+T1/K+Q/K**2
    M1 = 2*beta+J*B1+2*Q/K
    zero(B-(B0-B1*x))
    zero(C-(U-beta*eps))
    zero(D-(D0-D1*x-D2*x*x-D3*x**3-eps*(T0+M1*x)-Q*eta))
    h = (U0**2-B0*D0, -2*U0*U1+B0*D1+B1*D0,
         U1**2-2*U0*U2+B0*D2-B1*D1, 2*U1*U2+B0*D3-B1*D2, U2**2-B1*D3)
    e = (-2*beta*U0+B0*T0, 2*beta*U1+B0*M1-B1*T0, 2*beta*U2-B1*M1)
    H0 = sum(value*x**j for j, value in enumerate(h))
    Epoly = sum(value*x**j for j, value in enumerate(e))
    H = C*C-B*D
    zero(H-H0-eps*Epoly-beta**2*eps**2-B*Q*eta)
    Gamma = J*(1-B)+G*F*(2*yy-alpha+2*x*J/K)
    zero(sp.diff(H, eps)-(B*Gamma-2*beta*C))
    zero(sp.diff(H, eta)-B*Q)
    zero(Epoly-(B*Gamma-2*beta*U))
    zero(H0.subs(x, x0+shift)-sum(
        sum(sp.binomial(j, i)*h[j]*x0**(j-i) for j in range(i, 5))*shift**i for i in range(5)))

    # The (55i) constant floor R0 is positive at actual F,G,J,K. Its exact
    # propagation holds for every formal x>=x0, epsilon>=0, without B>0.
    Bbase = 1-G-F*d0
    c0, c1 = Bbase**2-A*(1-J), -2*Bbase*F*d1+A*(G+F*d0)
    c2 = F*F*d1*d1+A*F*(G/J+d1)
    y0, e0 = alpha+r*x0, r*x0*x0/K
    R0, R1 = c0+c1*y0+c2*y0*y0+A*G*F*e0, c1+2*c2*y0
    propagated = (R0+(r*R1+2*A*G*F*r*x0/K)*shift
                  +(c2*r*r+A*G*F*r/K)*shift**2+A*G*F*r*eps)
    zero((B*B-A*C).subs(x, x0+shift)-propagated)
    zero(e0.subs({r: rr, x0: K*K/(4*p)}, simultaneous=True)-a*J*J*K**3/(48*p**3))
    q = sp.symbols("q")
    c1_positive = G*(1-F-G)+F*d0*(5-F-2*q*G)+2*F**2*b*d0**2
    zero((c1-c1_positive).subs({d1: b*G*G/(2*p), d0: G*G/(2*p)})
         .subs(F, p/(2+b-(b+q)*G)).subs(q, 1-p))
    # Analytically c1>=G/3>0. On B>0, C<=B^2/A, use Fa(1-G)<p,
    # B<=1-G and 2J<=3A to get 2 beta C<BGJ<B Gamma. Hence Epoly>0
    # at epsilon=0, and both higher slacks can be lowered all the way to zero.

    # The killed reversible diagonal density is 4p sum c_j exp(-lambda_j t).
    # Weights omega_j=4p*c_j/lambda_j sum to one. For R=1/(1+lambda_j),
    # (K,x,v,w)=E[(1-R)*(1,R,R^2,R^3)]. Cauchy--Schwarz gives these signs.
    zero((1-K)*(1-K-x-v)-(1-K-x)**2-((1-K)*(x-v)-x*x))
    zero(x-x*x/(1-K)-x*x/K-x*(1-x/(K*(1-K))))
    zero(x*ww-vv**2-(x*eta-eps**2))
    rho = sp.symbols("rho")
    zero(x*H.subs(eta, eps**2/x+rho)
         -(x*H0+x*eps*Epoly+(x*beta**2+B*Q)*eps**2+x*B*Q*rho))

    # Fixed m, q->0: U2=O(q), B1=O(q^2), so only endpoint values enter.
    m = sp.symbols("m")
    rlead = sp.Rational(3, 16)*(m-3)
    B1lead = sp.Rational(1, 2)*sp.Rational(2, 9)*(m-2)*rlead
    U2lead = sp.Rational(1, 2)*sp.Rational(2, 3)*rlead/sp.Rational(4, 5)
    D3constant = sp.Rational(1, 4)/sp.Rational(4, 5)**2
    h4lead = sp.factor(U2lead**2-B1lead*D3constant)
    zero(h4lead+25*(m-3)*(m+1)/sp.Integer(12288))

    # Distinct simultaneous limit p_m=1-2/m. The return densities are
    # proportional to [1-2(1-theta)/m]^(-m-2), bounded by exp(6) for m>=4.
    # Their normalized limit is exp(-2theta)/Z0. Integral IBP gives Zj below.
    exponential, u = sp.symbols("E u")
    moments = [(1-1/exponential)/2]
    for j in range(1, 5):
        moments.append(sp.cancel(j*moments[-1]/2-1/(2*exponential)))
    ff, gg, jj, kk = [sp.cancel(moments[j]/moments[j-1]) for j in range(1, 5)]
    xx = kk*kk/4
    ypoint, vpoint, wpoint = jj*jj*(1+2*xx)/3, xx*xx/kk, xx**3/kk**2
    dpoint, zpoint = gg*gg*(1+2*ypoint)/2, ypoint**2/jj+2*jj*jj*vpoint/3
    point = {F: ff, G: gg, J: jj, K: kk, x: xx, y: ypoint,
             z2: zpoint, d: dpoint, v: vpoint, w: wpoint}
    bp, up, dp = [sp.cancel(value.subs(point, simultaneous=True))
                  for value in (Braw, Craw, 1-failure)]
    gap = sp.factor(up*up-bp*dp)
    numerator, denominator = sp.fraction(gap)
    expected = (11298771*exponential**7-479855421*exponential**6+8727946911*exponential**5
                -88141154049*exponential**4+533853269433*exponential**3
                -1939918623255*exponential**2+3917702255253*exponential-3393875313259)
    zero(numerator-expected)
    zero(denominator-82944*(exponential-5)**2*(exponential-3)**2*(exponential-1)**2*(3*exponential-19))

    def interval_rows(polynomial):
        """Integer binomial convolution for E=(369+370u)/(50(1+u))."""
        poly = sp.Poly(polynomial, exponential)
        degree = int(poly.degree())
        rows = [0]*(degree+1)
        for j in range(degree+1):
            value = int(poly.nth(j))
            assert value == poly.nth(j)
            for ell in range(j+1):
                left = value*comb(j, ell)*369**(j-ell)*370**ell*50**(degree-j)
                for hdegree in range(degree-j+1):
                    rows[ell+hdegree] += left*comb(degree-j, hdegree)
        return {"degree": degree, "positive_clearing_multiplier": str(50**degree),
                "ascending_coefficients": rows}

    negative = interval_rows(numerator)
    assert negative["ascending_coefficients"] == [
        -13193945604333667331, -85483057479080102960, -235533957083858054400,
        -357205201569817120000, -321354182948185600000, -170997904761676800000,
        -49625007177728000000, -6020929003520000000]
    assert all(value < 0 for value in negative["ascending_coefficients"])
    domain = {}
    for name, value, den in (("B", bp, 24*(exponential-3)*(exponential-1)),
                             ("U", up, 288*(exponential-5)*(exponential-3)*(exponential-1)),
                             ("A_squared_minus_B", (1-ff)**2-bp, 24*(exponential-3)*(exponential-1)**2)):
        num, actual_den = sp.fraction(sp.factor(value))
        zero(actual_den-den)
        rows = interval_rows(num)
        assert all(coefficient > 0 for coefficient in rows["ascending_coefficients"])
        domain[name] = {"formula": str(sp.factor(value)), "positive_interval_certificate": rows}
    lower = sum(sp.Rational(2)**j/sp.factorial(j) for j in range(9))
    upper = sum(sp.Rational(2)**j/sp.factorial(j) for j in range(7))
    upper += sp.Rational(2)**7/sp.factorial(7)/(1-sp.Rational(1, 4))
    assert lower == sp.Rational(2327, 315) > sp.Rational(369, 50)
    assert upper == sp.Rational(6983, 945) < sp.Rational(37, 5)
    assert lower > sp.Rational(155, 21) and sp.Rational(37, 5) < sp.Rational(67, 9)
    zero(kk-sp.Rational(8, 11)-2*(21*exponential-155)/(11*(3*exponential-19)))
    zero(sp.Rational(4, 5)-kk-2*(67-9*exponential)/(5*(3*exponential-19)))
    L = sp.symbols("L")
    fifth_floor = K*K/(4*p)*(1+c*L*L/(5*p))
    zero(fifth_floor-K*K/(4*p)-c*K*K*L*L/(20*p*p))
    limiting_L = sp.cancel((7-4/kk)/2)
    zero(limiting_L-(15*exponential-109)/(6*(exponential-7)))
    zero(limiting_L-sp.Rational(3, 4)-(21*exponential-155)/(12*(exponential-7)))

    return {
        "status": "audited reduction and relaxation obstructions; third survival gap OPEN",
        "scope": "Actual zeroth passage transforms; integer m>=4, 0<p<1; infinite horizon, boundary rate one",
        "target": "H=C^2-BD, (A,B,C,D)=(a_(m-1),a_(m-2),a_(m-3),a_(m-4))",
        "queue": {"first_failure_prefixes": prefixes, "prefix_counts": [1, 1, 2, 5],
                  "normalization": "jet_j=(-1)^j P^(j)(1)/j!=Pr(Poisson(T)=j)",
                  "four_passage_identity": "D=1-K-x(1-C)-v[J(1-B)+yGF]-wJGF"},
        "identities": {"local_first_step_jets": True, "quartic_and_higher_slack_coefficients": True,
                       "quartic_translation": True, "higher_slack_derivatives": True,
                       "strong_second_gap_propagation": True, "moment_bounds_algebra": True},
        "jet_relations": {"a": "q(m-3)", "b": "q(m-2)", "alpha": "J^2/(3p)",
                          "r": "a*alpha", "d0": "G^2/(2p)", "d1": "b*d0",
                          "y": "alpha+r*x", "d": "d0+d1*y", "z2": "y^2/J+r*v",
                          "v": "x^2/K+epsilon", "w": "x^3/K^2+2x*epsilon/K+eta"},
        "quartic": {"identity": "H=H0(x)+epsilon*E(x)+beta^2*epsilon^2+B(x)*Q*eta",
                    "degrees_in_x_epsilon_eta": [4, 2, 1],
                    "beta": "GF*r", "Q": "JGF", "B": "B0-B1*x", "U": "U0-U1*x-U2*x^2",
                    "constants": {name: str(value) for name, value in (
                        ("B0", B0), ("B1", B1), ("U0", U0), ("U1", U1), ("U2", U2),
                        ("D0", D0), ("D1", D1), ("D2", D2), ("D3", D3),
                        ("T0", T0), ("T1", T1), ("M1", M1))},
                    "h_ascending": ["U0^2-B0D0", "-2U0U1+B0D1+B1D0", "U1^2-2U0U2+B0D2-B1D1",
                                    "2U1U2+B0D3-B1D2", "U2^2-B1D3"],
                    "E_ascending": ["-2beta*U0+B0T0", "2beta*U1+B0M1-B1T0", "2beta*U2-B1M1"]},
        "zero_slack_propagation": {
            "dependency": "Proved report (55i): R0>0 at actual F,G,J,K, including m=4,5; not rerun",
            "identity": "B^2-AC=R0+(rR1+2AGFr*x0/K)s+(c2*r^2+AGFr/K)s^2+AGFr*epsilon, s=x-x0>=0",
            "positive_coefficients": "R1>0,c2>0; c1>=G/3>0 from actual last-passage equation",
            "epsilon_star": "0 for every formal x>=x0=K^2/(4p)",
            "monotonicity": "On B>0, H_eta=BQ>0 and H_epsilon=B*Gamma-2beta*C>0; use Fa(1-G)<p, B<=1-G, 2J<=3A",
            "conclusion": "H_actual>=H0(x_actual); equality at m=4, strict at m>=5",
            "intermediate_jets_need_not_be_physical": True},
        "moment_bounds": {
            "analytic_normalization": "Killed reversible diagonal density 4p*sum(c_j exp(-lambda_j t)); omega_j=4p*c_j/lambda_j sum to one",
            "mixture": "R=1/(1+lambda_j); (K,x,v,w)=E_omega[(1-R)(1,R,R^2,R^3)]",
            "bounds": ["x0<=x<=K(1-K)", "0<=epsilon<=x[1-x/(K(1-K))]", "eta>=epsilon^2/x"],
            "domain": "0<B(x)<=A^2 and U(x)>0 at actual x; not asserted on every formal jet"},
        "unbounded_shift_obstruction": {"fixed_m_q2_coefficient_of_h4": str(h4lead),
                                        "scope": "Every fixed m>=4, q->0; negative for small positive q"},
        "c2_endpoint_obstruction": {
            "scaling": "p_m=1-2/m, m->infinity; E=exp(2), formal x=x0=K^2/(4p), epsilon=eta=0",
            "limit_justification": "Unnormalized return density bounded by exp(6) for m>=4; dominated convergence of fixed moments",
            "limit_transforms": dict(zip(("F", "G", "J", "K"), map(str, (ff, gg, jj, kk)))),
            "N_E": str(numerator), "positive_denominator": str(denominator),
            "interval": ["369/50", "37/5"], "negative_interval_certificate": negative,
            "exp2_strict_bounds": [str(lower), str(upper)],
            "exponential_upper_tail_reason": "From degree seven onward successive term ratios <=1/4, with later strict inequality",
            "domain_certificates": domain, "sixth_passage_interval": "8/11<K<4/5",
            "finite_conclusion": "For all sufficiently large integer m: H0_m(x0_m)<0, x0_m<K_m(1-K_m), 0<B_m(x0_m)<A_m^2, U_m(x0_m)>0",
            "explicit_finite_threshold": None, "actual_first_derivative_substituted": False},
        "fifth_passage_floor": {"formula": "L=(1+4p+c-4p/K)/c, x>=x5=x0[1+cL^2/(5p)]>x0, c=q(m-4)>0",
                                "status": "Analytic first-step derivative bound for m>=5; refined quartic sign not proved or tested",
                                "limit_L": str(limiting_L), "limit_ratio": "x5/x0>49/40"},
        "third_gap_proved": False, "actual_LC_counterexample": False,
        "finite_time_LC_or_stationary_ULC_proved": False, "parameter_sweep": False}


def survival_third_fifth_floors_certificate():
    """Quantitative jet floors repair one exact simultaneous-limit discriminator.

    The generic derivative floors and interval width hold for m>=5. The
    actual limits p=1-c/m for c=2,3 are sign-certified with fifth and sixth
    passage floors respectively, and m=4..7 are checked at every rate.
    No all-size/all-rate theorem is claimed.
    """
    import hashlib
    from math import comb
    from sympy.polys.fields import field
    s = sp

    def values_from_fifth(p, c, L):
        """c=q(m-4); the formulas extend continuously to p=1,c>0."""
        q=1-p
        K=4*p/(1+4*p+c*(1-L))
        J=3*p/(1+3*p+(c+q)*(1-K))
        G=2*p/(1+2*p+(c+2*q)*(1-J))
        F=p/(1+p+(c+3*q)*(1-G))
        x5=K*K/(4*p)*(1+c*L*L/(5*p))
        eps5=c*K*K*L**3/(100*p**3)
        eta5=c*K*K*L**4/(500*p**4)
        return dict(F=F,G=G,J=J,K=K,p=p,a=c+q,b=c+2*q,
                    x5=x5,xplus=K*(1-K),eps5=eps5,eta5=eta5)

    def queue_gap(data, x, eps, eta):
        F,G,J,K,p,a,b=[data[k] for k in ('F','G','J','K','p','a','b')]
        y=J*J/(3*p)*(1+a*x)
        v=x*x/K+eps
        w=x**3/K**2+2*x*eps/K+eta
        d=G*G/(2*p)*(1+b*y)
        z2=y*y/J+a*J*J/(3*p)*v
        B=1-G-F*d
        C=1-J-y*G-y*d*F-z2*G*F
        D=1-K-x*(1-C)-v*(J*(1-B)+y*G*F)-w*J*G*F
        return C*C-B*D

    def quartic_coefficients(data):
        """Ascending coefficients of H(x,eps5,eta5), without expansion."""
        F,G,J,K,p,a,b=[data[k] for k in ('F','G','J','K','p','a','b')]
        eps,eta=data['eps5'],data['eta5']
        alpha=J*J/(3*p); r=a*alpha
        d0=G*G/(2*p); d1=b*d0; beta=G*F*r; Q=J*G*F
        B0=1-G-F*(d0+d1*alpha); B1=F*d1*r
        U0=1-J-(G+F*d0)*alpha-F*(d1+G/J)*alpha**2
        U1=r*(G+F*d0+2*F*(d1+G/J)*alpha)
        U2=r*r*F*(d1+G/J)+beta/K
        T0=J*(1-B0)+G*F*alpha; T1=J*B1+beta
        D0=1-K; D1=1-U0; D2=U1+T0/K; D3=U2+T1/K+Q/K**2
        M1=2*beta+J*B1+2*Q/K
        h=[U0*U0-B0*D0,-2*U0*U1+B0*D1+B1*D0,
           U1*U1-2*U0*U2+B0*D2-B1*D1,
           2*U1*U2+B0*D3-B1*D2,U2*U2-B1*D3]
        e=[-2*beta*U0+B0*T0,2*beta*U1+B0*M1-B1*T0,
           2*beta*U2-B1*M1]
        return [h[0]+eps*e[0]+beta*beta*eps*eps+B0*Q*eta,
                h[1]+eps*e[1]-B1*Q*eta,h[2]+eps*e[2],h[3],h[4]]

    def interval_coefficients(g, lower, upper):
        """Coefficients of (1+v)^4 P((lower+v*upper)/(1+v))."""
        return [sum(g[i]*comb(i,j)*lower**(i-j)*upper**j*comb(4-i,k-j)
                    for i in range(5)
                    for j in range(max(0,k-(4-i)),min(i,k)+1))
                for k in range(5)]

    def actual_small_bases():
        """Full passage jets at m=4..7, independently checked by killed resolvents."""
        Q,t=field('t',sp.QQ);p=t/(1+t);q=1/(1+t);out={}
        for m in (4,5,6,7):
         jets={m+1:[Q.zero]*4}
         for r in range(m,0,-1):
          down=(m-r)*q;later=jets[r+1]
          den=[1+r*p+down*(1-later[0]),-1-down*later[1],-down*later[2],-down*later[3]]
          row=[r*p/den[0]]
          for j in range(1,4):row.append(-sum(den[i]*row[j-i] for i in range(1,j+1))/den[0])
          jets[r]=row
         F=jets[1][0];G,d=jets[2][:2];J,y,zeta=jets[3][:3];K,x,v,w=jets[4]
         A=1-F;B=1-G-F*d;C=1-J-y*G-y*d*F-zeta*G*F
         D=1-K-x*(1-C)-v*(J*(1-B)+y*G*F)-w*J*G*F
         gap=C*C-B*D
         num,den=gap.numer,gap.denom
         assert all(a>0 for a in num.coeffs()) and all(a>0 for a in den.coeffs())
         rows=[[*mon,int(co)] for mon,co in sorted(num.items())]
         assert all(sp.Integer(int(co))==co for co in num.coeffs())
         probe=sp.Rational(2,3) # p=2/5
         def at(expr):return sp.cancel(expr.as_expr().subs(sp.Symbol('t'),probe))
         diagonal,_=boundary_profiles(m,sp.Rational(2,5))
         assert tuple(at(z) for z in (A,B,C,D))==tuple(diagonal[m-j] for j in (1,2,3,4))
         assert at(gap)==diagonal[m-3]**2-diagonal[m-2]*diagonal[m-4]>0
         out[str(m)]={'numerator_ascending_coefficients':[int(num.get((j,),0)) for j in range(num.degree()+1)],'denominator_ascending_coefficients':[int(den.get((j,),0)) for j in range(den.degree()+1)],'numerator_degree':num.degree(),'positive_numerator_terms':len(rows),'positive_denominator_terms':len(den),'numerator_sha256':hashlib.sha256(json.dumps(rows,separators=(',',':')).encode()).hexdigest(),'rational_gap_at_p_2_over_5':str(at(gap)),'all_four_queue_values_match_independent_boundary_resolvents':True}
        return out

    def sixth_floors(p,c,d,M):
        L=5*p/(1+5*p+d*(1-M))
        data=values_from_fifth(p,c,L)
        K=data['K']
        s1=M*M/(6*p);s2=M**3/(6*p)**2;s3=M**4/(6*p)**3
        ell1=L*L/(5*p)*(1+d*s1)
        eL=d*L*L/(5*p)*s2;hL=d*L*L/(5*p)*s3
        ell2=ell1*ell1/L+eL
        ell3=ell1**3/L**2+2*ell1*eL/L+hL
        xi=K*K/(4*p)*(1+c*ell1)
        eps=c*K*K/(4*p)*ell2;eta=c*K*K/(4*p)*ell3
        return data,{'M':M,'L':L,'ell1':ell1,'ell2':ell2,'ell3':ell3,
                     'xi6':xi,'epsilon6':eps,'eta6':eta}

    def actual_c3_refinement():
        Q,E,u=field('E,u',s.QQ)
        sums=[sum(s.Rational(3)**j/s.factorial(j) for j in range(r+1)) for r in range(7)]
        ratios={r:Q(s.Rational(r,3))*(E-sums[r])/(E-sums[r-1]) for r in range(1,7)}
        data,jets=sixth_floors(Q.one,3*Q.one,3*Q.one,ratios[6])
        assert jets['L']==ratios[5]
        assert all(data[k]==ratios[r] for k,r in [('F',1),('G',2),('J',3),('K',4)])
        EE,uu,z=s.symbols('E u z')
        low=sum(s.Rational(3)**j/s.factorial(j) for j in range(19))
        high=low+s.Rational(3)**19/s.factorial(19)/(1-s.Rational(3,20))
        assert low==s.Rational(19599925833041,975822848000)
        assert high==s.Rational(6330776047615183,315190779904000)
        assert low>20 and high<s.Rational(1291,64)

        def sign_certificate(poly):
            poly=s.Poly(poly,EE);degree=poly.degree()
            transformed=s.Poly(sum(coefficient*(low+high*z)**j*(1+z)**(degree-j)
                                  for (j,),coefficient in poly.terms()),z)
            multiplier,integer=transformed.clear_denoms()
            coefficients=[int(integer.nth(j)) for j in range(degree+1)]
            signs=sorted(set(1 if c>0 else -1 if c<0 else 0 for c in coefficients))
            return {'degree':degree,'clearing_multiplier':str(multiplier),
                    'ascending_coefficients':coefficients,'coefficient_signs':signs}

        def evaluate_package(lower,eps,eta,label):
            x=(lower+data['xplus']*u)/(1+u)
            H=queue_gap(data,x,eps,eta)*(1+u)**4
            N=s.Poly(H.numer.as_expr(),EE,uu)
            den=s.factor(H.denom.as_expr())
            shifted_den=s.Poly(H.denom.as_expr().subs(EE,20+z),z)
            assert all(c>0 for c in shifted_den.coeffs())
            assert N.degree_list()==(9,4) and len(N.terms())==50
            certificates=[]
            for k in range(5):
                poly=s.Poly(sum(c*EE**i for (i,j),c in N.terms() if j==k),EE)
                certificates.append({'u_degree':k,'numerator_ascending_coefficients':
                                      [int(poly.nth(i)) for i in range(10)],
                                      'interval_certificate':sign_certificate(poly)})
            rows=[[*powers,int(coefficient)] for powers,coefficient in sorted(N.terms())]
            return {'label':label,'lower_endpoint':str(lower.as_expr()),
                    'epsilon_floor':str(eps.as_expr()),'eta_floor':str(eta.as_expr()),
                    'positive_denominator':str(den),'numerator_degrees':[9,4],
                    'numerator_terms':50,'numerator_coefficient_rows':rows,
                    'numerator_sha256':hashlib.sha256(json.dumps(rows,separators=(',',':')).encode()).hexdigest(),
                    'Bernstein_certificates':certificates},H

        old,oldH=evaluate_package(data['x5'],data['eps5'],data['eta5'],'old fifth floors at actual L')
        assert [row['interval_certificate']['coefficient_signs'] for row in old['Bernstein_certificates']]==[[-1],[1],[1],[1],[1]]
        new,newH=evaluate_package(jets['xi6'],jets['epsilon6'],jets['eta6'],'sixth floors at actual L,M')
        assert all(row['interval_certificate']['coefficient_signs']==[1] for row in new['Bernstein_certificates'])
        coefficient_data=dict(data)
        coefficient_data.update(eps5=jets['epsilon6'],eta5=jets['eta6'])
        structured=interval_coefficients(quartic_coefficients(coefficient_data),
                                                 jets['xi6'],data['xplus'])
        assert newH==sum(structured[k]*u**k for k in range(5))

        # Exact width identity and physical M upper bound, checked without a grid.
        L,M,K=jets['L'],jets['M'],data['K']
        assert data['xplus']-jets['xi6']==3*K*K/4*((1-L)-jets['ell1'])
        assert (1-L)-jets['ell1']==L/5*((1-L)+3*((1-M)-L*M*M/6))
        assert Q(s.Rational(6,7))-M==5*(1291-64*E)/(56*(5*E-92))

        # General symbolic identity used in the width proof, not a finite case.
        p,d,ms=s.symbols('p d ms')
        ls=5*p/(1+5*p+d*(1-ms))
        ell1=ls*ls/(5*p)*(1+d*ms*ms/(6*p))
        assert s.cancel((1-ls)-ell1-ls/(5*p)*((1-ls)+d*((1-ms)-ls*ms*ms/(6*p))))==0

        F,G,J=[data[k] for k in ('F','G','J')]
        yplus=J*J/3*(1+3*data['xplus'])
        Bplus=1-G-F*G*G/2*(1+3*yplus)
        Bcert=sign_certificate(Bplus.numer.as_expr())
        assert Bcert['coefficient_signs']==[1]
        assert all(c>0 for c in s.Poly(Bplus.denom.as_expr().subs(EE,20+z),z).coeffs())

        record={
          'status':'PASS','scope':'only actual c=3 limit, p=1-3/m',
          'moment_ratio_formulas':{str(r):str(value.as_expr()) for r,value in ratios.items()},
          'exp3_lower_bound':str(low),'exp3_upper_bound':str(high),
          'old_fifth_floor_package':old,'new_sixth_floor_package':new,
          'old_floor_actual_L_lower_endpoint_strictly_negative':True,
          'new_floor_all_five_Bernstein_coefficients_strictly_positive':True,
          'independent_structured_quartic_identity':True,
          'width_identities_and_physical_M_range':True,
          'B_at_upper_endpoint_positive_certificate':Bcert,
          'B_at_upper_endpoint':str(Bplus.as_expr()),
          'finite_M_coordinate_depth8_candidate':'all-parameter signs OPEN; formulas in report',
          'finite_consequence':'Actual third survival gap positive for all sufficiently large m along p=1-3/m; no explicit threshold',
          'generic_sixth_floor_scope':'m>=6,0<p<1; positive reciprocal series propagates the first three sixth-passage jet floors; exact at m=6',
          'all_size_all_rate_gap_claim':False,'actual_gap_counterexample':False,
          'new_parameter_sweep':False,
        }
        return record

    pp,cc,LL=s.symbols('p c L',positive=True)
    KK=4*pp/(1+4*pp+cc*(1-LL))
    width=KK*(1-KK)-KK**2/(4*pp)*(1+cc*LL**2/(5*pp))
    assert s.cancel(width-cc*KK**2/(4*pp)*(1-LL-LL**2/(5*pp)))==0
    assert s.cancel((1-LL-LL**2/(5*pp)).subs(LL,5*pp/(1+5*pp))-1/(1+5*pp)**2)==0
    Q,E,u=field('E,u',s.QQ)
    L=(15*E-109)/(6*(E-7))
    data=values_from_fifth(Q.one,2*Q.one,L)
    moment_values={'F':(E-3)/(2*(E-1)),'G':(E-5)/(E-3),
                   'J':(3*E-19)/(2*(E-5)),'K':6*(E-7)/(3*E-19)}
    assert all(data[k]==v for k,v in moment_values.items())
    x5,xplus=data['x5'],data['xplus']
    K=data['K']
    assert xplus-x5==K*K/2*((1-L)-L*L/5)

    EE,uu,zz=s.symbols('E u zz')
    low=sum(s.Rational(2)**j/s.factorial(j) for j in range(13))
    high=low+s.Rational(2)**13/s.factorial(13)/(1-s.Rational(1,7))
    assert low==s.Rational(691283,93555) and high==s.Rational(134800213,18243225)
    assert low>7 and high<s.Rational(37,5)

    def sign_certificate(poly):
        poly=s.Poly(poly,EE); degree=poly.degree()
        transformed=s.Poly(sum(coefficient*(low+high*zz)**j*(1+zz)**(degree-j)
                              for (j,),coefficient in poly.terms()),zz)
        multiplier,integer_poly=transformed.clear_denoms()
        return {'degree':degree,'clearing_multiplier':str(multiplier),
                'ascending_coefficients':[int(integer_poly.nth(j)) for j in range(degree+1)]}

    bad=queue_gap(data,x5,0,0)
    bad_certificate=sign_certificate(bad.numer.as_expr())
    assert all(c<0 for c in bad_certificate['ascending_coefficients'])

    x=(x5+xplus*u)/(1+u)
    good=queue_gap(data,x,data['eps5'],data['eta5'])*(1+u)**4
    den=s.factor(good.denom.as_expr())
    expected_den=51840000*(EE-7)**2*(EE-5)**2*(EE-3)**2*(EE-1)**2*(3*EE-19)
    assert s.expand(den-expected_den)==0
    assert s.expand(bad.denom.as_expr()-expected_den)==0
    numerator=s.Poly(good.numer.as_expr(),EE,uu)
    assert numerator.degree_list()==(9,4) and len(numerator.terms())==50
    structured=interval_coefficients(quartic_coefficients(data),x5,xplus)
    assert good==sum(structured[k]*u**k for k in range(5))
    bernstein=[]
    for k in range(5):
        poly=s.Poly(sum(coefficient*EE**i for (i,j),coefficient in numerator.terms() if j==k),EE)
        cert=sign_certificate(poly)
        assert cert['degree']==9 and all(c>0 for c in cert['ascending_coefficients'])
        bernstein.append({'u_degree':k,'numerator_ascending_coefficients':
                          [int(poly.nth(i)) for i in range(10)],'positive_interval_certificate':cert})

    # B decreases with x. Its upper-endpoint value is positive on the same
    # single E interval, so the whole tested interval lies in B>0.
    F,G,J=[data[k] for k in ('F','G','J')]
    y=J*J/3*(1+2*xplus)
    Bplus=1-G-F*G*G/2*(1+2*y)
    Bcert=sign_certificate(Bplus.numer.as_expr())
    assert all(c>0 for c in Bcert['ascending_coefficients'])
    Bden_shift=s.Poly(Bplus.denom.as_expr().subs(EE,7+zz),zz)
    assert all(coefficient>0 for coefficient in Bden_shift.coeffs())
    rows=[[*powers,int(coefficient)] for powers,coefficient in sorted(numerator.terms())]
    result={
      'status':'positive c=2 limiting interval; all-size/all-rate third gap OPEN','tested_regime':'only the actual c=2 limit p=1-2/m',
      'fifth_value':str(L.as_expr()),'x5':str(x5.as_expr()),'xplus':str(xplus.as_expr()),
      'eps5':str(data['eps5'].as_expr()),'eta5':str(data['eta5'].as_expr()),
      'exp2_lower_bound':str(low),'exp2_upper_bound':str(high),
      'H0_x5_negative_numerator':str(bad.numer.as_expr()),
      'H0_x5_negative_certificate':bad_certificate,
      'common_positive_denominator':str(den),
      'restored_floors_numerator_degrees':[9,4],
      'restored_floors_numerator_terms':50,
      'restored_floors_ascending_coefficient_rows':rows,
      'restored_floors_coefficient_sha256':hashlib.sha256(json.dumps(rows,separators=(',',':')).encode()).hexdigest(),
      'positive_Bernstein_certificates':bernstein,
      'B_at_upper_endpoint':str(Bplus.as_expr()),'B_positive_certificate':Bcert,
      'B_denominator_positive_after_shift_E7':True,
      'five_structured_interval_coefficients_agree_with_direct_queue_gap':True,
      'strict_positive_full_c2_interval_with_eps5_eta5':True,
      'finite_all_size_signs':'OPEN; no finite-parameter expansion performed',
      'finite_consequence':'For every sufficiently large m, the actual third gap is positive along p=1-2/m; no explicit threshold',
      'fifth_jet_floors':{'domain':'m>=5,0<p<1','proof':'Reciprocal first-step series has nonnegative coefficients: ell_j>=L^(j+1)/(5p)^j',
          'formulas':['x>=K^2/(4p)*(1+c*L^2/(5p))','epsilon>=c*K^2*L^3/(100*p^3)','eta>=c*K^2*L^4/(500*p^4)'],
          'comparison':'H_actual>=H(x_actual,epsilon5,eta5) by the audited higher-slack monotonicity'},
      'all_size_all_rate_third_gap_claim':False,'new_parameter_sweep':False,
    }
    result['small_actual_bases']=actual_small_bases()
    result['sixth_passage_c3']=actual_c3_refinement()
    lower=s.S.Zero
    for birth in (8,7,6,5):
        lower=s.Rational(birth)/(1+birth+3*(1-lower))
    assert lower==s.Rational(115,153)
    formal=values_from_fifth(s.S.One,s.Integer(3),lower)
    formal_gap=s.cancel(queue_gap(formal,formal['x5'],formal['eps5'],formal['eta5']))
    assert formal_gap==s.Rational(-31670939693550051419413,29841111934881500548931031)<0
    assert formal['x5']<formal['xplus']
    result['depth8_relaxation_obstruction']={
        'path':'t->infinity,w=3t,u=0; p=t/(1+t),m=w+8',
        'L_lower_limit':str(lower),'gap_limit':str(formal_gap),
        'values':{key:str(formal[key]) for key in ('F','G','J','K','x5','xplus')},
        'conclusion':'Depth8 formal lower endpoint gives a negative restored-floor quartic; five global Bernstein signs on that enlarged domain are false',
        'actual_L_used':False,'actual_survival_counterexample':False}
    return result


def survival_third_lc_certificate():
    """All-size third infinite-horizon survival gap, boundary rate one.

    Sixth-passage derivative floors reduce the actual gap to a quartic
    on its physical first-jet interval. Retaining passages eight through
    six bounds the last free transform. All five Bernstein numerators
    then have positive coefficients in t>0,w,u>=0 (m=w+8,p=t/(1+t)).
    Exact actual bases m=4..7 complete the theorem. This does not assert
    full diagonal LC, finite-time LC, stationary ULC, or NA.

    Only this selected helper imports python-flint. Integer polynomials
    are staged before substitution; output retains hashes and statistics,
    not the 92,592 reproducible coefficient rows. Earlier helpers are
    neither modified nor rerun.
    """
    import hashlib
    import resource
    import time
    from contextlib import contextmanager
    from math import comb
    from sympy.polys.fields import field
    try:
        from flint import fmpz_mpoly_ctx
    except ImportError as error:
        raise SystemExit("--survival-third-lc requires optional python-flint") from error

    @contextmanager
    def bounded_memory():
        previous = resource.getrlimit(resource.RLIMIT_AS)
        limits = [1900*1024*1024]+[value for value in previous if value >= 0]
        limit = min(limits)
        resource.setrlimit(resource.RLIMIT_AS, (limit, previous[1]))
        try:
            yield limit
        finally:
            resource.setrlimit(resource.RLIMIT_AS, previous)

    def coefficient_rows(poly):
        return [[*map(int, poly.monomial(j)), int(poly.coefficient(j))]
                for j in range(len(poly)-1, -1, -1)]

    def digest(rows):
        return hashlib.sha256(json.dumps(rows, separators=(",", ":")).encode()).hexdigest()

    def polynomial_stats(poly):
        return {"terms": len(poly), "degrees": list(map(int, poly.degrees())),
                "degree_tw": max(sum(map(int, poly.monomial(j)[:2]))
                                 for j in range(len(poly))),
                "total_degree": int(poly.total_degree())}

    def actual_small_bases():
        """Only the four existing actual bases; no old limit audits."""
        ring, t = field("t", sp.QQ)
        p, q = t/(1+t), 1/(1+t)
        result = {}
        for m, expected_terms in zip((4, 5, 6, 7), (19, 26, 33, 40)):
            jets = {m+1: [ring.zero]*4}
            for r in range(m, 0, -1):
                down, later = (m-r)*q, jets[r+1]
                denominator = [1+r*p+down*(1-later[0]), -1-down*later[1],
                               -down*later[2], -down*later[3]]
                row = [r*p/denominator[0]]
                for j in range(1, 4):
                    row.append(-sum(denominator[i]*row[j-i]
                                    for i in range(1, j+1))/denominator[0])
                jets[r] = row
            F = jets[1][0]
            G, dG = jets[2][:2]
            J, y, zeta = jets[3][:3]
            K, x, second, third = jets[4]
            A = 1-F
            B = 1-G-F*dG
            C = 1-J-y*G-y*dG*F-zeta*G*F
            D = 1-K-x*(1-C)-second*(J*(1-B)+y*G*F)-third*J*G*F
            gap = C*C-B*D
            numerator, denominator = gap.numer, gap.denom
            assert all(value > 0 for value in numerator.coeffs())
            assert all(value > 0 for value in denominator.coeffs())
            assert len(numerator) == numerator.degree() == expected_terms
            assert numerator.get((1,), 0) > 0 and denominator.get((0,), 0) > 0
            num_rows = [[*powers, int(value)] for powers, value in sorted(numerator.items())]
            den_rows = [[*powers, int(value)] for powers, value in sorted(denominator.items())]
            assert all(sp.Integer(int(value)) == value
                       for value in numerator.coeffs()+denominator.coeffs())

            def at_probe(value):
                return sp.cancel(value.as_expr().subs(sp.Symbol("t"), sp.Rational(2, 3)))

            diagonal, _ = boundary_profiles(m, sp.Rational(2, 5))
            assert tuple(at_probe(value) for value in (A, B, C, D)) == tuple(
                diagonal[m-j] for j in (1, 2, 3, 4))
            assert at_probe(gap) == diagonal[m-3]**2-diagonal[m-2]*diagonal[m-4] > 0
            result[str(m)] = {
                "numerator_degree": numerator.degree(),
                "positive_numerator_terms": len(numerator),
                "positive_denominator_terms": len(denominator),
                "numerator_sha256": digest(num_rows),
                "denominator_sha256": digest(den_rows),
                "strict_t_anchor": int(numerator.get((1,), 0)),
                "positive_denominator_constant": int(denominator.get((0,), 0)),
                "rational_gap_at_p_2_over_5": str(at_probe(gap)),
                "all_four_queue_values_match_independent_boundary_resolvents": True,
            }
        return result

    started = time.monotonic()
    with bounded_memory() as memory_limit:
        ctx = fmpz_mpoly_ctx.get(["t", "w", "M"])
        t, w, M = ctx.gens()
        one, zero = ctx.constant(1), ctx.constant(0)

        class Rat:
            """Reduced fractions, cancelling before polynomial products."""
            def __init__(self, numerator, denominator=None, reduced=False):
                if isinstance(numerator, Rat):
                    self.n, self.d = numerator.n, numerator.d
                    return
                numerator = numerator if hasattr(numerator, "gcd") else ctx.constant(numerator)
                denominator = (one if denominator is None else denominator
                               if hasattr(denominator, "gcd") else ctx.constant(denominator))
                if denominator.is_zero():
                    raise ZeroDivisionError
                if numerator.is_zero():
                    self.n, self.d = zero, one
                    return
                if not reduced:
                    common = numerator.gcd(denominator)
                    numerator, denominator = numerator//common, denominator//common
                if denominator.coefficient(0) < 0:
                    numerator, denominator = -numerator, -denominator
                self.n, self.d = numerator, denominator

            def __add__(self, other):
                other = Rat(other)
                common = self.d.gcd(other.d)
                first, second = self.d//common, other.d//common
                numerator = self.n*second+other.n*first
                cancel = numerator.gcd(common)
                return Rat(numerator//cancel, first*(other.d//cancel), True)

            __radd__ = __add__

            def __neg__(self):
                return Rat(-self.n, self.d, True)

            def __sub__(self, other):
                return self+-Rat(other)

            def __rsub__(self, other):
                return Rat(other)+-self

            def __mul__(self, other):
                other = Rat(other)
                first, second = self.n.gcd(other.d), other.n.gcd(self.d)
                return Rat((self.n//first)*(other.n//second),
                           (self.d//second)*(other.d//first), True)

            __rmul__ = __mul__

            def __truediv__(self, other):
                other = Rat(other)
                if other.n.is_zero():
                    raise ZeroDivisionError
                first, second = self.n.gcd(other.n), self.d.gcd(other.d)
                return Rat((self.n//first)*(other.d//second),
                           (self.d//second)*(other.n//first), True)

            def __rtruediv__(self, other):
                return Rat(other)/self

            def __pow__(self, power):
                if power < 0:
                    return (Rat(1)/self)**(-power)
                return Rat(self.n**power, self.d**power, True)

            def __eq__(self, other):
                other = Rat(other)
                return self.n*other.d == other.n*self.d

        # Generic M coordinate, before the depth-eight interval substitution.
        kappa = {6: one, 5: 4+6*t+w-(w+3)*M}
        for birth in range(4, 0, -1):
            kappa[birth] = ((1+(birth+1)*t+w+8-birth)*kappa[birth+1]
                            -(w+8-birth)*(birth+1)*t*kappa[birth+2])
        p, q = Rat(t, 1+t), Rat(1, 1+t)
        c, down5 = Rat(w+4, 1+t), Rat(w+3, 1+t)
        F, G = Rat(t*kappa[2], kappa[1]), Rat(2*t*kappa[3], kappa[2])
        J, K = Rat(3*t*kappa[4], kappa[3]), Rat(4*t*kappa[5], kappa[4])
        L, MR = Rat(5*t, kappa[5]), Rat(M)
        ell1 = L**2/(5*p)*(1+down5*MR**2/(6*p))
        eL = down5*L**2/(5*p)*MR**3/(6*p)**2
        hL = down5*L**2/(5*p)*MR**4/(6*p)**3
        ell2 = ell1**2/L+eL
        ell3 = ell1**3/L**2+2*ell1*eL/L+hL
        lower, upper = K**2/(4*p)*(1+c*ell1), K*(1-K)
        epsilon, eta = c*K**2/(4*p)*ell2, c*K**2/(4*p)*ell3
        assert upper-lower == c*K**2/(4*p)*((1-L)-ell1)
        assert (1-L)-ell1 == L/(5*p)*((1-L)+down5*((1-MR)-L*MR**2/(6*p)))

        # Structured ascending quartic coefficients of H(X,epsilon6,eta6).
        alpha = J**2/(3*p)
        r = (c+q)*alpha
        d0, d1 = G**2/(2*p), (c+2*q)*G**2/(2*p)
        beta, Q = G*F*r, J*G*F
        B0, B1 = 1-G-F*(d0+d1*alpha), F*d1*r
        U0 = 1-J-(G+F*d0)*alpha-F*(d1+G/J)*alpha**2
        U1 = r*(G+F*d0+2*F*(d1+G/J)*alpha)
        U2 = r*r*F*(d1+G/J)+beta/K
        T0, T1 = J*(1-B0)+G*F*alpha, J*B1+beta
        D0, D1 = 1-K, 1-U0
        D2, D3 = U1+T0/K, U2+T1/K+Q/K**2
        gain = 2*beta+J*B1+2*Q/K
        h = [U0*U0-B0*D0, -2*U0*U1+B0*D1+B1*D0,
             U1*U1-2*U0*U2+B0*D2-B1*D1,
             2*U1*U2+B0*D3-B1*D2, U2*U2-B1*D3]
        e = [-2*beta*U0+B0*T0, 2*beta*U1+B0*gain-B1*T0,
             2*beta*U2-B1*gain]
        quartic = [h[0]+epsilon*e[0]+beta*beta*epsilon*epsilon+B0*Q*eta,
                   h[1]+epsilon*e[1]-B1*Q*eta, h[2]+epsilon*e[2], h[3], h[4]]

        # Independent polynomial arithmetic in X from the direct queue formula.
        def padd(first, second):
            return [(first[j] if j < len(first) else Rat(0))
                    +(second[j] if j < len(second) else Rat(0))
                    for j in range(max(len(first), len(second)))]

        def psub(first, second):
            return padd(first, [-value for value in second])

        def pmul(first, second):
            result = [Rat(0) for _ in range(len(first)+len(second)-1)]
            for i, value in enumerate(first):
                for j, other in enumerate(second):
                    result[i+j] = result[i+j]+value*other
            return result

        def pscale(poly, factor):
            return [value*factor for value in poly]

        Y, V = [alpha, r], [epsilon, Rat(0), 1/K]
        W = [eta, 2*epsilon/K, Rat(0), 1/K**2]
        GD = pscale(padd([Rat(1)], pscale(Y, c+2*q)), G**2/(2*p))
        Z = padd(pscale(pmul(Y, Y), 1/J), pscale(V, r))
        BP = psub([1-G], pscale(GD, F))
        CP = psub(psub(psub([1-J], pscale(Y, G)), pscale(pmul(Y, GD), F)),
                   pscale(Z, G*F))
        DP = psub(psub(psub([1-K], [Rat(0)]+psub([Rat(1)], CP)),
                        pmul(V, padd(pscale(psub([Rat(1)], BP), J), pscale(Y, G*F)))),
                   pscale(W, J*G*F))
        assert psub(pmul(CP, CP), pmul(BP, DP)) == quartic

        # Independent truncated continued fraction in h=1-z. Bottom M(h)
        # has the three exponential jet floors, with its constant M fixed.
        series = {7: [M, zero, zero, zero], 6: [6*t, -(1+t)*M, zero, zero]}
        for birth in range(5, 0, -1):
            later, earlier = series[birth+1], series[birth+2]
            series[birth] = [
                (1+(birth+1)*t+w+8-birth)*later[j]
                -((1+t)*later[j-1] if j else zero)
                -(w+8-birth)*(birth+1)*t*earlier[j] for j in range(4)]
            assert series[birth][0] == 6*t*kappa[birth]

        def passage_jets(birth, degree):
            denominator = series[birth]
            numerator = [birth*t*value for value in series[birth+1]]
            numerators, result = [], []
            for j in range(degree+1):
                value = numerator[j]*denominator[0]**j
                for ell in range(1, j+1):
                    value -= denominator[ell]*numerators[j-ell]*denominator[0]**(ell-1)
                numerators.append(value)
                result.append(Rat(value, denominator[0]**(j+1)))
            return result

        lower_y = alpha+r*lower
        lower_v = lower*lower/K+epsilon
        lower_w = lower**3/K**2+2*lower*epsilon/K+eta
        lower_z = lower_y**2/J+r*lower_v
        lower_d = G**2/(2*p)*(1+(c+2*q)*lower_y)
        assert passage_jets(4, 3) == [K, lower, lower_v, lower_w]
        assert passage_jets(3, 2) == [J, lower_y, lower_z]
        assert passage_jets(2, 1) == [G, lower_d]

        outctx = fmpz_mpoly_ctx.get(["t", "w", "u"])
        tt, ww, uu = outctx.gens()
        unit = outctx.constant(1)
        b8 = 1+9*tt+ww
        b7 = (2+8*tt+ww)*b8-8*tt*(ww+1)
        b6 = (3+7*tt+ww)*b7-7*tt*(ww+2)*b8
        NM = 6*tt*((1+7*tt)*b7+uu*b6)
        DM = (1+uu)*(1+7*tt)*b6
        # Separate numerator/denominator recurrence audits the lower bound.
        low_n, low_d = outctx.constant(0), unit
        for birth in (8, 7, 6):
            low_n, low_d = (birth*tt*low_d,
                           (1+(birth+1)*tt+ww+8-birth)*low_d-(ww+8-birth)*low_n)
        assert NM == low_n*(1+7*tt)+uu*6*tt*low_d
        assert DM == (1+uu)*(1+7*tt)*low_d
        for poly in (b8, b7, b6, NM, DM):
            assert all(poly.coefficient(j) > 0 for j in range(len(poly)))
        numerator_powers, denominator_powers = [unit], [unit]

        def homogeneous_substitution(rows, degree):
            while len(numerator_powers) <= degree:
                numerator_powers.append(numerator_powers[-1]*NM)
                denominator_powers.append(denominator_powers[-1]*DM)
            blocks = [{} for _ in range(degree+1)]
            for i, j, k, coefficient in rows:
                blocks[k][(i, j, 0)] = coefficient
            result = outctx.constant(0)
            for k, block in enumerate(blocks):
                result += outctx.from_dict(block)*numerator_powers[k]*denominator_powers[degree-k]
            # Independent Horner accumulation verifies the homogeneous map.
            horner, power = outctx.from_dict(blocks[degree]), unit
            for k in range(degree-1, -1, -1):
                power *= DM
                horner = horner*NM+outctx.from_dict(blocks[k])*power
            assert result == horner
            return result

        continuants = []
        for birth, constant in zip(range(1, 6), (40320, 5040, 720, 120, 24)):
            factor = homogeneous_substitution(coefficient_rows(kappa[birth]), 1)
            rows = coefficient_rows(factor)
            assert all(row[-1] > 0 for row in rows)
            assert factor(0, 0, 0) == constant
            continuants.append({"index": birth, **polynomial_stats(factor),
                                "constant": constant, "positive_terms": len(factor),
                                "negative_terms": 0, "coefficient_sha256": digest(rows)})

        expected_hashes = [
            "11469380c0c03823b735c6ab684af669f0b2bffd82553e0165c50a776c6640b7",
            "648aadf78a1f13abab77a3f1dc5a4ccfdcc10f68bf037d04a2a3c2f9f1387c4b",
            "cacb4a482f80683523c1c46e43d32dc18b51064a8079f1c1be45aa2a3002c854",
            "7ccb56fbd46cfb5f13be810106cc85476174d65d08e29e6b36214e361db35152",
            "ed04e9c10948e1018a7ac62f1e7c77fb968e3fc426aabdbfe86a2355638243a8",
        ]
        coefficients = []
        anchor0 = 1981728150031171584000000
        # The endpoint must pass before any other Bernstein index is built.
        for index in range(5):
            coefficient = sum(
                quartic[i]*comb(i, j)*lower**(i-j)*upper**j*comb(4-i, index-j)
                for i in range(5)
                for j in range(max(0, index-(4-i)), min(i, index)+1))
            remaining, powers = coefficient.d, []
            for birth in range(1, 6):
                power = 0
                while True:
                    quotient, remainder = divmod(remaining, kappa[birth])
                    if not remainder.is_zero():
                        break
                    remaining, power = quotient, power+1
                powers.append(power)
            assert remaining == 9
            assert powers == [2, 2, 2, 1, 0 if index == 0 else 2]
            degree = sum(powers)
            assert int(coefficient.n.degrees()[2]) == degree
            generic_rows = coefficient_rows(coefficient.n)
            generic_stats = polynomial_stats(coefficient.n)
            predicted_degree = generic_stats["degree_tw"]+4*degree
            dense_bound = (degree+1)*(predicted_degree+1)*(predicted_degree+2)//2
            numerator = homogeneous_substitution(generic_rows, degree)
            rows = coefficient_rows(numerator)
            assert len(rows) == (11592 if index == 0 else 20250)
            assert len(rows) <= dense_bound and all(row[-1] > 0 for row in rows)
            assert digest(rows) == expected_hashes[index]
            anchor = [row[-1] for row in rows if row[:3] == [1, 0, 0]]
            assert anchor == [anchor0*(1 if index == 0 else 576*comb(4, index))]
            coefficients.append({
                "index": index, "generic_numerator": generic_stats,
                "generic_numerator_sha256": digest(generic_rows),
                "denominator_constant": 9, "continuant_powers": powers,
                "numerator": polynomial_stats(numerator),
                "positive_terms": len(rows), "negative_terms": 0,
                "coefficient_sha256": digest(rows),
                "strict_t_anchor": anchor[0],
                "predicted_dense_monomial_bound": dense_bound,
                "independent_homogeneous_Horner_identity": True,
            })
            print("Third survival LC Bernstein", index, "passed:", len(rows),
                  "positive coefficients", flush=True)

        bases = actual_small_bases()
        assert sum(row["positive_terms"] for row in coefficients) == 92592
        result = {
            "status": "proved by sixth-passage floors and exact positive-coefficient certificates",
            "theorem": "a_(m-3)^2>a_(m-4)*a_(m-2) for every integer m>=4 and 0<p<1",
            "model": "actual infinite-horizon Ehrenfest survival diagonal; births p(m-y), deaths (1-p)y, boundary rate one",
            "all_size_all_rate_third_gap_claim": True,
            "full_survival_LC_claim": False, "finite_time_LC_claim": False,
            "stationary_ULC_or_NA_claim": False,
            "prior_input": "Strong second-gap certificate (55i) and its formal-jet propagation; higher-slack monotonicity (55j)",
            "parameterization": {"domain": "t>0,w,u>=0", "m": "w+8", "p": "t/(1+t)",
                "k8": "8t/(1+9t+w)", "k7": "7t/[1+8t+(w+1)(1-k8)]",
                "Mlo": "6t/[1+7t+(w+2)(1-k7)]", "Mhi": "6t/(1+7t)",
                "M": "(Mlo+u*Mhi)/(1+u)",
                "coverage": "Mlo<=actual P6(1)<Mhi by increasing first-step fractions; exact lower endpoint at m=8"},
            "sixth_floors": {
                "L": "5p/[1+5p+d(1-M)], d=(m-5)(1-p)",
                "ell1": "L^2/(5p)*(1+d*M^2/(6p))",
                "eL": "d*L^2*M^3/[(5p)(6p)^2]",
                "hL": "d*L^2*M^4/[(5p)(6p)^3]",
                "ell2": "ell1^2/L+eL", "ell3": "ell1^3/L^2+2*ell1*eL/L+hL",
                "xi6": "K^2/(4p)*(1+c*ell1), c=(m-4)(1-p)",
                "epsilon6": "c*K^2*ell2/(4p)", "eta6": "c*K^2*ell3/(4p)",
                "interval": "xi6<=xi_actual<=K(1-K)",
                "positive_width": "K(1-K)-xi6=c*K^2/(4p)*[(1-L)-ell1]>0",
                "width_proof": "(1-L)-ell1=L/(5p)*[(1-L)+d*((1-M)-L*M^2/(6p))]; M<=6p/(1+6p)",
                "comparison": "H_actual>=H(xi_actual,epsilon6,eta6) by positive higher-slack derivatives on the strong preceding-gap domain"},
            "Bernstein_convention": "[v^j](1+v)^4*P6((xi6+v*xplus)/(1+v)); normalized coefficient divides by binomial(4,j)",
            "strictness": "Each numerator has a positive t^1*w^0*u^0 anchor; denominators have positive constant terms; endpoints included",
            "positive_continuants": continuants,
            "Bernstein_coefficients": coefficients,
            "total_positive_terms": 92592, "total_negative_terms": 0,
            "exact_identity_checks": {
                "direct_queue_equals_structured_quartic": True,
                "independent_truncated_continued_fraction_jets": True,
                "both_physical_interval_width_identities": True,
                "independent_lower_continued_fraction_pair_construction": True,
                "all_five_homogeneous_Horner_identities": True,
                "all_five_audited_coefficient_hashes": True,
                "positive_continuant_denominators": True},
            "small_actual_bases": bases,
            "optional_dependency": "python-flint, imported locally by this helper only",
            "address_space_soft_cap_MiB": memory_limit//(1024*1024),
            "maxrss_MiB": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024,
            "elapsed_seconds": time.monotonic()-started,
            "raw_coefficient_rows_in_output": False,
            "runtime_scratch_file_dependencies": False, "new_parameter_sweep": False,
        }
    return result


def survival_high_slacks_certificate():
    """Bounded high-jet identities and the all-p m=7 lowering endpoint.

    The all-depth response theorem is an analytic input: its prefix proof
    is not inferred from checks at depths four/five. One univariate m=7
    certificate repairs the small-size gap in the formal lowering path.
    This proves no fourth survival LC gap or full stationary ULC statement.
    Earlier helpers, parameter grids, and scratch dependencies are absent.
    """
    import hashlib
    import math
    import resource
    import time
    from sympy.polys.fields import field

    def digest(rows):
        return hashlib.sha256(json.dumps(rows, separators=(",", ":")).encode()).hexdigest()

    def symbolic_depth_checks():
        jets = {(i, j): sp.Symbol("p%d_%d" % (i, j), positive=True)
                for i in range(1, 6) for j in range(i)}

        def failure_words(length):
            def extend(prefix, left):
                if len(prefix) == length:
                    if left == 0:
                        yield tuple(prefix)
                    return
                for count in range(left+1):
                    word = prefix+[count]
                    if len(word) < length and sum(word) < len(word):
                        continue
                    yield from extend(word, left-count)
            yield from extend([], length-1)

        def survival(depth):
            return 1-sp.Add(*(sp.Mul(*(jets[depth-j, count] for j, count in enumerate(word)))
                              for length in range(1, depth+1) for word in failure_words(length)))

        survivals = {r: survival(r) for r in range(1, 6)}
        counts = [len(list(failure_words(r))) for r in range(1, 6)]
        assert counts == [1, 1, 2, 5, 14]
        rho = sp.Symbol("rho", positive=True)
        high = []
        for depth in (4, 5):
            order = depth-2
            Q = sp.prod(jets[i, 0] for i in range(1, depth))
            gain = rho*Q/jets[depth-1, 0]
            x = jets[depth, 1]
            tail = Q*(1/jets[1, 0]+sum(jets[j, 1]/jets[j, 0] for j in range(2, depth)))
            loss = tail+2*x*Q/jets[depth, 0]+x*gain
            allowed = {(depth, order), (depth, order+1), (depth-1, order)}
            for (i, j), variable in jets.items():
                if i <= depth and j >= order and (i, j) not in allowed:
                    assert sp.diff(survivals[depth], variable) == 0
            response = -sp.diff(survivals[depth], jets[depth, order])
            response -= 2*x/jets[depth, 0]*sp.diff(survivals[depth], jets[depth, order+1])
            response -= rho*sp.diff(survivals[depth], jets[depth-1, order])
            assert sp.expand(response-loss) == 0
            assert sp.expand(-rho*sp.diff(survivals[depth-1], jets[depth-1, order])-gain) == 0
            assert sp.expand(-sp.diff(survivals[depth], jets[depth, depth-1])-Q) == 0
            high.append({"depth": depth, "highest_slack_response": True,
                         "next_highest_two_survival_responses": True,
                         "no_other_high_jet_support": True})

        F, G, J, K, L = sp.symbols("F G J K L", positive=True)
        dG, y, k1, x = sp.symbols("dG y k1 x", nonnegative=True)
        z0, k20, k30, ell30, ell40 = sp.symbols("z0 k20 k30 ell30 ell40", nonnegative=True)
        rho3, rho4 = sp.symbols("rho3 rho4", positive=True)
        sigma = sp.Symbol("sigma", nonnegative=True)
        u0 = x*x/L
        substitution = {
            jets[1, 0]: F, jets[2, 0]: G, jets[2, 1]: dG,
            jets[3, 0]: J, jets[3, 1]: y, jets[3, 2]: z0+rho3*rho4*sigma,
            jets[4, 0]: K, jets[4, 1]: k1, jets[4, 2]: k20+rho4*sigma,
            jets[4, 3]: k30+2*rho4*(k1/K+x/L)*sigma,
            jets[5, 0]: L, jets[5, 1]: x, jets[5, 2]: u0+sigma,
            jets[5, 3]: ell30+2*x*sigma/L,
            jets[5, 4]: ell40+3*x*x*sigma/L**2+sigma**2/L}
        B, C, D, E = [sp.expand(survivals[r].subs(substitution)) for r in (2, 3, 4, 5)]
        C0, D0, E0 = [value.subs(sigma, 0) for value in (C, D, E)]
        Q3, Q4 = F*G*J, F*G*J*K
        tail = J*(1-B)+y*G*F
        tail2 = K*(1-C)+k1*tail+(k20+rho4*sigma)*Q3
        tail3 = K*tail+k1*Q3
        queue_E = (1-L-x*(1-D)-(u0+sigma)*tail2
                   -(ell30+2*x*sigma/L)*tail3
                   -(ell40+3*x*x*sigma/L**2+sigma**2/L)*Q4)
        assert sp.expand(E-queue_E) == 0
        tail20 = K*(1-C0)+k1*tail+k20*Q3
        alpha = rho3*rho4*F*G
        beta = rho4*(tail+2*Q3*(k1/K+x/L))+k1*alpha
        gain2 = K*alpha+rho4*Q3
        gamma = x*beta+tail20+u0*gain2+2*x*tail3/L+3*x*x*Q4/L**2
        zeta = gain2+Q4/L
        assert sp.expand(C-(C0-alpha*sigma)) == 0
        assert sp.expand(D-(D0-beta*sigma)) == 0
        assert sp.expand(E-(E0-gamma*sigma-zeta*sigma**2)) == 0
        gap = D*D-C*E
        expected = (D0*D0-C0*E0+(C0*gamma+alpha*E0-2*D0*beta)*sigma
                    +(beta**2+C0*zeta-alpha*gamma)*sigma**2-alpha*zeta*sigma**3)
        assert sp.expand(gap-expected) == 0
        assert sp.expand(sp.diff(gap, sigma, 3)+6*alpha*zeta) == 0
        assert sp.expand(sp.diff(C*C-B*D, sigma, 2)-2*alpha*alpha) == 0
        low, upper, theta, h0, h1, scale = sp.symbols("low upper theta h0 h1 scale")
        derivative = lambda value: h0+h1*value-3*scale*value**2
        chord = (derivative((1-theta)*low+theta*upper)
                 -(1-theta)*derivative(low)-theta*derivative(upper))
        assert sp.expand(chord-3*scale*theta*(1-theta)*(upper-low)**2) == 0
        return {
            "scope": "symbolic identities only at queue depths four/five",
            "first_failure_prefix_counts": counts, "high_slack_checks": high,
            "depth_five": {
                "compact_fifth_queue_formula": True,
                "C_affine_response": True, "D_affine_response": True,
                "E_quadratic_response": True, "gap_cubic_identity": True,
                "strict_negative_cubic_coefficient": "-alpha*zeta; alpha=rho3*rho4*F*G>0, zeta=K*alpha+rho4*F*G*J+F*G*J*K/L>0",
                "derivative_concave_quadratic": True,
                "concave_quadratic_chord_identity": True,
                "preceding_gap_positive_quadratic_leading_term": "alpha^2",
                "physical_parameter_counterexample": False}}

    def m7_endpoint_certificate():
        ring, t = field("t", sp.QQ)
        p, q = t/(1+t), 1/(1+t)

        def predecessor(birth, later):
            down = (7-birth)*q
            denominator = [1+birth*p+down*(1-later[0]), -1-down*later[1],
                           -down*later[2], -down*later[3]]
            answer = [birth*p/denominator[0]]
            for j in range(1, 4):
                answer.append(-sum(denominator[k]*answer[j-k] for k in range(1, j+1))
                              /denominator[0])
            assert all(sum(denominator[k]*answer[j-k] for k in range(j+1))
                       == (birth*p if j == 0 else 0) for j in range(4))
            return answer

        jets = {8: [ring.zero]*4}
        for birth in range(7, 0, -1):
            jets[birth] = predecessor(birth, jets[birth+1])
        F = jets[1][0]
        G, dG = jets[2][:2]
        J, y, zeta = jets[3][:3]
        K, k1, k2, k3 = jets[4]
        L, x, l2, l3 = jets[5]
        M = jets[6][0]
        sigma2 = l2-x*x/L
        sigma3 = l3-x**3/L**2-2*x*sigma2/L
        assert sigma2 == 2*q*L*L*jets[6][2]/(5*p)
        assert sigma3 == 2*q*L*L*jets[6][3]/(5*p)
        hL = 2*q*L*L*M**4/(5*p*(6*p)**3)
        modified = {5: [L, x, l2, x**3/L**2+2*x*sigma2/L+hL]}
        for birth in range(4, 0, -1):
            modified[birth] = predecessor(birth, modified[birth+1])
            assert modified[birth][:3] == jets[birth][:3]
        B = 1-G-F*dG
        C = 1-J-y*G-y*dG*F-zeta*G*F
        Q3, Q4 = F*G*J, F*G*J*K
        tail = J*(1-B)+y*G*F
        Dactual = 1-K-k1*(1-C)-k2*tail-k3*Q3
        Dstar = 1-K-k1*(1-C)-k2*tail-modified[4][3]*Q3
        rho = 3*q*K*K/(4*p)
        Acal = rho*Q3
        Bcal = Q4*(1/F+dG/G+y/J+k1/K+2*x/L)+x*Acal
        assert Dstar-Dactual == Acal*(sigma3-hL)
        strong = C*C-B*Dstar
        endpoint = C*Bcal-2*Acal*Dstar
        assert endpoint == C*(Bcal-2*Acal*C/B)+2*Acal*strong/B

        def integer_pair(value):
            num = {k[0]: sp.Rational(c) for k, c in value.numer.items()}
            den = {k[0]: sp.Rational(c) for k, c in value.denom.items()}
            scale = sp.ilcm(*[c.q for c in list(num.values())+list(den.values())])
            ni = {k: int(c*scale) for k, c in num.items()}
            di = {k: int(c*scale) for k, c in den.items()}
            common = math.gcd(*list(ni.values()), *list(di.values()))
            ni, di = ({k: c//common for k, c in ni.items()},
                      {k: c//common for k, c in di.items()})
            if di[max(di)] < 0:
                ni, di = {k: -c for k, c in ni.items()}, {k: -c for k, c in di.items()}
            assert (ring.ring.from_dict({(k,): sp.QQ(c) for k, c in ni.items()})*value.denom
                    == ring.ring.from_dict({(k,): sp.QQ(c) for k, c in di.items()})*value.numer)
            return ni, di

        num, den = integer_pair(strong)
        nrows = [[k, c] for k, c in sorted(num.items())]
        drows = [[k, c] for k, c in sorted(den.items())]
        assert max(num) == max(den) == 44 and len(num) == 44 and len(den) == 45
        assert all(c > 0 for c in num.values()) and all(c > 0 for c in den.values())
        assert nrows[0] == [1, 3175200000] and den[0] == 23814000000
        assert digest(nrows) == "5932ae6dd97b373783d3ac4ed9b46b4774309ed631a372e6d44e064e42f90d62"
        assert digest(drows) == "af5d0d999cfee91327d0d9ed09b2fdf88d9962a4cfddb0a9cbfbfd3df424bf4d"
        # Independent prior actual m=7 base hashes, without calling its helper.
        actual_num, actual_den = integer_pair(C*C-B*Dactual)
        actual_num_hash = digest([[k, c] for k, c in sorted(actual_num.items())])
        actual_den_hash = digest([[k, c] for k, c in sorted(actual_den.items())])
        assert actual_num_hash == "c90cb02c9bce869705eeabba61b03fa37c16e0f51c106aeb30deb4fcaf0b8e6a"
        assert actual_den_hash == "739a74ec784d453768dc00d5e72e68ae462bfcc9f3d0c3529ebe5e1c21be04df"
        return {
            "scope": "m=7, p=t/(1+t), t>0, boundary rate one",
            "target": "C_actual^2-B_actual*Dstar>0",
            "modification": "L=P5; retain actual zeroth/first/second jets and lower reciprocal sigma3 to h_L",
            "h_L": "2*q*L^2*M^4/[(5p)(6p)^3], M=P6(1)",
            "numerator_degree": 44, "denominator_degree": 44,
            "numerator_ascending_coefficients": [num.get(k, 0) for k in range(45)],
            "denominator_ascending_coefficients": [den.get(k, 0) for k in range(45)],
            "positive_numerator_terms": 44, "positive_denominator_terms": 45,
            "negative_terms": 0, "strict_lowest_numerator_term": nrows[0],
            "positive_denominator_constant": den[0],
            "numerator_sha256": digest(nrows), "denominator_sha256": digest(drows),
            "actual_m7_base_hashes": {"numerator": actual_num_hash, "denominator": actual_den_hash},
            "exact_identity_checks": {
                "actual_and_modified_Taylor_convolutions": True,
                "all_zeroth_first_second_jets_preserved": True,
                "actual_reciprocal_slack_formulas": True,
                "Dstar_minus_Dactual_equals_Acal_times_removed_slack": True,
                "endpoint_response_decomposition": True,
                "integer_clearing_preserves_rational_function": True,
                "actual_m7_base_hashes_match": True},
            "endpoint_response_consequence": "C*Bcal-2*Acal*Dstar=C*(Bcal-2*Acal*C/B)+(2*Acal/B)*(C^2-B*Dstar)>3*Acal*C>0; use actual 0<C<=B and Bcal>5*Acal",
            "lowering_path_consequence": "D(sigma3)=Dstar-Acal*(sigma3-h_L), so the endpoint controls the whole actual-to-floor path",
            "m7_formal_scope_gap_closed": True}

    started = time.monotonic()
    previous_limits = resource.getrlimit(resource.RLIMIT_AS)
    limit = min([1900*1024*1024]+[value for value in previous_limits if value >= 0])
    resource.setrlimit(resource.RLIMIT_AS, (limit, previous_limits[1]))
    try:
        finite_checks = symbolic_depth_checks()
        endpoint_check = m7_endpoint_certificate()
        result = {
            "status": "PASS: depth-four/five identities and exact all-p m=7 endpoint",
            "model": "ordinary m-bit Ehrenfest passages, 0<p<1, independent boundary rate one",
            "analytic_all_depth_input": {
                "status": "proved analytically by the return-moment identity and prefix support argument; not extrapolated from finite checks",
                "range": "4<=r<=m, actual zeroth transforms, nonnegative formal first jets",
                "normalization": "S_i=a_(m-i); P_r(1-h)=f_r/[1-sum_(j>=1)sigma_j*h^j/f_r]",
                "return_moment_identity": "i*p*I_(i-1)-(m+1-i)*q*I_i=1, I_j=integral_0^1 theta^j*(p+q*theta)^(-m-2)dtheta",
                "zeroth_value_bound": "m*q*f1<p, hence delta_i*f1<p",
                "response_definitions": "Q=prod_(i=1)^(r-1)f_i; rho=delta_(r-1)*f_(r-1)^2/((r-1)*p); Acal=rho*Q/f_(r-1)",
                "response_Bcal": "Q*[1/f1+sum_(j=2)^(r-1)d_j/f_j+2*d_r/f_r]+d_r*Acal",
                "strict_response_gain": "Bcal/Acal>r-1+1/p>r",
                "actual_local_monotonicity": "Both highest reciprocal-slack partial derivatives of S_(r-1)^2-S_(r-2)*S_r are strictly positive at actual jets",
                "formal_path_hypotheses": "Hold lower reciprocal slacks actual; T=S_(r-3),U=S_(r-2) remain actual, 0<U<=T. Require U^2>=T*V along the lowering path, or the weaker V<=r*U/2, V=S_(r-1)",
                "actual_preceding_LC_alone_suffices_for_lowering": False},
            "finite_symbolic_checks": finite_checks,
            "m7_endpoint_certificate": endpoint_check,
            "fourth_gap_lowering_reduction": {
                "scope": "m>=5, actual x and sigma2; lower sigma3 to h_L and sigma4 to j_L",
                "comparison": "actual fourth gap>=H5(x_actual,sigma2_actual,h_L,j_L)",
                "dimension_coverage": {"m>=8": "prior formal sixth-floor third-gap certificate",
                                       "m=5,6": "floors already actual; equality",
                                       "m=7": "new exact endpoint certificate in this helper"},
                "h_L": "d*L^2*M^4/[(5p)(6p)^3], d=(m-5)*q",
                "j_L": "d*L^2*M^5/[(5p)(6p)^4]",
                "m5_convention": "omit M and set both floors to zero",
                "fourth_gap_LC_proved": False},
            "remaining_gap": {
                "second_slack_polynomial": "D(s)^2-C(s)*E(s) is cubic with coefficient -alpha*zeta<0; its derivative is strictly concave quadratic",
                "bounded_monotonicity_criterion": "On a nonempty [lo,hi], derivative nonnegative throughout iff it is nonnegative at both endpoints",
                "conditions_not_proved": [
                    "Sign both derivative endpoints on an interval covering the actual second-slack lowering path",
                    "Then sign the separate lower fourth-gap value H5(x,e_L,h_L,j_L)"],
                "interval_condition": "If using [e_L,smax(x)], require e_L<=smax(x)=x*(1-x/[L*(1-L)])",
                "response_route_conditions": "For the alternative response comparison retain B>0,C>=0,E>=0 and C^2>=B*D throughout; actual order C<=B is not automatic after formal lowering",
                "unbounded_obstruction_scope": "Derivative is eventually negative while preceding gap has positive quadratic leading term; C,D eventually negative, so this is not an actual positive-profile counterexample"},
            "full_survival_LC_proved": False, "stationary_ULC_or_NA_proved": False,
            "other_certificate_helpers_called": False, "runtime_scratch_dependencies": False,
            "new_size_or_rate_sweep": False, "large_coefficient_production": False,
            "address_space_soft_cap_MiB": limit//(1024*1024),
            "maxrss_MiB": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024,
            "elapsed_seconds": time.monotonic()-started}
    finally:
        resource.setrlimit(resource.RLIMIT_AS, previous_limits)
    return result


def survival_fourth_reduction_certificate():
    """Bounded fourth-gap reduction, c=3 interval, and four actual bases.

    Response growth and the m>=8 actual-first-jet reduction are analytic
    inputs, with their rational identities checked here. The fixed c=3
    interval and m=5,...,8 bases have exact coefficient certificates.
    The all-size, all-rate fourth survival gap remains open. No previous
    certificate helper is called; boundary_profiles is used only at the
    four prescribed rational normalization points.
    """
    import hashlib
    import math
    import resource
    import time
    from sympy.polys.fields import field
    from sympy.polys.rings import ring

    def digest(rows):
        return hashlib.sha256(json.dumps(rows, separators=(",", ":")).encode()).hexdigest()

    def response_growth_identity():
        F, G, J, K, L, p, theta, eta, xi = sp.symbols(
            "F G J K L p theta eta xi", positive=True)
        x, s = sp.symbols("x s", nonnegative=True)
        Q = F*G*J
        rho4, rho3 = theta*K, eta*J/K
        w0 = K/(4*p)
        w = w0+theta*x
        v0 = J/(3*p)+eta*w0
        v = v0+eta*theta*x
        u = G/(2*p)+xi*v
        U0 = 1/F+G/(2*p)+(1+xi)*v0
        h = eta*(1+xi)
        tau = 2+eta+h
        B0 = U0+(2+eta)*w0
        B1 = theta*tau+2/L
        k1, y, g1 = K*w, J*v, G*u
        k20 = k1**2/K+rho4*x**2/L
        j20 = y**2/J+rho3*k20
        B = 1-G-g1*F
        C0 = 1-J-y*(1-B)-j20*G*F
        Tstar = J*(1-B)+y*G*F
        T20 = K*(1-C0)+k1*Tstar+k20*Q
        T3 = K*Tstar+k1*Q
        alpha = rho3*rho4*F*G
        beta = rho4*(Tstar+2*Q*(k1/K+x/L))+k1*alpha
        Rstar = K*alpha+rho4*Q
        gamma = x*beta+T20+x**2/L*Rstar+2*x*T3/L+3*x**2*K*Q/L**2
        zeta = Rstar+K*Q/L
        a0 = (1/(F*G)+(v0+w0-2*theta)*U0+(1+eta)*w0**2
              -2*theta*(2+eta)*w0)/theta
        a1 = ((2+eta)*U0+(4+3*eta)*w0+h*(v0+w0)-2*theta*tau
              +2/(theta*L)*(U0+w0-2*theta))
        a2 = B1+theta*(1+eta)*(1+h)+2*(h+2+eta)/L+3/(theta*L**2)
        expected = a0+a1*x+a2*x**2+2*(1+eta+1/(theta*L))*s
        normalized = (gamma+2*zeta*s-2*beta)/(rho4*Q)
        assert sp.cancel(normalized-expected) == 0
        assert sp.cancel(beta/(rho4*Q)-(B0+B1*x)) == 0
        polynomial = sp.Poly(sp.cancel(normalized), x, s)
        assert set(polynomial.monoms()) == {(0, 0), (1, 0), (2, 0), (0, 1)}
        A0 = 1/F+G/(2*p)
        P0 = ((1-2*theta*G)/(F*G)+w0/F+G/(2*p)*w0+(1+eta)*w0**2
              -2*theta*G/(2*p)-2*theta*(2+eta)*w0)
        rearranged = (1+xi)*v0**2+(A0+(1+xi)*(w0-2*theta))*v0+P0
        assert sp.cancel(theta*a0-rearranged) == 0
        AG = 1/G+1-G
        positive_remainder = (4*(1-K)/K*AG
                              +2*(2-3*G)*(3+2*G)/(3*G)+2*(1-eta))
        assert sp.cancel(4/K*AG-4-2*eta-sp.Rational(4, 3)-positive_remainder) == 0
        return {
            "status": "PASS: exact response quadratic and positive-remainder identities",
            "definitions": {
                "zeroths": "F,G,J,K,L=P_1(1),...,P_5(1), actual Ehrenfest passages",
                "theta": "(m-4)*q*K/(4*p)", "eta": "(m-3)*q*J*K/(3*p)",
                "xi": "(m-2)*q*G*J/(2*p)", "x": "formal -P_5'(1)",
                "s": "second reciprocal slack of P_5", "w0": "K/(4*p)",
            },
            "normalized_identity": "(gamma+2*zeta*s-2*beta)/(rho4*F*G*J)=a0+a1*x+a2*x^2+2*(1+eta+1/(theta*L))*s",
            "coefficients": {"a0": str(a0), "a1": str(a1), "a2": str(a2)},
            "nonzero_monomials_x_s": [list(v) for v in polynomial.monoms()],
            "response_identity": True, "beta_affine_identity": True,
            "constant_coefficient_rearrangement": True,
            "constant_bound_remainder_decomposition": True,
            "analytic_sign_inputs": [
                "i*p*I_(i-1)-(m+1-i)*q*I_i=1 for actual return integrals",
                "theta*F<K/4, theta*G<K/2, eta<K<1, xi<J<1",
                "0<G<2/3 and w0=1-K-theta*(1-L)<1-K",
            ],
            "analytic_conclusion": "a0>4*w0/3, a1>2*theta, a2>0; response growth for all formal x,s>=0",
            "parameter_evaluations": 0, "root_calculations": 0,
        }

    def failure_words(length):
        def extend(prefix, left):
            if len(prefix) == length:
                if left == 0:
                    yield tuple(prefix)
                return
            for count in range(left+1):
                word = prefix+[count]
                if len(word) < length and sum(word) < len(word):
                    continue
                yield from extend(word, left-count)
        yield from extend([], length-1)

    prefixes = {k: list(failure_words(k)) for k in range(1, 6)}
    assert [len(prefixes[k]) for k in range(1, 6)] == [1, 1, 2, 5, 14]

    def queue_profiles(passage, one):
        return [one-sum((math.prod(passage[depth-j][count]
                                   for j, count in enumerate(word))
                         for length in range(1, depth+1)
                         for word in prefixes[length]), one*0)
                for depth in range(1, 6)]

    def c3_interval_certificate():
        Q, E = field("E", sp.QQ)
        R, X = ring("X", Q.to_domain())
        Z, z = ring("z", sp.QQ)
        order = 4
        integrals = [(E-1)/3]
        for i in range(1, 8):
            integrals.append((i*integrals[-1]-1)/3)
        f = {i: integrals[i]/integrals[i-1] for i in range(1, 8)}
        sums = [sum(sp.Rational(3)**j/sp.factorial(j) for j in range(i+1))
                for i in range(8)]
        for i in range(1, 8):
            assert f[i] == Q(sp.Rational(i, 3))*(E-sums[i])/(E-sums[i-1])
        for i in range(1, 7):
            assert f[i] == i/(1+i+3*(1-f[i+1]))
        L, M, N = f[5], f[6], f[7]

        def multiply(a, b, one):
            out = [one*0 for _ in range(order+1)]
            for i, ai in enumerate(a):
                for j, bj in enumerate(b):
                    if i+j <= order:
                        out[i+j] += ai*bj
            return out

        def geometric(f0, tail, one):
            assert tail[0] == 0
            power = [one]+[one*0]*order
            total = power[:]
            for _ in range(order):
                power = multiply(power, tail, one)
                total = [a+b for a, b in zip(total, power)]
            return [f0*c for c in total]

        def upstream(i, later, one):
            tail = [one*0]+[(one if k == 1 else one*0)*f[i]/i
                            +3*later[k]*f[i]/i for k in range(1, order+1)]
            row = geometric(one*f[i], tail, one)
            assert all(row[k] == one*f[i]/i*(row[k-1]+3*sum(
                (later[j]*row[k-j] for j in range(1, k+1)), one*0))
                for k in range(1, order+1))
            return row

        def reciprocal_slacks(series):
            reciprocal = [1/series[0]]
            for k in range(1, order+1):
                reciprocal.append(-sum(series[j]*reciprocal[k-j]
                                       for j in range(1, k+1))/series[0])
            return [-series[0]**2*reciprocal[k] for k in range(1, order+1)]

        sixth_base = [M**(k+1)/6**k for k in range(order+1)]
        x6, e6, h6, j6 = reciprocal_slacks(upstream(5, sixth_base, Q.one))
        assert x6 == L**2/5*(1+M**2/2)
        for k, slack in enumerate((e6, h6, j6), 2):
            assert slack == 3*L**2*M**(k+1)/(5*6**k)
        seventh_base = [N**(k+1)/7**k for k in range(order+1)]
        M_from_seventh = upstream(6, seventh_base, Q.one)
        L_from_seventh = upstream(5, M_from_seventh, Q.one)
        x7, e7, h7, j7 = reciprocal_slacks(L_from_seventh)
        assert M_from_seventh[0] == M and L_from_seventh[0] == L
        for k, slack in enumerate((e7, h7, j7), 2):
            assert slack == 3*L**2*M_from_seventh[k]/5

        def profiles_and_gap(slacks):
            tail = [R.zero, X*R(1/L)]+[R(s/L) for s in slacks]
            passage = {5: geometric(R(L), tail, R.one)}
            for i in range(4, 0, -1):
                passage[i] = upstream(i, passage[i+1], R.one)
            profiles = queue_profiles(passage, R.one)
            # Independent absorption recursion, without first-failure words.
            extinction = [R.one]
            for i in range(1, 6):
                nxt = [R.one]
                for queue in range(1, i+1):
                    nxt.append(sum((passage[i][n]*extinction[queue-1+n]
                                    for n in range(i-queue+1)), R.zero))
                extinction = nxt
                assert profiles[i-1] == R.one-extinction[1]
            return profiles, profiles[3]**2-profiles[2]*profiles[4]

        S6, H6 = profiles_and_gap((e6, h6, j6))
        S7, H7 = profiles_and_gap((e7, h7, j7))
        assert H7.degree() == 6
        upper = L*(1-L)
        low = sum(sp.Rational(3)**k/sp.factorial(k) for k in range(19))
        high = low+sp.Rational(3)**19/sp.factorial(19)/(1-sp.Rational(3, 20))
        assert (low, high) == (sp.Rational(19599925833041, 975822848000),
                               sp.Rational(6330776047615183, 315190779904000))

        def homogeneous_from_taylor(coefficients, lower, higher, degree):
            taylor = [(higher-lower)**k*sum(
                (coefficient*math.comb(d, k)*lower**(d-k)
                 for d, coefficient in coefficients.items() if d >= k), lower*0)
                for k in range(degree+1)]
            return [math.comb(degree, j)*sum(
                (taylor[k]*sp.Rational(math.comb(j, k), math.comb(degree, k))
                 for k in range(j+1)), lower*0) for j in range(degree+1)]

        polynomial_records = []

        def polynomial_packet(poly, expected_sign):
            degree = poly.degree()
            transformed = Z.zero
            for (k,), coefficient in poly.items():
                transformed += Z(coefficient)*(Z(low)+Z(high)*z)**k*(1+z)**(degree-k)
            independent = homogeneous_from_taylor(
                {k[0]: sp.Rational(c) for k, c in poly.items()}, low, high, degree)
            assert all(transformed[(k,)] == independent[k] for k in range(degree+1))
            assert len(transformed.items()) == degree+1
            assert all(expected_sign*c > 0 for c in transformed.coeffs())
            rows = [[k[0], str(c)] for k, c in sorted(poly.items())]
            transformed_rows = [[k[0], str(c)] for k, c in sorted(transformed.items())]
            record = {
                "degree": degree, "coefficients": rows,
                "interval_coefficient_signs": [expected_sign],
                "strict_transformed_coefficient_count": degree+1,
                "transformed_sha256": digest(transformed_rows),
                "affine_Taylor_conversion_matches": True,
            }
            polynomial_records.append(record)
            return record

        def packet(expression, numerator_sign, denominator_sign):
            return {"numerator": polynomial_packet(expression.numer, numerator_sign),
                    "denominator": polynomial_packet(expression.denom, denominator_sign)}

        sixth = {
            "scope": "c=3 limit, actual zeroths, sixth floors, x=x6; a formal obstruction only",
            "fourth_gap": packet(H6.evaluate(X, x6), -1, 1),
            "moment_envelope_margin": packet(x6*(1-x6/upper)-e6, -1, -1),
            "survivals": {str(i): packet(profile.evaluate(X, x6), 1, 1)
                          for i, profile in enumerate(S6, 1)},
            "actual_fourth_gap_counterexample": False,
        }
        seventh = {
            "x7_minus_x6": packet(x7-x6, 1, 1),
            "x7_old_higher_floors": packet(H6.evaluate(X, x7), -1, 1),
            "x7_all_seventh_floors": packet(H7.evaluate(X, x7), 1, 1),
            "x7_moment_envelope_new_e": packet(x7*(1-x7/upper)-e7, -1, -1),
        }
        degree = H7.degree()
        independent_bernstein = homogeneous_from_taylor(
            {k[0]: c for k, c in H7.items()}, x7, upper, degree)
        coefficients = []
        for j in range(degree+1):
            value = Q.zero
            for (k,), coefficient in H7.items():
                for ell in range(max(0, j-degree+k), min(j, k)+1):
                    value += (coefficient*math.comb(k, ell)*x7**(k-ell)*upper**ell
                              *math.comb(degree-k, j-ell))
            assert value == independent_bernstein[j]
            coefficients.append({"u_degree": j, **packet(value, 1, 1)})
        assert independent_bernstein[0] == H7.evaluate(X, x7)
        assert independent_bernstein[-1] == H7.evaluate(X, upper)
        interval = {"degree": degree, "bernstein_coefficients": coefficients,
                    "width": packet(upper-x7, 1, 1),
                    "lower_and_upper_endpoint_identities": True,
                    "independent_affine_Taylor_conversion_matches": True,
                    "all_seven_coefficients_strictly_positive": True}
        sign_count = sum(row["strict_transformed_coefficient_count"] for row in polynomial_records)
        assert len(polynomial_records) == 38 and sign_count == 320
        return {
            "status": "PASS: 38 polynomial records, 320 strict transformed coefficient signs",
            "parameter_scope": "only the c=3 limit p=1-3/m; E=exp(3)",
            "exp3_interval": [str(low), str(high)],
            "enclosure_proof": "sum through degree 18 < exp(3) < that sum + term19/(1-3/20)",
            "zeroth_construction": "I0=(E-1)/3; Ij=(j*I_(j-1)-1)/3; f_i=I_i/I_(i-1)",
            "polynomial_sign_conversion": "(1+z)^d*P((low+high*z)/(1+z)), z>=0; upper endpoint from leading coefficient",
            "hash_convention": "SHA256 compact JSON ascending nonzero [degree,rational-string] rows of the transformed polynomial",
            "first_jet_interval": "x7 <= x <= L*(1-L)",
            "homogeneous_Bernstein_convention": "(1+u)^6*H7((x7+L*(1-L)*u)/(1+u))",
            "geometric_series_equals_first_step_coefficients": True,
            "reciprocal_slack_identities": True,
            "queue_words_equal_independent_absorption_recursion": True,
            "sixth": sixth, "seventh": seventh, "interval": interval,
            "polynomial_certificate_count": len(polynomial_records),
            "coefficient_sign_count": sign_count,
            "actual_fourth_gap_counterexample": False,
        }

    def actual_small_bases():
        Q, t = field("t", sp.QQ)
        p, q = t/(1+t), 1/(1+t)
        tt = sp.Symbol("t")

        def primitive_pair(value):
            originals = [sp.Poly(poly.as_expr(), tt, domain=sp.QQ)
                         for poly in (value.numer, value.denom)]
            primitive, multipliers = [], []
            for poly in originals:
                scale, integer = poly.clear_denoms(convert=True)
                content, part = integer.primitive()
                multiplier = sp.Rational(content, scale)
                if part.LC() < 0:
                    part, multiplier = -part, -multiplier
                assert (poly-part*multiplier).is_zero
                primitive.append(part)
                multipliers.append(multiplier)
            numerator, denominator = primitive
            prefactor = multipliers[0]/multipliers[1]
            assert numerator.gcd(denominator).degree() == 0
            assert numerator.content() == denominator.content() == 1
            assert (originals[0]*denominator-prefactor*numerator*originals[1]).is_zero
            return prefactor, numerator, denominator

        def polynomial_packet(poly):
            rows = [[int(k[0]), int(c)] for k, c in sorted(poly.terms())]
            assert all(c > 0 for _, c in rows)
            return {
                "degree": int(poly.degree()), "nonzero_terms": len(rows),
                "positive_terms": len(rows), "negative_terms": 0,
                "ascending_coefficients": [int(poly.nth(k)) for k in range(poly.degree()+1)],
                "lowest_nonzero_term": rows[0], "coefficient_sha256": digest(rows),
                "primitive_integer_content": 1,
            }

        cases = {}
        for m in (5, 6, 7, 8):
            jets = {m+1: [Q.zero]*5}
            for birth in range(m, 0, -1):
                down, later = (m-birth)*q, jets[birth+1]
                denominator = [1+birth*p+down*(1-later[0]), -1-down*later[1]]
                denominator.extend(-down*later[j] for j in range(2, 5))
                row = [birth*p/denominator[0]]
                for j in range(1, 5):
                    row.append(-sum(denominator[k]*row[j-k]
                                    for k in range(1, j+1))/denominator[0])
                assert all(sum(denominator[k]*row[j-k] for k in range(j+1))
                           == (birth*p if j == 0 else 0) for j in range(5))
                jets[birth] = row
            survival = queue_profiles(jets, Q.one)
            F = jets[1][0]
            G, g1 = jets[2][:2]
            J, j1, j2 = jets[3][:3]
            K, k1, k2, k3 = jets[4][:4]
            L, l1, l2, l3, l4 = jets[5]
            A, B = 1-F, 1-G-F*g1
            C = 1-J-j1*G-j1*g1*F-j2*G*F
            Q3, Q4 = J*G*F, K*J*G*F
            tail = J*(1-B)+j1*G*F
            D = 1-K-k1*(1-C)-k2*tail-k3*Q3
            E = (1-L-l1*(1-D)-l2*(K*(1-C)+k1*tail+k2*Q3)
                 -l3*(K*tail+k1*Q3)-l4*Q4)
            assert (A, B, C, D, E) == tuple(survival)
            gap = D*D-C*E
            prefactor, numerator, denominator = primitive_pair(gap)
            ndata, ddata = polynomial_packet(numerator), polynomial_packet(denominator)
            assert prefactor > 0 and numerator.nth(1) > 0 and denominator.nth(0) == 1
            assert ndata["degree"] == ddata["degree"] == 9*m-16
            assert ndata["positive_terms"] == ndata["degree"]
            assert ddata["positive_terms"] == ddata["degree"]+1

            def at_probe(value):
                return sp.cancel(value.as_expr().subs(tt, sp.Rational(2, 3)))

            probe_values = [at_probe(value) for value in survival]
            diagonal, _ = boundary_profiles(m, sp.Rational(2, 5))
            assert probe_values == [diagonal[m-r] for r in range(1, 6)]
            gap_probe = at_probe(gap)
            assert gap_probe == diagonal[m-4]**2-diagonal[m-3]*diagonal[m-5] > 0
            cases[str(m)] = {
                "status": "strict all-p positive coefficient certificate",
                "target": "S4^2-S3*S5 = a_(m-4)^2-a_(m-3)*a_(m-5)",
                "normalization": "gap=positive_prefactor*primitive_numerator(t)/primitive_denominator(t)",
                "positive_prefactor": str(prefactor), "numerator": ndata, "denominator": ddata,
                "rational_pair_coprime": True, "primitive_integer_polynomials": True,
                "exact_clearing_identity": True, "strict_denominator_constant": 1,
                "normalized_jets_zero_through_four_convolution_checks": True,
                "queue_enumeration_equals_compact_five_survival_formula": True,
                "probe_p": "2/5", "probe_t": "2/3",
                "S1_through_S5_at_probe": list(map(str, probe_values)),
                "all_five_match_independent_boundary_profiles": True,
                "rational_gap_at_probe": str(gap_probe), "probe_gap_positive": True,
                "actual_negative_gap_found": False,
            }
        return {
            "status": "PASS: actual fourth-gap bases m=5,6,7,8, every 0<p<1",
            "parameterization": "p=t/(1+t), q=1/(1+t), t>0",
            "actual_passage_recursion": "P_r(z)=r*p/[z+r*p+(m-r)*q*(1-P_(r+1)(z))]; bottom P_m(z)=m*p/(z+m*p)",
            "jets": "P_r(1-h) coefficients through degree four, normalized by factorials",
            "first_failure_prefix_counts": [1, 1, 2, 5, 14],
            "coefficient_hash_convention": "SHA256 compact JSON ascending nonzero [degree,integer] rows of each primitive polynomial",
            "cases": cases, "all_four_bases_certified": True,
            "positive_numerator_coefficient_count": sum(row["numerator"]["positive_terms"] for row in cases.values()),
            "positive_denominator_coefficient_count": sum(row["denominator"]["positive_terms"] for row in cases.values()),
            "independent_normalization_cases": [{"m": m, "p": "2/5"} for m in (5, 6, 7, 8)],
        }

    started = time.monotonic()
    previous_limits = resource.getrlimit(resource.RLIMIT_AS)
    cap = min([1900*1024**2]+[v for v in previous_limits if v != resource.RLIM_INFINITY])
    resource.setrlimit(resource.RLIMIT_AS, (cap, previous_limits[1]))
    try:
        response = response_growth_identity()
        c3 = c3_interval_certificate()
        bases = actual_small_bases()
        result = {
            "status": "PASS: fourth-gap response reduction, c=3 interval, and four actual bases",
            "model": "actual infinite-horizon Ehrenfest survival; births p*(m-y), deaths (1-p)*y; independent boundary rate one",
            "target": "S4^2-S3*S5 = a_(m-4)^2-a_(m-3)*a_(m-5), S_r=a_(m-r)",
            "response_growth": response,
            "analytic_one_jet_reduction": {
                "status": "proved analytic input; exact response identities checked above",
                "scope": "integer m>=8, 0<p<1, actual zeroth transforms AND actual first L jet x",
                "floor_hypotheses": "e6<=e*<=sigma2_actual, h6<=h*<=sigma3_actual, j6<=j*<=sigma4_actual",
                "sixth_floor_definitions": "M=P6(1); rho5=(m-5)*q*L^2/(5*p); e6=rho5*M^3/(6*p)^2, h6=rho5*M^4/(6*p)^3, j6=rho5*M^5/(6*p)^4",
                "formal_profile_input": "0<C(s)<B, D(s)>0, E(s)>0, C(s)^2>B*D(s) throughout the actual lowering interval; uses strong preceding-gap certificates",
                "second_slack_derivative": "d_s(D^2-C*E)=alpha*E+2*beta*(C-D)+C*(gamma+2*zeta*s-2*beta)>0",
                "conclusion": "actual fourth gap >= H5(x_actual,e*,h*,j*)",
                "seventh_floors": "N=P7(1); Mhat(h)=6*p/[1-h+6*p+(m-6)*q*(1-N/(1-N*h/(7*p)))]; mu_j=[h^j]Mhat; x7=L^2/(5*p)+rho5*mu1; (e7,h7,j7)=rho5*(mu2,mu3,mu4)",
                "seventh_floor_bounds": "sixth floors <= seventh floors <= actual slacks; x7<=x_actual<=L*(1-L)",
                "arbitrary_formal_first_jet_lowering_proved": False,
                "all_size_fourth_gap_sign_follows_from_reduction_alone": False,
            },
            "c3_limit": c3,
            "analytic_c3_consequence": {
                "status": "proved by the analytic reduction and continuity of seven strict limiting Bernstein coefficients",
                "quantifiers": "there exists m0>=8 such that every integer m>=m0 with p=1-3/m has a_(m-4)^2>a_(m-3)*a_(m-5)",
                "zeroth_limit": "f_i -> I_i/I_(i-1), I_i=integral_0^1 theta^i*exp(3*(1-theta)) dtheta, by dominated convergence of return integrals",
                "continuity_inputs": "degree at most six; positive nonzero limiting denominators; all seven limiting coefficients and interval width strictly positive",
                "actual_first_jet_convergence_required": False,
                "uniform_positivity_over_moving_first_jet_interval": True,
                "explicit_threshold": False,
            },
            "actual_small_bases": bases,
            "remaining_gap": "all-size, all-rate fourth survival LC remains OPEN; no full survival LC or stationary ULC/NA conclusion",
            "all_size_fourth_gap_claim": False, "full_stationary_ULC_or_NA_claim": False,
            "old_high_slacks_metadata_modified": False,
            "large_generic_expansion_or_parameter_grid": False, "root_calculations": 0,
            "old_certificate_helpers_called": [], "boundary_profile_normalization_calls": 4,
            "effective_soft_limit_bytes": cap,
            "elapsed_seconds": round(time.monotonic()-started, 3),
            "peak_rss_MiB": round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024, 2),
        }
    finally:
        resource.setrlimit(resource.RLIMIT_AS, previous_limits)
    result["resource_limits_restored"] = resource.getrlimit(resource.RLIMIT_AS) == previous_limits
    return result


def survival_fourth_lc_certificate():
    """All-seven formal fourth-gap certificate and actual m=5,...,9 bases.

    The frozen producer's exact rational construction is retained here.
    Generic m>=10 uses seventh derivative floors, the depth-ten P7 value
    interval, and seven homogeneous Bernstein coefficients of a sextic.
    Independent reconstruction and coverage audits close the all-size
    conclusion; report (55q) records the analytical comparison and scope.
    Only this selected helper imports optional python-flint. Earlier
    certificate helpers are not called or modified.
    """
    import hashlib
    import math
    import resource
    import time
    from math import comb
    from sympy.polys.fields import field

    def digest(rows):
        return hashlib.sha256(json.dumps(rows, separators=(",", ":")).encode()).hexdigest()

    def failure_words(length):
        def extend(prefix, left):
            if len(prefix) == length:
                if left == 0:
                    yield tuple(prefix)
                return
            for count in range(left+1):
                word = prefix+[count]
                if len(word) < length and sum(word) < len(word):
                    continue
                yield from extend(word, left-count)
        yield from extend([], length-1)

    prefixes = {k: list(failure_words(k)) for k in range(1, 6)}

    assert [len(prefixes[k]) for k in range(1, 6)] == [1, 1, 2, 5, 14]

    def queue_profiles(passage, one):
        return [one-sum((math.prod(passage[depth-j][count]
                                   for j, count in enumerate(word))
                         for length in range(1, depth+1)
                         for word in prefixes[length]), one*0)
                for depth in range(1, 6)]

    def actual_small_bases():
        Q, t = field("t", sp.QQ)
        p, q = t/(1+t), 1/(1+t)
        tt = sp.Symbol("t")

        def primitive_pair(value):
            originals = [sp.Poly(poly.as_expr(), tt, domain=sp.QQ)
                         for poly in (value.numer, value.denom)]
            primitive, multipliers = [], []
            for poly in originals:
                scale, integer = poly.clear_denoms(convert=True)
                content, part = integer.primitive()
                multiplier = sp.Rational(content, scale)
                if part.LC() < 0:
                    part, multiplier = -part, -multiplier
                assert (poly-part*multiplier).is_zero
                primitive.append(part)
                multipliers.append(multiplier)
            numerator, denominator = primitive
            prefactor = multipliers[0]/multipliers[1]
            assert numerator.gcd(denominator).degree() == 0
            assert numerator.content() == denominator.content() == 1
            assert (originals[0]*denominator-prefactor*numerator*originals[1]).is_zero
            return prefactor, numerator, denominator

        def polynomial_packet(poly):
            rows = [[int(k[0]), int(c)] for k, c in sorted(poly.terms())]
            assert all(c > 0 for _, c in rows)
            return {
                "degree": int(poly.degree()), "nonzero_terms": len(rows),
                "positive_terms": len(rows), "negative_terms": 0,
                "ascending_coefficients": [int(poly.nth(k)) for k in range(poly.degree()+1)],
                "lowest_nonzero_term": rows[0], "coefficient_sha256": digest(rows),
                "primitive_integer_content": 1,
            }

        cases = {}
        for m in (5, 6, 7, 8, 9):
            jets = {m+1: [Q.zero]*5}
            for birth in range(m, 0, -1):
                down, later = (m-birth)*q, jets[birth+1]
                denominator = [1+birth*p+down*(1-later[0]), -1-down*later[1]]
                denominator.extend(-down*later[j] for j in range(2, 5))
                row = [birth*p/denominator[0]]
                for j in range(1, 5):
                    row.append(-sum(denominator[k]*row[j-k]
                                    for k in range(1, j+1))/denominator[0])
                assert all(sum(denominator[k]*row[j-k] for k in range(j+1))
                           == (birth*p if j == 0 else 0) for j in range(5))
                jets[birth] = row
            survival = queue_profiles(jets, Q.one)
            F = jets[1][0]
            G, g1 = jets[2][:2]
            J, j1, j2 = jets[3][:3]
            K, k1, k2, k3 = jets[4][:4]
            L, l1, l2, l3, l4 = jets[5]
            A, B = 1-F, 1-G-F*g1
            C = 1-J-j1*G-j1*g1*F-j2*G*F
            Q3, Q4 = J*G*F, K*J*G*F
            tail = J*(1-B)+j1*G*F
            D = 1-K-k1*(1-C)-k2*tail-k3*Q3
            E = (1-L-l1*(1-D)-l2*(K*(1-C)+k1*tail+k2*Q3)
                 -l3*(K*tail+k1*Q3)-l4*Q4)
            assert (A, B, C, D, E) == tuple(survival)
            gap = D*D-C*E
            prefactor, numerator, denominator = primitive_pair(gap)
            ndata, ddata = polynomial_packet(numerator), polynomial_packet(denominator)
            assert prefactor > 0 and numerator.nth(1) > 0 and denominator.nth(0) == 1
            assert ndata["degree"] == ddata["degree"] == 9*m-16
            assert ndata["positive_terms"] == ndata["degree"]
            assert ddata["positive_terms"] == ddata["degree"]+1

            def at_probe(value):
                return sp.cancel(value.as_expr().subs(tt, sp.Rational(2, 3)))

            probe_values = [at_probe(value) for value in survival]
            diagonal, _ = boundary_profiles(m, sp.Rational(2, 5))
            assert probe_values == [diagonal[m-r] for r in range(1, 6)]
            gap_probe = at_probe(gap)
            assert gap_probe == diagonal[m-4]**2-diagonal[m-3]*diagonal[m-5] > 0
            cases[str(m)] = {
                "status": "strict all-p positive coefficient certificate",
                "target": "S4^2-S3*S5 = a_(m-4)^2-a_(m-3)*a_(m-5)",
                "normalization": "gap=positive_prefactor*primitive_numerator(t)/primitive_denominator(t)",
                "positive_prefactor": str(prefactor), "numerator": ndata, "denominator": ddata,
                "rational_pair_coprime": True, "primitive_integer_polynomials": True,
                "exact_clearing_identity": True, "strict_denominator_constant": 1,
                "normalized_jets_zero_through_four_convolution_checks": True,
                "queue_enumeration_equals_compact_five_survival_formula": True,
                "probe_p": "2/5", "probe_t": "2/3",
                "S1_through_S5_at_probe": list(map(str, probe_values)),
                "all_five_match_independent_boundary_profiles": True,
                "rational_gap_at_probe": str(gap_probe), "probe_gap_positive": True,
                "actual_negative_gap_found": False,
            }
        return {
            "status": "PASS: actual fourth-gap bases m=5,6,7,8,9, every 0<p<1",
            "parameterization": "p=t/(1+t), q=1/(1+t), t>0",
            "actual_passage_recursion": "P_r(z)=r*p/[z+r*p+(m-r)*q*(1-P_(r+1)(z))]; bottom P_m(z)=m*p/(z+m*p)",
            "jets": "P_r(1-h) coefficients through degree four, normalized by factorials",
            "first_failure_prefix_counts": [1, 1, 2, 5, 14],
            "coefficient_hash_convention": "SHA256 compact JSON ascending nonzero [degree,integer] rows of each primitive polynomial",
            "cases": cases, "all_five_bases_certified": True,
            "positive_numerator_coefficient_count": sum(row["numerator"]["positive_terms"] for row in cases.values()),
            "positive_denominator_coefficient_count": sum(row["denominator"]["positive_terms"] for row in cases.values()),
            "independent_normalization_cases": [{"m": m, "p": "2/5"} for m in (5, 6, 7, 8, 9)],
        }

    def generic_certificate():
        STARTED = time.monotonic()
        STATE = {'status':'IN_PROGRESS',
                 'scope':'all seven degree-six Bernstein coefficients, analytic coverage pending audit',
                 'domain':'t>0,w>=0,u>=0; m=w+10,p=t/(1+t); seventh derivative floors',
                 'parameter_samples':0, 'retries':0, 'value_interval_depth':10,
                 'derivative_floor_depth':7, 'generic_old_w_substitution':'old_w=w+1',
                 'Bernstein_convention':'Cj=[v^j](1+v)^6 H((ell7+U*v)/(1+v)), U=L(1-L)',
                 'coefficients':[], 'completed_indices':[], 'fourth_gap_theorem':False}

        def checkpoint(phase, **values):
            # Retain producer progress state locally, without scratch I/O.
            STATE.update(values)
            STATE["phase"] = phase

        from flint import fmpz_mpoly_ctx
        ctx = fmpz_mpoly_ctx.get(['t','w','N'])
        t,w,N = ctx.gens()
        one,zero = ctx.constant(1),ctx.constant(0)

        class Rat:
            """The reviewed reduced-fraction arithmetic pattern, unchanged."""
            def __init__(self,numerator,denominator=None,reduced=False):
                if isinstance(numerator,Rat):
                    self.n,self.d=numerator.n,numerator.d
                    return
                numerator=numerator if hasattr(numerator,'gcd') else ctx.constant(numerator)
                denominator=(one if denominator is None else denominator
                             if hasattr(denominator,'gcd') else ctx.constant(denominator))
                if denominator.is_zero():raise ZeroDivisionError
                if numerator.is_zero():
                    self.n,self.d=zero,one
                    return
                if not reduced:
                    common=numerator.gcd(denominator)
                    numerator,denominator=numerator//common,denominator//common
                if denominator.coefficient(0)<0:numerator,denominator=-numerator,-denominator
                self.n,self.d=numerator,denominator
            def __add__(self,other):
                other=Rat(other)
                common=self.d.gcd(other.d)
                first,second=self.d//common,other.d//common
                numerator=self.n*second+other.n*first
                cancel=numerator.gcd(common)
                return Rat(numerator//cancel,first*(other.d//cancel),True)
            __radd__=__add__
            def __neg__(self):return Rat(-self.n,self.d,True)
            def __sub__(self,other):return self+-Rat(other)
            def __rsub__(self,other):return Rat(other)+-self
            def __mul__(self,other):
                other=Rat(other)
                first,second=self.n.gcd(other.d),other.n.gcd(self.d)
                return Rat((self.n//first)*(other.n//second),
                           (self.d//second)*(other.d//first),True)
            __rmul__=__mul__
            def __truediv__(self,other):
                other=Rat(other)
                if other.n.is_zero():raise ZeroDivisionError
                first,second=self.n.gcd(other.n),self.d.gcd(other.d)
                return Rat((self.n//first)*(other.d//second),
                           (self.d//second)*(other.n//first),True)
            def __rtruediv__(self,other):return Rat(other)/self
            def __pow__(self,power):
                if power<0:return (Rat(1)/self)**(-power)
                return Rat(self.n**power,self.d**power,True)
            def __eq__(self,other):
                other=Rat(other)
                return self.n*other.d==other.n*self.d

        def stats(poly):
            return {'terms':len(poly),'degrees':list(map(int,poly.degrees())) if poly else [0]*3,
                    'total_degree':int(poly.total_degree()) if poly else 0}

        def fingerprint(poly, stop_negative=False):
            sha=hashlib.sha256();sha.update(b'[')
            count=0;previous_powers=None;anchors=[]
            for index in range(len(poly)-1,-1,-1):
                powers=tuple(map(int,poly.monomial(index)))
                value=int(poly.coefficient(index))
                assert previous_powers is None or previous_powers<powers
                previous_powers=powers
                if stop_negative and value<0:
                    return {'first_negative':{'powers_t_w_u':list(powers),'coefficient':value},
                            'terms_examined_before_failure':count+1}
                if count:sha.update(b',')
                sha.update(json.dumps([*powers,value],separators=(',',':')).encode())
                if powers[1:]==(0,0):anchors.append([powers[0],value])
                count+=1
            sha.update(b']')
            return {'sha256':sha.hexdigest(),'terms_checked':count,'strict_w0_u0_anchors':anchors,
                    'first_negative':None}

        checkpoint('constructing continued-fraction series',effective_soft_limit_bytes=limit)
        # P7(1-h)=N/(1-N*h/(7p)); all four jets are the prescribed floors.
        series={8:[N,zero,zero,zero,zero],7:[7*t,-(1+t)*N,zero,zero,zero]}
        kappa={7:one,6:4+7*t+w-(w+3)*N}
        for birth in range(5,0,-1):
            kappa[birth]=((1+(birth+1)*t+w+9-birth)*kappa[birth+1]
                           -(w+9-birth)*(birth+1)*t*kappa[birth+2])
        for birth in range(6,0,-1):
            later,earlier=series[birth+1],series[birth+2]
            series[birth]=[(1+(birth+1)*t+w+9-birth)*later[j]
                           -((1+t)*later[j-1] if j else zero)
                           -(w+9-birth)*(birth+1)*t*earlier[j] for j in range(5)]
            assert series[birth][0]==7*t*kappa[birth]

        def passage_jets(birth,degree):
            denominator=series[birth]
            numerator=[birth*t*value for value in series[birth+1]]
            numerators=[];result=[]
            for j in range(degree+1):
                value=numerator[j]*denominator[0]**j
                for ell in range(1,j+1):
                    value-=denominator[ell]*numerators[j-ell]*denominator[0]**(ell-1)
                numerators.append(value)
                result.append(Rat(value,denominator[0]**(j+1)))
            return result

        base={i:passage_jets(i,min(4,i-1)) for i in range(1,6)}
        base[6]=passage_jets(6,4)
        p=Rat(t,1+t)
        delta={i:Rat(w+9-i,1+t) for i in range(1,7)}
        f={i:base[i][0] for i in range(1,7)}
        L,ell=base[5][:2]
        e,h,j=[delta[5]*L**2/(5*p)*base[6][k] for k in (2,3,4)]
        assert base[5][2:]==[ell**2/L+e,ell**3/L**2+2*ell*e/L+h,
                            ell**4/L**3+3*ell**2*e/L**2+e**2/L+2*ell*h/L+j]
        upper=L*(1-L)
        width=upper-ell
        assert width==delta[5]*L**2/(5*p)*(1-f[6]-base[6][1])
        checkpoint('seventh floors and symbolic interval ready')

        # Ascending polynomials in Y=X-ell, with exact rational coefficients.
        def trim(poly):
            while len(poly)>1 and poly[-1].n.is_zero():poly.pop()
            return poly
        def padd(a,b):
            return trim([(a[i] if i<len(a) else Rat(0))+
                         (b[i] if i<len(b) else Rat(0)) for i in range(max(len(a),len(b)))])
        def psub(a,b):return padd(a,[-v for v in b])
        def pscale(a,b):return trim([v*b for v in a])
        def pmul(a,b):
            result=[Rat(0) for _ in range(len(a)+len(b)-1)]
            for i,ai in enumerate(a):
                for k,bk in enumerate(b):result[i+k]=result[i+k]+ai*bk
            return trim(result)
        def psum(items):
            result=[Rat(0)]
            for item in items:result=padd(result,item)
            return result
        def peval(poly,value):
            result=Rat(0)
            for coeff in reversed(poly):result=result*value+coeff
            return result

        tail=[None,[ell/L,1/L],[e/L],[h/L],[j/L]]
        jets={5:[[L]]}
        for order in range(1,5):
            jets[5].append(psum(pmul(tail[k],jets[5][order-k]) for k in range(1,order+1)))
            assert jets[5][order][0]==base[5][order]
        for birth in range(4,0,-1):
            later=jets[birth+1];row=[[f[birth]]]
            for order in range(1,birth):
                row.append(pscale(padd(row[order-1],pscale(
                    psum(pmul(later[k],row[order-k]) for k in range(1,order+1)),
                    delta[birth])),f[birth]/(birth*p)))
                assert row[order][0]==base[birth][order]
            jets[birth]=row
            checkpoint('shifted first-step jets ready',birth=birth)

        def words(length):
            def recurse(prefix,total):
                position=len(prefix)
                if position==length-1:
                    yield tuple(prefix+[length-1-total]);return
                for count in range(length-total):
                    if total+count>=position+1:yield from recurse(prefix+[count],total+count)
            yield from recurse([],0)
        queue_words={i:list(words(i)) for i in range(1,6)}
        assert [len(queue_words[i]) for i in range(1,6)]==[1,1,2,5,14]
        def queue_survival(depth,table):
            result=[Rat(1)]
            for length in range(1,depth+1):
                for word in queue_words[length]:
                    term=[Rat(1)]
                    for offset,order in enumerate(word):term=pmul(term,table[depth-offset][order])
                    result=psub(result,term)
            return result
        S={}
        for depth in range(1,6):
            S[depth]=queue_survival(depth,jets)
            checkpoint('direct queue polynomial ready',depth=depth,Y_degree=len(S[depth])-1)
        sextic=psub(pmul(S[4],S[4]),pmul(S[3],S[5]))
        assert len(sextic)==7 and not sextic[-1].n.is_zero()

        # Independent absorption recursion: no first-failure word enumeration.
        extinction=[[Rat(1)]]
        for depth in range(1,6):
            nxt=[[Rat(1)]]
            for queue in range(1,depth+1):
                nxt.append(psum(pmul(jets[depth][n],extinction[queue-1+n])
                                for n in range(depth-queue+1)))
            extinction=nxt
            assert psub([Rat(1)],extinction[1])==S[depth]
        assert fingerprint(sextic[0].n)['sha256']=='773a59581c2ee9031259d2358e7f15662be7715d803c2f0f1eb41cd53fab4905'
        STATE['exact_checks']={'queue_words_equal_absorption_recursion':True,
            'seventh_floor_jets_match_continued_fraction':True,
            'lower_endpoint_generic_numerator_hash_matches':True,
            'symbolic_interval_width_identity':True,
            'inverse_Bernstein_reconstruction':False,
            'upper_endpoint_direct_queue_identity':False}
        STATE['shifted_sextic_coefficients']=[{'Y_degree':k,'numerator':stats(v.n),
            'denominator':stats(v.d),'numerator_sha256':fingerprint(v.n)['sha256'],
            'denominator_sha256':fingerprint(v.d)['sha256']} for k,v in enumerate(sextic)]
        checkpoint('degree-six H polynomial independently checked',sextic_degree=6)

        outctx=fmpz_mpoly_ctx.get(['t','w','u'])
        tt,ww,uu=outctx.gens();unit=outctx.constant(1)
        b10=1+11*tt+ww
        b9=(2+10*tt+ww)*b10-10*tt*(ww+1)
        b8=(3+9*tt+ww)*b9-9*tt*(ww+2)*b10
        b7=(4+8*tt+ww)*b8-8*tt*(ww+3)*b9
        NN=7*tt*((1+8*tt)*b8+uu*b7)
        DD=(1+uu)*(1+8*tt)*b7
        low_n,low_d=outctx.constant(0),unit
        for birth in (10,9,8,7):
            low_n,low_d=(birth*tt*low_d,
                (1+(birth+1)*tt+ww+10-birth)*low_d-(ww+10-birth)*low_n)
        assert NN==low_n*(1+8*tt)+uu*7*tt*low_d
        assert DD==(1+uu)*(1+8*tt)*low_d
        assert all(poly.coefficient(j)>0 for poly in (b10,b9,b8,b7,NN,DD) for j in range(len(poly)))

        def homogeneous(poly,degree,progress=False):
            blocks=[{} for _ in range(degree+1)]
            for index in range(len(poly)):
                i,j,k=map(int,poly.monomial(index))
                # Exactly compose old_w=w+1 before rational N substitution.
                from math import comb
                value=int(poly.coefficient(index))
                for power in range(j+1):
                    key=(i,power,0)
                    blocks[k][key]=blocks[k].get(key,0)+value*comb(j,power)
            value=outctx.from_dict(blocks[degree]);blocks[degree].clear()
            power=unit
            for k in range(degree-1,-1,-1):
                power*=DD
                part=outctx.from_dict(blocks[k]);blocks[k].clear()
                value=value*NN+part*power
                if progress:checkpoint('homogeneous Horner step',remaining_N_power=k,partial_terms=len(value))
            return value

        continuants=[]
        for birth in range(1,7):
            fact=homogeneous(kappa[birth],1)
            assert all(fact.coefficient(j)>0 for j in range(len(fact)))
            assert fact(0,0,0)>0
            continuants.append({'birth':birth,**stats(fact),'positive_constant':int(fact(0,0,0)),
                                'coefficient_sha256':fingerprint(fact)['sha256']})
        STATE['positive_continuants']=continuants
        scaled=[sextic[k]*width**k for k in range(7)]
        coefficients=[]
        for index in range(7):
            checkpoint('constructing Bernstein coefficient',Bernstein_index=index)
            # H(ell+width*z)=sum scaled[k] z^k, then z=v/(1+v).
            coefficient=sum((scaled[k]*comb(6-k,index-k) for k in range(index+1)),Rat(0))
            coefficients.append(coefficient)
            remaining=coefficient.d;powers=[]
            for birth in range(1,7):
                power=0
                while True:
                    quotient,remainder=divmod(remaining,kappa[birth])
                    if remainder:break
                    remaining=quotient;power+=1
                powers.append(power)
            extra={}
            for name,factor in [('t',t),('1+t',1+t)]:
                power=0
                while True:
                    quotient,remainder=divmod(remaining,factor)
                    if remainder:break
                    remaining=quotient;power+=1
                extra[name]=power
            assert remaining.total_degree()==0 and int(remaining.coefficient(0))>0
            num_degree=int(coefficient.n.degrees()[2])
            den_degree=int(coefficient.d.degrees()[2])
            assert den_degree==sum(powers)
            degree=max(num_degree,den_degree)
            generic={'index':index,'numerator':stats(coefficient.n),'denominator':stats(coefficient.d),
                'numerator_sha256':fingerprint(coefficient.n)['sha256'],
                'denominator_sha256':fingerprint(coefficient.d)['sha256'],
                'denominator_constant':int(remaining.coefficient(0)),
                'continuant_powers':powers,'extra_positive_powers':extra,
                'numerator_N_degree':num_degree,'denominator_N_degree':den_degree,
                'homogenizing_N_degree':degree,'extra_D_power_in_denominator':degree-den_degree}
            checkpoint('generic Bernstein coefficient ready',Bernstein_index=index,
                       generic_num_terms=len(coefficient.n),generic_N_degrees=[num_degree,den_degree],
                       denominator_continuant_powers=powers)
            numerator=homogeneous(coefficient.n,degree,True)
            certificate=fingerprint(numerator,True)
            if certificate['first_negative'] is not None:
                STATE['failed_coefficient']={'generic':generic,'numerator':stats(numerator),'certificate':certificate}
                checkpoint('first negative coefficient; stopped',status='NEGATIVE_COEFFICIENT_METHOD_OBSTRUCTION',
                           Bernstein_index=index,first_negative=certificate['first_negative'],
                           actual_counterexample=False)
                return
            assert certificate['strict_w0_u0_anchors']
            if index==0:
                assert certificate['sha256']=='325ebe4c94f70f792b7e846577f1ea79b78e1e63bc43584991230d36330d86de'
                STATE['exact_checks']['C0_homogeneous_hash_matches']=True
            # Reverse the binomial transform in every completed degree.
            reconstructed=sum((coefficients[k]*((-1)**(index-k))*comb(6-k,index-k)
                               for k in range(index+1)),Rat(0))
            assert reconstructed==scaled[index]
            record={'index':index,'generic':generic,'numerator':stats(numerator),
                    'certificate':certificate,'inverse_Bernstein_coefficient_check':True}
            STATE['coefficients'].append(record)
            STATE['completed_indices'].append(index)
            checkpoint('Bernstein coefficient passed',Bernstein_index=index,
                       positive_terms=len(numerator),coefficient_sha256=certificate['sha256'])

        # Endpoint U via a separate scalar first-step/queue calculation.
        upperjets={5:[L,upper,upper**2/L+e,upper**3/L**2+2*upper*e/L+h,
                      upper**4/L**3+3*upper**2*e/L**2+e**2/L+2*upper*h/L+j]}
        for birth in range(4,0,-1):
            later=upperjets[birth+1];row=[f[birth]]
            for order in range(1,birth):
                row.append(f[birth]/(birth*p)*(row[order-1]+delta[birth]*sum(
                    (later[k]*row[order-k] for k in range(1,order+1)),Rat(0))))
            upperjets[birth]=row
        direct={i:queue_survival(i,{k:[[v] for v in row] for k,row in upperjets.items()})[0]
                for i in (3,4,5)}
        uppergap=direct[4]**2-direct[3]*direct[5]
        assert coefficients[6]==uppergap==peval(sextic,width)
        STATE['exact_checks']['upper_endpoint_direct_queue_identity']=True
        STATE['exact_checks']['inverse_Bernstein_reconstruction']=True
        STATE['total_positive_terms']=sum(r['certificate']['terms_checked'] for r in STATE['coefficients'])
        checkpoint('all seven passed; pending independent and analytic audit',status='PASS_ALL7_FORMAL_CERTIFICATES',
                   total_positive_terms=STATE['total_positive_terms'],
                   all_seven_Bernstein_claim=True,fourth_gap_theorem=False)
        assert STATE["status"] == "PASS_ALL7_FORMAL_CERTIFICATES"
        assert STATE["completed_indices"] == list(range(7))
        assert STATE["total_positive_terms"] == 485204
        compact_coefficients = []
        for record in STATE["coefficients"]:
            record = dict(record)
            certificate = dict(record["certificate"])
            anchors = certificate.pop("strict_w0_u0_anchors")
            assert anchors and anchors[0][1] > 0
            certificate["strict_w0_u0_lowest_t_anchor"] = anchors[0]
            certificate["strict_w0_u0_anchor_count"] = len(anchors)
            record["certificate"] = certificate
            compact_coefficients.append(record)
        return {
            "status": "PASS_ALL7_FORMAL_CERTIFICATES",
            "domain": STATE["domain"],
            "free_value": "N=P7(1)",
            "generic_old_w_substitution": STATE["generic_old_w_substitution"],
            "value_interval_depth": 10, "derivative_floor_depth": 7,
            "Bernstein_convention": STATE["Bernstein_convention"],
            "sextic_degree": STATE["sextic_degree"],
            "completed_indices": STATE["completed_indices"],
            "coefficients": compact_coefficients,
            "total_positive_terms": STATE["total_positive_terms"],
            "positive_continuants": STATE["positive_continuants"],
            "shifted_sextic_coefficients": STATE["shifted_sextic_coefficients"],
            "exact_checks": STATE["exact_checks"],
            "coefficient_hash_convention": "SHA256 compact JSON lexicographically ascending nonzero [t-degree,w-degree,u-degree,integer] rows; generic hashes use N in place of u",
            "strict_anchor_convention": "lowest t-degree monomial at w=u=0; its positive coefficient gives strict positivity for t>0,w,u>=0",
            "full_coefficient_rows_saved": False,
            "parameter_samples": 0, "root_calculations": 0,
            "elapsed_seconds": round(time.monotonic()-STARTED, 3),
        }

    started = time.monotonic()
    previous = resource.getrlimit(resource.RLIMIT_AS)
    limit = min([1900*1024**2]+[v for v in previous if v != resource.RLIM_INFINITY])
    resource.setrlimit(resource.RLIMIT_AS, (limit, previous[1]))
    try:
        try:
            generic = generic_certificate()
        except ImportError as error:
            raise SystemExit("--survival-fourth-lc requires optional python-flint") from error
        assert generic is not None, "Formal coefficient certificate did not complete"
        bases = actual_small_bases()
        assert bases["positive_numerator_coefficient_count"] == 235
        assert bases["positive_denominator_coefficient_count"] == 240
        for case in bases["cases"].values():
            for part in ("numerator", "denominator"):
                # Retain the exact hashes, degrees and strict anchors, not rows.
                del case[part]["ascending_coefficients"]
        result = {
            "status": "PASS: all-seven formal Bernstein certificate and five actual all-p bases",
            "model": "actual infinite-horizon Ehrenfest survival, boundary rate one; births p*(m-y), deaths (1-p)*y",
            "target": "S4^2-S3*S5 = a_(m-4)^2-a_(m-3)*a_(m-5), S_r=a_(m-r)",
            "frozen_producer_source_sha256": "71fda83f0010158a7570d705e4116758d5342b42a07aebdb0f48e6f35f2befd0",
            "generic_m_ge_10": generic,
            "actual_small_bases": bases,
            "analytic_theorem_status": "PROVED: independent reconstruction, actual bases, and analytic coverage audits passed; see report (55q)",
            "conclusion": "Every integer m>=5 and every 0<p<1 satisfies the strict fourth infinite-horizon survival LC gap, with boundary rate one",
            "coverage_dependencies": [
                "actual-first-jet reduction at m>=8 using strong preceding formal gaps and response growth; seventh floors are coordinatewise admissible",
                "depth-ten lower and exponential upper bounds for actual P7(1), parameterized by u>=0 after m=w+10,p=t/(1+t)",
                "actual first jet lies in [ell7,L*(1-L)]; positive transformed denominators and seven strict Bernstein coefficients cover the full interval",
                "actual all-p bases m=5,6,7,8,9 cover the remaining sizes",
            ],
            "fourth_gap_theorem_promoted": True,
            "full_survival_LC_finite_time_LC_stationary_ULC_or_NA_claim": False,
            "old_certificate_helpers_called": [],
            "boundary_profile_normalization_calls": 5,
            "old_helper_metadata_modified": False,
            "parameter_grid_or_root_calls": False,
            "effective_soft_limit_bytes": limit,
            "elapsed_seconds": round(time.monotonic()-started, 3),
            "peak_rss_MiB": round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024, 2),
        }
    finally:
        resource.setrlimit(resource.RLIMIT_AS, previous)
    result["resource_limits_restored"] = resource.getrlimit(resource.RLIMIT_AS) == previous
    return result


def survival_slack_response_certificate():
    """One r=6 marked-word check of the all-depth slack-response formulas.

    The all-depth analytic lemma uses actual Ehrenfest zeroths and arbitrary
    nonnegative formal reciprocal slacks, propagated by first step. Its
    proof is separate from these finite symbolic identity checks. No gap
    polynomial, parameter sample, or earlier certificate helper is used.
    """
    import resource
    import time

    started = time.monotonic()
    previous = resource.getrlimit(resource.RLIMIT_AS)
    limit = min([1900*1024**2]+[v for v in previous if v != resource.RLIM_INFINITY])
    resource.setrlimit(resource.RLIMIT_AS, (limit, previous[1]))
    try:
        import sympy as sp
        r, n, k = 6, 5, 3
        f = dict(enumerate(sp.symbols('f1:7'),1))
        u = dict(enumerate(sp.symbols('u1:7'),1))
        v = dict(enumerate(sp.symbols('v1:7'),1))
        rho = dict(enumerate(sp.symbols('rho1:6'),1))
        def Q(j):
            return sp.prod(f[i] for i in range(1,j+1))
        def words(length):
            def recurse(prefix,left):
                if len(prefix)==length:
                    if left==0: yield tuple(prefix)
                    return
                for count in range(left+1):
                    row=prefix+[count]
                    if len(row)<length and sum(row)<len(row): continue
                    yield from recurse(row,left-count)
            yield from recurse([],length-1)
        def marked_jet(level,order):
            a=sp.prod(rho[j] for j in range(level,r))
            offset=order-k
            assert offset in (0,1,2)
            if offset==0: return a
            if offset==1: return 2*a*sum(u[j] for j in range(level,r+1))
            return a*(2*sum(v[j] for j in range(level,r+1))
                     +sum(u[j]**2 for j in range(level,r+1))
                     +4*sum(u[i]*u[j] for i in range(level,r+1)
                            for j in range(i+1,r+1)))
        def ordinary_jet(level,order):
            assert order in (0,1,2)
            return f[level]*(1 if order==0 else u[level] if order==1 else v[level])
        def response(depth):
            result=sp.Integer(0)
            marked_count=0
            for length in range(1,depth+1):
                for word in words(length):
                    for marked,order in enumerate(word):
                        if order<k: continue
                        term=marked_jet(depth-marked,order)
                        for pos,count in enumerate(word):
                            if pos!=marked:term*=ordinary_jet(depth-pos,count)
                        result+=term
                        marked_count+=1
            return sp.expand(result),marked_count
        theta=rho[n]/f[n]
        eta=rho[n-1]*f[n]/f[n-1]
        B0=rho[n]*Q(n-1)
        Z1=1/f[1]+sum(u[i] for i in range(2,n+1))
        Z2=(1/(f[1]*f[2])+sum(u[i] for i in range(3,n+1))/f[1]
            +sum(v[i] for i in range(3,n+1))
            +sum(u[i]*u[j] for i in range(2,n+1) for j in range(i+1,n+1)))
        alpha=rho[n-1]*rho[n]*Q(n-2)
        beta=B0*(Z1+(1+eta)*u[n]+2*u[r])
        gamma=(f[r]*u[r]*beta+f[r]*v[r]*B0*(1+eta)
               +Q(n)*(Z2+2*u[r]*Z1+2*v[r]+u[r]**2))
        records=[]
        for depth,formula in ((r-2,alpha),(r-1,beta),(r,gamma)):
            direct,count=response(depth)
            assert sp.expand(direct-formula)==0
            records.append({'depth':depth,'marked_word_terms':count,
                            'expanded_response_terms':len(sp.Add.make_args(direct)),
                            'formula_matches':True})
        th,et,a2,b2=sp.symbols('theta eta a2 b2')
        Btilde=Z1+(1+et)*u[n]
        constant=(1-2*th*f[2])/(f[1]*f[2])-2*th*a2
        remainder=(constant+(1/f[1]-2*th*(1+b2))*u[3]
                   +(1/f[1]-2*th)*sum(u[i] for i in range(4,n))
                   +(1/f[1]-2*th*(2+et))*u[n]
                   +sum(v[i] for i in range(3,n+1))
                   +sum(u[i]*u[j] for i in range(2,n+1) for j in range(i+1,n+1)))
        assert sp.expand((Z2-2*th*Btilde-remainder).subs(u[2],a2+b2*u[3]))==0
        beta_normal=Z1+(1+eta)*u[n]+2*u[r]
        response_decomposition=(B0/theta*(Z2-2*theta*(Z1+(1+eta)*u[n])
             +2*u[r]*(Z1-2*theta)+2*v[r]+u[r]**2
             +theta*f[r]*u[r]*beta_normal+theta*(1+eta)*f[r]*v[r]))
        assert sp.expand(gamma-2*beta-response_decomposition)==0
        result = {
            "status": "PASS: r=6 marked-word identities and two exact remainder identities",
            "scope": "one abstract r=6 response index check",
            "source_checker_sha256": "00b63807418fbf45fa529842d2aa72d528e1ce732ac2e5d62baf46f3aead2386",
            "checked_indices": {"r": r, "n": n, "k": k},
            "prefix_counts": [len(list(words(i))) for i in range(1, 7)],
            "response_checks": records,
            "linear_remainder_identity": True,
            "full_response_remainder_identity": True,
            "all_depth_conclusion_from_analytic_proof_not_this_test": True,
            "analytic_lemma": {
                "status": "PROVED: all-depth marked-word argument and strict margin passed independent analytic audit",
                "proof_location": "report (55r)",
                "hypotheses": "integers m>=r>=6; 0<p<1; q=1-p; boundary rate one; actual Ehrenfest zeroths f_i=P_i(1)",
                "formal_scope": "arbitrary nonnegative reciprocal slacks at passage r, with earlier passages propagated by first step",
                "reciprocal_convention": "P_r(1-h)=f_r/(1-sum_j sigma_j*h^j/f_r)",
                "definitions": {
                    "n": "r-1", "k": "r-3",
                    "S_j": "queue survival polynomial; at actual jets a_(m-j)",
                    "rho_i": "(m-i)*q*f_i^2/(i*p)",
                    "Q_j": "product_(i=1)^j f_i",
                    "u_i": "[h]P_i(1-h)/f_i", "v_i": "[h^2]P_i(1-h)/f_i",
                    "theta": "rho_n/f_n", "eta": "rho_(n-1)*f_n/f_(n-1)",
                    "B_star": "rho_n*Q_(n-1)",
                    "A": "-partial_(sigma_k) S_(r-2)",
                    "B": "-partial_(sigma_k) S_(r-1)",
                    "C": "-partial_(sigma_k) S_r",
                    "Z1": "1/f1+sum_(i=2)^n u_i",
                    "Z2": "1/(f1*f2)+sum_(i=3)^n u_i/f1+sum_(i=3)^n v_i+sum_(2<=i<j<=n)u_i*u_j",
                },
                "response_formulas": [
                    "A=rho_(n-1)*rho_n*Q_(n-2)",
                    "B=B_star*(Z1+(1+eta)*u_n+2*u_r)",
                    "C=f_r*u_r*B+f_r*v_r*B_star*(1+eta)+Q_n*(Z2+2*u_r*Z1+2*v_r+u_r^2)",
                ],
                "strict_comparison": "C-2*B > (B_star/p)*(1/f2-f2) > 0",
                "responses_independent_of_sigma_k": "r>=6 gives k>=3, while only normalized first and second jets occur",
                "actual_zeroth_sign_input": "i*p*I_(i-1)-(m+1-i)*q*I_i=1; in particular 1/f1>6*theta for n>=5",
                "conditional_gap_derivative": "For U=S_(r-2)>0, V=S_(r-1)<=U and W=S_r>=0: partial_(sigma_k)(V^2-U*W)=A*W+2*B*(U-V)+U*(C-2*B)>0",
                "remaining_profile_requirement": "A finite formal lowering path still requires positivity and V<=U throughout; actual LC alone does not supply them",
                "arbitrary_formal_profiles_asserted_positive_or_ordered": False,
                "all_depth_sign_inferred_from_r6_check": False,
            },
            "parameter_evaluations": 0, "gap_polynomials_constructed": 0,
            "size_or_rate_sweep": False, "fifth_gap_certificate": False,
            "full_survival_LC_stationary_ULC_or_NA_claim": False,
            "old_certificate_helpers_called": [],
            "effective_soft_limit_bytes": limit,
            "seconds": round(time.monotonic()-started, 4),
            "peak_rss_mib": round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024, 2),
        }
    finally:
        resource.setrlimit(resource.RLIMIT_AS, previous)
    result["limits_restored"] = resource.getrlimit(resource.RLIMIT_AS) == previous
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", help="Optional JSON certificate file")
    parser.add_argument("--symbolic-n", type=int,
                        help="Check every sleep rate and ULC input at this complete-graph size")
    parser.add_argument("--renewal", action="store_true",
                        help="Check the passage-count renewal formula and rank-flow obstruction")
    parser.add_argument("--reinjection-logcurvature", action="store_true",
                        help="Prove stationary ULC-ratio monotonicity at n=3,4")
    parser.add_argument("--reverse-boundary", action="store_true",
                        help="Audit the scalar boundary/sensitivity reverse-mean recursion")
    parser.add_argument("--rank-top-structure", action="store_true",
                        help="Check the positive top-rank kernel identity and one rational formula case")
    parser.add_argument("--rank-likelihood-order", action="store_true",
                        help="Prove strict stationary LR order at fixed n=5,p=1/2 for every feedback rate")
    parser.add_argument("--rank-top-coefficients", action="store_true",
                        help="Prove top-row coefficient TP2 through all-index symbolic mixed-minor tables")
    parser.add_argument("--rank-next-coefficients", action="store_true",
                        help="Check proved two-positive-support cases and low-degree next-row restoration")
    parser.add_argument("--rank-next-pointwise", action="store_true",
                        help="Prove pointwise next-row LR for arbitrary two positive support indices")
    parser.add_argument("--rank-next-pointwise-all", action="store_true",
                        help="Certify three/four-index pointwise LR closure (requires optional python-flint)")
    parser.add_argument("--rank-dimension-transfer", action="store_true",
                        help="Check the positive dimension identity and exact generic ULC obstruction at n=2 -> 3")
    parser.add_argument("--survival-top-lc", action="store_true",
                        help="Prove the top infinite-horizon survival LC gap at every size and rate")
    parser.add_argument("--survival-next-lc", action="store_true",
                        help="Prove the next infinite-horizon survival LC gap at every size and rate")
    parser.add_argument("--survival-third-reduction", action="store_true",
                        help="Check the third-gap jet reduction and exact endpoint obstruction")
    parser.add_argument("--survival-third-fifth-floors", action="store_true",
                        help="Check fifth-passage jet floors and the positive c=2 limiting interval")
    parser.add_argument("--survival-third-lc", action="store_true",
                        help="Prove the third infinite-horizon survival LC gap at every size and rate (requires optional python-flint)")
    parser.add_argument("--survival-high-slacks", action="store_true",
                        help="Check general high-jet identities and the all-rate m=7 formal endpoint")
    parser.add_argument("--survival-fourth-reduction", action="store_true",
                        help="Check fourth-gap response reduction, the c=3 interval, and actual m=5,...,8 bases")
    parser.add_argument("--survival-fourth-lc", action="store_true",
                        help="Check all seven generic fourth-gap coefficients and actual m=5,...,9 bases (requires optional python-flint)")
    parser.add_argument("--survival-slack-response", action="store_true",
                        help="Check r=6 marked-word responses and two exact remainder identities")
    args = parser.parse_args()
    if sum((args.renewal, args.reinjection_logcurvature, args.reverse_boundary,
            args.rank_top_structure, args.rank_likelihood_order, args.rank_top_coefficients,
            args.rank_next_coefficients, args.rank_next_pointwise, args.rank_next_pointwise_all,
            args.rank_dimension_transfer, args.survival_top_lc, args.survival_next_lc, args.survival_third_reduction, args.survival_third_fifth_floors,
            args.survival_third_lc, args.survival_high_slacks, args.survival_fourth_reduction, args.survival_fourth_lc,
            args.survival_slack_response,
            args.symbolic_n is not None)) > 1:
        parser.error("Choose only one specialized calculation")
    if args.survival_slack_response:
        result = survival_slack_response_certificate()
        if args.output:
            from pathlib import Path
            output = Path(args.output)
            saved = json.loads(output.read_text()) if output.exists() else {}
            assert isinstance(saved, dict), "Use the symbolic certificate object, not the rank-case list"
            saved["survival_slack_response"] = result
            output.write_text(json.dumps(saved, indent=2)+"\n")
        print("r=6 slack-response identities passed; the all-depth sign uses the separate analytical proof (55r)", flush=True)
        return
    if args.survival_fourth_lc:
        result = survival_fourth_lc_certificate()
        if args.output:
            from pathlib import Path
            output = Path(args.output)
            saved = json.loads(output.read_text()) if output.exists() else {}
            assert isinstance(saved, dict), "Use the symbolic certificate object, not the rank-case list"
            saved["survival_fourth_lc"] = result
            output.write_text(json.dumps(saved, indent=2)+"\n")
        print("Fourth survival gap: all seven uniform coefficients and five actual bases passed; analytical proof and independent audits in report (55q)", flush=True)
        return
    if args.survival_fourth_reduction:
        result = survival_fourth_reduction_certificate()
        if args.output:
            from pathlib import Path
            output = Path(args.output)
            saved = json.loads(output.read_text()) if output.exists() else {}
            assert isinstance(saved, dict), "Use the symbolic certificate object, not the rank-case list"
            saved["survival_fourth_reduction"] = result
            output.write_text(json.dumps(saved, indent=2)+"\n")
        print("Fourth-gap response, c=3 interval and four actual bases passed; all-size fourth LC remains open", flush=True)
        return
    if args.survival_high_slacks:
        result = survival_high_slacks_certificate()
        if args.output:
            from pathlib import Path
            output = Path(args.output)
            saved = json.loads(output.read_text()) if output.exists() else {}
            assert isinstance(saved, dict), "Use the symbolic certificate object, not the rank-case list"
            saved["survival_high_slacks"] = result
            output.write_text(json.dumps(saved, indent=2)+"\n")
        print("High-slack identities and all-p m=7 endpoint passed; fourth survival LC remains open", flush=True)
        return
    if args.survival_third_lc:
        result = survival_third_lc_certificate()
        if args.output:
            from pathlib import Path
            output = Path(args.output)
            saved = json.loads(output.read_text()) if output.exists() else {}
            assert isinstance(saved, dict), "Use the symbolic certificate object, not the rank-case list"
            saved["survival_third_lc"] = result
            output.write_text(json.dumps(saved, indent=2)+"\n")
        print("All-size third survival LC passed: 92,592 positive coefficients and four exact actual bases", flush=True)
        return
    if args.survival_third_fifth_floors:
        result = survival_third_fifth_floors_certificate()
        if args.output:
            from pathlib import Path
            output = Path(args.output)
            saved = json.loads(output.read_text()) if output.exists() else {}
            assert isinstance(saved, dict), "Use the symbolic certificate object, not the rank-case list"
            saved["survival_third_fifth_floors"] = result
            output.write_text(json.dumps(saved, indent=2)+"\n")
        print("Fifth/sixth passage floors, c=2,3 limiting intervals and bases passed; use --survival-third-lc for the all-size theorem", flush=True)
        return
    if args.survival_third_reduction:
        result = survival_third_reduction_certificate()
        if args.output:
            from pathlib import Path
            output = Path(args.output)
            saved = json.loads(output.read_text()) if output.exists() else {}
            assert isinstance(saved, dict), "Use the symbolic certificate object, not the rank-case list"
            saved["survival_third_reduction"] = result
            output.write_text(json.dumps(saved, indent=2)+"\n")
        print("Third survival-gap reduction and endpoint obstruction checks passed; use --survival-third-lc for the all-size theorem", flush=True)
        return
    if args.survival_next_lc:
        result = survival_next_lc_certificate()
        if args.output:
            from pathlib import Path
            output = Path(args.output)
            saved = json.loads(output.read_text()) if output.exists() else {}
            assert isinstance(saved, dict), "Use the symbolic certificate object, not the rank-case list"
            saved["survival_next_lc"] = result
            output.write_text(json.dumps(saved, indent=2)+"\n")
        print("All-size next survival LC certificate and independent resolvent audits passed", flush=True)
        return
    if args.survival_top_lc:
        result = survival_top_lc_certificate()
        if args.output:
            from pathlib import Path
            output = Path(args.output)
            saved = json.loads(output.read_text()) if output.exists() else {}
            assert isinstance(saved, dict), "Use the symbolic certificate object, not the rank-case list"
            saved["survival_top_lc"] = result
            output.write_text(json.dumps(saved, indent=2)+"\n")
        print("All-size top survival LC certificate and independent resolvent audits passed", flush=True)
        return
    if args.rank_dimension_transfer:
        result = rank_dimension_transfer_certificate()
        if args.output:
            from pathlib import Path
            output = Path(args.output)
            saved = json.loads(output.read_text()) if output.exists() else {}
            assert isinstance(saved, dict), "Use the symbolic certificate object, not the rank-case list"
            saved["rank_dimension_transfer"] = result
            output.write_text(json.dumps(saved, indent=2)+"\n")
        print("All 27 dimension-transfer and nonstationary ULC-obstruction checks passed", flush=True)
        return
    if args.rank_next_pointwise_all:
        result = rank_next_pointwise_all_certificate()
        if args.output:
            from pathlib import Path
            output = Path(args.output)
            saved = json.loads(output.read_text()) if output.exists() else {}
            assert isinstance(saved, dict), "Use the symbolic certificate object, not the rank-case list"
            saved["rank_next_pointwise_all"] = result
            output.write_text(json.dumps(saved, indent=2)+"\n")
        print("All three/four-index pointwise LR certificates and direct audits passed", flush=True)
        return
    if args.rank_next_pointwise:
        result = rank_next_pointwise_certificate()
        if args.output:
            from pathlib import Path
            output = Path(args.output)
            saved = json.loads(output.read_text()) if output.exists() else {}
            assert isinstance(saved, dict), "Use the symbolic certificate object, not the rank-case list"
            saved["rank_next_pointwise"] = result
            output.write_text(json.dumps(saved, indent=2)+"\n")
        print("All-index two-positive-support pointwise LR and single direct audit passed", flush=True)
        return
    if args.rank_next_coefficients:
        result = rank_next_coefficients_certificate()
        if args.output:
            from pathlib import Path
            output = Path(args.output)
            saved = json.loads(output.read_text()) if output.exists() else {}
            assert isinstance(saved, dict), "Use the symbolic certificate object, not the rank-case list"
            saved["rank_next_coefficients"] = result
            output.write_text(json.dumps(saved, indent=2)+"\n")
        print("Proved next-row support cases and 21 all-index low-degree restoration checks passed", flush=True)
        return
    if args.rank_top_coefficients:
        result = rank_top_coefficients_certificate()
        if args.output:
            from pathlib import Path
            output = Path(args.output)
            saved = json.loads(output.read_text()) if output.exists() else {}
            assert isinstance(saved, dict), "Use the symbolic certificate object, not the rank-case list"
            saved["rank_top_coefficients"] = result
            output.write_text(json.dumps(saved, indent=2)+"\n")
        print("All-size top-row symbolic coefficient kernels and single formula audit passed", flush=True)
        return
    if args.rank_likelihood_order:
        result = rank_likelihood_order_certificate()
        if args.output:
            from pathlib import Path
            output = Path(args.output)
            saved = json.loads(output.read_text()) if output.exists() else {}
            assert isinstance(saved, dict), "Use the symbolic certificate object, not the rank-case list"
            saved["rank_likelihood_order"] = result
            output.write_text(json.dumps(saved, indent=2)+"\n")
        print("Fixed n=5,p=1/2 all-feedback strict LR coefficient checks passed", flush=True)
        return
    if args.rank_top_structure:
        result = rank_top_structure_certificate()
        if args.output:
            from pathlib import Path
            output = Path(args.output)
            saved = json.loads(output.read_text()) if output.exists() else {}
            assert isinstance(saved, dict), "Use the symbolic certificate object, not the rank-case list"
            saved["rank_top_structure"] = result
            output.write_text(json.dumps(saved, indent=2)+"\n")
        print("Top-rank symbolic kernel and rational formula checks passed", flush=True)
        return
    if args.reverse_boundary:
        result = reverse_boundary_certificate()
        if args.output:
            from pathlib import Path
            output = Path(args.output)
            saved = json.loads(output.read_text()) if output.exists() else {}
            assert isinstance(saved, dict), "Use the symbolic certificate object, not the rank-case list"
            saved["reverse_boundary"] = result
            output.write_text(json.dumps(saved, indent=2)+"\n")
        print("Scalar reverse-boundary and conditional-moment checks passed", flush=True)
        return
    if args.reinjection_logcurvature:
        result = reinjection_logcurvature_certificate()
        if args.output:
            from pathlib import Path
            output = Path(args.output)
            saved = json.loads(output.read_text()) if output.exists() else {}
            assert isinstance(saved, dict), "Use the symbolic certificate object, not the rank-case list"
            saved["reinjection_logcurvature"] = result
            output.write_text(json.dumps(saved, indent=2)+"\n")
        print("Reinjection log-curvature and obstruction checks passed", flush=True)
        return
    if args.renewal:
        result = renewal_certificate()
        print(json.dumps(result, indent=2), flush=True)
        if args.output:
            from pathlib import Path
            output = Path(args.output)
            saved = json.loads(output.read_text()) if output.exists() else {}
            assert isinstance(saved, dict), "Use the symbolic certificate object, not the rank-case list"
            saved["renewal_representation"] = result
            output.write_text(json.dumps(saved, indent=2)+"\n")
        return
    if args.symbolic_n is not None:
        assert args.symbolic_n >= 2
        result = symbolic_lc_certificate(args.symbolic_n)
        print(args.symbolic_n, "all-rate, all-input ULC", result["all_rates_ULC_preservation_certified"],
              "failed coefficient tests", len(result["failed_coefficient_tests"]), flush=True)
        if args.output:
            from pathlib import Path
            output = Path(args.output)
            saved = json.loads(output.read_text()) if output.exists() else {"scope": "finite complete-graph sizes, all positive sleep rates", "cases": {}}
            assert saved.get("scope") == "finite complete-graph sizes, all positive sleep rates"
            saved["cases"][str(args.symbolic_n)] = result
            output.write_text(json.dumps(saved, indent=2) + "\n")
        return
    result = []
    # First match raw stabilization, then probe the larger symmetric cases.
    for n, rate in [(2, 1), (3, 1), (4, 1), (5, 1), (8, 1), (12, 1),
                    (8, sp.Rational(1, 64)), (8, 64)]:
        case = check(n, rate)
        result.append(case)
        print(n, "lambda", rate, "ULC", case["ULC"],
              "real roots", case["real_roots"], "all-input ULC",
              case["all_active_ULC_preservation"]["ULC_preservation_certified"], flush=True)
    if args.output:
        from pathlib import Path
        Path(args.output).write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()
