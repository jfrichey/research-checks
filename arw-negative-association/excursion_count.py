"""Exact diagnostics for the complete-graph count-excursion criterion.

The default mode preserves the original finite criterion checks. Use
--allsize-structure to check the count recurrences and determinant formulas,
the positive-rate obstruction to monotonicity of H, and the zero-rate
forward-tail log-concavity obstruction. Use --continuant-positivity for the
capture/top-entry identities and residuals supporting the analytic proof
of C[m-1,m]>0. Use --interior-pqd for the general-size polynomial identities
proving the r=1 stationary PQD cut and five saved reverse-kernel diagnostics.
Use --reverse-monotonicity for the all-rate reverse P/Q proof at m<=3.
Use --highrate-limit for the stationary vacancy coefficients, full excursion
limit, positive commutator, queue-word formulas and factorial-volume
diagnostics, all using the same independent m=3 count Taylor case.
Use --planar-capture for the bulk Green factorization, targeted birth-death
and gradient-killing obstructions, the Ehrenfest boundary identities, and
the finite-rate commutator remainder.
Positivity of the full C matrix remains open.
"""
import argparse
import json
from pathlib import Path

from itertools import combinations
from algebra import sp, graph_kernel, vertex_excursion_operators, trace_kernel
from exact import absorption, add, stable, stationary_law


def count_operators(n, rate):
    K = graph_kernel(tuple(tuple(j for j in range(n+1) if j != i)
                           for i in range(n)))
    O = vertex_excursion_operators(K, (1/sp.Rational(rate),)*n, 0)
    fullE = sp.eye(1 << (n-1))-O["R"]
    ranks = [[s for s in range(fullE.rows) if s.bit_count() == r]
             for r in range(n)]
    E, H = [sp.Matrix(n, n, lambda i,j:
                       sum(full[s,ranks[j][0]] for s in ranks[i]))
            for full in [fullE,O["H"]]]
    for full,lumped in [(fullE,E),(O["H"],H)]:
        for j,rank in enumerate(ranks):
            for col in rank:
                assert all(sum(full[s,col] for s in ranks[i]) == lumped[i,j]
                           for i in range(n))
    nu = sp.Matrix([sum(O["nu"][s] for s in rank) for rank in ranks])
    assert H*nu == nu and sum(nu) == 1
    assert all(value >= 0 for value in E)
    assert all(E[i,j] == 0 for i in range(n) for j in range(i))
    return K,E,H,nu


def count_criterion_certificate():
    records = []
    for n,rate in [(3,1),(3,64),(3,sp.Rational(1,64)),(4,1)]:
        K,E,H,nu = count_operators(n,rate)
        L = sp.Matrix(n,n,lambda i,j: 1 if i <= j else 0)
        Dnu = sp.diag(*nu)
        F = Dnu.inv()*E*Dnu
        C = L.inv()*F*L
        assert E*Dnu*L == Dnu*L*C
        assert all(C[i,i] == E[i,i] for i in range(n))
        c_fail = [(i,j,str(C[i,j])) for i in range(n) for j in range(n)
                  if C[i,j] < 0]
        h_fail = [(i,r,str(sum(H[j,i+1]-H[j,i] for j in range(r,n))))
                  for i in range(n-1) for r in range(1,n)
                  if sum(H[j,i+1]-H[j,i] for j in range(r,n)) < 0]
        tp2_fail = [(a,b,c,d,str(E[a,c]*E[b,d]-E[a,d]*E[b,c]))
                    for a,b in combinations(range(n),2)
                    for c,d in combinations(range(n),2)
                    if E[a,c]*E[b,d] < E[a,d]*E[b,c]]
        assert not c_fail and not h_fail
        for threshold in range(1,n):
            tail = sp.Matrix([[int(i >= threshold) for i in range(n)]])
            for f in [tail*(H-sp.eye(n)), (tail*nu)[0]*sp.ones(1,n)-tail*H]:
                a = f*Dnu*L
                assert a[-1] == 0 and all(value >= 0 for value in a)
        record = dict(n=n,rate=str(rate),C_nonnegative=True,H_monotone=True,
                      TP2_failures=tp2_fail,
                      E=[[str(x) for x in E.row(i)] for i in range(n)],
                      H=[[str(x) for x in H.row(i)] for i in range(n)],
                      nu=list(map(str,nu)),
                      C=[[str(x) for x in C.row(i)] for i in range(n)])
        if n in (3,4) and rate == 1:
            starts = [(1,)+stable((1 << r)-1,n-1) for r in range(n)]
            raw,nt = absorption(starts,K,(sp.Rational(1,2),)*n)
            rev,nr = absorption(starts,K,(sp.Rational(1,2),)*n,reverse=True)
            assert raw == rev
            A = sp.zeros(n)
            B = sp.zeros(n)
            for col in range(n):
                for mask in range(1 << n):
                    row = (mask >> 1).bit_count()
                    (B if mask & 1 else A)[row,col] += raw[col,mask]
            rawE = 2*sp.eye(n)-B.inv()
            rawH = A*B.inv()*(sp.eye(n)-rawE).inv()
            assert rawE == E and rawH == H
            rawnu,counts = stationary_law(trace_kernel(K,[0]),(sp.Rational(1,2),)*(n-1))
            assert [sum(rawnu[s] for s in range(1 << (n-1)) if s.bit_count()==i)
                    for i in range(n)] == list(nu)
            record["raw_toppling_order_counts"] = [nt,nr]
            record["raw_trace_stationarity_counts"] = counts
        records.append(record)
        print(n,rate,"C nonnegative; H monotone; TP2 failures",len(tp2_fail),flush=True)

    # Stationarity, count nonincrease, and all excursion algebra identities
    # alone do not imply the terminal-tail inequality.
    E = sp.Matrix([[sp.Rational(1,5),sp.Rational(1,100)],[0,sp.Rational(1,10)]])
    H = sp.Matrix([[sp.Rational(1,4),sp.Rational(3,4)],
                   [sp.Rational(3,4),sp.Rational(1,4)]])
    nu = sp.Matrix([sp.Rational(1,2)]*2)
    J = H*(sp.eye(2)-E)
    assert H*nu == nu and all(v >= 0 for v in J)
    assert sp.ones(1,2)*(E+J) == sp.ones(1,2)
    w = E*nu
    gap = sum(w)*nu[1]-(H*w)[1]
    assert gap == -sp.Rational(11,800)
    obstruction = dict(scope="abstract kernels, not an ARW counterexample",
                       E=str(E),H=str(H),J=str(J),nu=list(map(str,nu)),
                       terminal_tail_gap_at_k1=str(gap))

    # With larger downward jumps, one-step decreasing density is insufficient
    # even when H is monotone and every excursion algebra identity holds.
    E = sp.Matrix([[0,0,sp.Rational(2,5)],
                   [0,sp.Rational(7,20),0],
                   [0,0,sp.Rational(1,100)]])
    H = sp.Matrix([[sp.Rational(1,2),sp.Rational(29,100),sp.Rational(21,100)],
                   [sp.Rational(1,2),sp.Rational(13,50),sp.Rational(6,25)],
                   [0,sp.Rational(9,20),sp.Rational(11,20)]])
    nu = sp.Matrix([sp.Rational(1,3)]*3)
    J = H*(sp.eye(3)-E)
    assert H*nu == nu and all(v >= 0 for v in J)
    assert sp.ones(1,3)*(E+J) == sp.ones(1,3)
    assert all(sum(H[j,i+1]-H[j,i] for j in range(r,3)) >= 0
               for i in range(2) for r in range(1,3))
    h = [v/nu[i] for i,v in enumerate(E*nu)]
    assert h[0] >= h[1] >= h[2]
    tail = sp.Matrix([[0,1,1]])
    gap = (tail*(H-sp.eye(3))*E**2*nu)[0]
    assert gap == -sp.Rational(5591,500000)
    one_step_obstruction = dict(
        scope="abstract kernels, not an ARW counterexample",
        E=str(E),H=str(H),J=str(J),nu=list(map(str,nu)),
        initial_density=list(map(str,h)),
        continuation_gap_at_k2=str(gap),
    )
    certificate = dict(
        scope="specified finite parameters only; structural theorem is separate",
        cases=records, stationarity_only_obstruction=obstruction,
        one_step_density_obstruction=one_step_obstruction,
    )
    return certificate


def count_recursions(m, p):
    """Actual hub/capture laws R_b,U_b; p includes the self-loop conversion.

    R stabilizes with no capture site. U records outcomes that fill the
    initially empty removed site. Its extra capture clock has rate one.
    The affine correction removes the forbidden initial self-loop at v.
    """
    p = sp.sympify(p)
    q = 1-p
    R = [sp.zeros(1, m+1)]
    R[0][0, 0] = 1
    U = [sp.zeros(1, m+1)]
    for b in range(1, m+2):
        size = min(b, m+1)
        killed = sp.zeros(size)
        source_r = sp.zeros(size, m+1)
        source_u = sp.zeros(size, m+1)
        for y in range(size):
            birth, death = p*(m-y), q*y
            killed[y, y] = 1+birth+death
            if y+1 < size:
                killed[y, y+1] = -birth
            if y:
                killed[y, y-1] = -death
            source_r[y, :] = R[b-1][y, :]
            source_u[y, :] = U[b-1][y, :] + R[b-1][y, :]
        if b <= m:
            source_r[b-1, b] += p*(m-b+1)
        rb = killed.inv()*source_r
        ub = (killed+sp.eye(size)).inv()*source_u
        assert killed*rb == source_r
        assert (killed+sp.eye(size))*ub == source_u
        if b <= m:
            boundary = sp.zeros(1, m+1)
            boundary[0, b] = 1
            rb = rb.col_join(boundary)
            ub = ub.col_join(sp.zeros(1, m+1))
        R.append(rb)
        U.append(ub)
    H = sp.Matrix(m+1, m+1, lambda i,j: R[j+1][j, i])
    Ehat = sp.Matrix(m+1, m+1, lambda i,j: U[j+1][j, i])
    E = ((m+2)*Ehat-sp.eye(m+1))/(m+1)
    pi = sp.Matrix([sp.binomial(m,y)*p**y*q**(m-y)
                    for y in range(m+1)])
    nu = R[m].T*pi
    return R, U, H, E, nu


def count_continuant(m, p, b, a):
    """Positive binomial formula for det(a I-L_b), used at a=1,2."""
    return sp.factorial(b)*sum(
        sp.binomial(m+a, ell)*sp.binomial(a+b-ell-1, b-ell)
        *p**ell*(1-p)**(b-ell) for ell in range(b+1))


def check_count_formulas(m, p, H, E, nu):
    """Check the positive column recurrence and exact diagonal formula."""
    assert H*nu == nu and sum(nu) == 1
    assert sp.ones(1,m+1)*H == sp.ones(1,m+1)
    assert all(value >= 0 for value in H)
    assert all(value >= 0 for value in E)
    assert all(sum(E[:,j]) <= 1 for j in range(m+1))
    assert all(E[i,j] == 0 for i in range(m+1) for j in range(i))
    determinants = {}
    for a in (1,2):
        values = [count_continuant(m,p,b,a) for b in range(m+2)]
        assert values[0] == 1 and values[1] == a+p*m
        for b in range(1,m+1):
            assert values[b+1] == (
                (a+p*(m-b)+(1-p)*b)*values[b]
                -p*(m-b+1)*(1-p)*b*values[b-1])
        assert all(value > 0 for value in values)
        determinants[a] = values
    identity = sp.eye(m+1)
    for j in range(m+1):
        birth, death = p*(m-j), (1-p)*j
        rhs = identity[:,j]
        if j < m:
            rhs += birth*identity[:,j+1]
        if j:
            rhs += death*H*H[:,j-1]
        assert (1+birth+death)*H[:,j] == rhs
        pivot = 1+birth+death-(death*H[j,j-1] if j else 0)
        assert pivot == determinants[1][j+1]/determinants[1][j]
        diagonal = 1-sp.Rational(m+2,m+1)*(
            determinants[1][j+1]*determinants[2][j]
            /(determinants[1][j]*determinants[2][j+1]))
        assert E[j,j] == diagonal


def allsize_structure_certificate(saved_cases):
    """Only the new structural diagnostics; do not rerun the default suite."""
    expected = {(3,sp.S.One), (3,sp.Integer(64)),
                (3,sp.Rational(1,64)), (4,sp.S.One)}
    baseline = {(case["n"],sp.Rational(case["rate"])): case
                for case in saved_cases}
    assert expected <= baseline.keys()
    records = []
    for n,rate in sorted(expected):
        m = n-1
        p = n*rate/(n+1+n*rate)
        _, _, H, E, nu = count_recursions(m,p)
        case = baseline[n,rate]
        assert H == sp.Matrix(case["H"])
        assert E == sp.Matrix(case["E"])
        assert nu == sp.Matrix(case["nu"])
        check_count_formulas(m,p,H,E,nu)
        records.append(dict(
            n=n, rate=str(rate), p=str(p),
            E_H_nu_match_saved_criterion=True,
            positive_column_recurrence_verified=True,
            continuant_recurrence_verified=True,
            diagonal_determinant_formula_verified=True))
        print(n,rate,"count recurrences and diagonal formula verified",flush=True)

    # The zero-sleep limit predicts this positive-rate failure of count
    # monotonicity. Verify the actual traced ARW, not just its recurrence.
    n, rate = 4, sp.Rational(1,64)
    m = n-1
    p = n*rate/(n+1+n*rate)
    assert p == sp.Rational(1,81)
    _, _, H, E, nu = count_recursions(m,p)
    check_count_formulas(m,p,H,E,nu)
    original = graph_kernel(tuple(tuple(j for j in range(n+1) if j != i)
                                  for i in range(n)))
    traced = trace_kernel(original,[0])
    starts = [add(stable((1 << j)-1,m),v)
              for j in range(m+1) for v in range(m)]
    sleep = (rate/(1+rate),)*m
    raw, count = absorption(starts,traced,sleep)
    reverse, reverse_count = absorption(starts,traced,sleep,reverse=True)
    assert raw == reverse
    rawH = sp.eye(m+1)/n
    for j in range(m+1):
        for v in range(m):
            for mask in range(1 << m):
                rawH[mask.bit_count(),j] += raw[j*m+v,mask]/n
    assert rawH == H
    rawnu, stationary_counts = stationary_law(traced,sleep)
    assert nu == sp.Matrix([
        sum(rawnu[mask] for mask in range(1 << m) if mask.bit_count() == i)
        for i in range(m+1)])
    gap = sum(H[i,3]-H[i,2] for i in range(1,m+1))
    assert gap == -sp.Rational(26233100,452394123)
    L = sp.Matrix(m+1,m+1,lambda i,j: 1 if i <= j else 0)
    W = sp.diag(*nu)
    C = L.inv()*W.inv()*E*W*L
    pqd = [sum(nu[j]*sum(H[i,j] for i in range(t,m+1))
               for j in range(r+1,m+1))
           -sum(nu[j] for j in range(r+1,m+1))
            *sum(nu[i] for i in range(t,m+1))
           for r in range(m) for t in range(1,m+1)]
    assert all(value >= 0 for value in C)
    assert all(value >= 0 for value in pqd)
    positive_rate = dict(
        n=n, rate=str(rate), p=str(p),
        H=[[str(value) for value in H.row(i)] for i in range(m+1)],
        nu=list(map(str,nu)), threshold=1, input_counts=[2,3],
        tail_difference=str(gap),
        raw_toppling_order_counts=[count,reverse_count],
        raw_trace_stationarity_counts=stationary_counts,
        positive_column_recurrence_verified=True,
        diagonal_determinant_formula_verified=True,
        C_nonnegative=True, stationary_PQD_nonnegative=True)
    print(n,rate,"raw H monotonicity obstruction:",gap,flush=True)

    # Minimal saved forward-tail LC witness. At p=0 the independent
    # first-event recursion has no birth term and does not use resolvents.
    m, p = 7, sp.S.Zero
    _, _, H, E, nu = count_recursions(m,p)
    check_count_formulas(m,p,H,E,nu)
    pure = sp.zeros(m+1)
    pure[0,0] = 1
    for j in range(1,m+1):
        pure[:,j] = sp.Rational(j,j+1)*pure*pure[:,j-1]
        pure[j,j] += sp.Rational(1,j+1)
    assert pure == H
    h = [1-H[0,j] for j in range(m+1)]
    gap = h[6]**2-h[5]*h[7]
    assert gap == -sp.Rational(11737717,32112640000)
    zero_rate = dict(
        outside_sites=m, p="0", threshold=1, input_counts=[5,6,7],
        tail_probabilities=list(map(str,h[5:8])), LC_gap=str(gap),
        independent_zero_rate_recurrence_verified=True,
        scope="strict gap and continuity imply failure at small positive p")
    print("m=7,p=0 forward-tail LC obstruction:",gap,flush=True)
    return dict(
        scope="finite exact formula checks and specified obstructions; "
              "all-rate density-cone positivity and stationary PQD remain open",
        recurrence_checks=records,
        positive_rate_H_monotonicity_obstruction=positive_rate,
        zero_rate_forward_tail_LC_obstruction=zero_rate)


def continuant_polynomials(m, p, a):
    """Continuants polynomial in the discount a, so derivatives are exact."""
    values = [sp.S.One, a+p*m]
    for i in range(1,m+1):
        values.append(sp.expand(
            (a+p*(m-i)+(1-p)*i)*values[-1]
            -p*(m-i+1)*(1-p)*i*values[-2]))
    return values


def discounted_hub_kernel(m, p, a):
    """Positive column recurrence for E[exp(-(a-1)T); terminal count]."""
    values = continuant_polynomials(m,p,a)
    kernel = sp.zeros(m+1)
    identity = sp.eye(m+1)
    for j in range(m+1):
        column = identity[:,j]
        if j < m:
            column += p*(m-j)*identity[:,j+1]
        if j:
            column += (1-p)*j*sum(
                (kernel[k,j-1]*kernel[:,k] for k in range(j)),
                sp.zeros(m+1,1))
        kernel[:,j] = column/(values[j+1]/values[j])
    return kernel


def continuant_operator_check(m, p):
    """Compare capture/continuant formulas with the independent R_b,U_b solves."""
    _, _, H, E, nu = count_recursions(m,p)
    H2 = discounted_hub_kernel(m,p,sp.Integer(2))
    assert discounted_hub_kernel(m,p,sp.S.One) == H
    identity = sp.eye(m+1)
    Ehat = (identity+(m+1)*E)/(m+2)
    assert H2 == H*(identity-Ehat)
    death = sp.zeros(m+1)
    for j in range(1,m+1):
        death[j-1,j] = (1-p)*j
    A1 = sp.diag(*[1+p*(m-j)+(1-p)*j for j in range(m+1)])
    V1, V2 = A1-H*death, A1+identity-H2*death
    assert Ehat == identity-V1*V2.inv()

    a = sp.Symbol("a")
    values = continuant_polynomials(m,p,a)
    kappa = [sp.cancel(values[j+1]/values[j]) for j in range(m+1)]
    k1 = [value.subs(a,1) for value in kappa]
    k2 = [value.subs(a,2) for value in kappa]
    ell1 = [sp.cancel(sp.diff(value,a)/value).subs(a,1) for value in kappa]
    ell2 = [sp.cancel(sp.diff(value,a)/value).subs(a,2) for value in kappa]
    for i in range(m+1):
        assert H[i,i] == ell1[i] and H2[i,i] == ell2[i]

    # Positive capture-column recurrence; no U_b values enter this calculation.
    captured = sp.zeros(m+1)
    for j in range(m+1):
        column = identity[:,j]
        if j:
            column += (1-p)*j*(H*captured[:,j-1]+sum(
                (H2[k,j-1]*captured[:,k] for k in range(j)),
                sp.zeros(m+1,1)))
        captured[:,j] = column/k2[j]
    assert captured == Ehat

    L = sp.Matrix(m+1,m+1,lambda i,j: 1 if i <= j else 0)
    W = sp.diag(*nu)
    C = L.inv()*W.inv()*E*W*L
    alpha = sp.Rational(m+2,m+1)
    for i in range(m):
        entry = (1-p)*(i+1)/k2[i+1]*(ell1[i]-k1[i]/k2[i]*ell2[i])
        assert Ehat[i,i+1] == entry
        formula = alpha*(k1[i+1]/k2[i+1]-k1[i]/k2[i]
                         +entry*nu[i+1]/nu[i])
        assert C[i,i+1] == formula
    assert nu[m]/nu[m-1] == p/(k1[m-1]*(1-ell1[m]))
    top_kappa = alpha*(
        k1[m]*(k1[m]-sp.diff(kappa[m],a).subs(a,2))
        /(k2[m]*(k1[m]-sp.diff(kappa[m],a).subs(a,1)))
        -k1[m-1]/k2[m-1])
    assert C[m-1,m] == top_kappa

    q = 1-p
    g = lambda b: sum(sp.binomial(m,k)*p**(m-k)*q**k/sp.Integer(b+k)
                      for k in range(m+1))
    gp = lambda b: -sum(sp.binomial(m,k)*p**(m-k)*q**k/sp.Integer(b+k)**2
                        for k in range(m+1))
    h = lambda b: sum(sp.binomial(m-1,k)*p**(m-1-k)*q**k
                      /sp.Integer((b+k)*(b+k+1)) for k in range(m))
    for b in (1,2):
        assert kappa[m].subs(a,b) == 1/g(b)
        assert kappa[m-1].subs(a,b) == g(b)/h(b)
    numerator = h(1)*(g(2)**2+g(1)*gp(2))-g(1)*h(2)*(g(1)+gp(1))
    denominator = g(2)*h(1)*(g(1)+gp(1))
    assert denominator > 0
    assert alpha*numerator/denominator == C[m-1,m] > 0

    # One symbolic generating function verifies every return count for the
    # left comparison on the outside-full event, with no truncation in k.
    z = sp.Symbol("z")
    d, e, c = E[m-1,m-1], E[m,m], C[m-1,m]
    flux = H[m,m-1]*nu[m-1]
    assert flux == (1-H[m,m])*nu[m] > 0
    assert 0 <= d < 1 and 0 <= e < 1 and c > 0
    resolvent = (identity-z*E).inv()*nu
    left_series = ((H-identity)*resolvent)[m]
    expected = z*flux*c/((1-z*d)*(1-z*e))
    assert sp.cancel(left_series-expected) == 0
    occupied = (1-z)*resolvent
    outside_marginal = nu-H*occupied+occupied
    assert sp.cancel(nu[m]-outside_marginal[m]-(1-z)*expected) == 0
    return dict(
        outside_sites=m, p=str(p), discount_factorization_verified=True,
        positive_capture_recurrence_verified=True,
        logarithmic_diagonal_identity_verified=True,
        first_superdiagonal_formulas_verified=True,
        top_stationary_ratio_verified=True,
        top_continuant_and_integral_formulas_verified=True,
        top_C_entry=str(C[m-1,m]),
        all_k_full_event_left_resolvent_verified=True,
        outside_full_marginal_trace_resolvent_verified=True,
        full_event_resolvent_parameters=dict(
            d=str(d), e=str(e), A=str(flux), c=str(c)))


def top_gap_polynomial_check(m):
    """One symbolic polynomial case, with all factorization residuals saved."""
    t = sp.Symbol("t")
    Gfun = lambda a: sp.Poly(sum(sp.binomial(m,k)*t**k/sp.Integer(a+k)
                                for k in range(m+1)),t)
    Jfun = lambda a: sp.Poly(sum(sp.binomial(m,k)*t**k/sp.Integer(a+k)**2
                                for k in range(m+1)),t)
    Ffun = lambda a: sp.Poly(sum(sp.binomial(m-1,k)*t**k
                                /sp.Integer((a+k)*(a+k+1)) for k in range(m)),t)
    G, B, J, J2 = Gfun(1), Gfun(2), Jfun(1), Jfun(2)
    F, F2 = Ffun(1), Ffun(2)
    U = G-J
    S = J-G+B
    H = G*S-B
    D = sp.Poly(1,t)+sp.Poly(m*t-1,t)*G
    P = F*(B**2-G*J2)-F2*G*U
    residuals = {
        "mt_F": m*sp.Poly(t,t)*F-(G-1),
        "mt_F2": m*sp.Poly(t,t)*F2-(2*B-G),
        "J2_identity": (m+2)*sp.Poly(t,t)*J2-(U+sp.Poly(t,t)*B),
        "B_identity": (m+2)*sp.Poly(t,t)*B
                      -(1+sp.Poly((m+1)*t-1,t)*G),
        "positive_D_representation": D-m*sp.Poly(t,t)*(F+(m+2)*sp.Poly(t,t)*F2),
        "P_factorization": m*(m+2)*sp.Poly(t**2,t)*P-D*H,
    }
    assert all(value.is_zero for value in residuals.values())
    assert all(S.nth(k) == sp.binomial(m,k)/((k+1)**2*(k+2))
               for k in range(m+1))
    assert H.nth(0) == H.nth(1) == 0
    assert all(H.nth(k) > 0 for k in range(2,2*m+1))
    assert D.nth(0) == 0 and all(D.nth(k) > 0 for k in range(1,m+2))
    assert P.nth(0) == 0 and all(P.nth(k) > 0 for k in range(1,3*m))
    if m == 1:
        assert H == sp.Poly(t**2/24,t)
        assert D == sp.Poly(t*(1+t)/2,t)
        assert P == sp.Poly(t*(1+t)/144,t)
    return dict(
        outside_sites=m, residuals={name:str(value.as_expr())
                                   for name,value in residuals.items()},
        H_coefficients=[str(H.nth(k)) for k in range(2*m+1)],
        D_coefficients=[str(D.nth(k)) for k in range(m+2)],
        P_coefficients=[str(P.nth(k)) for k in range(3*m)],
        strict_positive_P_degrees=[1,3*m-1])


def continuant_positivity_certificate():
    """Audit residuals for the all-size analytic top-entry positivity proof.

    These checks accompany the proof, not a size/rate sweep. For 2<=ell<=m,
    retain G-degree ell with S-degree 0 and G-degree ell-1 with S-degree 1,
    where S=J-G+B. Their normalized sum exceeds B's coefficient because
    m/(m-ell+1)>1 and the displayed residual is nonnegative. For ell>m,
    B has no coefficient and the convolution is strictly positive.
    Thus H has positive coefficients in degrees 2..2m. The positive
    representation of D then proves positivity of P in degrees 1..3m-1.
    """
    m,t,G,B,J = sp.symbols("m t G B J")
    F = (G-1)/(m*t)
    F2 = (2*B-G)/(m*t)
    J2 = (G-J+t*B)/((m+2)*t)
    D = 1+(m*t-1)*G
    H = (G-1)*B-G*(G-J)
    P = F*(B**2-G*J2)-F2*G*(G-J)
    replacement = (1+((m+1)*t-1)*G)/((m+2)*t)
    residuals = {
        "general_P_factorization": sp.factor(
            (m*(m+2)*t**2*P-D*H).subs(B,replacement)),
        "general_positive_D_representation": sp.factor(
            (D-m*t*(F+(m+2)*t*F2)).subs(B,replacement)),
    }
    k,ell = sp.symbols("k ell")
    residuals["small_factor_coefficient"] = sp.factor(
        1/(k+1)**2-1/(k+1)+1/(k+2)-1/((k+1)**2*(k+2)))
    residuals["two_term_lower_bound"] = sp.factor(
        1/(2*(ell+1))+sp.Rational(1,12)-1/(ell+2)
        -(ell-1)*(ell-2)/(12*(ell+1)*(ell+2)))
    residuals["strict_excess"] = sp.factor(
        m/(m-ell+1)-1-(ell-1)/(m-ell+1))

    # The top-two-state identity is formal, including equal diagonal values.
    d,e,b,z = sp.symbols("d e b z")
    pair = sp.Matrix([[d,b],[0,e]])
    pair_series = (sp.Matrix([[1,-1]])*(sp.eye(2)-z*pair).inv()
                   *sp.ones(2,1))[0]
    residuals["general_all_k_left_resolvent"] = sp.cancel(
        pair_series-z*(d+b-e)/((1-z*d)*(1-z*e)))
    assert all(value == 0 for value in residuals.values())

    # Preserve the two targeted operator checks from the scratch derivation.
    operator_checks = [continuant_operator_check(2,sp.Rational(3,7)),
                       continuant_operator_check(3,sp.Rational(4,9))]
    # Explicit boundary m=1, and the single m=4 symbolic residual used earlier.
    edge = top_gap_polynomial_check(1)
    symbolic = top_gap_polynomial_check(4)
    print("Capture/top identities and general factorization residuals verified; "
          "m=1 edge and m=4 symbolic polynomial verified.",flush=True)
    return dict(
        scope="targeted exact checks accompanying the all-size analytic proof; "
              "the full C>=0 criterion and stationary PQD remain open",
        audited_conclusion="C[m-1,m]>0 for every m>=1 and 0<p<1",
        audited_full_event_consequence=dict(
            assumptions="m>=1, 0<p<1, and finite positive inserted rate; "
                        "s is its original sleep probability, q=1-s",
            all_k_left_gap="A*c*sum_{ell=0}^{k-1} d^(k-1-ell)*e^ell; zero at k=0",
            marginal_trace_gap="s*q*A*c/((1-q*d)*(1-q*e)) > 0",
            scope="outside-full marginal event only; no conditional/full NA "
                  "or right-half excursion comparison is claimed"),
        operator_checks=operator_checks,
        analytic_proof_residuals={name:str(value) for name,value in residuals.items()},
        coefficient_argument=dict(
            assumptions="m>=1 integer; for the lower bound, 2<=ell<=m",
            retained_terms="G_degree_ell*S_degree_0 and G_degree_(ell-1)*S_degree_1",
            normalized_lower_bound="1/(2*(ell+1))+m/(12*(m-ell+1))",
            normalized_B_coefficient="1/(ell+2)",
            nonnegative_base_difference="(ell-1)*(ell-2)/(12*(ell+1)*(ell+2))",
            strictly_positive_excess="(ell-1)/(12*(m-ell+1))",
            H_positive_degrees="2,...,2m; degrees 0 and 1 vanish",
            D_positive_degrees="1,...,m+1; degree 0 vanishes",
            P_positive_degrees="1,...,3m-1; degree 0 vanishes"),
        m1_edge=edge, symbolic_m4=symbolic)


def interior_pqd_polynomials():
    """Check the finite-degree algebra underlying the all-size r=1 proof.

    m remains an indeterminate. Only degrees affected by subtraction need
    explicit formulas; all higher coefficients are those of a product of
    positive polynomials and (1+z)^m. No bounded-size extrapolation is used.
    """
    m,z,power = sp.symbols("m z power", positive=True)
    q,p = 1/(1+z),z/(1+z)
    a,b = m*p,(m-1)*p
    d = (1+a)*(1+b)+q
    dc = sp.Poly(2+2*(m+1)*z+m*(m+1)*z**2,z)
    ac,bc = sp.Poly(1+(m+1)*z,z),sp.Poly(1+m*z,z)
    ec = sp.Poly((1+z)**2+m*z*dc.as_expr(),z)
    kc = sp.Poly(1+(2*m+1)*z+m*(m+1)*z**2,z)
    fc = sp.Poly(2+4*m*z+m*(m+1)*z**2,z)
    numerator0 = ac.as_expr()*dc.as_expr()*power-2*bc.as_expr()*ec.as_expr()
    numerator1 = dc.as_expr()*kc.as_expr()*power-fc.as_expr()*ec.as_expr()
    upper0 = 2*(1+b)/(power*(1+a)*d)
    upper1 = (1+z)*(d+2*b*q)/(power*d)
    bound0,bound1 = q/(q+a*d),(q+a*(2+b))/(q+a*d)
    residuals = {
        "R0_minus_U0": sp.cancel(bound0-upper0
            -(1+z)**2*numerator0/(power*ec.as_expr()*ac.as_expr()*dc.as_expr())),
        "R1_minus_U1": sp.cancel(bound1-upper1
            -(1+z)*numerator1/(power*ec.as_expr()*dc.as_expr())),
    }

    def choose(k):
        return sp.prod(m-j for j in range(k))/sp.factorial(k)

    def low_coefficients(left,right,degree):
        return [sp.factor(sum(left.nth(j)*choose(k-j)
                              for j in range(min(k,left.degree())+1))-right.nth(k))
                for k in range(degree+1)]

    coefficients0 = low_coefficients(ac*dc,2*bc*ec,4)
    coefficients1 = low_coefficients(dc*kc,fc*ec,5)
    expected0 = [0,0,0,m*(m-1)*(m+1)/3,
                 m*(m-1)*(15*m**2+5*m+2)/12]
    expected1 = [0,0,m*(m-1),7*m*(m-1)*(m+1)/3,
                 m*(m-1)*(43*m**2+33*m+14)/12,
                 m*(m-1)*(103*m**3+43*m**2-2*m-2)/30]
    assert all(sp.expand(a-b)==0 for a,b in zip(coefficients0,expected0))
    assert all(sp.expand(a-b)==0 for a,b in zip(coefficients1,expected1))
    assert (2*bc*ec).degree()==4 and (fc*ec).degree()==5
    shift = sp.Symbol("h",nonnegative=True)
    shifted = sp.expand((103*m**3+43*m**2-2*m-2).subs(m,shift+2))
    assert shifted==103*shift**3+661*shift**2+1406*shift+990

    # Formal operator identity, without assuming that H,A,D commute.
    hs,aa,dd = sp.symbols("H A D",commutative=False)
    residuals["stationary_factorization"] = sp.expand(
        (hs-1)*(aa-dd-hs*dd)-(hs*aa-hs**2*dd-aa+dd))
    assert all(value==0 for value in residuals.values())
    return dict(
        assumptions="integer m>=2 and z=p/(1-p)>0; power denotes (1+z)^m",
        residuals={name:str(value) for name,value in residuals.items()},
        polynomials={"Dcal":str(dc.as_expr()),"Acal":str(ac.as_expr()),
                     "Bcal":str(bc.as_expr()),"Ecal":str(ec.as_expr()),
                     "Kcal":str(kc.as_expr()),"Fcal":str(fc.as_expr())},
        definitions={"P0":"Acal*Dcal*(1+z)^m-2*Bcal*Ecal",
                     "P1":"Dcal*Kcal*(1+z)^m-Fcal*Ecal"},
        P0_coefficients_degrees_0_to_4=list(map(str,coefficients0)),
        P1_coefficients_degrees_0_to_5=list(map(str,coefficients1)),
        last_coefficient_positive_shift=str(shifted),
        higher_coefficients_argument="Subtracted degrees are 4 and 5; higher coefficients come only from positive products. Falling factorials include m=2,3.")


def reverse_continuation_diagnostic(case):
    """One preserved finite case; this is not an all-size order claim."""
    H,nu = sp.Matrix(case["H"]),sp.Matrix(case["nu"])
    size,rate = H.rows,sp.Rational(case["rate"])
    m = size-1
    p = size*rate/(size+1+size*rate)
    q = 1-p
    beta = [p*(m-j) for j in range(size)]
    delta = [q*j for j in range(size)]
    pi = sp.Matrix([sp.binomial(m,j)*p**j*q**(m-j) for j in range(size)])
    identity = sp.eye(size)
    birth,death = sp.zeros(size),sp.zeros(size)
    for j in range(size):
        if j<m: birth[j+1,j]=beta[j]
        if j: death[j-1,j]=delta[j]
    diagonal = sp.diag(*(1+beta[j]+delta[j] for j in range(size)))
    generator = birth+death-sp.diag(*(beta[j]+delta[j] for j in range(size)))
    F = diagonal-death-H*death
    assert H*diagonal==identity+birth+H**2*death
    assert (H-identity)*F==generator and generator*pi==sp.zeros(size,1)
    assert death*pi==sp.diag(*beta)*pi
    assert sum(F*pi)==1 and F*pi==nu and H*nu==nu
    assert nu==(identity+sp.diag(*delta))*pi-H*sp.diag(*beta)*pi

    a,b = m*p,(m-1)*p
    d = (1+a)*(1+b)+q
    u,v,w = 1/(1+a),q/((1+a)*d),b*(1+a)/d
    assert (H[0,0],H[0,1],H[2,1])==(u,v,w)
    x,y,s = nu[0],nu[1],nu[0]+nu[1]
    upper0 = 2*q**m*(1+b)/((1+a)*d)
    upper1 = q**(m-1)*(d+2*b*q)/d
    assert x<=upper0<q/(q+a*d)
    assert s<=upper1<(q+a*(2+b))/(q+a*d)
    r1_gaps = {1:x*u+y*v-x*s,2:s*(1-s)-y*w}
    for threshold in range(1,size):
        direct = (sum(nu[j]*sum(H[i,j] for i in range(threshold,size))
                      for j in range(2,size))-sum(nu[2:,0])*sum(nu[threshold:,0]))
        expected = r1_gaps.get(threshold,s*sum(nu[threshold:,0]))
        assert direct==expected and direct>0
        r1_gaps[threshold]=direct

    W = sp.diag(*nu)
    prefix = sp.Matrix(size,size,lambda i,r:int(bool(i<=r)))
    reverse = W.inv()*H*W
    V = diagonal-H*death
    assert all(V[i,j]==0 for i in range(size) for j in range(i))
    assert all(V[i,i]>0 for i in range(size))
    assert all(V[i,j]<=0 for i in range(size) for j in range(i+1,size))
    h = V.inv()*nu
    assert all(value>0 for value in h) and (identity+birth)*h==nu
    R = W.inv()*(identity+birth)*sp.diag(*h)
    Q = sp.diag(*h).inv()*V.inv()*W
    assert reverse==R*Q
    for kernel in (reverse,R,Q):
        assert kernel*sp.ones(size,1)==sp.ones(size,1)
        assert all(value>=0 for value in kernel)
    assert all((prefix.inv()*R*prefix)[i,j]>=0
               for i in range(size) for j in range(size))
    stop = []
    for i in range(size):
        s_i = nu[i]/(V[i,i]*h[i])
        coefficients = [delta[j]*H[i,j-1]*h[j]/(V[i,i]*h[i])
                        for j in range(i+1,size)]
        assert s_i>0 and all(value>=0 for value in coefficients)
        assert s_i+sum(coefficients)==1
        expected = s_i*identity[i,:]+sum(
            (coefficients[j-i-1]*Q[j,:] for j in range(i+1,size)),sp.zeros(1,size))
        assert expected==Q[i,:]
        stop.append(s_i)
    diagnostics = {}
    for name,kernel in (("reverse_P",reverse),("upward_Q",Q)):
        criterion = prefix.inv()*kernel*prefix
        negative = [(i,r,str(criterion[i,r])) for i in range(size)
                    for r in range(size) if criterion[i,r]<0]
        assert not negative
        diagnostics[name] = dict(monotone_at_this_parameter=True,
            criterion=[[str(value) for value in criterion.row(i)] for i in range(size)])
    return dict(m=m,p=str(p),rate=str(rate),
                stationary_binomial_identity_verified=True,
                stationary_signed_vector_normalization_verified=True,
                r1_gaps={str(t):str(value) for t,value in r1_gaps.items()},
                reverse_factorization_verified=True,
                positive_upward_row_recurrence_verified=True,
                upward_stopping_probabilities=list(map(str,stop)),
                diagnostics=diagnostics)


def interior_pqd_certificate(saved):
    """Symbolic all-size r=1 audit, plus exactly five saved finite checks."""
    polynomials = interior_pqd_polynomials()
    cases = list(saved["excursion_count_criterion"]["cases"])
    cases.append(saved["excursion_allsize_structure"]["positive_rate_H_monotonicity_obstruction"])
    indexed = {}
    for case in cases:
        size,rate = len(case["H"]),sp.Rational(case["rate"])
        indexed[(size-1,size*rate/(size+1+size*rate))]=case
    requested = [(2,sp.Rational(3,7)),(2,sp.Rational(48,49)),
                 (2,sp.Rational(3,259)),(3,sp.Rational(4,9)),(3,sp.Rational(1,81))]
    assert all(key in indexed for key in requested)
    records = [reverse_continuation_diagnostic(indexed[key]) for key in requested]
    print("General-m r=1 polynomial identities and five preserved reverse P/Q diagnostics verified.",flush=True)
    return dict(
        audited_conclusion="stationary PQD cut r=1 is strict for every m>=2, 0<p<1, and 1<=t<=m",
        proof_scope="analytic coefficient proof with m indeterminate; finite checks are not used to infer all-size positivity",
        stationary_identity="F=A-D-HD; (H-I)F=L; F*pi has mass one, hence F*pi=nu; D*pi=diag(beta)*pi",
        reverse_factorization="P=R Q, with R supported on {i-1,i} and Q an upward stopping kernel",
        reverse_order_scope="P and Q monotonicity checked only at the five preserved parameters; all-size monotonicity and the overshoot bound remain open",
        analytic_polynomial_checks=polynomials,saved_case_diagnostics=records)


def reverse_monotonicity_certificate():
    """Prove the inverse criterion at m<=3 with the rate indeterminate.

    T=Q^-1 has row sums one and is upper triangular. If L[i,r]=1[i<=r],
    the only nontrivial off-diagonals of L^-1 T L at m=3 are (0,1),
    (0,2),(1,2). Their negatives have positive rational factorizations.
    This fixed-size proof makes no assertion about m>=4 or the C matrix.
    """
    z = sp.Symbol("z",positive=True)
    p,q = z/(1+z),1/(1+z)

    # At m=2 the sole sign follows from a 2x2 determinant cancellation.
    a,b = 2*p,p
    k,d = 1+a,(1+a)*(1+b)+q
    u,w = 1/k,a/k
    v,t = q/(k*d),(k**2+q*a)/(k*d)
    assert sp.cancel(u*t-w*v-1/d)==0
    pi0,pi1 = q**2,a*q
    nu0 = pi0-a*pi0*u-b*pi1*v
    nu1 = (1+q)*pi1-a*pi0*w-b*pi1*t
    cross = sp.cancel(nu1*v-nu0*t)
    assert sp.cancel(cross+q**2/(k*d))==0

    # The first size at which the inverse criterion can be stronger than
    # Q monotonicity. No selected-rate substitution or size sweep is used.
    m,size = 3,4
    beta = [p*(m-j) for j in range(size)]
    delta = [q*j for j in range(size)]
    birth,death = sp.zeros(size),sp.zeros(size)
    for j in range(size):
        if j<m: birth[j+1,j]=beta[j]
        if j: death[j-1,j]=delta[j]
    diagonal = sp.diag(*(1+beta[j]+delta[j] for j in range(size)))
    identity,H = sp.eye(size),sp.zeros(size)
    for j in range(size):
        pivot = sp.cancel(diagonal[j,j]-(delta[j]*H[j,j-1] if j else 0))
        column = identity[:,j]+birth[:,j]
        if j:
            column += delta[j]*sum(
                (H[k,j-1]*H[:,k] for k in range(j)),sp.zeros(size,1))
        H[:,j] = column.applyfunc(lambda value:sp.cancel(value/pivot))
    pi = sp.Matrix([sp.binomial(m,j)*p**j*q**(m-j) for j in range(size)])
    nu = ((diagonal-death-H*death)*pi).applyfunc(sp.cancel)
    h = sp.zeros(size,1)
    for i in range(size):
        h[i] = sp.cancel(nu[i]-(beta[i-1]*h[i-1] if i else 0))
    V = diagonal-H*death
    assert all(sp.cancel(value)==0 for value in V*h-nu)
    assert all(sp.cancel(value)==0 for value in H*nu-nu)
    assert sp.cancel(sum(nu)-1)==0
    T = (sp.diag(*nu).inv()*V*sp.diag(*h)).applyfunc(sp.cancel)
    prefix = sp.Matrix(size,size,lambda i,r:int(bool(i<=r)))
    criterion = (prefix.inv()*T*prefix).applyfunc(sp.cancel)
    assert all(sp.cancel(value)==0 for value in T*sp.ones(size,1)-sp.ones(size,1))
    assert all(criterion[i,j]==0 for i in range(size) for j in range(i))
    assert all(criterion[i,m]==0 for i in range(m)) and criterion[m,m]==1
    records = []
    for i,r in [(0,1),(0,2),(1,2)]:
        numerator,denominator = sp.fraction(sp.cancel(-criterion[i,r]))
        npoly,dpoly = sp.Poly(numerator,z),sp.Poly(denominator,z)
        assert not npoly.is_zero and not dpoly.is_zero
        assert all(value>=0 for value in npoly.all_coeffs())
        assert all(value>=0 for value in dpoly.all_coeffs())
        records.append(dict(i=i,r=r,
            negative_offdiagonal=str(sp.factor(-criterion[i,r])),
            numerator_degree=int(npoly.degree()),denominator_degree=int(dpoly.degree()),
            numerator_nonnegative_coefficients=True,denominator_nonnegative_coefficients=True))
    print("All-rate reverse P/Q monotonicity at m<=3 verified by exact inverse-criterion signs.",flush=True)
    return dict(
        audited_conclusion="Q and stationary reverse P are stochastically monotone for every 0<p<1 when m<=3",
        proof_scope="m=0,1 automatic; exact m=2 cancellation and symbolic m=3 polynomial identities; no claim for m>=4 or C>=0",
        inverse_generator="T=diag(nu)^(-1) V diag(h)=I-G; G_ij=delta_j*H[i,j-1]*h_j/nu_i for j>i, G*1=0; Q=(I-G)^(-1)",
        sufficient_criterion="For 0<=i<r<m, sum_{j>r} delta_j*h_j*(H[i,j-1]/nu_i-H[i+1,j-1]/nu[i+1])<=0",
        implication="S=L^(-1) T L is upper triangular with positive diagonal and nonpositive off-diagonal; S^(-1)=L^(-1) Q L>=0; P=R Q with R monotone",
        m2_cross_identity="nu_1 H_01-nu_0 H_11=-q^2/((1+2p)*((1+2p)*(1+p)+q))<0",
        m2_cross_in_z=str(sp.factor(cross)),
        m3_parameter="z=p/(1-p)>0",
        m3_stationary_identities_verified=True,m3_inverse_boundary_column_verified=True,
        m3_nu=[str(sp.factor(value)) for value in nu],
        m3_h=[str(sp.factor(value)) for value in h],
        m3_inverse_offdiagonals=records)


def count_highrate_taylor(m, order):
    """Taylor coefficients in q=1-p of the actual count equations.

    This solves the R_b,U_b equations directly, without the limiting
    triangular basis or stationary vacancy formula. Return E,H,nu jets.
    The forced-launch correction is included in E.
    """
    R = [[sp.zeros(1, m+1) for _ in range(order+1)]]
    R[0][0][0, 0] = 1
    U = [[sp.zeros(1, m+1) for _ in range(order+1)]]
    for b in range(1, m+2):
        size = min(b, m+1)
        D0, D1 = sp.zeros(size), sp.zeros(size)
        for y in range(size):
            D0[y, y], D1[y, y] = 1+m-y, 2*y-m
            if y+1 < size:
                D0[y, y+1], D1[y, y+1] = -(m-y), m-y
            if y:
                D1[y, y-1] = -y
        G1, G2 = D0.inv(), (D0+sp.eye(size)).inv()
        rb, ub = [], []
        for k in range(order+1):
            sr = R[b-1][k].copy()
            su = U[b-1][k]+R[b-1][k]
            if b <= m:
                if k == 0:
                    sr[b-1, b] += m-b+1
                elif k == 1:
                    sr[b-1, b] -= m-b+1
            if k:
                sr -= D1*rb[k-1]
                su -= D1*ub[k-1]
            rb.append(G1*sr)
            ub.append(G2*su)
            assert D0*rb[k] == sr
            assert (D0+sp.eye(size))*ub[k] == su
        if b <= m:
            for k in range(order+1):
                boundary = sp.zeros(1, m+1)
                if k == 0:
                    boundary[0, b] = 1
                rb[k] = rb[k].col_join(boundary)
                ub[k] = ub[k].col_join(sp.zeros(1, m+1))
        R.append(rb)
        U.append(ub)
    E, H = [], []
    for k in range(order+1):
        H.append(sp.Matrix(m+1, m+1, lambda i,j: R[j+1][k][j, i]))
        capture = sp.Matrix(m+1, m+1, lambda i,j: U[j+1][k][j, i])
        E.append(sp.Rational(m+2, m+1)*capture)
        if k == 0:
            E[k] -= sp.eye(m+1)/(m+1)
    pi = [sp.Matrix([
        sp.binomial(m, y)*sp.binomial(y, k-(m-y))*(-1)**(k-(m-y))
        if m-y <= k <= m else sp.S.Zero for y in range(m+1)])
          for k in range(order+1)]
    nu = [sum((R[m][r].T*pi[k-r] for r in range(k+1)), sp.zeros(m+1, 1))
          for k in range(order+1)]
    return E, H, nu


def highrate_basis(size, a):
    """Universal vacancy basis A_a; indices are vacancy counts."""
    a = sp.sympify(a)
    return sp.Matrix(size, size, lambda ell,h:
                     sp.binomial(ell, h)*(-1/(a+ell))**(ell-h)
                     if ell >= h else 0)


def highrate_hub_limit(size, a):
    """A_a^-1 diag(1/(a+ell)) A_a, independent of ambient size."""
    a = sp.sympify(a)
    basis = highrate_basis(size, a)
    return basis.inv()*sp.diag(*[1/(a+ell) for ell in range(size)])*basis


def highrate_queue_tail_certificate(H, G, K, c):
    """Enumerate every capture mark in the existing m=3 queue paths.

    The ordered-time/FKG assertion has an analytic proof. These finite
    residuals check its word formula and normalization, not the open
    covering or arrival-count monotonicity hypotheses.
    """
    assert H.shape == G.shape == K.shape == (4, 4) and c.shape == (4, 1)

    def words(d, arrivals=0, services=0, prefix=()):
        if arrivals == d and services == d+1:
            yield prefix
            return
        if 1+arrivals-services == 0:
            return
        if arrivals < d:
            yield from words(d, arrivals+1, services, prefix+('A',))
        if services < d+1:
            yield from words(d, arrivals, services+1, prefix+('S',))

    cache = {}

    def coefficients(d, x):
        if (d, x) not in cache:
            F, M = sp.S.Zero, sp.S.Zero
            for word in words(d):
                arrivals, rates = 0, []
                for event in word:
                    rates.append(sp.Integer(x+arrivals))
                    arrivals += event == 'A'
                F += sp.prod(1/rate for rate in rates)
                for mark, event in enumerate(word):
                    if event == 'S':
                        M += sp.prod(1/(rate+int(index <= mark))
                                     for index, rate in enumerate(rates))
            cache[d, x] = F, M
        return cache[d, x]

    a = [c[ell]/sp.factorial(ell) for ell in range(4)]
    records = []
    for ell in range(4):
        Fmix, Mmix = [], []
        for h in range(ell+1):
            d, x = ell-h, h+1
            F, M = coefficients(d, x)
            Fshift, _ = coefficients(d, x+1)
            factor = sp.factorial(ell)/sp.factorial(h)
            assert factor*F == H[ell, h]
            assert factor*M == K[ell, h]
            assert M == x*F-(x+d)*Fshift
            discount = Fshift/F
            lower = sp.Rational(x*(x+1), (x+d+1)**2)
            upper = sp.Rational(x, x+d)*sp.Rational(x+d, x+d+1)**(d+1)
            assert lower <= discount <= upper
            Fmix.append(a[h]*F/a[ell])
            Mmix.append(a[h]*M/a[ell])
        assert sum(Fmix) == 1
        for r in range(ell+1):
            value = sum(Mmix[r:])
            assert value == sum(K[ell, h]*c[h]/c[ell] for h in range(r, ell+1))
            records.append(dict(vacancy_count=ell, cutoff=r, capture_tail=str(value)))

    g2 = sum(G[2, h]*c[h]/c[2] for h in range(1, 3))
    g3 = sum(G[3, h]*c[h]/c[3] for h in range(1, 4))
    assert g3-g2 == -sp.Rational(1199, 68000)
    return dict(
        scope="same m=3 limit only; no new size/rate case and no proof of F10",
        all_unmarked_word_coefficients_match=True,
        all_capture_positions_match=True,
        stationary_mixture_normalization=True,
        marked_sum_cancellation_identity=True,
        capture_tails=records,
        discounted_tail_shortcut_obstruction=dict(
            old=str(g2), new=str(g3), difference=str(g3-g2)),
        queue_formula="H[l,h]=(l!/h!) F_d(x), K[l,h]=(l!/h!) M_d(x), d=l-h, x=h+1; F sums all admissible queue words and M sums every service capture mark",
        stationary_tail_formula="a_l=c_l/l!; P_l(h)=a_h F_(l-h)(h+1)/a_l; T[l,r]=E_(h~P_l)[1_(h>=r) R_(l-h)(h+1)], R_d=M_d/F_d",
        ordered_time_formula=dict(
            domain="0<u_1<...<u_d<t, 0<v_1<...<v_d<t, u_j<v_j",
            density="exp(sum_j u_j-(x+d)t), normalized by F_d(x)",
            capture_reward="exp(-t)+sum_j exp(-v_j)",
            analytic_FKG_result="For fixed integer d>=0, R_d(x) is strictly increasing for real x>0: the ordered cone is a sublattice, its density is log-modular, and increasing x tilts by exp(-epsilon*t); the decreasing reward has positive covariance with this tilt",
            primary_FKG_source="https://link.springer.com/article/10.1007/BF01651330",
            validation_scope="FKG and the continuous grid/truncation passage are analytic inputs; the finite check verifies the equivalent word weights and capture normalization"),
        duration_discount_bounds="x(x+1)/(x+d+1)^2 <= F_d(x+1)/F_d(x) <= x/(x+d)*((x+d)/(x+d+1))^(d+1), d>=0, x>0",
        duration_bounds_checked_on_existing_entries=True,
        unresolved_conditions=[
            "P_l <=st Law(1+X), X~P_(l-1), is not proved",
            "R_(d+1)(x)>=R_d(x) for integer x>=1 is not proved"],
        full_tail_positivity_proved=False)


def highrate_factorial_volume_certificate(H, G, c):
    """Corrected covering ratios and the same m=3 parking volume.

    The n=4,k=1 scalar component uses only S_3; it is an obstruction
    to a termwise proof, not a new graph case or a counterexample to
    factorial log-convexity.
    """
    assert H.shape == G.shape == (4, 4) and c.shape == (4, 1)
    assert list(c) == [1, sp.Rational(1, 2), sp.Rational(2, 9), sp.Rational(17, 192)]
    assert H*c == c
    P = sp.diag(*c).inv()*H*sp.diag(*c)
    ratios, uncorrected = [], {}
    for ell in range(1, 4):
        for h in range(ell):
            b_ell = ell*c[ell]/c[ell-1]
            b_x = (h+1)*c[h+1]/c[h]
            actual = P[ell, h+1]/P[ell-1, h]
            uncorrected[ell, h] = b_x/b_ell*G[ell-1, h]/H[ell-1, h]
            corrected = sp.Rational(ell**2, (h+1)**2)*uncorrected[ell, h]
            assert sp.cancel(actual-corrected) == 0
            ratios.append([ell, h, str(actual)])
    actual_example = P[2, 1]/P[1, 0]
    assert actual_example == sp.Rational(1, 2)
    assert uncorrected[2, 0] == sp.Rational(1, 8)
    false_discount = G[2, 0]/H[2, 0]
    assert false_discount == sp.Rational(21, 160) > sp.Rational(1, 9)

    # Pitman--Stanley Theorem 1: the five Catalan compositions for S_3.
    w, v = sp.symbols('w v', nonnegative=True)
    x2, x3 = sp.Rational(1, 12), sp.Rational(1, 6)
    V3 = w*x2*x3+w*w*x2/2+w*x2*x2/2+w*w*x3/2+w**3/6
    assert sp.expand(V3) == 5*w/288+w*w/8+w**3/6
    assert V3.subs(w, sp.Rational(1, 4)) == c[3]/6
    pi = [sp.Rational(5, 17), sp.Rational(9, 17), sp.Rational(3, 17)]
    assert sum(pi) == 1
    a = [c[j]/sp.factorial(j) for j in range(4)]
    Q3 = V3.subs(w, sp.Rational(1, 4)-v)
    assert sp.expand(Q3-sum((-v)**j*a[3-j]/sp.factorial(j) for j in range(4))) == 0
    for j in range(4):
        moment = sum(pi[k-1]*sp.factorial(k)/sp.factorial(k-j)
                     for k in (1, 2, 3) if k >= j)
        assert moment == a[3-j]/(4**j*a[3])

    # This actual positive-basis component has negative individual gap;
    # the mixture criterion needs compensation between components.
    n, k = 4, 1
    A = (1-sp.Rational(1, n+1)**(k+1))/(n*(k+1))
    B = (1-sp.Rational(2, n+2)**(k+2))/(n*n*(k+1)*(k+2)) \
        -sp.Rational(1, n+1)**(k+1)/(n*(n+2)*(k+1))
    direct_A = sp.integrate((1-n*v)**k, (v, 0, sp.Rational(1, n+1)))
    direct_B = sp.integrate(v*(1-n*v)**k, (v, 0, sp.Rational(1, n+2))) \
        +sp.Rational(1, n+2)*sp.integrate(
            (1-n*v)**k, (v, sp.Rational(1, n+2), sp.Rational(1, n+1)))
    assert A == direct_A == sp.Rational(3, 25)
    assert B == direct_B == sp.Rational(149, 16200)
    gap = B-sp.Rational(n, n+1)**2*A*A
    assert gap == -sp.Rational(187, 10125000)
    return dict(
        scope="existing m=3 limit plus one formal parking component",
        corrected_covering_ratio_checks=ratios,
        missing_factor_example_actual=str(actual_example),
        missing_factor_example_uncorrected=str(uncorrected[2, 0]),
        false_discount_shortcut=dict(actual=str(false_discount), claimed_bound="1/9"),
        pitman_stanley_a3=str(a[3]), parking_pi4=list(map(str, pi)),
        scalar_integrals_checked=True, parking_factorial_moments_checked=True,
        componentwise_gap_n4_k1=str(gap),
        factorial_logconvexity_proved=False, covering_proved=False,
        corrected_covering_ratio="P_l(h+1)/P_(l-1)(h)=[l^2/(h+1)^2]*(b_(h+1)/b_l)*(G_(l-1,h)/H_(l-1,h)); b_l=l*c_l/c_(l-1)",
        correction_consequence="Monotonicity of b and G/H<=1 alone do not establish pointwise covering",
        primary_volume_source="https://live-statistics.pantheon.berkeley.edu/sites/default/files/tech-reports/560.pdf",
        volume_specialization="Pitman--Stanley Theorem 1: a_n=c_n/n!=V_n(1/(n+1),1/[n(n+1)],...,1/6); no log-convexity theorem is imported",
        parking_criterion=dict(
            distribution="K counts coordinates in [0,1/n] under uniform measure on S_(n-1); Q_(n-1)(v)/a_(n-1)=E(1-nv)^K",
            A_nk="[1-(n+1)^(-k-1)]/[n(k+1)]",
            B_nk="[1-(2/(n+2))^(k+2)]/[n^2(k+1)(k+2)]-(n+1)^(-k-1)/[n(n+2)(k+1)]",
            exact_target="b_(n+1)>=b_n iff E B_(n,K)>=[n/(n+1)]^2*(E A_(n,K))^2",
            factorial_moments="E(K)_j=a_(n-1-j)/(n^j*a_(n-1)), 0<=j<=n-1",
            exact_gap="E[B-gamma*A^2]+gamma*Var(A), gamma=[n/(n+1)]^2",
            obstruction_scope="The negative n=4,k=1 component rejects termwise positivity; the actual mixture inequality remains open"))


def highrate_limit_certificate():
    """One m=3 discriminator plus exact structural identity residuals.

    General positivity statements have analytic proofs; this is not a
    size/rate sweep or evidence promoted to an all-size C theorem.
    """
    def zero(matrix):
        assert all(sp.cancel(value) == 0 for value in matrix)

    def rows(matrix):
        return [[str(value) for value in row] for row in matrix.tolist()]

    def apply_jet(matrices, vectors):
        return [sum((matrices[j]*vectors[k-j] for j in range(k+1)),
                    sp.zeros(size, 1)) for k in range(order+1)]

    m, order = 3, 3
    size = m+1
    Ejet, Hjet, nujet = count_highrate_taylor(m, order)
    alpha = sp.Rational(m+2, m+1)
    c = highrate_basis(size, 1).inv()[:, 0]
    assert list(c) == [1, sp.Rational(1, 2), sp.Rational(2, 9), sp.Rational(17, 192)]
    for ell in range(size):
        i = m-ell
        assert all(nujet[k][i] == 0 for k in range(ell))
        assert nujet[ell][i] == sp.binomial(m, ell)*c[ell]
        for h in range(ell+1):
            assert all(Ejet[k][i, m-h] == 0 for k in range(ell-h))
            assert all(Hjet[k][i, m-h] == 0 for k in range(ell-h))
    for actual, expected in zip(apply_jet(Hjet, nujet), nujet):
        zero(actual-expected)
    for k, matrix in enumerate(Hjet):
        zero(sp.ones(1, size)*matrix-(sp.ones(1, size) if k == 0 else sp.zeros(1, size)))

    raw_E = sp.Matrix(size, size, lambda ell,h:
                      Ejet[ell-h][m-ell, m-h]*sp.binomial(m, h)/sp.binomial(m, ell)
                      if ell >= h else 0)
    raw_H = sp.Matrix(size, size, lambda ell,h:
                      Hjet[ell-h][m-ell, m-h]*sp.binomial(m, h)/sp.binomial(m, ell)
                      if ell >= h else 0)
    H, G = highrate_hub_limit(size, 1), highrate_hub_limit(size, 2)
    R = H.inv()*G
    K = sp.eye(size)-R
    Estar = sp.eye(size)-alpha*R
    zero(raw_E-Estar)
    zero(raw_H-H)
    zero(H*c-c)
    assert all(Estar[ell, ell] == sp.Rational(m-ell, (m+1)*(ell+2))
               for ell in range(size))

    P = sp.zeros(size)
    for i in range(size):
        P[m-i, i] = 1
    W = sp.diag(*c)
    Z = sp.Matrix(size, size, lambda ell,h: 1 if ell >= h else 0)
    F_count = P*W.inv()*Estar*W*P
    C_count = P*Z.inv()*W.inv()*Estar*W*Z*P
    expected_C = sp.Matrix([
        [0, sp.Rational(203, 1360), sp.Rational(7421, 54400), sp.Rational(145689, 1088000)],
        [0, sp.Rational(1, 16), sp.Rational(65, 384), sp.Rational(635, 4608)],
        [0, 0, sp.Rational(1, 6), sp.Rational(5, 36)],
        [0, 0, 0, sp.Rational(3, 8)]])
    zero(C_count-expected_C)
    assert all(C_count[i, j] > 0 for i in range(m) for j in range(i+1, size))
    first_super = []
    for i in range(m):
        ell = m-i
        predicted = alpha*(ell*(2*ell+3)*c[ell-1]/c[ell]-(ell+1)*(ell+2)) \
                    / ((ell+1)**2*(ell+2)**2)
        lower = alpha*sp.Rational(3*ell**2+3*ell-2, (ell+1)**2*(ell+2)**2)
        assert C_count[i, i+1] == predicted and predicted >= lower > 0
        first_super.append(dict(row=i, vacancy_defect=ell,
                                C_limit=str(predicted), rigorous_lower_bound=str(lower)))

    J = sp.zeros(size)
    for ell in range(m):
        J[ell+1, ell] = ell+1
    D1, D2 = sp.diag(*range(1, size+1)), sp.diag(*range(2, size+2))
    zero(H*D1-sp.eye(size)-H**2*J)
    zero(G*D2-sp.eye(size)-G**2*J)
    zero(K*D2-sp.eye(size)-H*K*J-K*G*J)
    zero(K-H*D1+D1*G)
    zero(H*J-J*G)
    assert all(K[ell, h] > 0 for ell in range(size) for h in range(ell+1))

    a = sp.symbols('a', positive=True)
    Da = sp.diag(*[a+ell for ell in range(size)])
    B = Da-J
    Pi = sp.Matrix(size, size, lambda ell,h: sp.binomial(ell, h) if ell >= h else 0)
    M = sp.Matrix.hstack(*[B**h*sp.eye(size)[:, h] for h in range(size)])
    zero(B*Pi-Pi*Da)
    zero(Pi.inv()*M-sp.diag(*[(a+ell)**ell for ell in range(size)])*highrate_basis(size, a))
    zero(B*M*Da-B**2*M-M*J)

    tau = sp.symbols('tau', nonnegative=True)
    occupied = (sp.eye(size)+tau*R).inv()*c
    empty = tau*G*occupied
    zero(D1*highrate_basis(size, 1)*empty+highrate_basis(size, 1)*occupied-sp.eye(size)[:, 0])
    zero(D2*highrate_basis(size, 2)*empty-tau*highrate_basis(size, 2)*occupied)
    assert sp.cancel(occupied[0]-2/(2+tau)) == 0
    assert sp.cancel(empty[0]-tau/(2+tau)) == 0
    c_next = highrate_basis(size+1, 1).inv()[1:, 0]
    zero(G*c-c_next)

    # Same m=3 matrices: structural commutator identity and the dominant
    # coefficients used in the separate analytic uniform-in-k argument.
    Q = K*H-H*K
    zero(K-G-H*K*J*G)
    zero(K-sum((H**r*G*(J*G)**r for r in range(m+1)), sp.zeros(size)))
    zero(Q-H*Q*J*H-H*K*H*K*J*H)
    zero(Q-sum((H**r*K*H*K*(J*H)**r for r in range(1, m+1)), sp.zeros(size)))
    assert all(Q[ell, h] > 0 for ell in range(size) for h in range(ell))
    assert all(Q[ell, ell-1] == sp.Rational(1, (ell+1)**2*(ell+2)**2)
               for ell in range(1, size))
    d0 = Estar[0, 0]
    perron = sp.ones(size, 1)
    for ell in range(1, size):
        perron[ell] = sum(Estar[ell, h]*perron[h] for h in range(ell)) \
                      / (d0-Estar[ell, ell])
    zero(Estar*perron-d0*perron)
    left_dominant = (sp.eye(size)-H)*perron
    right_dominant = H*perron-c
    left_k1 = (sp.eye(size)-H)*Estar*c
    right_k1 = H*Estar*c-d0*c
    zero(left_k1-alpha*Q*c)
    assert all(value > 0 for vector in (left_dominant, right_dominant, left_k1, right_k1)
               for value in vector[1:])

    # Independently check the sign and normalization of both leading tail
    # gaps at k=1 directly from the actual count Taylor jets.
    wjet = apply_jet(Ejet, nujet)
    hwjet = apply_jet(Hjet, wjet)
    gaps = []
    for ell in range(1, size):
        threshold = m-ell+1
        left = [sum(hwjet[k][i]-wjet[k][i] for i in range(threshold, size))
                for k in range(order+1)]
        right = [sum(sum(wjet[j])*sum(nujet[k-j][i] for i in range(threshold, size))
                     for j in range(k+1))-sum(hwjet[k][i] for i in range(threshold, size))
                 for k in range(order+1)]
        assert all(left[k] == right[k] == 0 for k in range(ell))
        assert left[ell] == sp.binomial(m, ell)*left_k1[ell] > 0
        assert right[ell] == sp.binomial(m, ell)*right_k1[ell] > 0
        gaps.append(dict(vacancy_cut=ell, occupation_threshold=threshold,
                         left_q_coefficient=str(left[ell]), right_q_coefficient=str(right[ell])))

    return dict(
        scope="single m=3 count Taylor check and structural algebra residuals; all-size C positivity remains open",
        outside_sites=m, Taylor_order=order, expansion_parameter="q=1-p",
        vacancy_coefficients=list(map(str, c)),
        first_superdiagonal=first_super,
        first_superdiagonal_theorem="For each fixed m, all first-superdiagonal C entries are positive sufficiently near p=1; c_l<=c_(l-1)/2 gives the displayed lower bound",
        basis="A_a[l,h]=binom(l,h)(-1/(a+l))^(l-h), l>=h; Hcal_a=A_a^-1 diag(1/(a+l)) A_a",
        full_limit="With D_q=diag(binom(m,l)q^l), D_q^-1 P E P D_q -> I-alpha Hcal_1^-1 Hcal_2; alpha=(m+2)/(m+1)",
        scaled_excursion_limit_vacancy_order=rows(Estar),
        weighted_excursion_limit_count_order=rows(F_count),
        C_limit_count_order=rows(C_count), positive_strict_upper_C_entries=6,
        unresolved_tail_criterion="T[l,r]=sum_(h=r)^l Kcal[l,h]c_h/c_l; require T[l,r]>T[l-1,r] for 0<=r<l; Kcal=I-Hcal_1^-1 Hcal_2",
        positive_capture_recurrence="Kcal D_2=I+Hcal_1 Kcal J+Kcal Hcal_2 J",
        finite_root_joint_law="tau=alpha/lambda_v, R=Hcal_1^-1 Hcal_2; b=(I+tau R)^-1 c, a=tau Hcal_2 b",
        next_volume_identity="Hcal_2 c=(c_1,...,c_(m+1))^T",
        residuals=dict(raw_E_all_leading_entries=True, raw_H_all_leading_entries=True,
                       raw_stationarity_through_order=True, quadratic_equations=True,
                       positive_capture_recurrence=True, symbolic_pascal_diagonalization=True,
                       finite_root_joint_harmonic_equations=True, next_volume_identity=True),
        termwise_weighted_capture_obstruction=dict(
            K10_c0_over_c1=str(K[1, 0]/c[1]), K20_c0_over_c2=str(K[2, 0]/c[2]),
            second_minus_first=str(K[2, 0]/c[2]-K[1, 0]/c[1])),
        highrate_sandwich_structure=dict(
            analytic_status="working proof: positive commutator plus a strictly separated dominant eigenvalue gives the full high-rate sandwich uniformly over k, without full C positivity",
            scope="fixed m, sufficiently large common outside rate, arbitrary positive finite inserted rate; one-coordinate NA, not full NA between blocks",
            commutator_identity="Q=Kcal Hcal_1-Hcal_1 Kcal=Hcal_1 Q J Hcal_1+Hcal_1 Kcal Hcal_1 Kcal J Hcal_1",
            finite_positive_expansion="Q=sum_(r=1)^m Hcal_1^r Kcal Hcal_1 Kcal (J Hcal_1)^r",
            commutator_vacancy_order=rows(Q), commutator_residuals_zero=True,
            dominant_eigenvalue=str(d0),
            left_dominant_coefficients=list(map(str, left_dominant[1:])),
            right_dominant_coefficients=list(map(str, right_dominant[1:])),
            raw_k1_tail_leading_coefficients=gaps),
        queue_tail=highrate_queue_tail_certificate(H, G, K, c),
        factorial_volume=highrate_factorial_volume_certificate(H, G, c))


def finite_capture_commutator_check():
    """One saved actual case: retain the finite-rate birth remainder."""
    m,p = 3,sp.Rational(4,9)
    q = 1-p
    _,_,H,E,nu = count_recursions(m,p)
    G = discounted_hub_kernel(m,p,2)
    birth,death = sp.zeros(m+1),sp.zeros(m+1)
    for j in range(m+1):
        if j<m: birth[j+1,j]=p*(m-j)
        if j: death[j-1,j]=q*j
    identity = sp.eye(m+1)
    R = identity+birth
    U,T = R.inv()*H,R.inv()*G
    K = identity-H.inv()*G
    assert K==(identity+(m+1)*E)/(m+2)
    assert min(U)>=0 and min(T)>=0
    Q = K*H-H*K
    A,S = T*R-R*T,death*G-H*death
    remainder = (identity+H*K*death)*A*U+H*K*S*U
    assert K==T+H*K*death*T
    assert Q==H*Q*death*U+H*K*H*K*death*U+remainder
    flux = [H[j+1,j]*K[j,j+1] for j in range(m)]
    assert min(flux)>0
    for j in range(m+1):
        assert Q[j,j]==(flux[j] if j<m else 0)-(flux[j-1] if j else 0)
    assert Q[m,m]==-sp.Rational(157757915,4045267824)
    assert (sp.ones(1,m+1)*Q*nu)[0]==0
    return dict(
        scope="one saved actual complete-graph case; no rate or size sweep",
        outside_sites=m,original_common_rate="1",refresh_parameter=str(p),
        definitions="R=I+B; U=R^-1 H; T=R^-1 G; K=I-H^-1 G; Q=KH-HK; A=TR-RT; S=DG-HD",
        finite_remainder="Q=H Q D U+H K H K D U+(I+H K D)A U+H K S U",
        positive_birth_flux=list(map(str,flux)),
        negative_full_count_diagonal=str(Q[m,m]),
        capture_identity_residual_zero=True,full_remainder_identity_residual_zero=True,
        diagonal_flux_identity_residual_zero=True,stationary_total_source_zero=True,
        consequence="Positive diagonal weighting or vacancy scaling cannot make the finite-rate commutator entrywise nonnegative; the excursion sandwich is not refuted.")


def capture_gradient_obstructions():
    """Two specified families separate gradient killing from capture order."""
    M,eps,a,theta = sp.symbols("M eps a theta",positive=True)

    def top_entry(v,u):
        u1,u2 = u.subs(a,1),u.subs(a,2)
        first = u1*(u1-sp.diff(u,a).subs(a,2))/(u2*(u1-sp.diff(u,a).subs(a,1)))
        return sp.factor(first-v.subs(a,1)/v.subs(a,2)),first

    v = a+2*M-M/(a+1)
    u = a+1-M/v
    unit,unit_first = top_entry(v,u)
    unit_expected = -(10*M**4+19*M**3-42*M**2-78*M-36) \
                    /(6*(2*M+3)*(5*M+6)*(M**2+2*M+2))
    assert sp.cancel(unit-unit_expected)==0
    assert sp.limit(unit,M,sp.oo)==-sp.Rational(1,6)
    assert sp.limit(unit_first,M,sp.oo)==sp.Rational(11,15)
    assert sp.limit(v.subs(a,1)/v.subs(a,2),M,sp.oo)==sp.Rational(9,10)
    assert sp.cancel(sp.limit(u,M,sp.oo)-2*a*(a+1)/(2*a+1))==0

    phase = sp.Matrix([[-1,1,0],[M,-2*M,M],[0,1,-1]])
    dual = sp.Matrix([[-M,M],[M,-M]])
    gradient = sp.Matrix([[-1,1,0],[0,-1,1]])
    assert gradient*phase==(dual-sp.eye(2))*gradient
    assert sp.Matrix([[-1,1]])*dual==-2*M*sp.Matrix([[-1,1]])
    average = sp.Matrix([[1-theta,theta,0],[0,1-theta,theta]])
    defect = dual*average-average*phase
    assert sp.expand(defect[0,0]-(1-theta-M))==0
    assert sp.expand(defect[1,2]-(theta-M))==0
    assert defect.subs({M:sp.Rational(1,2),theta:sp.Rational(1,2)})==sp.zeros(2,3)

    v = a+2-eps/(a+eps)
    u = a+eps-eps/v
    bounded,bounded_first = top_entry(v,u)
    bounded_expected = (3*eps**6+34*eps**5+151*eps**4+307*eps**3
                        +226*eps**2-81*eps-132) \
        /(2*(eps+1)*(eps+2)*(eps+4)*(3*eps+8)*(eps**2+4*eps+5))
    assert sp.cancel(bounded-bounded_expected)==0
    assert sp.limit(bounded,eps,0,dir="+")==-sp.Rational(33,160)
    assert sp.limit(bounded_first,eps,0,dir="+")==sp.Rational(87,160)
    assert sp.limit(v.subs(a,1)/v.subs(a,2),eps,0,dir="+")==sp.Rational(3,4)
    kappas = [a+eps,v,u]
    capture_diagonal = [1-k.subs(a,1)/k.subs(a,2) for k in kappas]
    assert [sp.limit(x,eps,0,dir="+") for x in capture_diagonal] \
           ==[sp.Rational(1,2),sp.Rational(1,4),sp.Rational(1,2)]
    middle = sp.factor((4*capture_diagonal[1]-1)/3)
    assert sp.cancel(middle-eps*(eps+5)/(3*(eps+1)*(3*eps+8)))==0

    def point(beta,delta):
        """Independent full six-state upward Green calculation."""
        beta,delta = list(map(sp.Rational,beta)),list(map(sp.Rational,delta))
        states = [(b,y) for b in range(1,4) for y in range(b)]
        indices = {state:i for i,state in enumerate(states)}
        Q,success = sp.zeros(6),sp.zeros(6,1)
        for (b,y),i in indices.items():
            Q[i,i]=-(1+beta[y]+delta[y])
            if y+1<b: Q[i,indices[b,y+1]]=beta[y]
            if y: Q[i,indices[b,y-1]]=delta[y]
            if b<3: Q[i,indices[b+1,y]]=1
            else: success[i]=1
        boundary = [indices[j+1,j] for j in range(3)]
        N1,N2 = (-Q).inv(),(sp.eye(6)-Q).inv()
        G1,G2 = N1.extract(boundary,boundary),N2.extract(boundary,boundary)
        pi = sp.Matrix([1,beta[0]/delta[1],beta[0]*beta[1]/(delta[1]*delta[2])])
        pi /= sum(pi)
        Pi = sp.diag(*pi)
        f = (N1*success).extract(boundary,[0])
        nu = Pi*G1.inv()*f
        birth = sp.zeros(3)
        for j in range(2): birth[j+1,j]=beta[j]
        H = (sp.eye(3)+birth)*Pi*G1*Pi.inv()
        H2 = (sp.eye(3)+birth)*Pi*G2*Pi.inv()
        assert H*nu==nu and sum(nu)==1 and min(nu)>0
        assert sp.ones(1,3)*H==sp.ones(1,3)
        Ehat = Pi*(sp.eye(3)-G1.inv()*G2)*Pi.inv()
        assert Ehat==sp.eye(3)-H.inv()*H2 and min(Ehat)>=0
        E = (4*Ehat-sp.eye(3))/3
        J = H*(sp.eye(3)-E)
        assert J==sp.Rational(4,3)*H2 and min(J)>=0
        assert sp.ones(1,3)*(E+J)==sp.ones(1,3)
        L = sp.Matrix(3,3,lambda i,j:int(bool(i<=j)))
        W = sp.diag(*nu)
        Chat = L.inv()*W.inv()*Ehat*W*L
        C = L.inv()*W.inv()*E*W*L
        left = ((H-sp.eye(3))*E*nu)[2]
        assert C[1,2]==sp.Rational(4,3)*Chat[1,2]<0
        assert left==H[2,1]*nu[1]*C[1,2]<0
        assert all(G1[i,j]/f[i]<=G1[i+1,j]/f[i+1]
                   for i in range(2) for j in range(i+1,3))
        masses = list(sp.ones(1,3)*E)
        return dict(
            beta=list(map(str,beta)),delta=list(map(str,delta)),
            gradient_killing=[str(beta[j]-beta[j+1]+delta[j+1]-delta[j]) for j in range(2)],
            per_count_rate_bounds=all(beta[j]<=2-j and delta[j]<=j for j in range(3)),
            nu=list(map(str,nu)),
            corrected_E=[[str(x) for x in E.row(i)] for i in range(3)],
            corrected_E_nonnegative=bool(min(E)>=0),
            corrected_E_column_masses=list(map(str,masses)),
            corrected_E_strictly_substochastic=bool(min(E)>=0 and max(masses)<1),
            J_nonnegative=True,E_plus_J_stochastic=True,stationarity_verified=True,
            planar_Green_ratio_verified=True,independent_six_state_check=True,
            negative_Chat_top=str(Chat[1,2]),negative_C_top=str(C[1,2]),
            negative_full_event_left_gap_at_k1=str(left))

    unit_point = point([1,10,0],[0,10,1])
    bounded_point = point([sp.Rational(1,10),1,0],[0,1,sp.Rational(1,10)])
    assert unit_point["gradient_killing"]==["1","1"]
    assert unit.subs(M,10)==sp.Rational(unit_point["negative_Chat_top"])
    assert not unit_point["corrected_E_nonnegative"]
    assert bounded_point["gradient_killing"]==["1/10","1/10"]
    assert bounded_point["per_count_rate_bounds"] and bounded_point["corrected_E_strictly_substochastic"]
    assert bounded.subs(eps,sp.Rational(1,10))==sp.Rational(bounded_point["negative_Chat_top"])
    return dict(
        scope="Two distinct targeted generic birth-death families; neither is an ARW/Ehrenfest counterexample.",
        unit_gradient_killing=dict(
            family="beta=(1,M,0), delta=(0,M,1), M>0",
            top_entry=str(unit),limit="-1/6",first_term_limit="11/15",v_ratio_limit="9/10",
            first_gradient_intertwining_verified=True,next_gradient_killing="2*M",
            stationary_Bernoulli_averaging_obstruction="Its two corner equations force theta=1-M and theta=M; closure occurs only at M=theta=1/2.",
            point=unit_point),
        bounded_rates_nonnegative_excursion=dict(
            family="beta=(eps,1,0), delta=(0,1,eps), eps>0",
            top_entry=str(bounded),limit="-33/160",first_term_limit="87/160",v_ratio_limit="3/4",
            capture_diagonal_limit=["1/2","1/4","1/2"],corrected_middle_diagonal=str(middle),
            small_eps_argument="Off-diagonal E entries are positive; E00=(2-eps)/(3*(2+eps)), E11=eps*(eps+5)/(3*(eps+1)*(3*eps+8)), and E22 tends to 1/3. Thus E>=0 for small positive eps, while Chat12 remains negative.",
            consequence="Even nonnegative gradient killing, count-rate bounds, nonnegative strictly substochastic E, J>=0, E+J stochastic, and Hnu=nu do not imply the generic excursion sandwich: its full-event left gap already fails at k=1.",
            point=bounded_point))


def ehrenfest_boundary_identities():
    """Audit the proved boundary identities at the single indexing case m=2."""
    def upward(m,p):
        states = [(b,y) for b in range(1,m+2) for y in range(b)]
        indices = {state:i for i,state in enumerate(states)}
        Q,success = sp.zeros(len(states)),sp.zeros(len(states),1)
        for (b,y),i in indices.items():
            birth,death = p*(m-y),(1-p)*y
            Q[i,i] = -(1+birth+death)
            if y+1<b: Q[i,indices[b,y+1]]=birth
            if y: Q[i,indices[b,y-1]]=death
            if b<m+1: Q[i,indices[b+1,y]]=1
            else: success[i]=1
        return states,indices,-Q,success

    def zero(matrix):
        assert all(sp.expand(x)==0 for x in matrix)

    m,p = 2,sp.Symbol("p")
    q = 1-p
    new,ni,A,c = upward(m,p)
    old,oi,Aplus,cplus = upward(m+1,p)
    gradient,average,boundary_source = (
        sp.zeros(len(new),len(old)) for _ in range(3))
    for (b,y),i in ni.items():
        lo,hi = oi[b+1,y],oi[b+1,y+1]
        gradient[i,hi],gradient[i,lo] = 1,-1
        average[i,lo],average[i,hi] = q,p
        if y==b-1: boundary_source[i,hi]=p*(m-y)
    zero(gradient*Aplus-(A+sp.eye(len(new)))*gradient-boundary_source)
    zero(A*average-average*Aplus-q*boundary_source)
    zero(gradient*cplus)
    zero(average*cplus-c)
    assert boundary_source[ni[m+1,m],:]==sp.zeros(1,len(old))

    old_injection,new_injection = sp.zeros(len(old),m+2),sp.zeros(len(new),m+1)
    shift = sp.zeros(m+1,m+2)
    for j in range(m+2): old_injection[oi[j+1,j],j]=1
    for j in range(m+1):
        new_injection[ni[j+1,j],j]=1
        shift[j,j+1]=1
    zero(gradient*old_injection-new_injection*shift)
    zero(average*old_injection-p*new_injection*shift)

    at = {p:sp.Rational(1,2)}
    A0,Aplus0 = A.subs(at),Aplus.subs(at)
    average0,source0 = average.subs(at),boundary_source.subs(at)
    N1,N2 = A0.inv(),(A0+sp.eye(len(new))).inv()
    old_success,new_success = Aplus0.inv()*cplus,N1*c
    zero(-gradient*old_success-N2*source0*old_success)
    zero(average0*old_success-new_success-sp.Rational(1,2)*N1*source0*old_success)
    boundary = [ni[i+1,i] for i in range(m+1)]
    shifted = [oi[i+2,i+1] for i in range(m+1)]
    f,t = new_success.extract(boundary,[0]),old_success.extract(shifted,[0])
    G1,G2 = N1.extract(boundary,boundary),N2.extract(boundary,boundary)
    birth = sp.diag(*[sp.Rational(m-i,2) for i in range(m+1)])
    zero(t-f-sp.Rational(1,2)*(G1-G2)*birth*t)
    assert all(x>0 for x in (t-f)[:-1]) and t[-1]==f[-1]==1
    g = G1.inv()*f
    assert g[-1]>0
    old_boundary = [oi[i+1,i] for i in range(m+2)]
    old_Green = Aplus0.inv().extract(old_boundary,old_boundary)
    Gplus = old_Green[1:,1:]
    zero(Gplus-(G1+G2)/2-(G1-G2)*birth*Gplus/2)
    old_g = old_Green.inv()*old_success.extract(old_boundary,[0])
    gshift = old_g[1:,0]
    capture = sp.eye(m+1)-G1.inv()*G2
    transfer = sp.eye(m+1)-capture/2
    zero(g-transfer*gshift)
    zero(f-(G1+G2)*gshift/2)
    for r in range(m+1):
        prefix = sp.diag(*[int(j<=r) for j in range(m+1)])
        correction = prefix*capture*(sp.eye(m+1)-prefix)*gshift/2
        zero(prefix*g-transfer*prefix*gshift+correction)
        assert min(correction)>=0
    print("Ehrenfest gradient, Bernoulli averaging, and prefix corrections verified at m=2.",flush=True)
    return dict(
        scope="All-size identities follow from the analytic generator calculation; one m=2 symbolic indexing audit and one p=1/2 Green check, not a proof of capture order.",
        symbolic_rate="p indeterminate",exact_rate="1/2",
        gradient_identity=True,stationary_bit_averaging_identity=True,
        zero_gradient_success_source=True,top_sloping_boundary_source_zero=True,
        f=list(map(str,f)),shifted_fplus=list(map(str,t)),
        cross_size_resolvent_identity=True,top_capture_source=str(g[-1]),
        arbitrary_boundary_source_identity=True,cross_size_Green_matrix_identity=True,
        capture_stationary_source_recursion=True,all_prefix_corrections_verified=True,
        prefix_correction="P_r g=(I-q Etilde)P_r g_plus-q P_r Etilde(I-P_r)g_plus; the correction is nonnegative and is subtracted.",
        success_flux=list(map(str,g)),shifted_success_flux=list(map(str,gshift)))


def planar_capture_certificate():
    """Six-state Green audit and targeted generic-BD obstructions.

    The generic rates are not Ehrenfest rates. Only Ehat is interpreted
    as capture; no physical validity of launch-corrected generic E is used.
    """
    beta = [sp.Integer(1),sp.Integer(10),sp.Integer(0)]
    delta = [sp.Integer(0),sp.Integer(1),sp.Integer(1)]
    size = 3
    birth,death = sp.zeros(size),sp.zeros(size)
    for j in range(size):
        if j+1<size: birth[j+1,j]=beta[j]
        if j: death[j-1,j]=delta[j]

    def discounted(a):
        diagonal = sp.diag(*(a+beta[j]+delta[j] for j in range(size)))
        H = sp.zeros(size)
        for j in range(size):
            pivot = diagonal[j,j]-(delta[j]*H[j,j-1] if j else 0)
            column = sp.eye(size)[:,j]+birth[:,j]
            if j:
                column += delta[j]*sum(
                    (H[k,j-1]*H[:,k] for k in range(j)),sp.zeros(size,1))
            H[:,j] = column/pivot
        V = diagonal-H*death
        assert H*V==sp.eye(size)+birth
        return H,V,V.inv()

    H,V,U = discounted(1)
    H2,_,U2 = discounted(2)
    pi = sp.Matrix([1,1,10])/12
    nu = (V-death)*pi
    assert H*nu==nu and sum(nu)==1 and min(nu)>0
    Ehat = sp.eye(size)-V*U2
    assert Ehat==sp.eye(size)-H.inv()*H2 and min(Ehat)>=0
    prefix = sp.Matrix(size,size,lambda i,j:int(bool(i<=j)))
    W,Pi = sp.diag(*nu),sp.diag(*pi)
    criterion = prefix.inv()*W.inv()*Ehat*W*prefix
    assert criterion[1,2]==-sp.Rational(1115,30324)

    # Independent full row generator, including all three interior states.
    states = [(b,y) for b in range(1,size+1) for y in range(b)]
    indices = {state:i for i,state in enumerate(states)}
    Q = sp.zeros(len(states))
    success = sp.zeros(len(states),1)
    for (b,y),row in indices.items():
        Q[row,row]=-(1+beta[y]+delta[y])
        if y+1<b: Q[row,indices[b,y+1]]=beta[y]
        if y: Q[row,indices[b,y-1]]=delta[y]
        if b<size: Q[row,indices[b+1,y]]=1
        else: success[row]=1
    N1,N2 = (-Q).inv(),(sp.eye(len(states))-Q).inv()
    boundary = [indices[j+1,j] for j in range(size)]
    interior = [i for i in range(len(states)) if i not in boundary]
    G1,G2 = N1.extract(boundary,boundary),N2.extract(boundary,boundary)
    assert G1==Pi.inv()*U*Pi and G2==Pi.inv()*U2*Pi
    assert N1-N2==N1*N2
    R = Q.extract(boundary,interior)*(-Q.extract(interior,interior)).inv()
    Phi = sp.zeros(size,len(states))
    for i in range(size):
        Phi[i,boundary[i]]=1
        for j,col in enumerate(interior): Phi[i,col]=R[i,j]
    g = Pi.inv()*nu
    assert Phi*success==g
    bulk = N2.extract(list(range(len(states))),boundary)
    assert Phi*bulk==Pi.inv()*Ehat*Pi
    for r in range(size):
        source = sp.Matrix([g[j] if j<=r else 0 for j in range(size)])
        assert sp.diag(*g).inv()*Phi*bulk*source==W.inv()*Ehat*W*prefix[:,r]
    h = U*nu
    assert all(U[i,j]/h[i]<=U[i+1,j]/h[i+1]
               for i in range(size-1) for j in range(i+1,size))

    # This limit selected the single rational example above.
    M,a = sp.symbols("M a",positive=True)
    v = a+M+1-1/(a+1)
    u = a+1-M/v
    u1,u2 = u.subs(a,1),u.subs(a,2)
    top = sp.factor(u1*(u1-sp.diff(u,a).subs(a,2))
                    /(u2*(u1-sp.diff(u,a).subs(a,1)))-v.subs(a,1)/v.subs(a,2))
    expected = (-4*M**3+9*M**2+69*M+180)/(6*(M+4)*(M+9)*(3*M+8))
    assert sp.cancel(top-expected)==0 and top.subs(M,10)==criterion[1,2]
    assert sp.limit(top,M,sp.oo)==-sp.Rational(2,9)

    def rows(matrix):
        return [[str(value) for value in matrix.row(i)] for i in range(matrix.rows)]

    finite = finite_capture_commutator_check()
    gradient = capture_gradient_obstructions()
    boundary_identities = ehrenfest_boundary_identities()
    print("Capture Green identities, generic-BD/gradient obstructions, and finite-rate commutator remainder verified.",flush=True)
    return dict(
        scope="Generic birth-death capture obstruction and one saved actual commutator check; neither refutes Ehrenfest capture positivity or the excursion sandwich.",
        bulk_factorization=dict(
            definitions="N_a=((a-1)I-Q)^-1; D is the diagonal boundary, I its interior complement; Phi=[I,Q_DI(-Q_II)^-1]",
            identities="Pi^-1 Ehat Pi=Phi (N_2)_(all,D); g=Pi^-1 nu=Phi c",
            missing_inequality="Phi_i(z_r) Phi_(i+1)(c)>=Phi_(i+1)(z_r) Phi_i(c), where z_r=(N_2)_(all,D)(g*1_[0,r])",
            full_transient_states=len(states),independent_Green_check=True,
            resolvent_identity_verified=True,all_prefix_reward_identities_verified=True),
        generic_birth_death_obstruction=dict(
            scope="Ehat only; rates are not Ehrenfest and launch-corrected generic E is not asserted physical",
            beta=list(map(str,beta)),delta=list(map(str,delta)),
            pi=list(map(str,pi)),nu=list(map(str,nu)),Ehat=rows(Ehat),criterion=rows(criterion),
            negative_top_entry=str(criterion[1,2]),planar_Green_ratio_still_holds=True,
            family="beta=(1,M,0), delta=(0,1,1), M>0",family_top_entry=str(top),family_limit="-2/9"),
        finite_rate_commutator=finite,gradient_killing_obstructions=gradient,
        ehrenfest_boundary_identities=boundary_identities)


def volume_moment_upper_endpoint_check():
    """Exact negative-branch cusp recursion; no moment calculation."""
    q, z, r, x, h, k = sp.symbols("q z r x h k", real=True)
    denominator = 1-z*q*(1-q)
    qv = q*(1-q)/denominator
    first = q/denominator
    second = q*(1-q)*(1-z*q*q)/denominator**3
    assert sp.cancel(sp.diff(-sp.log(1-q), q)*qv-first) == 0
    assert sp.cancel(sp.diff(first, q)*qv-second) == 0
    Dr = r*r+x*(1-r*r)/4
    velocity = -(1-r*r)/(2*Dr)
    assert sp.cancel((-2*qv).subs({q: (1-r)/2, z: 4-x})-velocity) == 0
    W = 1+h*(1-r*r)
    A1 = (1-r)/(2*W)
    assert sp.cancel(first.subs({q: (1-r)/2, z: 4-x})
                     .subs(x, 4*r*r*h)-r**(-2)*A1) == 0
    A = sp.Function("A")(r, h)
    bracket = (1-3*k)*A+r*sp.diff(A, r)-2*h*sp.diff(A, h)
    next_coefficient = (1-r*r)*((3*k-1)*A-r*sp.diff(A, r)
                               +2*h*sp.diff(A, h))/(2*W)
    assert sp.simplify(-(1-r*r)*bracket/(2*W)-next_coefficient) == 0
    A2 = (1-r*r)*(2*A1-r*sp.diff(A1, r)+2*h*sp.diff(A1, h))/(2*W)
    assert sp.cancel(second.subs({q: (1-r)/2, z: 4-x})
                     .subs(x, 4*r*r*h)-r**(-5)*A2) == 0
    assert A1.subs({r: 0, h: 0}) == A2.subs({r: 0, h: 0}) == sp.Rational(1, 2)
    return dict(
        scope="symbolic negative-branch derivatives and general cusp recursion only",
        first_and_second_log_s_derivatives=True, cusp_velocity=True,
        general_derivative_recursion=True, second_derivative_recursion_agrees=True,
        leading_constants="a_1=1/2; a_(k+1)=(3k-1)*a_k/2",
        conditional_analytic_consequence="every positive representing measure has all inverse moments of 4-Z finite; mu([4-h,4])=o(h^p) and m_n/4^n=o(n^-p) for every p>0",
        continuation_across_cusp_used=False,
        positive_measure_constructed_or_excluded=False)


def volume_moment_transform_certificate():
    """Symbolic transform and endpoint identities, without a moment sweep.

    The four moments are the already proved escape volumes. Existence of
    a positive representing measure is not assumed by this calculation.
    The conditional endpoint theorem uses the analytic proof in the report.
    """
    t, z = sp.symbols("t z")
    y, r, epsilon, parameter = sp.symbols("y r epsilon parameter", positive=True)
    moments = [sp.S.One, sp.Rational(1, 2), sp.Rational(4, 9), sp.Rational(17, 32)]
    series = sum(sp.laguerre(n, (n+1)*z)*t**(n+1)/(n+1) for n in range(4))
    average = sp.expand(sum(sp.expand(series).coeff(z, k)*moments[k] for k in range(4)))
    assert average == t
    w = (1-sp.exp(-series)).series(t, 0, 5).removeO()
    residual = (w-t*sp.exp(-z*w/(1-w))).series(t, 0, 5).removeO().expand()
    assert residual == 0

    tau = y*sp.exp(z*y)/(1+y)
    a = z*y*(1+y)
    denominator = 1+a
    first = y/(tau*denominator)
    second = y*(y-(y+2)*a-a*a)/(tau*tau*denominator**3)
    assert sp.simplify(sp.diff(sp.log(1+y), y)/sp.diff(tau, y)-first) == 0
    assert sp.simplify(sp.diff(first, y)/sp.diff(tau, y)-second) == 0
    dy_dz = -sp.diff(tau, z)/sp.diff(tau, y)
    assert sp.simplify(dy_dz/(1+y)+y*y/denominator) == 0
    terminal = sp.expand(y-(y+2)*(1+parameter)-(1+parameter)**2)
    assert sp.expand(terminal+3+(y+4)*parameter+parameter**2) == 0

    A = sp.log(1+r)+r/(1+r)-epsilon
    q = sp.log(1+1/r)
    dr = r/A
    qe = sp.diff(q, r)*dr
    qee = sp.diff(qe, epsilon)+sp.diff(qe, r)*dr
    assert sp.simplify(qe+1/((1+r)*A)) == 0
    assert sp.simplify(qee-(epsilon+r-sp.log(1+r))/((1+r)**2*A**3)) == 0
    boundary = epsilon/((1+r)**2*A**3)
    regular = 1/((1+r)*A)-(r-sp.log(1+r))/((1+r)**2*A**3)
    assert sp.simplify(qee+qe-boundary+regular) == 0
    kernel_mass = sp.integrate(epsilon/(epsilon**2+4*z)**sp.Rational(3, 2), (z, 0, sp.oo))
    stieltjes_mass = sp.integrate((parameter+z)**(-sp.Rational(3, 2)), (z, 0, sp.oo))
    assert kernel_mass == sp.Rational(1, 2)
    assert sp.simplify(stieltjes_mass-2/sp.sqrt(parameter)) == 0
    print("Volume moment transform, endpoint decomposition, and scaling identities passed.", flush=True)
    return dict(
        scope="four supplied moments and exact symbolic identities; no new moments or measure inversion",
        supplied_moments=list(map(str, moments)), laguerre_series_average=str(average),
        implicit_series_residual_through_t4=str(residual),
        t_and_z_derivative_identities=True, epsilon_derivative_identities=True,
        terminal_second_derivative_numerator=str(terminal),
        positive_boundary_minus_regular_identity=True,
        boundary_kernel_lebesgue_mass=str(kernel_mass),
        generalized_stieltjes_lebesgue_transform=str(stieltjes_mass),
        analytic_theorem=dict(
            status="proved in report conditional on positive representability; not inferred from the finite check",
            support="[0,4], essential supremum 4; representing measure unique if it exists",
            CDF_at_zero="F(z)~kappa*z, kappa=-2*integral Q_tt(1,z) mu(dz) in (0,6]",
            density_assumption_required=False),
        positive_measure_constructed=False,
        upper_endpoint=volume_moment_upper_endpoint_check())


def main():
    default_output = Path(__file__).with_name("algebra_certificates.json")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=default_output)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument(
        "--allsize-structure", action="store_true",
        help="run only the count-recurrence and obstruction diagnostics")
    mode.add_argument(
        "--continuant-positivity", action="store_true",
        help="run only capture/top-entry identities and positivity-proof residuals")
    mode.add_argument(
        "--interior-pqd", action="store_true",
        help="run the general-m r=1 proof identities and five saved reverse-kernel diagnostics")
    mode.add_argument(
        "--reverse-monotonicity", action="store_true",
        help="run the all-rate reverse P/Q monotonicity proof checks at m<=3")
    mode.add_argument(
        "--highrate-limit", action="store_true",
        help="check high-rate, queue and volume formulas using one m=3 count Taylor case")
    mode.add_argument(
        "--planar-capture", action="store_true",
        help="check capture Green/boundary identities and the targeted positivity obstructions")
    mode.add_argument(
        "--volume-moment-transform", action="store_true",
        help="check the Laguerre transform and conditional zero-endpoint identities")
    args = parser.parse_args()
    if args.allsize_structure:
        saved = json.loads(default_output.read_text())
        cases = saved["excursion_count_criterion"]["cases"]
        key = "excursion_allsize_structure"
        certificate = allsize_structure_certificate(cases)
    elif args.continuant_positivity:
        key = "excursion_continuant_positivity"
        certificate = continuant_positivity_certificate()
    elif args.interior_pqd:
        saved = json.loads(default_output.read_text())
        key = "excursion_interior_pqd"
        certificate = interior_pqd_certificate(saved)
    elif args.reverse_monotonicity:
        key = "excursion_reverse_monotonicity"
        certificate = reverse_monotonicity_certificate()
    elif args.highrate_limit:
        key = "excursion_highrate_limit"
        certificate = highrate_limit_certificate()
    elif args.planar_capture:
        key = "excursion_planar_capture"
        certificate = planar_capture_certificate()
    elif args.volume_moment_transform:
        key = "volume_moment_transform"
        certificate = volume_moment_transform_certificate()
    else:
        key = "excursion_count_criterion"
        certificate = count_criterion_certificate()
    # Reload immediately before writing, preserving other certificate entries.
    existing = json.loads(args.output.read_text()) if args.output.exists() else {}
    existing[key] = certificate
    args.output.write_text(json.dumps(existing, indent=2)+"\n")


if __name__ == "__main__":
    main()
