"""Exact marker/forest calculations for response-renewal §§15,17,19.

Four fixed examples only; standard library; no files written.  Run with
python3 -B research/arw/monotonicity/verify_marker_geometry.py

Checks the released-input marker law by two finite constructions: independent
latent subset plus exit forest, and loop-range-dependent thinning of that same
forest.  Includes a reversible three-site counterexample to strong Rayleigh
stability and an exact symbolic certificate for ordinary Rayleigh inequalities
on every three-site root domain.  One heterogeneous four-cycle checks the
arbitrary-cycle marker formula, including self-loops and two exit labels.
This does not establish an ARW occupation inequality or Rayleigh inequalities
in arbitrary larger strongly connected components.

AI-written research calculation, not a claim of novelty or user review.
"""

from collections import defaultdict
from fractions import Fraction as F
from itertools import combinations, product
from math import prod


def subsets(vertices):
    vertices = tuple(sorted(vertices))
    return [frozenset(c) for k in range(len(vertices) + 1)
            for c in combinations(vertices, k)]


def determinant(matrix):
    a = [list(map(F, row)) for row in matrix]
    result = F(1)
    for j in range(len(a)):
        pivot = next((i for i in range(j, len(a)) if a[i][j]), None)
        if pivot is None:
            return F(0)
        if pivot != j:
            a[j], a[pivot] = a[pivot], a[j]
            result = -result
        value = a[j][j]
        result *= value
        for i in range(j + 1, len(a)):
            ratio = a[i][j] / value
            for k in range(j + 1, len(a)):
                a[i][k] -= ratio * a[j][k]
    return result


def inverse(matrix):
    n = len(matrix)
    a = [[F(x) for x in row] + [F(i == j) for j in range(n)]
         for i, row in enumerate(matrix)]
    for j in range(n):
        pivot = next(i for i in range(j, n) if a[i][j])
        a[j], a[pivot] = a[pivot], a[j]
        value = a[j][j]
        a[j] = [x / value for x in a[j]]
        for i in range(n):
            if i != j:
                value = a[i][j]
                a[i] = [x - value * y for x, y in zip(a[i], a[j])]
    return [row[n:] for row in a]


def minor(matrix, indices):
    indices = sorted(indices)
    return determinant([[matrix[i][j] for j in indices] for i in indices])


def exit_forests(K):
    """Enumerate the at most 27 positive-arrow assignments in our examples."""
    n = len(K)
    choices = [[j for j in range(n + 1) if K[i][j]] for i in range(n)]
    for arrows in product(*choices):
        for start in range(n):
            seen, vertex = set(), start
            while vertex < n and vertex not in seen:
                seen.add(vertex)
                vertex = arrows[vertex]
            if vertex < n:
                break
        else:
            roots = frozenset(i for i, j in enumerate(arrows) if j == n)
            yield roots, prod(K[i][j] for i, j in enumerate(arrows))


def audit(K, p):
    n, Q = len(K), frozenset(range(len(K)))
    sets = subsets(Q)
    r = [(1 - x) / x for x in p]
    M = [[F(i == j) - K[i][j] for j in range(n)] for i in range(n)]
    d = {C: minor(M, C) for C in sets}
    assert all(x > 0 for x in d.values())
    Z = {R: d[R] * sum(prod(r[i] for i in R - C) / d[C]
                       for C in subsets(R)) for R in sets}
    nu = {C: d[Q] * prod(r[i] for i in Q - C) / (Z[Q] * d[C])
          for C in sets}
    assert sum(nu.values()) == 1
    G = inverse(M)
    H = [[G[i][j] * K[j][n] for j in range(n)] for i in range(n)]
    forests = list(exit_forests(K))
    assert sum(w for _, w in forests) == d[Q]
    forests = [(D, w / d[Q]) for D, w in forests]

    marker_mass = defaultdict(F)
    for C, cweight in nu.items():
        for D, fweight in forests:
            marker_mass[D - C] += cweight * fweight
    moments = {S: prod(r[i] * K[i][n] for i in S) * Z[Q - S] / Z[Q]
               for S in sets}
    for S in sets:
        assert sum(w for D, w in forests if S <= D) == minor(H, S)
        assert sum(w for E, w in marker_mass.items() if S <= E) == moments[S]

    # g(R) sums the weights p_support(alpha) w(alpha) of all balanced
    # words supported in R.  Inclusion-exclusion therefore gives the
    # exact range distribution without truncating the infinitely many words.
    g = {R: prod(p[i] for i in R) * Z[R] / d[R] for R in sets}
    loop_range = {L: sum((-1)**len(L - A) * g[A] for A in subsets(L)) / g[Q]
                  for L in sets}
    assert min(loop_range.values()) >= 0 and sum(loop_range.values()) == 1
    thinned_mass = defaultdict(F)
    for L, lweight in loop_range.items():
        conditional = defaultdict(F)
        for D, fweight in forests:
            eligible = D - L
            for E in subsets(eligible):
                weight = (fweight * prod(1 - p[i] for i in E)
                          * prod(p[i] for i in eligible - E))
                conditional[E] += weight
                thinned_mass[E] += lweight * weight
        assert sum(conditional.values()) == 1
        for S in sets:
            expected = (prod(1 - p[i] for i in S) * minor(H, S)
                        if not S & L else F(0))
            assert sum(w for E, w in conditional.items() if S <= E) == expected
    assert all(thinned_mass[S] == marker_mass[S] for S in sets)

    # Joint loss marks E_A and uncoloured forest roots D_B.  Formula applies
    # to disjoint A,B; arrows to a deleted vertex remain labelled exits.
    joint_checks = 0
    for A in sets:
        for B in subsets(Q - A):
            measured = sum(cw * fw for C, cw in nu.items() for D, fw in forests
                           if A <= D - C and B <= D)
            expected = (moments[A] * prod(K[i][n] for i in B)
                        * d[Q - A - B] / d[Q - A])
            assert measured == expected
            joint_checks += 1

    # Conditional pairwise negative correlation, including absent markers.
    conditional_checks = 0
    for assignment in product((-1, 0, 1), repeat=n):
        yes = {i for i, v in enumerate(assignment) if v == 1}
        no = {i for i, v in enumerate(assignment) if v == 0}
        free = [i for i, v in enumerate(assignment) if v == -1]
        rows = [(E, w) for E, w in marker_mass.items() if yes <= E and not no & E]
        total = sum(w for _, w in rows)
        if not total:
            continue
        for i, j in combinations(free, 2):
            mi = sum(w for E, w in rows if i in E)
            mj = sum(w for E, w in rows if j in E)
            mij = sum(w for E, w in rows if {i, j} <= E)
            assert total * mij <= mi * mj
            conditional_checks += 1
    return dict(Z=Z, mass=dict(marker_mass), moments=moments, nu=nu, H=H,
                loop_range=loop_range, forests=forests,
                joint_checks=joint_checks, conditional_checks=conditional_checks)


def cubic_discriminant(a, b, c, d):
    return b*b*c*c - 4*a*c**3 - 4*b**3*d - 27*a*a*d*d + 18*a*b*c*d


class Polynomial:
    """Small sparse exact polynomials used only for one symbolic identity."""

    def __init__(self, value=0):
        terms = value if isinstance(value, dict) else {(): F(value)}
        self.terms = {m: F(c) for m, c in terms.items() if c}

    @staticmethod
    def coerce(value):
        return value if isinstance(value, Polynomial) else Polynomial(value)

    def __add__(self, other):
        other = self.coerce(other)
        result = dict(self.terms)
        for m, c in other.terms.items():
            result[m] = result.get(m, F(0)) + c
        return Polynomial(result)

    __radd__ = __add__

    def __neg__(self):
        return Polynomial({m: -c for m, c in self.terms.items()})

    def __sub__(self, other):
        return self + -self.coerce(other)

    def __rsub__(self, other):
        return self.coerce(other) + -self

    def __mul__(self, other):
        other = self.coerce(other)
        result = defaultdict(F)
        for m, c in self.terms.items():
            for n, d in other.terms.items():
                result[tuple(sorted(m + n))] += c*d
        return Polynomial(dict(result))

    __rmul__ = __mul__

    def __eq__(self, other):
        return self.terms == self.coerce(other).terms


def symbolic_rayleigh_middle_coefficient():
    """Clear three pair determinants in the generic three-site certificate.

    Pair i,j is tested; k is the remaining vertex.  u,v,w are the mutual
    cycle weights on ij,ik,jk; t is the sum of directed triangle weights.
    kappa is K(k,s).  Every quantity in the displayed rational certificate
    is nonnegative for valid loopless kernels and nonnegative jump odds.
    """
    ri, rj, rk, u, v, w, t, kappa = [
        Polynomial({(name,): F(1)}) for name in
        ('ri', 'rj', 'rk', 'u', 'v', 'w', 't', 'kappa')]
    bi, bj, bk = 1+ri, 1+rj, 1+rk
    d = 1-u-v-w-t
    H = bi*bj-1
    Zij = bi*bj-u*H
    Zik = bi*bk-v*(bi*bk-1)
    Zjk = bj*bk-w*(bj*bk-1)
    den_u, den_v, den_w = 1-u, 1-v, 1-w
    common = den_u*den_v*den_w
    # Z_ijk = 1+d[(bi*bj*bk-1)+ri*w/(1-w)
    #                              +rj*v/(1-v)+rk*u/(1-u)].
    Z_cleared = ((1+d*(bi*bj*bk-1))*common
                 + d*(ri*w*den_u*den_v + rj*v*den_u*den_w
                      + rk*u*den_v*den_w))
    left = ((Zjk*bi+Zik*bj-Zij*bk-2*rk*kappa*u*H)*common
            - Z_cleared)
    right = ((2*u+t+rk*(2*u*(1-kappa)+t))*H+rk*t)*common
    right += rk*u*(v+w+t)*den_v*den_w
    right += ri*w*(u+v+t)*den_u*den_v
    right += rj*v*(u+w+t)*den_u*den_w
    assert left == right
    print('PASS: generic three-site Rayleigh middle coefficient, exact polynomial identity')


def three_site_rayleigh_audit(K, p, data):
    """Check every pair's polynomial, with no external-field sampling."""
    assert len(K) == 3
    Q = frozenset(range(3))
    # Removing jump self-loops preserves Z_R and a_i=r_i K(i,s).
    c = [1-K[i][i] for i in range(3)]
    jump = [[F(0) if i == j else K[i][j]/c[i] for j in range(4)]
            for i in range(3)]
    r = [c[i]*(1-p[i])/p[i] for i in range(3)]
    b = [1+x for x in r]
    a = [r[i]*jump[i][3] for i in range(3)]
    Z = data['Z'][Q]
    masses = data['mass']
    polynomials = []
    for i, j in combinations(range(3), 2):
        k = next(iter(Q-{i, j}))
        A = [masses.get(frozenset(), F(0)), masses.get(frozenset((k,)), F(0))]
        B = [masses.get(frozenset((i,)), F(0)), masses.get(frozenset((i, k)), F(0))]
        C = [masses.get(frozenset((j,)), F(0)), masses.get(frozenset((j, k)), F(0))]
        D = [masses.get(frozenset((i, j)), F(0)), masses.get(Q, F(0))]
        coefficients = [B[0]*C[0]-A[0]*D[0],
                        B[0]*C[1]+B[1]*C[0]-A[0]*D[1]-A[1]*D[0],
                        B[1]*C[1]-A[1]*D[1]]
        u = jump[i][j]*jump[j][i]
        v = jump[i][k]*jump[k][i]
        w = jump[j][k]*jump[k][j]
        t = (jump[i][j]*jump[j][k]*jump[k][i]
             + jump[i][k]*jump[k][j]*jump[j][i])
        H = b[i]*b[j]-1
        certificate = ((2*u+t+r[k]*(2*u*(1-jump[k][3])+t))*H
                       + r[k]*(t+u*(v+w+t)/(1-u))
                       + r[i]*w*(u+v+t)/(1-w)
                       + r[j]*v*(u+w+t)/(1-v))
        assert coefficients[1] == a[i]*a[j]*a[k]*certificate/(Z*Z)
        assert coefficients[2] == a[i]*a[j]*a[k]*a[k]*u*H/(Z*Z)
        assert min(coefficients) >= 0
        polynomials.append(coefficients)
    return polynomials


def labelled_cycle_audit():
    """One fixed cycle: original forest enumeration versus the closed PGF.

    Four heterogeneous rows, self-loops and two exit labels.  No truncation,
    field sampling or kernel sweep.  The arbitrary-length proof uses that
    every proper induced subgraph of a simple cycle is acyclic.
    """
    n = 4
    Q = frozenset(range(n))
    sets = subsets(Q)
    # Columns 0,...,3 are roots; columns 4,5 are distinct exit labels.
    K = [[F(1, 5), F(2, 5), F(0), F(0), F(1, 10), F(3, 10)],
         [F(0), F(1, 4), F(1, 2), F(0), F(1, 12), F(1, 6)],
         [F(0), F(0), F(1, 3), F(1, 3), F(1, 6), F(1, 6)],
         [F(1, 2), F(0), F(0), F(1, 6), F(1, 12), F(1, 4)]]
    p = [F(1, 3), F(2, 5), F(3, 7), F(4, 9)]
    assert all(sum(row) == 1 for row in K)
    r = [(1-x)/x for x in p]
    M = [[F(i == j)-K[i][j] for j in range(n)] for i in range(n)]
    d = {C: minor(M, C) for C in sets}
    Z = d[Q]*sum(prod(r[i] for i in Q-C)/d[C] for C in sets)
    nu = {C: d[Q]*prod(r[i] for i in Q-C)/(Z*d[C]) for C in sets}
    assert sum(nu.values()) == 1

    choices = [[j for j in range(n+2) if K[i][j]] for i in range(n)]
    forests = []
    for arrows in product(*choices):
        for start in range(n):
            seen, vertex = set(), start
            while vertex < n and vertex not in seen:
                seen.add(vertex)
                vertex = arrows[vertex]
            if vertex < n:
                break
        else:
            forests.append((arrows, prod(K[i][j] for i, j in enumerate(arrows))))
    assert len(forests) == 80
    assert sum(weight for _, weight in forests) == d[Q]
    measured = defaultdict(F)
    for C, cweight in nu.items():
        for arrows, weight in forests:
            marks = tuple(arrows[i]-n+1 if i not in C and arrows[i] >= n else 0
                          for i in range(n))
            measured[marks] += cweight*weight/d[Q]
    assert sum(measured.values()) == 1

    # Effective parameters after deleting self-loops.  Marker weights
    # r_i K(i,label) are invariant under this transformation.
    c = [1-K[i][i] for i in range(n)]
    effective_p = [1/(1+c[i]*r[i]) for i in range(n)]
    theta = prod(K[i][(i+1) % n]/c[i] for i in range(n))
    beta = theta*(1-prod(effective_p))
    u = [[r[i]*K[i][label]/(1+c[i]*r[i]) for label in (n, n+1)]
         for i in range(n)]
    assert theta == F(1, 10)
    assert beta == F(175453, 1840930)
    assert Z == (1-beta)/prod(effective_p)
    empty_product_mass = prod(1-sum(row) for row in u)
    assert empty_product_mass > beta
    for marks in product(range(3), repeat=n):
        independent = prod(u[i][mark-1] if mark else 1-sum(u[i])
                           for i, mark in enumerate(marks))
        expected = (independent-(beta if not any(marks) else 0))/(1-beta)
        assert measured[marks] == expected
    print('PASS: one heterogeneous four-cycle, self-loops and two exit labels;'
          ' all 81 marker masses match the product-minus-constant PGF')


def main():
    symbolic_rayleigh_middle_coefficient()
    K = [[F(0), F(1, 4), F(1, 4), F(1, 2)],
         [F(1, 4), F(0), F(1, 4), F(1, 2)],
         [F(1, 4), F(1, 4), F(0), F(1, 2)]]
    data = audit(K, [F(1, 2)] * 3)
    Q = frozenset(range(3))
    expected = [F(73, 212), F(33, 212), F(12, 212), F(4, 212)]
    assert all(data['mass'][S] == expected[len(S)] for S in subsets(Q))
    assert data['Z'][Q] == F(53, 8)
    assert data['H'] == [[F(3 if i == j else 1, 5) for j in range(3)]
                         for i in range(3)]
    ranges = {0: F(50, 53), 1: F(0), 2: F(5, 318), 3: F(1, 106)}
    assert all(data['loop_range'][L] == ranges[len(L)] for L in subsets(Q))
    assert cubic_discriminant(4, 36, 99, 73) == -15984
    t = -F(13, 3)
    assert 213 + 104*t + 12*t*t == -F(37, 3)
    print('reversible triangle: G=(73+33 e1+12 e2+4 e3)/212')
    print('count cubic discriminant = -15984: not strongly Rayleigh')
    print('loop ranges: empty=50/53, each pair=5/318, full=1/106')
    print('forest kernel H=(2I+J)/5; conditional thinning verified')
    triangle_rayleigh = three_site_rayleigh_audit(K, [F(1, 2)]*3, data)
    assert triangle_rayleigh == [[F(213, 212**2), F(104, 212**2), F(12, 212**2)]]*3
    print('triangle: every Rayleigh difference is (213+104 z+12 z^2)/212^2')

    directed = [[F(0), F(1, 2), F(0), F(1, 2)],
                [F(0), F(0), F(1, 2), F(1, 2)],
                [F(1, 2), F(0), F(0), F(1, 2)]]
    cycle = audit(directed, [F(1, 2)] * 3)
    assert all(cycle['mass'][S] == F(20 if not S else 3**(3 - len(S)), 57)
               for S in subsets(Q))
    # A real, non-Hermitian correlation kernel.  Its inverse has diagonal 4
    # and cyclic off-diagonal magnitudes (1,1,7), with product 7.
    A = [[F(4), F(-1), F(0)], [F(0), F(4), F(-1)],
         [F(-7), F(0), F(4)]]
    J = inverse(A)
    assert all(minor(J, S) == cycle['moments'][S] for S in subsets(Q))
    print('directed three-cycle: G=(product_i(3+z_i)-7)/57')
    print('directed example has a real nonsymmetric determinantal kernel')
    cycle_rayleigh = three_site_rayleigh_audit(directed, [F(1, 2)]*3, cycle)
    assert cycle_rayleigh == [[F(21, 57**2), F(7, 57**2), F(0)]]*3
    print('directed cycle: every Rayleigh difference is 7(3+z)/57^2')

    half = audit([[F(0), F(1, 2), F(1, 2)],
                  [F(1, 2), F(0), F(1, 2)]], [F(1, 2)] * 2)
    # Given C={1}, E_0=D_0.  The three forests are equiprobable.
    d0 = sum(w for D, w in half['forests'] if 0 in D)
    d1 = sum(w for D, w in half['forests'] if 1 in D)
    d01 = sum(w for D, w in half['forests'] if {0, 1} <= D)
    assert d01 - d0*d1 == -F(1, 9)
    # The single-active first-sleep matrix is (I-diag(1-p)K)^-1 diag(p).
    S11 = F(1, 2) / (1 - F(1, 4)**2)
    assert S11 == F(8, 15)
    assert d0*(1-d0)*(1-S11) == F(14, 135)
    print('given C={1}: Cov(E0,D1)=-1/9; Cov(E0,X1)=14/135')
    total = [data, cycle, half]
    print('PASS:', sum(x['joint_checks'] for x in total), 'joint root/marker',
          'identities;', sum(x['conditional_checks'] for x in total),
          'conditional marker-pair checks; no sampling or sweeps')
    labelled_cycle_audit()


if __name__ == '__main__':
    main()
