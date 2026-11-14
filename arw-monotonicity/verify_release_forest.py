#!/usr/bin/env python3
"""Exact resident-release forest certificates; standard library only.

Run from any directory:
    python3 -B /path/to/verify_release_forest.py
    python3 -B /path/to/verify_release_forest.py --full

The default checks the two-site formulas (including self-loops), two of the
nine prescribed kernels, and an independent raw four-site ARW calculation.
--full checks all nine prescribed kernels, with no parameter search.
All arithmetic is rational; no files are read or written.

Notation: P0 is binary stabilization with p_i=1, W is resident release,
and L=I-W. In i-vacant/i-present order,

    W = [[I,0],[B,Q]],  A=I-Q,  T=A^{-1}B.
    D(s)=det(I+sA)=sum_k d_k s^k,
    F_k=[s^k]adj(I+sA),  H_k=[[d_k I,0],[F_{k-1}B,F_k]].

Here F_{-1}=0, F_0=I, F_m=0, m=2^(n-1). The exact recurrence is
d_k=tr(A F_{k-1})/k, F_k=d_k I-A F_{k-1}. The scored quantity is
(delta_{C+x}-delta_C) P0 H_k 1_v, with i in C and x outside C.
In particular, the vacant-output contribution must not be discarded.

The endpoint audits check W A_i=A_i and W^infinity=A_i V_i, including
H_m=d_m A_i V_i. Here A_i is physical addition and V_i deletes i.
They also check P_c T P_c^{-1} equals ordinary IDLA addition for the
walk exiting i on the trace graph, where P_c is the occupied block of P0.
These identities do not imply A and T commute.

The generic two-site table is documented in two_site_prediction().
Finite certificates do not prove the general cyclic coefficient conjecture.
AI-written research code from Jacob Richey's mathematical specifications.
"""

import argparse
from collections import defaultdict
from dataclasses import dataclass
from fractions import Fraction as F
from itertools import product


def require(condition, message):
    """Keep scientific checks enabled even under python -O."""
    if not condition:
        raise AssertionError(message)


def identity(n):
    return [[F(a == b) for b in range(n)] for a in range(n)]


def zeros(n, width=None):
    return [[F(0)] * (n if width is None else width) for _ in range(n)]


def multiply(a, b):
    columns = list(zip(*b))
    return [[sum((x * y for x, y in zip(row, col) if x and y), F(0))
             for col in columns] for row in a]


def subtract(a, b):
    return [[x - y for x, y in zip(ar, br)] for ar, br in zip(a, b)]


def scale(c, a):
    return [[c * x for x in row] for row in a]


def inverse(matrix):
    n = len(matrix)
    rows = [list(row) + unit for row, unit in zip(matrix, identity(n))]
    for i in range(n):
        pivot = next((j for j in range(i, n) if rows[j][i]), None)
        require(pivot is not None, "singular rational matrix")
        rows[i], rows[pivot] = rows[pivot], rows[i]
        denominator = rows[i][i]
        rows[i] = [x / denominator for x in rows[i]]
        for j in range(n):
            if j != i and rows[j][i]:
                coefficient = rows[j][i]
                rows[j] = [x - coefficient * y
                           for x, y in zip(rows[j], rows[i])]
    return [row[n:] for row in rows]


def block(matrix, rows, cols):
    return [[matrix[a][b] for b in cols] for a in rows]


def stochastic(matrix):
    return all(all(q >= 0 for q in row) and sum(row) == 1
               for row in matrix)


def supports(n):
    return tuple(frozenset(i for i, bit in enumerate(bits) if bit)
                 for bits in product((0, 1), repeat=n))


class AdditionAlgebra:
    """Ordinary IDLA matrices and the complete-block inverse E=P^{-1}.

    D_j is IDLA addition in the binary active-input basis, J_j is its
    forced-jump mixture. Row C of E is
    delta_empty product_{j in C} ((1+r_j)D_j-r_j J_j).
    The physical sleeping-support addition is E D_j P.
    """

    def __init__(self, kernel):
        self.K = [[F(x) for x in row] for row in kernel]
        self.n = n = len(self.K)
        require(n > 0 and all(len(row) == n + 1 for row in self.K),
                "kernel must have one final sink column")
        require(stochastic(self.K), "kernel must be stochastic")
        reachable = {i for i in range(n) if self.K[i][n] > 0}
        while True:
            enlarged = reachable | {
                i for i in range(n)
                if any(self.K[i][j] > 0 for j in reachable)}
            if enlarged == reachable:
                break
            reachable = enlarged
        require(len(reachable) == n, "kernel must be sink-accessible")
        self.sets = sets = supports(n)
        self.pos = pos = {state: a for a, state in enumerate(sets)}
        size = len(sets)
        self.Id = identity(size)
        self.D = [zeros(size) for _ in range(n)]
        for a, state in enumerate(sets):
            occupied = sorted(state)
            green = inverse([[F(x == y) - self.K[x][y]
                              for y in occupied] for x in occupied])
            for x in range(n):
                if x not in state:
                    self.D[x][a][pos[state | {x}]] = F(1)
                    continue
                row = green[occupied.index(x)]
                self.D[x][a][a] = sum(
                    (q * self.K[z][n] for z, q in zip(occupied, row)), F(0))
                for y in range(n):
                    if y not in state:
                        self.D[x][a][pos[state | {y}]] = sum(
                            (q * self.K[z][y]
                             for z, q in zip(occupied, row)), F(0))
        require(all(stochastic(matrix) for matrix in self.D),
                "IDLA addition is not stochastic")
        self.J = [
            [[self.K[i][n] * self.Id[a][b] + sum(
                (self.K[i][j] * self.D[j][a][b] for j in range(n)), F(0))
              for b in range(size)] for a in range(size)]
            for i in range(n)]

    def inverse_basis(self, odds):
        require(len(odds) == self.n and all(r >= 0 for r in odds),
                "odds must be nonnegative")
        factors = [subtract(scale(1 + odds[j], self.D[j]),
                            scale(odds[j], self.J[j]))
                   for j in range(self.n)]
        result = []
        for state in self.sets:
            row = [F(not other) for other in self.sets]
            for j in sorted(state):
                row = multiply([row], factors[j])[0]
            result.append(row)
        return result


class ReleaseForest:
    """Exact release data and all finite forest-degree coefficients."""

    def __init__(self, algebra, odds, root, degrees=True):
        self.algebra, self.root = algebra, root
        self.odds = [F(x) for x in odds]
        require(0 <= root < algebra.n, "invalid release site")
        self.odds[root] = F(0)
        other = list(self.odds)
        other[root] = F(1)
        self.E0 = algebra.inverse_basis(self.odds)
        self.P0 = inverse(self.E0)
        self.L = multiply(subtract(algebra.inverse_basis(other), self.E0),
                          self.P0)
        self.W = subtract(algebra.Id, self.L)
        require(stochastic(self.P0) and stochastic(self.W),
                "P0 or W is not stochastic")
        self.vacant = [a for a, state in enumerate(algebra.sets)
                       if root not in state]
        self.present = [a for a, state in enumerate(algebra.sets)
                        if root in state]
        require(all(self.W[a] == algebra.Id[a] for a in self.vacant),
                "release must fix vacant supports")
        require(not any(self.P0[a][b] for a in self.present for b in self.vacant),
                "common root must remain occupied at p_i=1")
        self.Q = block(self.W, self.present, self.present)
        self.B = block(self.W, self.present, self.vacant)
        self.m = len(self.present)
        self.A = subtract(identity(self.m), self.Q)
        self.T = multiply(inverse(self.A), self.B)
        require(stochastic(self.T), "absorbing kernel is not stochastic")
        self.audit_endpoints()
        self.d, self.F = [F(1)], [identity(self.m)]
        if degrees:
            self.build_degrees()

    def audit_endpoints(self):
        """Audit W A_i=A_i, the absorbing projection, and trace addition."""
        alg, i = self.algebra, self.root
        Ai = multiply(multiply(self.E0, alg.D[i]), self.P0)
        Ji = multiply(multiply(self.E0, alg.J[i]), self.P0)
        require(stochastic(Ai) and stochastic(Ji), "physical additions")
        require(multiply(self.W, Ai) == Ai, "W A_i != A_i")
        require(multiply(Ji, Ai) == multiply(Ai, Ji), "addition commutation")
        deletion = zeros(len(alg.sets))
        for a, state in enumerate(alg.sets):
            deletion[a][alg.pos[state - {i}]] = F(1)
        physical_release = multiply(deletion, Ji)
        require(all(physical_release[a] == self.W[a] for a in self.present),
                "derivative release differs from remove + forced jump")
        self.projection = multiply(Ai, deletion)
        require(multiply(self.W, self.projection) == self.projection,
                "A_i V_i is not W-harmonic")
        require(multiply(self.projection, self.W) == self.projection,
                "A_i V_i is not supported on absorbing states")
        require(all(self.projection[a] == alg.Id[a] for a in self.vacant),
                "projection does not fix absorbing states")
        require(block(self.projection, self.present, self.vacant) == self.T
                and not any(self.projection[a][b]
                            for a in self.present for b in self.present),
                "A_i V_i differs from [[I,0],[A^{-1}B,0]]")

        # At a permanently occupied i, arrivals follow the first exit from i.
        remaining = [j for j in range(alg.n) if j != i]
        exit_sites = remaining + [alg.n]
        exit_law = [alg.K[i][j] / (1 - alg.K[i][i]) for j in exit_sites]
        trace_K = [[alg.K[j][k] + alg.K[j][i] * weight
                    for k, weight in zip(exit_sites, exit_law)]
                   for j in remaining]
        if not remaining:
            return
        trace = AdditionAlgebra(trace_K)
        trace_E = trace.inverse_basis([self.odds[j] for j in remaining])
        trace_P = inverse(trace_E)
        Pc = block(self.P0, self.present, self.present)
        require(Pc == trace_P, "occupied P0 block is not trace stabilization")
        trace_D = [
            [exit_law[-1] * trace.Id[a][b] + sum(
                (exit_law[j] * trace.D[j][a][b]
                 for j in range(len(remaining))), F(0))
             for b in range(self.m)] for a in range(self.m)]
        require(multiply(multiply(Pc, self.T), trace_E) == trace_D,
                "P_c T P_c^{-1} is not trace IDLA exit addition")

    def build_degrees(self):
        unit = identity(self.m)
        for k in range(1, self.m + 1):
            AF = multiply(self.A, self.F[-1])
            dk = sum((AF[a][a] for a in range(self.m)), F(0)) / k
            Fk = subtract(scale(dk, unit), AF)
            require(dk > 0 and all(x >= 0 for row in Fk for x in row),
                    f"invalid forest coefficient k={k}")
            self.d.append(dk)
            self.F.append(Fk)
            # B1=A1, so this verifies row normalization of the full H_k.
            Fprev_B1 = multiply(self.F[k - 1],
                               [[sum(row)] for row in self.B])
            require(all(sum(row) + Fprev_B1[a][0] == dk
                        for a, row in enumerate(Fk)), f"H_{k} row sums")
        require(not any(x for row in self.F[-1] for x in row), "F_m != 0")
        require(self.H(self.m) == scale(self.d[-1], self.projection),
                "H_m != d_m A_i V_i")

    def H(self, k):
        """Full coefficient in original support order, not just F_k."""
        require(0 <= k < len(self.d), "coefficient not constructed")
        result = zeros(len(self.algebra.sets))
        for a in self.vacant:
            result[a][a] = self.d[k]
        for a, row in zip(self.present, self.F[k]):
            for b, value in zip(self.present, row):
                result[a][b] = value
        if k:
            bottom_left = multiply(self.F[k - 1], self.B)
            for a, row in zip(self.present, bottom_left):
                for b, value in zip(self.vacant, row):
                    result[a][b] = value
        return result

    def target_coefficients(self):
        """All P0 H_k 1_v columns, retaining the vacant-block term."""
        fv0 = [[F(v in self.algebra.sets[a]) for v in range(self.algebra.n)]
               for a in self.vacant]
        fv1 = [[F(v in self.algebra.sets[a]) for v in range(self.algebra.n)]
               for a in self.present]
        Bf = multiply(self.B, fv0)
        result = []
        for k, dk in enumerate(self.d):
            occupied = multiply(self.F[k], fv1)
            if k:
                correction = multiply(self.F[k - 1], Bf)
                occupied = [[x + y for x, y in zip(ar, br)]
                            for ar, br in zip(occupied, correction)]
            full = zeros(len(self.algebra.sets), self.algebra.n)
            for a, row in zip(self.vacant, fv0):
                full[a] = [dk * x for x in row]
            for a, row in zip(self.present, occupied):
                full[a] = row
            result.append(multiply(self.P0, full))
        return result


def two_site_prediction(kernel, p_other):
    """Generic two-site formulas, with common root 0 and added site 1.

    After self-loop removal, routing probabilities are a=K01/(1-K00),
    b=K10/(1-K11), and effective sleep t=p/[p+(1-p)(1-K11)].
    Write q=1-t, alpha=t/(1-abq), e=1-abq, f=1-ab alpha.
    The physical row difference is alpha(delta_{01}-delta_0).

          k          target 0                   target 1
          0             0                        alpha
          1       alpha ab(1-bq)             alpha(2-a-ab)
          2             0                   alpha(1-a)e f

    All displayed factors are nonnegative for a,b in [0,1], ab<1,
    0<t<1. A root self-loop multiplies degree k by (1-K00)^k.
    The returned W is obtained by first-jump decomposition, independently
    of the inverse-basis derivative used in ReleaseForest.
    """
    K = [[F(x) for x in row] for row in kernel]
    require(len(K) == 2 and all(len(row) == 3 for row in K)
            and stochastic(K), "two-site stochastic kernel required")
    p = F(p_other)
    require(0 < p < 1 and K[0][0] < 1 and K[1][1] < 1,
            "invalid sleep or self-loop")
    c = 1 - K[0][0]
    a, b = K[0][1] / c, K[1][0] / (1 - K[1][1])
    require(a * b < 1, "two-site kernel must be accessible")
    t = p / (p + (1 - p) * (1 - K[1][1]))
    q = 1 - t
    alpha = t / (1 - a * b * q)
    e, f = 1 - a * b * q, 1 - a * b * alpha
    table = [[F(0), alpha],
             [alpha * a * b * (1 - b * q), alpha * (2 - a - a * b)],
             [F(0), alpha * (1 - a) * e * f]]
    table = [[c**k * x for x in row] for k, row in enumerate(table)]
    require(all(x >= 0 for row in table for x in row), "two-site table sign")
    # Support order: empty, {1}, {0}, {0,1}.
    W = identity(4)
    W[2] = [1 - a * b * q - a * t, a * t, a * b * q, F(0)]
    W[3] = [a * (1 - b)**2 * q, 1 - a + a * (1 - b) * t,
            a * b * (1 - alpha + q * (1 - b)), a * b * alpha]
    W = [[c * x + (1 - c) * F(r == col)
          for col, x in enumerate(row)] for r, row in enumerate(W)]
    return table, [F(1), c * (e + f), c * c * e * f], W


def audit_two_site():
    fixtures = [
        ("loop-free",
         [[F(0), F(2, 5), F(3, 5)], [F(3, 4), F(0), F(1, 4)]],
         F(5, 8),
         [[F(0), F(50, 71)], [F(345, 2272), F(65, 71)],
          [F(0), F(21, 71)]]),
        ("self-loops",
         [[F(1, 5), F(2, 5), F(2, 5)], [F(1, 2), F(1, 3), F(1, 6)]],
         F(2, 5),
         [[F(0), F(8, 13)], [F(3, 26), F(36, 65)],
          [F(0), F(8, 65)]]),
    ]
    for name, K, p, expected in fixtures:
        alg = AdditionAlgebra(K)
        data = ReleaseForest(alg, [F(0), (1 - p) / p], 0)
        table, polynomial, W = two_site_prediction(K, p)
        require(table == expected, f"{name}: changed two-site certificate")
        require(data.W == W and data.d == polynomial,
                f"{name}: two-site kernel or determinant formula")
        coeff = data.target_coefficients()
        low, high = alg.pos[frozenset({0})], alg.pos[frozenset({0, 1})]
        actual = [[matrix[high][v] - matrix[low][v] for v in range(2)]
                  for matrix in coeff]
        require(actual == table, f"{name}: two-site coefficient formula")
        print(f"PASS two-site {name}: 6 coefficients and full release formula")


@dataclass(frozen=True)
class Case:
    name: str
    kernel: list
    odds: tuple
    roots: tuple
    counts: tuple  # coefficients, positive, zero


def prescribed_cases():
    """The nine fixed certificates from the report; no generated search."""
    cases = [
        Case("single_leak_two",
             [[0, 1, 0], [F(1, 2), 0, F(1, 2)]],
             (F(1, 7), F(7)), (0, 1), (12, 7, 5)),
        Case("single_leak_three",
             [[0, 1, 0, 0], [0, 0, 1, 0], [F(2, 3), 0, 0, F(1, 3)]],
             (F(2), F(1, 10), F(10)), (0, 1, 2), (180, 150, 30)),
        Case("merge_feedback_three",
             [[0, F(1, 2), F(1, 2), 0], [0, 0, 1, 0],
              [F(1, 2), 0, 0, F(1, 2)]],
             (F(1, 10), F(1), F(10)), (0, 1, 2), (180, 152, 28)),
        Case("bidirected_four", bidirected_cycle(4),
             (F(1, 10), F(10), F(1, 3), F(3)), (0, 1, 2, 3),
             (1728, 1632, 96)),
        Case("merge_return_four",
             [[0, 1, 0, 0, 0], [0, 0, F(1, 2), F(1, 2), 0],
              [F(1, 2), 0, 0, 0, F(1, 2)], [0, 0, F(1, 2), 0, F(1, 2)]],
             (F(10), F(1, 10), F(5), F(1, 5)), (0, 1, 2, 3),
             (1728, 1629, 99)),
        Case("complete_asymmetric_four",
             [[F(x, 10) for x in row] for row in
              ((0, 1, 2, 6, 1), (4, 0, 1, 2, 3),
               (1, 5, 0, 3, 1), (2, 1, 1, 0, 6))],
             (F(1, 20), F(20), F(1, 5), F(5)), (0, 1, 2, 3),
             (1728, 1632, 96)),
    ]
    K = zeros(5, 6)
    for i in range(4):
        K[i][i + 1] = F(1)
    K[4][0], K[4][5] = F(3, 4), F(1, 4)
    cases.append(Case("single_leak_five", K,
                      (F(1), F(1, 5), F(5), F(1, 10), F(10)),
                      (0, 2), (5440, 5288, 152)))
    K = zeros(5, 6)
    K[0][1] = F(1)
    K[1][2] = K[1][3] = F(1, 2)
    K[2][4] = K[3][4] = F(1)
    K[4][0] = K[4][5] = F(1, 2)
    cases.append(Case("five_feedback_diamond", K,
                      (F(1, 10), F(10), F(1, 10), F(10), F(1)),
                      (0, 2), (5440, 5300, 140)))
    cases.append(Case("bidirected_five", bidirected_cycle(5),
                      (F(1, 100), F(100), F(1, 10), F(10), F(1)),
                      (0, 4), (5440, 5312, 128)))
    return cases


def bidirected_cycle(n):
    require(n >= 3, "two distinct neighbours required")
    return [[F(1, 3) if j in ((i + 1) % n, (i - 1) % n, n) else F(0)
             for j in range(n + 1)] for i in range(n)]


def audit_case(case):
    alg = AdditionAlgebra(case.kernel)
    checked = positive = zero = 0
    for root in case.roots:
        data = ReleaseForest(alg, case.odds, root)
        coeff = data.target_coefficients()
        for state in alg.sets:
            if root not in state:
                continue
            for x in range(alg.n):
                if x in state:
                    continue
                low, high = alg.pos[state], alg.pos[state | {x}]
                for k, matrix in enumerate(coeff):
                    for v in range(alg.n):
                        value = matrix[high][v] - matrix[low][v]
                        require(value >= 0,
                                f"{case.name}: root={root}, C={sorted(state)}, "
                                f"x={x}, v={v}, k={k}, coefficient={value}")
                        checked += 1
                        positive += int(value > 0)
                        zero += int(value == 0)
    counts = checked, positive, zero
    require(counts == case.counts, f"{case.name}: certificate counts {counts}")
    print(f"PASS {case.name}: {checked} coefficients "
          f"({positive} positive, {zero} zero); "
          f"{len(case.roots)} endpoint/trace audits")
    return counts


def raw_step(K, p, state):
    """Highest-active-site ARW; skip failed sleeps during collisions.

    Counts -1 denote sleeping residents, 0 vacancies, positive counts active
    particles. At count >=2 use the next jump directly, also when p_i=1.
    This is the continuous p_i=1 boundary used by the release identity.
    """
    n = len(K)
    site = max(i for i, count in enumerate(state) if count > 0)
    result = defaultdict(F)
    jump_weight = F(1)
    if state[site] == 1:
        asleep = list(state)
        asleep[site] = -1
        result[tuple(asleep)] += p[site]
        jump_weight -= p[site]
    for destination, probability in enumerate(K[site]):
        if not probability or not jump_weight:
            continue
        after = list(state)
        after[site] -= 1
        if destination < n:
            after[destination] = (2 if after[destination] == -1
                                  else after[destination] + 1)
        result[tuple(after)] += jump_weight * probability
    return dict(result)


def solve_transient(rows, rhs):
    """Sparse exact elimination of a transient I-P, then back substitution."""
    rows, rhs = [dict(row) for row in rows], [list(row) for row in rhs]
    size = len(rows)
    width = len(rhs[0]) if size else 0
    for i in range(size):
        pivot = rows[i].get(i, F(0))
        require(pivot > 0, "raw transient elimination has nonpositive pivot")
        tail = {j: value / pivot for j, value in rows[i].items()
                if j > i and value}
        rows[i] = tail
        rhs[i] = [x / pivot for x in rhs[i]]
        for j in range(i + 1, size):
            coefficient = rows[j].pop(i, F(0))
            if not coefficient:
                continue
            for k, value in tail.items():
                updated = rows[j].get(k, F(0)) - coefficient * value
                if updated:
                    rows[j][k] = updated
                else:
                    rows[j].pop(k, None)
            rhs[j] = [x - coefficient * y for x, y in zip(rhs[j], rhs[i])]
    solution = zeros(size, width)
    for i in range(size - 1, -1, -1):
        solution[i] = [
            rhs[i][v] - sum((q * solution[j][v] for j, q in rows[i].items()),
                            F(0)) for v in range(width)]
    return solution


def raw_absorption(initial, K, p, sets):
    """Explore raw configurations without using the addition algebra."""
    initial = list(dict.fromkeys(initial))
    active = lambda state: any(x > 0 for x in state)
    states = [state for state in initial if active(state)]
    index = {state: j for j, state in enumerate(states)}
    pos = {state: j for j, state in enumerate(sets)}
    rows, rhs = [], []
    for state in states:
        row = {index[state]: F(1)}
        out = [F(0)] * len(sets)
        transitions = raw_step(K, p, state)
        require(sum(transitions.values()) == 1, "raw transition normalization")
        for after, weight in transitions.items():
            if not weight:
                continue
            if active(after):
                if after not in index:
                    index[after] = len(states)
                    states.append(after)
                j = index[after]
                row[j] = row.get(j, F(0)) - weight
            else:
                support = frozenset(j for j, count in enumerate(after)
                                    if count == -1)
                out[pos[support]] += weight
        rows.append(row)
        rhs.append(out)
    solution = solve_transient(rows, rhs)
    require(stochastic(solution), "raw absorption laws are not stochastic")
    laws = {}
    for state in initial:
        if active(state):
            laws[state] = solution[index[state]]
        else:
            support = frozenset(j for j, count in enumerate(state) if count == -1)
            laws[state] = [F(other == support) for other in sets]
    return laws, len(states)


def audit_raw_four():
    """Prescribed bidirected four-cycle, root 0, other sleeps 1/2."""
    alg = AdditionAlgebra(bidirected_cycle(4))
    data = ReleaseForest(alg, (F(0), F(1), F(1), F(1)), 0, degrees=False)
    active = [tuple(int(j in state) for j in range(4)) for state in alg.sets]
    release_starts = {}
    for a in data.present:
        state = alg.sets[a]
        law = defaultdict(F)
        for destination, probability in enumerate(alg.K[0]):
            if not probability:
                continue
            after = [-int(j in state and j != 0) for j in range(4)]
            if destination < 4:
                after[destination] = 2 if after[destination] == -1 else 1
            law[tuple(after)] += probability
        release_starts[a] = dict(law)
    starts = active + [state for law in release_starts.values() for state in law]
    laws, raw_count = raw_absorption(
        starts, alg.K, [F(1), F(1, 2), F(1, 2), F(1, 2)], alg.sets)
    require([laws[state] for state in active] == data.P0,
            "raw four-site P0 mismatch")
    for a, law in release_starts.items():
        row = [sum((weight * laws[state][b] for state, weight in law.items()),
                   F(0)) for b in range(len(alg.sets))]
        require(row == data.W[a], f"raw four-site release row {a} mismatch")
    require(raw_count == 127, f"raw transient-state count changed: {raw_count}")
    print("PASS raw four-site ARW: 127 transient states, "
          "256 P0 entries, 128 nontrivial W entries")


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--full", action="store_true",
                        help="run all nine prescribed cases (no parameter search)")
    args = parser.parse_args()
    audit_two_site()
    cases = prescribed_cases()
    if not args.full:
        cases = [cases[0], cases[3]]
    totals = [0, 0, 0]
    for case in cases:
        for j, value in enumerate(audit_case(case)):
            totals[j] += value
    expected = (21876, 21102, 774) if args.full else (1740, 1639, 101)
    require(tuple(totals) == expected, f"unexpected overall totals {totals}")
    audit_raw_four()
    print(f"PASS {'full' if args.full else 'focused'} suite: "
          f"{len(cases)} prescribed cases; {totals[0]} coefficients "
          f"({totals[1]} positive, {totals[2]} zero).")
    print("General cyclic coefficient positivity remains OPEN.")


if __name__ == "__main__":
    main()
