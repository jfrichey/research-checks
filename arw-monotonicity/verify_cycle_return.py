#!/usr/bin/env python3
"""Exact killed-return certificate for every directed three-cycle.

Run from any directory with Python and SymPy:
    python3 -B /path/to/verify_cycle_return.py

Verified with SymPy 1.14.0. The program reads no data files, writes no
files, and imports no project or scratch modules. All calculations are
symbolic or rational; no parameter grid or floating-point tolerance is
used. Checks remain enabled under python -O.

Vertices are 0 -> 1 -> 2 -> 0, with forward probabilities (a,b,c) and
the remaining probabilities going to the sink. Assume 0<=a,b,c<=1,
abc<1, p_0=1, and (p_1,p_2)=(p,q) in (0,1)^2. Failed sleeps at piles
are skipped, including the continuous-extension convention at p_0=1.

Q(B,C) is the probability that releasing the resident at 0 from sleeping
support {0} union B leaves final support {0} union C. It is the occupied
block of RN J53/J113's resident-release kernel, equivalently the old N3
killed first-chance kernel. Missing mass is killed.

The theorem certified here is Qg(B)<=Qg(B') whenever B subset B' and
g is nonnegative and increasing. No entrance-law order or all-length
theorem is asserted.

Why six polynomials suffice:
  * Particle conservation gives Q(B,C)=0 for |C|>|B|.
  * Q1(B) is a times the probability of a positive first path return.
    Common output coins at busy buffers and ordered routing trials
    order this returned count when B grows. Thus survival is increasing.
  * The empty row is supported at empty. The only remaining comparisons
    are B=10,01 against 11, for the upper sets U1, U2, U+ (containing
    site 1, containing site 2, or nonempty).

Put k=abc, E=1-k(1-p)(1-q), F=1-k+kpq. E,F>0 on the domain; k=0
gives Q=0. For those six comparisons the polynomials
    H(B,U) = E F [Q(11,U)-Q(B,U)] / (k w(U)),
    w(U1)=p, w(U2)=q, w(U+)=1
have nonnegative Bernstein coefficients on [0,1]^5. The exact expected
degrees/counts/minima are below. There are 2,232 coefficients: 2,034
positive and 198 zero. Every polynomial is reconstructed from its full
Bernstein expansion and checked identically.

The raw one-instruction calculation is an independent finite validation
of Q's construction, not the parameter-uniform proof. A final exact
check preserves the common-first-pass cylinder obstruction: a cylinder
of probability 1/32 stops the lower process at {1}, but leaves the upper
at queue 2 and background {1}; its continuation retains 1 with probability
107/195. Thus conditional same-cylinder coupling fails despite the
whole-row theorem.

AI-written research code from Jacob Richey's mathematical specifications.
"""

import sys

sys.dont_write_bytecode = True

from collections import defaultdict, deque
from fractions import Fraction as F
from itertools import product
from math import comb

import sympy as S


VARIABLES = S.symbols("a b c p q", nonnegative=True)
a, b, c, p, q = VARIABLES
BITS = ((0, 0), (1, 0), (0, 1), (1, 1))
UPPER_SETS = ((1, 3), (2, 3), (1, 2, 3))
LOW = tuple((1, bits) for bits in BITS)
HIGH = ((2, (0, 0)), (2, (1, 0)), (2, (0, 1)), (3, (0, 0)))
STATES = LOW + HIGH
POSITION = {state: i for i, state in enumerate(STATES)}
# (multidegree in a,b,c,p,q, positive count, zero count, minimum positive)
EXPECTED = (
    ((2, 2, 3, 2, 3), 404, 28, S.Rational(1, 9)),
    ((1, 2, 2, 2, 1), 98, 10, S.Rational(1, 4)),
    ((2, 2, 3, 3, 3), 495, 81, S.Rational(1, 27)),
    ((2, 2, 3, 2, 3), 424, 8, S.Rational(1, 6)),
    ((1, 2, 2, 2, 1), 98, 10, S.Rational(1, 4)),
    ((2, 2, 3, 3, 3), 515, 61, S.Rational(1, 18)),
)


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def path_pass(residents, incoming, sleeps, forwards):
    """Exact cut-path law of (terminal bits, returned batch)."""
    row = {(incoming, ()): S.Integer(1)}
    for old, sleep, forward in zip(residents, sleeps, forwards):
        out = defaultdict(lambda: S.Integer(0))
        for (queue, bits), weight in row.items():
            coins = ((old, S.Integer(1)),) if queue == 0 else (
                (0, 1 - sleep), (1, sleep))
            for last, probability in coins:
                emitted = queue + old - last
                for returned in range(emitted + 1):
                    routing = (comb(emitted, returned) * forward**returned
                               * (1 - forward)**(emitted - returned))
                    out[returned, bits + (last,)] += weight * probability * routing
        row = dict(out)
    result = {(bits, returned): S.factor(weight)
              for (returned, bits), weight in row.items() if weight != 0}
    require(S.factor(sum(result.values()) - 1) == 0, "path normalization")
    return result


def symbolic_return():
    """Solve the four high-state equations; no time truncation."""
    transition = S.zeros(8)
    passes = {bits: path_pass(bits, 1, (p, q), (b, c)) for bits in BITS}
    for i, (queue, old) in enumerate(STATES):
        if queue > 1:
            transition[i, POSITION[queue - 1, old]] += 1 - a
        for (new, returned), probability in passes[old].items():
            next_queue = queue - 1 + returned
            if next_queue > 0:
                transition[i, POSITION[next_queue, new]] += a * probability
    transition = transition.applyfunc(S.factor)
    ll, lh = transition[:4, :4], transition[:4, 4:]
    hl, hh = transition[4:, :4], transition[4:, 4:]
    descent = ((S.eye(4) - hh).inv(method="DM") * hl).applyfunc(S.factor)
    require(((S.eye(4) - hh) * descent - hl).applyfunc(S.factor) == S.zeros(4),
            "high-state equations")
    require(all(S.factor(sum(descent.row(i)) - 1) == 0 for i in range(4)),
            "high-state absorption normalization")
    return (ll + lh * descent).applyfunc(S.factor), descent


def bernstein_certificate(expression):
    """Compute every exact coefficient and reconstruct the polynomial."""
    polynomial = S.Poly(expression, *VARIABLES)
    degrees = polynomial.degree_list()
    coefficients = []
    for index in product(*(range(degree + 1) for degree in degrees)):
        value = sum(
            coefficient * S.prod(
                S.binomial(i, j) / S.binomial(degree, j)
                for i, j, degree in zip(index, monomial, degrees))
            for monomial, coefficient in polynomial.terms()
            if all(i >= j for i, j in zip(index, monomial))
        )
        require(value >= 0, f"negative Bernstein coefficient {index}: {value}")
        coefficients.append((index, value))
    reconstructed = S.Poly(0, *VARIABLES)
    for index, value in coefficients:
        if value:
            term = value * S.prod(
                S.binomial(degree, i) * x**i * (1 - x)**(degree - i)
                for x, degree, i in zip(VARIABLES, degrees, index))
            reconstructed += S.Poly(term, *VARIABLES)
    require(reconstructed == polynomial, "reconstructed Bernstein identity")
    positive = [value for _, value in coefficients if value]
    return (degrees, len(positive), len(coefficients) - len(positive),
            min(positive))


def fraction_inverse(matrix):
    """Independent rational Gaussian elimination for the raw ARW audit."""
    size = len(matrix)
    rows = [list(row) + [F(i == j) for j in range(size)]
            for i, row in enumerate(matrix)]
    for i in range(size):
        pivot = next((j for j in range(i, size) if rows[j][i]), None)
        require(pivot is not None, "singular raw transient matrix")
        rows[i], rows[pivot] = rows[pivot], rows[i]
        denominator = rows[i][i]
        rows[i] = [x / denominator for x in rows[i]]
        for j in range(size):
            if j != i and rows[j][i]:
                weight = rows[j][i]
                rows[j] = [x - weight * y for x, y in zip(rows[j], rows[i])]
    return [row[size:] for row in rows]


def raw_return(parameters):
    """Raw legal topplings: -1 sleeper, 0 vacant, positive active count."""
    aa, bb, cc, pp, qq = map(F, parameters)
    forwards, sleeps = (aa, bb, cc), (F(1), pp, qq)
    # Conditional on the forced first jump reaching site 1.
    starts = [(0, 2 if x else 1, -y) for x, y in BITS]

    def step(state):
        vertex = next(i for i, count in enumerate(state) if count > 0)
        out = defaultdict(F)
        # Skip ineffective sleeps at piles, including the p_0=1 boundary.
        chance = sleeps[vertex] if state[vertex] == 1 else F(0)
        if chance:
            changed = list(state)
            changed[vertex] = -1
            out[tuple(changed)] += chance
        routes = (((vertex + 1) % 3, forwards[vertex]),
                  (3, 1 - forwards[vertex]))
        for destination, probability in routes:
            weight = (1 - chance) * probability
            if not weight:
                continue
            changed = list(state)
            changed[vertex] -= 1
            if destination < 3:
                changed[destination] = (2 if changed[destination] == -1
                                        else changed[destination] + 1)
            out[tuple(changed)] += weight
        require(sum(out.values()) == 1, "raw transition normalization")
        return out

    transient, seen, rows = [], set(), {}
    pending = deque(starts)
    while pending:
        state = pending.popleft()
        if state in seen:
            continue
        seen.add(state)
        if not any(count > 0 for count in state):
            continue
        transient.append(state)
        rows[state] = step(state)
        pending.extend(rows[state])
    index = {state: i for i, state in enumerate(transient)}
    size = len(transient)
    matrix = [[F(i == j) for j in range(size)] for i in range(size)]
    reward = [[F(0) for _ in BITS] for _ in transient]
    for i, state in enumerate(transient):
        for end, weight in rows[state].items():
            if end in index:
                matrix[i][index[end]] -= weight
            elif end[0] == -1:
                bits = (int(end[1] == -1), int(end[2] == -1))
                reward[i][BITS.index(bits)] += weight
    inverse = fraction_inverse(matrix)
    result = [[aa * sum(inverse[index[state]][j] * reward[j][k]
                        for j in range(size)) for k in range(4)]
              for state in starts]
    return result, size


def cylinder_obstruction(descent):
    """Check an actual common first-pass cylinder, not arbitrary queues."""
    half = S.Rational(1, 2)
    parameters = dict(zip(VARIABLES, (1, half, half, half, half)))
    # Common outputs (1,0), forwarded trial at 1, two forwarded trials at 2.
    probability = (p * b * (1 - q) * c**2).subs(parameters)
    require(probability == S.Rational(1, 32), "cylinder probability")
    returns = []
    for old in ((1, 0), (1, 1)):
        at_second = 1 + old[0] - 1  # output 1 at site 1; all arrows forward
        returned = at_second + old[1]  # output 0 at site 2
        returns.append(returned)
    require(returns == [1, 2], "cylinder returned batches")
    high = HIGH.index((2, (1, 0)))
    retained = S.factor(descent[high, 1] + descent[high, 3]).subs(parameters)
    require(retained == S.Rational(107, 195), "conditional continuation")
    print(f"cylinder probability={probability}; upper continuation retains "
          f"site 1 with probability={retained} < 1")


def main():
    print(f"SymPy {S.__version__}; certificate developed with 1.14.0", flush=True)
    kernel, descent = symbolic_return()
    kappa = a * b * c
    e = 1 - kappa * (1 - p) * (1 - q)
    f = 1 - kappa + kappa * p * q
    require(S.factor(kernel[0, 0] - kappa * (1 - p) * (1 - q)) == 0,
            "empty diagonal formula")
    require(all(kernel[0, j] == 0 for j in (1, 2, 3)), "empty-row support")
    require(kernel[1, 3] == kernel[2, 3] == 0, "particle-count support")
    require(S.factor(kernel[3, 3] - kappa * p * q / f) == 0,
            "full diagonal formula")

    polynomials, summaries = [], []
    for lower in (1, 2):
        for upper_set in UPPER_SETS:
            prefactor = p if upper_set == (1, 3) else (
                q if upper_set == (2, 3) else S.Integer(1))
            gain = S.factor(sum(kernel[3, j] - kernel[lower, j]
                                for j in upper_set))
            polynomial = S.cancel(e * f * gain / (kappa * prefactor))
            require(S.fraction(polynomial)[1] == 1, "polynomial denominator")
            summary = bernstein_certificate(polynomial)
            require(summary == EXPECTED[len(summaries)], "certificate counts")
            polynomials.append(polynomial)
            summaries.append(summary)
            print(f"B={BITS[lower]}, U={[BITS[j] for j in upper_set]}: "
                  f"degree={summary[0]}, positive={summary[1]}, "
                  f"zero={summary[2]}, min_positive={summary[3]}", flush=True)
    require(sum(item[1] for item in summaries) == 2034
            and sum(item[2] for item in summaries) == 198, "total counts")

    t, u = 1 - p, 1 - q
    compact_10 = ((1 - kappa) * (p + t * (1 - b) + t * (1 - b * c))
                  + kappa * p * (1 - b * t + t * (1 - u) * (1 - b * c)))
    compact_01 = ((1 - kappa) * (1 - b * c * t)
                  + kappa * p * (1 - b * t * (c + (1 - c) * u)))
    require(S.expand(polynomials[1] - compact_10) == 0, "first compact formula")
    require(S.expand(polynomials[4] - compact_01) == 0, "second compact formula")

    parameters = (F(2, 3), F(3, 5), F(4, 7), F(2, 5), F(3, 7))
    values = dict(zip(VARIABLES, map(S.Rational, parameters)))
    raw, transient_count = raw_return(parameters)
    require(transient_count == 28, "raw transient-state count")
    require(all(S.Rational(raw[i][j]) == kernel[i, j].subs(values)
                for i in range(4) for j in range(4)), "raw ARW comparison")
    print("raw ARW: 28 transient states; all 16 entries agree exactly")
    cylinder_obstruction(descent)
    print("PASS: six reconstructed symbolic identities; "
          "2,232 nonnegative Bernstein coefficients; no files written")


if __name__ == "__main__":
    main()
