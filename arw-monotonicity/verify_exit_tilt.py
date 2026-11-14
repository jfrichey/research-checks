"""Exact harmonic exit-tilt identities on one two-port ARW domain.

Run with Python's standard library:
    python3 -B research/arw/monotonicity/verify_exit_tilt.py

AI-written exploratory certificate. The general proof is in response-renewal.md.
The raw instruction-chain solver is shared with verify_boundary_response.py;
the tilted kernel and predicted likelihood ratios are constructed separately.
This fixed example checks all stable backgrounds, binary active inputs, and
one three-particle input. It is not a monotonicity search.
"""

from fractions import Fraction as F
from itertools import product
from math import prod
from collections import defaultdict, deque
import json

from verify_boundary_response import multiply, raw_laws, raw_response


def harmonic_extension(vertices, routes, boundary):
    """Solve h=K_D h+K_boundary z by rational elimination."""
    n = len(vertices)
    matrix = [
        [F(x == y) - routes[x].get(y, F(0)) for y in vertices]
        + [sum(w * boundary[y] for y, w in routes[x].items()
               if y not in vertices)]
        for x in vertices
    ]
    for col in range(n):
        pivot = next(i for i in range(col, n) if matrix[i][col])
        matrix[col], matrix[pivot] = matrix[pivot], matrix[col]
        d = matrix[col][col]
        matrix[col] = [w / d for w in matrix[col]]
        for i in range(n):
            if i != col:
                d = matrix[i][col]
                matrix[i] = [a - d * b for a, b in zip(matrix[i], matrix[col])]
    return {x: matrix[i][-1] for i, x in enumerate(vertices)}


def raw_jump_moments(vertices, routes, sleep, cost, configurations):
    """Solve a raw Feynman--Kac system without a current truncation."""
    vi = {x: i for i, x in enumerate(vertices)}
    states = set(configurations.values())
    waiting = deque(states)
    transitions = {}
    while waiting:
        state = waiting.popleft()
        if all(a <= 0 for a in state):
            continue
        i = min(i for i, count in enumerate(state) if count > 0)
        x = vertices[i]
        row = defaultdict(F)
        asleep = list(state)
        if asleep[i] == 1:
            asleep[i] = -1
        row[tuple(asleep)] += sleep[x]
        for y, weight in routes[x].items():
            end = list(state)
            end[i] -= 1
            if y in vi:
                j = vi[y]
                end[j] = 2 if end[j] == -1 else end[j] + 1
            row[tuple(end)] += (1 - sleep[x]) * weight * cost.get((x, y), F(1))
        transitions[state] = row
        for end in row:
            if end not in states:
                states.add(end)
                waiting.append(end)
    transient = sorted(transitions)
    terminal = sorted(states - set(transient))
    ti, ei = ({s: i for i, s in enumerate(v)} for v in (transient, terminal))
    n, m = len(transient), len(terminal)
    matrix = [[F(i == j) for j in range(n)] + [F(0)] * m for i in range(n)]
    for state, row in transitions.items():
        for end, weight in row.items():
            if end in ti:
                matrix[ti[state]][ti[end]] -= weight
            else:
                matrix[ti[state]][n + ei[end]] += weight
    for col in range(n):
        pivot = next(i for i in range(col, n) if matrix[i][col])
        matrix[col], matrix[pivot] = matrix[pivot], matrix[col]
        d = matrix[col][col]
        matrix[col] = [w / d for w in matrix[col]]
        for i in range(n):
            if i != col:
                d = matrix[i][col]
                if d:
                    matrix[i] = [a - d * b for a, b in zip(matrix[i], matrix[col])]
    rows = {}
    for key, start in configurations.items():
        law = {start: F(1)} if start in ei else dict(zip(terminal, matrix[ti[start]][n:]))
        rows[key] = {tuple(int(a == -1) for a in end): weight
                     for end, weight in law.items() if weight}
    return rows, n


def edge_current_audit(vertices, routes, sleep):
    cost = {(0, 1): F(3, 2)}
    weighted = {x: {y: w * cost.get((x, y), F(1)) for y, w in row.items()}
                for x, row in routes.items()}
    boundary = {-1: F(1), -2: F(1)}
    h = harmonic_extension(vertices, weighted, boundary)
    assert h == {0: F(7, 5), 1: F(6, 5)}
    tilted = {x: {y: w * (h[y] if y in h else F(1)) / h[x]
                  for y, w in row.items()} for x, row in weighted.items()}
    assert all(sum(row.values()) == 1 for row in tilted.values())
    configurations = {}
    for bits in product((0, 1), repeat=len(vertices)):
        for i, x in enumerate(vertices):
            state = [-b for b in bits]
            state[i] = 2 if bits[i] else 1
            configurations[x, bits] = tuple(state)
    rows, n = raw_jump_moments(vertices, routes, sleep, cost, configurations)
    ordinary = raw_response(vertices, tilted, sleep)
    entries = 0
    for (x, bits), row in rows.items():
        original_weight = h[x] * prod(h[y] ** b for y, b in zip(vertices, bits))
        predicted = {end: w * prod(h[y] ** b for y, b in zip(vertices, end))
                     / original_weight for end, w in row.items()}
        actual = defaultdict(F)
        for (end, _), w in ordinary.rows[x, bits].items():
            actual[end] += w
        assert predicted == actual
        assert sum(predicted.values()) == 1
        entries += len(predicted)
    eigenvalue = F(2, 3) * h[0] + F(1, 3) * h[1]
    assert eigenvalue == F(4, 3)
    states = tuple(product((0, 1), repeat=len(vertices)))
    matrix = [[sum((F(2, 3) if x == 0 else F(1, 3))
                   * rows[x, bits].get(end, F(0)) for x in vertices)
               for end in states] for bits in states]
    spectrum = [F(0), F(1, 3), F(1, 6), eigenvalue]
    check_spectrum(matrix, spectrum)
    return {"raw_transient_states": n, "terminal_entries": entries,
            "harmonic_weights": {str(x): str(w) for x, w in h.items()},
            "iid_drive_eigenvalue": str(eigenvalue),
            "full_spectrum": list(map(str, spectrum))}


def check_spectrum(matrix, roots):
    """Newton identities: traces through dimension determine the polynomial."""
    power = matrix
    for k in range(1, len(matrix) + 1):
        assert sum(power[i][i] for i in range(len(matrix))) == sum(x ** k for x in roots)
        if k < len(matrix):
            power = multiply(power, matrix)


def audit():
    vertices = (0, 1)
    # Two distinct killed exits, asymmetric feedback, and one self-loop.
    routes = {
        0: {0: F(1, 5), 1: F(2, 5), -1: F(1, 3), -2: F(1, 15)},
        1: {0: F(1, 2), -1: F(1, 6), -2: F(1, 3)},
    }
    sleep = {0: F(1, 3), 1: F(2, 5)}
    boundary = {-1: F(2), -2: F(1, 3)}
    h = harmonic_extension(vertices, routes, boundary)
    assert h == {0: F(13, 9), 1: F(7, 6)}
    tilted = {
        x: {y: w * (h[y] if y in h else boundary[y]) / h[x]
            for y, w in routes[x].items()}
        for x in vertices
    }
    assert all(sum(row.values()) == 1 for row in tilted.values())

    original = raw_response(vertices, routes, sleep)
    transformed = raw_response(vertices, tilted, sleep)
    assert original.ports == transformed.ports
    ports = original.ports

    def chi(bits):
        return prod(h[x] ** bit for x, bit in zip(vertices, bits))

    def exit_weight(counts):
        return prod(boundary[o] ** count for o, count in zip(ports, counts))

    addition_entries = 0
    for (source, bits), law in original.rows.items():
        initial_weight = h[source] * chi(bits)
        predicted = {
            (end, counts): w * exit_weight(counts) * chi(end) / initial_weight
            for (end, counts), w in law.items()
        }
        assert predicted == transformed.rows[source, bits]
        assert sum(predicted.values()) == 1
        addition_entries += len(predicted)

    inputs = {bits: bits for bits in product((0, 1), repeat=2)}
    inputs[(2, 1)] = (2, 1)
    p0, laws, n0 = raw_laws(vertices, routes, sleep, inputs)
    p1, tilted_laws, n1 = raw_laws(vertices, tilted, sleep, inputs)
    assert p0 == p1 == ports
    active_entries = 0
    for counts, law in laws.items():
        initial_weight = chi(counts)
        predicted = {
            (end, output): w * exit_weight(output) * chi(end) / initial_weight
            for (end, output), w in law.items()
        }
        assert predicted == tilted_laws[counts]
        assert sum(predicted.values()) == 1
        active_entries += len(predicted)

    # A positive eigenvector of each polynomial addition response, and of
    # their fixed probabilistic mixture, gives an exactly known eigenvalue.
    rho = {0: F(2, 3), 1: F(1, 3)}
    eigenvalue = sum(rho[x] * h[x] for x in vertices)
    assert eigenvalue == F(73, 54)
    for bits in product((0, 1), repeat=2):
        weighted = sum(
            rho[x] * w * exit_weight(output) * chi(end)
            for x in vertices
            for (end, output), w in original.rows[x, bits].items()
        )
        assert weighted == eigenvalue * chi(bits)

    states = tuple(product((0, 1), repeat=2))
    matrix = [[sum(rho[x] * w * exit_weight(output)
                   for x in vertices
                   for (end, output), w in original.rows[x, bits].items()
                   if end == target)
               for target in states] for bits in states]
    spectrum = [F(0), F(31, 54), F(4, 27), eigenvalue]
    check_spectrum(matrix, spectrum)

    # Complete the resident walks virtually: the two-port PGF, including
    # the retained particles, agrees with independent walk endpoints.
    # Evaluate at a second point to avoid testing a normalization only.
    second_boundary = {-1: F(3, 2), -2: F(4, 5)}
    second_h = harmonic_extension(vertices, routes, second_boundary)
    for counts, law in laws.items():
        completed = sum(
            w * prod(second_boundary[o] ** count
                     for o, count in zip(ports, output))
            * prod(second_h[x] ** bit for x, bit in zip(vertices, end))
            for (end, output), w in law.items()
        )
        assert completed == prod(second_h[x] ** count
                                 for x, count in zip(vertices, counts))

    return {
        "status": "PASS",
        "harmonic_weights": {str(x): str(w) for x, w in h.items()},
        "addition_joint_entries": addition_entries,
        "active_joint_entries": active_entries,
        "addition_transient_states": original.transient_states,
        "active_transient_states": n0,
        "tilted_active_transient_states": n1,
        "iid_drive_eigenvalue": str(eigenvalue),
        "full_spectrum": list(map(str, spectrum)),
        "input_cases": len(inputs),
        "edge_current": edge_current_audit(vertices, routes, sleep),
    }


if __name__ == "__main__":
    print(json.dumps(audit(), indent=2))
