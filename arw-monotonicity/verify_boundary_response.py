"""Exact finite boundary-response audits; Python standard library only.

Run: python3 -B research/arw/monotonicity/verify_boundary_response.py

AI-written exploratory mathematics from Jacob Richey's assigned question.
This checks three fixed motifs, not monotonicity or a parameter sweep.
No existing verifier, checkpoint, or certificate is imported or modified.
All probabilities are fractions. Raw ARW uses largest-site priority; the
queue composition uses smallest-site priority and complete local responses.
"""

from collections import defaultdict, deque
from dataclasses import dataclass
from fractions import Fraction as F
from itertools import product
import argparse
import json


SINK = -1


def absorption(starts, step, terminal):
    """Enumerate the reachable finite chain and solve (I-Q)H=S exactly."""
    states = list(dict.fromkeys(starts))
    seen = set(states)
    waiting = deque(states)
    rows = {}
    while waiting:
        state = waiting.popleft()
        if terminal(state):
            continue
        row = {end: w for end, w in step(state).items() if w}
        assert sum(row.values()) == 1 and min(row.values()) > 0
        rows[state] = row
        for end in row:
            if end not in seen:
                seen.add(end)
                states.append(end)
                waiting.append(end)
        assert len(states) < 2000, "The audit is deliberately bounded."

    transient = [s for s in states if not terminal(s)]
    ends = [s for s in states if terminal(s)]
    ti = {s: i for i, s in enumerate(transient)}
    ei = {s: i for i, s in enumerate(ends)}
    n, m = len(transient), len(ends)
    matrix = [[F(i == j) for j in range(n)] + [F(0)] * m
              for i in range(n)]
    for s, row in rows.items():
        for end, w in row.items():
            if end in ti:
                matrix[ti[s]][ti[end]] -= w
            else:
                matrix[ti[s]][n + ei[end]] += w

    for col in range(n):
        pivot = next(i for i in range(col, n) if matrix[i][col])
        matrix[col], matrix[pivot] = matrix[pivot], matrix[col]
        d = matrix[col][col]
        matrix[col][col:] = [x / d for x in matrix[col][col:]]
        for row in range(n):
            d = matrix[row][col]
            if row == col or not d:
                continue
            matrix[row][col:] = [a - d * b for a, b in
                                 zip(matrix[row][col:], matrix[col][col:])]

    laws = {s: {s: F(1)} for s in ends}
    for s, row in zip(transient, matrix):
        law = {end: w for end, w in zip(ends, row[n:]) if w}
        assert sum(law.values()) == 1 and min(law.values()) > 0
        laws[s] = law
    return {s: laws[s] for s in starts}, n


@dataclass
class Response:
    """One active input: (source, sleeping bits) -> (new bits, exit counts)."""

    vertices: tuple
    ports: tuple
    rows: dict
    transient_states: int


def raw_laws(vertices, routes, sleep, configurations):
    """Raw ARW for arbitrary initial activity, with labeled killed departures."""
    vertices = tuple(sorted(vertices))
    vi = {v: i for i, v in enumerate(vertices)}
    ports = tuple(sorted({y for x in vertices for y in routes[x] if y not in vi}))
    oi = {v: i for i, v in enumerate(ports)}
    zero = (0,) * len(ports)
    starts = {key: (tuple(state), zero) for key, state in configurations.items()}

    def step(state):
        active, output = state
        i = max(i for i, count in enumerate(active) if count > 0)
        x = vertices[i]
        row = defaultdict(F)
        asleep = list(active)
        if asleep[i] == 1:
            asleep[i] = -1
        row[tuple(asleep), output] += sleep[x]
        for y, weight in routes[x].items():
            end, out = list(active), list(output)
            end[i] -= 1
            if y in vi:
                j = vi[y]
                end[j] = 2 if end[j] == -1 else end[j] + 1
            else:
                out[oi[y]] += 1
            row[tuple(end), tuple(out)] += (1 - sleep[x]) * weight
        return row

    laws, n = absorption(list(starts.values()), step,
                         lambda s: all(a <= 0 for a in s[0]))
    rows = {}
    for key, start in starts.items():
        rows[key] = {(tuple(int(a == -1) for a in active), output): w
                     for (active, output), w in laws[start].items()}
        mass = sum(abs(a) for a in configurations[key])
        assert all(sum(bits) + sum(output) == mass for bits, output in rows[key])
    return ports, rows, n


def raw_response(vertices, routes, sleep):
    """Kill every departure from vertices and retain its destination label."""
    vertices = tuple(sorted(vertices))
    configurations = {}
    for bits in product((0, 1), repeat=len(vertices)):
        for i, x in enumerate(vertices):
            state = [-b for b in bits]
            state[i] = 2 if bits[i] else 1
            configurations[x, bits] = tuple(state)
    ports, rows, n = raw_laws(vertices, routes, sleep, configurations)
    return Response(vertices, ports, rows, n)


def glue(*modules):
    """Resolve all mutual ports; retain stable bits and pending arrival counts."""
    vertices = tuple(sorted(v for module in modules for v in module.vertices))
    assert len(set(vertices)) == len(vertices)
    vi = {v: i for i, v in enumerate(vertices)}
    owner = {v: module for module in modules for v in module.vertices}
    ports = tuple(sorted({v for module in modules for v in module.ports
                          if v not in vi}))
    oi = {v: i for i, v in enumerate(ports)}
    starts = {}
    for bits in product((0, 1), repeat=len(vertices)):
        for x in vertices:
            pending = tuple(int(v == x) for v in vertices)
            starts[x, bits] = bits, pending, (0,) * len(ports)

    def step(state):
        bits, pending, output = state
        i = min(i for i, count in enumerate(pending) if count)
        x = vertices[i]
        module = owner[x]
        local = tuple(bits[vi[v]] for v in module.vertices)
        row = defaultdict(F)
        for (new_local, departures), w in module.rows[x, local].items():
            new_bits, queue, out = list(bits), list(pending), list(output)
            queue[i] -= 1
            for v, b in zip(module.vertices, new_local):
                new_bits[vi[v]] = b
            for y, count in zip(module.ports, departures):
                if y in vi:
                    queue[vi[y]] += count
                else:
                    out[oi[y]] += count
            assert sum(new_bits) + sum(queue) + sum(out) == (
                sum(bits) + sum(pending) + sum(output))
            row[tuple(new_bits), tuple(queue), tuple(out)] += w
        return row

    laws, n = absorption(list(starts.values()), step, lambda s: not any(s[1]))
    rows = {}
    for key, start in starts.items():
        rows[key] = {(bits, output): w for (bits, _, output), w in laws[start].items()}
    return Response(vertices, ports, rows, n)


def same_response(a, b):
    assert a.vertices == b.vertices and a.ports == b.ports
    assert a.rows == b.rows
    return sum(len(row) for row in a.rows.values())


def input_word(module, word, bits=None):
    """Sequentially add a word of inputs, with all external ports absorbing."""
    if bits is None:
        bits = (0,) * len(module.vertices)
    law = {(bits, (0,) * len(module.ports)): F(1)}
    for x in word:
        updated = defaultdict(F)
        for (old, counts), a in law.items():
            for (new, extra), b in module.rows[x, old].items():
                output = tuple(u + v for u, v in zip(counts, extra))
                updated[new, output] += a * b
        law = dict(updated)
    assert sum(law.values()) == 1
    return law


def multiply(a, b):
    return [[sum(x * y for x, y in zip(row, col)) for col in zip(*b)] for row in a]


def buffer_audit():
    # One input, one feedback portal, and the original sink.
    routes = {0: {1: F(1, 2), SINK: F(1, 2)}}
    module = raw_response((0,), routes, {0: F(1, 2)})
    portal = module.ports.index(1)
    coeff = [[[F(0)] * 2 for _ in range(2)] for _ in range(3)]
    for bit in (0, 1):
        for (end, counts), w in module.rows[0, (bit,)].items():
            coeff[counts[portal]][bit][end[0]] += w
    expected = [
        [[F(1, 4), F(1, 2)], [F(1, 8), F(1, 4)]],
        [[F(1, 4), F(0)], [F(1, 4), F(1, 4)]],
        [[F(0), F(0)], [F(1, 8), F(0)]],
    ]
    assert coeff == expected
    xy, yx = multiply(coeff[0], coeff[1]), multiply(coeff[1], coeff[0])
    commutator = [[a - b for a, b in zip(u, v)] for u, v in zip(xy, yx)]
    assert commutator == [[F(1, 8), F(0)], [F(0), F(-1, 8)]]
    assert xy[0][1] / sum(xy[0]) == F(2, 5)
    assert yx[0][1] / sum(yx[0]) == F(2, 3)
    # Coefficients of A(z)^2 = ((1+z)/2) A(z), including zero entries.
    zero = [[F(0)] * 2 for _ in range(2)]
    for h in range(5):
        lhs = [[F(0)] * 2 for _ in range(2)]
        for j in range(3):
            if 0 <= h - j < 3:
                term = multiply(coeff[j], coeff[h - j])
                lhs = [[a + b for a, b in zip(u, v)] for u, v in zip(lhs, term)]
        right0 = coeff[h] if h < 3 else zero
        right1 = coeff[h - 1] if 0 <= h - 1 < 3 else zero
        rhs = [[(a + b) / 2 for a, b in zip(u, v)] for u, v in zip(right0, right1)]
        assert lhs == rhs
    return {
        "count_matrices": coeff,
        "A0_A1_minus_A1_A0": commutator,
        "occupation_given_counts_0_then_1": F(2, 5),
        "occupation_given_counts_1_then_0": F(2, 3),
        "polynomial_identity": "A(z)^2=((1+z)/2)A(z)",
    }


def pair_audit():
    routes = {0: {1: F(1)}, 1: {0: F(1, 2), SINK: F(1, 2)}}
    sleep = {0: F(1, 2), 1: F(1, 2)}
    a, b = [raw_response((i,), routes, sleep) for i in range(2)]
    composed = glue(a, b)
    direct = raw_response((0, 1), routes, sleep)
    entries = same_response(composed, direct)
    full = input_word(composed, (0, 1))
    assert full == input_word(direct, (1, 0))
    _, simultaneous, _ = raw_laws((0, 1), routes, sleep, {"full": (1, 1)})
    assert full == simultaneous["full"]
    law = {str(bits): w for (bits, _), w in full.items()}
    assert law == {"(0, 0)": F(4, 35), "(0, 1)": F(8, 35),
                   "(1, 0)": F(9, 35), "(1, 1)": F(14, 35)}
    return {"response_outcomes_checked": entries,
            "queue_transient_states": composed.transient_states,
            "raw_transient_states": direct.transient_states,
            "all_active_terminal_law": law}


def associativity_audit():
    # One fixed three-cycle with a single sink leak; no enumeration of kernels.
    routes = {0: {1: F(1)}, 1: {2: F(1)},
              2: {0: F(1, 2), SINK: F(1, 2)}}
    sleep = dict.fromkeys(range(3), F(1, 2))
    a, b, c = [raw_response((i,), routes, sleep) for i in range(3)]
    ab, bc = glue(a, b), glue(b, c)
    left, right = glue(ab, c), glue(a, bc)
    direct = raw_response((0, 1, 2), routes, sleep)
    entries = same_response(left, right)
    same_response(left, direct)
    # The open two-site response has two possible inputs. Check the entire
    # joint law of terminal bits and total exits for the two input orders.
    for bits in product((0, 1), repeat=2):
        assert input_word(ab, (0, 1), bits) == input_word(ab, (1, 0), bits)
    full = input_word(left, (0, 1, 2))
    assert full == input_word(direct, (2, 1, 0))
    _, simultaneous, _ = raw_laws((0, 1, 2), routes, sleep, {"full": (1, 1, 1)})
    assert full == simultaneous["full"]
    return {"response_outcomes_checked": entries,
            "left_queue_transient_states": left.transient_states,
            "right_queue_transient_states": right.transient_states,
            "raw_transient_states": direct.transient_states,
            "all_active_terminal_law": {str(bits): w for (bits, _), w in full.items()}}


def processor_dictionary_audit():
    """One fixed raw-book cylinder; no random sampling or book enumeration."""
    # Raw sleep symbols are strings; jump symbols are destination integers.
    # Every unspecified instruction is a sink jump, ensuring termination.
    words = {0: ["S", "S", 0, "S", "S", 1, SINK, "S"],
             1: [0, "S", SINK, "S", SINK]}

    def symbol(x, index):
        return words[x][index] if index < len(words[x]) else SINK

    def parse(x, initial_bit):
        bits, arrows, clocks, skipped = [initial_bit], [], [0], []
        offset = 0
        for m in range(1, 13):
            g = 0
            if bits[-1]:
                while symbol(x, offset) == "S":
                    offset += 1
                    g += 1
                arrows.append(symbol(x, offset))
                offset += 1
            trial = symbol(x, offset)
            offset += 1
            bits.append(int(trial == "S"))
            if trial != "S":
                arrows.append(trial)
            skipped.append(g)
            clocks.append(offset)
            assert len(arrows) == m + initial_bit - bits[m]
        return bits, arrows, clocks, skipped

    initial_bits = (0, 0)
    books = [parse(x, initial_bits[x]) for x in range(2)]

    def processor(initial, largest):
        pending, used = list(initial), [0, 0]
        while any(pending):
            choices = [x for x in range(2) if pending[x]]
            x = max(choices) if largest else min(choices)
            bits, arrows, _, _ = books[x]
            old = used[x]
            pending[x] -= 1
            used[x] += 1
            before = old + initial_bits[x] - bits[old]
            after = used[x] + initial_bits[x] - bits[used[x]]
            for y in arrows[before:after]:
                if y != SINK:
                    pending[y] += 1
            assert sum(used) < 20
        end = tuple(books[x][0][used[x]] for x in range(2))
        clocks = tuple(books[x][2][used[x]] for x in range(2))
        return end, tuple(used), clocks

    def raw(initial, largest):
        active, used = list(initial), [0, 0]
        while any(a > 0 for a in active):
            choices = [x for x in range(2) if active[x] > 0]
            x = max(choices) if largest else min(choices)
            instruction = symbol(x, used[x])
            used[x] += 1
            if instruction == "S":
                if active[x] == 1:
                    active[x] = -1
            else:
                active[x] -= 1
                if instruction != SINK:
                    y = instruction
                    active[y] = 2 if active[y] == -1 else active[y] + 1
            assert sum(used) < 40
        return tuple(int(a == -1) for a in active), tuple(used)

    witnesses = {}
    for initial in ((1, 0), (2, 0), (2, 1)):
        low = processor(initial, False)
        assert low == processor(initial, True)
        assert (low[0], low[2]) == raw(initial, False) == raw(initial, True)
        witnesses[str(initial)] = dict(terminal_bits=low[0],
                                      processed_inputs=low[1], raw_topplings=low[2])
    assert books[0][0][:5] == [0, 1, 1, 0, 1]
    assert books[0][3][:4] == [0, 1, 1, 0]

    # The last reset is biased even in one site with a self-loop.
    routes = {0: {0: F(1, 2), SINK: F(1, 2)}}
    _, laws, _ = raw_laws((0,), routes, {0: F(1, 2)}, {"one": (1,)})
    occupied = sum(w for (bits, _), w in laws["one"].items() if bits[0])
    assert occupied == F(2, 3)
    return {"raw_book_witnesses": witnesses,
            "site_0_first_four_reset_bits": books[0][0][1:5],
            "site_0_first_four_skipped_sleep_counts": books[0][3][:4],
            "self_loop_last_reset_occupation": occupied,
            "scope": "fixed-background pathwise dictionary; one cylinder and one exact law"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reset-only", action="store_true",
                        help="only the bounded raw-stack/reset-processor audit")
    args = parser.parse_args()
    if args.reset_only:
        print(json.dumps(processor_dictionary_audit(), indent=2, default=str))
        return
    results = {"buffer": buffer_audit(), "pair_feedback": pair_audit(),
               "three_module_associativity": associativity_audit()}
    print(json.dumps(results, indent=2, default=str))


if __name__ == "__main__":
    main()
