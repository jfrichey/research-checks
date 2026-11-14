"""Exact finite-book audits for SX's common path/stack probability space.

The general law and inverse arguments are in stack-exploration.md. These
checks retain every bit and arrow, including unused tails; no files are written.
Research implementation written by Codex, not independently user-validated.
"""

from collections import Counter, deque
from fractions import Fraction as F
from functools import lru_cache
from random import Random
import json


def execute(bits, roots, arrows=None, prescribed_paths=None):
    """Ordinary labelled buffer stabilization, followed by a complete drain.

    The chosen contact is the first co-location of label zero in each copy.
    Their roots are distinct in the audit. A label's position changes when its
    arrow is emitted, before its queued arrival is processed.
    """
    queues = [deque(range(len(row))) for row in roots]
    residents = [{}, {}]
    walks = [[[z] for z in row] for row in roots]
    offsets, prefix = Counter(), {}
    contact = None

    def emit(copy, label, destination):
        walks[copy][label].append(destination)
        if destination != -1:
            queues[copy].append(label)

    def deliver(copy, sleeping):
        label = queues[copy].popleft()
        site = walks[copy][label][-1]
        offsets[copy, site] += 1
        key = copy, site, offsets[copy, site]
        if prescribed_paths is None:
            arrow = arrows[key]
        else:
            path = prescribed_paths[copy][label]
            assert path[len(walks[copy][label])-1] == site
            arrow = path[len(walks[copy][label])]
        prefix[key] = arrow
        if site in residents[copy]:
            old_label, old_key = residents[copy].pop(site)
            emit(copy, old_label, prefix[old_key])
        if sleeping and bits[key]:
            residents[copy][site] = label, key
        else:
            emit(copy, label, arrow)

    while any(queues):
        for copy in (0, 1):
            if queues[copy]:
                deliver(copy, True)
                left, right = walks[0][0][-1], walks[1][0][-1]
                if contact is None and left == right and left != -1:
                    contact = (left, len(walks[0][0])-1, len(walks[1][0])-1)
        assert sum(offsets.values()) < 1000
    supports = [set(row) for row in residents]
    ordinary_offsets = offsets.copy()

    for copy in (0, 1):
        for site in sorted(residents[copy]):
            label, key = residents[copy][site]
            emit(copy, label, prefix[key])
        residents[copy].clear()
        while queues[copy]:
            deliver(copy, False)
            assert sum(offsets.values()) < 1000
    assert all(path[-1] == -1 for row in walks for path in row)
    if prescribed_paths is not None:
        assert walks == prescribed_paths
    return dict(paths=walks, prefix=prefix, offsets=offsets,
                ordinary_offsets=ordinary_offsets,
                contact=contact, supports=supports)


def decompose(bits, arrows, roots):
    data = execute(bits, roots, arrows=arrows)
    data['tails'] = {(copy, site, index-data['offsets'][copy, site]): arrow
                     for (copy, site, index), arrow in arrows.items()
                     if index > data['offsets'][copy, site]}
    return data


def reconstruct(bits, roots, paths, tails):
    data = execute(bits, roots, prescribed_paths=paths)
    arrows = dict(data['prefix'])
    for (copy, site, index), arrow in tails.items():
        arrows[copy, site, index+data['offsets'][copy, site]] = arrow
    return arrows, data


def exchange_at_contact(bits, arrows, roots):
    data = decompose(bits, arrows, roots)
    if data['contact'] is None:
        return dict(arrows), None
    site, left, right = data['contact']
    paths = [[list(path) for path in row] for row in data['paths']]
    a, b = paths[0][0], paths[1][0]
    assert a[left] == b[right] == site
    paths[0][0] = a[:left]+b[right:]
    paths[1][0] = b[:right]+a[left:]
    changed, decoded = reconstruct(bits, roots, paths, data['tails'])
    assert decoded['contact'] == data['contact']
    return changed, data['contact']


def forced_drain_counts(arrows, roots):
    offsets = Counter()
    for copy, row in enumerate(roots):
        queue = deque(row)
        while queue:
            site = queue.popleft()
            offsets[copy, site] += 1
            destination = arrows[copy, site, offsets[copy, site]]
            if destination != -1:
                queue.append(destination)
            assert sum(offsets.values()) < 1000
    return offsets


def terminal_coin_exchange(bits, arrows, roots, target):
    """Partial involution inside each DI pool; reject changed terminal indices."""
    data = decompose(bits, arrows, roots)
    cut = tuple(data['ordinary_offsets'][copy, target] for copy in (0, 1))
    occupied = tuple(target in row for row in data['supports'])
    if not all(cut) or occupied[0] == occupied[1]:
        return dict(bits), None
    changed = dict(bits)
    for copy, index in enumerate(cut):
        changed[copy, target, index] = 1-changed[copy, target, index]
    image = decompose(changed, arrows, roots)
    recovered = tuple(image['ordinary_offsets'][copy, target] for copy in (0, 1))
    if recovered != cut:
        return dict(bits), None
    assert tuple(target in row for row in image['supports']) == occupied[::-1]
    assert image['offsets'] == data['offsets']
    inventory = lambda field: Counter(site for (copy, site, index), bit in field.items()
                                     if bit and index <= data['offsets'][copy, site])
    assert inventory(changed) == inventory(bits)
    return changed, cut


def audit():
    rng = Random(20261011)
    roots = [[0, 2], [1, 0, 2]]
    sleeps = [(2, 5), (3, 7), (5, 11), (7, 13)]
    counts = Counter()
    for sample in range(256):
        bits, arrows = {}, {}
        for copy in (0, 1):
            for site, (num, den) in enumerate(sleeps):
                for index in range(1, 129):
                    bits[copy, site, index] = int(rng.randrange(den) < num)
                    arrows[copy, site, index] = rng.choice([0, 1, 2, 3, -1]) if index <= 64 else -1
        data = decompose(bits, arrows, roots)
        rebuilt, _ = reconstruct(bits, roots, data['paths'], data['tails'])
        assert rebuilt == arrows
        assert data['offsets'] == forced_drain_counts(arrows, roots)
        changed_bits = {key: 1-value for key, value in bits.items()}
        assert decompose(changed_bits, arrows, roots)['offsets'] == data['offsets']
        changed, cut = exchange_at_contact(bits, arrows, roots)
        restored, decoded = exchange_at_contact(bits, changed, roots)
        assert restored == arrows and cut == decoded
        counts['samples_with_contact'] += cut is not None
        counts['changed_arrow_fields'] += changed != arrows
        if cut is not None:
            counts['contacts_after_both_have_jumped'] += cut[1] > 0 and cut[2] > 0
        inventory = lambda field: Counter((site, arrow) for (_, site, _), arrow in field.items())
        assert inventory(changed) == inventory(arrows)
        new_data = decompose(bits, changed, roots)
        assert new_data['tails'] == data['tails']
        hit = lambda paths: any(3 in path for row in paths for path in row)
        assert hit(data['paths']) == hit(new_data['paths'])
        counts['joint_target_3_successes'] += hit(data['paths'])
        counts['joint_target_3_failures'] += not hit(data['paths'])
        for target in (0, 3):
            coin_image, coin_cut = terminal_coin_exchange(bits, arrows, roots, target)
            coin_restored, recovered_cut = terminal_coin_exchange(coin_image, arrows, roots, target)
            assert coin_restored == bits and recovered_cut == coin_cut
            if coin_cut is not None:
                counts['terminal_coin_exchanges'] += 1
                counts['terminal_coin_bad_pairs_repaired'] += target in data['supports'][0]
    assert counts['contacts_after_both_have_jumped']
    assert counts['joint_target_3_successes'] and counts['joint_target_3_failures']
    return dict(seed=20261011, sampled_cylinders=256, **counts,
                all_full_inverses=True, all_arrow_inventories=True,
                all_residual_tails_recovered=True, drain_independent_of_bits=True,
                contact_bits_unchanged=True, joint_success_invariant=True,
                terminal_coin_inverses_and_drain_inventories=True,
                scope='Representation and partial exchanges only; occupation sign open')


def fixed_path_occupation(paths, sleeps, target):
    """Exact coin average with prescribed labelled paths and FIFO arrivals."""
    @lru_cache(None)
    def solve(positions, queue, residents):
        if not queue:
            return F(residents[target] >= 0)
        label = queue[0]
        site = paths[label][positions[label]]
        pending, steps, resting = list(queue[1:]), list(positions), list(residents)

        def emit(particle, indices, packets):
            indices[particle] += 1
            if paths[particle][indices[particle]] != -1:
                packets.append(particle)

        if resting[site] >= 0:
            emit(resting[site], steps, pending)
            resting[site] = -1
        sleeping = list(resting)
        sleeping[site] = label
        stay = solve(tuple(steps), tuple(pending), tuple(sleeping))
        emit(label, steps, pending)
        jump = solve(tuple(steps), tuple(pending), tuple(resting))
        return sleeps[site]*stay+(1-sleeps[site])*jump

    return solve((0,)*len(paths), tuple(range(len(paths))), (-1,)*len(sleeps))


def frozen_path_obstruction():
    paths = ((0, 2, 0, -1), (1, 0, -1))
    sleeps = (F(1, 2),)*3
    lower = fixed_path_occupation(paths[:1], sleeps, 0)
    upper = fixed_path_occupation(paths, sleeps, 0)
    assert lower == F(5, 8) and upper == F(9, 16)
    return dict(common_path=paths[0], extra_path=paths[1],
                lower=str(lower), upper=str(upper), difference=str(upper-lower),
                verified=True,
                scope='Fixed shared labelled paths; not an unconditional ARW counterexample')


def shared_common_path_obstruction():
    """Average the extra's entire walk; fixing the common FIFO path still fails."""
    from collections import defaultdict
    from verify_response_renewal import stopped, arw_step

    path, sleeps = (0, 2, 0, -1), (F(1, 2),)*3
    kernel = {0: {2: F(1, 2), -1: F(1, 2)}, 1: {0: F(1)}, 2: {0: F(1)}}
    start = (0, 1, (0, 1), (-1,)*3)

    def emit(state, label):
        index, extra, queue, residents = state
        choices = {path[index+1]: F(1)} if label == 0 else kernel[extra]
        out = {}
        for destination, weight in choices.items():
            i, e = (index+1, extra) if label == 0 else (index, destination)
            pending = queue if destination == -1 else queue+(label,)
            out[i, e, pending, residents] = weight
        return out

    def transition(state):
        index, extra, queue, residents = state
        label = queue[0]
        site = path[index] if label == 0 else extra
        resting = list(residents)
        old = resting[site]
        resting[site] = -1
        base = (index, extra, queue[1:], tuple(resting))
        states = emit(base, old) if old >= 0 else {base: F(1)}
        out = defaultdict(F)
        for current, weight in states.items():
            row = list(current[3])
            row[site] = label
            out[current[:3]+(tuple(row),)] += weight*sleeps[site]
            for final, probability in emit(current, label).items():
                out[final] += weight*(1-sleeps[site])*probability
        return out

    laws, states = stopped([start], transition, lambda state: not state[2])
    upper = sum(weight for final, weight in laws[start].items() if final[3][0] >= 0)
    lower = fixed_path_occupation((path,), sleeps, 0)
    assert lower == F(5, 8) and upper == F(269, 448)
    K = [[F(0), F(0), F(1, 2), F(1, 2)],
         [F(1), F(0), F(0), F(0)], [F(1), F(0), F(0), F(0)]]
    starts = [(1, 0, 0), (1, 1, 0)]
    raw, _ = stopped(starts, lambda state: arw_step(K, sleeps, state),
                     lambda state: all(a <= 0 for a in state))
    ordinary = [sum(weight for final, weight in raw[s].items() if final[0] == -1)
                for s in starts]
    assert ordinary == [F(4, 7), F(41, 70)]
    return dict(common_path=path, conditional_path_probability='1/4',
                lower=str(lower), upper=str(upper), difference=str(upper-lower),
                exact_conditional_states=states, unconditional=list(map(str, ordinary)),
                verified=True, scope='Shared common full path in the specified FIFO representation')


if __name__ == '__main__':
    print(json.dumps(dict(representation=audit(), frozen_paths=frozen_path_obstruction(),
                          shared_common_path=shared_common_path_obstruction()), indent=2))
