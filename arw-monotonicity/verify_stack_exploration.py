"""Exact checks for the local constructions in stack-exploration.md.

Research code written by Codex from the project's saved record conventions.
These finite checks do not establish a global add-one coupling. General local
proofs and the remaining domain restrictions are in the maintained reports.
Run with Python 3; standard library only. No files are written.
"""

from collections import Counter, deque
from fractions import Fraction as F
from itertools import product
from random import Random
import json


def stabilize(records, pending, initial_residents=None):
    """Potential-arrow buffer model; sink is -1, keys are (site, arrival)."""
    residents = dict(initial_residents or {})
    arrivals = {}
    queue = deque(pending)
    trace, used_bits, used_arrows = [], set(), set()
    while queue:
        site = queue.popleft()
        index = arrivals.get(site, 0) + 1
        arrivals[site] = index
        key = site, index
        bit, arrow = records[key]  # Missing records are an error, not truncation.
        used_bits.add(key)
        outgoing = []
        if site in residents:
            old_key, old_arrow = residents.pop(site)
            used_arrows.add(old_key)
            outgoing.append(old_arrow)
        if bit:
            residents[site] = key, arrow
        else:
            used_arrows.add(key)
            outgoing.append(arrow)
        trace.append((site, index, bit, outgoing))
        queue.extend(destination for destination in outgoing if destination != -1)
        assert len(trace) < 1000, "Unexpected loop in a specified finite witness"
    return set(residents), used_bits, used_arrows, trace


def switch(records, active_extra):
    """A FIXED permutation, with contact index known in advance."""
    result = dict(records)
    contact_index = 2 if active_extra else 1
    pairs = [((0, 0, 1), (1, 0, 1)),
             ((0, 1, 1), (1, 1, contact_index))]
    for left, right in pairs:
        result[left], result[right] = result[right], result[left]
    return result


def record_weight(records, kernel, sleeps):
    weight = F(1)
    for (_, site, _), (bit, arrow) in records.items():
        weight *= sleeps[site] if bit else 1 - sleeps[site]
        weight *= kernel[site][arrow]
    return weight


def run_pair(records, active_extra, resident_arrow):
    outcomes = []
    for copy in (0, 1):
        book = {(site, index): record for (side, site, index), record
                in records.items() if side == copy}
        pending = [1, 0] if copy and active_extra else [0]
        residents = {1: ((1, 0), resident_arrow)} if copy and not active_extra else {}
        outcomes.append(stabilize(book, pending, residents))
    return outcomes


def check_record_witness(kind, active_extra):
    # The two original examples use symmetric jump rows. Sleep rates are
    # deliberately heterogeneous; none of the permutations exchanges sites.
    if kind == "BS":
        kernel = {0: {1: F(1, 2), -1: F(1, 2)},
                  1: {0: F(1, 2), -1: F(1, 2)}}
        contact_record, resident_arrow = (0, -1), -1
    else:
        kernel = {0: {1: F(1, 2), -1: F(1, 2)},
                  1: {0: F(1, 2), 2: F(1, 2)},
                  2: {1: F(1, 2), -1: F(1, 2)}}
        contact_record, resident_arrow = (1, 2), 2
    sleeps = {0: F(2, 5), 1: F(3, 7), 2: F(5, 11)}
    contact_index = 2 if active_extra else 1
    count, mass = 0, F(0)
    for root_padding, spare_bit, spare_arrow in product(
            kernel[0], (0, 1), kernel[1]):
        records = {
            (0, 0, 1): (1, root_padding),
            (1, 0, 1): (0, 1),
            (0, 1, 1): (spare_bit, spare_arrow),
            (1, 1, contact_index): contact_record,
        }
        if active_extra:
            records[1, 1, 1] = 1, resident_arrow
        if kind == "R":
            records[1, 2, 1] = 0, -1
        changed = switch(records, active_extra)
        assert switch(changed, active_extra) == records
        weight = record_weight(records, kernel, sleeps)
        assert weight == record_weight(changed, kernel, sleeps)
        if not active_extra:
            weight *= kernel[1][resident_arrow]  # Initial sleeper's padding.
        before = run_pair(records, active_extra, resident_arrow)
        after = run_pair(changed, active_extra, resident_arrow)
        assert (0 in before[0][0], 0 in before[1][0]) == (True, False)
        assert (0 in after[0][0], 0 in after[1][0]) == (False, True)
        if kind == "R":
            assert (2, 1) in before[1][1] and (2, 1) in before[1][2]
            assert (2, 1) not in after[1][1] and (2, 1) not in after[1][2]
        resident_key = (1, 1) if active_extra else (1, 0)
        assert resident_key in before[1][2]
        assert resident_key not in after[1][2]
        mass += weight
        count += 1
    p, q = sleeps[0], sleeps[1]
    expected = p * (1-p) * kernel[0][1] * kernel[1][resident_arrow]
    expected *= (q if contact_record[0] else 1-q) * kernel[1][contact_record[1]]
    if active_extra:
        expected *= q
    if kind == "R":
        expected *= (1-sleeps[2]) * kernel[2][-1]
    assert mass == expected  # The spare-record probabilities summed to one.
    return dict(example=kind, active_extra=active_extra,
                padded_cylinders=count, bad_mass=str(mass), image_mass=str(mass))


def colour_law(hits, split=None):
    raw = {}
    for bits in product((0, 1), repeat=len(hits)):
        allowed = any(bits) if split is None else any(bits[:split]) and any(bits[split:])
        if not allowed:
            continue
        weight = F(1)
        for bit, probability in zip(bits, hits):
            weight *= probability if bit else 1-probability
        if weight:
            raw[bits] = weight
    total = sum(raw.values(), F(0))
    return {bits: weight / total for bits, weight in raw.items()}


def check_long_excursion(length):
    """One fixed switch, with an arbitrarily long untouched creator block."""
    kernel = {0: {1: F(1, 2), -1: F(1, 2)},
              1: {0: F(1, 2), 2: F(1, 2)},
              2: {2: F(1, 2), -1: F(1, 2)}}
    sleeps = {0: F(2, 5), 1: F(3, 7), 2: F(5, 11)}
    records = {(0, 0, 1): (1, -1), (1, 0, 1): (0, 1),
               (0, 1, 1): (0, 2), (1, 1, 1): (1, 2), (1, 1, 2): (1, 2)}
    for index in range(1, length + 1):
        records[1, 2, index] = 0, (2 if index < length else -1)
    changed = switch(records, True)
    assert switch(changed, True) == records
    assert record_weight(records, kernel, sleeps) == record_weight(changed, kernel, sleeps)
    before = run_pair(records, True, 2)
    after = run_pair(changed, True, 2)
    assert (0 in before[0][0], 0 in before[1][0]) == (True, False)
    assert (0 in after[0][0], 0 in after[1][0]) == (False, True)
    for index in range(1, length + 1):
        assert (2, index) in before[1][1]
        assert (2, index) not in after[1][1]
    return dict(excursion_length=length, all_excursion_records_become_unused=True)


def stopped_walk(records, copy, source, offsets):
    """Read one explorer to its first sleep/sink, retaining its full records.

    Existing resident arrows would be held as pending packets. Thus the
    explorer's word is an ordinary one-particle word even with residents.
    """
    words, path = {}, []
    site = source
    while site != -1:
        word = words.setdefault(site, [])
        index = offsets.get((copy, site), 0) + len(word) + 1
        record = records[copy, site, index]
        word.append(record)
        bit, arrow = record
        path.append((site, index, bit, arrow))
        if bit:
            return words, path, site
        site = arrow
        assert len(path) < 100, "Unexpected loop in a finite witness"
    return words, path, -1


def exchange_record_blocks(records, words, offsets):
    """Exchange specified local prefixes below fixed offsets, retaining tails.

    For each site the two suffixes W_- T_-, W_+ T_+ become
    W_+ T_-, W_- T_+. Per-copy finite padding lengths may change.
    """
    changed = {}
    sites = {site for _, site, _ in records}
    for copy in (0, 1):
        for site in sites:
            book = [records[copy, site, index]
                    for index in range(1, 1 + sum(side == copy and z == site
                                                 for side, z, _ in records))]
            cut = offsets.get((copy, site), 0)
            old_word = words[copy].get(site, [])
            new_word = words[1-copy].get(site, [])
            replacement = book[:cut] + new_word + book[cut + len(old_word):]
            for index, record in enumerate(replacement, 1):
                changed[copy, site, index] = record
    return changed


def switch_walk_prefixes(records, source, offsets):
    """Exchange variable-length stopped words, preserving each copy's tail."""
    walks = [stopped_walk(records, copy, source, offsets) for copy in (0, 1)]
    return exchange_record_blocks(records, [walk[0] for walk in walks], offsets), walks


def check_common_wakeup():
    """Two common particles; the suppressed avalanche returns to both sites."""
    kernel = {0: {1: F(1, 2), -1: F(1, 2)},
              1: {2: F(1, 2), -1: F(1, 2)},
              2: {0: F(1, 3), 1: F(1, 3), -1: F(1, 3)}}
    sleeps = {0: F(2, 5), 1: F(3, 7), 2: F(5, 11)}
    records = {
        (0, 2, 1): (1, -1), (1, 2, 1): (1, 0),
        (0, 0, 1): (1, -1), (1, 0, 1): (0, 1),
        (0, 1, 1): (0, -1), (1, 1, 1): (1, 2),
        (1, 1, 2): (1, 2), (1, 2, 2): (0, 1),
        (1, 0, 2): (0, -1), (1, 1, 3): (0, -1),
        (1, 2, 3): (0, -1),
    }
    changed = switch(records, True)
    assert switch(changed, True) == records
    assert record_weight(records, kernel, sleeps) == record_weight(changed, kernel, sleeps)
    endpoints = []
    for field in (records, changed):
        pair = []
        for copy in (0, 1):
            book = {(z, i): r for (side, z, i), r in field.items() if side == copy}
            pair.append(stabilize(book, [2, 1, 0] if copy else [2, 0]))
        endpoints.append([sorted(result[0]) for result in pair])
    assert endpoints == [[[0, 2], []], [[1, 2], [0, 1, 2]]]
    return dict(before=endpoints[0], after=endpoints[1],
                common_sleeper_woken_before_switch=True,
                original_avalanche_returns_to_both_switch_sites=True)


def check_variable_walks():
    """Recover selected addresses despite returns and unequal local lengths."""
    kernel = {0: {0: F(1, 4), 1: F(1, 4), 3: F(1, 4), -1: F(1, 4)},
              1: {2: F(1, 2), -1: F(1, 2)},
              2: {0: F(1, 3), 1: F(1, 3), -1: F(1, 3)},
              3: {0: F(1, 2), -1: F(1, 2)}}
    sleeps = {0: F(2, 5), 1: F(3, 7), 2: F(5, 11), 3: F(7, 13)}
    checks = []
    # Common creator at 2; extra creator at 1 in the upper copy.
    offsets = {(0, 2): 1, (1, 2): 1, (1, 1): 1}
    for lower_loops, upper_loops in ((0, 0), (1, 3), (4, 1)):
        for padding in (0, 1):
            records = {(copy, z, i): (padding, -1)
                       for copy in (0, 1) for z in range(4) for i in range(1, 18)}
            records[0, 2, 1] = 1, -1
            records[1, 2, 1] = 1, 0
            records[1, 1, 1] = 1, 2
            for i in range(1, lower_loops + 1):
                records[0, 0, i] = 0, 0
            records[0, 0, lower_loops + 1] = 1, -1
            for i in range(1, upper_loops + 1):
                records[1, 0, i] = 0, 0
            records[1, 0, upper_loops + 1] = 0, 3
            records[1, 3, 1] = 0, 0
            records[1, 0, upper_loops + 2] = 0, 1
            records[1, 1, 2] = 1, 2
            # The extra wakes the common sleeper; both return through 0/1.
            records[1, 2, 2] = 0, 1
            records[1, 0, upper_loops + 3] = 0, -1
            records[1, 1, 3] = 0, -1
            records[1, 2, 3] = 0, -1
            changed, walks = switch_walk_prefixes(records, 0, offsets)
            restored, reverse_walks = switch_walk_prefixes(changed, 0, offsets)
            assert restored == records
            assert [w[2] for w in walks] == [0, 1]
            assert [w[2] for w in reverse_walks] == [1, 0]
            assert [w[0] for w in reverse_walks] == [w[0] for w in reversed(walks)]
            assert record_weight(records, kernel, sleeps) == record_weight(changed, kernel, sleeps)
            endpoints = []
            for field in (records, changed):
                pair = []
                for copy in (0, 1):
                    book = {(z, i): r for (side, z, i), r in field.items() if side == copy}
                    pair.append(stabilize(book, [2, 1, 0] if copy else [2, 0]))
                endpoints.append([sorted(result[0]) for result in pair])
            assert endpoints == [[[0, 2], []], [[1, 2], [0, 1, 2]]]
            checks.append(dict(lower_loops=lower_loops, upper_loops=upper_loops,
                               padding=padding, inverse=True, equal_weight=True))
    return checks


def check_terminal_overlap(return_fields=False):
    """Two distinct bad cylinders have the same terminal-transposition image."""
    kernel = {0: {1: F(1, 2), -1: F(1, 2)},
              1: {0: F(1, 2), 2: F(1, 2)},
              2: {1: F(1, 2), -1: F(1, 2)}}
    sleeps = {0: F(2, 5), 1: F(3, 7), 2: F(5, 11)}
    image = {
        (0, 0, 1): (0, 1), (0, 0, 2): (0, -1),
        (0, 1, 1): (0, 2), (0, 1, 2): (0, 0),
        (0, 2, 1): (0, 1), (0, 2, 2): (1, -1),
        (1, 0, 1): (1, 1),
        (1, 1, 1): (0, 2), (1, 1, 2): (1, 2),
        (1, 2, 1): (0, -1), (1, 2, 2): (1, 1),
    }

    def run(field):
        outcomes = []
        for copy in (0, 1):
            book = {(z, i): r for (side, z, i), r in field.items() if side == copy}
            outcomes.append(stabilize(book, [0, 1, 2] if copy else [0, 2]))
        return outcomes

    def exchange(field, n, m):
        result = dict(field)
        a, b = (0, 0, n), (1, 0, m)
        result[a], result[b] = result[b], result[a]
        return result

    image_out = run(image)
    assert tuple(0 in outcome[0] for outcome in image_out) == (False, True)
    sources = []
    for n in (1, 2):
        source = exchange(image, n, 1)
        out = run(source)
        assert tuple(0 in outcome[0] for outcome in out) == (True, False)
        cuts = [max(i for z, i in outcome[1] if z == 0) for outcome in out]
        assert cuts == [n, 1]
        assert exchange(source, *cuts) == image
        assert record_weight(source, kernel, sleeps) == record_weight(image, kernel, sleeps)
        sources.append(source)
    assert sources[0] != sources[1]
    if return_fields:
        return kernel, sleeps, image, sources
    return dict(lower_input=[0, 2], upper_input=[0, 1, 2],
                distinct_equal_weight_bad_sources=2, common_good_image=True,
                source_terminal_indices=[[1, 1], [2, 1]],
                cylinder_weight=str(record_weight(image, kernel, sleeps)))


def creator_contact_offsets(records):
    """Upper creator at 0, immediate extra sleep at 1, root 2 first hits 0."""
    words, _, endpoint = stopped_walk(records, 1, 0, {})
    assert endpoint == 0
    offsets = {(1, z): len(word) for z, word in words.items()}
    extra_index = offsets.get((1, 1), 0) + 1
    assert records[1, 1, extra_index][0] == 1
    offsets[1, 1] = extra_index
    site, prefix = 2, []
    while site != 0:
        index = offsets.get((1, site), 0) + 1
        bit, arrow = records[1, site, index]
        assert bit == 0 and arrow != -1  # This witness reaches the contact.
        offsets[1, site] = index
        prefix.append((site, index, bit, arrow))
        site = arrow
        assert len(prefix) < 100
    return offsets, prefix


def staged_creator_pair(records, copy):
    """Run root 0, upper extra 1, then root 2; hold woken residents' arrows."""
    counts, residents, pending, trace = {}, {}, deque(), []

    def visit(site):
        index = counts.get(site, 0) + 1
        counts[site] = index
        bit, arrow = records[copy, site, index]
        if site in residents:
            old_arrow = residents.pop(site)
            if old_arrow != -1:
                pending.append(old_arrow)
        if bit:
            residents[site] = arrow
        trace.append((site, index, bit, arrow))
        return -1 if bit else arrow

    def walk(site):
        while site != -1:
            site = visit(site)

    walk(0)
    if copy:
        assert 1 not in residents
        assert visit(1) == -1 and 1 in residents
    walk(2)
    residual_start = len(trace)
    while pending:
        destination = visit(pending.popleft())
        if destination != -1:
            pending.append(destination)
    return set(residents), trace[residual_start:]


def check_creator_contact_switch():
    """Exchange a lower creator word with an upper word after first contact.

    The upper explorer wakes BOTH the extra at 1 and the common creator at 0
    before it first sleeps. Thus this lies outside the earlier no-contact lemma.
    """
    kernel = {0: {z: F(1, 5) for z in (0, 1, 2, 3, -1)},
              1: {0: F(1, 2), 2: F(1, 2)},
              2: {1: F(1, 2), -1: F(1, 2)},
              3: {0: F(1, 2), -1: F(1, 2)}}
    sleeps = {0: F(2, 5), 1: F(3, 7), 2: F(5, 11), 3: F(7, 13)}
    checks = []
    for lower_loops, creator_loops, explorer_loops in ((0, 0, 0), (2, 1, 3), (5, 3, 1)):
        records = {(c, z, i): (0, 2 if z == 1 else -1)
                   for c in (0, 1) for z in range(4) for i in range(1, 31)}
        for i in range(1, lower_loops + 1):
            records[0, 0, i] = 0, 0
        records[0, 0, lower_loops + 1] = 1, 3
        for i in range(1, creator_loops + 1):
            records[1, 0, i] = 0, 0
        records[1, 0, creator_loops + 1] = 1, 2
        records[1, 1, 1] = 1, 2
        records[1, 2, 1] = 0, 1
        records[1, 1, 2] = 0, 0
        for i in range(creator_loops + 2, creator_loops + 2 + explorer_loops):
            records[1, 0, i] = 0, 0
        records[1, 0, creator_loops + 2 + explorer_loops] = 0, 3
        records[1, 3, 1] = 1, -1
        records[1, 2, 2] = 1, -1
        records[1, 2, 3] = 0, -1
        offsets, prefix = creator_contact_offsets(records)
        changed, _ = switch_walk_prefixes(records, 0, offsets)
        reverse_offsets, reverse_prefix = creator_contact_offsets(changed)
        assert (reverse_offsets, reverse_prefix) == (offsets, prefix)
        restored, _ = switch_walk_prefixes(changed, 0, reverse_offsets)
        assert restored == records
        assert record_weight(records, kernel, sleeps) == record_weight(changed, kernel, sleeps) > 0
        before = [staged_creator_pair(records, copy) for copy in (0, 1)]
        after = [staged_creator_pair(changed, copy) for copy in (0, 1)]
        assert [result[0] for result in before] == [{0}, {3}]
        assert [result[0] for result in after] == [{3}, {0}]
        assert all(z != 0 for z, _, _, _ in before[1][1])
        for field, expected in ((records, before), (changed, after)):
            for copy in (0, 1):
                book = {(z, i): r for (side, z, i), r in field.items() if side == copy}
                fifo = stabilize(book, [0, 1, 2] if copy else [0, 2])
                assert fifo[0] == expected[copy][0]
        checks.append(dict(lower_loops=lower_loops, creator_loops=creator_loops,
                           explorer_loops=explorer_loops, recovered_first_contact=True,
                           before=[[0], [3]], after=[[3], [0]], inverse=True))
    return checks


def word_vertices(word):
    """Sleep ends at a decorated absorbing port; its potential arrow is kept."""
    vertices = [word[0][0]]
    for site, bit, arrow in word:
        assert vertices[-1] == site
        vertices.append(("sleep", site, arrow) if bit else arrow)
    return vertices


def fomin_switch(first, second, include_two_sleeps=False):
    """Fomin's switch with an invariant restriction on the endpoint pair.

    The retained LE indices are incoming-edge indices, NOT last visit indices.
    Source: https://arxiv.org/pdf/math/0004083, proof of Theorem 6.1.
    """
    first_vertices, second_vertices = word_vertices(first), word_vertices(second)
    endpoints = [vertices[-1] for vertices in (first_vertices, second_vertices)]
    one_exit = sum(endpoint == -1 for endpoint in endpoints) == 1
    distinct_sleeps = (include_two_sleeps and all(isinstance(e, tuple) for e in endpoints)
                       and endpoints[0][1] != endpoints[1][1])
    if not (one_exit or distinct_sleeps):
        return first, second, None
    loop_erased, positions = [], {}
    for index, site in enumerate(first_vertices):
        if site in positions:
            loop_erased = loop_erased[:positions[site] + 1]
            positions = {z: i for i, (z, _) in enumerate(loop_erased)}
        else:
            positions[site] = len(loop_erased)
            loop_erased.append((site, index))
    second_sites = {z for z in second_vertices if isinstance(z, int) and z >= 0}
    for site, first_cut in loop_erased:
        if isinstance(site, int) and site >= 0 and site in second_sites:
            second_cut = second_vertices.index(site)
            return (first[:first_cut] + second[second_cut:],
                    second[:second_cut] + first[first_cut:],
                    (site, first_cut, second_cut))
    return first, second, None


def lower_fomin_switch(records, include_two_sleeps=False):
    """Read the lower roots 0,2 with resident emissions held, and rewrite words."""
    creator, first_trace, _ = stopped_walk(records, 0, 0, {})
    offsets = {(0, z): len(word) for z, word in creator.items()}
    _, second_trace, _ = stopped_walk(records, 0, 2, offsets)
    first = [(z, b, d) for z, _, b, d in first_trace]
    second = [(z, b, d) for z, _, b, d in second_trace]
    changed_first, changed_second, cut = fomin_switch(first, second, include_two_sleeps)
    restored_first, restored_second, _ = fomin_switch(changed_first, changed_second,
                                                     include_two_sleeps)
    assert (restored_first, restored_second) == (first, second)
    changed, counts = dict(records), {}
    for site, bit, arrow in changed_first + changed_second:
        index = counts.get(site, 0) + 1
        counts[site] = index
        changed[0, site, index] = bit, arrow
    assert counts == dict(Counter(z for z, _, _ in first + second))
    return changed, cut


def reorder_lower_words(records, reverse=False):
    """A bijection between fixed reading orders (0,2) and (2,0).

    Retain each root's own stopped word. The inverse reads the other fixed
    order; no history-dependent order selection is being asserted.
    """
    read_order = (2, 0) if reverse else (0, 2)
    offsets, words = {}, {}
    for root in read_order:
        local, trace, _ = stopped_walk(records, 0, root, offsets)
        words[root] = [(z, bit, arrow) for z, _, bit, arrow in trace]
        for site, word in local.items():
            offsets[0, site] = offsets.get((0, site), 0) + len(word)
    changed, counts = dict(records), {}
    for root in reversed(read_order):
        for site, bit, arrow in words[root]:
            index = counts.get(site, 0) + 1
            counts[site] = index
            changed[0, site, index] = bit, arrow
    assert counts == {site: count for (copy, site), count in offsets.items()}
    return changed


def lower_cycle_switch(records, reverse=False):
    """Remove one cross-endpoint return cycle, retaining it in unused books.

    The inverse reads the two new root words, then the unused jump-only
    cycle from the first endpoint through its first return there.
    """
    offsets, words = {}, []
    for root in (0, 2):
        local, trace, _ = stopped_walk(records, 0, root, offsets)
        words.append([(z, bit, arrow) for z, _, bit, arrow in trace])
        for site, word in local.items():
            offsets[0, site] = offsets.get((0, site), 0) + len(word)
    first, second = words
    assert first[-1][1] == second[-1][1] == 1
    if not reverse:
        w, r = first[-1][0], second[-1][0]
        assert r != w and all(z != r for z, _, _ in second[:-1])
        first_cut = max(i for i, (z, _, _) in enumerate(first) if z == r)
        second_cut = next(i for i, (z, _, _) in enumerate(second) if z == w)
        assert all(z != w for z, _, _ in first[first_cut:-1])
        cycle = first[first_cut:-1] + second[second_cut:-1]
        assert all(bit == 0 for _, bit, _ in cycle)
        assert cycle[0][0] == r and cycle[-1][2] == r
        replacement = first[:first_cut] + second[-1:] + second[:second_cut] + first[-1:] + cycle
        before = first + second
    else:
        r, w = first[-1][0], second[-1][0]
        assert r != w and all(z != r and z != w for z, _, _ in second[:-1])
        cycle, site = [], r
        while True:
            index = offsets.get((0, site), 0) + 1
            offsets[0, site] = index
            bit, arrow = records[0, site, index]
            assert bit == 0 and arrow != -1
            cycle.append((site, bit, arrow))
            site = arrow
            if site == r:
                break
            assert len(cycle) < 100
        split = next(i for i, (z, _, _) in enumerate(cycle) if z == w)
        replacement = first[:-1] + cycle[:split] + second[-1:] + second[:-1] + cycle[split:] + first[-1:]
        before = first + second + cycle
    assert Counter(replacement) == Counter(before)
    changed, counts = dict(records), {}
    for site, bit, arrow in replacement:
        counts[site] = counts.get(site, 0) + 1
        changed[0, site, counts[site]] = bit, arrow
    return changed, cycle


def lower_cycle_involution(records):
    """Pair the cycle-removal domain with its disjoint, decodable image."""
    offsets, words = {}, []
    for root in (0, 2):
        local, trace, _ = stopped_walk(records, 0, root, offsets)
        words.append([(z, bit, arrow) for z, _, bit, arrow in trace])
        for site, word in local.items():
            offsets[0, site] = offsets.get((0, site), 0) + len(word)
    first, second = words
    if first[-1][1] != 1 or second[-1][1] != 1 or first[-1][0] == second[-1][0]:
        return dict(records), "fixed"
    w, r = first[-1][0], second[-1][0]
    cuts = [i for i, (z, _, _) in enumerate(first) if z == r]
    if (cuts and all(z != r for z, _, _ in second[:-1])
            and any(z == w for z, _, _ in second)
            and all(z != w for z, _, _ in first[cuts[-1]:-1])):
        return lower_cycle_switch(records)[0], "erase"
    # In the image the second path avoids the first endpoint and reaches
    # its own endpoint for the first time at its sleep.
    r, w = first[-1][0], second[-1][0]
    if any(z == r or z == w for z, _, _ in second[:-1]):
        return dict(records), "fixed"
    site, visited, steps = r, set(), 0
    while True:
        visited.add(site)
        index = offsets.get((0, site), 0)+1
        offsets[0, site] = index
        bit, arrow = records[0, site, index]
        if bit or arrow == -1:
            return dict(records), "fixed"
        site = arrow
        steps += 1
        assert steps < 100
        if site == r:
            if w in visited:
                return lower_cycle_switch(records, reverse=True)[0], "insert"
            return dict(records), "fixed"


def check_fomin_selector(include_two_sleeps=False):
    """Exhaust short recurrent paths to audit the incoming-edge cut convention."""
    families = []
    for source in (0, 2):
        paths = []
        for length in range(1, 5):
            for interior in product(range(3), repeat=length - 1):
                sites = (source,) + interior
                for sink in (False, True):
                    word = [(sites[j], 0, sites[j+1]) for j in range(length-1)]
                    word.append((sites[-1], 0, -1) if sink
                                else (sites[-1], 1, sites[-1]))
                    paths.append(word)
        families.append(paths)
    checked = switched = 0
    for first in families[0]:
        for second in families[1]:
            a, b, cut = fomin_switch(first, second, include_two_sleeps)
            c, d, _ = fomin_switch(a, b, include_two_sleeps)
            assert (c, d) == (first, second)
            assert Counter(a + b) == Counter(first + second)
            checked += 1
            switched += cut is not None
    return dict(path_pairs=checked, nontrivial_switches=switched,
                include_two_sleeps=include_two_sleeps,
                inverse_and_record_inventory=True)


def check_cycle_selector():
    """Exhaust short word pairs, retaining the entire newly unused cycle."""
    families = []
    for root in (0, 2):
        words = []
        for length in range(1, 5):
            for interior in product(range(3), repeat=length-1):
                sites = (root,) + interior
                prefix = [(z, 0, t) for z, t in zip(sites, sites[1:])]
                words.extend((prefix+[(sites[-1], 1, 0)], prefix+[(sites[-1], 0, -1)]))
        families.append(words)
    checked = erased = 0
    for first in families[0]:
        for second in families[1]:
            records = {(0, z, i): (1, -1) for z in range(3) for i in range(1, 12)}
            counts = {}
            for site, bit, arrow in first+second:
                counts[site] = counts.get(site, 0)+1
                records[0, site, counts[site]] = bit, arrow
            changed, mode = lower_cycle_involution(records)
            restored, reverse_mode = lower_cycle_involution(changed)
            assert restored == records
            if mode == "erase":
                assert reverse_mode == "insert"
                book = {(z, i): record for (copy, z, i), record in changed.items()}
                assert stabilize(book, [0, 2])[0] == {first[-1][0], second[-1][0]}
                erased += 1
            else:
                assert mode == reverse_mode == "fixed"
            original_inventory = Counter((z, bit, arrow)
                for (copy, z, i), (bit, arrow) in records.items())
            new_inventory = Counter((z, bit, arrow)
                for (copy, z, i), (bit, arrow) in changed.items())
            assert original_inventory == new_inventory
            checked += 1
    return dict(path_pairs=checked, erased_and_reinserted=erased,
                complete_inverse=True, site_record_inventory=True)


def check_contact_then_fomin(two_sleep_endpoints=False):
    """The first exchange gives (1,1); the loop-erased second exchange repairs it."""
    kernel = {0: {z: F(1, 5) for z in (0, 1, 2, 3, -1)},
              1: {0: F(1, 2), 2: F(1, 2)},
              2: {z: F(1, 3) for z in (1, 3, -1)},
              3: {z: F(1, 3) for z in (0, 3, -1)}}
    sleeps = {0: F(2, 5), 1: F(3, 7), 2: F(5, 11), 3: F(7, 13)}
    if two_sleep_endpoints:
        kernel[3] = {z: F(1, 4) for z in (0, 3, 4, -1)}
        kernel[4] = {-1: F(1)}
        sleeps[4] = F(11, 17)
    checks = []
    for lower_loops, creator_loops, explorer_loops, return_loops in (
            (0, 0, 0, 0), (2, 1, 3, 2), (5, 3, 1, 4)):
        records = {(c, z, i): (0, 2 if z == 1 else -1)
                   for c in (0, 1) for z in kernel for i in range(1, 41)}
        for i in range(1, lower_loops + 1):
            records[0, 0, i] = 0, 0
        records[0, 0, lower_loops + 1] = 1, 3
        records[0, 0, lower_loops + 2] = 1, -1
        records[0, 2, 1] = 0, 3
        for i in range(1, return_loops + 1):
            records[0, 3, i] = 0, 3
        records[0, 3, return_loops + 1] = 0, 4 if two_sleep_endpoints else -1
        if two_sleep_endpoints:
            records[0, 4, 1] = 1, -1
        for i in range(1, creator_loops + 1):
            records[1, 0, i] = 0, 0
        records[1, 0, creator_loops + 1] = 1, 2
        records[1, 1, 1] = 1, 2
        records[1, 2, 1] = 0, 1
        records[1, 1, 2] = 0, 0
        for i in range(creator_loops + 2, creator_loops + 2 + explorer_loops):
            records[1, 0, i] = 0, 0
        records[1, 0, creator_loops + 2 + explorer_loops] = 0, 3
        records[1, 3, 1] = 1, 0
        records[1, 2, 2] = 1, -1
        records[1, 2, 3] = 0, -1
        offsets, _ = creator_contact_offsets(records)
        intermediate, _ = switch_walk_prefixes(records, 0, offsets)
        changed, cut = lower_fomin_switch(intermediate, two_sleep_endpoints)
        recovered_intermediate, _ = lower_fomin_switch(changed, two_sleep_endpoints)
        assert recovered_intermediate == intermediate
        reverse_offsets, _ = creator_contact_offsets(recovered_intermediate)
        restored, _ = switch_walk_prefixes(recovered_intermediate, 0, reverse_offsets)
        assert restored == records
        endpoints = [[sorted(staged_creator_pair(field, copy)[0]) for copy in (0, 1)]
                     for field in (records, intermediate, changed)]
        lower_before = [0, 4] if two_sleep_endpoints else [0]
        lower_after = [3, 4] if two_sleep_endpoints else [3]
        assert endpoints == [[lower_before, [3]], [lower_before, [0]], [lower_after, [0]]]
        assert record_weight(records, kernel, sleeps) == record_weight(intermediate, kernel, sleeps)
        assert record_weight(records, kernel, sleeps) == record_weight(changed, kernel, sleeps) > 0
        for field, expected in zip((records, intermediate, changed), endpoints):
            for copy in (0, 1):
                book = {(z, i): r for (side, z, i), r in field.items() if side == copy}
                assert sorted(stabilize(book, [0, 1, 2] if copy else [0, 2])[0]) == expected[copy]
        reordered = reorder_lower_words(intermediate)
        assert lower_cycle_involution(intermediate)[0] == intermediate
        assert reorder_lower_words(reordered, reverse=True) == intermediate
        assert record_weight(reordered, kernel, sleeps) == record_weight(records, kernel, sleeps)
        for copy in (0, 1):
            book = {(z, i): r for (side, z, i), r in reordered.items() if side == copy}
            assert sorted(stabilize(book, [0, 1, 2] if copy else [0, 2])[0]) == endpoints[-1][copy]
        checks.append(dict(loops=[lower_loops, creator_loops, explorer_loops, return_loops],
                           loop_erased_cut=cut, endpoints=endpoints,
                           inverse_in_reverse_order=True, equal_weight=True,
                           fixed_root_order_repair=True, fixed_root_order_inverse=True))
    return checks


def check_two_sleep_obstruction():
    """A return to the second endpoint can survive the loop-erased repair."""
    kernel = {z: {d: F(1, 7) for d in (*range(6), -1)} for z in range(6)}
    sleeps = dict(enumerate(map(F, ("2/5", "3/7", "5/11", "7/13", "11/17", "13/19"))))
    records = {(0, z, i): (0, -1) for z in range(6) for i in range(1, 8)}
    words = {0: [(0, 5), (1, -1)], 2: [(0, 4)],
             3: [(1, 0), (0, 4)], 4: [(0, 5), (1, 0)],
             5: [(0, 3), (0, 3)]}
    for site, word in words.items():
        for index, record in enumerate(word, 1):
            records[0, site, index] = record
    changed, cut = lower_fomin_switch(records, True)
    restored, _ = lower_fomin_switch(changed, True)
    assert restored == records
    assert record_weight(records, kernel, sleeps) == record_weight(changed, kernel, sleeps) > 0
    endpoints = []
    for field in (records, changed):
        book = {(z, i): record for (side, z, i), record in field.items()}
        endpoints.append(sorted(stabilize(book, [0, 2])[0]))
    assert endpoints == [[0, 4], [0, 3]]
    reordered = reorder_lower_words(records)
    assert reorder_lower_words(reordered, reverse=True) == records
    assert record_weight(reordered, kernel, sleeps) == record_weight(records, kernel, sleeps)
    book = {(z, i): record for (side, z, i), record in reordered.items()}
    assert stabilize(book, [0, 2])[0] == {3, 4}
    return dict(first_word_sites=[0, 5, 3], second_word_sites=[2, 4, 5, 3, 4],
                loop_erased_cut=cut, final_supports=endpoints,
                target_still_occupied=True, inverse=True, equal_weight=True,
                fixed_root_order_repairs_to=[3, 4])


def check_mutual_sleep_dependencies():
    """Mutual endpoint visits: one erased-loop case and one remaining cycle."""
    kernel = {z: {d: F(1, 7) for d in (*range(6), -1)} for z in range(6)}
    sleeps = dict(enumerate(map(F, ("2/5", "3/7", "5/11", "7/13", "11/17", "13/19"))))
    cases = [
        ([(0, 0, 4), (4, 0, 0), (0, 0, 5), (5, 0, 3), (3, 1, 0)],
         [(2, 0, 5), (5, 0, 3), (3, 0, 4), (4, 1, 0)], {3, 4}),
        ([(0, 0, 4), (4, 0, 3), (3, 1, 0)],
         [(2, 0, 3), (3, 0, 4), (4, 1, 0)], {0, 3}),
    ]
    result = []
    for first, second, fomin_endpoint in cases:
        records = {(0, z, i): (0, -1) for z in range(6) for i in range(1, 12)}
        counts = {}
        for site, bit, arrow in first + second:
            counts[site] = counts.get(site, 0) + 1
            records[0, site, counts[site]] = bit, arrow
        records[0, 0, counts[0] + 1] = 1, -1
        changed, cut = lower_fomin_switch(records, True)
        reordered = reorder_lower_words(records)
        assert lower_fomin_switch(changed, True)[0] == records
        assert reorder_lower_words(reordered, reverse=True) == records
        assert record_weight(records, kernel, sleeps) == record_weight(changed, kernel, sleeps)
        assert record_weight(records, kernel, sleeps) == record_weight(reordered, kernel, sleeps)
        endpoints = []
        for field in (records, changed, reordered):
            book = {(z, i): record for (side, z, i), record in field.items()}
            endpoints.append(stabilize(book, [0, 2])[0])
        assert endpoints == [{0, 4}, fomin_endpoint, {0, 3}]
        result.append(dict(first_sites=[z for z, _, _ in first],
                           second_sites=[z for z, _, _ in second], cut=cut,
                           original_fomin_reordered=[sorted(s) for s in endpoints],
                           full_inverses=True, equal_weights=True))
    return result


def check_contact_then_cycle():
    """Full bad-to-good pairs with a retained, inverse-readable unused cycle."""
    kernel = {z: {d: F(1, 7) for d in (*range(6), -1)} for z in range(6)}
    sleeps = dict(enumerate(map(F, ("2/5", "3/7", "5/11", "7/13", "11/17", "13/19"))))
    checks = []
    for loops in (0, 2, 5):
        W_sites = [0, 4] + ([0, 4] if loops == 5 else []) + [5]*(loops+1) + [3]
        R_sites = [2, 3] + [5, 3]*loops + [4]
        W = [(z, 0, t) for z, t in zip(W_sites, W_sites[1:])] + [(3, 1, 0)]
        R = [(z, 0, t) for z, t in zip(R_sites, R_sites[1:])] + [(4, 1, 0)]
        records = {(c, z, i): (0, -1)
                   for c in (0, 1) for z in range(6) for i in range(1, 45)}
        records[0, 0, 1] = 1, -1
        records[0, 0, 2] = 1, -1
        counts = {0: 1}
        for site, bit, arrow in R:
            counts[site] = counts.get(site, 0)+1
            records[0, site, counts[site]] = bit, arrow
        records[1, 0, 1] = 1, -1
        records[1, 1, 1] = 1, -1
        records[1, 2, 1] = 0, 0
        counts = {0: 1, 1: 1, 2: 1}
        for site, bit, arrow in W:
            counts[site] = counts.get(site, 0)+1
            records[1, site, counts[site]] = bit, arrow
        offsets, _ = creator_contact_offsets(records)
        intermediate, _ = switch_walk_prefixes(records, 0, offsets)
        changed, cycle = lower_cycle_switch(intermediate)
        assert lower_cycle_involution(intermediate) == (changed, "erase")
        assert lower_cycle_involution(changed) == (intermediate, "insert")
        recovered, decoded_cycle = lower_cycle_switch(changed, reverse=True)
        assert recovered == intermediate and decoded_cycle == cycle
        reverse_offsets, _ = creator_contact_offsets(recovered)
        restored, _ = switch_walk_prefixes(recovered, 0, reverse_offsets)
        assert restored == records
        assert record_weight(records, kernel, sleeps) == record_weight(intermediate, kernel, sleeps)
        assert record_weight(records, kernel, sleeps) == record_weight(changed, kernel, sleeps) > 0
        endpoints = [[sorted(staged_creator_pair(field, copy)[0]) for copy in (0, 1)]
                     for field in (records, intermediate, changed)]
        assert endpoints == [[[0, 4], [1, 3]], [[0, 4], [0, 1]], [[3, 4], [0, 1]]]
        for field, expected in zip((records, intermediate, changed), endpoints):
            for copy in (0, 1):
                book = {(z, i): record for (c, z, i), record in field.items() if c == copy}
                assert sorted(stabilize(book, [0, 1, 2] if copy else [0, 2])[0]) == expected[copy]
        combined = reorder_lower_words(changed)
        assert lower_cycle_involution(reorder_lower_words(combined, reverse=True))[0] == intermediate
        assert record_weight(combined, kernel, sleeps) == record_weight(records, kernel, sleeps)
        for copy in (0, 1):
            book = {(z, i): record for (c, z, i), record in combined.items() if c == copy}
            assert sorted(stabilize(book, [0, 1, 2] if copy else [0, 2])[0]) == endpoints[-1][copy]
        checks.append(dict(loops=loops, endpoints=endpoints, retained_cycle=cycle,
                           decoded_first_return=True, complete_inverse=True, equal_weights=True,
                           fixed_cycle_then_root_order_repairs=True))
    return checks


def last_target_offsets(records):
    """Recover the upper FIFO prefix before its final arrival at target 0."""
    book = {(z, i): record for (copy, z, i), record in records.items() if copy == 1}
    trace = stabilize(book, [0, 1, 2])[3]
    cut = max(i for i, (site, _, _, _) in enumerate(trace) if site == 0)
    offsets, residents = {}, {}
    for site, index, bit, _ in trace[:cut]:
        offsets[1, site] = index
        residents.pop(site, None)
        if bit:
            residents[site] = records[1, site, index][1]
    return offsets, residents, trace[cut][1]


def check_last_contact_creator():
    """A final-arrival cut repairs pairs with earlier upper target returns."""
    kernel = {z: {d: F(1, 6) for d in (*range(5), -1)} for z in range(5)}
    sleeps = dict(enumerate(map(F, ("2/5", "3/7", "5/11", "7/13", "11/17"))))
    checks = []
    for loops in (0, 2, 5):
        records = {(copy, z, i): (0, -1)
                   for copy in (0, 1) for z in range(5) for i in range(1, 35)}
        records[0, 0, 1] = records[0, 0, 2] = 1, -1
        records[0, 2, 1] = 0, 4
        records[1, 0, 1] = 1, 0
        records[1, 1, 1] = 1, -1
        records[1, 2, 1] = 0, 0
        for index in range(2, loops+2):
            records[1, 0, index] = 0, 0
        records[1, 0, loops+2] = 0, 3
        records[1, 0, loops+3] = 0, 4
        records[1, 3, 1] = 1, -1
        records[1, 4, 1] = 1, 0
        offsets, residents, last_index = last_target_offsets(records)
        _, word, endpoint = stopped_walk(records, 1, 0, offsets)
        assert endpoint == 4 and not ({z for z, _, _, _ in word} & (set(residents)-{0}))
        intermediate, _ = switch_walk_prefixes(records, 0, offsets)
        changed = reorder_lower_words(intermediate)
        restored_intermediate = reorder_lower_words(changed, reverse=True)
        reverse_offsets, _, reverse_index = last_target_offsets(restored_intermediate)
        assert (reverse_offsets, reverse_index) == (offsets, last_index)
        restored, _ = switch_walk_prefixes(restored_intermediate, 0, reverse_offsets)
        assert restored == records
        early_offsets, _ = creator_contact_offsets(records)
        early, _ = switch_walk_prefixes(records, 0, early_offsets)
        endpoints = []
        for field in (records, early, intermediate, changed):
            endpoints.append([])
            for copy in (0, 1):
                book = {(z, i): record for (c, z, i), record in field.items() if c == copy}
                endpoints[-1].append(sorted(stabilize(book, [0, 1, 2] if copy else [0, 2])[0]))
            assert record_weight(field, kernel, sleeps) == record_weight(records, kernel, sleeps) > 0
        assert endpoints == [[[0], [1, 3, 4]], [[3], [1, 4]],
                             [[0], [0, 1, 3]], [[4], [0, 1, 3]]]
        checks.append(dict(loops=loops, final_upper_target_index=last_index,
                           original_early_last_reordered=endpoints,
                           recovered_last_cut=True, complete_inverse=True, equal_weights=True))
    return checks


def marked_target_offsets(records):
    """The sleep immediately preceding the final target record marks a cut."""
    book = {(z, i): record for (copy, z, i), record in records.items() if copy == 1}
    trace = stabilize(book, [0, 1, 2])[3]
    last_position = max(i for i, (site, _, _, _) in enumerate(trace) if site == 0)
    previous_sleep = max((index for site, index, bit, _ in trace[:last_position]
                          if site == 0 and bit), default=0)
    cut_index = previous_sleep+1
    cut = next(i for i, (site, index, _, _) in enumerate(trace)
               if site == 0 and index == cut_index)
    offsets, residents = {}, {}
    for site, index, bit, _ in trace[:cut]:
        offsets[1, site] = index
        residents.pop(site, None)
        if bit:
            residents[site] = records[1, site, index][1]
    return offsets, residents, cut_index, trace[last_position][1]


def check_marked_returning_creator():
    """A preceding sleep recovers the cut even when lower L returns to 0."""
    kernel = {z: {d: F(1, 6) for d in (*range(5), -1)} for z in range(5)}
    sleeps = dict(enumerate(map(F, ("2/5", "3/7", "5/11", "7/13", "11/17"))))
    checks = []
    for last_index, lower_loops in ((3, 2), (5, 1), (8, 4)):
        records = {(copy, z, i): (0, -1)
                   for copy in (0, 1) for z in range(5) for i in range(1, 35)}
        for index in range(1, lower_loops+1):
            records[0, 0, index] = 0, 0
        if last_index == 5:
            records[0, 0, 1] = 0, 4
            records[0, 4, 1] = 0, 0
        records[0, 0, lower_loops+1] = records[0, 0, lower_loops+2] = 1, -1
        records[0, 2, 1] = 0, 4
        for index in range(1, last_index-1):
            records[1, 0, index] = 1, 0
        records[1, 0, last_index-1] = 1, 3
        records[1, 0, last_index] = 0, 4
        records[1, 1, 1] = 1, -1
        records[1, 2, 1] = 0, 0
        records[1, 3, 1] = 1, -1
        records[1, 4, 1] = 1, 0
        offsets, residents, cut_index, old_final_index = marked_target_offsets(records)
        assert cut_index == old_final_index == last_index and 0 in residents
        intermediate, _ = switch_walk_prefixes(records, 0, offsets)
        changed = reorder_lower_words(intermediate)
        recovered = reorder_lower_words(changed, reverse=True)
        reverse_offsets, _, reverse_cut, new_final_index = marked_target_offsets(recovered)
        assert (reverse_offsets, reverse_cut) == (offsets, cut_index)
        assert new_final_index == last_index+lower_loops
        restored, _ = switch_walk_prefixes(recovered, 0, reverse_offsets)
        assert restored == records
        endpoints = []
        for field in (records, intermediate, changed):
            endpoints.append([])
            for copy in (0, 1):
                book = {(z, i): record for (c, z, i), record in field.items() if c == copy}
                endpoints[-1].append(sorted(stabilize(book, [0, 1, 2] if copy else [0, 2])[0]))
            assert record_weight(field, kernel, sleeps) == record_weight(records, kernel, sleeps) > 0
        assert endpoints == [[[0], [1, 3, 4]], [[0], [0, 1, 3]], [[4], [0, 1, 3]]]
        checks.append(dict(lower_creator_returns=lower_loops, old_last_index=last_index,
                           new_last_index=new_final_index, recovered_cut=reverse_cut,
                           off_site_return=last_index == 5,
                           endpoints=endpoints, complete_inverse=True, equal_weights=True))
    return checks


def creation_prefix(records, copy, roots):
    """Stable support, read offsets, and terminal resident-arrow coordinates."""
    book = {(z, i): record for (c, z, i), record in records.items() if c == copy}
    support, _, _, trace = stabilize(book, roots)
    offsets = {}
    for site, index, _, _ in trace:
        offsets[copy, site] = index
    resident_keys = {z: (copy, z, offsets[copy, z]) for z in support}
    return support, offsets, resident_keys


def avalanche_block(records, copy, offsets, resident_keys, source):
    """Read a complete avalanche; initial resident arrows are separate data."""
    book = {(z, i-offsets.get((copy, z), 0)): record
            for (c, z, i), record in records.items()
            if c == copy and i > offsets.get((copy, z), 0)}
    residents = {z: ((z, 0), records[key][1]) for z, key in resident_keys.items()}
    support, _, _, trace = stabilize(book, [] if source == -1 else [source], residents)
    counts = {}
    for site, index, _, _ in trace:
        counts[site] = index
    words = {z: [book[z, i] for i in range(1, count+1)] for z, count in counts.items()}
    return words, support, trace


def matched_terminal_avalanche(records, lower_roots, upper_roots, n, m,
                               reserve_missing=False):
    """The counterfactual lower sleep recovers the full creation background."""
    lower_target, upper_target = (0, 0, n), (1, 0, m)
    assert records[lower_target][0] == 0 and records[upper_target][0] == 1
    source = records[lower_target][1]
    counterfactual = dict(records)
    counterfactual[lower_target] = records[upper_target]
    stable, lower_offsets, lower_keys = creation_prefix(counterfactual, 0, lower_roots)
    assert 0 in stable and lower_offsets[0, 0] == n
    B = stable-{0}
    sigma, upper_offsets, upper_keys = creation_prefix(records, 1, upper_roots)
    assert 0 in sigma
    C = B & sigma
    if reserve_missing:
        for z in B-sigma:
            j = upper_offsets.get((1, z), 0) + 1
            while records[1, z, j][0] != 1:
                j += 1
            upper_keys[z] = (1, z, j)
            upper_offsets[1, z] = j
        C = B
    lower_keys = {z: lower_keys[z] for z in C}
    selected_upper_keys = {z: upper_keys[z] for z in C}
    offsets = lower_offsets | upper_offsets
    lower = avalanche_block(records, 0, offsets, lower_keys, source)
    upper = avalanche_block(records, 1, offsets, selected_upper_keys, source)
    changed = exchange_record_blocks(records, [lower[0], upper[0]], offsets)
    for z in C:
        lo, hi = lower_keys[z], selected_upper_keys[z]
        lower_bit, lower_arrow = changed[lo]
        upper_bit, upper_arrow = changed[hi]
        assert lower_bit == upper_bit == 1
        changed[lo], changed[hi] = (1, upper_arrow), (1, lower_arrow)
    return changed, B, sigma, lower, upper


def matched_avalanche_switch(records, common_roots):
    """Special case with a fixed lower target index one."""
    _, _, m = last_target_offsets_general(records, common_roots)
    return matched_terminal_avalanche(records, [0]+list(common_roots),
                                      [0, 1]+list(common_roots), 1, m)


def terminal_pair_repair(records, lower_roots, upper_roots, reverse=False,
                         reserve_missing=False):
    """Both terminal indices are decoded from the final image, without marks."""
    lower, lower_offsets, _ = creation_prefix(records, 0, lower_roots)
    upper, upper_offsets, _ = creation_prefix(records, 1, upper_roots)
    assert (0 in lower, 0 in upper) == ((False, True) if reverse else (True, False))
    n, m = lower_offsets[0, 0], upper_offsets[1, 0]
    lo, hi = (0, 0, n), (1, 0, m)
    if reverse:
        changed, B, sigma, _, _ = matched_terminal_avalanche(
            records, lower_roots, upper_roots, n, m, reserve_missing)
        changed[lo], changed[hi] = changed[hi], changed[lo]
    else:
        intermediate = dict(records)
        intermediate[lo], intermediate[hi] = intermediate[hi], intermediate[lo]
        changed, B, sigma, _, upper_aux = matched_terminal_avalanche(
            intermediate, lower_roots, upper_roots, n, m, reserve_missing)
        forbidden = {0} if reserve_missing else ({0} | (B-sigma))
        assert all(z not in forbidden for z, _, _, _ in upper_aux[2])
        assert creation_prefix(changed, 0, lower_roots)[1][0, 0] == n
        assert creation_prefix(changed, 1, upper_roots)[1][1, 0] == m
    return changed, dict(n=n, m=m, background=sorted(B), upper_stable=sorted(sigma))


def check_repaired_terminal_overlap():
    """Whole avalanches separate the two old overlapping terminal-swap images."""
    kernel, sleeps, old_image, sources = check_terminal_overlap(return_fields=True)
    images, checks = [], []
    for source in sources:
        changed, cuts = terminal_pair_repair(source, [0, 2], [0, 1, 2])
        restored, decoded = terminal_pair_repair(changed, [0, 2], [0, 1, 2], reverse=True)
        assert restored == source and decoded == cuts
        assert record_weight(changed, kernel, sleeps) == record_weight(source, kernel, sleeps)
        lower, _, _ = creation_prefix(changed, 0, [0, 2])
        upper, _, _ = creation_prefix(changed, 1, [0, 1, 2])
        assert 0 not in lower and 0 in upper
        images.append(changed)
        checks.append(dict(cuts=cuts, final_supports=[sorted(lower), sorted(upper)],
                           full_inverse=True, equal_weight=True))
    assert images[0] != images[1]
    return dict(old_terminal_swap_images_coincide=True, new_images_are_distinct=True,
                source_and_decoded_indices=[[1, 1], [2, 1]], checks=checks)


def check_matched_avalanches():
    """Complete wake cascades, including equal endpoints and two sleepers."""
    kernel = {z: {d: F(1, 7) for d in (*range(6), -1)} for z in range(6)}
    sleeps = dict(enumerate(map(F, ("2/5", "3/7", "5/11", "7/13", "11/17", "13/19"))))
    checks = []
    for common_roots in ((2,), (2, 4)):
        records = {(copy, z, i): (0, -1)
                   for copy in (0, 1) for z in range(6) for i in range(1, 45)}
        for index in range(1, 8):
            records[0, 0, index] = 1, -1
        for z in common_roots:
            records[0, z, 1] = 1, 0
            records[1, z, 1] = 1, -1
        records[0, 3, 1] = records[1, 3, 1] = 0, 2
        if len(common_roots) == 1:
            records[0, 2, 2] = records[1, 2, 2] = 1, 0
            outside_sleeper = 4
        else:
            records[0, 2, 2] = records[1, 2, 2] = 0, 4
            records[0, 4, 2] = records[1, 4, 2] = 1, 0
            outside_sleeper = 5
        records[1, 0, 1] = 0, outside_sleeper
        records[1, 1, 1] = 0, 0
        records[1, 0, 2] = 0, 3
        records[1, outside_sleeper, 1] = 1, -1
        # First repair the upper last target bit, at fixed lower index one.
        _, _, last_index = last_target_offsets_general(records, common_roots)
        intermediate = dict(records)
        lo, hi = (0, 0, 1), (1, 0, last_index)
        intermediate[lo], intermediate[hi] = intermediate[hi], intermediate[lo]
        changed, B, sigma, lower, upper = matched_avalanche_switch(intermediate, common_roots)
        assert all(z != 0 for z, _, _, _ in upper[2])
        # The upper auxiliary avalanche is dominated by adding the same
        # packet to its entire actual stable background.
        _, upper_offsets, upper_keys = creation_prefix(intermediate, 1, [0, 1]+list(common_roots))
        full_upper = avalanche_block(intermediate, 1, upper_offsets, upper_keys, 3)
        assert all(z != 0 for z, _, _, _ in full_upper[2])
        assert {z for z, _, _, _ in upper[2]} <= {z for z, _, _, _ in full_upper[2]}
        restored, *_ = matched_avalanche_switch(changed, common_roots)
        assert restored == intermediate
        _, _, image_last = last_target_offsets_general(restored, common_roots)
        assert image_last == last_index
        restored[lo], restored[hi] = restored[hi], restored[lo]
        assert restored == records
        endpoints = []
        for field in (records, intermediate, changed):
            endpoints.append([])
            for copy in (0, 1):
                roots = ([0, 1] if copy else [0])+list(common_roots)
                book = {(z, i): record for (c, z, i), record in field.items() if c == copy}
                endpoints[-1].append(sorted(stabilize(book, roots)[0]))
            assert record_weight(field, kernel, sleeps) == record_weight(records, kernel, sleeps) > 0
        assert 0 in endpoints[0][0] and 0 not in endpoints[0][1]
        assert 0 in endpoints[1][0] and 0 in endpoints[1][1]
        assert 0 not in endpoints[2][0] and 0 in endpoints[2][1]
        checks.append(dict(common_roots=common_roots, matched_background=sorted(B),
                           upper_stable=sorted(sigma), endpoints=endpoints,
                           lower_avalanche_before=sorted(lower[1]),
                           transplanted_avalanche=sorted(upper[1]),
                           unused_resident_arrows_exchanged=True,
                           complete_inverse=True, equal_weights=True))
    return checks


def last_target_offsets_general(records, common_roots):
    """Only the final target index is needed for the general witness inputs."""
    book = {(z, i): record for (c, z, i), record in records.items() if c == 1}
    trace = stabilize(book, [0, 1]+list(common_roots))[3]
    last_index = max(index for z, index, _, _ in trace if z == 0)
    return None, None, last_index


def check_completed_word_contacts():
    """Both the final word and its residual can visit an old common sleeper."""
    kernel = {z: {d: F(1, 6) for d in (*range(5), -1)} for z in range(5)}
    sleeps = dict(enumerate(map(F, ("2/5", "3/7", "5/11", "7/13", "11/17"))))
    checks = []
    for returning in (False, True):
        records = {(copy, z, i): (0, -1)
                   for copy in (0, 1) for z in range(5) for i in range(1, 30)}
        records[0, 0, 1] = 1, -1
        records[0, 2, 1] = 1, 0
        if returning:
            records[0, 0, 1] = 0, 2
            records[0, 0, 2] = 1, -1
            records[0, 2, 1] = 0, 0
            records[0, 2, 2] = 1, 0
        records[1, 0, 1] = 1, 3
        records[1, 1, 1] = 0, 0
        records[1, 2, 1] = 1, -1
        records[1, 0, 2] = 0, 3
        records[1, 3, 1] = 0, 2
        records[1, 2, 2] = 0, 4
        records[1, 4, 1] = 1, -1
        records[1, 3, 2] = 0, 2
        records[1, 2, 3] = 1, -1
        offsets, residents, last_index = last_target_offsets(records)
        assert set(residents) == {0, 2} and last_index == 2
        changed, walks = switch_walk_prefixes(records, 0, offsets)
        touched = [{z for z, _, _, _ in walk[1]} & (set(residents)-{0}) for walk in walks]
        assert touched == ([{2}, {2}] if returning else [set(), {2}])
        if returning:
            reverse_offsets, _, cut, final_index = marked_target_offsets(changed)
            assert cut == 2 and final_index == 3
        else:
            reverse_offsets, _, final_index = last_target_offsets(changed)
            assert final_index == 2
        assert reverse_offsets == offsets
        restored, _ = switch_walk_prefixes(changed, 0, reverse_offsets)
        assert restored == records
        assert record_weight(changed, kernel, sleeps) == record_weight(records, kernel, sleeps) > 0
        endpoints = []
        for field in (records, changed):
            endpoints.append([])
            for copy in (0, 1):
                book = {(z, i): record for (c, z, i), record in field.items() if c == copy}
                endpoints[-1].append(sorted(stabilize(book, [0, 1, 2] if copy else [0, 2])[0]))
        assert endpoints == [[[0, 2], [2, 4]], [[2, 4], [0, 2]]]
        checks.append(dict(returning_lower_creator=returning, endpoints=endpoints,
                           contacted_background=[sorted(s) for s in touched],
                           complete_inverse=True, equal_weights=True))
    return checks


def check_terminal_repair_samples():
    """Fixed-seed finite cylinders audit the general inverse beyond hand cases.

    Every site has 64 sampled records followed by eight sink jumps, so these
    executions terminate without discarding long histories or truncating them.
    Site-specific combined record inventories certify equal product weights.
    """
    rng = Random(20261010)
    p = [F(2, 5), F(3, 7), F(5, 11), F(7, 13)]
    destinations = [0, 1, 2, 3, -1]
    bad = matched = extended = lower_returns = upper_returns = 0
    backgrounds = Counter()
    for sample in range(512):
        lower_roots = [0, 2] if sample % 2 else [0, 2, 3]
        upper_roots = [0, 1]+lower_roots[1:]
        records = {}
        for copy in (0, 1):
            for site in range(4):
                for index in range(1, 73):
                    records[copy, site, index] = (
                        (int(rng.randrange(p[site].denominator) < p[site].numerator),
                         rng.choice(destinations)) if index <= 64 else (0, -1))
        lower, lower_offsets, _ = creation_prefix(records, 0, lower_roots)
        upper, upper_offsets, _ = creation_prefix(records, 1, upper_roots)
        if not (0 in lower and 0 not in upper):
            continue
        bad += 1
        n, m = lower_offsets[0, 0], upper_offsets[1, 0]
        probe = dict(records)
        lo, hi = (0, 0, n), (1, 0, m)
        probe[lo], probe[hi] = probe[hi], probe[lo]
        sigma, offsets, keys = creation_prefix(probe, 1, upper_roots)
        B = lower-{0}
        selected = {z: keys[z] for z in B & sigma}
        upper_aux = avalanche_block(probe, 1, offsets, selected, probe[lo][1])
        assert all(z != 0 for z, _, _, _ in upper_aux[2])
        if any(z in B-sigma for z, _, _, _ in upper_aux[2]):
            continue
        changed, cuts = terminal_pair_repair(records, lower_roots, upper_roots)
        restored, decoded = terminal_pair_repair(changed, lower_roots, upper_roots, reverse=True)
        assert restored == records and cuts == decoded
        original_inventory = Counter((z, bit, arrow)
            for (copy, z, index), (bit, arrow) in records.items())
        changed_inventory = Counter((z, bit, arrow)
            for (copy, z, index), (bit, arrow) in changed.items())
        assert original_inventory == changed_inventory
        matched += B <= sigma
        extended += not B <= sigma
        lower_returns += n > 1
        upper_returns += m > 1
        backgrounds[len(lower)-1] += 1
    assert matched and lower_returns and upper_returns and backgrounds[2]
    return dict(seed=20261010, sampled_cylinders=512, bad_pairs=bad,
                matched_pairs_repaired=matched, unmatched_pairs_repaired=extended,
                outside_extended_domain=bad-matched-extended,
                lower_noninitial_terminal=lower_returns, upper_noninitial_terminal=upper_returns,
                background_sizes=dict(backgrounds), all_full_inverses=True,
                all_site_record_inventories=True, no_discarded_long_histories=True)


def check_reserved_terminal_samples():
    """A fixed first-unused-sleep reservation permits actual unmatched wakes."""
    rng = Random(20261010)
    p = [F(2, 5), F(3, 7), F(5, 11), F(7, 13)]
    counts = Counter()
    for sample in range(512):
        lower_roots = [0, 2] if sample % 2 else [0, 2, 3]
        upper_roots = [0, 1]+lower_roots[1:]
        records = {}
        for copy in (0, 1):
            for site in range(4):
                for index in range(1, 82):
                    records[copy, site, index] = (
                        (int(rng.randrange(p[site].denominator) < p[site].numerator),
                         rng.choice([0, 1, 2, 3, -1])) if index <= 64 else
                        (int(index == 73), -1))
        lower, low_offsets, _ = creation_prefix(records, 0, lower_roots)
        upper, upp_offsets, _ = creation_prefix(records, 1, upper_roots)
        if not (0 in lower and 0 not in upper):
            continue
        counts['bad_pairs'] += 1
        n, m = low_offsets[0, 0], upp_offsets[1, 0]
        probe = dict(records)
        lo, hi = (0, 0, n), (1, 0, m)
        probe[lo], probe[hi] = probe[hi], probe[lo]
        _, B, sigma, _, aux = matched_terminal_avalanche(
            probe, lower_roots, upper_roots, n, m, reserve_missing=True)
        if any(z == 0 for z, _, _, _ in aux[2]):
            counts['auxiliary_target_hits'] += 1
            continue
        changed, cuts = terminal_pair_repair(
            records, lower_roots, upper_roots, reserve_missing=True)
        restored, decoded = terminal_pair_repair(
            changed, lower_roots, upper_roots, reverse=True, reserve_missing=True)
        assert restored == records and cuts == decoded
        inventory = lambda field: Counter((z, bit, arrow)
            for (_, z, _), (bit, arrow) in field.items())
        assert inventory(changed) == inventory(records)
        counts['repaired_pairs'] += 1
        counts['unmatched_backgrounds'] += bool(B-sigma)
        counts['unmatched_residents_woken'] += any(z in B-sigma for z, _, _, _ in aux[2])
    assert counts['unmatched_residents_woken']
    return dict(seed=20261010, sampled_cylinders=512, **counts,
                all_full_inverses=True, all_record_inventories=True,
                no_discarded_long_histories=True)


def check_forced_auxiliary_return():
    """A missing sleeper with its only arrow to v defeats every reservation."""
    records = {(copy, site, index): (0, 0 if site == 1 else -1)
               for copy in (0, 1) for site in range(3) for index in range(1, 12)}
    records[0, 0, 1] = (1, -1)
    records[0, 1, 1] = (1, 0)
    records[1, 0, 2] = (0, 1)
    records[1, 1, 2] = (1, 0)
    lr, ur = [0, 1], [0, 1, 2]
    lower, loff, _ = creation_prefix(records, 0, lr)
    upper, uoff, _ = creation_prefix(records, 1, ur)
    assert lower == {0, 1} and upper == {1}
    n, m = loff[0, 0], uoff[1, 0]
    probe = dict(records)
    probe[0, 0, n], probe[1, 0, m] = probe[1, 0, m], probe[0, 0, n]
    _, B, sigma, _, auxiliary = matched_terminal_avalanche(
        probe, lr, ur, n, m, reserve_missing=True)
    assert B == {1} and sigma == {0}
    assert auxiliary[2][0][0] == 1 and 0 in auxiliary[2][0][3]
    assert any(site == 0 for site, _, _, _ in auxiliary[2])
    kernel = {0: {1: F(1, 2), -1: F(1, 2)}, 1: {0: F(1)},
              2: {0: F(1, 2), -1: F(1, 2)}}
    weight = record_weight(records, kernel, {z: F(1, 2) for z in range(3)})
    assert weight > 0
    return dict(original_supports=[sorted(lower), sorted(upper)],
                repaired_upper=sorted(sigma), missing_background=sorted(B-sigma),
                target_indices=[n, m], first_auxiliary_arrival=auxiliary[2][0],
                positive_cylinder_weight=str(weight), forced_target_arrow=True,
                scope='Obstruction to background-preserving target avoidance only')


def target_history_component(records, lower_roots, upper_roots, target=0):
    """Union of emitted finite-vertex edges; sink exits do not join components."""
    neighbors = {z: set() for _, z, _ in records}
    for copy, roots in enumerate((lower_roots, upper_roots)):
        book = {(z, i): record for (c, z, i), record in records.items() if c == copy}
        trace = stabilize(book, roots)[3]
        for site, _, _, outgoing in trace:
            for destination in outgoing:
                if destination != -1:
                    neighbors[site].add(destination)
                    neighbors[destination].add(site)
    component, queue = {target}, [target]
    while queue:
        site = queue.pop()
        for destination in neighbors[site]-component:
            component.add(destination)
            queue.append(destination)
    return component


def switch_history_component(records, lower_roots, upper_roots, target=0):
    component = target_history_component(records, lower_roots, upper_roots, target)
    extra = Counter(upper_roots)-Counter(lower_roots)
    assert sum(extra.values()) == 1 and not (set(extra) & component)
    changed = {(1-copy if z in component else copy, z, i): record
               for (copy, z, i), record in records.items()}
    return changed, component


def check_history_components():
    rng = Random(20261010)
    p = [F(2, 5), F(3, 7), F(5, 11), F(7, 13)]
    counts = Counter()
    for sample in range(1024):
        lower_roots = [0, 2] if sample % 2 else [2]
        upper_roots = lower_roots+[1]
        records = {(copy, site, index):
                   ((int(rng.randrange(p[site].denominator) < p[site].numerator),
                     rng.choice([0, 1, 2, 3, -1])) if index <= 64 else (0, -1))
                   for copy in (0, 1) for site in range(4) for index in range(1, 73)}
        component = target_history_component(records, lower_roots, upper_roots)
        if 1 in component:
            continue
        counts['outside_source_component'] += 1
        changed, copied = switch_history_component(records, lower_roots, upper_roots)
        restored, decoded = switch_history_component(changed, lower_roots, upper_roots)
        assert restored == records and component == copied == decoded
        before = [creation_prefix(records, c, roots)[0]
                  for c, roots in enumerate((lower_roots, upper_roots))]
        after = [creation_prefix(changed, c, roots)[0]
                 for c, roots in enumerate((lower_roots, upper_roots))]
        for c in (0, 1):
            assert after[c] & component == before[1-c] & component
            assert after[c]-component == before[c]-component
        counts['bad_pairs_repaired'] += 0 in before[0] and 0 not in before[1]
        counts['bad_upper_never_visited_target'] += (
            0 in before[0] and not creation_prefix(records, 1, upper_roots)[1].get((1, 0), 0))
    assert counts['bad_upper_never_visited_target']
    return dict(seed=20261010, sampled_cylinders=1024, **counts,
                all_full_inverses=True, component_invariant=True,
                all_inside_supports_exchanged=True, all_outside_supports_unchanged=True)


def swap_colours(bits, i, j):
    changed = list(bits)
    changed[i], changed[j] = changed[j], changed[i]
    return tuple(changed)


def check_colours():
    result = []
    for lower, upper, expected in [(1, 2, F(1, 3)), (2, 3, F(4, 21))]:
        hits = [F(1, 2)] * (lower + upper)
        separate = colour_law(hits, split=lower)
        escaped = sum((weight for bits, weight in separate.items()
                       if swap_colours(bits, 0, lower) not in separate), F(0))
        assert escaped == expected
        joint = colour_law(hits)
        assert all(weight == joint[swap_colours(bits, 0, lower)]
                   for bits, weight in joint.items())
        result.append(dict(lower=lower, upper=upper,
                           separate_conditioning_escape=str(escaped),
                           joint_conditioning_invariant=True))
    # Unequal positions away from contact; only the exchanged h-values agree.
    hits = [F(2, 5), F(1, 7), F(2, 5), F(3, 4), F(1, 3)]
    joint = colour_law(hits)
    assert all(weight == joint[swap_colours(bits, 0, 2)]
               for bits, weight in joint.items())
    return result


if __name__ == "__main__":
    print(json.dumps({
        "record_switches": [check_record_witness(kind, active)
                            for kind in ("BS", "R") for active in (False, True)],
        "long_excursions": [check_long_excursion(length) for length in (1, 2, 6)],
        "common_wakeup": check_common_wakeup(),
        "variable_stopped_walks": check_variable_walks(),
        "terminal_swap_overlap": check_terminal_overlap(),
        "creator_contact_switch": check_creator_contact_switch(),
        "fomin_selector": check_fomin_selector(),
        "fomin_selector_two_sleeps": check_fomin_selector(True),
        "cycle_selector": check_cycle_selector(),
        "contact_then_fomin": check_contact_then_fomin(),
        "contact_then_fomin_two_sleeps": check_contact_then_fomin(True),
        "two_sleep_obstruction": check_two_sleep_obstruction(),
        "mutual_sleep_dependencies": check_mutual_sleep_dependencies(),
        "contact_then_cycle": check_contact_then_cycle(),
        "last_contact_creator": check_last_contact_creator(),
        "marked_returning_creator": check_marked_returning_creator(),
        "matched_avalanches": check_matched_avalanches(),
        "repaired_terminal_overlap": check_repaired_terminal_overlap(),
        "completed_word_contacts": check_completed_word_contacts(),
        "general_terminal_repair_samples": check_terminal_repair_samples(),
        "reserved_terminal_repair_samples": check_reserved_terminal_samples(),
        "forced_auxiliary_return": check_forced_auxiliary_return(),
        "history_components": check_history_components(),
        "colours": check_colours(),
        "status": "Exact local checks only; global selection/inverse and monotonicity open",
    }, indent=2))
