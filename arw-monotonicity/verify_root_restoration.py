"""Finite-book audits of the root-restoration theorem in response-renewal.md.

Written by Codex during the ongoing ARW research task. These are exact
combinatorial/rational checks, not a proof of occupation monotonicity.
Standard library only; no files are written when this script runs.
"""

from collections import Counter
from fractions import Fraction as F
from itertools import product
from random import Random
from verify_stack_exploration import stabilize


def counts(result):
    return Counter(z for z, _ in result[2])


def restore(records, initial, roots, largest=False):
    """Restore missing initial roots; discard edits outside the final prefix."""
    assert set(roots) <= set(initial)
    changed = dict(records)
    previous = None
    edits = []
    while True:
        result = stabilize(changed, initial)
        emitted = counts(result)
        if previous is not None:
            assert all(emitted[z] <= previous[z] for z in set(emitted) | set(previous))
            assert sum(emitted.values()) < sum(previous.values())
        missing = sorted(set(roots) - result[0])
        if not missing:
            canonical = dict(records)
            for z in roots:
                index = max(i for zz, i in result[1] if zz == z)
                canonical[z, index] = 1, canonical[z, index][1]
            audited = stabilize(canonical, initial)
            assert audited[:3] == result[:3]
            return canonical, audited, edits
        z = missing[-1 if largest else 0]
        index = max(i for zz, i in result[1] if zz == z)
        assert changed[z, index][0] == 0
        changed[z, index] = 1, changed[z, index][1]
        previous = emitted
        edits.append((z, index))


def terminal_subsets(records, initial, roots):
    """All inverse candidates and the downward-closed admissibility family."""
    root_list = sorted(roots)
    result = stabilize(records, initial)
    assert set(roots) <= result[0]
    terminal = {z: max(i for zz, i in result[1] if zz == z) for z in roots}
    good = []
    for mask in range(1 << len(root_list)):
        inverse = dict(records)
        for j, z in enumerate(root_list):
            if mask >> j & 1:
                key = z, terminal[z]
                inverse[key] = 0, inverse[key][1]
        canonical, _, _ = restore(inverse, initial, roots)
        if canonical == records:
            good.append(mask)
    for mask in good:
        for submask in range(1 << len(root_list)):
            if submask & mask == submask:
                assert submask in good
    return good


def finite_product_audit():
    """Exhaustive two-site product law with deterministic sleeping sink tails."""
    keys = [(z, i) for z in range(2) for i in (1, 2)]
    options = list(product((0, 1), (-1, 0, 1)))
    p = (F(2, 5), F(3, 7))
    ell = [(1 - t) / t for t in p]
    weights, pushforward, books = {}, {}, {}
    comparisons = 0
    for row in product(options, repeat=len(keys)):
        records = {(z, i): (1, -1) for z in range(2) for i in range(1, 12)}
        records.update(zip(keys, row))
        weight = F(1)
        for (z, _), (bit, _) in zip(keys, row):
            weight *= (p[z] if bit else 1 - p[z]) / 3
        weights[row] = weight
        books[row] = records
        for roots in ({0}, {1}, {0, 1}):
            canonical, result, _ = restore(records, [0, 1], roots)
            other, other_result, _ = restore(records, [0, 1], roots, largest=True)
            assert canonical == other and result[:3] == other_result[:3]
            lower = counts(result)
            for x in range(2):
                _, upper, _ = restore(records, [0, 1, x], roots)
                assert all(counts(upper)[z] >= lower[z] for z in range(2))
                comparisons += 1
            if roots == {0, 1}:
                image = tuple(canonical[key] for key in keys)
                pushforward[image] = pushforward.get(image, F(0)) + weight
        # Decreasing a sleep bit can only increase restored emission counts.
        for key in keys:
            if records[key][0]:
                less_sleep = dict(records)
                less_sleep[key] = 0, records[key][1]
                _, base, _ = restore(records, [0, 1], {0, 1})
                _, larger, _ = restore(less_sleep, [0, 1], {0, 1})
                assert all(counts(larger)[z] >= counts(base)[z] for z in range(2))
                comparisons += 1
    for image, mass in pushforward.items():
        good = terminal_subsets(books[image], [0, 1], {0, 1})
        final = stabilize(books[image], [0, 1])
        terminal_odds = [ell[z] if max(i for zz, i in final[1] if zz == z) <= 2
                         else F(0) for z in range(2)]
        density = sum((terminal_odds[0] if mask & 1 else 1) *
                      (terminal_odds[1] if mask & 2 else 1) for mask in good)
        assert mass == weights[image] * density
    assert sum(pushforward.values(), F(0)) == 1
    return dict(books=len(books), image_books=len(pushforward),
                monotone_comparisons=comparisons, density_identity=True)


def nonfactorizing_fibre():
    records = {(z, i): (1, -1) for z in range(3) for i in range(1, 12)}
    records.update({(0, 1): (1, 1), (2, 1): (1, 1),
                    (1, 1): (1, 0), (1, 2): (0, 2)})
    good = terminal_subsets(records, [0, 2], {0, 2})
    assert good == [0, 1, 2]  # Each singleton retracts; the pair does not.
    return dict(roots=[0, 2], admissible_masks=good)


def joint_fibre_obstruction():
    """A two-common-root fibre has signed target weight (1-ell_1)(1+ell_2)."""
    records = {(z, i): (1, -1) for z in range(3) for i in range(1, 12)}
    records.update({(0, 1): (1, -1), (0, 2): (0, 0),
                    (0, 3): (0, 1), (1, 1): (1, 0),
                    (1, 2): (0, -1)})
    lower, upper, roots = [1, 2], [0, 1, 2], {1, 2}
    def joint(book):
        result, lo, _ = restore(book, lower, roots)
        result, hi, _ = restore(result, upper, roots)
        assert stabilize(result, lower)[:3] == lo[:3]
        return result, lo, hi
    canonical, lo, hi = joint(records)
    assert canonical == records
    keys = sorted({(z, max(i for zz, i in rr[1] if zz == z))
                   for z in roots for rr in (lo, hi)})
    assert keys == [(1, 1), (2, 1)]
    coefficients = {}
    for mask in range(4):
        inverse = dict(records)
        for j, key in enumerate(keys):
            if mask >> j & 1:
                inverse[key] = 0, inverse[key][1]
        assert joint(inverse)[0] == records
        lo = stabilize(inverse, lower)[0]
        hi = stabilize(inverse, upper)[0]
        coefficients[mask] = int(0 in hi) - int(0 in lo)
    assert coefficients == {0: 1, 1: -1, 2: 1, 3: -1}
    return dict(joint_fibre_coefficients=coefficients,
                signed_weight_at_odds_3_1=-4)


def joint_fibre_all_parameters():
    """A fixed joint fibre has weight -ell_1 for every positive jump odds."""
    records = {(z, i): (1, -1) for z in range(2) for i in range(1, 12)}
    records.update({(0, 1): (0, 1), (0, 2): (0, -1),
                    (1, 1): (0, 1), (1, 2): (1, 0), (1, 3): (1, -1)})
    lower, upper, roots = [1], [0, 1], {1}

    def joint(book):
        result, lo, _ = restore(book, lower, roots)
        result, hi, _ = restore(result, upper, roots)
        assert stabilize(result, lower)[:3] == lo[:3]
        return result, lo, hi

    canonical, lo, hi = joint(records)
    assert canonical == records and lo[0] == hi[0] == {1}
    keys = [(1, 2), (1, 3)]
    coefficients = {}
    for mask in range(4):
        inverse = dict(records)
        for j, key in enumerate(keys):
            if mask >> j & 1:
                inverse[key] = 0, inverse[key][1]
        if joint(inverse)[0] == records:
            lo = stabilize(inverse, lower)[0]
            hi = stabilize(inverse, upper)[0]
            coefficients[mask] = int(1 in hi)-int(1 in lo)
    assert coefficients == {0: 0, 2: -1}
    return dict(all_parameters_fibre=coefficients, signed_weight="-ell_1")


def seeded_audit():
    rng = Random(20261015)
    checked = 0
    for n in (2, 3, 4):
        for _ in range(200):
            records = {(z, i): (rng.randrange(2), rng.randrange(-1, n))
                       for z in range(n) for i in range(1, 41)}
            records.update({(z, i): (1, -1) for z in range(n) for i in range(41, 70)})
            initial = rng.sample(range(n), rng.randint(1, n))
            roots = set(rng.sample(initial, rng.randint(1, len(initial))))
            lo, lr, _ = restore(records, initial, roots)
            hi, hr, _ = restore(records, initial, roots, largest=True)
            assert lo == hi and lr[:3] == hr[:3]
            _, upper, _ = restore(records, initial + [rng.randrange(n)], roots)
            assert all(counts(upper)[z] >= counts(lr)[z] for z in range(n))
            checked += 1
    return dict(seeded_books=checked)


def shifted_book(records, result):
    """Remove emitted prefixes; occupied sites retain their unread terminal arrow."""
    emitted = counts(result)
    return {(z, i-emitted[z]): record for (z, i), record in records.items()
            if i > emitted[z]}


def continuation_fibre_audit():
    rng = Random(2026101603)
    checked, inverse_checks = 0, 0
    for n in (2, 3, 4):
        for _ in range(200):
            records = {(z, i): (rng.randrange(2), rng.randrange(-1, n))
                       for z in range(n) for i in range(1, 45)}
            records.update({(z, i): (1, -1) for z in range(n) for i in range(45, 100)})
            initial = [rng.randrange(n) for _ in range(rng.randint(1, n+2))]
            root_list = rng.sample(sorted(set(initial)), rng.randint(1, len(set(initial))))
            roots = set(root_list)
            image, result, _ = restore(records, initial, roots)
            B = sorted(result[0])
            fresh = shifted_book(image, result)
            assert stabilize(fresh, B)[0] == set(B)
            good = terminal_subsets(image, initial, roots)
            fresh_good = terminal_subsets(fresh, B, roots)
            assert good == fresh_good
            for mask in range(1 << len(roots)):
                old_inverse, new_inverse = dict(image), dict(fresh)
                for j, z in enumerate(sorted(roots)):
                    if mask >> j & 1:
                        old_key = z, counts(result)[z]+1
                        old_inverse[old_key] = 0, old_inverse[old_key][1]
                        new_inverse[z, 1] = 0, new_inverse[z, 1][1]
                old_result = stabilize(old_inverse, initial)
                new_result = stabilize(new_inverse, B)
                assert old_result[0] == new_result[0]
                assert all(counts(old_result)[z]-counts(result)[z] == counts(new_result)[z]
                           for z in range(n))
                inverse_checks += 1
            checked += 1
    return dict(continuation_fibres=checked, inverse_replays=inverse_checks)


def fresh_boundary_book(base, roots, boundary):
    """Couple all fresh U_B from U_R: shift bits at B minus R, never arrows.

    Root first bits are already one and their later records stay fixed.
    At each occupied nonroot insert a leading one in its bit sequence only.
    Independent base bits/arrows give the RN16 fresh marginal for every B;
    these are not tails sampled under the biased canonical-image law.
    """
    assert roots <= boundary
    assert all(base[z, 1][0] == 1 for z in roots)
    return {(z, i): ((1 if i == 1 else base[z, i-1][0])
                     if z in boundary-roots else bit, arrow)
            for (z, i), (bit, arrow) in base.items()}


def _fresh_boundary_comparisons(base, roots):
    """Compare every nested B in one simultaneous finite-book coupling."""
    sites = {z for z, _ in base}
    outside = sorted(sites-roots)
    boundaries = [frozenset(roots | {z for j, z in enumerate(outside)
                                   if mask >> j & 1})
                  for mask in range(1 << len(outside))]
    books = {B: fresh_boundary_book(base, roots, B) for B in boundaries}
    emissions, edits, families = {}, {}, {}
    stats = Counter(coupled_boundaries=len(boundaries))
    for B, book in books.items():
        unflipped = stabilize(book, sorted(B))
        assert unflipped[0] == set(B) and not counts(unflipped)
        assert all(book[key][1] == base[key][1] for key in base)
        assert all(book[z, i] == base[z, i] for z, i in base if z in roots)
        families[B] = []
        for mask in range(1 << len(roots)):
            inverse = dict(book)
            for j, z in enumerate(sorted(roots)):
                if mask >> j & 1:
                    inverse[z, 1] = 0, inverse[z, 1][1]
            canonical, result, changes = restore(inverse, sorted(B), roots)
            emissions[B, mask] = counts(result)
            edits[B, mask] = changes
            admissible = canonical == book
            assert admissible == (not emissions[B, mask])
            if admissible:
                families[B].append(mask)
        assert 0 in families[B]
    for B in boundaries:
        for C in boundaries:
            if not B < C:
                continue
            lower, upper = books[B], books[C]
            for (z, i), (bit, arrow) in lower.items():
                if z in C-B:
                    assert upper[z, i] == (1 if i == 1 else lower[z, i-1][0], arrow)
                else:
                    assert upper[z, i] == (bit, arrow)
            for mask in range(1 << len(roots)):
                lo, hi = emissions[B, mask], emissions[C, mask]
                assert all(hi[z] >= lo[z] for z in sites), (B, C, mask, lo, hi)
                # Every upper restoration edit lies at/after the lower
                # canonical terminal index. This is the proof's induction.
                assert all(i >= lo[z]+1 for z, i in edits[C, mask])
                stats['canonical_emission_comparisons'] += 1
                stats['restoration_edit_barriers'] += len(edits[C, mask])
            assert set(families[C]) <= set(families[B]), (B, C, families)
            stats['nested_boundary_pairs'] += 1
            stats['strict_fibre_inclusions'] += int(set(families[C]) < set(families[B]))
    return stats, families


def fresh_boundary_coupling_audit():
    """Bounded checks of decreasing fresh fibres; no parameter search.

    All boundaries containing {0,2} share a base book in each of 300
    seeded RN12 four-cycle trials. Deterministic sleeping/sink tails make
    these finite-book checks, not samples of the infinite ARW fibre law.
    Also retain the nonfactorizing and later-root-bit witnesses exactly.
    """
    seed, trials = 20261016091, 300
    rng = Random(seed)
    p, roots = (F(1, 2), F(1, 2), F(1, 8), F(1, 8)), {0, 2}
    total = Counter()
    for _ in range(trials):
        base = {(z, i): (int(rng.randrange(p[z].denominator) < p[z].numerator),
                         (z+1) % 4 if rng.randrange(32) < 31 else -1)
                for z in range(4) for i in range(1, 51)}
        base.update({(z, i): (1, -1) for z in range(4) for i in range(51, 151)})
        for z in roots:
            base[z, 1] = 1, base[z, 1][1]
        stats, _ = _fresh_boundary_comparisons(base, roots)
        total.update(stats)

    witnesses = []
    base = {(z, i): (1, -1) for z in range(3) for i in range(1, 60)}
    base.update({(0, 1): (1, 1), (2, 1): (1, 1),
                 (1, 1): (1, 0), (1, 2): (0, 2)})
    stats, families = _fresh_boundary_comparisons(base, {0, 2})
    pair = [families[frozenset({0, 2})], families[frozenset(range(3))]]
    assert pair == [[0, 1, 2], [0]]
    total.update(stats)
    witnesses.append(dict(name='RN11_nonfactorizing', lower_upper_masks=pair))
    for bit in (1, 0):
        base = {(z, i): (1, -1) for z in range(3) for i in range(1, 60)}
        base.update({(0, 1): (1, 1), (0, 2): (0, 2),
                     (1, 1): (1, 1), (1, 2): (0, 2), (1, 3): (bit, 0),
                     (2, 1): (1, 2), (2, 2): (0, 0), (2, 3): (0, 1)})
        stats, families = _fresh_boundary_comparisons(base, {0, 1})
        pair = [families[frozenset({0, 1})], families[frozenset(range(3))]]
        assert pair == [[0, 1] if bit else [0], [0]]
        total.update(stats)
        witnesses.append(dict(name='RN16_later_root_bit', bit=bit, lower_upper_masks=pair))
    return dict(seed=seed, fixed_kernel='RN12_four_cycle', seeded_base_books=trials,
                **total, witnesses=witnesses)


def fresh_boundary_non_supermodularity_witness():
    """A fixed-book mixed difference is -r_0; no averaged claim is made."""
    base = {(z, i): (1, -1) for z in range(3) for i in range(1, 12)}
    base.update({(0, 1): (1, 1), (1, 1): (1, 2), (2, 1): (1, 0)})
    _, families = _fresh_boundary_comparisons(base, {0})
    boundaries = list(map(frozenset, ({0}, {0, 1}, {0, 2}, {0, 1, 2})))
    masks = [families[B] for B in boundaries]
    assert masks == [[0, 1], [0, 1], [0, 1], [0]]
    # Coefficients of 1 and r_0 in Z(B); retain the exact formal sign.
    coefficients = [(int(0 in family), int(1 in family)) for family in masks]
    mixed = [coefficients[3][j]-coefficients[1][j]-coefficients[2][j]
             +coefficients[0][j] for j in range(2)]
    assert mixed == [0, -1]
    return dict(non_supermodular_boundaries=list(map(sorted, boundaries)),
                admissible_masks=masks, mixed_difference='-r_0',
                scope='fixed-book fibre weight only; no claim about averaged m_R')


def archived_entrance_mass_audit():
    """Exact m_R monotonicity and normalized rooted bounds on saved matrices.

    Test all nonempty R, all boundary enlargements and all binary inputs
    D containing R. RN16's input-independent fresh G=Gamma gives the same
    normalization for arbitrary all-active eta>=1_R; that general identity
    is not inferred from these finite binary checks.
    """
    import verify_response_renewal as rn

    audits = []
    for case_number, case in enumerate(rn.CASES):
        supports = [frozenset(z for z, bit in enumerate(bits) if bit)
                    for bits in case['subsets']]
        pos = {B: i for i, B in enumerate(supports)}
        P = [[F(q) for q in row] for row in case['stabilization_matrix']]
        assert all(sum(row) == 1 and min(row) >= 0 for row in P)
        stats = Counter()
        for roots in supports:
            if not roots:
                continue
            boundary = [B for B in supports if roots <= B]
            Q = [[P[pos[D]][pos[B]] for B in boundary] for D in boundary]
            inverse = rn.inverse(Q)
            gamma = rn.multiply(inverse, [P[pos[D]] for D in boundary])
            assert all(min(row) >= 0 for row in gamma)
            masses = dict(zip(boundary, map(sum, gamma)))
            assert list(masses.values()) == list(map(sum, inverse))  # Q^-1 1.
            for B, row in zip(boundary, gamma):
                assert all(row[pos[C]] == F(B == C) for C in boundary)
                assert all(not row[pos[C]] for C in supports if len(C) > len(B))
                assert sum(row[pos[C]] for C in boundary)/masses[B] == 1/masses[B]
            baseline = P[pos[roots]][pos[roots]]
            assert masses[roots] == 1/baseline
            for B in boundary:
                for C in boundary:
                    if B < C:
                        assert masses[C] <= masses[B], (case_number, roots, B, C)
                        stats['nested_mass_inequalities'] += 1
                        stats['strict_mass_inequalities'] += int(masses[C] < masses[B])
            for D in boundary:
                weights = [P[pos[D]][pos[B]] for B in boundary]
                # Check the whole entrance expansion, including its exact
                # normalization, across every original binary input D.
                assert rn.multiply([weights], gamma)[0] == P[pos[D]]
                assert sum(w*masses[B] for B, w in zip(boundary, weights)) == 1
                joint = sum(weights)
                slack = sum(w*(masses[roots]-masses[B])
                            for B, w in zip(boundary, weights))
                assert slack == masses[roots]*joint-1 >= 0
                assert joint >= baseline
                stats['normalized_input_expansions'] += 1
            stats['root_sets'] += 1
        audits.append(dict(case=case_number, n=len(case['K']), **stats))
    return dict(archived_entrance_mass_audits=audits)


def scc_entrance_factorization_audit():
    """One fixed graph: SCC product, self-loops, and strict source reachability.

    Its condensation is {3}->{0,1}->{2}; roots are {0,2}, so U={1}.
    Vertex 3 is a root-free upstream component. Sources 0 and 2 test
    multiplicities at roots, with strict gain only for the former.
    """
    import verify_response_renewal as rn
    from verify_feedback_cone import binary_stabilization

    K = [[F(1,4), F(1,2), F(0), F(0), F(1,4)],
         [F(1,2), F(0), F(1,4), F(0), F(1,4)],
         [F(0), F(0), F(1,3), F(0), F(2,3)],
         [F(1,2), F(0), F(0), F(1,4), F(1,4)]]
    p = [F(1,2), F(2,3), F(1,3), F(3,5)]
    roots = frozenset({0, 2})
    components = list(map(frozenset, ({3}, {0, 1}, {2})))

    def masses(supports, P, R):
        pos = {B: i for i, B in enumerate(supports)}
        boundary = [B for B in supports if R <= B]
        Q = [[P[pos[D]][pos[B]] for B in boundary] for D in boundary]
        return dict(zip(boundary, map(sum, rn.inverse(Q))))

    supports, P, _, _ = binary_stabilization(K, p)
    global_m = masses(supports, P, roots)
    local_m = []
    for D in components:
        vertices = sorted(D)
        killed_K = [[K[z][w] for w in vertices] +
                    [sum(K[z][w] for w in range(5) if w not in D)]
                    for z in vertices]
        states, local_P, _, _ = binary_stabilization(killed_K, [p[z] for z in vertices])
        states = [frozenset(vertices[z] for z in C) for C in states]
        local_m.append(masses(states, local_P, roots & D))
    for B, value in global_m.items():
        expected = F(1)
        for D, factor in zip(components, local_m):
            expected *= factor[B & D]
        assert value == expected
        for y in set(range(4))-B:
            difference = value-global_m[B | {y}]
            assert difference >= 0 and (difference > 0) == (y == 1)
    assert global_m[roots] == global_m[roots | {3}] == F(35,9)
    assert global_m[roots | {1}] == global_m[roots | {1,3}] == F(245,72)

    def family(book, B, R):
        return {frozenset(z for j, z in enumerate(sorted(R)) if mask >> j & 1)
                for mask in terminal_subsets(book, sorted(B), R)}

    rng = Random(2026101618)
    routes = ((0, 1, 1, -1), (0, 0, 2, -1), (2, -1, -1), (0, 0, 3, -1))
    checked = 0
    for _ in range(30):
        base = {(z, i): (int(rng.randrange(p[z].denominator) < p[z].numerator),
                         rng.choice(routes[z]))
                for z in range(4) for i in range(1, 31)}
        base.update({(z, i): (1, -1) for z in range(4) for i in range(31, 100)})
        for z in roots:
            base[z, 1] = 1, base[z, 1][1]
        for B in global_m:
            book = fresh_boundary_book(base, roots, B)
            factors = []
            for D in components:
                killed = {(z, i): (bit, arrow if arrow in D else -1)
                          for (z, i), (bit, arrow) in book.items() if z in D}
                factors.append(family(killed, B & D, roots & D))
            expected = {frozenset().union(*parts) for parts in product(*factors)}
            assert family(book, B, roots) == expected
            checked += 1

    starts = {x: tuple(int(z in roots)+int(z == x) for z in range(4)) for x in range(4)}
    raw, states = rn.stopped(starts.values(), lambda s: rn.arw_step(K, p, s),
                             lambda s: all(q <= 0 for q in s))
    gains = [sum(q for end, q in raw[starts[x]].items() if all(end[z] == -1 for z in roots))
             -1/global_m[roots] for x in range(4)]
    assert gains == [F(3,175), F(9,350), F(0), F(2,525)]
    assert [x for x, gain in enumerate(gains) if gain > 0] == [0, 1, 3]
    return dict(sccs=list(map(sorted, components)), roots=sorted(roots),
                factorized_finite_book_families=checked, exact_boundary_products=len(global_m),
                boundary_masses=[(sorted(B), str(q)) for B, q in global_m.items()],
                joint_gains_by_source=list(map(str, gains)), raw_transient_states=states)


def pop_arrow_cycles(records, n):
    """Ordinary directed cycle popping; returns prefix counts and popped range."""
    emitted = [0]*n
    while True:
        arrows = {z: records[z, emitted[z]+1][1] for z in range(n)}
        cycle = None
        for start in range(n):
            path, location = [], {}
            z = start
            while z != -1 and z not in location:
                location[z] = len(path)
                path.append(z)
                z = arrows[z]
            if z != -1:
                cycle = path[location[z]:]
                break
        if cycle is None:
            return emitted
        for z in cycle:
            emitted[z] += 1
        assert sum(emitted) < 1000, "Unexpected cycle-popping loop in a finite book"


def full_root_cycle_audit():
    rng = Random(2026101602)
    checked = 0
    for n in (2, 3, 4):
        for _ in range(300):
            records = {(z, i): (rng.randrange(2), rng.randrange(-1, n))
                       for z in range(n) for i in range(1, 30)}
            records.update({(z, i): (1, -1) for z in range(n) for i in range(30, 100)})
            initial = list(range(n))
            image, result, _ = restore(records, initial, set(initial))
            fresh = shifted_book(image, result)
            emitted = pop_arrow_cycles(fresh, n)
            allowed = sum(1 << z for z in range(n) if emitted[z] == 0)
            expected = [mask for mask in range(1 << n) if mask & allowed == mask]
            assert terminal_subsets(image, initial, set(initial)) == expected
            checked += 1
    return dict(full_root_cycle_fibres=checked)


def partial_root_bit_dependence():
    records = {(z, i): (1, -1) for z in range(3) for i in range(1, 20)}
    records.update({(0, 1): (1, 1), (0, 2): (0, 2),
                    (1, 1): (1, 1), (1, 2): (0, 2), (1, 3): (1, 0),
                    (2, 1): (1, 2), (2, 2): (0, 0), (2, 3): (0, 1)})
    before = terminal_subsets(records, [0, 1], {0, 1})
    records[1, 3] = 0, 0
    after = terminal_subsets(records, [0, 1], {0, 1})
    assert before == [0, 1] and after == [0]
    return dict(partial_root_families=[before, after], changed_bit=(1, 3))


def killed_root_replay_audit():
    """Killing at root-set exits permits an arrow-only creation representation."""
    rng = Random(2026101606)
    checked, inverse_checks = 0, 0
    for n in (3, 4, 5):
        for _ in range(200):
            R = set(rng.sample(range(n), rng.randint(1, min(3, n-1))))
            roots = sorted(R)
            book = {(z, i): (rng.randrange(2), rng.randrange(-1, n))
                    for z in range(n) for i in range(1, 40)}
            book.update({(z, i): (1, -1) for z in range(n) for i in range(40, 100)})
            killed = {key: (bit, arrow if key[0] in R and arrow in R else -1)
                      for key, (bit, arrow) in book.items()}
            killed_image, result, _ = restore(killed, roots, R)
            image = {key: (killed_image[key][0], arrow)
                     for key, (_, arrow) in book.items()}
            assert stabilize(image, roots)[:3] == result[:3]
            fresh = shifted_book(image, result)
            shifted_killed = shifted_book(killed_image, result)
            popped = pop_arrow_cycles(shifted_killed, n)
            allowed = sum(1 << j for j, z in enumerate(roots) if popped[z] == 0)
            good = terminal_subsets(killed_image, roots, R)
            assert good == [mask for mask in range(1 << len(roots)) if mask & allowed == mask]
            for mask in good:
                old, new = dict(image), dict(fresh)
                for j, z in enumerate(roots):
                    if mask >> j & 1:
                        key = z, counts(result)[z]+1
                        old[key] = 0, old[key][1]
                        new[z, 1] = 0, new[z, 1][1]
                old_end = stabilize(old, roots)
                new_end = stabilize(new, roots)
                assert old_end[0] == new_end[0]
                assert all(counts(old_end)[z]-counts(result)[z] == counts(new_end)[z]
                           for z in range(n))
                inverse_checks += 1
            checked += 1
    return dict(killed_root_replays=checked, full_graph_inverses=inverse_checks)



def active_cycle_pop(records, roots, largest=False):
    roots = set(roots)
    popped = Counter()
    while True:
        arrows = {z: records[z, popped[z]+1][1] for z in roots}
        cycles = set()
        for start in sorted(roots):
            path, at = [], {}
            z = start
            while z in roots and z not in at:
                at[z] = len(path)
                path.append(z)
                z = arrows[z]
            if z in at:
                cycle = path[at[z]:]
                pivot = cycle.index(min(cycle))
                cycle = tuple(cycle[pivot:]+cycle[:pivot])
                if any(records[y, popped[y]+1][0] == 0 for y in cycle):
                    cycles.add(cycle)
        if not cycles:
            shifted = {(z, i-popped[z]): record for (z, i), record in records.items()
                       if i > popped[z]}
            return shifted, popped
        cycle = sorted(cycles)[-1 if largest else 0]
        for z in cycle:
            popped[z] += 1
        assert sum(popped.values()) < 1000


def active_cycle_replay_audit():
    rng = Random(2026101608)
    checked = 0
    for n in (2, 3, 4, 5):
        for _ in range(200):
            roots = set(rng.sample(range(n), rng.randint(1, n)))
            book = {(z, i): (rng.randrange(2), rng.randrange(-1, n))
                    for z in range(n) for i in range(1, 35)}
            book.update({(z, i): (1, -1) for z in range(n) for i in range(35, 100)})
            shifted, popped = active_cycle_pop(book, roots)
            other, other_popped = active_cycle_pop(book, roots, largest=True)
            assert other == shifted and other_popped == popped
            for initial in (sorted(roots), sorted(roots)+[rng.randrange(n)]):
                old = stabilize(book, initial)
                new = stabilize(shifted, initial)
                assert old[0] == new[0]
                assert all(counts(old)[z] == popped[z]+counts(new)[z] for z in range(n))
                checked += 1
    return dict(active_cycle_replays=checked)

def active_cycle_lattice_obstruction():
    x = F(12, 25)
    h = {frozenset(C): F(0) for C in ((0,), (0, 1), (0, 2), (0, 1, 2))}
    K = [[F(1, 25) if y == 3 else (x if y != z else F(0))
          for y in range(4)] for z in range(3)]
    for arrows in product(range(4), repeat=3):
        weight = F(1)
        for z, y in enumerate(arrows):
            weight *= K[z][y]
        cyclic = set()
        for start in range(3):
            path, at = [], {}
            z = start
            while z < 3 and z not in at:
                at[z] = len(path)
                path.append(z)
                z = arrows[z]
            if z in at:
                cyclic.update(path[at[z]:])
        for C in h:
            if cyclic <= C:
                h[C] += weight
    d = h[frozenset({0})]
    pair = h[frozenset({0, 1})]
    assert h[frozenset({0, 2})] == pair and h[frozenset(range(3))] == 1
    gap = pair*pair-d
    assert d == F(1369, 15625) and pair == F(4969, 15625)
    assert gap == F(3300336, 244140625)
    return dict(active_cycle_lattice_gap=str(gap))



def single_hole_cycle_audit():
    """Exact boundary states after erasing a cycle through a missing input."""
    rng = Random(2026101609)
    checked = 0
    for k in range(2, 6):
        n = k+2
        for bits in product((0, 1), repeat=k):
            if all(bits):
                continue
            for trial in range(20):
                U = {(z, i): (rng.randrange(2), rng.randrange(-1, n))
                     for z in range(n) for i in range(1, 40)}
                U.update({(z, i): (0, -1) for z in range(n) for i in range(40, 110)})
                old = {(z, i+1 if z < k else i): record for (z, i), record in U.items()}
                old.update({(z, 1): (bits[z], (z+1) % k) for z in range(k)})
                extras = [rng.randrange(n) for _ in range(trial % 4)]
                lower = stabilize(old, list(range(1, k))+extras)[0]
                boundary, residents = U, {}
                if bits[0] == 0 and 0 in bits[1:]:
                    pending = list(range(1, k))+extras
                elif bits[0] == 1:
                    j = next(i for i in range(1, k) if bits[i] == 0)
                    pending = list(range(j+1, k))+extras
                    residents = {i: ((i, 0), i+1) for i in range(j)}
                else:
                    pending = extras
                    residents = {i: ((i, 0), (i+1) % k) for i in range(1, k)}
                    boundary = {(z, i+1 if z == 0 else i): record
                                for (z, i), record in U.items()}
                    boundary[0, 1] = 0, 1
                assert lower == stabilize(boundary, pending, residents)[0]
                assert stabilize(old, list(range(k))+extras)[0] == stabilize(
                    U, list(range(k))+extras)[0]
                checked += 1
    # An unused zero at the hole cannot be discarded in a dormant component.
    U = {(z, i): (0, -1) for z in range(2) for i in range(1, 20)}
    U[0, 1] = 1, -1
    kept = {(z, i+1 if z == 0 else i): record for (z, i), record in U.items()}
    kept[0, 1] = 0, 1
    resident = {1: ((1, 0), 0)}
    assert stabilize(kept, [], resident)[0] == {1}
    assert stabilize(kept, [0], resident)[0] == {0}
    assert stabilize(U, [0], resident)[0] == {0, 1}
    # With two reserved sleepers, one-sleeper boundaries cannot reproduce the
    # no-arrival result on a zero/sink suffix.
    book = {(z, i): (0, -1) for z in range(3) for i in range(1, 20)}
    book.update({(0, 1): (0, 1), (1, 1): (1, 2), (2, 1): (1, 0)})
    assert stabilize(book, [1, 2])[0] == {1, 2}
    self_loop_checks = 0
    for trial in range(200):
        n = 4
        U = {(z, i): (rng.randrange(2), rng.randrange(-1, n))
             for z in range(n) for i in range(1, 35)}
        U.update({(z, i): (0, -1) for z in range(n) for i in range(35, 100)})
        old = {(z, i+1 if z == 0 else i): record for (z, i), record in U.items()}
        old[0, 1] = 0, 0
        extras = [rng.randrange(n) for _ in range(trial % 5)]
        assert stabilize(old, extras)[0] == stabilize(U, extras)[0]
        self_loop_checks += 1
    return dict(single_hole_joint_checks=checked, missing_input_self_loops=self_loop_checks)



def forest_promotion_audit():
    """Staged popping and the reverse-local-word closed-walk bijection."""
    rng = Random(2026101610)
    checked = 0
    for n in (2, 3, 4, 5):
        for _ in range(100):
            T = set(rng.sample(range(n), rng.randint(1, n)))
            x = rng.choice(sorted(T))
            R = T-{x}
            records = {(z, i): (0, rng.randrange(-1, n))
                       for z in range(n) for i in range(1, 35)}
            records.update({(z, i): (0, -1)
                            for z in range(n) for i in range(35, 100)})
            after_R, old = active_cycle_pop(records, R)
            after_T, extra = active_cycle_pop(after_R, T)
            direct, total = active_cycle_pop(records, T)
            assert direct == after_T
            assert all(total[z] == old[z]+extra[z] for z in range(n))
            beta = {z: [after_R[z, i][1] for i in range(1, extra[z]+1)] for z in T}
            support = {z for z in T if beta[z]}
            if support:
                assert x in support
                for start in support-{x}:
                    z, seen = start, set()
                    while z != x and z not in seen:
                        seen.add(z)
                        z = beta[z][0]
                    assert z == x
                reversed_words = {z: list(reversed(beta[z])) for z in T}
                positions = Counter()
                departures = {z: [] for z in T}
                z = x
                while positions[z] < len(reversed_words[z]):
                    y = reversed_words[z][positions[z]]
                    positions[z] += 1
                    departures[z].append(y)
                    z = y
                assert z == x and all(positions[z] == len(beta[z]) for z in T)
                assert all(list(reversed(departures[z])) == beta[z] for z in T)
            for z in R:
                restored = beta[z][0] if beta[z] else after_T[z, 1][1]
                assert restored == after_R[z, 1][1]
            checked += 1
    return dict(staged_forests_and_closed_words=checked)



def record_prefix_projection_audit():
    """Finite paired-word checks of the greatest M_R left divisor."""
    from itertools import combinations
    n = 3
    domains = [frozenset(C) for k in range(n+1)
               for C in combinations(range(n), k)]

    def project(h, R):
        choices = [[k for k in range(len(row)+1)
                    if (z in R or k == 0)
                    and (k == len(row) or row[k][0] == 1)]
                   for z, row in enumerate(h)]
        feasible = []
        for cut in product(*choices):
            incoming = [0]*n
            for row, k in zip(h, cut):
                for _, y in row[:k]:
                    incoming[y] += 1
            if tuple(incoming) == cut:
                feasible.append(cut)
        greatest = tuple(max(c[z] for c in feasible) for z in range(n))
        assert greatest in feasible
        return (tuple(row[:k] for row, k in zip(h, greatest)),
                tuple(row[k:] for row, k in zip(h, greatest)))

    counts = checks = 0
    G = {R: [F(0)]*5 for R in domains}
    quotient = {(R, x): [F(0)]*5 for R in domains for x in range(n) if x not in R}
    p = [F(1, 3), F(2, 5), F(3, 7)]
    for degree in range(5):
        lengths = [a for a in product(range(degree+1), repeat=n) if sum(a) == degree]
        for sizes in lengths:
            for arrows in product(range(n), repeat=degree):
                incoming = [arrows.count(z) for z in range(n)]
                if tuple(incoming) != sizes:
                    continue
                rows, offset = [], 0
                for k in sizes:
                    rows.append(arrows[offset:offset+k])
                    offset += k
                free = degree-sum(k > 0 for k in sizes)
                for later in product((0, 1), repeat=free):
                    h, offset, weight = [], 0, F(1)
                    for z, row in enumerate(rows):
                        bits = ((1,)+later[offset:offset+len(row)-1]) if row else ()
                        offset += max(0, len(row)-1)
                        h.append(tuple(zip(bits, row)))
                        for bit in bits:
                            weight *= (p[z] if bit else 1-p[z])/5
                    h = tuple(h)
                    L = frozenset(z for z, row in enumerate(h) if row)
                    projections = {}
                    for R in domains:
                        alpha, beta = project(h, R)
                        projections[R] = alpha
                        assert all(a+b == row for a, b, row in zip(alpha, beta, h))
                        assert not any(project(beta, R)[0])
                        if L <= R:
                            G[R][degree] += weight
                        for x in range(n):
                            if x not in R and L <= R | {x} and not any(alpha):
                                quotient[R, x][degree] += weight
                    for R in domains:
                        for S in domains:
                            if R <= S:
                                assert project(projections[S], R)[0] == projections[R]
                                checks += 1
                    counts += 1
    for (R, x), series in quotient.items():
        for degree in range(5):
            assert sum(G[R][k]*series[degree-k] for k in range(degree+1)) == G[R | {x}][degree]
    return dict(paired_words=counts, nested_prefix_checks=checks,
                promotion_series_degree=4)


def surplus_forest_path_audit():
    """Consumed balanced excursions, forest path, and remaining U-prefix."""
    rng = Random(2026101614)
    checked = 0
    for n in (2, 3, 4, 5):
        for _ in range(100):
            R = set(rng.sample(range(n), rng.randint(1, n)))
            source = rng.choice(sorted(R))
            book = {(z, i): (0, rng.randrange(-1, n))
                    for z in range(n) for i in range(1, 35)}
            book.update({(z, i): (0, -1) for z in range(n) for i in range(35, 150)})
            after, popped = active_cycle_pop(book, R)
            forest_path, z = [], source
            while z in R:
                assert z not in forest_path
                forest_path.append(z)
                z = after[z, 1][1]
            exit_label = z
            P, used, visited, le = set(forest_path), Counter(), set(), [source]
            z = source
            while z in R:
                visited.add(z)
                used[z] += 1
                y = book[z, used[z]][1]
                if y in R:
                    if y in le:
                        le = le[:le.index(y)+1]
                    else:
                        le.append(y)
                z = y
            assert z == exit_label and le == forest_path
            eps = {z: used[z]-int(z in P) for z in R}
            delta = {z: popped[z]-eps[z] for z in R}
            assert all(eps[z] >= 0 and delta[z] >= 0 for z in R)
            assert all(delta[z] == 0 for z in P)
            assert all(used[z] <= popped[z]+int(z in P) for z in R)
            for cuts, start in ((eps, {z: 0 for z in R}), (delta, eps)):
                incoming = Counter()
                for z in R:
                    for i in range(start[z]+1, start[z]+cuts[z]+1):
                        incoming[book[z, i][1]] += 1
                assert all(incoming[z] == cuts[z] for z in R)
                assert not any(incoming[z] for z in incoming if z not in R)
            L = {z for z in R if eps[z]}
            assert visited == P | L
            for a in L-P:
                z, seen = a, set()
                while z not in P and z not in seen:
                    seen.add(z)
                    z = book[z, eps[z]][1]
                assert z in P
            shifted = {(z, i): book[z, i+used[z]]
                       for z in range(n) for i in range(1, 100)}
            _, rest = active_cycle_pop(shifted, R-P)
            assert all(rest[z] == delta[z] for z in R)
            checked += 1
    return dict(surplus_forest_path_checks=checked)


if __name__ == "__main__":
    print(record_prefix_projection_audit(), flush=True)
    print(surplus_forest_path_audit(), flush=True)

    print(forest_promotion_audit(), flush=True)
    print(single_hole_cycle_audit(), flush=True)
    print(active_cycle_replay_audit(), flush=True)
    print(active_cycle_lattice_obstruction(), flush=True)
    print(finite_product_audit(), flush=True)
    print(nonfactorizing_fibre(), flush=True)
    print(joint_fibre_obstruction(), flush=True)
    print(joint_fibre_all_parameters(), flush=True)
    print(seeded_audit(), flush=True)
    print(continuation_fibre_audit(), flush=True)
    print(fresh_boundary_coupling_audit(), flush=True)
    print(fresh_boundary_non_supermodularity_witness(), flush=True)
    print(archived_entrance_mass_audit(), flush=True)
    print(scc_entrance_factorization_audit(), flush=True)
    print(full_root_cycle_audit(), flush=True)
    print(partial_root_bit_dependence(), flush=True)
    print(killed_root_replay_audit(), flush=True)
