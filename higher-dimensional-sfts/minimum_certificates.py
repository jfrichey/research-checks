"""Exact local-lemma certificates for strict rectangular entropy minimality.

Every comparison is rational. This certifies a finite set of inequalities on
an infinite conditional dependency graph by checking the finite perturbation
of a constant assignment. The research note gives the mathematical reduction
and the separate large-support proof. Standard library only.
"""

from fractions import Fraction as F
from itertools import permutations
import json
from math import factorial
from pathlib import Path

from ordering_checks import agreement, binary_representatives, displacements, rectangle


def certificate(h, w, dx, dy, q):
    s = h*w
    degree = (2*h-1)*(2*w-1)-1
    p = F(1, q**s)
    base = 2*p
    steps = [(a, b) for a in range(1-h, h) for b in range(1-w, w)
             if (a, b) != (0, 0)]

    def allowed(a, b):
        return abs(a) >= h or abs(b) >= w

    # C consists of anchors whose supports miss the first fixed copy.
    # J consists of the anchors in C which meet the second fixed copy.
    near = {(dx+a, dy+b) for a, b in steps if allowed(dx+a, dy+b)}
    checked = near | {(a+x, b+y) for a, b in near for x, y in steps
                      if allowed(a+x, b+y)}

    def probability(a, b):
        overlap = max(0, h-abs(a-dx))*max(0, w-abs(b-dy))
        return F(1, q**(s-overlap))

    def value(a, b):
        return 2*probability(a, b)

    assert 0 < base < 1
    assert all(0 < value(a, b) < 1 for a, b in near)
    # Far from J, degree is at most the full overlap degree and x=2p.
    minimum = (1-base)**degree
    worst = None
    for a, b in sorted(checked):
        product = F(1)
        for x, y in steps:
            if allowed(a+x, b+y):
                product *= 1-value(a+x, b+y)
        if product < minimum:
            minimum, worst = product, (a, b)
    # Since x_k=2 p'_k, the local-lemma inequality is exactly product >=1/2.
    assert 2*minimum >= 1, (h, w, dx, dy, q, minimum)
    avoid_factor = F(1)
    for a, b in near:
        avoid_factor *= 1-value(a, b)
    overlap = (h-abs(dx))*(w-abs(dy))
    pair_probability = F(1, q**(2*s-overlap))
    return dict(shape=[h, w], displacement=[dx, dy], q=q,
                multiplier=2, changed_probability_sites=len(near),
                explicitly_checked_sites=len(checked), worst_anchor=worst,
                minimum_neighbor_product=str(minimum),
                displayed_decimal=float(minimum),
                conditional_avoidance_factor=str(avoid_factor),
                pair_probability=str(pair_probability),
                defect_lower_bound=str(pair_probability*avoid_factor),
                exact_test_passed=True)


def binary_square_minimum():
    """Separate certificate for the exceptional binary 3x3 minimum.

    The lower comparison is the positive formal Shearer partition function,
    not an overlap-free binary pattern (none exists). All quantities below
    are rational bounds on the analytic estimates proved in the note.
    """
    p, rho, a = F(1, 512), F(1, 256), F(1, 16)
    exp_a_upper = F(16, 15)  # exp(a) <= 1/(1-a)
    kp = 25*p*exp_a_upper/(1-rho*exp_a_upper)**2
    assert kp <= a
    upper = p*rho*exp_a_upper**2/(1-rho**2*exp_a_upper**2)
    assert upper == F(1, 115198)
    lower = {}
    x = F(1, 100)
    assert x*(1-x)**24 >= F(1, 128)
    lower['maximal_overlap_2'] = F(1, 2**16)*(1-x)**24
    x = F(1, 200)
    assert x*(1-x)**24 >= F(1, 256)
    lower['two_corner_agreement_pairs'] = F(2, 2**17)*(1-x)**25
    local = [certificate(3, 3, dx, dy, 2)
             for dx, dy in [(0, 2), (1, 1), (0, 1)]]
    for r, record in zip((3, 4, 6), local):
        lower[f'maximal_overlap_{r}'] = F(record['defect_lower_bound'])
    gaps = {name: c/(1+c)-upper for name, c in lower.items()}
    assert min(gaps.values()) > 0
    return dict(status='Exact rational constants for the working proof of '
                       'strict binary 3x3 minimality at every negative '
                       'coupling and at the hard endpoint.',
                kp_parameter=str(a), kp_sum_bound=str(kp),
                single_corner_upper_bound=str(upper),
                defect_bounds={k: str(v) for k, v in lower.items()},
                pressure_gap_coefficients={k: str(v) for k, v in gaps.items()},
                minimum_displayed_gap=float(min(gaps.values())),
                additional_local_certificates=local)


def exp_upper(x, degree=10):
    """Rational Taylor upper bound, with a geometric bound on the tail."""
    assert 0 <= x < degree+2
    return (sum(x**j/F(factorial(j)) for j in range(degree+1))
            + x**(degree+1)/F(factorial(degree+1))/(1-x/F(degree+2)))


def conditioned_pattern_certificate(h, w, g, fixed, extra, multiplier, single):
    """Local lemma after fixing compatible copies of a specified G class.

    For one fixed copy, use all other anchors. For a pair, use lexicographically
    forward anchors disjoint from the first copy, plus specified extra events
    to avoid. Impossible conditional events receive probability zero.
    """
    support = set(rectangle(h, w))
    ds = displacements(support)
    steps = set(ds) | {(-a, -b) for a, b in ds}
    good = {d for i, d in enumerate(ds) if g >> i & 1}
    good |= {(-a, -b) for a, b in good}
    union = set().union(*({(a+x, b+y) for a, b in support} for x, y in fixed))
    p, degree = F(1, 2**(h*w)), len(steps)

    def allowed(z):
        if single:
            return z != (0, 0)
        return z > (0, 0) and (z not in steps or z in extra)

    near = {(a+x, b+y) for x, y in fixed for a, b in steps | {(0, 0)}
            if allowed((a+x, b+y))}
    checked = near | {(a+x, b+y) for a, b in near for x, y in steps
                      if allowed((a+x, b+y))}

    def value(z):
        differences = [(z[0]-a, z[1]-b) for a, b in fixed]
        if any(d in steps and d not in good for d in differences):
            return F(0)
        free = len({(a+z[0], b+z[1]) for a, b in support} - union)
        return multiplier*F(1, 2**free)

    values = {z: value(z) for z in near}
    base = multiplier*p
    if not 0 < base < 1 or any(x >= 1 for x in values.values()):
        return None
    minimum = (1-base)**degree
    for a, b in checked:
        product = F(1)
        for x, y in steps:
            z = (a+x, b+y)
            if allowed(z):
                product *= 1-values.get(z, base)
        minimum = min(minimum, product)
    if minimum < 1/multiplier:
        return None
    factor = F(1)
    for x in values.values():
        factor *= 1-x
    return dict(fixed_anchors=fixed, extra_avoided=extra,
                domain='all other anchors' if single else 'forward nonneighbors plus extras',
                multiplier=str(multiplier),
                active_changed_sites=[dict(anchor=z, x=str(x))
                                      for z, x in sorted(values.items()) if x],
                explicitly_checked_sites=len(checked),
                minimum_neighbor_product=str(minimum),
                conditional_avoidance_factor=str(factor),
                exact_test_passed=True)


def binary_two_by_three_minimum():
    h, w = 2, 3
    p, rho, alpha = F(1, 64), F(1, 32), F(1, 40)
    kp = p*exp_upper(15*alpha)/(1-rho*exp_upper(9*alpha))**2
    upper = p*rho*exp_upper(24*alpha)/(1-rho**2*exp_upper(18*alpha))
    assert kp <= alpha
    ds = displacements(rectangle(h, w))
    groups = binary_representatives(h, w)
    assert len(groups) == 18
    multipliers = [F(7, 5), F(3, 2), F(7, 4), F(2), F(5, 2)]
    records, minima = [], []
    for g, word in sorted(groups.items()):
        assert agreement(rectangle(h, w), word) == g
        good = [d for i, d in enumerate(ds) if g >> i & 1]
        if len(good) == 1 and (h-abs(good[0][0]))*(w-abs(good[0][1])) == 1:
            minima.append(g)
            continue
        b, one = F(0), None
        for multiplier in multipliers:
            record = conditioned_pattern_certificate(
                h, w, g, ((0, 0),), (), multiplier, single=True)
            if record:
                candidate = p*F(record['conditional_avoidance_factor'])
                if candidate > b:
                    b, one = candidate, record
        best, best_records, chosen_best = F(0), None, None
        for length in (1, 2, 3):
            for chosen in permutations(good, length):
                total, witnesses = F(0), []
                for j, d in enumerate(chosen):
                    pair_best, pair_record = F(0), None
                    for multiplier in multipliers:
                        record = conditioned_pattern_certificate(
                            h, w, g, ((0, 0), d), chosen[:j], multiplier, single=False)
                        if record:
                            assert all(F(row['x']) >= b for row in record['active_changed_sites'])
                            r = (h-abs(d[0]))*(w-abs(d[1]))
                            probability = F(1, 2**(2*h*w-r))
                            factor = F(record['conditional_avoidance_factor'])
                            c = min(probability, probability*factor/(1-b)**7)
                            if c > pair_best:
                                pair_best = c
                                pair_record = dict(**record, pair_probability=str(probability),
                                                   uniform_defect=str(c))
                    if pair_record is None:
                        break
                    total += pair_best
                    witnesses.append(pair_record)
                else:
                    if total > best:
                        best, best_records, chosen_best = total, witnesses, chosen
        gap = best/(1+best)-upper
        assert gap > 0, (word, gap)
        records.append(dict(word=word, agreement_mask=g, conditional_rate_lower_bound=str(b),
                            single_copy_certificate=one, ordered_neighbors=chosen_best,
                            witness_certificates=best_records, total_uniform_defect=str(best),
                            pressure_gap_coefficient=str(gap)))
    assert len(minima) == 2 and len(records) == 16
    return dict(status='Exact complete binary 2x3 pattern-class certificate. '
                       'The analytic reduction in the note gives strict planar '
                       'pressure minimality for every beta<0 and strict entropy minimality.',
                minimizing_agreement_masks=minima, kp_parameter=str(alpha),
                kp_sum_bound=str(kp), single_corner_upper_bound=str(upper),
                minimum_displayed_gap=float(min(F(row['pressure_gap_coefficient'])
                                                for row in records)), records=records)


def canonical(h, w, dx, dy):
    dx, dy = abs(dx), abs(dy)
    if h == w and dx > dy:
        dx, dy = dy, dx
    return h, w, dx, dy


def main():
    specs = [(2, 4, 1, 0, 2), (2, 4, 1, 1, 2), (2, 4, 0, 2, 2),
             (2, 5, 0, 3, 2), (2, 5, 1, 1, 2),
             (2, 5, 1, 0, 2), (2, 5, 0, 2, 2),
             (2, 6, 0, 3, 2), (2, 6, 1, 0, 2),
             (3, 4, 0, 2, 2), (3, 4, 1, 1, 2),
             (2, 2, 0, 1, 3), (2, 3, 1, 0, 3)]
    records = [certificate(*spec) for spec in specs]
    available = set(specs)
    uniform_cases, certified_cases = 0, 0
    # All nondegenerate rectangles of area below 16, up to transpose.
    for h in range(2, 4):
        for w in range(h, 16//h+1):
            if h*w >= 16:
                continue
            for q in (2, 3):
                if q == 2 and (h, w) in [(2, 2), (2, 3), (3, 3)]:
                    continue  # empty G is not realizable on these binary shapes
                for dx in range(h):
                    for dy in range(1-w, w):
                        if (dx, dy) <= (0, 0):
                            continue
                        overlap = (h-dx)*(w-abs(dy))
                        free = max(overlap, h*w-overlap)
                        # e<11/4, so this rational condition implies the
                        # uniform conditional local lemma for a maximal pair.
                        if q**free >= F(11, 4)*(2*h-1)*(2*w-1):
                            uniform_cases += 1
                        else:
                            key = canonical(h, w, dx, dy)+(q,)
                            assert key in available, key
                            certified_cases += 1
    # For s>=16 and rectangles, D+1<=4s. A two-step induction extends these.
    for s in (16, 17):
        assert 2**((s+1)//2) >= 11*s
    # For arbitrary supports, D+1<=s(s-1)+1; same two-step induction.
    for s in (19, 21, 22):
        assert 2**((s+1)//2) >= F(11, 4)*(s*(s-1)+1)
    assert 2**10 >= F(11, 4)*363  # s=20, maximal overlap r=10
    output = dict(status='Exact rational certificates; all-volume consequence '
                        'uses the maximal-overlap/local-lemma/Shearer proof '
                        'in entropy-ordering.md. No floating spectral input.',
                  small_shape_uniform_displacements=uniform_cases,
                  small_shape_certified_displacements=certified_cases,
                  rectangle_large_area_bases=[16, 17],
                  arbitrary_support_large_size_bases=[19, 20, 21, 22],
                  certificates=records,
                  exceptional_binary_square=binary_square_minimum(),
                  exceptional_binary_two_by_three=binary_two_by_three_minimum())
    Path(__file__).with_name('minimum-certificates.json').write_text(
        json.dumps(output, indent=2) + '\n')
    print('All 13 rational certificates passed; small-shape coverage:',
          uniform_cases, 'uniform cases and', certified_cases, 'certificate cases.')
    print('Binary 3x3 corner minimum certified; uniform gap coefficient >',
          output['exceptional_binary_square']['minimum_displayed_gap'])
    print('Binary 2x3 corner minimum certified; uniform gap coefficient >',
          output['exceptional_binary_two_by_three']['minimum_displayed_gap'])


if __name__ == '__main__':
    main()
