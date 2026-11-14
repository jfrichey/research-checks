#!/usr/bin/env python3
"""Exact checks of the large-sleep vacancy expansion in report.md.

Concrete questions: does the escape-probability recurrence agree with raw
symbolic stabilization, and can duplicated conditional-correlation motifs
violate its disjoint cylinder inequalities? Region volumes supply a second
check of the recurrence. The proof, rather than these cases, establishes
the general coefficient inequalities and the stated large-rate NA theorem.
For endpoint paths, the first correction is also checked: its proposed
consecutive-site interaction formula must agree with every vacancy moment,
including unequal sleep rates, before using it for a large-rate NA proof.
"""

import argparse
import itertools
import json
from pathlib import Path

from exact import graph_kernel, increasing_events, sp


THREE_SITE = {
    "path3": ((1, 3), (0, 2), (1, 3)),
    "triangle3": ((1, 2, 3), (0, 2, 3), (0, 1, 3)),
    "path3_sink_endpoint": ((1, 3), (0, 2), (1,)),
    "path3_sink_center": ((1,), (0, 2, 3), (1,)),
    "path3_sink_endpoint_center": ((1, 3), (0, 2, 3), (1,)),
    "path3_sink_all": ((1, 3), (0, 2, 3), (1, 3)),
    "triangle3_sink_one": ((1, 2, 3), (0, 2), (0, 1)),
    "triangle3_sink_two": ((1, 2, 3), (0, 2, 3), (0, 1)),
}


def coefficients(neighbors):
    n = len(neighbors)
    laplacian = sp.eye(n)
    for i, row in enumerate(graph_kernel(neighbors)):
        for j, prob in row:
            if j < n:
                laplacian[i, j] -= prob
    green = laplacian.inv()
    escape, coeff = {}, {0: sp.S.One}
    for mask in range(1, 1 << n):
        vertices = [i for i in range(n) if mask >> i & 1]
        minor = green.extract(vertices, vertices)
        eq = minor.inv()*sp.ones(len(vertices), 1)
        assert minor*eq == sp.ones(len(vertices), 1)
        escape[mask] = dict(zip(vertices, eq))
        assert all(v >= 0 for v in eq) and any(v > 0 for v in eq)
        value = 0
        sub = (mask-1) & mask
        while True:
            diff = [i for i in vertices if not sub >> i & 1]
            value += (-1)**(len(diff)+1)*coeff[sub]*sp.prod(escape[mask][i] for i in diff)
            if sub == 0:
                break
            sub = (sub-1) & mask
        coeff[mask] = sp.factor(value)
        assert coeff[mask] > 0
    return coeff, escape


def region_volume(mask, escape):
    """Independent exact integration over boxes cut out by all thresholds."""
    vertices = list(escape[mask])
    subsets = [b for b in escape if b & mask == b]
    cuts = {v: sorted({sp.S.Zero} | {escape[b][v] for b in subsets if v in escape[b]})
            for v in vertices}
    intervals = [list(zip(cuts[v][:-1], cuts[v][1:])) for v in vertices]
    volume = sp.S.Zero
    for box in itertools.product(*intervals):
        point = {v: (lo+hi)/2 for v, (lo, hi) in zip(vertices, box)}
        if all(any(point[v] < threshold for v, threshold in escape[b].items())
               for b in subsets):
            volume += sp.prod(hi-lo for lo, hi in box)
    return volume


def check(neighbors, symbolic=None, integrate=False):
    n = len(neighbors)
    coeff, escape = coefficients(neighbors)
    if symbolic is not None:
        p = sp.Symbol("p")
        law = tuple(sp.sympify(v) for v in symbolic)
        for mask in range(1, 1 << n):
            joint = sum(law[state] for state in range(1 << n) if state & mask == 0)
            leading = sp.cancel(joint/(1-p)**mask.bit_count()).subs(p, 1)
            assert leading == coeff[mask]
    if integrate:
        for mask in range(1, 1 << n):
            assert region_volume(mask, escape) == coeff[mask]
    gaps = {}
    for left in range(1, 1 << n):
        for right in range(left+1, 1 << n):
            if left & right:
                continue
            gap = coeff[left]*coeff[right]-coeff[left | right]
            assert gap >= 0
            gaps[f"{left},{right}"] = gap
    pair_escape = all(escape[(1 << i) | (1 << j)][i] > 0
                      and escape[(1 << i) | (1 << j)][j] > 0
                      for i in range(n) for j in range(i+1, n))
    # These cases have connected internal graphs. The theorem additionally
    # needs irreducibility when inferring strictness from pair escape.
    if pair_escape:
        assert all(value > 0 for value in gaps.values())
    return {"neighbors": neighbors, "symbolic_stabilization_check": symbolic is not None,
            "independent_region_volume_check": integrate, "pair_escape_condition": pair_escape,
            "coefficients": {str(k): str(v) for k, v in coeff.items()},
            "disjoint_gaps": {k: str(v) for k, v in gaps.items()}}


def adjacency(n, edges, sinks):
    return tuple(tuple([b if a == v else a for a, b in edges if v in (a, b)]
                       + ([n] if v in sinks else [])) for v in range(n))


def check_path_correction(n, alpha=None, symbolic=None, conductances=None, holding=None):
    from algebra import high_rate_series, cylinder_gap_series

    alpha = tuple(map(sp.Rational, alpha or (1,)*n))
    conductances = tuple(map(sp.Rational, conductances or (1,)*n))
    holding = tuple(map(sp.Rational, holding or (0,)*n))
    resistance = [sum(1/g for g in conductances[:i+1]) for i in range(n)]
    mass = [conductances[i]+(conductances[i+1] if i+1 < n else 0)+holding[i]
            for i in range(n)]
    kernel = tuple(tuple([(i-1 if i else n, conductances[i]/mass[i])]
                         + ([(i+1, conductances[i+1]/mass[i])] if i+1 < n else [])
                         + ([(i, holding[i]/mass[i])] if holding[i] else []))
                   for i in range(n))
    c = [1/(mass[i]*resistance[i]) for i in range(n)]
    law = high_rate_series(kernel, n+1, alpha)
    relative, leading = {0: sp.S.Zero}, {0: sp.S.One}
    for mask in range(1, 1 << n):
        vertices = [i for i in range(n) if mask >> i & 1]
        leading[mask] = sp.prod(c[i]/alpha[i] for i in vertices)
        relative[mask] = sum(
            -c[i]/alpha[i] + sum(resistance[i]*c[w]/(alpha[w]*(resistance[w]-resistance[i]))
                                 for w in range(i+1, nxt))
            for i, nxt in zip(vertices, vertices[1:]+[n]))
        moment = [sum(row[state] for state in range(1 << n) if state & mask == 0)
                  for row in law]
        order = len(vertices)
        assert all(value == 0 for value in moment[:order])
        assert moment[order] == leading[mask]
        assert moment[order+1] == leading[mask]*relative[mask]
    checked = 0
    for left in range(1, 1 << n):
        for right in range(left+1, 1 << n):
            if left & right:
                continue
            gap = cylinder_gap_series(law, left, right)
            order = (left | right).bit_count()+1
            predicted = leading[left | right]*(relative[left]+relative[right]
                                               -relative[left | right])
            assert predicted > 0
            assert all(value == 0 for value in gap[:order])
            assert gap[order] == predicted
            checked += 1
    if symbolic is not None:
        assert alpha == (sp.S.One,)*n
        assert conductances == (sp.S.One,)*n and not any(holding)
        p, t = sp.symbols("p t")
        raw = [sp.sympify(v) for v in symbolic]
        for mask in range(1, 1 << n):
            moment = sum(raw[state] for state in range(1 << n) if state & mask == 0)
            normalized = sp.cancel(moment.subs(p, 1/(1+t))/t**mask.bit_count())
            numerator, denominator = map(lambda x: sp.Poly(x, t), sp.fraction(normalized))
            first = numerator.nth(0)/denominator.nth(0)
            second = (numerator.nth(1)-first*denominator.nth(1))/denominator.nth(0)
            assert first == leading[mask]
            assert second == leading[mask]*relative[mask]
    interactions = {f"{i},{j}": resistance[i]*sum(c[w]/(alpha[w]*(resistance[w]-resistance[i]))
                                                 for w in range(j, n))
                    for i in range(n) for j in range(i+1, n)}
    return {"n": n, "alpha": list(map(str, alpha)),
            "conductances": list(map(str, conductances)), "holding": list(map(str, holding)),
            "all_moments_checked": (1 << n)-1, "strict_disjoint_gaps_checked": checked,
            "independent_symbolic_check": symbolic is not None,
            "relative_corrections": {str(k): str(v) for k, v in relative.items()},
            "consecutive_interactions": {k: str(v) for k, v in interactions.items()}}


def check_cut_vertex_corrections(neighbors, alpha):
    """Test the proposed general correction on branching/cyclic cut graphs.

    The new formula must match the independent algebraic series for every
    dominated subset, and resolve all zero leading gaps. Full Boolean-event
    checks also test the next-order cancellation used to deduce NA.
    """
    from algebra import high_rate_series

    n = len(neighbors)
    alpha = tuple(map(sp.Rational, alpha))
    kernel = graph_kernel(neighbors)
    coeff, escape = coefficients(neighbors)
    laplacian = sp.eye(n)
    for v, row in enumerate(kernel):
        for u, prob in row:
            if u < n:
                laplacian[v, u] -= prob
    green = laplacian.inv()
    inverses = {}
    for mask in range(1, 1 << n):
        sites = [i for i in range(n) if mask >> i & 1]
        inverses[mask] = (sites, green.extract(sites, sites).inv())

    def hit(v, u, taboo=0):
        sites, inverse = inverses[(1 << v) | (1 << u) | taboo]
        return -inverse[sites.index(v), sites.index(u)]

    law = high_rate_series(kernel, n+1, alpha)
    relative = {0: sp.S.Zero}
    for mask in range(1, 1 << n):
        moment = [sum(row[state] for state in range(1 << n) if state & mask == 0)
                  for row in law]
        leading = coeff[mask]/sp.prod(alpha[i] for i in range(n) if mask >> i & 1)
        assert moment[mask.bit_count()] == leading
        relative[mask] = moment[mask.bit_count()+1]/leading
    penalties, monotone_checks, zero_gaps = {}, 0, 0
    for v in range(n):
        descendants = sum(1 << u for u in range(n) if u != v and green[u, v] == green[v, v])
        cv = coeff[1 << v]
        for block in range(1, 1 << n):
            if block & descendants != block:
                continue
            terms = []
            for u in range(n):
                if u == v:
                    continue
                cu = coeff[1 << u]
                if not descendants >> u & 1:
                    hh = green[u, v]/(green[v, v]-green[u, v])
                    ev = escape[(1 << v) | (1 << u)][v]
                    term = hh*(cv-ev)*(cu-coeff[block | (1 << u)]/coeff[block])/(cv*alpha[u])
                elif block >> u & 1:
                    term = cu*hit(v, u)/(cv*alpha[u])
                else:
                    term = (cu*hit(v, u)-coeff[block | (1 << u)]/coeff[block]
                            *hit(v, u, block))/(cv*alpha[u])
                assert term >= 0
                terms.append(term)
            penalty = sum(terms)
            assert penalty == relative[1 << v]+relative[block]-relative[block | (1 << v)]
            assert penalty > 0
            penalties[v, block] = penalty
        for block in range(1, 1 << n):
            for other in range(block+1, 1 << n):
                if block & other or (block | other) & descendants != block | other:
                    continue
                if coeff[block]*coeff[other] == coeff[block | other]:
                    assert penalties[v, block | other] >= penalties[v, block]
                    assert penalties[v, block | other] >= penalties[v, other]
                    monotone_checks += 2
    for left in range(1, 1 << n):
        for right in range(left+1, 1 << n):
            if not left & right and coeff[left]*coeff[right] == coeff[left | right]:
                assert relative[left]+relative[right]-relative[left | right] > 0
                zero_gaps += 1
    event_cache, comparisons, orders = {}, 0, {}

    def event_series(event):
        if event not in event_cache:
            # Event bits encode vacancy states, while law encodes sleepers.
            event_cache[event] = tuple(sum(row[state] for state in range(1 << n)
                                           if event >> (((1 << n)-1) ^ state) & 1)
                                       for row in law)
        return event_cache[event]

    lifted = {}
    for mask in range(1, 1 << n):
        if mask.bit_count() == n:
            continue
        sites = tuple(i for i in range(n) if mask >> i & 1)
        lifted[mask] = [sum(1 << state for state in range(1 << n)
                            if sum(((state >> i) & 1) << j for j, i in enumerate(sites)) in event)
                        for event in increasing_events(sites)]
    for left, events in lifted.items():
        for right, other_events in lifted.items():
            if right <= left or left & right:
                continue
            for event in events:
                aa = event_series(event)
                for other in other_events:
                    bb, joint = event_series(other), event_series(event & other)
                    gaps = [sum(aa[j]*bb[k-j] for j in range(k+1))-joint[k] for k in range(n+2)]
                    first = next((k, value) for k, value in enumerate(gaps) if value)
                    assert first[1] > 0
                    orders[first[0]] = orders.get(first[0], 0)+1
                    comparisons += 1
    return {"neighbors": neighbors, "alpha": list(map(str, alpha)),
            "dominated_block_identities_checked": len(penalties),
            "monotonicity_checks": monotone_checks, "zero_leading_cylinder_gaps_resolved": zero_gaps,
            "full_NA_event_comparisons": comparisons, "positive_event_gap_orders": orders,
            "relative_corrections": {str(k): str(v) for k, v in relative.items()},
            "dominated_block_penalties": {f"{v},{mask}": str(value)
                                          for (v, mask), value in penalties.items()}}


def conditional_tree_obstruction():
    """Disprove a Rayleigh/CNA route even for a tree with one sink edge.

    The leading conditional covariance uses the proved vacancy theorem;
    region integration and the IDLA series independently verify its inputs.
    The finite check also tests full NA, which conditioning does not settle.
    """
    from algebra import high_rate_series, stationary
    from exact import na_gaps

    neighbors = adjacency(5, [(0, 1), (1, 2), (1, 3), (3, 4)], {0})
    leading, escape = coefficients(neighbors)
    masks = (4, 12, 20, 28)  # {2}, {2,3}, {2,4}, {2,3,4}
    expected = (sp.Rational(1, 3), sp.Rational(7, 150),
                sp.Rational(7, 96), sp.Rational(17, 1600))
    assert tuple(leading[mask] for mask in masks) == expected
    series = high_rate_series(graph_kernel(neighbors), 3)
    for mask in masks:
        assert region_volume(mask, escape) == leading[mask]
        assert sum(series[mask.bit_count()][state] for state in range(32)
                   if state & mask == 0) == leading[mask]
    coefficient = leading[28]/leading[4]-leading[12]*leading[20]/leading[4]**2
    assert coefficient == sp.Rational(1, 800)
    law = stationary(graph_kernel(neighbors), (sp.Rational(1, 64),)*5)
    moments = {mask: sum(law[state] for state in range(32) if state & mask == 0)
               for mask in masks}
    covariance = moments[28]/moments[4]-moments[12]*moments[20]/moments[4]**2
    assert covariance > 0
    gaps = list(na_gaps(law))
    assert all(record[-1] >= 0 for record in gaps)
    return {"neighbors": neighbors, "condition_empty": [2], "pair": [3, 4],
            "vacancy_leading_coefficients": {str(mask): str(leading[mask]) for mask in masks},
            "region_and_IDLA_series_agree": True,
            "conditional_covariance_leading_order": 2,
            "conditional_covariance_leading_coefficient": str(coefficient),
            "finite_witness_lambda": 64, "finite_conditional_covariance": str(covariance),
            "finite_unconditioned_NA_comparisons": len(gaps),
            "finite_unconditioned_NA_passes": True}


def check_vacancy_bk_series(neighbors):
    """Check the first nonzero BK gap directly from independent IDLA series.

    This validates the sparse BK argument on all increasing-event pairs
    at at most four sites, including overlapping events. Only its first
    two possible orders are used, as in the abstract criterion.
    """
    from algebra import high_rate_series

    n = len(neighbors)
    assert n <= 4
    series = high_rate_series(graph_kernel(neighbors), n+1)
    denominator = sp.ilcm(*(entry.q for row in series for entry in row))
    scale = int(denominator)
    rows = [[int(denominator*row[((1 << n)-1) ^ mask]) for mask in range(1 << n)]
            for row in series]
    events = increasing_events(tuple(range(n)))
    triggers = [tuple(mask for mask in event
                      if not any(other != mask and other & mask == other for other in event))
                for event in events]
    minimum = [min(mask.bit_count() for mask in choices) for choices in triggers]
    values = [[sum(row[mask] for mask in event) for row in rows] for event in events]
    cylinder_bits = [sum(1 << state for state in range(1 << n) if state & mask == mask)
                     for mask in range(1 << n)]
    cache = {}
    orders, checked = {}, 0
    for left, aa in enumerate(triggers):
        a = minimum[left]
        for right in range(left, len(triggers)):
            b = minimum[right]
            bits = 0
            for first in aa:
                for second in triggers[right]:
                    if not first & second:
                        bits |= cylinder_bits[first | second]
            if bits not in cache:
                cache[bits] = [sum(row[mask] for mask in range(1 << n) if bits >> mask & 1)
                               for row in rows]
            joint = cache[bits]
            assert not any(joint[:a+b])
            leading = values[left][a]*values[right][b]
            if a+b <= n+1:
                leading -= scale*joint[a+b]
            assert leading >= 0
            order = a+b
            if leading == 0:
                correction = (values[left][a+1]*values[right][b]
                              +values[left][a]*values[right][b+1])
                if a+b+1 <= n+1:
                    correction -= scale*joint[a+b+1]
                assert correction > 0
                order += 1
            orders[str(order)] = orders.get(str(order), 0)+1
            checked += 1
    return {"neighbors": neighbors, "indicators": "vacancies", "parameter": "t=1/lambda",
            "all_increasing_event_pairs": checked, "positive_gap_orders": orders,
            "first_two_orders_suffice": True}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", help="Optional JSON certificate file")
    args = parser.parse_args()
    known = json.loads(Path(__file__).with_name("certificates.json").read_text())
    result = {}
    for name, neighbors in THREE_SITE.items():
        result[name] = check(neighbors, known[name]["law_in_p"], integrate=True)
    motifs = [
        ("single_sink_path", 4, [(0, 1), (1, 2), (2, 3)], {0}),
        ("two_sink_edges_path", 4, [(0, 1), (1, 2), (2, 3)], {0, 3}),
        ("duplicated_gate", 4, [(0, 1), (1, 2), (1, 3)], {0, 1, 3}),
        ("connected_gate_twins", 4, [(0, 1), (1, 2), (1, 3), (0, 3)], {0, 1, 3}),
        ("shared_gate", 5, [(0, 1), (1, 2), (0, 3), (3, 4)], {0, 1, 3}),
        ("bridged_motifs", 6, [(0, 1), (1, 2), (3, 4), (4, 5), (0, 3)], {0, 1, 3, 4}),
    ]
    for name, n, edges, sinks in motifs:
        result[name] = check(adjacency(n, edges, sinks), integrate=n <= 4)
    for name, case in result.items():
        gaps = [sp.Rational(v) for v in case["disjoint_gaps"].values()]
        print(name, "strict gaps", sum(bool(v > 0) for v in gaps),
              "zero gaps", sum(bool(v == 0) for v in gaps), flush=True)
    result["endpoint_path_first_corrections"] = []
    for n, alpha in [(n, None) for n in range(2, 7)] + [(4, (1, 2, 3, 5))]:
        symbolic = known["path3_sink_endpoint"]["law_in_p"] if n == 3 else None
        case = check_path_correction(n, alpha, symbolic)
        result["endpoint_path_first_corrections"].append(case)
        print("endpoint path", n, "alpha", case["alpha"], "first correction verified", flush=True)
    case = check_path_correction(4, (1, 2, 3, 5), conductances=(2, 3, 5, 7), holding=(1, 2, 1, 3))
    result["endpoint_path_first_corrections"].append(case)
    print("weighted endpoint path with holding: first correction verified", flush=True)
    result["cut_vertex_first_corrections"] = {}
    for name, edges, alpha in [
        ("rooted_fork", [(0, 1), (0, 2), (1, 3), (2, 4)], (1,)*5),
        ("cycle_below_cut", [(0, 1), (1, 2), (1, 3), (2, 3), (3, 4)], (1, 2, 3, 5, 7)),
    ]:
        case = check_cut_vertex_corrections(adjacency(5, edges, {0}), alpha)
        result["cut_vertex_first_corrections"][name] = case
        print(name, "cut-vertex formula verified; full event gaps", case["full_NA_event_comparisons"], flush=True)
    result["conditional_tree_obstruction"] = conditional_tree_obstruction()
    print("five-site tree: conditional covariance is positive; all unconditioned NA gaps pass at lambda=64", flush=True)
    result["vacancy_BK_series"] = {}
    for name, neighbors in (
            ("endpoint_path4", ((1, 4), (0, 2), (1, 3), (2,))),
            ("stemmed_star4", ((1, 4), (0, 2, 3), (1,), (1,)))):
        result["vacancy_BK_series"][name] = check_vacancy_bk_series(neighbors)
        print(name, "all 13,861 vacancy BK event series gaps pass", flush=True)
    if args.output:
        Path(args.output).write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()
