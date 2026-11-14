#!/usr/bin/env python3
"""Exact BK and NA certificates on conditional-correlation tree obstructions.

Question: does stationary BK survive where conditional negative correlation
fails? Check every pair of increasing events, including overlapping pairs,
at two specified rates and for both sleepers and vacancies. These finite
checks do not establish BK at other rates or on arbitrary graphs. Separate
NA certificates cover every positive common rate on the specified six-site
tree with two pendant length-two arms.
"""

import argparse
from functools import lru_cache
import json
from pathlib import Path

from algebra import stationary
from exact import graph_kernel, na_gaps, sp


@lru_cache(None)
def upsets(n):
    """Truth tables of all increasing events, including the two constants."""
    if n == 0:
        return (0, 1)
    previous = upsets(n-1)
    width = 1 << (n-1)
    return tuple(low | (high << width) for low in previous for high in previous
                 if low & ~high == 0)


def minimal_witnesses(bits, n):
    return tuple(state for state in range(1 << n) if bits >> state & 1
                 and not any(bits >> (state ^ (1 << i)) & 1
                             for i in range(n) if state >> i & 1))


def occurrence(left, right_bits, zero_bits):
    """Truth table of disjoint occurrence from one event's minimal witnesses.

    Right occurs outside witness I iff its truth table holds with I zeroed.
    Shifting the zero-I states by integer I sets precisely those bits to one.
    """
    value = 0
    for witness in left:
        value |= (right_bits & zero_bits[witness]) << witness
    return value


def certificate(law):
    n = len(law).bit_length()-1
    assert n <= 5 and len(law) == 1 << n
    law = tuple(map(sp.Rational, law))
    assert sum(law) == 1 and min(law) >= 0
    denominator = int(sp.ilcm(*(value.q for value in law)))
    weights = tuple(int(denominator*value) for value in law)
    events = upsets(n)
    probabilities = {bits: sum(weights[state] for state in range(1 << n)
                              if bits >> state & 1) for bits in events}
    scaled = {bits: denominator*value for bits, value in probabilities.items()}
    zero_bits = [sum(1 << state for state in range(1 << n) if state & mask == 0)
                 for mask in range(1 << n)]
    tails = [denominator*sum(weights[state] for state in range(1 << n)
                             if state.bit_count() >= k) for k in range(2*n+1)]
    full = (1 << (1 << n))-1
    records = []
    for bits in events:
        if bits in (0, full):
            continue
        witnesses = minimal_witnesses(bits, n)
        records.append((bits, probabilities[bits], witnesses,
                        min(mask.bit_count() for mask in witnesses)))
    checked, calculated = 0, 0
    result = {"sites": n, "law": list(map(str, law)),
              "nonconstant_increasing_events": len(records)}
    for index, (aa, pa, wa, ma) in enumerate(records):
        for bb, pb, wb, mb in records[index:]:
            checked += 1
            product = pa*pb
            # Every disjoint witness has at least ma+mb occupied sites.
            if product >= tails[ma+mb]:
                continue
            calculated += 1
            joint = (occurrence(wa, bb, zero_bits) if len(wa) <= len(wb)
                     else occurrence(wb, aa, zero_bits))
            if product < scaled[joint]:
                result.update({"BK_passes": False,
                               "pairs_checked_before_failure": checked,
                               "left_minimal_witnesses": wa,
                               "right_minimal_witnesses": wb,
                               "left_probability": str(sp.Rational(pa, denominator)),
                               "right_probability": str(sp.Rational(pb, denominator)),
                               "disjoint_probability": str(sp.Rational(probabilities[joint], denominator)),
                               "gap": str(sp.Rational(product-scaled[joint], denominator**2))})
                return result
    assert checked == len(records)*(len(records)+1)//2
    result.update({"BK_passes": True, "all_event_pairs": checked,
                   "pairs_using_exact_rank_tail_bound": checked-calculated,
                   "pairs_using_disjoint_occurrence": calculated})
    return result


def na_certificate(law):
    """Check NA through all bipartitions; unused coordinates may be padded.

    This integer truth-table check permits six sites, where only one side
    of a bipartition can require the 7,581 increasing five-site events.
    Constants are omitted. A failure records both local event truth tables.
    """
    n = len(law).bit_length()-1
    assert len(law) == 1 << n and 2 <= n <= 6
    law = tuple(map(sp.Rational, law))
    assert sum(law) == 1 and min(law) >= 0
    denominator = int(sp.ilcm(*(value.q for value in law)))
    weights = [int(denominator*value) for value in law]
    full, checked = (1 << n)-1, 0

    def probability(bits):
        value = 0
        while bits:
            low = bits & -bits
            value += weights[low.bit_length()-1]
            bits -= low
        return value

    for left in range(1, full):
        right = full ^ left
        if left > right:
            continue
        sides = []
        for side in (left, right):
            coordinates = [i for i in range(n) if side >> i & 1]
            projections = [sum(((mask >> i) & 1) << j for j, i in enumerate(coordinates))
                           for mask in range(1 << n)]
            events = []
            for bits in upsets(len(coordinates)):
                if bits in (0, (1 << (1 << len(coordinates)))-1):
                    continue
                lifted = sum(1 << mask for mask, local in enumerate(projections) if bits >> local & 1)
                events.append((lifted, probability(lifted), bits))
            sides.append(events)
        for aa, pa, event_a in sides[0]:
            for bb, pb, event_b in sides[1]:
                checked += 1
                gap = pa*pb-denominator*probability(aa & bb)
                if gap < 0:
                    return {"NA_passes": False, "pairs_checked": checked,
                            "left_support_mask": left, "right_support_mask": right,
                            "left_event_bits": event_a, "right_event_bits": event_b,
                            "gap": str(sp.Rational(gap, denominator**2))}
    return {"NA_passes": True, "bipartition_event_pairs": checked}


def balanced_tree_na_certificate():
    """Can two pendant length-two arms amplify conditional dependence?

    Check all disjoint increasing events on this one specified six-site
    tree at rates 1 and 64. These finite checks are not an all-rate proof.
    """
    assert na_certificate(tuple(sp.Rational(x, 1000)
                                for x in (200, 240, 240, 294, 6, 10, 10, 0)))["NA_passes"]
    assert not na_certificate(tuple(map(sp.Rational,
                                      ("5/78", "19/390", "68/195", "7/13"))))["NA_passes"]
    neighbors = ((1, 6), (0, 2, 4), (1, 3), (2,), (1, 5), (4,))
    records = []
    for rate in (1, 64):
        law = stationary(graph_kernel(neighbors), (sp.Rational(1, rate),)*6)
        result = na_certificate(law)
        assert result["NA_passes"] and result["bipartition_event_pairs"] == 58674
        records.append({"lambda": rate, "law": list(map(str, law)), **result})
    return {"neighbors": neighbors,
            "scope": "specified tree and two common rates only", "cases": records,
            "NA_checker_accepts_known_NA_law_and_detects_insertion_counterexample": True}


def symbolic_na_gaps(weights, denominator):
    """Coefficientwise NA on all bipartitions, with polynomial subset sums.

    For each event on the smaller block, compute the signed gap measure
    on the larger block once. Removing the least true state of an upset
    leaves another upset, so every event gap then needs one polynomial
    addition. Negative coefficients only fail this sufficient criterion.
    """
    n, variable = len(weights).bit_length()-1, denominator.gens[0]
    assert 2 <= n <= 6 and len(weights) == 1 << n
    zero = sp.Poly(0, variable)
    assert sum(weights, zero) == denominator
    full, checked, zeros = (1 << n)-1, 0, 0
    positive_denominator = denominator != zero and all(value >= 0 for value in denominator.all_coeffs())
    for left in range(1, full):
        right = full ^ left
        if left > right:
            continue
        aa, bb = sorted((left, right), key=lambda side: (side.bit_count(), side))
        coordinates = [[i for i in range(n) if side >> i & 1] for side in (aa, bb)]
        a_size, b_size = (1 << len(co) for co in coordinates)
        joint = [[zero]*b_size for _ in range(a_size)]
        for mask, weight in enumerate(weights):
            a, b = [sum(((mask >> i) & 1) << j for j, i in enumerate(co)) for co in coordinates]
            joint[a][b] = weight
        marginal_b = [sum((joint[a][b] for a in range(a_size)), zero) for b in range(b_size)]
        for event_a in upsets(len(coordinates[0])):
            if event_a in (0, (1 << a_size)-1):
                continue
            conditional = [sum((joint[a][b] for a in range(a_size) if event_a >> a & 1), zero)
                           for b in range(b_size)]
            probability_a = sum(conditional, zero)
            signed = [probability_a*marginal_b[b]-denominator*conditional[b] for b in range(b_size)]
            assert sum(signed, zero).is_zero
            gaps = {0: zero}
            for event_b in sorted(upsets(len(coordinates[1]))):
                if event_b == 0:
                    continue
                low = event_b & -event_b
                gap = gaps[event_b ^ low]+signed[low.bit_length()-1]
                gaps[event_b] = gap
                if event_b == (1 << b_size)-1:
                    assert gap.is_zero
                    continue
                checked += 1
                zeros += gap.is_zero
                negatives = sum(1 for value in gap.all_coeffs() if value < 0)
                if negatives:
                    return {"coefficientwise_NA_certified": False, "pairs_checked": checked,
                            "left_support_mask": aa, "right_support_mask": bb,
                            "left_event_bits": event_a, "right_event_bits": event_b,
                            "negative_coefficients": negatives,
                            "scope": "failure of a sufficient coefficient criterion, not an NA counterexample"}
    return {"coefficientwise_NA_certified": positive_denominator,
            "positive_coefficient_denominator": positive_denominator,
            "bipartition_event_pairs": checked, "identically_zero_gaps": zeros}


def balanced_tree_symbolic_na_certificate():
    """All-rate check on the six-site tree with two pendant length-two arms."""
    from algebra import stationary_orbit_polynomials

    t = sp.Symbol("t")
    product = symbolic_na_gaps([sp.Poly(t**mask.bit_count(), t) for mask in range(8)],
                              sp.Poly((1+t)**3, t))
    assert product["coefficientwise_NA_certified"] and product["identically_zero_gaps"] == 12
    counter = symbolic_na_gaps([sp.Poly(value, t) for value in (25, 19, 136, 210)], sp.Poly(390, t))
    assert not counter["coefficientwise_NA_certified"]
    neighbors = ((1, 6), (0, 2, 4), (1, 3), (2,), (1, 5), (4,))
    permutation = (0, 1, 4, 5, 2, 3)
    orbits = {tuple(sorted({mask, sum(((mask >> i) & 1) << permutation[i] for i in range(6))}))
              for mask in range(64)}
    weights, denominator = stationary_orbit_polynomials(graph_kernel(neighbors), orbits, t)
    print("six-site symbolic law: full basis residual verified; denominator degree",
          denominator.degree(), flush=True)
    for rate in (1, 64):
        point = sp.Rational(1, rate)
        finite = stationary(graph_kernel(neighbors), (point,)*6)
        assert tuple(weight.eval(point)/denominator.eval(point) for weight in weights) == finite
    result = symbolic_na_gaps(weights, denominator)
    return {"neighbors": neighbors, "parameter": "t=1/lambda", "symmetry_permutation": permutation,
            "stable_state_orbits": len(orbits), "full_basis_residual_verified": True,
            "normalization_verified": True, "rational_laws_at_rates_1_and_64_verified": True,
            "denominator_degree": denominator.degree(),
            "denominator_descending_coefficients": list(map(str, denominator.all_coeffs())),
            "orbit_numerators": [{"masks": orbit,
                                  "per_state_descending_coefficients": list(map(str, weights[orbit[0]].all_coeffs()))}
                                 for orbit in sorted(orbits)],
            "gap_certificate": result,
            "polynomial_checker_accepts_product_law_and_detects_insertion_counterexample": True}


def complement_defect_check(n=4):
    """Check the strict-containment lemma used for sparse BK complementation.

    G(Y)=not E(1-Y), H(Y)=not F(1-Y). Every state in
    (G intersect H) minus (G square H) must strictly contain a state in
    the complement-image of (E intersect F) minus (E square F).
    The report proves this for arbitrary n; enumeration checks the encoding.
    """
    size = 1 << n
    full = (1 << size)-1
    zeros = [sum(1 << state for state in range(size) if state & mask == 0)
             for mask in range(size)]
    strict_cylinders = [sum(1 << state for state in range(size)
                            if state != mask and state & mask == mask)
                        for mask in range(size)]
    events = tuple(bits for bits in upsets(n) if bits not in (0, full))
    witnesses = {bits: minimal_witnesses(bits, n) for bits in events}
    dual = {bits: sum(1 << state for state in range(size)
                     if not (bits >> ((size-1) ^ state) & 1)) for bits in events}
    checked, nonempty = 0, 0
    for index, aa in enumerate(events):
        for bb in events[index:]:
            gg, hh = dual[aa], dual[bb]
            defect = aa & bb & ~occurrence(witnesses[aa], bb, zeros)
            qq = sum(1 << state for state in range(size)
                     if defect >> ((size-1) ^ state) & 1)
            dd = gg & hh & ~occurrence(witnesses[gg], hh, zeros)
            closure = 0
            for state in range(size):
                if qq >> state & 1:
                    closure |= strict_cylinders[state]
            assert dd & ~closure == 0
            checked += 1
            nonempty += bool(dd)
    return {"sites": n, "all_event_pairs": checked,
            "pairs_with_nonempty_dual_defect": nonempty,
            "strict_containment_always_holds": True}


def symbolic_certificate(neighbors):
    """Test coefficientwise BK positivity at every common sleep rate.

    This stronger algebraic sufficient condition is tested only on the
    specified graph. A negative coefficient alone would not refute BK.
    """
    t = sp.Symbol("t")
    n = len(neighbors)
    assert n <= 4
    law = stationary(graph_kernel(neighbors), (t,)*n)
    fractions = [tuple(sp.Poly(part, t) for part in sp.fraction(mass)) for mass in law]
    denominator = sp.Poly(1, t)
    for _, value in fractions:
        denominator = sp.lcm(denominator, value)
    if denominator.LC() < 0:
        denominator = -denominator
    weights = [num*denominator.exquo(den) for num, den in fractions]
    zero = sp.Poly(0, t)
    assert sum(weights, zero) == denominator
    positive = all(value >= 0 for value in denominator.all_coeffs())
    events = upsets(n)
    full = (1 << (1 << n))-1
    chosen = tuple(bits for bits in events if bits not in (0, full))
    zeros = [sum(1 << state for state in range(1 << n) if state & mask == 0)
             for mask in range(1 << n)]
    witnesses = {bits: minimal_witnesses(bits, n) for bits in events}
    result = {"neighbors": neighbors, "parameter": "t=1/lambda",
              "denominator": str(sp.factor(denominator.as_expr())),
              "sleeper_numerators": [str(sp.factor(value.as_expr())) for value in weights],
              "denominator_positive_coefficients": positive, "indicators": {}}
    for orientation, ww in (("sleepers", weights), ("vacancies", weights[::-1])):
        values = {bits: sum((ww[state] for state in range(1 << n) if bits >> state & 1), zero)
                  for bits in events}
        scaled = {bits: denominator*value for bits, value in values.items()}
        checked, failures = 0, []
        for index, aa in enumerate(chosen):
            for bb in chosen[index:]:
                joint = occurrence(witnesses[aa], bb, zeros)
                gap = values[aa]*values[bb]-scaled[joint]
                checked += 1
                if any(value < 0 for value in gap.all_coeffs()):
                    failures.append({"left_minimal_witnesses": witnesses[aa],
                                     "right_minimal_witnesses": witnesses[bb],
                                     "gap_numerator": str(sp.factor(gap.as_expr()))})
        result["indicators"][orientation] = {"all_event_pairs": checked,
                                              "failed_coefficient_tests": failures,
                                              "all_rates_BK_certified": positive and not failures}
    return result


def product_neighborhood_certificate():
    """Check strict product BK for all four-site pairs at unequal parameters.

    Overlapping essential supports have gap at least the square of the
    smallest atom. The remaining pairs have exactly zero gap. This is a
    finite check of the lemma used for the NA neighborhood corollary.
    """
    n = 4
    parameters = (sp.Rational(1, 2), sp.Rational(2, 3),
                  sp.Rational(3, 4), sp.Rational(5, 6))
    law = tuple(sp.prod(parameters[i] if state >> i & 1 else 1-parameters[i]
                        for i in range(n)) for state in range(1 << n))
    denominator = int(sp.ilcm(*(value.q for value in law)))
    weights = tuple(int(value*denominator) for value in law)
    minimum = min(weights)
    zeros = [sum(1 << state for state in range(1 << n) if state & mask == 0)
             for mask in range(1 << n)]
    events = upsets(n)
    probabilities = {bits: sum(weights[state] for state in range(1 << n) if bits >> state & 1)
                     for bits in events}
    records = []
    for bits in events:
        if bits in (0, (1 << (1 << n))-1):
            continue
        witnesses = minimal_witnesses(bits, n)
        support = 0
        for witness in witnesses:
            support |= witness
        records.append((bits, witnesses, support))
    overlap, disjoint = 0, 0
    for index, (aa, wa, sa) in enumerate(records):
        for bb, wb, sb in records[index:]:
            joint = occurrence(wa, bb, zeros)
            gap = probabilities[aa]*probabilities[bb]-denominator*probabilities[joint]
            if sa & sb:
                assert gap >= minimum**2
                overlap += 1
            else:
                assert gap == 0
                disjoint += 1
    assert overlap+disjoint == 13861
    return {"parameters": list(map(str, parameters)), "all_event_pairs": overlap+disjoint,
            "overlapping_essential_support_pairs": overlap,
            "disjoint_essential_support_pairs": disjoint,
            "minimum_atom": str(sp.Rational(minimum, denominator)),
            "overlap_gaps_at_least_minimum_atom_squared": True,
            "zero_gaps_exactly_disjoint_essential_supports": True}


def two_copy_witness_obstruction():
    """A saturated ARW law violates disjoint-witness transfer between copies.

    This tests a precise proposed switching target, stronger than ordinary
    BK. The complete graph with three non-sink vertices discriminates:
    its large-rate vacancy covariance predicts failure of the strengthening,
    while the ordinary BK gap for the same events remains positive.
    """
    from exact import stationary_law

    kernel = graph_kernel(((1, 2, 3), (0, 2, 3), (0, 1, 3)))
    law, counts = stationary_law(kernel, (sp.Rational(2, 3),)*3)
    expected = tuple(sp.Rational(weight, 78064) for weight in
                     (1051, 3747, 3747, 11592, 3747, 11592, 11592, 30996))
    assert law == expected == stationary(kernel, (sp.Rational(1, 2),)*3)
    left, right = (1, 4), (2, 4)

    def holds(state, witnesses):
        return any(state & witness == witness for witness in witnesses)

    def disjoint(first, second):
        return any(first & a == a and second & b == b and not a & b
                   for a in left for b in right)

    pa = sum(law[s] for s in range(8) if holds(s, left))
    pb = sum(law[s] for s in range(8) if holds(s, right))
    one_copy = sum(law[s] for s in range(8) if disjoint(s, s))
    two_copy = sum(law[s]*law[t] for s in range(8) for t in range(8)
                   if disjoint(s, t))
    assert one_copy == sum(law[s] for s in range(8) if s.bit_count() >= 2)
    strong_gap, bk_gap = two_copy-one_copy, pa*pb-one_copy
    assert strong_gap == -sp.Rational(1803573, 6093988096)
    assert bk_gap == sp.Rational(58370337, 1523497024)
    vacancy = lambda mask: sum(law[s] for s in range(8) if not s & mask)
    va, vb, vab, vac, vbc, vabc = map(vacancy, (1, 2, 3, 5, 6, 7))
    assert strong_gap == vab-va*vb+va*vbc+vb*vac-2*vabc
    green = (sp.eye(3)-(sp.ones(3)-sp.eye(3))/3).inv()
    singleton = 1/green[0, 0]
    escapes = green.extract([0, 1], [0, 1]).inv()*sp.ones(2, 1)
    pair = singleton*sum(escapes)-escapes[0]*escapes[1]
    assert pair-singleton**2 == -sp.Rational(4, 81)
    return {
        "scope": "counterexample to disjoint-witness transfer; not to BK or NA",
        "neighbors": ((1, 2, 3), (0, 2, 3), (0, 1, 3)),
        "sleep_rate": 2, "stationary_law": list(map(str, law)),
        "left_minimal_witnesses": left, "right_minimal_witnesses": right,
        "left_probability": str(pa), "right_probability": str(pb),
        "one_copy_disjoint_probability": str(one_copy),
        "two_copy_disjoint_witness_probability": str(two_copy),
        "strengthened_gap": str(strong_gap), "ordinary_BK_gap": str(bk_gap),
        "large_rate_leading_coefficient": "-4/81",
        "raw_transient_state_counts": counts,
        "raw_law_both_orders_every_addition_and_IDLA_law_checked": True}


def validate():
    """Compare the bit operation with its definition and detect a known failure."""
    n = 3
    zeros = [sum(1 << state for state in range(8) if state & mask == 0)
             for mask in range(8)]
    cylinders = [sum(1 << state for state in range(8) if state & mask == mask)
                 for mask in range(8)]
    pairs = 0
    for aa in upsets(n):
        wa = minimal_witnesses(aa, n)
        for bb in upsets(n):
            wb = minimal_witnesses(bb, n)
            direct = 0
            for first in wa:
                for second in wb:
                    if not first & second:
                        direct |= cylinders[first | second]
            assert occurrence(wa, bb, zeros) == direct
            pairs += 1
    # Markstrom, Example 2.1: NA but not BK.
    law = tuple(sp.Rational(value, 1000) for value in (200, 240, 240, 294, 6, 10, 10, 0))
    assert all(record[-1] >= 0 for record in na_gaps(law))
    left = sum(1 << state for state in range(8) if state & 5)
    right = sum(1 << state for state in range(8) if state & 6)
    joint = occurrence(minimal_witnesses(left, n), right, zeros)
    prob = lambda bits: sum(law[state] for state in range(8) if bits >> state & 1)
    gap = prob(left)*prob(right)-prob(joint)
    assert gap == -sp.Rational(1, 2500)
    assert not certificate(law)["BK_passes"]
    assert certificate((sp.Rational(1, 8),)*8)["BK_passes"]
    return {"bit_operation_matches_definition_on_all_three_site_pairs": pairs,
            "Markstrom_NA_example_BK_gap": str(gap),
            "source": "https://emis.de/ft/41869", "uniform_product_passes": True,
            "complement_defect_lemma": complement_defect_check(),
            "product_neighborhood_lemma": product_neighborhood_certificate()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--symbolic-case", choices=("conditional_base3", "endpoint_path4"))
    parser.add_argument("--balanced-tree-na", action="store_true")
    parser.add_argument("--balanced-tree-na-symbolic", action="store_true")
    parser.add_argument("--validate-only", action="store_true")
    parser.add_argument("--two-copy-obstruction", action="store_true")
    args = parser.parse_args()
    saved = json.loads(args.output.read_text()) if args.output and args.output.exists() else {}
    if args.two_copy_obstruction:
        saved["two_copy_witness_obstruction"] = two_copy_witness_obstruction()
        print(json.dumps(saved["two_copy_witness_obstruction"], indent=2), flush=True)
        if args.output:
            args.output.write_text(json.dumps(saved, indent=2)+"\n")
        return
    if args.validate_only:
        saved["validation"] = validate()
        print(json.dumps(saved["validation"], indent=2), flush=True)
        if args.output:
            args.output.write_text(json.dumps(saved, indent=2)+"\n")
        return
    if args.balanced_tree_na_symbolic:
        saved["balanced_tree_NA_symbolic"] = balanced_tree_symbolic_na_certificate()
        print(saved["balanced_tree_NA_symbolic"]["gap_certificate"], flush=True)
        if args.output:
            args.output.write_text(json.dumps(saved, indent=2)+"\n")
        return
    if args.balanced_tree_na:
        saved["balanced_tree_NA"] = balanced_tree_na_certificate()
        print("balanced six-site tree: all 58,674 NA pairs pass at rates 1 and 64", flush=True)
        if args.output:
            args.output.write_text(json.dumps(saved, indent=2)+"\n")
        return
    if args.symbolic_case:
        examples = {"conditional_base3": ((1, 3), (0, 2, 3), (1,)),
                    "endpoint_path4": ((1, 4), (0, 2), (1, 3), (2,))}
        case = symbolic_certificate(examples[args.symbolic_case])
        saved["validation"] = validate()
        saved.setdefault("symbolic_cases", {})[args.symbolic_case] = case
        print(args.symbolic_case, case["indicators"], flush=True)
        if args.output:
            args.output.write_text(json.dumps(saved, indent=2)+"\n")
        return
    neighbors = ((1, 5), (0, 2, 3), (1,), (1, 4), (3,))
    result = {"validation": validate(), "neighbors": neighbors, "cases": {}}
    if "symbolic_cases" in saved:
        result["symbolic_cases"] = saved["symbolic_cases"]
    if "balanced_tree_NA" in saved:
        result["balanced_tree_NA"] = saved["balanced_tree_NA"]
    if "balanced_tree_NA_symbolic" in saved:
        result["balanced_tree_NA_symbolic"] = saved["balanced_tree_NA_symbolic"]
    if "two_copy_witness_obstruction" in saved:
        result["two_copy_witness_obstruction"] = saved["two_copy_witness_obstruction"]
    for rate in (1, 64):
        law = stationary(graph_kernel(neighbors), (sp.Rational(1, rate),)*5)
        for indicators, values in (("sleepers", law), ("vacancies", tuple(reversed(law)))):
            key = f"lambda_{rate}_{indicators}"
            result["cases"][key] = certificate(values)
            print(key, {k: v for k, v in result["cases"][key].items() if k != "law"}, flush=True)
    if args.output:
        args.output.write_text(json.dumps(result, indent=2)+"\n")


if __name__ == "__main__":
    main()
