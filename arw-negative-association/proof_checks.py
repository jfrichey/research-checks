#!/usr/bin/env python3
"""Symbolic checks supporting report.md; no numerical sign decisions.

Purpose: verify the Green-minor expansion on the three-site path and test
whether its small-sleep NA result extends to every common sleep rate.
Coefficientwise sign certificates are sufficient, not necessary.
Run with python3 -B research/arw/negative-association/proof_checks.py.
"""

import json
from exact import (
    DomainMatrix, absorption, graph_kernel, na_check, na_gaps, polynomial,
    rayleigh_differences, sp, stable_mask, transitions,
)


def symbolic_law(kernel, sleep):
    """Same finite chain as exact.py, over a rational-function field.

    Kernel/sleep parameters must be interpreted in their probabilistic
    ranges by the caller. Check equations and normalization as field
    identities; do not attempt symbolic positivity comparisons.
    """
    n = len(kernel)
    states = [(1,) * n]
    index = {states[0]: 0}
    rows = []
    for state in states:
        row = transitions(state, kernel, sleep)
        rows.append(row)
        for nxt in row:
            if any(v > 0 for v in nxt) and nxt not in index:
                index[nxt] = len(states)
                states.append(nxt)
    left = sp.eye(len(states))
    right = sp.zeros(len(states), 1 << n)
    for i, row in enumerate(rows):
        for nxt, prob in row.items():
            if nxt in index:
                left[i, index[nxt]] -= prob
            else:
                right[i, stable_mask(nxt)] += prob
    a = DomainMatrix.from_Matrix(left).to_field()
    b = DomainMatrix.from_Matrix(right).convert_to(a.domain)
    solution = a.inv().matmul(b)
    assert all(v == 0 for v in (a.matmul(solution) - b).to_Matrix())
    ones = DomainMatrix.from_Matrix(sp.ones(1 << n, 1)).convert_to(a.domain)
    assert all(v == 1 for v in solution.matmul(ones).to_Matrix())
    return tuple(sp.factor(v) for v in solution.to_Matrix()[0, :])


def positive_rational(expr, variable):
    """A checkable sufficient certificate of nonnegativity for variable>0."""
    num, den = sp.fraction(sp.cancel(expr))
    n, d = sp.Poly(num, variable), sp.Poly(den, variable)
    if d.LC() < 0:
        n, d = -n, -d
    return {
        "certified": all(c >= 0 for c in n.all_coeffs())
                     and all(c >= 0 for c in d.all_coeffs()),
        "numerator": str(sp.factor(n.as_expr())),
        "denominator": str(sp.factor(d.as_expr())),
    }


def two_site_check():
    p, r, a, b = sp.symbols("p r a b")
    kernel = (((1, a), (2, 1-a)), ((0, b), (2, 1-b)))
    mu = symbolic_law(kernel, (p, r))
    q, s = 1-p, 1-r
    d, h = 1-a*b+a*b*p*r, 1-a*b*q*s
    m11 = p*r/d
    m10 = p*s*(1-a*b+a*b*r*(1-b)*q)/(d*h)
    m01 = r*q*(1-a*b+a*b*p*(1-a)*s)/(d*h)
    claimed = (1-m11-m10-m01, m10, m01, m11)
    assert all(sp.cancel(x-y) == 0 for x, y in zip(mu, claimed))
    det = a*b*p*r*q*s*(1-a+a*(1-b)*q)*(1-b+b*(1-a)*s)/(d*d*h*h)
    assert sp.cancel(mu[1]*mu[2]-mu[0]*mu[3]-det) == 0
    return "all-parameter law and nonnegative determinant verified"


def preservation_check():
    kernel = graph_kernel(((1, 2), (0, 2)))
    laws, _ = absorption([(1, -1), (0, 2)], kernel, (sp.Rational(1, 2),) * 2)
    point = tuple(laws.row(0))
    uniform = tuple((laws[0, j]+laws[1, j])/2 for j in range(4))
    assert point == tuple(sp.Rational(v, 390) for v in (25, 19, 136, 210))
    assert uniform == tuple(sp.Rational(v, 156) for v in (25, 19, 58, 54))
    assert na_check(point)[1] == sp.Rational(1333, 76050)
    assert na_check(uniform)[1] == sp.Rational(31, 3042)
    return "point and uniform-drive positive covariances verified"


def three_site_check(neighbors):
    kernel = graph_kernel(neighbors)
    p, rate = sp.symbols("p rate")
    mu = symbolic_law(kernel, (p,) * 3)
    laplacian = sp.eye(3)
    for i, row in enumerate(kernel):
        for j, prob in row:
            if j < 3:
                laplacian[i, j] -= prob
    green = laplacian.inv()
    coefficients = {}
    for mask in range(1, 8):
        subset = [i for i in range(3) if mask >> i & 1]
        joint = sum(mu[s] for s in range(8) if s & mask == mask)
        lead = sp.cancel(joint / p**len(subset)).subs(p, 0)
        minor = green.extract(subset, subset).det()
        assert lead == minor
        coefficients[str(subset)] = str(lead)
    rate_mu = tuple(sp.factor(v.subs(p, rate/(1+rate))) for v in mu)
    gaps = [positive_rational(gap, rate) for _, gap in na_gaps(rate_mu)]
    stability = {}
    z, _ = polynomial(rate_mu)
    for pair, difference in rayleigh_differences(rate_mu).items():
        other = next(i for i in range(3) if i not in pair)
        quad = sp.Poly(difference, z[other])
        aa, bb, cc = (quad.nth(j) for j in (2, 1, 0))
        stability[str(pair)] = {
            "quadratic": str(sp.factor(difference)),
            "a": positive_rational(aa, rate),
            "c": positive_rational(cc, rate),
            "4ac-bb": positive_rational(4*aa*cc-bb*bb, rate),
        }
    return {
        "law_in_p": [str(v) for v in mu],
        "leading_Green_minors": coefficients,
        "all_rates_NA_certified": all(g["certified"] for g in gaps),
        "NA_gap_certificates": gaps,
        "strong_Rayleigh_certificates": stability,
    }


def main():
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", help="Optional JSON certificate file")
    args = parser.parse_args()
    result = {"two_site": two_site_check(), "preservation": preservation_check()}
    for name, neighbors in {
        "path3": ((1, 3), (0, 2), (1, 3)),
        "triangle3": ((1, 2, 3), (0, 2, 3), (0, 1, 3)),
        "path3_sink_endpoint": ((1, 3), (0, 2), (1,)),
        "path3_sink_center": ((1,), (0, 2, 3), (1,)),
        "path3_sink_endpoint_center": ((1, 3), (0, 2, 3), (1,)),
        "path3_sink_all": ((1, 3), (0, 2, 3), (1, 3)),
        "triangle3_sink_one": ((1, 2, 3), (0, 2), (0, 1)),
        "triangle3_sink_two": ((1, 2, 3), (0, 2, 3), (0, 1)),
    }.items():
        result[name] = three_site_check(neighbors)
        assert result[name]["all_rates_NA_certified"]
        strong = all(c["certified"] for pair in result[name]["strong_Rayleigh_certificates"].values()
                     for key, c in pair.items() if key != "quadratic")
        result[name]["all_rates_strong_Rayleigh_certified"] = strong
        print(name, "all-rate NA certified", result[name]["all_rates_NA_certified"],
              "all-rate strong Rayleigh certified", strong, flush=True)
    result["stationary_counterexamples"] = stationary_counterexamples(result)
    result["targeted_NA_diagnostics"] = targeted_na_checks()
    if args.output:
        from pathlib import Path
        Path(args.output).write_text(json.dumps(result, indent=2) + "\n")
    print(result["two_site"], flush=True)
    print(result["preservation"], flush=True)


def stationary_counterexamples(cases):
    """Separate rational witnesses from failure of a sufficient sign test."""
    p, rate = sp.symbols("p rate")
    mu = tuple(sp.sympify(v).subs(p, sp.Rational(64, 65))
               for v in cases["path3_sink_endpoint_center"]["law_in_p"])
    from exact import stationary_law
    direct, _ = stationary_law(graph_kernel(((1, 3), (0, 2, 3), (1,))),
                              (sp.Rational(64, 65),) * 3)
    assert mu == direct
    d = rayleigh_differences(mu)[1, 2].subs(sp.Symbol("z0"), sp.Rational(1, 100))
    assert d < 0 and na_check(mu)[1] < 0
    rayleigh = {
        "neighbors": [[1, 3], [0, 2, 3], [1]],
        "lambda": 64,
        "law": [str(v) for v in mu],
        "pair": [1, 2],
        "external_field": ["1/100", "1", "1"],
        "Rayleigh_difference": str(d),
    }
    rate_mu = tuple(sp.factor(sp.sympify(v).subs(p, rate/(1+rate)))
                    for v in cases["path3_sink_endpoint_center"]["law_in_p"])
    prob0 = sum(rate_mu[s] for s in (0, 2, 4, 6))
    delta0 = rate_mu[2]*rate_mu[4]-rate_mu[0]*rate_mu[6]
    assert sp.limit(rate*prob0, rate, sp.oo) == sp.Rational(3, 4)
    assert sp.limit(rate**4*delta0, rate, sp.oo) == -sp.Rational(1, 540)
    assert sp.limit(-rate**2*delta0/prob0**2, rate, sp.oo) == sp.Rational(4, 1215)
    mu = tuple(sp.sympify(v).subs(p, sp.Rational(3, 4))
               for v in cases["path3_sink_endpoint"]["law_in_p"])
    direct, _ = stationary_law(graph_kernel(((1, 3), (0, 2), (1,))),
                              (sp.Rational(3, 4),) * 3)
    assert mu == direct
    d = rayleigh_differences(mu)[0, 1].subs(sp.Symbol("z2"), -sp.Rational(1, 12))
    assert d == -sp.Rational(68136427, 519532701427600)
    return {
        "Rayleigh_failure": rayleigh,
        "strong_Rayleigh_failure": {
            "neighbors": [[1, 3], [0, 2], [1]],
            "lambda": 3,
            "law": [str(v) for v in mu],
            "pair": [0, 1],
            "remaining_field": "-1/12",
            "Rayleigh_difference": str(d),
        },
        "large_lambda_conditioned_covariance_coefficient": "4/1215",
    }


def targeted_na_checks():
    from exact import stationary_law
    results = []
    for attach in (0, 1, 2):
        neighbors = [[1, 4], [0, 2, 4], [1], []]
        neighbors[attach].append(3)
        neighbors[3].append(attach)
        mu, _ = stationary_law(graph_kernel(neighbors), (sp.Rational(64, 65),)*4)
        count, cov, witness = na_check(mu)
        assert count == 174 and cov < 0
        results.append({"neighbors": neighbors, "rates": ["64"]*4,
                        "tests": count, "largest_covariance": str(cov), "witness": witness})
    for rate0 in (sp.Rational(1, 64), 1, 64, 4096):
        rates = tuple(map(sp.Rational, (rate0, 64, 64)))
        neighbors = ((1, 3), (0, 2, 3), (1,))
        mu, _ = stationary_law(graph_kernel(neighbors), tuple(r/(1+r) for r in rates))
        count, cov, witness = na_check(mu)
        assert count == 15 and cov < 0
        results.append({"neighbors": neighbors, "rates": list(map(str, rates)),
                        "tests": count, "largest_covariance": str(cov), "witness": witness})
    return results


if __name__ == "__main__":
    main()
