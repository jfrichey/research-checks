#!/usr/bin/env python3
"""Rational certificate for the early-label cutoff error in the noise bound.

AI-written research, 2026-10-10. This implements the analytic bound in
balloon-process.tex, Proposition `prop:certified-cutoff`. It certifies an
approximation error, not the value of the limiting noise probability.
All inequalities and outward rounding use fractions.Fraction.
"""
from fractions import Fraction as F
from math import ceil, factorial
from pathlib import Path
import argparse
import json
import platform


def exponential_lower(x, degree):
    """Positive Taylor partial sum <= exp(x), for rational x >= 0."""
    assert x >= 0 and degree >= 0
    term = total = F(1)
    for k in range(1, degree + 1):
        term *= x / k
        total += term
    return total


def exponential_upper(x, degree):
    """Bound the omitted positive Taylor terms by a geometric series."""
    assert 0 <= x < degree + 2
    lower = exponential_lower(x, degree)
    first_omitted = x**(degree + 1) / factorial(degree + 1)
    return lower + first_omitted / (1 - x / (degree + 2))


def negative_exponential_upper(x):
    assert x >= 0
    degree = ceil(2 * x) + 40
    return 1 / exponential_lower(x, degree), degree


def decimal_upper(x, places=12):
    """Decimal string rounded towards +infinity, without floating point."""
    unit = 10**places
    numerator = ceil(x * unit)
    whole, remainder = divmod(numerator, unit)
    return f'{whole}.{remainder:0{places}d}'


def gamma_certificate():
    # gamma < H_n - log(n), since H_n-log(n) decreases to gamma.
    # log(2) = 2 sum_{j>=0} (1/3)^(2j+1)/(2j+1); every term is positive.
    count = 1024
    log_terms = 40
    harmonic = sum((F(1, k) for k in range(1, count + 1)), F(0))
    log_two_lower = 2 * sum((F(1, (2*j+1) * 3**(2*j+1))
                             for j in range(log_terms)), F(0))
    gamma_upper = harmonic - 10 * log_two_lower
    rounded_gamma_upper = F(ceil(gamma_upper * 10**9), 10**9)
    exp_upper = exponential_upper(rounded_gamma_upper, 40)
    convenient_upper = F(891, 500)
    assert exp_upper < convenient_upper
    return convenient_upper, dict(
        harmonic_count=count, log_two_positive_series_terms=log_terms,
        gamma_upper=str(rounded_gamma_upper), exponential_upper_degree=40,
        exponential_upper_decimal=decimal_upper(exp_upper),
        euler_constant_exponential_upper=str(convenient_upper))


def cutoff_certificate(reciprocal_cutoff, depth, exp_gamma_upper, target):
    cutoff = F(1, reciprocal_cutoff)
    v = F(reciprocal_cutoff - depth)
    assert depth >= 1 and v >= 8
    second_moment_upper = 4 * exp_gamma_upper * (1 + exp_gamma_upper)
    main = cutoff * (2 * (1 + exp_gamma_upper)
                     + F(3, 4) * second_moment_upper / factorial(depth))
    empty_exp, empty_degree = negative_exponential_upper(F(reciprocal_cutoff-1, 4))
    moment_exp, moment_degree = negative_exponential_upper(v / 4)
    empty_tail = 7 * empty_exp
    moment_tail = cutoff * 7 * v**2 * moment_exp
    total = main + empty_tail + moment_tail
    assert total < target
    return dict(
        cutoff=str(cutoff), factorial_index=depth,
        second_moment_upper=str(second_moment_upper), exponential_moment_upper=7,
        empty_tail_taylor_degree=empty_degree,
        second_moment_tail_taylor_degree=moment_degree,
        main_upper_decimal=decimal_upper(main),
        empty_tail_upper_decimal=decimal_upper(empty_tail),
        second_moment_tail_upper_decimal=decimal_upper(moment_tail),
        total_upper_decimal=decimal_upper(total),
        certified_strict_upper=str(target))


def certificate():
    exp_gamma_upper, gamma = gamma_certificate()
    inputs = ((100, F(55641, 10**6)), (200, F(27821, 10**6)),
              (1000, F(5565, 10**6)))
    return dict(
        status='rigorous rational evaluation of the proved analytic cutoff bound',
        scope='early-label agreement error; not an enclosure of p_infinity',
        python_version=platform.python_version(),
        arithmetic='fractions.Fraction; all decimal output rounded upwards',
        gamma_certificate=gamma,
        cutoff_certificates=[cutoff_certificate(r, 10, exp_gamma_upper, target)
                             for r, target in inputs])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    report = json.dumps(certificate(), indent=2) + '\n'
    if args.output:
        args.output.write_text(report)
    print(report, end='')
