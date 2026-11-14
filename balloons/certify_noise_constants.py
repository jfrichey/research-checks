#!/usr/bin/env python3
"""Outward integer-interval certificate for the exact noise bound constants.

AI-written research, 2026-10-15. No binary floating-point operations are used.
The certificate evaluates c_star, c_dagger, c_ddagger, and u_star;
it does not identify p_infinity or certify the historical candidate near .658.
"""
import argparse
import json
from pathlib import Path
import platform

SCALE = 10**50
LOG_TERMS = 60
EXP_DEGREE = 128
E1_DEGREE = 128
GAMMA_COUNT = 2**18
HORIZON = 16


def ceildiv(a, b):
    assert b > 0
    return -((-a)//b)


def rational(a, b=1):
    return a*SCALE//b, ceildiv(a*SCALE, b)


def add(x, y):
    return x[0]+y[0], x[1]+y[1]


def sub(x, y):
    return x[0]-y[1], x[1]-y[0]


def mul_positive(x, y):
    assert x[0] >= 0 and y[0] >= 0
    return x[0]*y[0]//SCALE, ceildiv(x[1]*y[1], SCALE)


def div_positive(x, y):
    assert x[0] >= 0 and y[0] > 0
    return x[0]*SCALE//y[1], ceildiv(x[1]*SCALE, y[0])


def div_int(x, n):
    return x[0]//n, ceildiv(x[1], n)


def times_int(x, n):
    return (n*x[0], n*x[1]) if n >= 0 else (n*x[1], n*x[0])


def log_reduced(num, den):
    """log(num/den), with 1 <= num/den <= 2, via positive atanh series."""
    assert den <= num <= 2*den
    u = rational(num-den, num+den)
    u2 = mul_positive(u, u)
    power, total = u, (0, 0)
    for j in range(LOG_TERMS):
        total = add(total, div_int(power, 2*j+1))
        power = mul_positive(power, u2)
    # Remaining powers have denominators >= 2*LOG_TERMS+1.
    tail = div_positive(div_int(power, 2*LOG_TERMS+1),
                        sub((SCALE, SCALE), u2))
    return 2*total[0], 2*(total[1]+tail[1])


LOG_TWO = log_reduced(2, 1)


def gamma_interval():
    lower = sum(SCALE//j for j in range(1, GAMMA_COUNT+1))
    upper = sum(ceildiv(SCALE, j) for j in range(1, GAMMA_COUNT+1))
    center = sub((lower, upper), times_int(LOG_TWO, 18))
    # 1/(2(n+1)) < H_n-log(n)-gamma < 1/(2n).
    return (center[0]-rational(1, 2*GAMMA_COUNT)[1],
            center[1]-rational(1, 2*(GAMMA_COUNT+1))[0])


def exp_interval(x):
    assert 0 <= x[0] <= x[1] <= 16*SCALE
    term = total = (SCALE, SCALE)
    for j in range(1, EXP_DEGREE+1):
        term = div_int(mul_positive(term, x), j)
        total = add(total, term)
    first = div_int(mul_positive(term, x), EXP_DEGREE+1)
    denominator = SCALE-ceildiv(x[1], EXP_DEGREE+2)
    tail_upper = ceildiv(first[1]*SCALE, denominator)
    return total[0], total[1]+tail_upper


def laplace_at_grid(j, grid_power, gamma):
    """F(s)=tanh(E1(s)/2) for s=j/2**grid_power; F(0)=1."""
    if j == 0:
        return SCALE, SCALE
    binary_power = j.bit_length()-1
    log_s = add(log_reduced(j, 2**binary_power),
                times_int(LOG_TWO, binary_power-grid_power))
    s = rational(j, 2**grid_power)
    term, series = (SCALE, SCALE), (0, 0)
    for n in range(1, E1_DEGREE+1):
        term = div_int(mul_positive(term, s), n)
        contribution = div_int(term, n)
        series = add(series, contribution) if n % 2 == 0 else sub(series, contribution)
    # The alternating tail decreases in magnitude after this truncation.
    assert E1_DEGREE+1 > HORIZON
    first = div_int(div_int(mul_positive(term, s), E1_DEGREE+1),
                    E1_DEGREE+1)
    if E1_DEGREE % 2 == 0:
        series = (series[0]-first[1], series[1])
    else:
        series = (series[0], series[1]+first[1])
    e1 = sub(sub((0, 0), gamma), add(log_s, series))
    assert e1[1] > 0
    e1 = (max(0, e1[0]), e1[1])
    exponential = exp_interval(e1)
    return ((exponential[0]-SCALE)*SCALE//(exponential[0]+SCALE),
            ceildiv((exponential[1]-SCALE)*SCALE, exponential[1]+SCALE))


def decimal(x, places, upper=False):
    scaled = ceildiv(x*10**places, SCALE) if upper else x*10**places//SCALE
    sign = '-' if scaled < 0 else ''
    whole, rest = divmod(abs(scaled), 10**places)
    return f'{sign}{whole}.{rest:0{places}d}'


def report_interval(x):
    assert x[0] <= x[1]
    return dict(lower_numerator=str(x[0]), upper_numerator=str(x[1]),
                denominator=str(SCALE),
                outward_decimal_12=[decimal(x[0], 12), decimal(x[1], 12, True)],
                outward_decimal_9=[decimal(x[0], 9), decimal(x[1], 9, True)])


def decimal_checks():
    """Independent standard-library checks; not an input to the certificate."""
    from decimal import Decimal, localcontext
    checks = 0
    with localcontext() as context:
        context.prec = 100
        scale = Decimal(SCALE)
        for j in range(257):
            interval = exp_interval(rational(j, 16))
            reference = (Decimal(j)/16).exp()
            assert Decimal(interval[0])/scale <= reference <= Decimal(interval[1])/scale
            checks += 1
        for den in (1, 3, 7, 31, 127):
            for num in range(den, 2*den+1):
                interval = log_reduced(num, den)
                reference = (Decimal(num)/den).ln()
                assert Decimal(interval[0])/scale <= reference <= Decimal(interval[1])/scale
                checks += 1
    return dict(exp_log_comparisons=checks, decimal_precision=100,
                status='implementation validation only; certificate uses integer bounds')


def certificate(grid_power=12):
    assert grid_power >= 8
    denominator = 2**grid_power
    count = HORIZON*denominator
    assert count % 2 == 0
    gamma = gamma_interval()
    step_exp = exp_interval(rational(1, denominator))
    step_decay = div_positive((SCALE, SCALE), step_exp)
    decay = (SCALE, SCALE)
    sum_g = sum_h = (0, 0)
    for j in range(count+1):
        f = laplace_at_grid(j, grid_power, gamma)
        decay2 = mul_positive(decay, decay)
        g = mul_positive(decay2, mul_positive(f, f))
        h = mul_positive(rational(j, denominator), mul_positive(decay, f))
        weight = 1 if j in (0, count) else 4 if j % 2 else 2
        sum_g = add(sum_g, times_int(g, weight))
        sum_h = add(sum_h, times_int(h, weight))
        if j != count:
            decay = mul_positive(decay, step_decay)
        if j and j % 8192 == 0:
            print(f'Certified endpoint intervals: {j}/{count}', flush=True)

    integral_g = div_int(sum_g, 3*denominator)
    integral_h = div_int(sum_h, 3*denominator)
    # G(s)=E exp(-s(2+Y1+Y2)); |G''''| <= 82*4!*4**4=503808.
    # H(s)=E s exp(-s(1+Y)); |H''''| <= 5*E(1+Y)**3 < 17920.
    # These use E exp(Y/4)<7, exp(1/2)<5/3, exp(1/4)<4/3.
    assert exp_interval(rational(1, 2))[1] < rational(5, 3)[0]
    assert exp_interval(rational(1, 4))[1] < rational(4, 3)[0]
    error_g = ceildiv(HORIZON*503808*SCALE, 180*denominator**4)
    error_h = ceildiv(HORIZON*17920*SCALE, 180*denominator**4)
    # Y >= 1 gives tails G <= exp(-4s), H <= s exp(-2s).
    tail_g = div_int(mul_positive(decay2, decay2), 4)[1]
    tail_h = mul_positive(rational(2*HORIZON+1, 4), decay2)[1]
    integral_g = (integral_g[0]-error_g, integral_g[1]+error_g+tail_g)
    integral_h = (integral_h[0]-error_h, integral_h[1]+error_h+tail_h)
    c_star = sub(rational(1, 4), div_int(integral_g, 2))
    log_32, log_85 = log_reduced(3, 2), log_reduced(8, 5)
    eta = div_positive(log_85, sub(rational(12), log_85))
    extra = sub(add(integral_h, div_int(log_32, 2)), rational(1, 12))
    c_dagger = add(c_star, mul_positive(eta, extra))
    c_balanced = div_int(log_85, 24)
    c_extra = sub(add(rational(1, 48), div_int(LOG_TWO, 24)),
                  div_int(log_reduced(5, 3), 12))
    eta_plus = div_positive(div_int(log_32, 6),
                           sub(sub(rational(1), times_int(c_balanced, 2)),
                               mul_positive(rational(2, 5), c_extra)))
    c_ddagger = add(c_star, mul_positive(eta_plus, extra))
    assert c_ddagger[0] > c_dagger[1]
    u_star = sub(sub(rational(23, 24), div_int(log_32, 4)),
                 div_int(integral_h, 2))
    return dict(
        status='rigorous outward integer-interval evaluation of proved exact bounds',
        scope='c_star, c_dagger, c_ddagger, u_star; does not certify the candidate p_infinity near .658',
        python_version=platform.python_version(), arithmetic='integer intervals / 10**50; no floats',
        parameters=dict(horizon=HORIZON, grid_power=grid_power, simpson_intervals=count,
                        gamma_harmonic_count=GAMMA_COUNT, log_terms=LOG_TERMS,
                        exponential_degree=EXP_DEGREE, exponential_integral_degree=E1_DEGREE),
        gamma=report_interval(gamma),
        implementation_checks=decimal_checks(),
        error_bounds=dict(simpson_g=decimal(error_g, 15, True),
                          simpson_h=decimal(error_h, 15, True),
                          tail_g=decimal(tail_g, 15, True), tail_h=decimal(tail_h, 15, True)),
        constants={name: report_interval(value) for name, value in
                   [('c_star', c_star), ('laplace_second_moment_integral', integral_h),
                    ('c_dagger', c_dagger), ('half_interval_extra', c_extra),
                    ('eta_plus', eta_plus), ('c_ddagger', c_ddagger), ('u_star', u_star)]},
        certified_p_infinity_enclosure=[decimal(c_ddagger[0], 10), decimal(u_star[1], 10, True)])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--grid-power', type=int, default=12)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = certificate(args.grid_power)
    content = json.dumps(result, indent=2)+'\n'
    if args.output:
        args.output.write_text(content)
    print(content, end='')
