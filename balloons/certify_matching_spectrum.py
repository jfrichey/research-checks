#!/usr/bin/env python3
"""Certify the exact quadratic endpoint coefficient of the matching spectrum.

AI-written research, 2026-10-15. This is pi_2 for the matching spectrum,
not the one-coordinate survivor-noise probability p_infinity.
Uses the outward integer arithmetic in certify_noise_constants.py.
"""
import argparse
import json
from pathlib import Path
import platform

from certify_noise_constants import (
    SCALE, add, ceildiv, div_int, div_positive, exp_interval,
    gamma_interval, mul_positive, rational, report_interval, sub, times_int,
)

HORIZON = 32
EIN_DEGREE = 128


def ein_interval(u):
    """Alternating entire series for Ein(u), including its first omitted term."""
    assert 0 <= u[0] <= u[1] <= HORIZON*SCALE
    term, total = rational(1), rational(0)
    for j in range(1, EIN_DEGREE+1):
        term = div_int(mul_positive(term, u), j)
        value = div_int(term, j)
        total = add(total, value) if j % 2 else sub(total, value)
    assert EIN_DEGREE % 2 == 0 and EIN_DEGREE > HORIZON
    first = div_int(div_int(mul_positive(term, u), EIN_DEGREE+1),
                    EIN_DEGREE+1)
    return total[0], total[1]+first[1]


def integrand_interval(u, decay):
    if u == (0, 0):
        return rational(1)
    ein = ein_interval(u)
    inverse_exp = div_positive(rational(1), exp_interval(times_int(ein, 2)))
    quotient = div_positive(sub(rational(1), decay), u)
    return mul_positive(inverse_exp, mul_positive(quotient, quotient))


def decimal_checks():
    """Independent implementation checks, not analytic error justification."""
    from decimal import Decimal, localcontext
    checks = 0
    with localcontext() as ctx:
        ctx.prec = 100
        scale = Decimal(SCALE)
        for j in range(129):
            u = rational(j, 4)
            half_decay = div_positive(rational(1), exp_interval(div_int(u, 2)))
            decay = mul_positive(half_decay, half_decay)
            enclosure = integrand_interval(u, decay)
            if not j:
                reference = Decimal(1)
            else:
                x = Decimal(j)/4
                term, ein = Decimal(1), Decimal(0)
                for k in range(1, 257):
                    term *= x/k
                    ein += term/k if k % 2 else -term/k
                reference = (-2*ein).exp()*((1-(-x).exp())/x)**2
            assert Decimal(enclosure[0])/scale <= reference <= Decimal(enclosure[1])/scale
            checks += 1
    return dict(integrand_comparisons=checks, decimal_precision=100,
                status='implementation validation only')


def certificate(grid_power=11):
    assert grid_power >= 7
    denominator = 2**grid_power
    count = HORIZON*denominator
    decay = rational(1)
    step_decay = div_positive(rational(1), exp_interval(rational(1, denominator)))
    total = rational(0)
    for j in range(count+1):
        value = integrand_interval(rational(j, denominator), decay)
        weight = 1 if j in (0, count) else 4 if j % 2 else 2
        total = add(total, times_int(value, weight))
        if j != count:
            decay = mul_positive(decay, step_decay)
        if j and j % 8192 == 0:
            print(f'Certified matching endpoints: {j}/{count}', flush=True)
    simpson = div_int(total, 3*denominator)
    # Integrand = E exp(-u T), T=V1+V2+U1+U2, Dickman V and uniform U.
    # Cumulants 3,7/6,2/3,29/60 give E T^4=4697/30 < 157.
    error = ceildiv(HORIZON*157*SCALE, 180*denominator**4)
    gamma = gamma_interval()
    e_minus_two_gamma = div_positive(rational(1), exp_interval(times_int(gamma, 2)))
    tail_scale = div_int(e_minus_two_gamma, 3*HORIZON**3)
    # Ein(u)=gamma+log(u)+E1(u), and 0<E1(u)<=exp(-u)/u.
    e1_upper = div_int(decay, HORIZON)
    tail_factor = div_positive(rational(1), exp_interval(times_int(e1_upper, 2)))
    one_minus_decay = sub(rational(1), decay)
    tail_factor = mul_positive(tail_factor, mul_positive(one_minus_decay, one_minus_decay))
    tail_lower = mul_positive(tail_scale, tail_factor)[0]
    tail = tail_lower, tail_scale[1]
    result = simpson[0]-error+tail[0], simpson[1]+error+tail[1]
    return dict(
        status='rigorous outward integer-interval evaluation of the exact matching coefficient',
        scope='pi_2 for the integer matching spectrum; not p_infinity or the survivor spectrum',
        python_version=platform.python_version(), arithmetic='integer intervals / 10**50; no floats',
        parameters=dict(horizon=HORIZON, grid_power=grid_power, simpson_intervals=count,
                        ein_series_degree=EIN_DEGREE, fourth_derivative_bound=157),
        fourth_moment_exact=dict(numerator=4697, denominator=30),
        gamma=report_interval(gamma),
        simpson_value=report_interval(simpson),
        simpson_error=report_interval((0, error)),
        tail=report_interval(tail),
        pi_2=report_interval(result),
        implementation_checks=decimal_checks(),
    )


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--grid-power', type=int, default=11)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = certificate(args.grid_power)
    content = json.dumps(result, indent=2)+'\n'
    if args.output:
        args.output.write_text(content)
    print(content, end='')
