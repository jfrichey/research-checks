"""Exact polynomial certificate for three sites with a permanent extra.

See drain-inventory.md. This proves a parameter-uniform occupation inequality,
not the stronger arbitrary three-site drain-inventory coefficient conjecture.
All polynomial arithmetic is integer arithmetic; no symbolic dependency or
floating-point positivity test is used. Research implementation by Codex,
without independent user validation.
"""

from collections import defaultdict
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
from math import comb, factorial
from random import Random
import hashlib
import json

from verify_inventory_corners import stopped, boundary_step


class Poly:
    """Integer polynomials in p0,p1,a,b,d,e,g,h."""
    def __init__(self, terms):
        self.terms = {key: value for key, value in terms.items() if value}

    @staticmethod
    def constant(value):
        return Poly({(0,)*8: value})

    @staticmethod
    def variable(index):
        key = [0]*8
        key[index] = 1
        return Poly({tuple(key): 1})

    def __add__(self, other):
        if isinstance(other, int):
            other = Poly.constant(other)
        out = defaultdict(int, self.terms)
        for key, value in other.terms.items():
            out[key] += value
        return Poly(out)

    __radd__ = __add__

    def __neg__(self):
        return Poly({key: -value for key, value in self.terms.items()})

    def __sub__(self, other):
        return self + (-other)

    def __rsub__(self, other):
        return (-self) + other

    def __mul__(self, other):
        if isinstance(other, int):
            other = Poly.constant(other)
        out = defaultdict(int)
        for a, c in self.terms.items():
            for b, d in other.terms.items():
                out[tuple(x+y for x, y in zip(a, b))] += c*d
        return Poly(out)

    __rmul__ = __mul__

    def evaluate(self, values):
        total = F(0)
        for key, coefficient in self.terms.items():
            term = F(coefficient)
            for value, exponent in zip(values, key):
                term *= value**exponent
            total += term
        return total


def reduced_numerator():
    p0,p1,a,b,d,e,g,h = [Poly.variable(i) for i in range(8)]
    q0,q1 = 1-p0,1-p1
    s0,s1 = b*g,e*h
    u,w = a+b*h,d+e*g
    A0,A1,B0,B1 = 1-q0*s0,1-q1*s1,1-s0,1-s1
    G = A0*A1-q0*q1*u*w
    D0,E0 = 1-q0*q1*a*d,1-(1-p0*p1)*a*d
    M = A0*A1*B0*B1-(A0*A1-p0*p1)*u*w
    S = 1-q1*d*(a+b)-q0*q1*a*d*(e+p1*d)
    R = q1*(p0*d*b*A1+q0*(a*e+b)*w)
    Z = E0*D0*A1-S*G-R*D0
    numerator = M*Z+E0*D0*p1*q0*u*w*B0*(A1-q1*w)

    # Exact division by 1-p0, checked by multiplication afterwards.
    by_remaining = defaultdict(dict)
    for key, value in numerator.terms.items():
        by_remaining[key[1:]][key[0]] = value
    quotient = {}
    for rest, coefficients in by_remaining.items():
        running = 0
        degree = max(coefficients)
        for i in range(degree):
            running += coefficients.get(i, 0)
            if running:
                quotient[(i,)+rest] = running
        assert coefficients.get(degree, 0) == -running
    value = Poly(quotient)
    assert (q0*value).terms == numerator.terms
    return value


def certificate():
    poly = reduced_numerator()
    assert len(poly.terms) == 1603
    degrees = tuple(max(key[2+2*i]+key[3+2*i] for key in poly.terms)
                    for i in range(3))
    assert degrees == (4,5,4)
    assert max(key[0] for key in poly.terms) == 2
    assert max(key[1] for key in poly.terms) == 4

    @lru_cache(None)
    def row_elevation(i,j,r):
        unused = r-i-j
        out = []
        for a in range(unused+1):
            for b in range(unused-a+1):
                c = unused-a-b
                coefficient = factorial(unused)//(factorial(a)*factorial(b)*factorial(c))
                out.append(((i+a,j+b),coefficient))
        return out

    homogeneous = defaultdict(int)
    for key, value in poly.terms.items():
        i,j,a,b,d,e,g,h = key
        choices = (row_elevation(a,b,degrees[0]),row_elevation(d,e,degrees[1]),
                   row_elevation(g,h,degrees[2]))
        for (r0,c0),(r1,c1),(r2,c2) in product(*choices):
            homogeneous[(i,j)+r0+r1+r2] += value*c0*c1*c2

    coefficients = defaultdict(int)
    for key, value in homogeneous.items():
        if not value:
            continue
        i,j,*rows = key
        for k0 in range(i,3):
            for k1 in range(j,5):
                coefficients[(k0,k1)+tuple(rows)] += value*comb(2-i,k0-i)*comb(4-j,k1-j)
    nonzero = {key: value for key, value in coefficients.items() if value}
    assert len(nonzero) == 47808
    assert all(value > 0 for value in nonzero.values())
    # Independently evaluate the fully homogeneous basis at integer values.
    # The corresponding normalized row and sleep arguments need not share
    # the same denominators; this also audits every implicit sink exponent.
    reconstructed = 0
    for key, value in nonzero.items():
        k0,k1,a,b,d,e,g,h = key
        reconstructed += (value * 2**k0 * 3**(2-k0) * 3**k1 * 4**(4-k1)
                          * 2**b * 3**(4-a-b)
                          * 2**d * 3**e
                          * 3**g * 2**(4-g-h))
    normalized = (F(2,5),F(3,7),F(1,6),F(2,6),F(2,6),F(3,6),F(3,6),F(1,6))
    assert poly.evaluate(normalized) * 5**2 * 7**4 * 6**13 == reconstructed
    digest = hashlib.sha256(json.dumps(sorted(nonzero.items()),separators=(',',':')).encode()).hexdigest()
    return dict(check='three_site_permanent_extra_polynomial_certificate',
                sleep_degrees=[2,4],row_degrees=list(degrees),
                positive_integer_coefficients=len(nonzero),
                zero_coefficients=3*5*15*21*15-len(nonzero),negative_coefficients=0,
                coefficient_sha256=digest,min_positive=min(nonzero.values()),
                max_positive=max(nonzero.values()),verified=True), poly


def occupation_audit(poly):
    rng = Random(41082)
    for trial in range(24):
        rows = []
        for z in range(3):
            weights = [rng.randrange(1,5) if d != z else 0 for d in range(4)]
            total = sum(weights)
            rows.append(tuple(F(value,total) for value in weights))
        K = tuple(rows)
        p0,p1 = [rng.choice([F(1,3),F(1,2),F(2,3)]) for _ in range(2)]
        q0,q1 = 1-p0,1-p1
        a,b = K[0][1:3]
        d,e = K[1][0],K[1][2]
        g,h = K[2][:2]
        s0,s1,u,w = b*g,e*h,a+b*h,d+e*g
        A0,A1,B0,B1 = 1-q0*s0,1-q1*s1,1-s0,1-s1
        G = A0*A1-q0*q1*u*w
        D0,E0 = 1-q0*q1*a*d,1-(1-p0*p1)*a*d
        M = A0*A1*B0*B1-(A0*A1-p0*p1)*u*w
        S = 1-q1*d*(a+b)-q0*q1*a*d*(e+p1*d)
        R = q1*(p0*d*b*A1+q0*(a*e+b)*w)
        lower = p0*(S*G+R*D0)/(E0*G*D0)
        upper = p0*(A1*M+p1*q0*u*w*B0*(A1-q1*w))/(G*M)
        predicted = p0*q0*poly.evaluate((p0,p1,a,b,d,e,g,h))/(E0*D0*G*M)
        assert min(E0,D0,G,M) > 0
        assert upper-lower == predicted >= 0
        laws,_ = stopped(((1,1,0),(1,1,1)),lambda state: boundary_step(K,(p0,p1,F(1)),state),
                         lambda state: all(a<=0 for a in state))
        actual = [sum(weight for state,weight in laws[eta].items() if state[0] == -1)
                  for eta in ((1,1,0),(1,1,1))]
        assert actual == [lower,upper]
    return dict(check='general_three_site_occupation_formulas',rational_cases=24,
                all_formulas_equal_raw_absorption=True,verified=True)


if __name__ == '__main__':
    result,poly = certificate()
    print(json.dumps(result),flush=True)
    print(json.dumps(occupation_audit(poly)),flush=True)
