"""Exact Egyptian-fraction constructions and finite, certified length searches.

All fractions are reduced with integer arithmetic. ``shortest`` searches by
length, with the first denominator bounded by floor(k / remainder); it never
uses an arbitrary denominator cutoff. Its cost can grow rapidly.
"""
from __future__ import annotations

from functools import lru_cache
from math import gcd


def reduced(a: int, b: int) -> tuple[int, int]:
    g = gcd(a, b)
    return a // g, b // g


def greedy(a: int, b: int, *, bit_limit: int | None = None) -> tuple[int, ...]:
    if not 0 < a < b:
        raise ValueError('Require 0 < a < b.')
    a, b = reduced(a, b)
    terms = []
    while a:
        d = (b + a - 1) // a
        terms.append(d)
        a, b = reduced(a * d - b, b * d)
        if bit_limit is not None and b.bit_length() > bit_limit:
            raise ArithmeticError('Exact computation exceeded the declared bit limit.')
    return tuple(terms)


def binary_remainder(a: int, b: int) -> tuple[int, ...]:
    if not 0 < a < b:
        raise ValueError('Require 0 < a < b.')
    m = (b - 1).bit_length()
    q = 1 << m
    s, r = divmod(a * q, b)
    first = [q // (1 << j) for j in range(m) if (s >> j) & 1]
    second = [b * q // (1 << j) for j in range(m) if (r >> j) & 1]
    return tuple(sorted(first + second))


@lru_cache(maxsize=100_000)
def _square_divisors(b: int) -> tuple[int, ...]:
    from sympy import factorint
    divs = [1]
    for p, e in factorint(b).items():
        divs = [d * int(p)**j for d in divs for j in range(2*e + 1)]
    return tuple(sorted(d for d in divs if d <= b))


def _two(a: int, b: int, lower: int) -> tuple[int, ...] | None:
    if a == 1 and b >= lower:
        return (b,)
    for d in _square_divisors(b):
        if (b + d) % a:
            continue
        x = (b + d) // a
        if x < lower:
            continue
        other = b*b // d
        if (b + other) % a == 0:
            y = (b + other) // a
            if x < y:
                return (x, y)
    return None


def _search(a: int, b: int, slots: int, lower: int) -> tuple[int, ...] | None:
    if a == 1 and b >= lower:
        return (b,)
    if slots == 1 or a * lower >= slots * b:
        return None
    if slots == 2:
        return _two(a, b, lower)
    for d in range(max(lower, (b + a - 1)//a), slots*b//a + 1):
        c, q = reduced(a*d-b, b*d)
        if c == 0:
            return (d,)
        tail = _search(c, q, slots-1, d+1)
        if tail is not None:
            return (d,) + tail
    return None


def shortest(a: int, b: int) -> tuple[int, ...]:
    """Return a shortest representation, proving failure of all shorter lengths."""
    if not 0 < a < b:
        raise ValueError('Require 0 < a < b.')
    a, b = reduced(a, b)
    known = min(greedy(a, b), binary_remainder(a, b), key=len)
    for k in range(1, len(known)):
        found = _search(a, b, k, 2)
        if found is not None:
            return found
    return known


def four_greedy_length(n: int) -> int:
    if n <= 4:
        raise ValueError('Require n > 4.')
    if n % 4 == 0:
        return 1
    if n % 24 in (1, 17):
        return 4
    if n % 24 in (5, 13):
        return 3
    return 2


def four_needs_three(n: int) -> bool:
    """Exact criterion for minimum length >= 3; makes no Erdos--Straus assumption."""
    if n <= 4:
        raise ValueError('Require n > 4.')
    if n % 2 == 0:
        return False
    from sympy import factorint
    return all(int(p) % 4 == 1 for p in factorint(n))
