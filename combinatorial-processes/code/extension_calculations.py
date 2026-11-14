"""Deterministic calculations for the November 3 manuscript extensions."""
from __future__ import annotations

import numpy as np


def adsorption_distributions(max_n: int) -> list[np.ndarray]:
    """Evaluate the exact PGF splitting recurrence in double precision.

    C_-1 = C_0 = 1 and C_n = z/n sum_{i=-1}^{n-2} C_i C_{n-3-i}.
    The returned list is indexed by n >= 0; no simulation is used.
    """
    pmfs = [np.array([1.0])]
    empty = pmfs[0]
    for n in range(1, max_n + 1):
        out = np.zeros((n+1)//2 + 1)
        for i in range(-1, n-1):
            j = n-3-i
            if i > j:
                break
            left, right = pmfs[i] if i >= 0 else empty, pmfs[j] if j >= 0 else empty
            conv = np.convolve(left, right)
            out[1:len(conv)+1] += (1 if i == j else 2)*conv/n
        pmfs.append(out)
    return pmfs


def standardized_pmf(pmf: np.ndarray):
    from scipy.special import ndtr
    x = np.arange(len(pmf), dtype=np.longdouble)
    p = pmf.astype(np.longdouble)
    mean = np.sum(p*x)
    variance = np.sum(p*(x-mean)**2)
    third = np.sum(p*(x-mean)**3)
    z = np.asarray((x-mean)/np.sqrt(variance), dtype=float)
    cdf = np.cumsum(pmf)
    normal = ndtr(z)
    distance = max(np.max(np.abs(cdf-normal)), np.max(np.abs(cdf-pmf-normal)))
    return z, cdf, float(mean), float(variance), float(third), float(distance)


def four_exception_sieve(limit: int) -> np.ndarray:
    """Boolean n: odd, with every prime divisor congruent to 1 modulo 4.

    Includes 1; callers studying 4/n exclude n <= 4.
    """
    is_prime = np.ones(limit+1, dtype=bool)
    is_prime[:2] = False
    for p in range(2, int(limit**.5)+1):
        if is_prime[p]:
            is_prime[p*p::p] = False
    exceptional = np.ones(limit+1, dtype=bool)
    exceptional[::2] = False
    for p in np.flatnonzero(is_prime):
        if p % 4 == 3:
            exceptional[p::p] = False
    return exceptional
