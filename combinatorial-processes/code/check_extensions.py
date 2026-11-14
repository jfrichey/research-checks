"""Substantive finite and analytic checks for the November 3 extensions.

These checks audit constants, enumerations, search certificates, and source/data
consistency. They do not constitute independent proofs of the limit theorems.
"""
from collections import Counter
from fractions import Fraction as F
from itertools import permutations
from math import factorial, gcd, exp
import hashlib
import json
from pathlib import Path

import mpmath as mp
import numpy as np
import sympy as sp
from scipy.integrate import quad

from egyptian_fractions import greedy, binary_remainder, shortest, four_greedy_length, four_needs_three
from extension_calculations import adsorption_distributions, standardized_pmf, four_exception_sieve

ROOT = Path(__file__).resolve().parents[1]


def matchings(items):
    if not items:
        yield ()
        return
    a = items[0]
    for i in range(1, len(items)):
        b = items[i]
        for rest in matchings(items[1:i]+items[i+1:]):
            yield ((a, b),)+rest


def sock_paths():
    for n in range(1, 7):
        counts, maxima = Counter(), Counter()
        for matching in matchings(tuple(range(2*n))):
            openings = {a for a, b in matching}
            path = [0]
            for j in range(2*n):
                path.append(path[-1]+(1 if j in openings else -1))
            counts[tuple(path)] += 1
            maxima[max(path)] += 1
        for path, count in counts.items():
            weight = 1
            for a, b in zip(path, path[1:]):
                if b < a:
                    weight *= a
            assert weight == count
        assert sum(counts.values()) == factorial(2*n)//(2**n*factorial(n))
        for h in range(n+1):
            matrix = np.zeros((h+1, h+1), dtype=object)
            for j in range(h):
                matrix[j, j+1] = 1
                matrix[j+1, j] = j+1
            assert np.linalg.matrix_power(matrix, 2*n)[0, 0] == sum(v for k, v in maxima.items() if k <= h)
    print('PASS: weighted Dyck paths and exact maximum CDF versus all matchings through six pairs')


def adsorption():
    pmfs = adsorption_distributions(640)
    for n in range(1, 8):
        counts = Counter()
        for order in permutations(range(n)):
            occupied = set()
            for j in order:
                if all(abs(j-k)>1 for k in occupied):
                    occupied.add(j)
            counts[len(occupied)] += 1
        for k, p in enumerate(pmfs[n]):
            assert abs(p-counts[k]/factorial(n)) < 1e-14
    s = sp.symbols('s')
    c1, c2, c3 = sp.symbols('c1 c2 c3')
    y = 1+c1*s+c2*s**2+c3*s**3
    eq = sp.series((y-1)*sp.exp(2*y)-sp.tanh(s/4)*(y+1), s, 0, 4).removeO()
    sol = {}
    for k, unknown in enumerate((c1,c2,c3), 1):
        sol[unknown] = sp.solve(eq.coeff(s, k).subs(sol), unknown)[0]
    lam = sp.series(s/2-sp.log(y.subs(sol)), s, 0, 4).removeO()
    expected = [(1-sp.exp(-2))/2, sp.exp(-4), (sp.exp(-2)-61*sp.exp(-6))/16]
    for k, value in enumerate(expected, 1):
        assert sp.simplify(lam.coeff(s,k)*factorial(k)-value) == 0
    a = standardized_pmf(pmfs[160]); b = standardized_pmf(pmfs[640])
    for idx, value in zip((2,3,4), expected):
        assert abs((b[idx]-a[idx])/480-float(value)) < 1e-10
    print('PASS: full adsorption PMFs by permutation enumeration; pole derivatives and cumulant slopes')


def uniform_renewal():
    a, y = sp.symbols('a y', real=True)
    low, high = sp.exp(a), sp.exp(a)-sp.exp(a+y-1)
    mass = sp.integrate(low, (y,0,1-a))+sp.integrate(high,(y,1-a,1))
    invariant = sp.integrate(2*a*low, (a,0,1-y))+sp.integrate(2*a*high, (a,1-y,1))
    assert sp.simplify(mass-1) == 0
    assert sp.simplify(invariant-2*(1-y)) == 0
    for av in np.linspace(.01,1,80):
        for yv in np.linspace(.01,.99,80):
            density = exp(av) if av+yv <= 1 else exp(av)-exp(av+yv-1)
            assert density >= 1-yv-1e-14
    mp.mp.dps = 80
    for x in [mp.mpf(j)/4 for j in range(1,161)]:
        m = mp.fsum((-1)**j/mp.factorial(j)*(x-j)**j*mp.exp(x-j) for j in range(int(x)+1))
        assert abs(m-2*x-mp.mpf(2)/3) <= mp.power(2,1-int(x))
    print('PASS: residual kernel normalization, invariant density, minorization, and exponential mean bound')


def rapid_records():
    mp.mp.dps = 180
    masses = [mp.exp(-1)/mp.factorial(k) for k in range(170)]
    tails = [mp.fsum(masses[k:]) for k in range(170)]
    def unique(n):
        return mp.fsum(n*masses[k]*mp.exp((n-1)*mp.log1p(-tails[k]))
                       for k in range(1,165))  # k=0 has F_-1=0
    target_errors, troughs = [], []
    for j in (12,30,60):
        n = int(1/masses[j])
        target_errors.append(abs(unique(n)-1/mp.e))
        n = int(1/mp.sqrt(masses[j]*masses[j+1]))
        troughs.append(unique(n))
    assert target_errors[0] > target_errors[1] > target_errors[2]
    assert target_errors[-1] < mp.mpf('.02')
    assert troughs[0] > troughs[1] > troughs[2]
    print('PASS: high-precision Poisson unique-maximum subsequences and their distinct limits')


def brute_exists(r, slots, lower=2):
    """Small-case oracle using Fraction arithmetic and direct denominator loops."""
    if r.numerator == 1 and r.denominator >= lower:
        return True
    if slots == 1:
        return False
    start = max(lower, (r.denominator+r.numerator-1)//r.numerator)
    stop = slots*r.denominator//r.numerator
    for d in range(start, stop+1):
        tail = r-F(1,d)
        if tail > 0 and brute_exists(tail, slots-1, d+1):
            return True
    return False


def egyptian_checks():
    for n in range(2,101):
        for k in range(1,n):
            for ds in (greedy(k,n), binary_remainder(k,n)):
                assert sum((F(1,d) for d in ds), F()) == F(k,n)
                assert list(ds) == sorted(set(ds))
            assert len(greedy(k,n)) <= k//gcd(k,n)
            assert len(binary_remainder(k,n)) <= 2*(n-1).bit_length()
            if n <= 18:
                best = shortest(k,n)
                assert sum((F(1,d) for d in best), F()) == F(k,n)
                assert all(not brute_exists(F(k,n), h) for h in range(1,len(best)))
    exception = four_exception_sieve(10_000)
    for n in range(5,10_001):
        assert len(greedy(4,n)) == four_greedy_length(n)
        assert bool(exception[n]) == four_needs_three(n)
        if n <= 500:
            assert (len(shortest(4,n)) >= 3) == bool(exception[n])
    assert Counter(four_greedy_length(n) for n in range(5,29)) == {1:6,2:14,3:2,4:2}
    for a in range(1,8):
        period = factorial(a)
        first = [len(greedy(a,n)) for n in range(a+1,a+period+1)]
        for shift in (period, 2*period):
            assert first == [len(greedy(a,n+shift)) for n in range(a+1,a+period+1)]
    d, b, r = sp.symbols('d b r', positive=True)
    ratio = sp.integrate(d*r-1,(r,1/d,b))/sp.integrate(r,(r,1/d,b))
    assert sp.simplify(ratio-(d-2/(b+1/d))) == 0
    assert sp.simplify(ratio.subs(b,1/(d-1))-d/(2*d-1)) == 0
    data = json.loads((ROOT/'computations'/'egyptian-lengths.json').read_text())
    for n, row in data['optimal'].items():
        n = int(n)
        assert len(row['representations']) == n-1
        for k, ds in enumerate(row['representations'],1):
            assert sum((F(1,d) for d in ds),F()) == F(k,n)
            assert ds == sorted(set(ds))
    for n, row in data['greedy'].items():
        assert sum(row['counts'].values()) == int(n)-1 == row['completed']
        assert row['censored'] == 0
    assert abs(data['optimal']['127']['mean']-487/126) < 1e-14
    assert abs(data['optimal']['128']['mean']-386/127) < 1e-14
    print('PASS: exact Egyptian constructions, minimum-search oracle, residue periods, and all saved certificates')


def manifests():
    data = json.loads((ROOT/'figures'/'extensions-manifest.json').read_text())
    for name, row in data['figures'].items():
        assert hashlib.sha256((ROOT/'figures'/(name+'.pdf')).read_bytes()).hexdigest() == row['sha256']
    for folder, key in [('computations','data_sha256'),('code','code_sha256')]:
        for name, digest in data[key].items():
            assert hashlib.sha256((ROOT/folder/name).read_bytes()).hexdigest() == digest
    print('PASS: extension figure, data, and generating-code hashes')


if __name__ == '__main__':
    for check in (sock_paths, adsorption, uniform_renewal, rapid_records, egyptian_checks, manifests):
        check()
