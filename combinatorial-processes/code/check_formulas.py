#!/usr/bin/env python3
"""Independent finite-model and analytic checks for notes.tex.

Generated with ChatGPT, 2026-11-02. These checks support the manuscript audit;
they do not certify asymptotic proofs. Requires sympy, scipy, numpy, mpmath.
Run from any directory; prints results and writes no files.
"""
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from itertools import permutations, product
from math import comb, factorial, exp

import mpmath as mp
import numpy as np
import sympy as sp
from scipy.integrate import quad
from scipy.stats import irwinhall


def falling(n, r):
    if n < r:
        return 0
    return factorial(n) // factorial(n-r)


def socks():
    for n, l in [(2, 2), (3, 2), (4, 2), (2, 3), (2, 4)]:
        N = n*l
        totals = np.zeros(N+1, dtype=np.int64)
        cross = np.zeros((N+1, N+1), dtype=np.int64)
        for order in permutations(range(N)):
            counts = [0]*n
            ys = [0]
            completed = 0
            for obj in order:
                counts[obj//l] += 1
                completed += counts[obj//l] == l
                ys.append(completed)
            y = np.array(ys, dtype=np.int64)
            totals += y
            cross += np.outer(y, y)
        size = factorial(N)
        for k in range(N+1):
            mean = F(n*falling(k, l), falling(N, l))
            assert F(int(totals[k]), size) == mean
            second = mean + F(n*(n-1)*falling(k, 2*l), falling(N, 2*l))
            assert F(int(cross[k, k]), size) == second
            for m in range(k, N+1):
                mixed = F(0) if m < l else mean + (
                    F(n*(n-1)*falling(k, l), falling(N, l))
                    * F(falling(m-l, l), falling(N-l, l)))
                assert F(int(cross[k, m]), size) == mixed
        if l == 2:
            for k in range(N-2):
                mean_k = F(k*(N-k), N-1)
                mean_next = F((k+1)*(N-k-1), N-1)
                for x in range(N+1):
                    prob = F(x, N-k)
                    cond = ((1-prob)*(x+1-mean_next)
                            + prob*(x-1-mean_next)) / ((N-k-1)*(N-k-2))
                    assert cond == (x-mean_k)/((N-k)*(N-k-1))
    print('PASS: sock means, second moments, mixed moments, and pair martingale')


def permutations_and_seats():
    a = [F(0)]*9
    b = [F(0)]*9
    for n in range(1, 9):
        masks = Counter()
        seats = Counter()
        for order in permutations(range(n)):
            mask = sum(1 << j for j in range(n-1) if order[j+1] == order[j]+1)
            masks[mask] += 1
            occupied = 0
            for i in order:
                if not occupied & ((1 << i) | (1 << (i+1)) | ((1 << (i-1)) if i else 0)):
                    occupied |= 1 << i
            seats[occupied.bit_count()] += 1
        for subset in range(1 << (n-1)):
            count = sum(v for mask, v in masks.items() if mask & subset == subset)
            assert count == factorial(n-subset.bit_count())
        a[n] = 1 + F(2, n)*sum(a[1:n-1])
        b[n] = (2*a[n]-1 + F(2, n)*sum(a[i]*a[n-3-i] for i in range(1, n-3))
                + F(2, n)*sum(b[1:n-1]))
        assert F(sum(k*v for k, v in seats.items()), factorial(n)) == a[n]
        assert F(sum(k*k*v for k, v in seats.items()), factorial(n)) == b[n]
    x, z = sp.symbols('x z')
    f = (1-sp.exp(-2*x))/(2*(1-x)**2)
    g = ((1+x)*(1+sp.exp(-4*x))-2*(1-x+2*x*x)*sp.exp(-2*x))/(4*(1-x)**3)
    h = -(x*sp.sqrt(z)*sp.cosh(x*sp.sqrt(z))+(x*z+z-1)*sp.sinh(x*sp.sqrt(z)))/(
        (x-1)*sp.sqrt(z)*sp.cosh(x*sp.sqrt(z))+(x*z-1)*sp.sinh(x*sp.sqrt(z)))
    assert sp.simplify(sp.diff(g,x)-(2*sp.diff(f,x)-1/(1-x)**2+2*x*x*f*f+2*x*g/(1-x))) == 0
    assert sp.simplify(sp.diff(h,x)-z*(1+x+x*h)**2) == 0
    y = sp.symbols('y', real=True)
    hy = sp.exp(-y)-(1-y)/sp.E
    hz = sp.exp(-z)-(1-z)/sp.E
    four = sp.integrate(2*(1-y)*hy*sp.integrate(hz,(z,y,1)),(y,0,1))
    assert sp.simplify(four-(sp.Rational(1,2)-sp.Rational(79,30)*sp.exp(-2))) == 0
    print('PASS: all permutations/seatings through n=8; seating ODEs and four-point integral')


def coupon_finite():
    for n in [2, 3, 4]:
        length = 6
        total = np.zeros(length+1, dtype=np.int64)
        cross = np.zeros((length+1, length+1), dtype=np.int64)
        for draws in product(range(n), repeat=length):
            seen = set()
            x = [0]
            for v in draws:
                seen.add(v)
                x.append(len(seen))
            x = np.array(x)
            total += x
            cross += np.outer(x, x)
        size = n**length
        for k in range(length+1):
            ak = F(n-1,n)**k
            assert F(int(total[k]),size) == n*(1-ak)
            for m in range(k,length+1):
                am = F(n-1,n)**m
                cov = n*am*(1-ak)+n*(n-1)*(F(n-2,n)**k*F(n-1,n)**(m-k)-ak*am)
                actual = F(int(cross[k,m]),size)-F(int(total[k]*total[m]),size*size)
                assert cov == actual
    print('PASS: coupon covariance against every six-draw history for n=2,3,4')


def hitting_times():
    mp.mp.dps = 60
    for x in [.2, .858, 1., 1.5, 2., 2.3, 4.2, 10.5]:
        xx = mp.mpf(str(x))
        mean = sum((-1)**j*(xx-j)**j*mp.exp(xx-j)/mp.factorial(j)
                   for j in range(int(mp.floor(xx))+1))
        tails = [1.0] + [irwinhall(k).cdf(x) for k in range(1, 80)]
        mean_ind = sum(tails)
        second_ind = sum((2*k+1)*v for k, v in enumerate(tails))
        assert abs(float(mean)-mean_ind) < 1e-10
        assert abs(float((2*xx+1)*mean)-second_ind) < 1e-9
    for k in range(1,11):
        for p in [F(1,3), F(1,2), F(3,4)]:
            q = 1-p
            for s in range(1,k+1):
                # Inversion via the negative-binomial law for S_s.
                inverse = sum(F(comb(t-1,s-1))*p**s*q**(t-s) for t in range(s,k))
                binomial = sum(F(comb(k-1,r))*p**r*q**(k-1-r) for r in range(s,k))
                assert inverse == binomial
    print('PASS: uniform-sum moments from independent CDFs; geometric inversion')


@lru_cache(None)
def survivors(x, y):
    if x*y == 0:
        return {x+y:F(1)}
    result = Counter()
    for state in [(x-1,y),(x,y-1)]:
        for r,p in survivors(*state).items():
            result[r] += p/2
    return dict(result)


def arrows():
    for k in range(1,17):
        dist = Counter()
        for b in range(k+1):
            for r,p in survivors(b,k-b).items():
                dist[r] += F(comb(k,b),2**k)*p
        assert sum(dist.values()) == 1
        assert sum((r*r+r)*p for r,p in dist.items()) == 2*k
    for k in range(2,11):
        counts = Counter()
        for signs in product([-1,1],repeat=k):
            z = sum(signs[i] != signs[(i+1)%k] for i in range(k))
            counts[z] += 1
            for i in range(k):
                if signs[i] == signs[(i+1)%k]:
                    continue
                after = signs[:i]+signs[i+1:]
                zz = sum(after[j] != after[(j+1)%len(after)] for j in range(len(after)))
                assert zz-z in [0,-2]
        assert counts == Counter({r:2*comb(k,r) for r in range(0,k+1,2)})
    print('PASS: arrow absorption identity through k=16 and all cyclic signs through k=10')


def circle():
    for n in range(1,16):
        for theta in [.03,.2,.4,.5,.7,.9]:
            q = lambda d:max(0.,1-theta-min(d,theta))**n
            integral = 2*quad(q,0,.5,points=sorted(set([min(theta,.5),min(1-theta,.5)])),epsabs=1e-12)[0]
            formula = 2/(n+1)*(1-theta)**(n+1)+(n-1)/(n+1)*max(1-2*theta,0)**(n+1)
            assert abs(integral-formula)<1e-11
            if n==1:
                assert abs(formula-(1-theta)**2)<1e-12
    for theta in [F(1,3),F(1,2),F(3,5),F(9,10)]:
        val=sum((-1)**j*comb(2,j)*max(F(0),1-j*theta) for j in range(3))
        assert val==max(F(0),2*theta-1)
    print('PASS: circle second moments from geometric integrals, including one-arc variance zero')


def remaining():
    x,h,y,q,t = sp.symbols('x h y q t', positive=True)
    assert sp.simplify(sp.exp(-sp.exp(-y-h))*sp.exp(-(1-sp.exp(-h))*sp.exp(-y))-sp.exp(-sp.exp(-y))) == 0
    # Uniform waiting times, timeout t in (0,1): independently integrate the success law.
    for timeout in [.1,.5,1.]:
        for ss in [0.,.2,2.]:
            direct=quad(lambda u:exp(-ss*u),0,timeout)[0]/(1-(1-timeout)*exp(-ss*timeout))
            rewritten=1-ss*quad(lambda u:exp(-ss*u)*(1-u),0,timeout)[0]/(1-(1-timeout)*exp(-ss*timeout))
            assert abs(direct-rewritten)<1e-12
    local=(1-q)**-2*(1-2*q*(1-q)/(1-q*q)+q*q*(1-q)/(1-q**3))
    assert sp.simplify(local-(1+q*q)/((1-q*q)*(1-q**3)))==0
    log_local=sum(q**(j+1)*(1-q**j)/(j*(1-q**(j+1))) for j in range(1,5))
    assert sp.series(log_local,q,0,4).removeO() == q*q-q**3/2
    # Uniform perpetuity, independently from the Poisson-integral cumulants.
    # Its r-th cumulant is integral_0^1 u^r du/u = 1/r.
    moments=[F(1)]
    for r in range(1,7):
        moments.append(sum(F(comb(r-1,j-1),j)*moments[r-j] for j in range(1,r+1)))
        fixed=F(1,r)*sum(F(comb(r,j))*moments[j] for j in range(r))
        assert fixed == moments[r]
    print('PASS: Gumbel invariant law, fisherman transform, divisor factors, perpetuity moments')
    print('Divisor second moment:',mp.zeta(2)**2*mp.zeta(3)/mp.zeta(4))


if __name__=='__main__':
    for check in [socks,permutations_and_seats,coupon_finite,hitting_times,arrows,circle,remaining]:
        check()
    print('All finite-model and formula checks passed. Asymptotic proofs require mathematical review.')
