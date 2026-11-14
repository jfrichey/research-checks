#!/usr/bin/env python3
"""Reproducible checks for word-potentials.tex (2026-10-11).

Exact certificates and floating-point conjecture scans are reported separately.
Run: python3 research/sft/gibbs/word_potentials.py --scan --output PATH
Initial implementation developed collaboratively with ChatGPT by OpenAI from
Jacob Richey's mathematical direction. Research code; no independent audit.
"""
import argparse
import itertools
import json
import math
from fractions import Fraction
from pathlib import Path

import numpy as np
import sympy as sp
from scipy.linalg import eig
from scipy.optimize import brentq


def borders(word):
    return tuple(k for k in range(1, len(word) + 1) if word[:k] == word[-k:])


def prefix_matrix(word, q, t):
    """Counting transfer matrix; derivative is with respect to log(t)."""
    m = len(word)
    matrix = np.zeros((m, m))
    derivative = np.zeros_like(matrix)
    for i in range(m):
        for letter in range(q):
            suffix = word[:i] + (letter,)
            hit = suffix[-m:] == word
            j = max(k for k in range(m) if suffix[len(suffix)-k:] == word[:k])
            matrix[i, j] += t if hit else 1
            if hit:
                derivative[i, j] += t
    return matrix, derivative


def pressure_density(word, q, t):
    matrix, derivative = prefix_matrix(word, q, t)
    values, left, right = eig(matrix, left=True, right=True)
    j = np.argmax(values.real)
    lam = float(values[j].real)
    l, r = left[:, j].real, right[:, j].real
    density = float(l @ derivative @ r / (lam * (l @ r)))
    return math.log(lam), density


def characteristic(word, q, x, t):
    m = len(word)
    proper = sum(x**(b-1) for b in borders(word) if b < m)
    return sp.expand((x-q)*x**(m-1) - (t-1)*(1+(x-q)*proper))


class Interval:
    """Closed rational intervals; elementary operations use exact arithmetic."""
    def __init__(self, a, b=None):
        self.a = Fraction(a)
        self.b = Fraction(a if b is None else b)

    def __add__(self, other):
        other = other if isinstance(other, Interval) else Interval(other)
        return Interval(self.a+other.a, self.b+other.b)

    __radd__ = __add__

    def __mul__(self, other):
        other = other if isinstance(other, Interval) else Interval(other)
        products = [self.a*other.a, self.a*other.b, self.b*other.a, self.b*other.b]
        return Interval(min(products), max(products))

    __rmul__ = __mul__

    def __truediv__(self, other):
        other = other if isinstance(other, Interval) else Interval(other)
        assert other.a*other.b > 0, 'Division interval contains zero'
        return self * Interval(1/other.b, 1/other.a)


def evaluate_polynomial(poly, interval):
    result = Interval(0)
    for coefficient in poly.all_coeffs():
        result = result*interval + Fraction(coefficient)
    return result


def exact_density_certificate(word, t_value):
    x, t = sp.symbols('x t')
    polynomial = characteristic(word, 2, x, t)
    specialized = sp.Poly(polynomial.subs(t, t_value), x)
    roots = specialized.intervals(eps=sp.Rational(1, 10**20))
    (a, b), multiplicity = roots[-1]
    assert a > 2 and multiplicity == 1
    root = Interval(Fraction(a), Fraction(b))
    numerator = sp.Poly((-t*sp.diff(polynomial, t)).subs(t, t_value), x)
    denominator = sp.Poly((x*sp.diff(polynomial, x)).subs(t, t_value), x)
    density = evaluate_polynomial(numerator, root)/evaluate_polynomial(denominator, root)
    scale = 10**9
    low = (density.a.numerator*scale)//density.a.denominator
    high = -((-density.b.numerator*scale)//density.b.denominator)
    assert Fraction(low, scale) <= density.a <= density.b <= Fraction(high, scale)
    return {
        'word': ''.join(map(str, word)), 't': t_value,
        'root_interval': [str(a), str(b)],
        'density_interval': [f'{low/scale:.9f}', f'{high/scale:.9f}'],
        'method': 'exact rational root isolation and rational interval evaluation',
    }, density


def full_matrix(word, t):
    m = len(word)
    states = list(itertools.product(range(2), repeat=m-1))
    index = {state: i for i, state in enumerate(states)}
    matrix = sp.zeros(len(states))
    for i, state in enumerate(states):
        for a in range(2):
            edge = state+(a,)
            matrix[i, index[edge[1:]]] += t if edge == word else 1
    return matrix


def scc_sizes(matrix):
    n = matrix.rows
    reach = [[i == j or matrix[i,j] != 0 for j in range(n)] for i in range(n)]
    for k in range(n):
        for i in range(n):
            for j in range(n):
                reach[i][j] = reach[i][j] or (reach[i][k] and reach[k][j])
    unseen = set(range(n)); sizes = []
    while unseen:
        i = min(unseen)
        component = {j for j in unseen if reach[i][j] and reach[j][i]}
        unseen -= component; sizes.append(len(component))
    return sorted(sizes)


def validate():
    checks = 0
    x, t = sp.symbols('x t')
    for m in range(2, 6):
        for word in itertools.product(range(2), repeat=m):
            n = m+2
            expected = sum(2**sum(s[i:i+m] == word for i in range(n-m+1))
                           for s in itertools.product(range(2), repeat=n))
            matrix, _ = prefix_matrix(word, 2, 2)
            actual = (np.linalg.matrix_power(matrix, n) @ np.ones(m))[0]
            assert abs(actual-expected) < 1e-8
            polynomial = characteristic(word, 2, x, t)
            for tv in [.5, 2]:
                pressure, density = pressure_density(word, 2, tv)
                lam = math.exp(pressure)
                residual = float(polynomial.subs({x: lam, t: tv}))
                f_x = float(sp.diff(polynomial, x).subs({x: lam, t: tv}))
                f_t = float(sp.diff(polynomial, t).subs({x: lam, t: tv}))
                assert abs(residual) < 1e-9
                assert abs(density+tv*f_t/(lam*f_x)) < 1e-10
            checks += 1
    certificates = []
    for tv in [2, 64]:
        u, iu = exact_density_certificate(tuple(map(int, '01001')), tv)
        v, iv = exact_density_certificate(tuple(map(int, '00100')), tv)
        assert (iv.a > iu.b) if tv == 2 else (iu.a > iv.b)
        certificates.extend([u, v])
    z = sp.symbols('z')
    obstruction = []
    for word in ['0001', '0011']:
        matrix = full_matrix(tuple(map(int, word)), t)
        determinant = sp.expand((sp.eye(matrix.rows)-z*matrix).det())
        assert sp.expand(determinant - (1-2*z+(1-t)*z**4)) == 0
        sizes = scc_sizes(matrix.subs(t, 0))
        assert sizes == ([1, 7] if word == '0001' else [8])
        obstruction.append({'word': word, 'det_I_minus_z_M': str(determinant),
                            'avoidance_scc_sizes': sizes})
    return {'validation_words': checks, 'density_reversal_certificates': certificates,
            'conjugacy_obstruction_checks': obstruction}


def representatives(q, m):
    reps = {}
    for word in itertools.product(range(q), repeat=m):
        reps.setdefault(borders(word), word)
    profiles = sorted(reps, key=lambda b: sum(q**(k-1) for k in b))
    return profiles, [reps[b] for b in profiles]


def border_profiles(maximum):
    """Enumerate all realizable profiles over arbitrary finite alphabets.

    If b is the longest proper border, the smaller borders are exactly the
    borders of that b-letter prefix. For every recursively generated candidate,
    impose its letter equalities and check that their equivalence classes force
    no extra border. Giving each class its own letter realizes every survivor.
    Thus every binary profile is included, without a binary-realization input.
    """
    rows = {0: {(): ()}}
    for m in range(1, maximum+1):
        profiles = {}
        for b in range(m):
            for prefix in rows[b]:
                candidate = prefix+(m,)
                parent = list(range(m))

                def root(i):
                    while parent[i] != i:
                        parent[i] = parent[parent[i]]
                        i = parent[i]
                    return i

                for border in prefix:
                    for i in range(border):
                        parent[root(i)] = root(m-border+i)
                classes = {}; word = []
                for i in range(m):
                    word.append(classes.setdefault(root(i), len(classes)))
                word = tuple(word)
                if borders(word) == candidate:
                    profiles[candidate] = word
        rows[m] = profiles
    return rows


def pair_crossing_certificates(u, v, m, q=2):
    """Exact gap and conditional-mean criteria; u,v are proper period sets."""
    new = v-u
    a = min(new)
    offsets = [n-a for n in new]
    scale = max(offsets)
    ell = a-max(n for n in u if n < a)
    j = min([n for n in u if n > a]+[m])-a
    gap_numerator = sum((n+ell)*(n-j)*q**(scale-n) for n in offsets)
    if gap_numerator <= 0:
        return 'gap', {'strict_gap': gap_numerator < 0}
    common = u & v
    removed = u-v
    max_common = max(common)
    numerator = sum(p*p*q**(max_common-p) for p in common)
    denominator = sum(p*q**(max_common-p) for p in common)
    ell_numerator = a*denominator-numerator
    j = min(removed | {m})-a
    mean_numerator = sum((n*denominator+ell_numerator)*(n-j)*q**(scale-n)
                         for n in offsets)
    details = {'gap_criterion': str(Fraction(gap_numerator, q**scale)),
               'common_period_mean': str(Fraction(numerator, denominator)),
               'first_new_period': a, 'upper_gap': j,
               'mean_criterion': str(Fraction(mean_numerator, denominator*q**scale))}
    kind = 'common_mean' if ell_numerator > 0 and mean_numerator <= 0 else 'uncovered'
    return kind, details


def binary_proof_checks():
    """Exact arithmetic for the finite and infinite bounds in the binary proof."""
    affine = [(Fraction(5, 9), Fraction(56, 27)),
              (Fraction(7, 16), Fraction(3, 2)),
              (Fraction(7, 16), Fraction(17, 8)),
              (Fraction(25, 64), Fraction(59, 32)),
              (Fraction(3, 8), Fraction(17, 8)),
              (Fraction(5, 16), Fraction(2))]
    for slope, intercept in affine:
        assert slope < 1 and 5*slope+intercept < 5

    def multiples(step, offset):
        z = Fraction(1, 2**step)
        return step*step*z*(1+z)/(1-z)**3-offset*step*z/(1-z)**2

    period_two = 4*(multiples(2, 3)-multiples(6, 3))
    assert period_two == Fraction(-1376, 27783)
    rows = []
    expected = {3: [-3, -3, Fraction(-7, 4), Fraction(-3, 2),
                    Fraction(-3, 16), Fraction(-1, 2), Fraction(-51, 64)],
                4: [-4, -4, Fraction(-5, 2), Fraction(-3, 4),
                    Fraction(-3, 4), Fraction(-11, 8), Fraction(-19, 32)]}
    for k, last, tail in [(3, 10, Fraction(15, 32)), (4, 11, Fraction(83, 64))]:
        bounds = []
        for d in range(k+1, last+1):
            support = set(range(k, d, k)) | set(range(max(k+1, d-k+2), d))
            value = sum(Fraction(s*(s-k-1), 2**(s-k)) for s in support)
            assert value < 0
            bounds.append(value)
        assert bounds == expected[k]
        infinite = 2**k*multiples(k, k+1)
        assert infinite+tail < 0
        rows.append({'least_period': k, 'finite_lengths': [k+1, last],
                     'finite_upper_bounds': list(map(str, bounds)),
                     'infinite_multiples_sum': str(infinite),
                     'additional_tail_bound': str(tail),
                     'large_length_upper_bound': str(infinite+tail)})
    return {'status': 'exact arithmetic checks for the working binary proof',
            'checked_on': '2026-10-11',
            'nonconstant_moment_affine_bounds': [
                {'slope': str(a), 'intercept': str(b), 'slack_at_five': str(5-5*a-b)}
                for a, b in affine],
            'period_two_excluding_multiples_three': str(period_two),
            'small_period_bounds': rows}


def marker_swap_checks():
    """Exhaustive local certificate for 1100011 <-> 1101011."""
    u = tuple(map(int, '1100011'))
    v = tuple(map(int, '1101011'))

    def rule(window):
        trigger = window[:3] == (1, 1, 0) and window[4:] == (0, 1, 1)
        return 1-window[3] if trigger else window[3]

    count = 0
    for window in itertools.product(range(2), repeat=13):
        image = tuple(rule(window[i:i+7]) for i in range(7))
        assert rule(image) == window[6]
        central = window[3:10]
        assert (central == u) == (image == v)
        assert (central == v) == (image == u)
        assert (window[6:8] == (1, 1)) == (image[3:5] == (1, 1))
        count += 1
    assert borders(u) == borders(v) == (1, 2, 7)
    return {'status': 'exhaustive finite-window proof of the local identities',
            'checked_on': '2026-10-11', 'u': '1100011', 'v': '1101011',
            'border_lengths': [1, 2, 7], 'radius': 3, 'input_window_length': 13,
            'windows_checked': count,
            'identities': ['h squared is identity', 'u indicator equals v after h',
                           'v indicator equals u after h', '11-marker field unchanged']}


def increment_lift_checks():
    """Exact partition checks and an exhaustive parity-preserving lifted code."""
    finite_cases = 0
    for q, maximum in [(2, 5), (3, 3)]:
        for m in range(2, maximum+1):
            n = m+2
            for word in itertools.product(range(q), repeat=m):
                orbit = {tuple((a+c) % q for a in word) for c in range(q)}
                difference = tuple((word[i+1]-word[i]) % q for i in range(m-1))
                lhs = sum(2**sum(x[i:i+m] in orbit for i in range(n-m+1))
                          for x in itertools.product(range(q), repeat=n))
                rhs = q*sum(2**sum(y[i:i+m-1] == difference for i in range(n-m+1))
                            for y in itertools.product(range(q), repeat=n-1))
                assert lhs == rhs
                finite_cases += 1

    def code(word):
        return sum(int(a) << i for i, a in enumerate(word))

    u, v = '0100000010', '0100110010'
    complement = lambda w: ''.join(str(1-int(a)) for a in w)
    f_words, g_words = [u, complement(u)], [v, complement(v)]
    f, g = list(map(code, f_words)), list(map(code, g_words))

    def selected(words):
        result = np.zeros(words.shape, dtype=bool)
        for value in f+g:
            result |= words == value
        return result

    def lift(window):
        return ((window >> 5) & 1) ^ (selected(window & 1023) |
                                      selected((window >> 1) & 1023))

    # Bit i represents the i-th symbol, with no floating-point computation.
    windows = np.arange(1 << 21, dtype=np.uint32)
    images = np.zeros(windows.shape, dtype=np.uint16)
    for i in range(11):
        images |= np.asarray(lift((windows >> i) & 2047), dtype=np.uint16) << i
    assert np.all(lift(images) == ((windows >> 10) & 1))
    before, after = (windows >> 5) & 1023, images & 1023
    assert np.all(((before == f[0]) | (before == f[1])) ==
                  ((after == g[0]) | (after == g[1])))
    assert np.all(((before == g[0]) | (before == g[1])) ==
                  ((after == f[0]) | (after == f[1])))

    short = np.arange(1 << 12, dtype=np.uint16)
    differences = (short ^ (short >> 1)) & 2047
    trigger = np.zeros(short.shape, dtype=bool)
    for value in map(code, ['110000011', '110101011']):
        trigger |= ((differences & 511) == value) | (((differences >> 2) & 511) == value)
    base_image = ((differences >> 5) & 1) ^ trigger
    assert np.all((lift(short & 2047) ^ lift((short >> 1) & 2047)) == base_image)

    t = sp.symbols('t')
    same, different = sp.Matrix([[t, 1], [1, t]]), sp.Matrix([[1, t], [t, 1]])
    assert same*sp.ones(2, 1) == different*sp.ones(2, 1) == (t+1)*sp.ones(2, 1)
    assert sp.trace(same)-sp.trace(different) == 2*t-2
    return {'status': 'exact finite enumeration and integer-vectorized local proof',
            'checked_on': '2026-10-11', 'partition_function_cases': finite_cases,
            'partition_checks': {'q_2_lengths': [2, 5], 'q_3_lengths': [2, 3],
                                 'n': 'm+2', 'fugacity': 2},
            'f_words': f_words, 'g_words': g_words, 'radius': 5,
            'involution_and_indicator_window_length': 21,
            'involution_and_indicator_windows': 1 << 21,
            'factor_identity_window_length': 12, 'factor_identity_windows': 1 << 12,
            'two_letter_obstruction': {'common_row_sum': str(t+1),
                                      'formal_trace_difference': str(2*t-2)}}


def many_crossings_certificate(N):
 """Rational parameters and convex-secant bounds for the finite construction."""
 F = Fraction
 assert isinstance(N, int) and N >= 1
 L=2*(N+2);D=2**L;B=4*(N+2)**2;T=2*N+3
 E=2*L*(B+T)+T+1+(64*L*D-1).bit_length()
 k=[L//2+i for i in range(N+2)]
 cycles={}
 for i,ki in enumerate(k):
  word=(1,)*ki+(0,)*(L-ki)
  for complement in [False,True]:
   cycle_word=tuple(1-a for a in word) if complement else word
   rotations={cycle_word[j:]+cycle_word[:j] for j in range(L)}
   assert len(rotations)==L
   for state in rotations:
    assert state not in cycles or (i==0 and cycles[state]==0)
    cycles[state]=i
 assert len(cycles)==(2*N+3)*L
 backgrounds=[]
 for parity in [1,0]:
  intercepts=[-B]*(N+2);intercepts[0]=0;previous=0
  for i in range(1,N+2):
   if i%2==parity:
    intercepts[i]=intercepts[previous]-(k[i]-k[previous])*(2*i-1)
    previous=i
  assert all(-B<=a<=0 for a in intercepts)
  backgrounds.append(intercepts)
 assert 2**(E-2*L*(B+T)-T-1)>=64*L*D
 samples=[]
 for j in range(1,N+2):
  expected=[]
  for parity,intercepts in zip([1,0],backgrounds):
   index=max([0]+[i for i in range(1,j+1) if i%2==parity])
   slope=F(k[index],L);expected.append(slope)
   lines=[(a,b) for a,ki in zip(intercepts,k) for b in [ki,L-ki]]
   for beta_units in [F(2*j)-F(1,2),F(2*j),F(2*j)+F(1,2)]:
    largest=max(a+b*beta_units for a,b in lines)
    slopes={F(b,L) for a,b in lines if a+b*beta_units==largest}
    assert slopes=={slope}
  assert expected[0]-expected[1]==F(1 if j%2 else -1,L)
  error=F(1,16*L)
  intervals=[(r-error,r+error) for r in expected]
  assert (intervals[0][0]>intervals[1][1]) if j%2 else (intervals[1][0]>intervals[0][1])
  samples.append({'coupling':f'{2*j}*log(2)','limiting_frequencies':list(map(str,expected)),'certified_frequency_intervals':[[str(a),str(b)] for a,b in intervals]})
 return {'N':N,'cycle_length':L,'vertices':D,'local_background_range':L+1,'used_cycle_states':len(cycles),'observed_word':'1','cycle_mark_counts':k,'cycle_weight_log2':backgrounds,'epsilon':f'2^(-{E})','pressure_error_bound':str(F(1,64*L)),'frequency_error_bound':str(F(1,16*L)),'reference_frequencies':['1/2','1/2'],'samples':samples}

def many_crossing_checks():
    for N in range(1, 13):
        many_crossings_certificate(N)
    return {'status': 'exact parameter, envelope, and comparison-bound certificate',
            'checked_on': '2026-10-11', 'parameter_checks_N': [1, 12],
            'graph': 'Binary de Bruijn graph on L-letter states. Use cycles of '
                     '1^k 0^(L-k) and their complements, sharing the central cycle. '
                     'Put the displayed weight on one distinguished edge per '
                     'noncentral cycle, and 1 on its other edges. Add epsilon to '
                     'the base weight of every de Bruijn edge. Tilt by the first '
                     'bit of the current state.',
            'bound_method': 'positive-vector Perron comparison and convex secants; '
                            'no floating-point eigenvalues',
            'example': many_crossings_certificate(3)}


def weighted_overlap_data(word, probabilities):
    eta = math.prod(probabilities[a] for a in word)
    periods = {len(word)-b for b in borders(word) if b < len(word)}
    theta = tuple(math.prod(probabilities[a] for a in word[-s:]) if s in periods
                  else Fraction(0) for s in range(1, len(word)))
    return eta, theta


def weighted_profile_checks():
    """Transfer-matrix consistency plus an exact biased pressure-crossing example."""
    comparisons = 0
    for probabilities in [(Fraction(7, 20), Fraction(13, 20)),
                          (Fraction(1, 5), Fraction(4, 5))]:
        for m in range(2, 6):
            for word in itertools.product(range(2), repeat=m):
                eta, theta = weighted_overlap_data(word, probabilities)
                for t in [.5, 2]:
                    matrix = np.zeros((m, m)); derivative = np.zeros_like(matrix)
                    for i in range(m):
                        for a in range(2):
                            suffix = word[:i]+(a,)
                            hit = suffix[-m:] == word
                            j = max(k for k in range(m)
                                    if suffix[len(suffix)-k:] == word[:k])
                            weight = float(probabilities[a])*(t if hit else 1)
                            matrix[i, j] += weight
                            if hit:
                                derivative[i, j] += weight
                    values, left, right = eig(matrix, left=True, right=True)
                    j = np.argmax(values.real); x = float(values[j].real)
                    l, r = left[:, j].real, right[:, j].real
                    frequency = float(l @ derivative @ r/(x*(l @ r)))
                    L = m-1
                    d = sum(float(c)*x**(L-s) for s, c in enumerate(theta, 1))
                    dp = sum((L-s)*float(c)*x**(L-s-1)
                             for s, c in enumerate(theta, 1) if s < L)
                    residual = (x-1)*x**L-(t-1)*(float(eta)+(x-1)*d)
                    f_x = x**L+(x-1)*L*x**(L-1)-(t-1)*(d+(x-1)*dp)
                    predicted = t*(float(eta)+(x-1)*d)/(x*f_x)
                    assert abs(residual) < 1e-10
                    assert abs(frequency-predicted) < 1e-9
                    comparisons += 1
    x = sp.symbols('x')
    probabilities = (Fraction(7, 20), Fraction(13, 20))
    polynomials = []; masses = []
    for word in ['1110000111', '0011111100']:
        eta, theta = weighted_overlap_data(tuple(map(int, word)), probabilities)
        masses.append(eta)
        polynomials.append(sum(sp.Rational(c.numerator, c.denominator)*x**(9-s)
                               for s, c in enumerate(theta, 1)))
    assert masses[0] == masses[1]
    q = 9800*x*x-15600*x-3549
    factor = sp.Rational(40, 107653)*sp.Rational(masses[0].numerator,
                                               masses[0].denominator)
    assert sp.expand(polynomials[0]-polynomials[1]-factor*q) == 0
    intervals = sp.Poly(q, x).intervals(eps=sp.Rational(1, 10**15))
    assert len(intervals) == 2 and intervals[0][0][1] < 0 < 1 < intervals[1][0][0]
    return {'checked_on': '2026-10-11',
            'transfer_validation': {'status': 'floating-point formula consistency',
                                    'comparisons': comparisons, 'length_range': [2, 5],
                                    'one_probabilities': ['13/20', '4/5'],
                                    'fugacities': ['1/2', '2']},
            'pressure_crossing': {'status': 'exact rational polynomial certificate',
                                  'u': '1110000111', 'v': '0011111100',
                                  'one_probability': '13/20',
                                  'common_word_probability': str(masses[0]),
                                  'normalized_D_u_minus_D_v': str(q),
                                  'positive_root_interval': list(map(str, intervals[1][0])),
                                  'polynomial_value_at_one': -9349}}


def polynomial_density(coefficients, q, t):
    """Algebraic root branch, also for binary polynomials not realized by words."""
    L = len(coefficients)

    def values(x):
        d = sum(c*x**j for j, c in enumerate(coefficients))
        dp = sum(j*c*x**(j-1) for j, c in enumerate(coefficients) if j)
        return d, dp

    def characteristic_at(x):
        d, _ = values(x)
        return (x-q)*x**L-(t-1)*(1+(x-q)*d)

    bracket = (q*L/(L+1), q) if t < 1 else (q, q+t)
    x = brentq(characteristic_at, *bracket, xtol=1e-13)
    d, dp = values(x)
    f_x = x**L+(x-q)*L*x**(L-1)-(t-1)*(d+(x-q)*dp)
    return t*(1+(x-q)*d)/(x*f_x)


def polynomial_scan():
    """Adversarial finite check of both new orders without word realizability."""
    negative_t = [.0001, .02, .2, .8, .99]
    positive_beta = [.001, .03, .1, .3, 1, 2, 3, 4, 5, 6, 8, 12, 20]
    result = {'status': 'floating-point consistency check, not the proof',
              'negative': [], 'positive': []}
    for q in [2, 3, 4]:
        failures = []; count = 0
        for L in range(1, 9):
            profiles = sorted(itertools.product([0, 1], repeat=L),
                              key=lambda p: sum(c*2**j for j, c in enumerate(p)))
            count += len(profiles)
            for t in negative_t:
                frequencies = [polynomial_density(p, q, t) for p in profiles]
                for i, (a, b) in enumerate(zip(frequencies, frequencies[1:])):
                    if a < b-1e-12:
                        failures.append({'L': L, 't': t,
                                         'u': profiles[i], 'v': profiles[i+1]})
        result['negative'].append({'q': q, 'L_range': [1, 8],
            'profiles_total': count, 't_grid': negative_t,
            'tolerance': 1e-12, 'failures': failures})
    failures = []; count = 0
    for L in range(1, 12):
        groups = {}
        for p in itertools.product([0, 1], repeat=L):
            degree = max([j for j, c in enumerate(p) if c]+[-1])
            groups.setdefault(degree, []).append(p)
            count += 1
        for beta in positive_beta:
            extrema = []
            for degree, profiles in sorted(groups.items()):
                values = [(polynomial_density(p, 2, math.exp(beta)), p)
                          for p in profiles]
                extrema.append((degree, min(values), max(values)))
            for a, b in zip(extrema, extrema[1:]):
                if a[2][0] > b[1][0]+1e-11:
                    failures.append({'L': L, 'beta': beta,
                                     'lower_degree_maximum': a,
                                     'higher_degree_minimum': b})
    result['positive'].append({'q': 2, 'L_range': [1, 11],
        'profiles_total': count, 'beta_grid': positive_beta,
        'tolerance': 1e-11, 'failures': failures})
    return result


def crossing_certificates(maximum=70):
    """Exact hypothesis checks on exhaustively generated border profiles.

    Coverage contains every binary word. The check uses integer-scaled rational
    arithmetic and the proved theorems, without sampling any coupling.
    """
    all_profiles = border_profiles(maximum)
    validation = []
    for q, limit in [(2, min(16, maximum)), (3, min(9, maximum))]:
        for m in range(1, limit+1):
            direct = {borders(w) for w in itertools.product(range(q), repeat=m)}
            assert set(all_profiles[m]) == direct
        validation.append({'q': q, 'length_range': [1, limit],
                           'method': 'equality with exhaustive word enumeration'})
    rows = []
    for m in range(2, maximum+1):
        groups = {}
        for profile in sorted(all_profiles[m], key=lambda b: sum(2**k for k in b)):
            if len(profile) > 1:
                groups.setdefault(m-profile[-2], []).append(profile)
        counts = {'gap': 0, 'common_mean': 0, 'uncovered': 0}
        inclusions = noninclusions = 0
        uncovered = []; resolved = []
        for group in groups.values():
            for i, u in enumerate(group):
                pu = {m-b for b in u[:-1]}
                max_u = max(pu)
                mean_numerator = sum(p*p*2**(max_u-p) for p in pu)
                mean_denominator = sum(p*2**(max_u-p) for p in pu)
                for v in group[i+1:]:
                    pv = {m-b for b in v[:-1]}
                    kind, details = pair_crossing_certificates(pu, pv, m)
                    counts[kind] += 1
                    if pu <= pv:
                        inclusions += 1
                        assert mean_numerator < min(pv-pu)*mean_denominator
                    else:
                        noninclusions += 1
                        assert kind == 'gap' and details['strict_gap']
                    if kind == 'uncovered':
                        uncovered.append({'u_borders': u, 'v_borders': v, **details})
                    elif kind == 'common_mean' and m <= 20:
                        resolved.append({'u_borders': u, 'v_borders': v, **details})
        rows.append({'m': m, 'realized_border_classes': len(all_profiles[m]),
                     'same_period_pairs': sum(counts.values()),
                     'rational_gap_criterion': counts['gap'],
                     'common_mean_criterion': counts['common_mean'],
                     'inclusion_mean_assertions': inclusions,
                     'noninclusion_strict_gap_assertions': noninclusions,
                     'small_common_mean_examples': resolved, 'uncovered': uncovered})
    return {'status': 'exact enumeration and rational theorem-hypothesis checks',
            'q': 2, 'length_range': [2, maximum],
            'enumeration': 'recursive longest-border candidates and equality closure; '
                           'all finite-alphabet profiles, containing every binary profile',
            'enumerator_validation': validation, 'rows': rows}


def biased_overlap_scan():
    """Finite tests of the weighted-overlap conjectures; all cases are retained."""
    def profile(w):
     return tuple((b,sum(w[:b])) for b in range(1,len(w)) if w[:b]==w[-b:])

    def transitions(w):
     m=len(w); out=[]
     for i in range(m):
      for a in range(2):
       suffix=w[:i]+(a,);hit=suffix[-m:]==w
       j=max(k for k in range(m) if suffix[len(suffix)-k:]==w[:k])
       out.append((i,j,a,hit))
     return out

    def density(w,edges,bias,beta):
     m=len(w);t=math.exp(beta);A=np.zeros((m,m));D=np.zeros((m,m))
     for i,j,a,hit in edges:
      z=bias if a else 1-bias
      A[i,j]+=z*(t if hit else 1)
      if hit:D[i,j]+=z*t
     values,l,r=eig(A,left=True,right=True);idx=np.argmax(values.real)
     left,right=l[:,idx].real,r[:,idx].real;lam=float(values[idx].real)
     return math.log(lam),float(left@D@right/(lam*(left@right)))

    result=[]
    negative=[-12,-8,-4,-2,-1,-.3,-.05]
    positive=list(np.linspace(.05,8,80))+list(np.linspace(8.2,32,100))
    for m in range(3,13):
     groups={}
     for w in itertools.product(range(2),repeat=m):
      if 0<sum(w)<m:
       groups.setdefault(sum(w),{}).setdefault(profile(w),w)
     cases=[];pairs=[]
     for number,group in groups.items():
      base=len(cases)
      for pr,w in group.items():cases.append((w,pr,transitions(w)))
      items=list(group)
      for i,u in enumerate(items):
       for j,v in enumerate(items):
        if u!=v and set(dict(u))<=set(dict(v)) and all(dict(v)[b]<=k for b,k in u): pairs.append((base+i,base+j))
     for bias in [.65,.8,.95]:
      values=np.array([[density(w,edges,bias,float(beta)) for w,pr,edges in cases] for beta in negative+positive])
      fails=[]
      for i,j in pairs:
       diff=values[:,j,1]-values[:,i,1]
       pressure_diff=values[:,j,0]-values[:,i,0]
       assert np.min(pressure_diff)>-1e-9
       if max(diff[:len(negative)])>1e-11:
        k=int(np.argmax(diff[:len(negative)]));fails.append({'type':'negative_reversal','beta':negative[k],'delta':float(diff[k]),'u':cases[i][0],'v':cases[j][0],'profiles':[cases[i][1],cases[j][1]]})
       if tuple(b for b,k in cases[i][1])[-1:]!=tuple(b for b,k in cases[j][1])[-1:]:
        pos=diff[len(negative):]
        if min(pos)<-1e-10:
         k=int(np.argmin(pos));fails.append({'type':'distinct_period_positive_reversal','beta':float(positive[k]),'delta':float(pos[k]),'u':cases[i][0],'v':cases[j][0],'profiles':[cases[i][1],cases[j][1]]})
       else:
        pos=diff[len(negative):];sig=np.sign(pos[np.abs(pos)>1e-10]);changes=int(np.count_nonzero(sig[1:]!=sig[:-1]))
        if changes>1:fails.append({'type':'multiple_positive_crossings','changes':changes,'u':cases[i][0],'v':cases[j][0],'profiles':[cases[i][1],cases[j][1]]})
      result.append({'m':m,'bias':bias,'classes':len(cases),'pairs':len(pairs),'failures':fails})

    return {'status': 'floating-point evidence only, not a universal assertion',
            'checked_on': '2026-10-11', 'q': 2, 'length_range': [3,12],
            'bias_means_probability_of_one': [.65,.8,.95],
            'hypothesis': 'equal word composition and coordinatewise weighted overlap domination',
            'negative_beta': negative,
            'positive_beta_grid': [[.05,8,80],[8.2,32,100]],
            'negative_tolerance': 1e-11, 'positive_tolerance': 1e-10,
            'rows': result}


def scan():
    result = {'status': 'floating-point evidence, not universal inequalities',
              'comparison_tolerance': 1e-10, 'order_scans': [], 'crossing_scans': []}
    negative = [-12, -6, -3, -1, -.2, -.02]
    positive = [.02, .2, 1, 3, 6, 12, 24]
    for q, maximum in [(2, 16), (3, 9), (4, 7)]:
        total = 0; failures = []
        for m in range(2, maximum+1):
            profiles, words = representatives(q, m); total += len(words)
            for beta in negative+positive:
                frequencies = [pressure_density(w, q, math.exp(beta))[1] for w in words]
                for i in range(len(words)):
                    for j in range(i+1, len(words)):
                        if beta < 0:
                            bad = frequencies[j] > frequencies[i]+1e-10
                        else:
                            different_period = profiles[i][-2:-1] != profiles[j][-2:-1]
                            bad = different_period and frequencies[j] < frequencies[i]-1e-10
                        if bad:
                            failures.append({'u': words[i], 'v': words[j], 'beta': beta})
        result['order_scans'].append({'q': q, 'lengths': [2, maximum],
              'border_classes_total': total, 'negative_beta': negative,
              'positive_beta': positive, 'failures': failures})
    grid = np.concatenate([np.linspace(.02, 8, 140), np.linspace(8.2, 30, 90)])
    for m in [5, 8, 12, 16]:
        profiles, words = representatives(2, m)
        frequencies = np.array([[pressure_density(w, 2, math.exp(float(beta)))[1]
                                  for w in words] for beta in grid])
        maximum_changes = 0; failures = []
        for i in range(len(words)):
            for j in range(i+1, len(words)):
                difference = frequencies[:,j]-frequencies[:,i]
                signs = np.sign(difference[np.abs(difference) > 1e-9])
                changes = int(np.count_nonzero(signs[1:] != signs[:-1]))
                maximum_changes = max(maximum_changes, changes)
                if changes > 1:
                    failures.append({'u': words[i], 'v': words[j], 'changes': changes})
        result['crossing_scans'].append({'q': 2, 'm': m, 'classes': len(words),
            'beta_grid': {'segments': [[.02, 8, 140], [8.2, 30, 90]]},
            'sign_tolerance': 1e-9, 'maximum_observed_sign_changes': maximum_changes,
            'multiple_crossings': failures})
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--scan', action='store_true')
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    report = {'date': '2026-10-11', 'exact_and_validation': validate(),
              'binary_proof_checks': binary_proof_checks(),
              'marker_swap_checks': marker_swap_checks(),
              'increment_lift_checks': increment_lift_checks(),
              'many_crossing_checks': many_crossing_checks(),
              'weighted_profile_checks': weighted_profile_checks()}
    if args.scan:
        report['numerical_scans'] = scan()
        report['polynomial_scans'] = polynomial_scan()
        report['crossing_certificates'] = crossing_certificates()
        report['biased_overlap_scan'] = biased_overlap_scan()
    payload = json.dumps(report, indent=2)+'\n'
    if args.output:
        args.output.write_text(payload)
        print(f'Checks passed; evidence saved to {args.output}')
    else:
        print(payload, end='')
