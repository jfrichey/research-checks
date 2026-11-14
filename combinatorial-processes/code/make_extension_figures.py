"""Build three deterministic extension figures and their numerical records.

No random samples, capped optimal searches, or silently omitted greedy inputs.
PDFs and reusable data belong to the project; all previews stay under /tmp.
"""
from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime, timezone
from fractions import Fraction
import hashlib
import json
from pathlib import Path

from make_figures import plt, panels, TEAL, ORANGE, PURPLE, GRAY
import numpy as np
from scipy.special import ndtr
from sympy import prevprime

from egyptian_fractions import greedy, shortest
from extension_calculations import adsorption_distributions, standardized_pmf, four_exception_sieve

ROOT = Path(__file__).resolve().parents[1]
RECORD = {'date': '2026-11-03', 'method': 'deterministic, no random seeds',
          'software': {'numpy': np.__version__}, 'figures': {}}


def save_json(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, sort_keys=True) + '\n')


def finish(fig, name, preview, **parameters):
    pdf = ROOT / 'figures' / (name+'.pdf')
    fig.savefig(pdf, metadata={'Creator': 'make_extension_figures.py', 'Author': 'ChatGPT',
        'CreationDate': datetime(2026, 11, 3, tzinfo=timezone.utc),
        'ModDate': datetime(2026, 11, 3, tzinfo=timezone.utc)})
    fig.savefig(preview/(name+'.png'), dpi=180)
    plt.close(fig)
    RECORD['figures'][name] = parameters | {'sha256': hashlib.sha256(pdf.read_bytes()).hexdigest()}


def adsorption(preview):
    pmfs = adsorption_distributions(640)
    data = {}
    ns = list(range(20, 641, 20))
    for n in ns:
        z, cdf, mean, var, third, dist = standardized_pmf(pmfs[n])
        assert abs(pmfs[n].sum()-1) < 2e-12
        assert pmfs[n].min() >= 0
        data[n] = {'mean': mean, 'variance': var, 'third_centered': third,
                   'kolmogorov': dist, 'sqrt_n_kolmogorov': np.sqrt(n)*dist}
    fig, ax = panels()
    for n, color in zip((40, 160, 640), (TEAL, PURPLE, ORANGE)):
        z, cdf, *_ = standardized_pmf(pmfs[n])
        ax[0].step(z, cdf, where='post', color=color, lw=1.15, label=f'n = {n}')
    xx = np.linspace(-3, 3, 400)
    ax[0].plot(xx, ndtr(xx), color=GRAY, lw=1.2, ls='--', label='normal')
    ax[0].set(xlim=(-3, 3), ylim=(-.015, 1.015), xlabel='standardized occupancy',
              ylabel='cumulative probability', title='The jumps shrink on the normal scale')
    ax[0].legend(loc='lower right')
    ax[1].plot(ns, [data[n]['sqrt_n_kolmogorov'] for n in ns], 'o-', ms=3, color=TEAL)
    ax[1].set(xlabel='number of sites n', ylabel=r'$\sqrt{n}\,\sup_x |F_n(x)-\Phi(x)|$',
              title='A finite-size check of the rate')
    ax[1].set_ylim(bottom=0)
    finish(fig, '04-adsorption-clt', preview, displayed_n=[40, 160, 640],
           distance_n=ns, arithmetic='double precision exact-recurrence evaluation')
    save_json(ROOT/'computations'/'adsorption-clt.json', {
        'method': 'PGF splitting recurrence, double precision; CDF includes both sides of jumps',
        'statistics': data,
        'pmfs': {n: pmfs[n].tolist() for n in (40, 160, 640)}})
    print('Built adsorption distributions through n=640.', flush=True)


def egyptian(preview):
    optimal = {}
    for n in (31, 32, 61, 64, 127, 128):
        reps = [shortest(k, n) for k in range(1, n)]
        for k, ds in enumerate(reps, 1):
            assert sum((Fraction(1, d) for d in ds), Fraction()) == Fraction(k, n)
            assert list(ds) == sorted(set(ds))
        counts = Counter(map(len, reps))
        optimal[n] = {'counts': dict(counts), 'mean': sum(map(len, reps))/(n-1),
                      'representations': reps}
    denominators = sorted({n for j in range(5, 16) for n in (2**j, int(prevprime(2**j)))})
    census = {}
    for n in denominators:
        counts = Counter(len(greedy(k, n, bit_limit=200_000)) for k in range(1, n))
        census[n] = {'counts': dict(counts), 'mean': sum(k*v for k, v in counts.items())/(n-1),
                     'max_length': max(counts), 'completed': sum(counts.values()), 'censored': 0}
        print(f'Full greedy census: N={n}, mean={census[n]["mean"]:.5f}.', flush=True)
    fig, ax = panels()
    lengths = np.arange(1, 8)
    ax[0].bar(lengths-.18, [optimal[127]['counts'].get(int(k), 0)/126 for k in lengths],
              width=.36, color=TEAL, label='minimum')
    ax[0].bar(lengths+.18, [census[127]['counts'].get(int(k), 0)/126 for k in lengths],
              width=.36, color=ORANGE, label='greedy')
    ax[0].set(xlabel='number of unit fractions', ylabel='proportion of numerators',
              title='The complete grid with N = 127', xticks=lengths)
    ax[0].legend()
    for values, color, label in [([2**j for j in range(5,16)], TEAL, 'powers of two'),
                               ([int(prevprime(2**j)) for j in range(5,16)], ORANGE, 'nearby primes')]:
        ax[1].plot(values, [census[n]['mean'] for n in values], 'o-', color=color, ms=3, label=label)
    ax[1].set_xscale('log', base=2)
    ax[1].set(xlabel='denominator N (log scale)', ylabel='mean greedy length',
              title='Arithmetic remains visible')
    ax[1].legend()
    finish(fig, '16-egyptian-lengths', preview, histogram_N=127, greedy_denominators=denominators,
           arithmetic='exact Python integers', bit_limit=200_000, censored_inputs=0)
    save_json(ROOT/'computations'/'egyptian-lengths.json', {
        'date': '2026-11-03', 'definition': 'uniform numerator k=1,...,N-1; distinct unit fractions',
        'optimal_method': 'complete search by increasing length, finite first-denominator bound; exact certificates',
        'greedy_bit_limit': 200_000, 'censored_inputs': 0, 'optimal': optimal, 'greedy': census})


def fixed_four(preview):
    limit = 1_000_000
    exceptional = four_exception_sieve(limit)
    exceptional[:5] = False
    counts = np.cumsum(exceptional)
    ns = np.unique(np.logspace(1, 6, 200).astype(int))
    fig, ax = panels()
    x = np.arange(1, 5)
    ax[0].bar(x-.18, [1/4, 3/4, 0, 0], width=.36, color=TEAL, label='minimum')
    ax[0].bar(x+.18, [1/4, 7/12, 1/12, 1/12], width=.36, color=ORANGE, label='greedy')
    ax[0].set(xticks=x, xlabel='number of unit fractions', ylabel='limiting probability',
              ylim=(0,.85), title='Uniform large denominators')
    ax[0].legend()
    ax[1].semilogx(ns, counts[ns]/(ns-4), color=PURPLE)
    ax[1].set(xlabel='cutoff N', ylabel=r'proportion of $5\leq n\leq N$',
              title='Minimum length at least three', ylim=(0,.23))
    finish(fig, '16-four-over-n', preview, sieve_limit=limit,
           exact_criterion='n odd, all prime divisors 1 mod 4', plot_cutoffs=ns.tolist())
    save_json(ROOT/'computations'/'four-over-n.json', {
        'method': 'exact integer sieve; excludes n <= 4',
        'greedy_period': 24, 'greedy_period_counts': [6,14,2,2],
        'optimal_limiting_probabilities': ['1/4','3/4','0','0'],
        'exception_counts': {int(n): int(counts[n]) for n in ns}})


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--preview-dir', type=Path, default=Path('/tmp/mathlab-process-nov3/figure-previews'))
    args = parser.parse_args()
    args.preview_dir.mkdir(parents=True, exist_ok=True)
    adsorption(args.preview_dir)
    egyptian(args.preview_dir)
    fixed_four(args.preview_dir)
    RECORD['data_sha256'] = {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                             for p in sorted((ROOT/'computations').glob('*.json'))}
    RECORD['code_sha256'] = {name: hashlib.sha256((ROOT/'code'/name).read_bytes()).hexdigest()
        for name in ('make_extension_figures.py', 'extension_calculations.py', 'egyptian_fractions.py')}
    save_json(ROOT/'figures'/'extensions-manifest.json', RECORD)
    print('Built and recorded all three extension figures.', flush=True)
