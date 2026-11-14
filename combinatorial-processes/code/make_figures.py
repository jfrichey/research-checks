"""Rebuild the illustrated notebook's figures from exact formulas and simulations.

Run from any directory. Vector PDFs go to ../figures; optional PNG previews and
the machine-readable numerical record go to --preview-dir (default /tmp).
No imported image is modified. Each topic has its own recorded NumPy seed.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import itertools
import json
import math
import os
from pathlib import Path

os.environ.setdefault('MPLCONFIGDIR', '/tmp/mathlab-combinatorial-matplotlib')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Rectangle
import mpmath as mp
import numpy as np
from scipy import integrate, optimize, special, stats

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'figures'
TEAL, ORANGE, PURPLE, GRAY = '#167D8D', '#D66B36', '#7662A5', '#64717A'
COLORS = [TEAL, ORANGE, PURPLE, GRAY]
SEEDS = {i: 2026110200 + i for i in range(1, 16)}
RECORD = {'seeds': SEEDS, 'software': {'numpy': np.__version__, 'matplotlib': matplotlib.__version__}, 'figures': {}}
plt.rcParams.update({
    'font.family': 'serif', 'font.serif': ['STIXGeneral'],
    'mathtext.fontset': 'stix', 'font.size': 10,
    'axes.titlesize': 11, 'axes.labelsize': 10,
    'xtick.labelsize': 9, 'ytick.labelsize': 9,
    'axes.spines.top': False, 'axes.spines.right': False,
    'axes.edgecolor': GRAY, 'axes.linewidth': .65,
    'lines.linewidth': 1.6, 'legend.frameon': False,
    'legend.fontsize': 9, 'pdf.fonttype': 42,
    'savefig.facecolor': 'white', 'figure.facecolor': 'white',
})


def panels(n=2, height=2.9):
    return plt.subplots(1, n, figsize=(7.2, height), layout='constrained')


def finish(fig, name, preview, **params):
    pdf = OUT / f'{name}.pdf'
    fig.savefig(pdf, metadata={'Creator': 'make_figures.py', 'Author': 'ChatGPT',
        'CreationDate': datetime(2026, 11, 2, tzinfo=timezone.utc),
        'ModDate': datetime(2026, 11, 2, tzinfo=timezone.utc)})
    fig.savefig(preview / f'{name}.png', dpi=180)
    plt.close(fig)
    RECORD['figures'][name] = params | {'sha256': hashlib.sha256(pdf.read_bytes()).hexdigest()}


def socks(preview):
    rng = np.random.default_rng(SEEDS[1])
    n = 600
    draw = rng.permutation(np.repeat(np.arange(n), 2))
    seen = np.zeros(n, dtype=int)
    path = [0]
    for sock in draw:
        seen[sock] += 1
        path.append(path[-1] + (1 if seen[sock] == 1 else -1))
    k = np.arange(2*n+1)
    t = k/(2*n)
    mean = k - k*(k-1)/(2*n-1)
    centered = (np.asarray(path)-mean)/(2*np.sqrt(n))
    fig, ax = panels()
    ax[0].plot(t, path, color=TEAL, lw=1, label='sample')
    ax[0].plot(t, mean, color=ORANGE, label='exact mean')
    ax[0].set(xlabel='fraction drawn', ylabel='unmatched socks', title='A pile builds, then disappears')
    ax[0].legend(loc='upper left')
    ax[1].fill_between(t, -2*t*(1-t), 2*t*(1-t), color=TEAL, alpha=.10,
                       label=r'$\pm 2$ limiting standard deviations')
    ax[1].plot(t, centered, color=TEAL, lw=.9)
    ax[1].axhline(0, color=GRAY, lw=.6)
    ax[1].set(xlabel='fraction drawn', ylabel=r'$(X_k-\mathbb{E} X_k)/(2\sqrt{n})$', title='Fluctuations remember the endpoint')
    finish(fig, '01-socks', preview, n=n, seed=SEEDS[1], band='pointwise asymptotic scale, not simultaneous coverage')


def dishes(preview):
    rng = np.random.default_rng(SEEDS[2])
    n = 450
    fig, ax = panels(3)
    models = [('Empty the sink', .5, True), ('Clean one: p = 1/2', .5, False),
              ('Clean one: p = 0.35', .35, False)]
    for a, (title, p, reset), color in zip(ax, models, COLORS):
        path = [0]
        for _ in range(n):
            x = path[-1]
            if reset:
                x = 0 if rng.random() < p else x+1
            else:
                x = 1 if x == 0 else x+(-1 if rng.random()<p else 1)
            path.append(x)
        a.plot(np.arange(n+1), path, color=color, lw=1)
        a.set(title=title, xlabel='step', ylabel='dishes')
    finish(fig, '02-dishes', preview, steps=n, seed=SEEDS[2], models=models)


def coupons(preview):
    rng = np.random.default_rng(SEEDS[3])
    n = 200
    draws = rng.integers(n, size=1600)
    seen, path = set(), [0]
    for x in draws:
        seen.add(int(x)); path.append(len(seen))
    k = np.arange(len(path))
    fig, ax = panels()
    ax[0].plot(k/n, np.asarray(path)/n, color=TEAL, lw=1, label='sample')
    ax[0].plot(k/n, 1-(1-1/n)**k, color=ORANGE, label='exact mean')
    ax[0].set(xlabel=r'draws / $n$', ylabel='fraction collected', title='Most types arrive early', ylim=(0,1.04))
    ax[0].legend(loc='lower right')
    n2, reps = 1000, 6000
    times = rng.geometric(np.arange(1,n2+1)/n2, size=(reps,n2)).sum(axis=1)
    values = times/n2-np.log(n2)
    ax[1].hist(values, bins=np.linspace(-2.5,7,45), density=True, color=TEAL, alpha=.35,
               label='simulated collection times')
    x = np.linspace(-2.5,7,400)
    ax[1].plot(x, np.exp(-x-np.exp(-x)), color=ORANGE, label='Gumbel density')
    ax[1].set(xlabel=r'$\tau/n-\log n$', ylabel='density', title='The last few set the clock')
    ax[1].legend(fontsize=8)
    finish(fig, '03-coupons', preview, path_n=n, completion_n=n2, repetitions=reps, seed=SEEDS[3])


def greedy(priorities):
    taken = np.zeros(len(priorities), dtype=bool)
    for j in np.argsort(priorities):
        if (j == 0 or not taken[j-1]) and (j+1 == len(taken) or not taken[j+1]):
            taken[j] = True
    return taken


def seating(preview):
    rng = np.random.default_rng(SEEDS[4])
    fig, ax = panels()
    n = 24
    w = rng.random(n)
    occupied = greedy(w)
    ax[0].vlines(np.arange(n), 0, w, color=GRAY, lw=.8, alpha=.65)
    ax[0].scatter(np.arange(n), w, c=np.where(occupied, TEAL, '#B7BEC3'), s=32, zorder=3)
    ax[0].scatter(np.where(occupied)[0], np.full(occupied.sum(),-.09), color=TEAL, s=38, marker='s')
    ax[0].set(xlabel='seat', ylabel='arrival priority', ylim=(-.17,1.05), title='Earlier arrivals block their neighbors')
    n2, reps = 160, 5000
    counts = np.array([greedy(rng.random(n2)).sum() for _ in range(reps)])
    mean = .5*(1-np.exp(-2))*n2 + .5*(1-3*np.exp(-2))
    variance = np.exp(-4)*(n2+3)
    ax[1].hist(counts, bins=np.arange(counts.min()-.5,counts.max()+1.5), density=True,
               color=TEAL, alpha=.4)
    x = np.linspace(counts.min()-1, counts.max()+1, 250)
    ax[1].plot(x, stats.norm.pdf(x,mean,np.sqrt(variance)), color=ORANGE, label='normal approximation')
    ax[1].set(xlabel='occupied seats', ylabel='probability per seat count', title='A narrow range around 43.2%')
    ax[1].legend(fontsize=8)
    finish(fig, '04-seating', preview, diagram_n=n, histogram_n=n2, repetitions=reps, seed=SEEDS[4])


def permutations(preview):
    fig, ax = panels()
    blocks = [[3,4,5],[8],[1,2],[6,7]]
    x = 0
    for block, color in zip(blocks, COLORS):
        ax[0].add_patch(Rectangle((x-.33,.5),len(block)-.14,.55,facecolor=color,alpha=.15,edgecolor=color))
        for v in block:
            ax[0].text(x,.77,str(v),ha='center',va='center',fontsize=17,color=color); x += 1
        x += .25
    ax[0].text(3.6,.1,'Four prescribed successions form four blocks',ha='center',fontsize=10)
    ax[0].set(xlim=(-.6,8.7),ylim=(-.15,1.5), title='Glue each required succession into a block')
    ax[0].axis('off')
    support = np.arange(7)
    for n,color,marker in [(6,TEAL,'o'),(20,PURPLE,'s')]:
        probabilities = [sum((1-r/n)*(-1)**(r-j)/math.factorial(r-j)
                             for r in range(j,n))/math.factorial(j) if j<n else 0 for j in support]
        ax[1].plot(support, probabilities, marker=marker, color=color, label=f'n = {n}', lw=1)
    ax[1].plot(support, stats.poisson.pmf(support,1), color=ORANGE, ls='--', label='Poisson(1)')
    ax[1].set(xlabel='number of successions',ylabel='probability',title='The Poisson limit appears quickly')
    ax[1].legend()
    finish(fig, '05-permutations', preview, exact_pmf_n=[6,20], illustration=[v for b in blocks for v in b])


def uniform_mean(x):
    return math.fsum((-1)**j * (x-j)**j / math.factorial(j) * math.exp(x-j)
                     for j in range(math.floor(x)+1))


def sums(preview):
    x = np.linspace(.001,8,700)
    mean = np.array([uniform_mean(float(u)) for u in x])
    variance = mean*(2*x+1-mean)
    fig, ax = panels()
    ax[0].plot(x, mean-2*x, color=TEAL)
    ax[0].axhline(2/3,color=ORANGE,ls='--',label=r'$2/3$')
    ax[0].set(xlabel='threshold',ylabel=r'$\mathbb{E}\tau_x-2x$', title='The overshoot leaves a constant correction')
    ax[0].legend()
    ax[1].plot(x,variance,color=TEAL,label='exact variance')
    ax[1].plot(x,2*x/3,color=ORANGE,ls='--',label=r'$2x/3$')
    ax[1].set(xlabel='threshold',ylabel='variance',title='One mean formula gives the second moment')
    ax[1].legend()
    finish(fig, '06-sums', preview, method='finite exponential sums', threshold_range=[.001,8])


def pursuit(preview):
    fig, ax = panels()
    k, radius = 5, 1.
    radii = np.geomspace(1,.003,800)
    theta = np.pi/k
    phase = np.cos(theta)/np.sin(theta)*np.log(1/radii)
    for j in range(k):
        angle = phase+2*np.pi*j/k
        ax[0].plot(radii*np.cos(angle),radii*np.sin(angle),color=COLORS[j%4],lw=1.1)
        ax[0].scatter([np.cos(2*np.pi*j/k)],[np.sin(2*np.pi*j/k)],color=COLORS[j%4],s=25)
    ax[0].set(title='Euclidean pursuit: a shrinking pentagon',aspect='equal',xlim=(-1.1,1.1),ylim=(-1.1,1.1))
    ax[0].axis('off')
    ax[1].add_patch(Circle((0,0),1,fill=False,ec='#B7BEC3',lw=2))
    positions = np.array([0,.32,.64])*2*np.pi
    for j,angle in enumerate(positions):
        ax[1].scatter([np.cos(angle)],[np.sin(angle)],color=TEAL,s=45)
        ax[1].annotate('',xy=(np.cos(angle+.22),np.sin(angle+.22)),xytext=(np.cos(angle),np.sin(angle)),
                        arrowprops={'arrowstyle':'->','color':TEAL,'lw':1.8})
    ax[1].text(0,-1.28,'All three shortest displacements point the same way',ha='center',fontsize=9)
    ax[1].set(title='Circle pursuit: a translating configuration',aspect='equal',xlim=(-1.25,1.25),ylim=(-1.45,1.15))
    ax[1].axis('off')
    finish(fig, '07-pursuit', preview, regular_polygon=k, exact_formula=True, circle_positions=[0,.32,.64])


def arrow_path(n, rng):
    # Deliberately direct implementation for figures; no dependency on the
    # separately maintained fast arrow-process implementation.
    signs = rng.integers(0,2,n).tolist()
    runs = [1]
    for j in range(1,n):
        if signs[j] == signs[j-1]: runs[-1] += 1
        else: runs.append(1)
    if len(runs)>1 and signs[0]==signs[-1]: runs[0] += runs.pop()
    counts = [0 if len(runs)==1 else len(runs)]
    while len(runs)>1:
        j = int(rng.integers(len(runs)))
        runs[j] -= 1
        if runs[j] == 0:
            if len(runs)==2:
                runs = [sum(runs)]
            else:
                runs = runs[j:]+runs[:j]
                runs = [runs[-1]+runs[1]]+runs[2:-1]
        counts.append(0 if len(runs)==1 else len(runs))
    return np.array(counts), runs[0]


def arrows(preview):
    rng = np.random.default_rng(SEEDS[7])
    fig, ax = panels()
    ns = [1000,4000,16000]
    for n,color in zip(ns,COLORS):
        counts,R = arrow_path(n,rng)
        assert len(counts)-1 == n-R
        assert counts[-1]==0 and np.all(np.isin(np.diff(counts),[-2,0]))
        ax[0].plot(np.arange(len(counts))/n,counts/n,color=color,lw=1,label=f'k = {n:,}')
    ax[0].set(xlabel=r'deletions / $k$',ylabel=r'disagreements / $k$',title='A candidate macroscopic curve')
    ax[0].legend(fontsize=8)
    k,reps = 2500,15000
    plus = rng.binomial(k,.5,size=reps)
    minus = k-plus
    active = (plus>0)&(minus>0)
    while active.any():
        ids = np.flatnonzero(active)
        coin = rng.integers(0,2,len(ids))
        plus[ids] -= coin
        minus[ids] -= 1-coin
        active = (plus>0)&(minus>0)
    R = plus+minus
    ax[1].hist(R/np.sqrt(k),bins=np.linspace(0,5,38),density=True,color=TEAL,alpha=.35)
    x=np.linspace(0,5,300)
    ax[1].plot(x,np.exp(-x*x/4)/np.sqrt(np.pi),color=ORANGE,label=r'$|N(0,2)|$')
    ax[1].set(xlabel=r'survivors / $\sqrt{k}$',ylabel='density',title='The terminal count has an explicit limit')
    ax[1].legend()
    finish(fig, '07-arrows', preview, path_sizes=ns, terminal_k=k, repetitions=reps, seed=SEEDS[7],
           terminal_mean=float(R.mean()), second_identity_estimate=float(np.mean(R.astype(float)**2+R)),
           second_identity_exact=2*k)


def coverage(preview):
    rng=np.random.default_rng(SEEDS[8])
    fig,ax=panels()
    theta=.16
    centers=np.sort(rng.random(12))
    angles=np.linspace(0,1,3000,endpoint=False)
    distance=np.abs(angles[:,None]-centers)
    covered=(np.minimum(distance,1-distance)<=theta/2).any(axis=1)
    ax[0].plot(np.cos(2*np.pi*angles),np.sin(2*np.pi*angles),color='#DCE2E5',lw=8)
    xx=np.where(covered,np.cos(2*np.pi*angles),np.nan)
    yy=np.where(covered,np.sin(2*np.pi*angles),np.nan)
    ax[0].plot(xx,yy,color=TEAL,lw=8,solid_capstyle='butt')
    ax[0].scatter(1.12*np.cos(2*np.pi*centers),1.12*np.sin(2*np.pi*centers),color=ORANGE,s=15)
    ax[0].set(title='Coverage is decided by the largest spacing',aspect='equal',xlim=(-1.25,1.25),ylim=(-1.25,1.25))
    ax[0].axis('off')
    reps=5000
    for n,color in zip([20,100,500],COLORS):
        e=rng.exponential(size=(reps,n))
        maxima=(e/e.sum(axis=1)[:,None]).max(axis=1)
        values=np.sort(n*maxima-np.log(n))
        ax[1].plot(values,np.arange(1,reps+1)/reps,color=color,lw=1,label=f'n = {n}')
    x=np.linspace(-2.5,7,300)
    ax[1].plot(x,np.exp(-np.exp(-x)),color=GRAY,ls='--',label='Gumbel CDF')
    ax[1].set(xlabel=r'$n\max_i D_i-\log n$',ylabel='cumulative probability',xlim=(-2.5,6),title='A second route to the Gumbel law')
    ax[1].legend(fontsize=8,loc='lower right')
    finish(fig,'08-coverage',preview,arc_length=theta,centers=centers.tolist(),spacing_n=[20,100,500],repetitions=reps,seed=SEEDS[8])


def inverse(preview):
    fig,ax=panels()
    path=np.array([0,1,1,2,5,5,5,6,8,8,10])
    time=np.arange(len(path))
    ax[0].step(time,path,where='post',color=TEAL)
    ax[0].axhline(4,color=ORANGE,ls='--',lw=1)
    ax[0].vlines(4,4,5,color=PURPLE,lw=4)
    ax[0].scatter([4],[5],color=PURPLE,s=30)
    ax[0].annotate('overshoot',xy=(4,4.6),xytext=(5.2,3),arrowprops={'arrowstyle':'->','color':GRAY})
    ax[0].annotate('one plateau, three observation times',xy=(5.5,5),xytext=(.7,9),
                    arrowprops={'arrowstyle':'->','color':GRAY},fontsize=9)
    ax[0].set(xlabel='time',ylabel='height',title='A nondecreasing sample path')
    levels=np.linspace(.01,10,1000)
    hitting=np.searchsorted(path,levels,side='left')
    ax[1].plot(levels,hitting,color=ORANGE)
    ax[1].scatter([4],[4],color=PURPLE,s=30)
    ax[1].set(xlabel='target height',ylabel='first hitting time',title='Its generalized inverse')
    finish(fig,'09-inverse',preview,path=path.tolist(),threshold=4)


def maxima(preview):
    mp.mp.dps=40
    phases=np.linspace(0,1,301)
    a=mp.mpf(2)
    values=np.array([float(sum(a**(mp.mpf(str(r))-j)*mp.exp(-a**(mp.mpf(str(r))-j))
                               for j in range(-12,140))-1/mp.log(2)) for r in phases])
    fig,ax=panels()
    ax[0].plot(phases,values,color=TEAL)
    ax[0].axhline(0,color=GRAY,lw=.6)
    ax[0].ticklabel_format(axis='y',style='sci',scilimits=(0,0))
    ax[0].set(xlabel=r'phase $\{\log_2 n\}$',ylabel=r'$H(r)-1/\log 2$',title='A tiny oscillation that never vanishes')
    lam=np.linspace(0,7,500)
    j=np.arange(100)
    expectation=np.array([100*np.sum(stats.poisson.pmf(j,l)*stats.poisson.cdf(j,l)**99) for l in lam])
    ax[1].plot(lam,expectation,color=ORANGE)
    ax[1].set(xlabel=r'Poisson mean $\lambda$',ylabel='expected number tied at the maximum',yscale='log',title='Discrete maxima retain multiplicities')
    finish(fig,'10-maxima',preview,geometric_p=.5,periodic_sum_j=[-12,139],poisson_n=100,poisson_sum_j=[0,99],
           periodic_range=[float(values.min()),float(values.max())])


def gumbel_cov(h):
    if h==0: return np.pi**2/6
    r=np.exp(-h)
    return special.spence(1-r)-h*np.log1p(-r)


def gumbel(preview):
    rng=np.random.default_rng(SEEDS[11])
    fig,ax=panels()
    t,u=0.,-np.log(rng.exponential())
    ts,us=[t],[u]
    while t<24:
        wait=np.log1p(np.exp(u)*rng.exponential())
        if t+wait>24:
            ts.append(24);us.append(u-(24-t));break
        t+=wait;u-=wait;ts.append(t);us.append(u)
        u+=rng.exponential();ts.append(t);us.append(u)
    ax[0].plot(ts,us,color=TEAL,lw=1.2)
    ax[0].axhline(np.euler_gamma,color=GRAY,ls='--',lw=.8,label='stationary mean')
    ax[0].set(xlabel='logarithmic time',ylabel=r'$U(t)$',title='Linear drift, exponential upward jumps')
    ax[0].legend(fontsize=8)
    h=np.linspace(0,7,300)
    ax[1].plot(h,[gumbel_cov(float(v))/(np.pi**2/6) for v in h],color=ORANGE,label='exact correlation')
    reps=150000
    lags=np.array([.25,.5,1,2,4,6])
    estimates,errors=[],[]
    for lag in lags:
        r=np.exp(-lag);v=rng.exponential(size=reps);e=rng.exponential(size=reps)
        u0=-np.log(v);uh=-np.log(np.minimum(v/r,e/(1-r)))
        product=(u0-np.euler_gamma)*(uh-np.euler_gamma)/(np.pi**2/6)
        estimates.append(product.mean());errors.append(1.96*product.std(ddof=1)/np.sqrt(reps))
    ax[1].errorbar(lags,estimates,yerr=errors,fmt='o',ms=3,color=TEAL,lw=.8,label='independent pairs; 95% error bars')
    ax[1].set(xlabel='time separation',ylabel='correlation',title='An exact measure of memory',ylim=(-.03,1.03))
    ax[1].legend(fontsize=8)
    finish(fig,'11-gumbel',preview,path_duration=24,pair_repetitions=reps,seed=SEEDS[11],lags=lags.tolist(),
           correlations=estimates,error_bars=errors)


def ruler(preview):
    fig,ax=panels()
    x=np.linspace(-1,1,600);f=1-np.abs(x);height=.55
    ax[0].plot(x,f,color=TEAL)
    ax[0].fill_between(x,0,f,where=f<=height,color=TEAL,alpha=.3)
    ax[0].axhline(height,color=ORANGE,ls='--',label=r'height $h$')
    ax[0].set(xlabel='location',ylabel=r'$f(y)$',title='Sample a location with density f',ylim=(0,1.1))
    ax[0].legend()
    z=np.linspace(0,1,400)
    ax[1].plot(z,2*z,color=TEAL)
    ax[1].fill_between(z,0,2*z,where=z<=height,color=TEAL,alpha=.3)
    ax[1].axvline(height,color=ORANGE,ls='--')
    ax[1].set(xlabel='density height',ylabel=r'$g(h)=2h$',title='The sampled height favors taller values',ylim=(0,2.15))
    finish(fig,'12-ruler',preview,density='f(y)=(1-|y|) on [-1,1]',height=height,shaded_probability=height**2)


def products(preview):
    rng=np.random.default_rng(SEEDS[13])
    fig,ax=panels()
    n=40
    logproducts=np.cumsum(np.log(rng.random(n)))
    ax[0].plot(np.arange(1,n+1),logproducts,color=TEAL,marker='.',ms=3,label='one product sequence')
    ax[0].plot(np.arange(1,n+1),-np.arange(1,n+1),color=ORANGE,ls='--',label=r'$-n$')
    ax[0].set(xlabel='number of factors',ylabel='log product',title='Multiplication becomes Poisson arrival times')
    ax[0].legend(fontsize=8)
    reps,terms=90000,60
    product=np.ones(reps);total=np.zeros(reps)
    for _ in range(terms):
        product*=rng.random(reps);total+=product
    ax[1].hist(total,bins=np.linspace(0,4.5,70),density=True,color=TEAL,alpha=.35)
    s=np.linspace(.001,2,300)
    density=np.exp(-np.euler_gamma)*np.where(s<=1,1,1-np.log(s))
    ax[1].plot(s,density,color=ORANGE,label='exact density on (0, 2)')
    ax[1].set(xlabel=r'$S=\sum_n Y_n$',ylabel='density',title='The sum has the Dickman distribution')
    ax[1].legend(fontsize=8)
    finish(fig,'13-products',preview,path_terms=n,sum_terms=terms,repetitions=reps,seed=SEEDS[13],
           expected_omitted_sum=2.**(-terms),mean=float(total.mean()),variance=float(total.var()))


def fishing(preview):
    w,fast,slow=.7,2.,.1
    def survival(t):return w*np.exp(-fast*t)+(1-w)*np.exp(-slow*t)
    def mean(t,c):return (c+w*(-np.expm1(-fast*t))/fast+(1-w)*(-np.expm1(-slow*t))/slow)/(1-survival(t))
    fig,ax=panels()
    t=np.linspace(.001,10,500)
    ax[0].plot(t,survival(t),color=TEAL,label='survival probability')
    ax[0].plot(t,w*np.exp(-fast*t),color=GRAY,ls=':',label='fast component')
    ax[0].plot(t,(1-w)*np.exp(-slow*t),color=PURPLE,ls='--',label='slow component')
    ax[0].set(xlabel='time without a bite',ylabel='probability',title='A mixture of fast and slow casts')
    ax[0].legend(fontsize=8)
    optima={}
    for c,color in zip([0,.2,1.],COLORS):
        ax[1].plot(t,mean(t,c),color=color,label=f'setup cost c = {c:g}')
        if c:
            result=optimize.minimize_scalar(lambda u:mean(u,c),bounds=(.001,100),method='bounded')
            optima[str(c)]={'timeout':float(result.x),'mean':float(result.fun)}
            ax[1].scatter(result.x,result.fun,color=color,s=25)
    ax[1].set(xlabel='timeout',ylabel='expected total catch time',title='Setup cost selects a finite timeout',ylim=(.4,7.5))
    ax[1].legend(fontsize=8)
    finish(fig,'14-fishing',preview,mixture_weight=w,rates=[fast,slow],costs=[0,.2,1.],optima=optima)


def divisors(preview):
    rng=np.random.default_rng(SEEDS[15])
    N=100000
    sigma=np.zeros(N+1)
    for d in range(1,N+1): sigma[d::d]+=1/d
    fig,ax=panels()
    ax[0].vlines(np.arange(1,181),1,sigma[1:181],color=TEAL,lw=.8)
    ax[0].axhline(2,color=ORANGE,ls='--',label='abundance threshold')
    ax[0].set(xlabel='integer',ylabel=r'$\sigma_1(n)/n$',title='Arithmetic spikes at individual integers')
    ax[0].legend(fontsize=8)
    for m,color in [(1000,PURPLE),(N,TEAL)]:
        values=np.sort(sigma[1:m+1])
        ax[1].plot(values,np.arange(1,m+1)/m,color=color,lw=1,label=f'uniform integer up to {m:,}')
    primes=[p for p in range(2,102) if all(p%d for d in range(2,math.isqrt(p)+1))]
    reps=100000
    product=np.ones(reps)
    for p in primes:
        valuation=rng.geometric(1-1/p,size=reps)-1
        product*=(1-np.power(float(p),-valuation-1))/(1-1/p)
    values=np.sort(product)
    ax[1].plot(values,np.arange(1,reps+1)/reps,color=ORANGE,ls='--',lw=1,label='independent primes up to 101')
    ax[1].set(xlabel='normalized divisor sum',ylabel='cumulative probability',title='A stable law after averaging',xlim=(1,4.3))
    ax[1].legend(fontsize=8,loc='lower right')
    finish(fig,'15-divisors',preview,integer_cutoffs=[1000,N],prime_cutoff=101,repetitions=reps,seed=SEEDS[15],
           mean_arithmetic=float(sigma[1:].mean()),mean_product=float(product.mean()),log_tail_bound=1/101)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--preview-dir',type=Path,default=Path('/tmp/mathlab-process-illustrated/figure-previews'))
    args=parser.parse_args()
    OUT.mkdir(exist_ok=True);args.preview_dir.mkdir(parents=True,exist_ok=True)
    for function in [socks,dishes,coupons,seating,permutations,sums,pursuit,arrows,coverage,inverse,maxima,gumbel,ruler,products,fishing,divisors]:
        function(args.preview_dir)
        print(f'Built {function.__name__}',flush=True)
    # The durable manifest contains scientific parameters and output hashes.
    # Previews are disposable renderings and stay in /tmp.
    (OUT/'manifest.json').write_text(json.dumps(RECORD,indent=2)+'\n')


if __name__=='__main__': main()
