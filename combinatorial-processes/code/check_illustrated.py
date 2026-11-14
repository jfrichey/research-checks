"""Checks for the illustrated revision's added calculations; writes no files."""
from collections import Counter
from fractions import Fraction as F
from itertools import permutations, product
from math import comb, factorial, gcd, lcm
import mpmath as mp
import numpy as np

def dish_balance():
    for p in [F(3,5),F(2,3),F(4,5)]:
        pi=lambda j: (2*p-1)/(2*p) if j==0 else (2*p-1)/(2*p*p)*((1-p)/p)**(j-1)
        assert pi(0)==p*pi(1)
        for j in range(1,20): assert (1-p)*pi(j)==p*pi(j+1)
        assert pi(0)+pi(1)/(1-(1-p)/p)==1
    print('PASS: reflected dish-chain invariant probabilities')

def permutation_masses():
    for n in range(1,9):
        actual=Counter(sum(b==a+1 for a,b in zip(order,order[1:])) for order in permutations(range(n)))
        for j in range(n):
            predicted=sum((1-F(r,n))*F((-1)**(r-j),factorial(r-j)) for r in range(j,n))/factorial(j)
            assert predicted==F(actual[j],factorial(n))
    print('PASS: full permutation mass functions against enumeration through n=8')

def circle_targets():
    for theta in [F(1,5),F(1,2),F(3,4)]:
        for d in [F(0),F(1,10),F(3,10),F(1,2)]:
            cuts=sorted({F(0),F(1),(-theta/2)%1,(theta/2)%1,(d-theta/2)%1,(d+theta/2)%1})
            categories=Counter()
            for a,b in zip(cuts,cuts[1:]):
                x=(a+b)/2
                hit=[]
                for center in [F(0),d]:
                    delta=abs(x-center)
                    hit.append(min(delta,1-delta)<theta/2)
                categories[tuple(hit)]+=b-a
            both=categories[True,True]
            only=categories[True,False]+categories[False,True]
            none=categories[False,False]
            from_recursion=(1+only/theta)/(1-none)
            u=min(F(1),theta+min(d,theta))
            assert from_recursion==2/theta-1/u
    print('PASS: two-target coverage means from exact arc partitions')

def arrow_drift():
    for k in range(2,10):
        for signs in product([0,1],repeat=k):
            edges=[i for i in range(k) if signs[i]!=signs[(i+1)%k]]
            z=len(edges)
            if not z: continue
            singleton=sum(signs[(i-1)%k]!=signs[i]!=signs[(i+1)%k] for i in range(k))
            change=0
            for i in edges:
                new=signs[:i]+signs[i+1:]
                zz=sum(new[j]!=new[(j+1)%len(new)] for j in range(len(new)))
                change+=zz-z
            assert F(change,z)==F(-2*singleton,z)
    print('PASS: singleton-run boundary drift on every circle through k=9')

def pursuit():
    for k in [3,5,9]:
        theta=np.pi/k
        for t in [.01,.3,.7]:
            def points(t):
                r=1-t*np.sin(theta)
                phase=np.cos(theta)/np.sin(theta)*np.log(1/r)
                angles=2*np.pi*np.arange(k)/k+phase
                return r*np.exp(1j*angles)
            x=points(t)
            target=np.roll(x,-1)-x
            target/=np.abs(target)
            velocity=(points(t+1e-6)-points(t-1e-6))/(2e-6)
            assert np.max(np.abs(velocity-target))<1e-8
    print('PASS: exact regular-polygon trajectories satisfy unit-speed pursuit')

def geometric_joint():
    mp.mp.dps=45
    for p in [mp.mpf('.5'),mp.mpf('.8')]:
        q=1-p;a=1/q
        for phase in [mp.mpf(0),mp.mpf('.3'),mp.mpf('.7')]:
            mass=mp.mpf(0);mean=mp.mpf(0)
            for j in range(-15,200):
                A=a**(phase-j)
                mass+=mp.exp(-q*A)-mp.exp(-A)
                mean+=p*A*mp.exp(-q*A)
            H=p/q*sum(a**(phase-j)*mp.exp(-a**(phase-j)) for j in range(-15,200))
            assert abs(mass-1)<mp.mpf('1e-35')
            assert abs(mean-H)<mp.mpf('1e-35')
    n=2**20;m=20;p=mp.mpf('.5')
    for j in range(-2,4):
        pk=p**(m+j+1);lower=1-p**(m+j)
        A=mp.mpf(2)**(-j)
        for b in range(1,5):
            finite=mp.mpf(comb(n,b))*pk**b*lower**(n-b)
            limit=mp.exp(-A)*(p*A)**b/mp.factorial(b)
            assert abs(finite-limit)<mp.mpf('2e-6')
    print('PASS: joint geometric height/tie limits normalize and match finite samples')

def gumbel_covariance():
    mp.mp.dps=40
    for h in map(mp.mpf,['.1','.5','1','2','5']):
        a=mp.expm1(h)
        # Conditional expectation under the independent-exponential coupling.
        direct=mp.quad(lambda x:mp.exp(-x)*(-mp.log(x)-mp.euler)*mp.e1(a*x),[0,1,mp.inf])
        integral=mp.quad(lambda x:x/mp.expm1(x),[h,mp.inf])
        closed=mp.polylog(2,mp.exp(-h))-h*mp.log1p(-mp.exp(-h))
        assert abs(direct-integral)<mp.mpf('1e-30')
        assert abs(integral-closed)<mp.mpf('1e-30')
    assert abs(mp.quad(lambda x:2*x*x/mp.expm1(x),[0,1,mp.inf])-4*mp.zeta(3))<mp.mpf('1e-30')
    print('PASS: Gumbel covariance by independent conditional integration; 4*zeta(3) coefficient')

def divisor_pairs():
    D=8;N=lcm(*range(1,D+1))
    actual=sum(sum(F(1,d) for d in range(1,D+1) if n%d==0)**2 for n in range(1,N+1))/N
    formula=sum(F(gcd(a,b),a*a*b*b) for a in range(1,D+1) for b in range(1,D+1))
    assert actual==formula
    print('PASS: divisor-pair second moment from a complete residue period')

if __name__=='__main__':
    for check in [dish_balance,permutation_masses,circle_targets,arrow_drift,pursuit,geometric_joint,gumbel_covariance,divisor_pairs]:
        check()
    print('All added-calculation checks passed.')

