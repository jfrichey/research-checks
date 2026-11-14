"""Compact exact audit of the four RN25 boundary coefficient identities.

Recorded verification status, 2026-10-16: only checks 1-2 were collected
from the fresh agent execution. Checks 3-4 remain uncollected; no full
passing run is asserted. See RN section 26.

Run:
  python3 -B \
      research/arw/monotonicity/verify_three_site_boundary.py

Standalone: no repository imports, cached polynomials, numerical samples,
or generated files. Requires SymPy. All four checks reconstruct RN25's
compact matrices; H is defined by equation (122), not by an artifact.

Indexing (after tracing onto the three sites and removing self-loops):
  0 = occupation target v, 1 = added/missing site c, 2 = other root a.
  Lower input = {0,2}; upper input = {0,1,2}.
  r_i = (1-p_i)/p_i, so p_1 -> 0 means r_1 -> infinity.
  K = [[0,a,b],[c,0,d],[e,f,0]], k_i = 1-sum_j K_ij.
  d_i = det(I-K_{V\\{i}}), delta = det(I-K).

Important notation distinction:
  H_boundary = [r_1^3] H.
  RN25's endpoint H_0,H_4 instead mean [r_2^0]H,[r_2^4]H.
  Our fourth identity is [r_1^3 r_2^0]H.
"""

import time
import sympy as s


def data():
    a, b, c, d, e, f = edges = s.symbols("a b c d e f")
    r0, r1, r2 = r = s.symbols("r0 r1 r2")
    K = s.Matrix([[0, a, b], [c, 0, d], [e, f, 0]])
    k = [1 - sum(K[i, j] for j in range(3)) for i in range(3)]
    ds = [1 - d*f, 1 - b*e, 1 - a*c]
    delta = s.expand((s.eye(3) - K).det())
    Pd = s.prod(ds)
    Z = []
    N = s.zeros(3)
    bbar = []
    for i in range(3):
        j, l = [z for z in range(3) if z != i]
        Z.append(1 + ds[i]*(r[j] + r[l] + r[j]*r[l]))
        N[i, i] = Z[i]
        for j in range(3):
            if i != j:
                l = 3-i-j
                N[i, j] = -r[i]*(K[i, j]*(1+ds[i]*r[l])
                                  + K[i, l]*K[l, j])
        j, l = [z for z in range(3) if z != i]
        bbar.append(r[i]*(delta*Pd/ds[i]
                     + Pd*ds[i]*(k[i]*(r[j]+r[l]+r[j]*r[l])
                         + r[l]*K[i, j]*k[j]+r[j]*K[i, l]*k[l])))
    Zbar = Pd + delta*(Pd*(r0*r1+r0*r2+r1*r2+r0*r1*r2)
                      + sum(r[i]*Pd/ds[i] for i in range(3)))
    M = s.eye(3) + s.diag(*r)*(s.eye(3)-K)
    D_one = s.expand(M.det())
    u = M.adjugate()[:, 0]
    L = s.zeros(3)
    cv = s.zeros(3, 1)
    for i in range(3):
        j, l = [z for z in range(3) if z != i]
        L[i, i] = r[j]*r[l]*(K[j, i]*k[l]+K[l, i]*k[j]
                                       + K[j, i]*K[l, i]*k[i])
        cv[i] = r[j]*r[l]*(r[i]*k[j]*k[l]
                           + (k[j]+K[j, i]*k[i])*(k[l]+K[l, i]*k[i]))
        for j in range(3):
            if i != j:
                l = 3-i-j
                L[i, j] = -r[l]*((1+r[j])*k[l]+K[l, j]*k[j])
    return dict(edges=edges, r=r, K=K, k=k, ds=ds, delta=delta,
                Pd=Pd, Z=Z, N=N, bbar=s.Matrix(bbar), Zbar=Zbar,
                D_one=D_one, u=u, L=L, cv=cv)


def rn122_rhs(D, substitutions):
    """Build r0*d1*H from RN25 (122), specializing before multiplication."""
    sub = lambda obj: obj.subs(substitutions, simultaneous=True)
    N = sub(D['N'])
    u = sub(D['u'])
    done = sub(D['D_one'])
    xi = done*s.Matrix([0, 1, 1]) - sub(D['L'])*u
    source = sub(D['bbar']) - sub(D['Zbar'])*s.Matrix([0, 1, 0])
    return (sub(D['Pd'])*N.det()*(done-(sub(D['cv']).T*u)[0])
            + (xi.T*N.adjugate()*source)[0])


def coeff(expr, variable, power):
    return s.Poly(s.expand(expr), variable).coeff_monomial(variable**power)


def assert_zero(expr):
    assert s.cancel(expr) == 0


def four_coefficient_audit(D):
    r0, r1, r2 = D['r']
    d0, d1, d2 = D['ds']
    delta = D['delta']
    K, k = D['K'], D['k']
    Dstar = 1+d2*r0+d0*r2+delta*r0*r2
    Lam = 1+r0+r2+d1*r0*r2
    Ztrace = 1+delta*(r0/d0+r2/d2+r0*r2)
    assert_zero(coeff(D['D_one'], r1, 1)-Dstar)
    print('PASS 1: [r1] D_one = 1+d2*r0+d0*r2+delta*r0*r2', flush=True)
    assert_zero(coeff(D['Zbar'], r1, 1)-d0*d2*delta*D['Z'][1])
    print('PASS 2: [r1] Zbar = d0*d2*delta*Z1', flush=True)
    assert_zero(coeff(D['N'].det(), r1, 2)-d0*d2*Lam*Ztrace)
    print('PASS 3: [r1^2] detN = d0*d2*Lam*Ztrace', flush=True)
    g01 = K[0, 1]+K[0, 2]*K[2, 1]
    g10 = K[1, 0]+K[1, 2]*K[2, 0]
    ell1 = k[1]+K[1, 2]*k[2]
    expected = d2**2*g01*g10*ell1*(1+r0)*(1+d2*r0)
    rhs_at_r2_zero = rn122_rhs(D, {r2: 0})
    assert_zero(coeff(rhs_at_r2_zero, r1, 3)-r0*d1*expected)
    assert s.Poly(s.expand(rhs_at_r2_zero), r1).degree() == 3
    print('PASS 4: [r1^3 r2^0] H = d2^2*g01*g10*ell1*(1+r0)*(1+d2*r0)',
          flush=True)


if __name__ == '__main__':
    started = time.monotonic()
    four_coefficient_audit(data())
    print('Four exact identities passed in %.2fs.' % (time.monotonic()-started))
