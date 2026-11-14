"""Exact three-site rare-sleep field formula and minimal-support audit.

Requires SymPy (checked with 1.14.0). Uses small scalar polynomials and
walk Green matrices; no parameter sampling or output files.
AI-written research code. See response-renewal.md Section 22.
"""
import sympy as s
from collections import defaultdict
from itertools import combinations
from math import factorial

a,b,c,d,e,f=s.symbols('a b c d e f')
edges=(a,b,c,d,e,f)
K=s.Matrix([[0,a,b],[c,0,d],[e,f,0]])
k0,k1,k2=1-a-b,1-c-d,1-e-f
d0,d1,d2=1-d*f,1-b*e,1-a*c
det=1-a*c-b*e-d*f-a*d*e-b*c*f
H1,H2=k1+d*k2,k2+f*k1
T1,T2=a+b*f,b+a*d

# The slope is core/(d0^2*d1*d2*det). This polynomial is the cleared
# Green/avoidance formula, not a fitted or interpolated expression.
core=d1*d2*(det*d0*d0-k0*d0**3
    -a*(k1*k1*d0**3+H2*(H1*H1-k1*k1*d0*d0))
    -b*(k2*k2*d0**3+H1*(H2*H2-k2*k2*d0*d0)))
core-=d0*d1*(d2*k1+k0*c)*T2*H2*(d0-H2)
core-=d0*d2*(d1*k2+k0*e)*T1*H1*(d0-H1)
core=s.Poly(core,*edges)
slope=core.as_expr()/(d0*d0*d1*d2*det)

def row_transform(original,row,degree,sign=1):
    out=defaultdict(int)
    x,y=2*row,2*row+1
    for mon,coefficient in original.items():
        n=degree-mon[x]-mon[y]
        for i in range(n+1):
            for j in range(n-i+1):
                m=list(mon);m[x]+=i;m[y]+=j
                out[tuple(m)]+=coefficient*sign**(i+j)*factorial(n)//(
                    factorial(i)*factorial(j)*factorial(n-i-j))
    return {m:c for m,c in out.items() if c}

def support_audit():
    arcs=[(0,1),(0,2),(1,0),(1,2),(2,0),(2,1)]
    def strongly_connected(ids):
        reach={i:{i} for i in range(3)}
        for index in ids:
            x,y=arcs[index];reach[x].add(y)
        for _ in range(3):
            for i in range(3):reach[i]|=set().union(*(reach[j] for j in list(reach[i])))
        return all(len(v)==3 for v in reach.values())
    unsigned=0;exceptional=[]
    for size in (3,4):
        for ids in combinations(range(6),size):
            if not strongly_connected(ids):continue
            P=s.Poly(core.as_expr().subs({edges[i]:0 for i in range(6) if i not in ids}),*edges)
            terms={m:int(c) for m,c in P.terms()}
            degrees=[max(m[2*i]+m[2*i+1] for m in terms) for i in range(3)]
            for row in range(3):terms=row_transform(terms,row,degrees[row])
            if any(c<0 for c in terms.values()):exceptional.append(tuple(arcs[i] for i in ids))
            else:unsigned+=1
    assert unsigned==9
    assert exceptional==[((0,1),(1,0),(1,2),(2,0)),
                         ((0,2),(1,0),(2,0),(2,1))]
    print('VERIFIED two cycle supports and seven other four-edge supports have nonnegative slopes',flush=True)
    print('Exceptional supports (negative weights exist):',exceptional,flush=True)

def negative_support_audit():
    original={m:int(c) for m,c in core.terms()}
    terms=original
    degrees=(3,5,5)
    for row in range(3):terms=row_transform(terms,row,degrees[row])
    negative=[m for m,c in terms.items() if c<0]
    assert len(terms)==2700 and len(negative)==58
    for m in negative:
        sinks=[degrees[i]-m[2*i]-m[2*i+1] for i in range(3)]
        assert m[2]>0 and m[4]>0 and sinks[0]>=2
        assert ((m[0]>0 and m[3]>0 and sinks[2]>0)
                or (m[1]>0 and m[5]>0 and sinks[1]>0))
    restored=terms
    for row in (2,1,0):restored=row_transform(restored,row,degrees[row],-1)
    assert restored==original
    print('VERIFIED 58 negative terms out of 2700 all require the stated shortcut-cycle/sink motif',flush=True)

def family_audit():
    B,C,E,F=s.symbols('B C E F')
    U=1+B+E+F+B*F
    D=(1+C)*(1+B+E+F)+B*F
    J=B*(C+1)-F
    R=J*E*((C+1)*E+C*F)+(B+1)*(F+1)*(C+F+1)*(B*(C+1)*E+F*(B*(C+1)+1))
    expected=B*C*R/((1+B)*(1+C)**2*(1+E+F)*U*D)
    terms={m:int(c) for m,c in core.terms() if m[0]==m[3]==0}
    degrees=[max(m[2*i]+m[2*i+1] for m in terms) for i in range(3)]
    # Homogenize only the surviving edge in rows 0 and 1. Terms involving
    # an absent edge are set to zero immediately after each row transform.
    for row in range(3):
        terms=row_transform(terms,row,degrees[row])
        terms={m:c for m,c in terms.items() if m[0]==m[3]==0}
    rate=s.Poly.from_dict({(m[1],m[2],m[4],m[5]):c for m,c in terms.items()},(B,C,E,F))
    assert all(n>=3 for n in degrees)
    rhs=B*C*R*(1+B)**(degrees[0]-3)*(1+C)**(degrees[1]-3)*(1+E+F)**(degrees[2]-3)
    assert rate==s.Poly(rhs,B,C,E,F)
    simple=s.factor(expected.subs({B:1,C:1,F:3}))
    assert simple==s.factor((-2*E*E+77*E+360)/(8*(E+4)*(E+8)*(2*E+13)))
    assert simple.subs(E,60)==-s.Rational(555,1157632)
    print('VERIFIED four-edge family sign formula and witness slope:',simple.subs(E,60),flush=True)

if __name__=='__main__':
    assert len(core.terms())==231
    print('VERIFIED cleared rare-sleep numerator has 231 terms',flush=True)
    kap1,kap2=k1+c*k0,k2+e*k0
    assert s.Poly(det-d1*kap1-d2*kap2+kap1*kap2,*edges).is_zero
    # Leading y term: k-Gamma=k0(I-K)e0, so G(k-Gamma)=k0e0.
    rhs=s.Matrix([k0,k1-kap1,k2-kap2])
    assert rhs == k0*(s.eye(3)-K)[:,0]
    # Exact common-odds renewal split N=t^2 D(I-K)+t(3D-adj(I-K))+I.
    t=s.symbols('t'); D=s.diag(d0,d1,d2); minors=(d0,d1,d2)
    N=s.zeros(3)
    for i in range(3):
        N[i,i]=1+minors[i]*(2*t+t*t)
        for j in range(3):
            if i!=j:
                l=3-i-j
                N[i,j]=-t*(K[i,j]*(1+minors[i]*t)+K[i,l]*K[l,j])
    expected=t*t*D*(s.eye(3)-K)+t*(3*D-(s.eye(3)-K).adjugate())+s.eye(3)
    assert all(s.Poly(q,t,*edges).is_zero for q in N-expected)
    # Audit the Green identity behind the first correction to the renewal.
    adj=(s.eye(3)-K).adjugate()
    for i,j,minor in [(1,2,d2),(2,1,d1)]:
        assert s.Poly(adj[i,0]*adj[j,j]-adj[i,j]*adj[j,0]-det*K[i,0],*edges).is_zero
    print('VERIFIED killed-Green form of the two escape weights',flush=True)
    family_audit()
    support_audit()
    negative_support_audit()
