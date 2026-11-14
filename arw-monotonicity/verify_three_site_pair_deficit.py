"""Exact three-site pair-deficit certificate; see response-renewal.md Section 25.

Run with python3 -B and SymPy available (audited with 1.14.0).
Rebuilds and checks the saved coefficients; writes no files. AI-written.
All arithmetic is integer polynomial arithmetic. No parameter sampling.
"""
import sympy as s
import time
import gzip
import json
import hashlib
from pathlib import Path
from collections import defaultdict
from math import factorial

START=time.monotonic()
def say(*args): print(round(time.monotonic()-START,2),*args,flush=True)

def build():
    a,b,c,d,e,f=edges=s.symbols('a b c d e f')
    r0,r1,r2=r=s.symbols('r0 r1 r2')
    gens=edges+r
    K=s.Matrix([[0,a,b],[c,0,d],[e,f,0]])
    k=[1-a-b,1-c-d,1-e-f]
    ds=[1-d*f,1-b*e,1-a*c]
    detK=s.expand((s.eye(3)-K).det())
    Z=[1+ds[i]*(r[j]+r[l]+r[j]*r[l])
       for i in range(3) for j,l in [[z for z in range(3) if z!=i]]]
    N=s.zeros(3)
    for i in range(3):
        N[i,i]=Z[i]
        for j in range(3):
            if i!=j:
                l=3-i-j
                N[i,j]=-r[i]*(K[i,j]*(1+ds[i]*r[l])+K[i,l]*K[l,j])
    P=s.prod(ds)
    Zbar=s.expand(P+detK*(P*(r0*r1+r0*r2+r1*r2+r0*r1*r2)
                         +sum(r[i]*s.prod(ds[j] for j in range(3) if j!=i) for i in range(3))))
    bbar=[]
    for i in range(3):
        j,l=[z for z in range(3) if z!=i]
        bbar.append(s.expand(r[i]*(detK*s.prod(ds[z] for z in range(3) if z!=i)
                     +P*ds[i]*(k[i]*(r[j]+r[l]+r[j]*r[l])
                        +r[l]*K[i,j]*k[j]+r[j]*K[i,l]*k[l]))))
    data=dict(gens=gens,r=r,K=K,k=k,ds=ds,detK=detK,Z=Z,N=N,P=P,Zbar=Zbar,bbar=bbar,
              D1=s.expand((s.eye(3)+s.diag(*r)*(s.eye(3)-K)).det()))
    gens=data['gens']; r=data['r']; K=data['K']; ds=data['ds']
    poly=lambda q:s.Poly(q,*gens,domain=s.ZZ)
    N=data['N']; adj=N.adjugate()
    detN=poly(N.det())
    assert poly(data['P'])==poly(s.prod(ds))
    # Independently reconstruct the cleared renewal system from its formulas.
    for i in range(3):
        j,k=[j for j in range(3) if j!=i]
        assert ds[i]==1-K[j,k]*K[k,j]
        assert poly(N[i,i])==poly(1+ds[i]*(r[j]+r[k]+r[j]*r[k]))
        for j in range(3):
            if i==j: continue
            k=3-i-j
            assert poly(N[i,j])==poly(-r[i]*(K[i,j]*(1+ds[i]*r[k])+K[i,k]*K[k,j]))
    source=[poly(data['bbar'][i]-(data['Zbar'] if i==1 else 0)) for i in range(3)]
    adjp=[[poly(adj[i,j]) for j in range(3)] for i in range(3)]
    currents=[sum((adjp[i][j]*source[j] for j in range(3)),poly(0)) for i in range(3)]
    # Independent column-replacement determinant checks the adjugate orientation.
    from itertools import permutations
    def determinant(matrix):
        total=poly(0)
        for perm in permutations(range(3)):
            parity=(-1)**sum(perm[i]>perm[j] for i in range(3) for j in range(i+1,3))
            product=poly(parity)
            for i in range(3):product*=matrix[i][perm[i]]
            total+=product
        return total
    for col in range(3):
        replaced=[[source[i] if j==col else poly(N[i,j]) for j in range(3)] for i in range(3)]
        assert currents[col]==determinant(replaced)
    say('adjugate source agrees with all three column-replacement determinants')
    candidates=[('d1',poly(ds[1]))]
    results={}
    for J in [(1,),(1,2),(0,1,2)]:
        raw=poly(data['P'])*detN+sum((currents[i] for i in J),poly(0))
        H=raw; cancelled=[]
        say('family',J,'raw monomials',len(H.terms()))
        for name,divisor in candidates:
            multiplicity=0
            while not H.is_zero:
                quotient,remainder=H.div(divisor)
                if not remainder.is_zero: break
                H=quotient; multiplicity+=1
            if multiplicity: cancelled.append((name,multiplicity,divisor))
        assert raw==H*s.prod(divisor**m for _,m,divisor in cancelled)
        say('family',J,'cancelled',[(name,m) for name,m,_ in cancelled],
            'remaining terms',len(H.terms()),'sleep blocks',len({m[-3:] for m in H.monoms()}))
        results[J]={'raw':raw,'H':H,'cancelled':cancelled}
    data2={'gens':gens,'r':r,'families':results}
    return data2

def expand_rates(data):
    cache={}
    def trinomial(n):
        if n not in cache:
            cache[n]=[(i,j,factorial(n)//(factorial(i)*factorial(j)*factorial(n-i-j)))
                      for i in range(n+1) for j in range(n-i+1)]
        return cache[n]
    def transform(poly,row,degree,inverse=False):
        result=defaultdict(int)
        a,b=2*row,2*row+1
        for mon,c in poly.items():
            n=degree-mon[a]-mon[b]
            assert n>=0
            for i,j,weight in trinomial(n):
                key=list(mon); key[a]+=i; key[b]+=j
                result[tuple(key)]+=c*weight*((-1)**(i+j) if inverse else 1)
        return {m:c for m,c in result.items() if c}
    for J,entry in data['families'].items():
        blocks=defaultdict(dict)
        for mon,c in entry['H'].terms():blocks[mon[-3:]][mon[:6]]=int(c)
        certificate={}; negative=[]
        for sleep,original in sorted(blocks.items()):
            degrees=tuple(max(m[2*i]+m[2*i+1] for m in original) for i in range(3))
            result=original
            for row in range(3):result=transform(result,row,degrees[row])
            neg={m:c for m,c in result.items() if c<0}
            certificate[sleep]=(degrees,result)
            if neg: negative.append((sleep,len(neg),list(neg.items())[:2]))
            restored=result
            for row in (2,1,0):restored=transform(restored,row,degrees[row],inverse=True)
            assert restored==original,(J,sleep)
        entry['certificate']=certificate;entry['negative']=negative
        say('family',J,'blocks',len(certificate),'rate monomials',sum(len(v[1]) for v in certificate.values()),
            'negative blocks',len(negative),'negative witnesses',negative[:3])
    return data

def artifact(data):
    result={
        'vertices':[0,1,2], 'source_vertex':1,
        'jump_variables':['K01','K02','K10','K12','K20','K21'],
        'sleep_odds_variables':['r0','r1','r2'],
        'cancelled_factor':'d1=1-K02*K20',
        'identity':'pi(V)+sum_{i in J}(beta_i-gamma_i)=d1*H_J/(Zbar*det(N))',
        'coefficient_interpretation':'For each block, multiply jump monomial m by product_i k_i^(row_degrees[i]-m[2*i]-m[2*i+1]).',
        'families':[]}
    for J,entry in data['families'].items():
        assert not entry['negative']
        assert [(name,m) for name,m,_ in entry['cancelled']]==[('d1',1)]
        blocks=[]
        for sleep,(degrees,terms) in sorted(entry['certificate'].items()):
            assert all(isinstance(c,int) and c>0 for c in terms.values())
            assert all(all(m[2*i]+m[2*i+1]<=degrees[i] for i in range(3)) for m in terms)
            blocks.append({'sleep_exponents':list(sleep),'row_degrees':list(degrees),
                           'terms':[[list(m),c] for m,c in sorted(terms.items())]})
        result['families'].append({'J':list(J),'raw_monomials':len(entry['raw'].terms()),
            'reduced_monomials':len(entry['H'].terms()),'sleep_blocks':len(blocks),
            'rate_monomials':sum(len(block['terms']) for block in blocks),'blocks':blocks})
    path=Path(__file__).with_name("three-site-pair-deficit-certificate.json.gz")
    encoded=gzip.compress(json.dumps(result,sort_keys=True,separators=(',',':')).encode(),mtime=0)
    assert path.exists(),"The saved coefficient certificate is required."
    assert json.loads(gzip.decompress(path.read_bytes()))==result
    say('saved coefficient artifact verified',str(path),'SHA256',hashlib.sha256(path.read_bytes()).hexdigest())
    return result

def structural_audit(data):
    """Small Schur, boundary, and leading-coefficient identities from RN25."""
    a,b,c,d,e,f=data['gens'][:6]
    r0,r1,r2=r=data['r']; edges=(a,b,c,d,e,f)
    K=s.Matrix([[0,a,b],[c,0,d],[e,f,0]])
    k=[1-a-b,1-c-d,1-e-f]
    ds=[1-d*f,1-b*e,1-a*c]
    delta=s.expand((s.eye(3)-K).det())
    Z=[1+ds[i]*(r[j]+r[l]+r[j]*r[l])
       for i in range(3) for j,l in [[v for v in range(3) if v!=i]]]
    blocks={}
    for J,entry in data['families'].items():
        coeff=defaultdict(dict)
        for mon,w in entry['H'].terms():coeff[mon[-3:]][mon[:6]]=w
        blocks[J]={alpha:s.Poly.from_dict(poly,edges) for alpha,poly in coeff.items()}
    g01=a+b*f;g10=c+d*e;g21=f+e*a;g12=d+c*b
    L0=k[0]+b*k[2];L2=k[2]+e*k[0]
    Q0=a*(1-c)+b*f+c*a*b*e
    Q2=f*(1-d)+e*a+d*f*e*b
    expected={
        ((1,),(1,0,0)):ds[2]*g01*g10,
        ((1,),(0,0,1)):ds[0]*g21*g12,
        ((1,2),(1,0,0)):ds[2]*g01*g10,
        ((1,2),(0,0,1)):ds[0]*ds[2]*L2,
        ((0,1,2),(1,0,0)):ds[0]*ds[2]*L0,
        ((0,1,2),(0,0,1)):ds[0]*ds[2]*L2,
        ((1,),(2,2,2)):ds[0]*ds[2]*delta*(ds[0]*c*Q0+ds[2]*d*Q2),
        ((1,2),(2,2,3)):ds[0]**2*ds[2]**2*delta*L2,
        ((0,1,2),(3,2,2)):ds[0]**2*ds[2]**2*delta*L0,
        ((0,1,2),(2,2,3)):ds[0]**2*ds[2]**2*delta*L2}
    for (J,alpha),value in expected.items():
        assert blocks[J][alpha]==s.Poly(value,*edges),(J,alpha)
    assert max(map(sum,blocks[(1,)]))==6
    assert max(map(sum,blocks[(1,2)]))==max(map(sum,blocks[(0,1,2)]))==7
    def boundary(J,index):
        return s.Poly(sum((v.as_expr()*s.prod(r[i]**alpha[i] for i in range(3))
                           for alpha,v in blocks[J].items() if alpha[index]==0),s.Integer(0)),
                      *(edges+r))
    assert boundary((1,),2)==s.Poly(r0*ds[2]*g01*g10*(1+r1*(k[1]+d*k[2]))*Z[2],*(edges+r))
    assert boundary((1,),0)==s.Poly(r2*ds[0]*g21*g12*(1+r1*(k[1]+c*k[0]))*Z[0],*(edges+r))
    assert boundary((1,2),2)==boundary((1,),2)
    assert boundary((0,1,2),0)==boundary((1,2),0)
    qstar=k[1]+c*L0/(ds[1]*ds[2])+d*L2/(ds[1]*ds[0])
    assert s.cancel(1-qstar-c*Q0/(ds[1]*ds[2])-d*Q2/(ds[1]*ds[0]))==0

    # Generic M-matrix Schur cut formula, independent of the ARW entries.
    z0,z1,z2=s.symbols('Z0 Z1 Z2')
    A01,A02,A10,A12,A20,A21=s.symbols('A01 A02 A10 A12 A20 A21')
    b0,b1,b2,zv=s.symbols('b0 b1 b2 ZV')
    N=s.Matrix([[z0,-A01,-A02],[-A10,z1,-A12],[-A20,-A21,z2]])
    D=z0*z2-A02*A20
    V0=z2*A01+A02*A21;V2=A20*A01+z0*A21
    W0=z2*b0+A02*b2;W2=A20*b0+z0*b2
    U=z1*D-A10*V0-A12*V2
    L=(zv-b1)*D-A10*W0-A12*W2
    assert s.expand(U-N.det())==0
    adj_source=N.adjugate()*s.Matrix([-b0,zv-b1,-b2])
    assert s.expand(adj_source[1]-L)==0
    for B in [(),(0,),(2,),(0,2)]:
        Vs=sum({0:V0,2:V2}[i] for i in B)
        Ws=sum({0:W0,2:W2}[i] for i in B)
        direct=U-sum(adj_source[i] for i in (1,)+B)
        assert s.expand(D*direct-(U*(D+Ws)-L*(D+Vs)))==0
    R=s.Matrix([[0,s.Rational(4,5),0],[0,0,0],[0,s.Rational(4,5),0]])
    force=s.Matrix([s.Rational(1,10),s.Rational(1,10),s.Rational(11,10)])
    excess=(s.eye(3)-R).inv()*(s.Rational(17,20)*s.eye(3)[:,1]-force)
    assert excess==s.Matrix([s.Rational(1,2),s.Rational(3,4),-s.Rational(1,2)])
    assert excess[1]<1 and sum(excess)<1 and excess[0]+excess[1]>1
    say('structural audit: 10 coefficient identities, boundary factors, escape limit, all four Schur cuts')

if __name__=='__main__':
    data=build()
    expand_rates(data)
    artifact(data)
    structural_audit(data)
