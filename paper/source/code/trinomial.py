#!/usr/bin/env python3
"""Exact rational certificates for the finite-degree sign table (Python 3)."""
if not __debug__:
    raise SystemExit("Enable assertions for the rational enclosures.")

from fractions import Fraction as F
from dataclasses import dataclass

@dataclass(frozen=True)
class Box:
    lo: F
    hi: F
    def __add__(s,t):
        return Box(s.lo+t.lo,s.hi+t.hi)
    def __neg__(s):
        return Box(-s.hi,-s.lo)
    def __sub__(s,t):
        return s+-t
    def __mul__(s,t):
        p=[x*y for x in (s.lo,s.hi) for y in (t.lo,t.hi)]
        return Box(min(p),max(p))
    def power(s,n):
        assert s.lo>=0 and n>=0
        return Box(s.lo**n,s.hi**n)
    def reciprocal(s):
        assert s.lo>0
        return Box(1/s.hi,1/s.lo)

def point(x):
    return Box(F(x),F(x))

def root(f,lo,hi):
    lo,hi=F(lo),F(hi)
    assert f(lo)<0<f(hi)
    for _ in range(64):
        m=(lo+hi)/2
        v=f(m)
        if v==0:
            return point(m)
        if v<0:
            lo=m
        else:
            hi=m
    return Box(lo,hi)

def certify(d,c,L,U):
    a=Box(F(c)-F(1,10**9),F(c)+F(1,10**9))
    tau,rho=F(d-2,d-1),F(d-1,d)
    assert 0<a.lo<a.hi<tau**(d-2)/(d-1)
    s=[]; r=[]; t=[]
    for A in (a.lo,a.hi):
        s.append(root(lambda x:x**(d-2)*(d*x+d-1)-A,0,1))
        r.append(root(lambda x:x**(d-2)*(d*x-d+1)+A,tau,rho))
        t.append(root(lambda x:x**(d-1)-A,0,1))
    def f(x,A,z):
        return x.power(d)+point(z)*(x.power(d-1)-point(A)*x)
    assert (f(s[0],a.lo,1)-f(r[0],a.lo,-1)).lo>0
    assert (f(s[1],a.hi,1)-f(r[1],a.hi,-1)).hi<0
    s=Box(s[0].lo,s[1].hi)
    r=Box(r[1].lo,r[0].hi)
    t=Box(t[0].lo,t[1].hi)
    assert s.hi<t.lo and r.hi<1
    def Q(x,z):
        return (point(F(d*d,2*d-2))*x.power(2*d-2)
              + point(F(z*2*d*(d-1),2*d-3))*x.power(2*d-3)
              + point(F((d-1)**2,2*d-4))*x.power(2*d-4)
              - point(F(z*2*d,d-1))*a*x.power(d-1)
              - point(F(2*(d-1),d-2))*a*x.power(d-2))
    h=t*r*s.reciprocal()
    log_h=Box(1-1/h.lo,h.hi-1)
    J=point(1000*d)*(Q(t,1)-Q(s,1)-Q(point(1),-1)+Q(r,-1)+a.power(2)*log_h)
    if not L<J.lo<=J.hi<U:
        raise ArithmeticError(f"certificate failed for d={d}: outside ({L}, {U})")
    print(f"d={d:2}: {L:4} < 1000*d*I_d < {U:3} (exact)")

for row in [(6,'0.062500000',130,135),(8,'0.040486856',25,27),
            (10,'0.029550858',-24,-22),(12,'0.023104700',-51,-50),
            (14,'0.018889538',-69,-67),(16,'0.015933561',-81,-79),
            (18,'0.013753521',-89,-88)]:
    certify(*row)

assert F(19,20)**19 < F(19,50)
assert 20*F(5,3)*F(2,3)**19 < F(1,60)
assert 20*F(3,4)**19 < F(1,6)
assert (F(19,50)+F(1,60))/(F(18,19)+F(2,3)) < F(1,4)
assert F(245,1296)-F(1,5) == -F(71,6480)
print("Degree >= 20: all rational constants in the analytic bound verified.")
