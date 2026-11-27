"""Exact checks for the printed degree-ten seed intervals in Article, Lemma 3.1."""
import json
from fractions import Fraction as F
from pathlib import Path

from sympy import expand, symbols

from baseline_checks import RationalInterval, pt


def cell(n, scale):
    return RationalInterval(F(n, scale), F(n + 1, scale))


def strictly_inside(interval, lower, upper):
    assert lower < interval.lo <= interval.hi < upper


def certify():
    data = json.loads((Path(__file__).resolve().parents[1] / "inputs/degree-ten-seed.json").read_text())
    x, b, a, sign = symbols("x b a sign")
    poly = x**14 - 9*x**10 - 8*x**9 - 8*x**5 - 9*x**4 + 1
    factor = (x+1)**4 * (x**4-x**3+x**2-x+1) * ((x**2-x+1)**3-3*x**3)
    assert expand(poly-factor) == 0
    v = cell(data["radius"]["left_integer"], 10**data["radius"]["digits"])
    assert poly.subs(x, v.lo) < 0 < poly.subs(x, v.hi)
    b0 = pt(F(10, 9)) * (v.power(5)+pt(1)) * (v.power(4)-pt(1)).reciprocal()
    a0 = (pt(9)*b0+pt(10))*pt(F(1, 5))

    def profile(r, s):
        return r.power(10) + pt(s)*(b0*r.power(9)-a0*r.power(5))

    def derivative_factor(r):
        return pt(10)*r.power(5)-pt(9)*b0*r.power(4)+pt(5)*a0

    # The same six terms give the symbolic primitive and its interval extension.
    terms = [(F(100,18),18,0,0,0), (F(180,17),17,1,0,1),
             (F(81,16),16,2,0,0), (-F(100,13),13,0,1,1),
             (-F(90,12),12,1,1,0), (F(25,8),8,0,2,0)]
    primitive = sum(c*x**i*b**j*a**k*sign**m for c,i,j,k,m in terms)
    for s in (1, -1):
        q = x**10+s*(b*x**9-a*x**5)
        assert expand(x*primitive.subs(sign,s).diff(x)-q.diff(x)**2) == 0

    def action(r, s):
        return sum((pt(c*s**m)*r.power(i)*b0.power(j)*a0.power(k)
                    for c,i,j,k,m in terms), pt(0))

    peak = cell(data["peak"]["left_integer"], 10**data["peak"]["digits"])
    assert derivative_factor(pt(peak.lo)).lo > 0 > derivative_factor(pt(peak.hi)).hi
    assert peak.hi < F(18,25)*b0.lo < F(18,25)*b0.hi < v.lo
    strictly_inside(profile(peak, -1), 0, 14)
    strictly_inside(pt(50)+pt(36)*b0, 130, 145)
    strictly_inside(v.power(7)*(pt(50)*v-pt(36)*b0), 850, 1100)
    assert 145*v.hi**2 < 850
    print("PASS degree-ten factorization, isolating cell, peak, curvature and primitive identities.")
    for h, plus, minus, lower, upper in data["tests"]:
        tp, tm = cell(plus, 10**data["turn_digits"]), cell(minus, 10**data["turn_digits"])
        for t, s in [(tp,1),(tm,-1)]:
            assert profile(pt(t.lo),s).hi < h < profile(pt(t.hi),s).lo
        gap = action(tp,1)-action(pt(1),1)-action(tm,-1)+action(v,-1)
        strictly_inside(gap, lower, upper)
        print(f"PASS h={h}: outer radius cells {plus},{minus} at 1e-9; {lower} < action gap < {upper}.")


if __name__ == "__main__":
    certify()
