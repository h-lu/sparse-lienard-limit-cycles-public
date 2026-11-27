"""Rational ascent signs and wall-error bounds for the eight-term support in the supplement."""
import argparse
from fractions import Fraction as F


C = {0: F(1), 1: -F(19, 10), 2: F(1)}
G = {0: F(1), 3: -F(1)}


def value(poly, x):
    return sum(c * x**i for i, c in poly.items())


def wronskian(even, odd, gap, stride):
    result = {}
    for i, a in even.items():
        for j, b in odd.items():
            n = i + j
            result[n] = result.get(n, 0) + (gap + stride*(j-i)) * a * b
    return result


def integrated(even, odd, gap, stride, exponent_sum):
    h = wronskian(even, odd, gap, stride)
    return {i: c / (exponent_sum + stride*i) for i, c in h.items()}


def fixed_eight():
    m = 2**38
    scale = 8
    odd = {**G, 16: F(511, scale**16)}
    integral = integrated(C, odd, -1, 2, 11)
    slope = {6+2*i: 8*c*scale**i for i, c in C.items()}
    even_part = {i+2: c/(i+2) for i, c in slope.items()}
    odd_part = {5+2*i: c*scale**i for i, c in odd.items()}
    assert sorted([3, m, *even_part, *odd_part]) == [3, 5, 8, 10, 11, 12, 37, 274877906944]
    assert sum(slope.values()) == F(1992, 5) and sum(even_part.values()) == F(2363, 75)
    assert sum(odd_part.values()) == 0 and sum(i*c for i, c in odd_part.items()) == 13280
    nodes = tuple(map(F, ("1/5", "1/2", "1", "2", "11/2", "7", "15/2")))
    for x, direction in zip(nodes, (-1, 1, -1, 1, -1, 1, 1)):
        rho = x/scale
        normalized = value(integral, x) / value(C, x)
        assert direction*normalized > 0 and rho**5 * normalized**2 > F(1, 400000)**2
        assert rho <= F(15, 16) and 0 < sum(c*rho**(i//2) for i, c in even_part.items()) < 21
    assert 6 + 4/(1-F(19, 20)) == 86 < m-2
    assert sum(abs(c) for c in odd_part.values()) == 1024
    assert F(400*1024, 1)/F(2, 5) < 2**20
    assert F(15, 16)**16 < F(1, 2) and (m-8)//2 > 16*60 and 2**40 > 10**9
    assert F(1, 400000)-F(1, 10**9) > F(1, 410000)
    assert sum(even_part.values())-sum(slope.values())/m > 31
    assert F(1, 1)-F(37*40, m) > F(3, 4)
    assert 2*sum(map(abs, even_part.values())) < 114 and 2*sum(map(abs, slope.values())) < 1300
    assert 400*3**16 < m and 40-F(800, m) > 39 and 2*sum(slope.values()) > 114+1
    assert 4*1300 < 10*m
    curvature_bound = sum(i*(i-1)*abs(c) for i, c in odd_part.items())
    assert curvature_bound == 736992 and F(2*curvature_bound*40, m) < 1
    assert 4*14000 < 16*13000 and (40+4)*14000*440000 < m
    print("Fourteen-cycle support verified: seven ascent tests and the wall bounds for M=2^38 pass.")


if __name__ == "__main__":
    if not __debug__:
        raise SystemExit("Enable assertions to check the enclosures.")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--terms", type=int, choices=(8,), default=8)
    parser.parse_args()
    fixed_eight()
