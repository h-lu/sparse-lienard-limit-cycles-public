"""Whole-domain interval checks for the degree-22 valley in Corollary 5.2."""
from flint import arb, ctx
from sympy import Poly, Rational, diff, integrate, symbols

ctx.prec = 512
ETA = arb(1) / 2**20
QAMP = arb(1) / 2**10


def even(r):
    return r**10 * ((1 - r**2)**2 + ETA * r**2)


def odd(r):
    return r**9 * ((1 - r**2)**2 - QAMP * r**2)


def positive_tests():
    # The zeros of E'(r)/r^9 are the two roots in z=r^2.
    middle = 12 * (2 - ETA)
    discriminant = middle**2 - 4 * 14 * 10
    crest = ((middle - discriminant.sqrt()) / 28).sqrt()
    valley = ((middle + discriminant.sqrt()) / 28).sqrt()
    low, high = arb(0), arb(1) / 2
    minimum = even(valley)
    assert 0 < minimum and high < crest < valley < 1
    assert even(low) < minimum < even(high)
    for _ in range(128):
        mid = (low + high) / 2
        if even(mid) < minimum:
            low = mid
        elif even(mid) > minimum:
            high = mid
        else:
            raise ArithmeticError("Unresolved branch comparison")
    a = low.union(high)
    outer = valley + arb(1) / 2**30
    assert 0 < even(a / 2) < minimum < even(outer) < even(crest) < even(arb(2))
    r = symbols("r")
    e = r**10 * ((1-r**2)**2 + Rational(1, 2**20)*r**2)
    q = r**9 * ((1-r**2)**2 - Rational(1, 2**10)*r**2)
    f = (diff(e, r)/r).expand()
    u = integrate(f*diff(q, r) - diff(f, r)*q, r)
    assert u.subs(r, 0) == 0

    def value(poly, at):
        result = arb(0)
        for coeff in Poly(poly, r).all_coeffs():
            num, den = coeff.as_numer_denom()
            result = result * at + arb(int(num)) / int(den)
        return result

    edge = value(u, a) + value(f, a)*odd(valley)
    last = value(u, arb(2)) + edge - value(u, valley)
    first = value(u, a / 2)
    for label, inner in (("V_o", edge), ("V_i", edge + value(u, outer) - value(u, valley))):
        assert first > 0 and inner < 0 and last > 0
        print(f"{label} at E(a/2), E(c+2^-30), E(2):", *(v.str(9) for v in (first, inner, last)))


def negative_tests():
    # Fold [1,infinity) at r=1/s. The resulting density on [0,1]
    # is s^12*B/(H*s^14+C)^2. The other density is s^10*B/(H+s^10*C)^2.
    # B=(1-s^2)^2-q*s^2 and C=(1-s^2)^2+eta*s^2.
    cells = 2**17
    step = arb(1) / cells
    depths = (ETA**2, ETA, arb(1))
    integrals = [arb(0) for _ in depths]
    for i in range(cells):
        lo, hi = i*step, (i+1)*step
        square = ((1-hi**2)**2).union((1-lo**2)**2)
        z = (lo**2).union(hi**2)
        base = square + ETA * z
        numerator = square - QAMP * z
        s10 = (lo**10).union(hi**10)
        s12 = (lo**12).union(hi**12)
        s14 = (lo**14).union(hi**14)
        for j, depth in enumerate(depths):
            direct = depth + s10 * base
            inverse = depth * s14 + base
            assert direct > 0 and inverse > 0
            integrals[j] += step * numerator * (s10 / direct**2 + s12 / inverse**2)
    print(f"Reciprocal substitution on [1,infinity); {cells} equal full cells on [0,1].")
    bounds = ((6990000000, 7020000000), (-265000, -260000), ("0.0537", "0.0538"))
    for name, value, (lower, upper) in zip(("eta^2", "eta", "1"), integrals, bounds, strict=True):
        assert arb(lower) < value < arb(upper)
        print(f"J({name}) endpoints:", value.lower().str(12), value.upper().str(12))


if __name__ == "__main__":
    if not __debug__:
        raise SystemExit("Run with assertions enabled.")
    print("Degree-22 finite profile: d=10, p=2, eta=2^-20, q=2^-10; Arb, 512 bits.")
    positive_tests()
    negative_tests()
