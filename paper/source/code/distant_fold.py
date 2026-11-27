"""Exact rational bounds for the degree-ten seed curvature and degree-eighteen fold cancellation."""

from fractions import Fraction as F


def evaluate(poly, x):
    return sum(coef * x**power for power, coef in poly.items())


def polynomial_bounds(poly, lo, hi):
    """Termwise enclosure for nonnegative arguments."""
    assert 0 <= lo <= hi
    lower = upper = F(0)
    for power, coef in poly.items():
        ends = (coef * lo**power, coef * hi**power)
        lower += min(ends)
        upper += max(ends)
    return lower, upper


def degree_ten_minima():
    polynomial = {14: 1, 10: -9, 9: -8, 5: -8, 4: -9, 0: 1}
    lo, hi = F(19219, 10000), F(19220, 10000)
    assert evaluate(polynomial, lo) < 0 < evaluate(polynomial, hi)
    b_lo = 10 * (lo**5 + 1) / (9 * (hi**4 - 1))
    b_hi = 10 * (hi**5 + 1) / (9 * (lo**4 - 1))
    # At the two critical radii, use q''_+(1)=50+36b and
    # q''_-(v)=v**7*(50v-36b).
    assert 130 < 50 + 36 * b_lo <= 50 + 36 * b_hi < 145
    assert 850 < lo**7 * (50 * lo - 36 * b_hi)
    assert hi**7 * (50 * hi - 36 * b_lo) < 1100
    assert hi**2 / 850 < F(1, 145)
    print("The incoming integral has the required sign near the common seed minimum.")


def terminal_cancellation():
    # (9/5)**(1/8) is the descending zero in degree eighteen.
    lo, hi = F(1), F(3, 2)
    for _ in range(90):
        mid = (lo + hi) / 2
        if mid**8 < F(9, 5):
            lo = mid
        else:
            hi = mid
    for rho, sign in [(F(6901267, 10**6), -1), (F(6901268, 10**6), 1)]:
        action = {}
        for e, f in {8: 10, 16: -10}.items():
            for p, q in {9: 1, 11: -1 - rho, 13: rho}.items():
                coefficient = f * q * F(p - e, p + e)
                action[p + e] = action.get(p + e, 0) + coefficient
        at_crest = evaluate(action, F(1))
        assert at_crest < 0
        interval = polynomial_bounds(action, lo, hi)
        assert min(sign * (value - at_crest) for value in interval) > 0
    print("The degree-eighteen terminal cancellation lies in (6.901267, 6.901268).")


if __name__ == "__main__":
    if not __debug__:
        raise SystemExit("Run without -O to check all inequalities.")
    degree_ten_minima()
    terminal_cancellation()
