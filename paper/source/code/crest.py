"""Algebraic checks for the origin and single-crest construction.

Run with Python 3 and mpmath. The coefficient calculation is rational.
All finite-wall signs use directed interval arithmetic and certified brackets.
"""
import argparse
from fractions import Fraction

from mpmath import iv


def evaluate(coefficients, z):
    value = 0
    for coefficient in reversed(coefficients):
        value = coefficient + z * value
    return value


def rational_block(m):
    e, q = 8, 9
    roots = [Fraction(i, 64 * m) for i in range(1, m)] + [Fraction(1)]
    action = [Fraction(q - e, q + e)]
    for root in roots:
        new = action + [Fraction(0)]
        for j, coefficient in enumerate(action):
            new[j + 1] -= coefficient / root
        action = new
    odd = [coefficient * Fraction(q + 2*j + e, q + 2*j - e)
           for j, coefficient in enumerate(action)]
    assert odd[0] == 1 and all(odd)
    assert all(evaluate(action, root) == 0 for root in roots)
    before_last = (-1)**(m - 1)
    assert before_last * evaluate(odd, 1) > 0
    assert before_last * odd[-1] < 0
    right = 2
    while before_last * evaluate(odd, right) >= 0:
        right *= 2
    print(f"Rational block with {m} added odd powers: {m} simple action roots;"
          f" final primitive root has 1 < r^2 < {right}.")


def number(x):
    x = Fraction(x)
    return iv.mpf(x.numerator) / x.denominator


def sign(x):
    if x > 0:
        return 1
    if x < 0:
        return -1
    raise ArithmeticError(f"interval straddles zero: {x}")


def check_enclosure(value, bounds, denominator):
    lo, hi = (Fraction(numerator, denominator) for numerator in bounds)
    assert number(lo) < value < number(hi)
    return f"({lo}, {hi})"


def check_wall(odd_power, wall, tests, depth_bounds):
    iv.dps = 95
    d = 10
    print(f"Five-term support 5, 9, 10, {odd_power}, {wall}.")

    def profile(r):
        return r**d - number(Fraction(d, wall))*r**wall

    def slope(r):
        return d*(r**(d - 2) - r**(wall - 2))

    def action(r, start):
        total = iv.mpf(0)
        for e, weight in ((d - 2, d), (wall - 2, -d)):
            for p, odd in ((9, 1), (odd_power, -1)):
                coefficient = number(Fraction(weight * odd * (p - e), p + e))
                total += coefficient * (r**(p + e) - start**(p + e))
        return total

    def descending(level):
        # Use the regular coordinate u = M(r - 1), with exact dyadic midpoints.
        lo, hi = Fraction(0), Fraction(32)

        def residual(u):
            return profile(1 + number(u) / wall) - level

        assert sign(residual(lo)) == 1 and sign(residual(hi)) == -1
        for _ in range(125):
            mid = (lo + hi) / 2
            if sign(residual(mid)) == 1:
                lo = mid
            else:
                hi = mid
        return 1 + iv.mpf([number(lo).a, number(hi).b]) / wall

    for a, left_bounds, right_bounds, denominator in tests:
        radius = number(a)
        z = descending(profile(radius))
        left = action(radius, 0) / slope(radius)
        right = left - action(z, 1) / slope(z)
        assert sign(slope(radius)) == 1 and sign(slope(z)) == -1
        left_interval = check_enclosure(left, left_bounds, denominator)
        right_interval = check_enclosure(right, right_bounds, denominator)
        print(f"a={a}: U/f in {left_interval}, B_R/f in {right_interval}")
    for s, bounds in zip((0, Fraction(1, 2), 1), depth_bounds):
        terminal = action(descending(-number(s)), 1)
        result = check_enclosure(terminal, bounds, 1)
        print(f"depth {s}: descending action in {result}")


if __name__ == "__main__":
    if not __debug__:
        raise SystemExit("Run without -O to enable the certificate assertions.")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--new-odd", type=int, choices=(1, 3), default=1,
                        help="number of added odd powers in the rational block")
    args = parser.parse_args()
    rational_block(args.new_odd)
    if args.new_odd == 1:
        check_wall(11, 20_000_000, [
            (Fraction(1, 2), (37, 39), (36, 37), 1_000_000),
            (Fraction(3, 4), (-23, -22), (-23, -22), 10_000),
            (Fraction(4, 5), (-57, -56), (-57, -56), 10_000),
        ], [(-26, -25), (-39, -38), (-53, -52)])
        check_wall(29, 16_384, [
            (Fraction(5, 6), (8, 9), (2, 3), 1_000),
            (Fraction(19, 20), (-92, -91), (-97, -96), 1_000),
            (Fraction(39, 40), (-226, -225), (-231, -230), 1_000),
        ], [(-110, -109), (-176, -175), (-247, -246)])
