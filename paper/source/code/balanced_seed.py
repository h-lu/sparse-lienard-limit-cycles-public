"""Degree-ten balanced seed: exact rational action enclosures at heights 20, 100 and 200."""
import argparse
from fractions import Fraction as Q


def isolate(f, lo, hi, steps=84):
    """Enclose an isolated negative-to-positive crossing on the given interval."""
    lo, hi = Q(lo), Q(hi)
    assert f(lo) < 0 < f(hi)
    for _ in range(steps):
        mid = (lo + hi) / 2
        val = f(mid)
        if val == 0:
            return mid, mid
        if val < 0:
            lo = mid
        else:
            hi = mid
    return lo, hi


def seed_parameters():
    def well_polynomial(z):
        return z**14 - 9*z**10 - 8*z**9 - 8*z**5 - 9*z**4 + 1

    zl, zu = isolate(well_polynomial, Q(19, 10), 2)
    bl = Q(10, 9) * (zl**5 + 1) / (zu**4 - 1)
    bu = Q(10, 9) * (zu**5 + 1) / (zl**4 - 1)
    al, au = (9*bl + 10)/5, (9*bu + 10)/5
    assert 1 + bu - al < -2
    c, e = Q(1460, 1000), Q(1462, 1000)
    assert e < 18*bl/25 < 18*bu/25 < zl
    assert 10*c**5 - 9*bu*c**4 + 5*al > 0
    assert 10*e**5 - 9*bl*e**4 + 5*au < 0
    assert e**10 - bl*c**9 + au*e**5 < 20
    return (al, au), (bl, bu), (zl, zu)


def seed_turn(height, side, a, b):
    al, au = a
    bl, bu = b

    def profile(x, aa, bb):
        return x**10 + side * (bb*x**9 - aa*x**5) - height

    if side == 1:
        start, lower_coeffs, upper_coeffs = Q(5, 4), (au, bl), (al, bu)
        assert 10*start**5 + 9*bl*start**4 - 5*au > 0
    else:
        start, lower_coeffs, upper_coeffs = Q(2), (al, bu), (au, bl)
        assert start > 18*bu/25 and 10*start**5 - 9*bu*start**4 + 5*al > 0
    lo, _ = isolate(lambda t: profile(t, *upper_coeffs), start, 4)
    _, hi = isolate(lambda t: profile(t, *lower_coeffs), start, 4)
    return lo, hi


def seed_primitive(a, b, t, side):
    # Each unsigned factor a^i b^j t^k is positive and increasing.
    terms = [(Q(100,18), 0, 0, 18), (side*Q(180,17), 0, 1, 17),
             (Q(81,16), 0, 2, 16), (-side*Q(100,13), 1, 0, 13),
             (-Q(90,12), 1, 1, 12), (Q(25,8), 2, 0, 8)]
    lo, hi = Q(0), Q(0)
    for weight, i, j, k in terms:
        ends = [weight * a[m]**i * b[m]**j * t[m]**k for m in (0, 1)]
        lo += min(ends)
        hi += max(ends)
    return lo, hi


def check_seed():
    a, b, z = seed_parameters()
    for h, lower, upper in [(20, -89, -88), (100, 1293, 1295), (200, 6884, 6886)]:
        actions = []
        for side, well in [(1, (Q(1), Q(1))), (-1, z)]:
            turn = seed_turn(h, side, a, b)
            tl, tu = seed_primitive(a, b, turn, side)
            wl, wu = seed_primitive(a, b, well, side)
            actions.append((tl - wu, tu - wl))
        (pl, pu), (ml, mu) = actions
        assert lower < pl - mu <= pu - ml < upper
        print(f"seed h={h}: action difference in ({lower},{upper})")


if __name__ == "__main__":
    if not __debug__:
        raise SystemExit("Run without -O so all rational certificate checks are enabled.")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--degree", type=int, choices=(10,), default=10,
                        help="degree of the balanced seed")
    parser.parse_args()
    check_seed()
