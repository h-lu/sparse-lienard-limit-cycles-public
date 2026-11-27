"""Exact endpoint and action inequalities for the two constructions."""

from dataclasses import dataclass
from fractions import Fraction as Q


@dataclass(frozen=True)
class Interval:
    lo: Q
    hi: Q

    @classmethod
    def point(cls, value):
        return cls(Q(value), Q(value))

    def __add__(self, other):
        return Interval(self.lo + other.lo, self.hi + other.hi)

    def __neg__(self):
        return Interval(-self.hi, -self.lo)

    def __sub__(self, other):
        return self + (-other)

    def __mul__(self, other):
        products = [a * b for a in (self.lo, self.hi) for b in (other.lo, other.hi)]
        return Interval(min(products), max(products))

    def reciprocal(self):
        assert self.lo > 0
        return Interval(1 / self.hi, 1 / self.lo)

    def power(self, exponent):
        assert self.lo >= 0 and exponent >= 0
        return Interval(self.lo**exponent, self.hi**exponent)


def value(poly, x):
    return sum((c * x.power(n) for n, c in poly.items()), Interval.point(0))


def root_bracket(poly, height, lower):
    """The caller establishes increase on the outer branch analytically."""
    lower = Q(lower)
    upper = lower + Q(3, 10**14)
    assert value(poly, Interval.point(lower)).hi < height
    assert value(poly, Interval.point(upper)).lo > height
    return Interval(lower, upper)


def action(poly, lower, upper):
    derivative = [(n, Interval.point(n) * c) for n, c in poly.items()]
    total = Interval.point(0)
    for i, a in derivative:
        for j, b in derivative:
            n = i + j - 2
            assert n > 0
            total += a * b * (upper.power(n) - lower.power(n)) * Interval.point(Q(1, n))
    return total


def verify_table(name, profiles, lower_limits, rows):
    for height, plus_root, minus_root, limits in rows:
        height = Q(height)
        roots = [root_bracket(poly, height, left)
                 for poly, left in zip(profiles, (plus_root, minus_root))]
        assert all(lower.hi < upper.lo for lower, upper in zip(lower_limits, roots))
        actions = [action(poly, lower, upper)
                   for poly, lower, upper in zip(profiles, lower_limits, roots)]
        diff = actions[0] - actions[1]
        low, high = map(Q, limits)
        assert low < diff.lo <= diff.hi < high
        print(f"{name}   h = {height}   {low} < S_plus - S_minus < {high}")


def outer_minimum():
    one = Interval.point(1)
    profiles = [{10: one, 11: Interval.point(sign * 2), 12: one} for sign in (-1, 1)]
    assert Q(5, 6)**10 * Q(1, 6)**2 < Q(1, 100)
    verify_table("outer minimum", profiles, (one, Interval.point(0)), [
        ("1/100", "1.07097474629757", ".57608191657593", ("3/10000", "4/10000")),
        ("1", "1.28519903324534", ".88127146163356", ("-3/2", "-1")),
        ("2", "1.33433799530645", ".93884595326092", ("-5", "-4")),
    ])


def degree_ten():
    one = Interval.point(1)
    r = Interval(Q("1.92194269368393019"), Q("1.92194269368393020"))

    def well_equation(x):
        return 9 * (x**10 - 1) * (x**4 - 1) - 8 * (x**5 + 1) * (1 + x**9)

    assert well_equation(r.lo) < 0 < well_equation(r.hi)
    b = Interval.point(Q(10, 9)) * (r.power(5) + one) * (r.power(4) - one).reciprocal()
    a = (Interval.point(9) * b + Interval.point(10)) * Interval.point(Q(1, 5))
    profiles = [{5: -a, 9: b, 10: one}, {5: a, 9: -b, 10: one}]
    assert value(profiles[0], one).hi < -2

    # The derivative of V divided by x^4 has its sole minimum at 18 B / 25.
    inner = Interval(Q("1.46100899849472"), Q("1.46100899849473"))
    critical = Interval.point(Q(18, 25)) * b
    assert inner.hi < critical.lo and critical.hi < r.lo
    reduced = {0: Interval.point(5) * a, 4: Interval.point(-9) * b,
               5: Interval.point(10)}
    assert value(reduced, Interval.point(inner.lo)).lo > 0
    assert value(reduced, Interval.point(inner.hi)).hi < 0
    assert value(profiles[1], inner).hi < 14

    verify_table("degree ten", profiles, (one, r), [
        ("20", "1.32047248100449", "2.08589292441446", ("-89", "-88")),
        ("100", "1.49773775614208", "2.20694152078300", ("1293", "1295")),
        ("200", "1.59405002861009", "2.27907463266905", ("6884", "6886")),
    ])


if __name__ == "__main__":
    if not __debug__:
        raise SystemExit("Enable assertions for the rational enclosures.")
    outer_minimum()
    degree_ten()
