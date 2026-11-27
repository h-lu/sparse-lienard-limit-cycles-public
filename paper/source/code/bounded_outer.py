"""Exact finite certificates for the local bounded- and exterior-return construction."""

from fractions import Fraction as F


def cubic_start():
    separation = F(15, 8)
    lift = F(1, 2**30)
    even = {0: 1, 2: -3 + lift, 3: 2}
    odd = {0: 1, 1: -F(3, 2) / separation + lift,
           3: 1 / (2 * separation**3)}
    coefficients = [F(0)] * 7
    for i, c in even.items():
        for j, g in odd.items():
            exponent = i + j
            factor = F(1 + 1024 * (j - i), 17 + 1024 * exponent)
            coefficients[exponent] += c * g * factor
    nodes = [F(1, 16), F(1, 6), 1, separation, 4, 32, 4096]
    scaled_floors = []
    for index, z in enumerate(nodes):
        value = sum(coef * z**i for i, coef in enumerate(coefficients))
        assert (-1)**index * value > F(1, 1000)
        scaled_floors.append((10**6 * value) // (1 + z**6))
    print("At s=17 the seven signed cubic tests exceed 1/1000 in magnitude.")
    print("Floors of 10^6 J_17(z)/(1+z^6):", *scaled_floors)


def evaluate(poly, radius):
    return sum(coef * radius**degree for degree, coef in poly.items())


def descending(level, even):
    lo, hi = F(1), F(3, 2)
    assert evaluate(even, lo) > level > evaluate(even, hi)
    for _ in range(86):
        mid = (lo + hi) / 2
        at_mid = evaluate(even, mid)
        if at_mid == level:
            return mid, mid
        if at_mid > level:
            lo = mid
        else:
            hi = mid
    return lo, hi


def on_interval(poly, interval):
    """Termwise rational bounds for a polynomial on a positive interval."""
    lo, hi = interval
    low = high = F(0)
    for power, coef in poly.items():
        ends = coef * lo**power, coef * hi**power
        low += min(ends)
        high += max(ends)
    return low, high


def divide(num, denom):
    assert denom[1] < 0 or denom[0] > 0
    quotients = [x / y for x in num for y in denom]
    return min(quotients), max(quotients)


def enclosed(interval, lower, upper):
    assert lower < interval[0] <= interval[1] < upper


def actions(slope, rho):
    odd = {9: F(1), 11: -1 - rho, 13: rho}
    action = {}
    for e, a in slope.items():
        for p, b in odd.items():
            power = e + p
            action[power] = action.get(power, 0) + a * b * F(p - e, p + e)
    tail = dict(action)
    tail[0] = -evaluate(action, F(1))
    return action, tail


def small_support():
    rho = F(6901267, 10**6)
    even = {10: F(1), 18: -F(5, 9)}
    slope = {power - 2: power * coef for power, coef in even.items()}
    action, tail = actions(slope, rho)
    rows = [(F(1, 5), 10**9, (5, 6), (4, 5)),
            (F(1, 2), 10**6, (-296, -295), (-306, -305)),
            (F(3, 4), 10**4, (-130, -129), (-135, -134))]
    for a, unit, left_bounds, right_bounds in rows:
        z = descending(evaluate(even, a), even)
        left = evaluate(action, a) / evaluate(slope, a)
        corr_lo, corr_hi = divide(on_interval(tail, z), on_interval(slope, z))
        enclosed((left * unit, left * unit), *left_bounds)
        enclosed(((left - corr_hi) * unit, (left - corr_lo) * unit), *right_bounds)
    levels = [(F(0), 10**9, (-10, -9)),
              (F(1, 10**8), 10**9, (-9, -8)),
              (F(1, 2), 100, (18, 19)),
              (F(1), 100, (61, 62))]
    for depth, unit, bounds in levels:
        interval = on_interval(tail, descending(-depth, even))
        enclosed(tuple(v * unit for v in interval), *bounds)
    _, other_tail = actions(slope, rho + F(1, 10**6))
    other_end = on_interval(other_tail, descending(0, even))
    enclosed(tuple(v * 10**9 for v in other_end), 15, 16)
    print("Degree 18: both actions at three rising radii and the exterior action "
          "at zero and three positive depths pass exact rational enclosures.")


if __name__ == "__main__":
    if not __debug__:
        raise SystemExit("Run without -O to check all inequalities.")
    cubic_start()
    small_support()
