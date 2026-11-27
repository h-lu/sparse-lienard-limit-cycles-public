"""Exact algebra for the moving endpoints, the binomial example and the outer oval."""
import sympy as sp


def elimination_identity():
    """Denominator-cleared differentiation in Article, Proposition 3.2."""
    x, mu = sp.symbols("x mu")
    g, h = (sp.Function(name)(x) for name in ("g", "h"))
    bracket = lambda a, b: a * sp.diff(b, x) - sp.diff(a, x) * b
    g1, h1 = sp.diff(g, x) - mu*g, sp.diff(h, x) - mu*h
    wr = bracket(g, h)
    left = g1 * (sp.diff(wr, x) - mu*wr) - sp.diff(g1, x)*wr
    assert sp.expand(left - g * bracket(g1, h1)) == 0
    print("Symbolic Euler elimination identity verified.")


def crest_variations():
    r, t, delta, speed = sp.symbols("r t delta speed", real=True)
    even = r**2 / 2 - r**4 / 4
    odd = r**3 - r**5
    profile = even + delta * odd
    slope = sp.cancel(even.diff(r) / r)
    rise, fall = sp.Rational(1, 5), sp.Rational(7, 5)
    assert even.subs(r, rise) == even.subs(r, fall)
    assert slope.subs(r, rise) > 0 > slope.subs(r, fall)
    assert even.diff(r).subs(r, 1) == 0 == odd.subs(r, 1)

    def primitive(integrand, start):
        return sp.integrate(sp.expand(integrand), (r, start, t))

    def action(start):
        return primitive(sp.cancel(profile.diff(r)**2 / r), start)

    def moving_cost(cost, endpoint, height_speed):
        # E(R(delta)) + delta Q(R(delta)) = height + delta * height_speed.
        turn_speed = (height_speed - odd.subs(r, endpoint)) / even.diff(r).subs(r, endpoint)
        change = cost.diff(delta) + turn_speed * cost.diff(t)
        return change.subs({delta: 0, t: endpoint})

    variations = {}
    for start, endpoint in [(0, rise), (1, fall)]:
        cost = action(start)
        expected = primitive(slope * odd.diff(r) - slope.diff(r) * odd, start)
        variations[endpoint] = expected.subs(t, endpoint)
        assert moving_cost(cost, endpoint, 0) == variations[endpoint]
        height_derivative = cost.diff(t).subs(delta, 0) / even.diff(r).subs(r, t)
        assert sp.simplify(height_derivative - slope.subs(r, t)) == 0

    # Equal heights at a fixed plus turn move the minus turn to first order.
    minus_profile = even - delta * odd
    equality = minus_profile.subs(r, rise + delta * speed) - profile.subs(r, rise)
    turn_velocity, = sp.solve(equality.diff(delta).subs(delta, 0), speed)
    assert turn_velocity == 2 * odd.subs(r, rise) / even.diff(r).subs(r, rise)

    # Equal outgoing costs on two perturbed descending branches, by reflection.
    outgoing = action(1)
    level_velocity, = sp.solve(moving_cost(outgoing, fall, speed), speed)
    assert level_velocity == -variations[fall] / slope.subs(r, fall)
    incoming_change = moving_cost(action(0), rise, level_velocity)
    reflected = variations[rise] - slope.subs(r, rise) * variations[fall] / slope.subs(r, fall)
    assert incoming_change == reflected and reflected != 0
    print("Exact fixed-level and equal-cost derivatives agree on ascending and descending inverse branches.")
    print(f"At equal-height radii {rise} and {fall}, f(descending)={slope.subs(r, fall)}, reflected variation={reflected}.")


def binomial_and_lift():
    radius, angle = sp.symbols("R t", real=True)
    x = radius * sp.cos(angle)
    drift = sp.integrate(x**2 - x**4, (angle, 0, 2 * sp.pi))
    root = 2 / sp.sqrt(3)
    assert sp.expand(drift - sp.pi * radius**2 * (1 - 3 * radius**2 / 4)) == 0
    assert drift.subs(radius, root) == 0 and drift.diff(radius).subs(radius, root) < 0

    x, y, primitive, v = sp.symbols("x y P v", real=True)
    vprime = 1 / (1 + primitive**2)
    along_flow = (x - y * vprime) * (y - primitive) - (y - v) * x
    completed = -x * (primitive - v) - vprime * (y - primitive / 2)**2 + vprime * primitive**2 / 4
    assert sp.simplify(along_flow - completed) == 0
    print("The averaged binomial energy drift has a simple zero at 2/sqrt(3); the outer-oval derivative identity holds symbolically.")


if __name__ == "__main__":
    if not __debug__:
        raise SystemExit("Enable assertions for the symbolic checks.")
    elimination_identity()
    crest_variations()
    binomial_and_lift()
