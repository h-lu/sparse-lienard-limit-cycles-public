"""Arb enclosures for the fold and finite valley profile."""
from fractions import Fraction
from importlib.metadata import version

from flint import arb, ctx


def ball(value):
    return arb(str(value))


def sign(value):
    if value > 0:
        return "+"
    if value < 0:
        return "-"
    raise ArithmeticError(f"Undetermined sign {value}")


def action_terms(slope, odd):
    """Coefficients of the primitive of f Q' - f' Q."""
    terms = {}
    for i, a in slope.items():
        for j, b in odd.items():
            power = i + j
            terms[power] = terms.get(power, 0) + a * b * Fraction(j - i, power)
    return terms


def polynomial(terms, at):
    return sum(ball(c) * at**p for p, c in terms.items())


def fold_directions():
    print("\nThree-equation fold, finite algebra (no finite Miranda root computed)")
    for degree in (10, 14, 18, 67080):
        g = (1 + arb(8) / degree).root(8)
        f = {degree - 2: degree, degree + 6: -degree}

        def variation(odd):
            terms = action_terms(f, odd)
            return sum(ball(c) * (g**p - 1) for p, c in terms.items())

        # The first row of the two-column matrix is (1,1).
        det = variation({degree + 1: -1, degree + 3: 1})
        intercept = variation({degree - 1: 1, degree + 1: -1})
        assert det > 0 and intercept < 0
        rho = -intercept / det
        if degree == 10:
            assert rho > ball('6.901267') and rho < ball('6.901268')
        print(f"d={degree}: det={det.str(9)}, rho={rho.str(10)}")


def exact_crest_hits():
    """Piecewise affine even profile in u, phi=u then 2-u, joined at u=1."""
    print("\nActual outgoing legs on a piecewise affine crest")
    for n in (16, 32, 64, 128):
        eps = arb(1) / n

        def upper_square(turn):
            # Solve b_v=eps/b-1 to the crest. Then b_v=eps/b+1 to the axis.
            y = -(-(-1 - (turn - 1) / eps).exp()).lambertw()
            crest_gap = eps * (1 - y)
            x1 = 1 + crest_gap / eps
            x0 = -(-x1 * (-x1 - 1 / eps).exp()).lambertw(-1)
            hit = eps * (x0 - 1)
            return hit**2 / 2

        late, t0, t1 = arb(2), ball("5/4"), ball("3/2")
        upper_late = upper_square(late)
        a0 = upper_late - upper_square(t0)
        a1 = upper_late - upper_square(t1)
        assert 0 < a1 and a1 < a0
        error = (a1 / a0).log() + (t1 - t0) / eps
        assert abs(error) < ball('1/100')
        print(f"eps=1/{n}: log(A(3/2)/A(5/4))+1/(4 eps)={error.str(9)}")


def valley_tests():
    print("\nStride-four family, d=10, eta=2^-20, q=2^-10")
    degree, stride = 10, 4
    eta = Fraction(1, 2**20)
    q = Fraction(1, 2**10)
    even = {degree: 1, degree + stride: eta - 2, degree + 2 * stride: 1}
    odd = {degree - 1: 1, degree + stride - 1: -2 - q, degree + 2 * stride - 1: 1}
    slope = {p - 2: p * c for p, c in even.items()}
    u = action_terms(slope, odd)
    lead = degree + 2 * stride
    middle = ball((degree + stride) * (2 - eta))
    disc = (middle**2 - 4 * lead * degree).sqrt()
    crest = ((middle - disc) / (2 * lead)).root(stride)
    minimum = ((middle + disc) / (2 * lead)).root(stride)
    m = polynomial(even, minimum)
    lo, hi = arb(0), ball('1/2')
    assert hi < crest and polynomial(even, hi) > m
    for _ in range(130):
        mid = (lo + hi) / 2
        if polynomial(even, mid) > m:
            hi = mid
        elif polynomial(even, mid) < m:
            lo = mid
        else:
            raise ArithmeticError("Bisection needs more precision")
    a = lo.union(hi)
    after = minimum + ball(Fraction(1, 2**30))
    edge = polynomial(u, a) + polynomial(slope, a) * polynomial(odd, minimum)
    h0 = polynomial(even, a / 2)
    h1 = polynomial(even, after)
    h2 = polynomial(even, arb(2))
    assert h0 < m and m < h1 and h1 < polynomial(even, crest) < h2
    shift = edge - polynomial(u, minimum)
    vo = [polynomial(u, a / 2), edge, polynomial(u, arb(2)) + shift]
    vi = [vo[0], edge + polynomial(u, after) - polynomial(u, minimum), vo[2]]
    for label, values in [('Vo', vo), ('Vi', vi)]:
        assert [sign(v) for v in values] == ['+', '-', '+']
        print(label, ' '.join(v.str(7) for v in values))

    # Enclose the integral over each full cell, including changes of sign.
    # For r>=2, E>=r^(d+2p)/4, and |rQ|<=3r^(d+2p).
    depths = [ball(eta**2), ball(eta), arb(1)]
    integrals = [arb(0) for _ in depths]
    smaller_odd = arb(0)
    extra_power = arb(0)
    panels = 2**16
    width = arb(2) / panels
    eta_ball, q_ball = ball(eta), ball(q)
    for i in range(panels):
        left, right = i * width, (i + 1) * width
        assert right <= 1 or left >= 1
        r = left.union(right)
        z = r**stride
        # Monotonicity on either side of 1 encloses the square without cancellation.
        shape = ((1 - left**stride)**2).union((1 - right**stride)**2)
        rd = r**degree
        e = rd * (shape + eta_ball * z)
        numerator = rd * (shape - q_ball * z)
        for j, h in enumerate(depths):
            assert h + e > 0
            integrals[j] += width * numerator / (h + e)**2
        middle_weight = width * rd / (eta_ball + e)**2
        smaller_odd += middle_weight * (shape - eta_ball * z)
        extra_power += middle_weight * (shape - q_ball * z + 2 * q_ball * r**2)
    tail = ball(Fraction(48, (degree + 2 * stride - 1) * 2**(degree + 2 * stride - 1)))
    assert 0 < tail < ball('0.000022')
    published = [("6900000000", "7100000000"), ("-140000", "-100000"), ("0.044", "0.045")]
    for j, h in enumerate(depths):
        integrals[j] += (-tail).union(tail)
        low, high = published[j]
        assert ball(low) < integrals[j] < ball(high)
        assert sign(integrals[j]) == ['+', '-', '+'][j]
        print(f"H={h.str(8)}: J lower {integrals[j].lower().str(10)}, "
              f"upper {integrals[j].upper().str(10)} in ({low}, {high})")
    print(f"{panels} whole-cell enclosures on [0,2], absolute tail <= {tail.str(7)}")
    for label, value in [('q=eta', smaller_odd), ('one additional r^11', extra_power)]:
        enclosure = value + (-tail).union(tail)
        assert enclosure > 0
        print(f"H=eta, {label}: J enclosure {enclosure.str(7)}")
    for scale, expected in [(2, '-'), (8, '+')]:
        at_one = eta - Fraction(18, 26 * scale**8)
        assert sign(ball(at_one)) == expected
        print(f"R={scale}, even valley value at r=1 after adding the wall: {at_one}")


if __name__ == "__main__":
    if not __debug__:
        raise SystemExit("Assertions must be enabled.")
    ctx.prec = 512
    print("python-flint", version("python-flint"), "Arb precision", ctx.prec)
    fold_directions()
    exact_crest_hits()
    valley_tests()
    print("\nPASS independent finite checks")
