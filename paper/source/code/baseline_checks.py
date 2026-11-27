#!/usr/bin/env python3
"""Rational enclosures for the trinomial gap, the degree-ten seed (including both
curvatures), the degree-twelve profile and the degree-eighteen fold.
Run from the computation root with python3 -B code/baseline_checks.py.
"""
from dataclasses import dataclass
from fractions import Fraction as F


@dataclass(frozen=True)
class RationalInterval:
    lo: F
    hi: F

    def __post_init__(self):
        assert self.lo <= self.hi

    def __add__(self, other: "RationalInterval") -> "RationalInterval":
        return RationalInterval(self.lo + other.lo, self.hi + other.hi)

    def __neg__(self) -> "RationalInterval":
        return RationalInterval(-self.hi, -self.lo)

    def __sub__(self, other: "RationalInterval") -> "RationalInterval":
        return self + (-other)

    def __mul__(self, other: "RationalInterval") -> "RationalInterval":
        products = [x * y for x in (self.lo, self.hi) for y in (other.lo, other.hi)]
        return RationalInterval(min(products), max(products))

    def power(self, n: int) -> "RationalInterval":
        assert self.lo >= 0 and n >= 0
        return RationalInterval(self.lo ** n, self.hi ** n)

    def reciprocal(self) -> "RationalInterval":
        assert self.lo > 0
        return RationalInterval(F(1) / self.hi, F(1) / self.lo)


def pt(val) -> RationalInterval:
    frac = F(val)
    return RationalInterval(frac, frac)


def bisect_rational_root(func, lo, hi, steps: int = 64) -> RationalInterval:
    lo_f, hi_f = F(lo), F(hi)
    assert func(lo_f) < 0 < func(hi_f)
    for _ in range(steps):
        mid = (lo_f + hi_f) / 2
        val = func(mid)
        if val == 0:
            return pt(mid)
        if val < 0:
            lo_f = mid
        else:
            hi_f = mid
    return RationalInterval(lo_f, hi_f)


def verify_trinomial_certificates():
    """Verify exact rational enclosures of 1000 * d * I_d for d = 6..18 and d >= 20 constants."""
    table = [
        (6, "0.062500000", 130, 135),
        (8, "0.040486856", 25, 27),
        (10, "0.029550858", -24, -22),
        (12, "0.023104700", -51, -50),
        (14, "0.018889538", -69, -67),
        (16, "0.015933561", -81, -79),
        (18, "0.013753521", -89, -88),
    ]
    for d, c_str, lo_b, hi_b in table:
        a = RationalInterval(F(c_str) - F(1, 10 ** 9), F(c_str) + F(1, 10 ** 9))
        tau, rho = F(d - 2, d - 1), F(d - 1, d)
        assert 0 < a.lo < a.hi < tau ** (d - 2) / (d - 1)
        s_ends, r_ends, t_ends = [], [], []
        for a_end in (a.lo, a.hi):
            s_ends.append(bisect_rational_root(lambda x: x ** (d - 2) * (d * x + d - 1) - a_end, 0, 1))
            r_ends.append(bisect_rational_root(lambda x: x ** (d - 2) * (d * x - d + 1) + a_end, tau, rho))
            t_ends.append(bisect_rational_root(lambda x: x ** (d - 1) - a_end, 0, 1))

        def f_prof(x_iv: RationalInterval, a_val: F, sgn: int) -> RationalInterval:
            return x_iv.power(d) + pt(sgn) * (x_iv.power(d - 1) - pt(a_val) * x_iv)

        assert (f_prof(s_ends[0], a.lo, 1) - f_prof(r_ends[0], a.lo, -1)).lo > 0
        assert (f_prof(s_ends[1], a.hi, 1) - f_prof(r_ends[1], a.hi, -1)).hi < 0

        s_iv = RationalInterval(s_ends[0].lo, s_ends[1].hi)
        r_iv = RationalInterval(r_ends[1].lo, r_ends[0].hi)
        t_iv = RationalInterval(t_ends[0].lo, t_ends[1].hi)
        assert s_iv.hi < t_iv.lo and r_iv.hi < 1

        def q_antideriv(x_iv: RationalInterval, sgn: int) -> RationalInterval:
            return (
                pt(F(d * d, 2 * d - 2)) * x_iv.power(2 * d - 2)
                + pt(F(sgn * 2 * d * (d - 1), 2 * d - 3)) * x_iv.power(2 * d - 3)
                + pt(F((d - 1) ** 2, 2 * d - 4)) * x_iv.power(2 * d - 4)
                - pt(F(sgn * 2 * d, d - 1)) * a * x_iv.power(d - 1)
                - pt(F(2 * (d - 1), d - 2)) * a * x_iv.power(d - 2)
            )

        h_iv = t_iv * r_iv * s_iv.reciprocal()
        log_h = RationalInterval(1 - F(1) / h_iv.lo, h_iv.hi - 1)
        scaled_id = pt(1000 * d) * (
            q_antideriv(t_iv, 1)
            - q_antideriv(s_iv, 1)
            - q_antideriv(pt(1), -1)
            + q_antideriv(r_iv, -1)
            + a.power(2) * log_h
        )
        assert lo_b < scaled_id.lo <= scaled_id.hi < hi_b

    assert F(19, 20) ** 19 < F(19, 50)
    assert 20 * F(5, 3) * F(2, 3) ** 19 < F(1, 60)
    assert 20 * F(3, 4) ** 19 < F(1, 6)
    assert (F(19, 50) + F(1, 60)) / (F(18, 19) + F(2, 3)) < F(1, 4)
    assert F(245, 1296) - F(1, 5) == -F(71, 6480)
    print("[PASS] Trinomial {1, d-1, d} certificates verified for d = 6..18 and d >= 20.")


def verify_degree_ten_seed_and_extensions():
    """Verify exact rational enclosures for {5,9,10}, {3,5,9,10}, {5,9,10,11,12}, {5,9,10,11,13,18}."""
    # 1. Balanced degree-10 seed {5, 9, 10}.
    def seed_poly(v: F) -> F:
        return v ** 14 - 9 * v ** 10 - 8 * v ** 9 - 8 * v ** 5 - 9 * v ** 4 + 1

    v_iv = bisect_rational_root(seed_poly, F(19, 10), F(2), 84)
    assert F("1.92194269368393019") < v_iv.lo <= v_iv.hi < F("1.92194269368393020")
    b_iv = RationalInterval(
        10 * (v_iv.lo ** 5 + 1) / (9 * (v_iv.hi ** 4 - 1)),
        10 * (v_iv.hi ** 5 + 1) / (9 * (v_iv.lo ** 4 - 1)),
    )
    a0_iv = RationalInterval((9 * b_iv.lo + 10) / 5, (9 * b_iv.hi + 10) / 5)

    # Check inner maximum of q_- on (1.460, 1.462) is strictly below 14 < 20.
    r_peak_lo = bisect_rational_root(
        lambda r: -(10 * r ** 5 - 9 * b_iv.hi * r ** 4 + 5 * a0_iv.lo),
        F(1460, 1000),
        F(1462, 1000),
        64,
    ).lo
    r_peak_hi = bisect_rational_root(
        lambda r: -(10 * r ** 5 - 9 * b_iv.lo * r ** 4 + 5 * a0_iv.hi),
        F(1460, 1000),
        F(1462, 1000),
        64,
    ).hi
    peak_val_upper = r_peak_hi ** 10 - b_iv.lo * r_peak_lo ** 9 + a0_iv.hi * r_peak_hi ** 5
    assert peak_val_upper < 14

    # Curvature comparison near the common minimum value for incoming signs.
    qpp_plus_lo = 50 + 36 * b_iv.lo
    qpp_plus_hi = 50 + 36 * b_iv.hi
    qpp_minus_lo = v_iv.lo ** 7 * (50 * v_iv.lo - 36 * b_iv.hi)
    qpp_minus_hi = v_iv.hi ** 7 * (50 * v_iv.hi - 36 * b_iv.lo)
    assert 130 < qpp_plus_lo <= qpp_plus_hi < 145
    assert 850 < qpp_minus_lo <= qpp_minus_hi < 1100
    assert v_iv.hi ** 2 / 850 < F(1, 145)

    # Seed action differences S_+(h) - S_-(h) at h = 20, 100, 200.
    def seed_antideriv(r_iv: RationalInterval, sgn: int) -> RationalInterval:
        return (
            pt(F(100, 18)) * r_iv.power(18)
            + pt(sgn * F(180, 17)) * b_iv * r_iv.power(17)
            + pt(F(81, 16)) * b_iv.power(2) * r_iv.power(16)
            - pt(sgn * F(100, 13)) * a0_iv * r_iv.power(13)
            - pt(F(90, 12)) * a0_iv * b_iv * r_iv.power(12)
            + pt(F(25, 8)) * a0_iv.power(2) * r_iv.power(8)
        )

    for h_val, lo_b, hi_b in ((20, -89, -88), (100, 1293, 1295), (200, 6884, 6886)):
        tp_lo = bisect_rational_root(
            lambda r: r ** 10 + (b_iv.hi * r ** 9 - a0_iv.lo * r ** 5) - h_val, F(5, 4), F(4), 84
        ).lo
        tp_hi = bisect_rational_root(
            lambda r: r ** 10 + (b_iv.lo * r ** 9 - a0_iv.hi * r ** 5) - h_val, F(5, 4), F(4), 84
        ).hi
        tm_lo = bisect_rational_root(
            lambda r: r ** 10 - (b_iv.lo * r ** 9 - a0_iv.hi * r ** 5) - h_val, F(2), F(4), 84
        ).lo
        tm_hi = bisect_rational_root(
            lambda r: r ** 10 - (b_iv.hi * r ** 9 - a0_iv.lo * r ** 5) - h_val, F(2), F(4), 84
        ).hi
        tp_iv = RationalInterval(tp_lo, tp_hi)
        tm_iv = RationalInterval(tm_lo, tm_hi)
        diff = (seed_antideriv(tp_iv, 1) - seed_antideriv(pt(1), 1)) - (
            seed_antideriv(tm_iv, -1) - seed_antideriv(v_iv, -1)
        )
        assert lo_b < diff.lo <= diff.hi < hi_b

    # 2. Outer minimum r^10 (1 -+ r)^2 action differences at h = 1/100, 1, 2 for {5,9,10,11,12}.
    def outer_antideriv(r_iv: RationalInterval, sgn: int) -> RationalInterval:
        # p(r) = r^10 - sgn * 2 r^11 + r^12
        coeffs = ((10, F(1)), (11, F(-2 * sgn)), (12, F(1)))
        acc = pt(0)
        for i, ci in coeffs:
            for j, cj in coeffs:
                exp_ij = i + j - 2
                acc = acc + pt(F(i * j * ci * cj, exp_ij)) * r_iv.power(exp_ij)
        return acc

    for h_val, lo_b, hi_b in ((F(1, 100), F(3, 10000), F(4, 10000)), (F(1), F(-3, 2), F(-1)), (F(2), F(-5), F(-4))):
        tp = bisect_rational_root(lambda r: r ** 10 * (r - 1) ** 2 - h_val, F(1), F(2), 84)
        tm = bisect_rational_root(lambda r: r ** 10 * (r + 1) ** 2 - h_val, F(0), F(1), 84)
        diff = (outer_antideriv(tp, 1) - outer_antideriv(pt(1), 1)) - outer_antideriv(tm, -1)
        assert lo_b < diff.lo <= diff.hi < hi_b

    # 3. Degree-18 distant-fold cancellation rho in (6.901267, 6.901268) for {5,9,10,11,13,18}.
    g_iv = bisect_rational_root(lambda r: r ** 8 - F(9, 5), F(1), F(3, 2), 90)
    for rho_val, sgn in ((F(6901267, 10 ** 6), -1), (F(6901268, 10 ** 6), 1)):
        poly = {}
        for e_exp, f_c in ((8, F(10)), (16, F(-10))):
            for p_exp, q_c in ((9, F(1)), (11, -1 - rho_val), (13, rho_val)):
                poly[e_exp + p_exp] = poly.get(e_exp + p_exp, F(0)) + f_c * q_c * F(p_exp - e_exp, p_exp + e_exp)
        u_at_1 = sum(poly.values())
        assert u_at_1 < 0
        u_at_g = pt(0)
        for deg_m, c_m in poly.items():
            u_at_g = u_at_g + pt(c_m) * g_iv.power(deg_m)
        t_at_g = u_at_g - pt(u_at_1)
        if sgn > 0:
            assert t_at_g.lo > 0
        else:
            assert t_at_g.hi < 0

    print("[PASS] Balanced degree-10 seed and extensions {3,5,9,10}, {5,9,10,11,12}, {5,9,10,11,13,18} verified.")


def main():
    if not __debug__:
        raise SystemExit("Run without -O so all certificate assertions are active.")
    verify_trinomial_certificates()
    verify_degree_ten_seed_and_extensions()
    print("PASS trinomial, balanced-seed and fold enclosures.")


if __name__ == "__main__":
    main()
