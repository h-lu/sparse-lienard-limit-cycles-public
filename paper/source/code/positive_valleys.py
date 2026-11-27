#!/usr/bin/env python3
"""Rigorous rational and Arb ball-arithmetic certificates for Theorem 1.1.

This script verifies all finite-degree algebraic and integral enclosures
supporting the positive-valley sparse Lienard construction:
  1. Exact rational signs of the limiting Wronskian W_D(z) at the six positive
     nodes z in {1/8192, 1/8, 1/2, 3/2, 16, 8192} for stride D = 258, and the
     three exponential base inequalities for incoming heights -H at depths H = w^n, with
     w in {1/2, 2/3, 7/8} with k_0 = 2/3, s = 3/4.
  2. Exact rational root enclosures for all four initial odd-prefix classes
     j in {0, 1, 2, 3} (fractional offset gamma = 2j/258) at the six positive
     nodes and five depth nodes.
  3. Validated 768-bit Arb enclosures for the balanced compact seed on
     {m_0, 2m_0 - 1, 2m_0} at m_0 = 33153, d_0 = 66306:
       0.81 < H_b / d_0 < 0.82,
       4 < S_+(d_0) / m_0^3 < 5,
       16 < S_-(d_0) / m_0^3 < 17,
       q_-''(v) / (v^2 q_+''(1)) > 10^5.
  4. Validated 768-bit Arb enclosures for the isolated positive-valley group at
     order n = 257 (d = 66306, eta = (2/3)^257, q = (3/4)^128):
       - enclosures of the local maximum z_p, local minimum z_c, and pre-crest
         equal-height inverse z_a = a^D;
       - strict signs [+, -, +, -, +, -, +] for both positive-height variation
         lists V_o and V_i, including the moving-minimum term f(a)Q(c) and the
         post-valley incoming test z_v = z_c + 2^-200;
       - independent central-difference cross-check at perturbation tau = 2^-240
         by bisecting the perturbed extrema and level inverses and integrating
         the exact squared radial slopes directly;
       - rigorous signed Riemann-sum + analytic-tail enclosures of n H D J(H)
         over (0, infinity) at H = w^257 for w in {1/2, 2/3, 7/8, 3/2, 4},
         certifying the alternating signs [+, -, +, -, +].
  5. Exact support collection, parity separation, degree formula, odd-power
     count b = 2L + j + 4, explicit ten-power support
       {33153, 66305, 66306, 66563, 66822, 67079, 67080, 67081, 67083, 67088},
     and the count arithmetic 16L + 10 + 3j = 4k - 14 - j.

From the supplied source directory or the standalone computation root, run
  python3 -B code/positive_valleys.py
"""

from fractions import Fraction as F
from flint import arb, ctx


STRIDE = 258
ORDER = 257
PRE_NODES = tuple(map(F, ("1/8192", "1/8", "1/2")))
POST_NODES = tuple(map(F, ("3/2", "16", "8192")))
DEPTH_NODES = tuple(map(F, ("1/2", "2/3", "7/8", "3/2", "4")))


def ball(rational) -> arb:
    """Convert an exact Fraction, int, or decimal/rational string to an Arb ball."""
    return arb(str(rational))


def certified_sign(value: arb) -> int:
    """Return +1 or -1 when an Arb ball lies strictly on one side of zero."""
    if value > 0:
        return 1
    if value < 0:
        return -1
    raise ArithmeticError(f"Unresolved sign in interval enclosure: {value}")


def horner(coefficients, x):
    """Evaluate a polynomial given in ascending power order."""
    acc = 0
    for c in reversed(coefficients):
        acc = acc * x + c
    return acc


def action_polynomial_coeffs(slope_coeffs, odd_coeffs, degree: int):
    """Exact rational coefficients in z = r^D of U(r) / (d * r^(2d - 3))."""
    out = [F(0)] * (len(slope_coeffs) + len(odd_coeffs) - 1)
    for i, left in enumerate(slope_coeffs):
        for j, right in enumerate(odd_coeffs):
            out[i + j] += (
                left * right * (1 + STRIDE * (j - i))
                / (2 * degree - 3 + STRIDE * (i + j))
            )
    return out


def verify_limit_signs_and_exponents():
    """Verify rational limit signs of W_D(z) and base inequalities for incoming heights -H."""
    expected = {
        PRE_NODES[0]: 1, PRE_NODES[1]: -1, PRE_NODES[2]: 1,
        POST_NODES[0]: 1, POST_NODES[1]: -1, POST_NODES[2]: 1,
    }
    for z in PRE_NODES + POST_NODES:
        c_val = (1 - z) ** 2 * (1 + 2 * z)
        b_val = (z - 1) * (z * z + z - 8) / 8
        wronskian_z = c_val * (3 * z * z - 9) / 8 - 6 * z * (z - 1) * b_val
        assert wronskian_z == -(z - 1) ** 2 * (3 * z * z - 30 * z + 9) / 8
        w_d = c_val * b_val + STRIDE * z * wronskian_z
        assert (1 if w_d > 0 else -1) == expected[z]

    w_small, k0, w_recov = DEPTH_NODES[:3]
    s = F(3, 4)
    # Check ((1/2)*sqrt(s)/k_0^(3/2))^2 < 1, s/k_0 > 1, s/(7/8) < 1.
    assert w_small ** 2 * s < k0 ** 3
    assert k0 < s < w_recov
    print("[PASS] Limiting W_D(z) signs and exponential bases for negative incoming heights verified.")


def enclose_fractional_power(z: F, prefix_j: int):
    """Enclose z^(-2*prefix_j / D) strictly inside a rational dyadic interval."""
    if prefix_j == 0:
        return F(1), F(1)
    target = z ** (-2 * prefix_j)
    lo, hi = F(1, 2), F(2)
    assert lo ** STRIDE < target < hi ** STRIDE
    for _ in range(50):
        mid = (lo + hi) / 2
        if mid ** STRIDE < target:
            lo = mid
        else:
            hi = mid
    return lo, hi


def verify_prefix_residue_classes():
    """Verify exact rational sign margins for all four prefix lengths j = 0,1,2,3."""
    wanted_pos = (1, -1, 1, 1, -1, 1)
    margins = (F(9, 10), F(10), F(20), F(300), F(10 ** 7), F(10 ** 22))
    wanted_depth = (1, 1, 1, -1, 1)
    for j in range(4):
        gamma = F(2 * j, STRIDE)
        for z, sgn, margin in zip(PRE_NODES + POST_NODES, wanted_pos, margins):
            c_val = (1 - z) ** 2 * (1 + 2 * z)
            dc_val = 6 * z * (z - 1)
            for factor in enclose_fractional_power(z, j):
                b_val = 1 - F(9, 8) * z + z ** 3 * factor / 8
                db_val = -F(9, 8) + (3 - gamma) * z * z * factor / 8
                w_val = (1 + 2 * j) * c_val * b_val + STRIDE * z * (c_val * db_val - dc_val * b_val)
                assert sgn * w_val > margin
        for z, sgn in zip(DEPTH_NODES, wanted_depth):
            for factor in enclose_fractional_power(z, j):
                b_val = 1 - F(9, 8) * z + z ** 3 * factor / 8
                assert sgn * b_val > 0
    print("[PASS] All four odd-prefix classes j in {0, 1, 2, 3} pass exact rational checks.")


def bisect_crossing(func, lo: F, hi: F, steps: int) -> arb:
    """Enclose a simple sign change from negative to positive using exact dyadic midpoints."""
    assert certified_sign(func(ball(lo))) == -1
    assert certified_sign(func(ball(hi))) == 1
    for _ in range(steps):
        mid = (lo + hi) / 2
        if certified_sign(func(ball(mid))) < 0:
            lo = mid
        else:
            hi = mid
    return ball(lo).union(ball(hi))


def verify_compact_seed_33153():
    """Certify the compact seed {m_0, 2m_0 - 1, 2m_0} at m_0 = 33153, d_0 = 66306."""
    d0 = STRIDE * ORDER
    m0 = d0 // 2
    assert m0 == 33153 and d0 == 66306 and m0 % 2 == 1

    def seed_defining_poly(t: arb) -> arb:
        r = t.root(m0)
        return (
            t ** 3 / r
            - (d0 - 1) * (t * t + t / r)
            - 2 * (m0 - 1) * (t * t / r + t)
            + 1
        )

    tv = bisect_crossing(seed_defining_poly, F(2 * m0), F(8 * m0), 430)
    v = tv.root(m0)
    b0 = d0 * (tv + 1) / ((d0 - 1) * (tv / v - 1))
    a0 = 2 + (d0 - 1) * b0 / m0

    def profile_in_t(t: arb, side: int) -> arb:
        return t * t + side * (b0 * t * t / t.root(m0) - a0 * t)

    def minus_reduced_deriv(t: arb) -> arb:
        return d0 * t - (d0 - 1) * b0 * t / t.root(m0) + m0 * a0

    peak_t = bisect_crossing(lambda t: -minus_reduced_deriv(t), F(1), F(2 * m0), 240)
    height_b = profile_in_t(peak_t, -1)
    assert height_b < d0 and peak_t < tv

    turn_plus = bisect_crossing(lambda t: profile_in_t(t, 1) - d0, F(1), F(d0), 240).root(m0)
    turn_minus = bisect_crossing(lambda t: profile_in_t(t, -1) - d0, F(2 * m0), F(16 * m0), 240).root(m0)
    assert turn_plus > 1 and turn_minus > v

    def exact_slope_integral(side: int, start: arb, end: arb) -> arb:
        terms = ((d0, arb(1)), (d0 - 1, side * b0), (m0, -side * a0))
        total = arb(0)
        for deg_i, ci in terms:
            for deg_j, cj in terms:
                exp_ij = deg_i + deg_j - 2
                total += ci * cj * deg_i * deg_j * (end ** exp_ij - start ** exp_ij) / exp_ij
        return total

    s_plus = exact_slope_integral(1, arb(1), turn_plus)
    s_minus = exact_slope_integral(-1, v, turn_minus)
    curv_plus = 2 * m0 * m0 + (d0 - 1) * (m0 - 1) * b0
    curv_minus = v ** (d0 - 3) * (2 * m0 * m0 * v - (d0 - 1) * (m0 - 1) * b0)
    curv_ratio = curv_minus / (v * v * curv_plus)

    assert ball("81/100") < height_b / d0 < ball("82/100")
    assert 4 < s_plus / m0 ** 3 < 5
    assert 16 < s_minus / m0 ** 3 < 17
    assert curv_ratio > 10 ** 5
    print(
        f"[PASS] Compact seed m_0={m0}: H_b/d_0 in {(height_b / d0).str(8)}, "
        f"S_+/m_0^3 in {(s_plus / m0 ** 3).str(8)}, "
        f"S_-/m_0^3 in {(s_minus / m0 ** 3).str(8)}, "
        f"curvature ratio in {curv_ratio.str(8)}."
    )


def even_shape(z, eta):
    return (1 - z) ** 2 * (1 + 2 * z) + eta * z * z


def odd_shape(z, q):
    return (z - 1) * (z * z + z - 8) / 8 - q * z


def verify_finite_valley_positive_lists(eta_frac: F, q_frac: F):
    """Certify both 7-sign positive lists at n = 257 and cross-check by perturbed endpoints."""
    n = ORDER
    d = n * STRIDE
    eta, q = ball(eta_frac), ball(q_frac)
    even_c = [F(1), F(0), eta_frac - 3, F(2)]
    slope_c = [c * F(n + i, n) for i, c in enumerate(even_c)]
    odd_c = [F(1), -F(9, 8) - q_frac, F(0), F(1, 8)]
    act_coeffs = list(map(ball, action_polynomial_coeffs(slope_c, odd_c, d)))

    def energy(z: arb) -> arb:
        return z ** n * even_shape(z, eta)

    def slope_factor(z: arb) -> arb:
        return even_shape(z, eta) + z * z * (6 * (z - 1) + 2 * eta) / n

    def action_u(z: arb) -> arb:
        return d * z ** (2 * n) * horner(act_coeffs, z) / z.root(STRIDE) ** 3

    def slope_f(z: arb) -> arb:
        return d * z ** n * slope_factor(z) / z.root(STRIDE) ** 2

    def odd_q(z: arb) -> arb:
        return z ** n * odd_shape(z, q) / z.root(STRIDE)

    z_crest = bisect_crossing(lambda z: -slope_factor(z), F(1, 2), F(n - 1, n), 320)
    z_valley = bisect_crossing(slope_factor, F(n - 1, n), F(1), 380)
    m_val = energy(z_valley)
    crest_val = energy(z_crest)
    assert ball("7/8") < z_crest
    assert ball("99229/100000") < z_crest < ball("99230/100000")
    assert ball("239/100") * arb(10) ** -44 < 1 - z_valley < ball("240/100") * arb(10) ** -44
    assert ball("555/100") * arb(10) ** -46 < m_val < ball("556/100") * arb(10) ** -46
    assert ball("242/100") * arb(10) ** -5 < crest_val < ball("243/100") * arb(10) ** -5
    z_pre = bisect_crossing(lambda z: energy(z) - m_val, F(1, 2), F(7, 8), 240)
    assert ball("67022/100000") < z_pre < ball("67023/100000")

    v_ord_mid = action_u(z_pre) + slope_f(z_pre) * odd_q(z_valley)
    assert certified_sign(v_ord_mid) == -1

    z_post_valley = z_valley + arb(2) ** -200
    assert m_val < energy(z_post_valley) < crest_val
    v_inc_mid = v_ord_mid + action_u(z_post_valley) - action_u(z_valley)
    assert certified_sign(v_inc_mid) == -1

    assert all(energy(ball(z)) < m_val for z in PRE_NODES)
    assert all(energy(ball(z)) > crest_val for z in POST_NODES)

    vals_pre = [action_u(ball(z)) for z in PRE_NODES]
    vals_post = [action_u(ball(z)) + v_ord_mid - action_u(z_valley) for z in POST_NODES]
    pos_boxes = (
        (vals_pre[0], "176/100", "177/100", -2012),
        (vals_pre[1], "-511/100", "-510/100", -464),
        (vals_pre[2], "207/100", "208/100", -154),
        (v_ord_mid, "-371/100", "-370/100", -57),
        (v_inc_mid, "-371/100", "-370/100", -57),
        (vals_post[0], "571/100", "572/100", 92),
        (vals_post[1], "-122/100", "-121/100", 626),
        (vals_post[2], "984/100", "985/100", 2033),
    )
    for val, lo_s, hi_s, exp10 in pos_boxes:
        scale = arb(10) ** exp10
        assert ball(lo_s) * scale < val < ball(hi_s) * scale
    signs_pre = [certified_sign(val) for val in vals_pre]
    signs_post = [certified_sign(val) for val in vals_post]
    assert signs_pre + [certified_sign(v_ord_mid)] + signs_post == [1, -1, 1, -1, 1, -1, 1]
    assert signs_pre + [certified_sign(v_inc_mid)] + signs_post == [1, -1, 1, -1, 1, -1, 1]

    # Independent cross-check of V_o = U(a) + f(a)Q(c) by perturbing both profiles
    # at tau = 2^-240, relocating the perturbed minimum and pre-crest inverse, and
    # integrating the squared radial slope directly.
    tau = arb(2) ** -240

    def perturbed_plateau_integral(side: int) -> arb:
        def p_energy(z: arb) -> arb:
            return energy(z) + side * tau * odd_q(z)

        def p_deriv_factor(z: arb) -> arb:
            db = (3 * z * z - 9) / 8 - q
            odd_part = ((n - arb(1) / STRIDE) * odd_shape(z, q) + z * db) / z.root(STRIDE)
            return n * slope_factor(z) + side * tau * odd_part

        p_bottom = bisect_crossing(p_deriv_factor, F(n - 1, n), F(1), 500)
        p_min = p_energy(p_bottom)
        p_early = bisect_crossing(lambda z: p_energy(z) - p_min, F(1, 2), F(7, 8), 380)
        rad = p_early.root(STRIDE)
        terms = [(d + i * STRIDE, ball(c)) for i, c in enumerate(even_c) if c]
        terms.extend((d - 1 + i * STRIDE, side * tau * ball(c)) for i, c in enumerate(odd_c) if c)
        return sum(
            ci * cj * deg_i * deg_j * rad ** (deg_i + deg_j - 2) / (deg_i + deg_j - 2)
            for deg_i, ci in terms
            for deg_j, cj in terms
        )

    diff_quotient = (perturbed_plateau_integral(1) - perturbed_plateau_integral(-1)) / (2 * tau)
    assert diff_quotient < 0
    assert ball("999/1000") < diff_quotient / v_ord_mid < ball("1001/1000")
    print(f"[PASS] n={n}: both 7-sign positive variation lists and direct perturbed cross-check verified.")


def dyadic_valley_panels():
    """Partition [1/8, 8] into dyadic subintervals graded geometrically toward z = 1."""
    for left, right in (
        (F(1, 8), F(1, 4)),
        (F(1, 4), F(1, 2)),
        (F(3, 2), F(2)),
        (F(2), F(4)),
        (F(4), F(8)),
    ):
        yield left, right, 1024
    for k in range(1, 100):
        far, near = F(1, 2 ** k), F(1, 2 ** (k + 1))
        yield 1 - far, 1 - near, 512
        yield 1 + near, 1 + far, 512
    rad = F(1, 2 ** 100)
    yield 1 - rad, F(1), 1
    yield F(1), 1 + rad, 1


def verify_finite_valley_depth_integrals(eta_frac: F, q_frac: F):
    """Enclose n H D J(H) on (0, infinity) at all five depth nodes for n = 257."""
    n = ORDER
    eta, q = ball(eta_frac), ball(q_frac)
    heights = [ball(w) ** n for w in DEPTH_NODES]
    accum = [arb(0) for _ in heights]

    partition = sorted(dyadic_valley_panels())
    assert partition[0][0] == F(1, 8) and partition[-1][1] == F(8)
    assert all(p1[1] == p2[0] for p1, p2 in zip(partition, partition[1:]))

    assert sum(count for _, _, count in partition) == 106498
    for left, right, count in partition:
        assert right <= 1 or left >= 1
        step = ball((right - left) / count)
        a_pt = ball(left)
        for _ in range(count):
            b_pt = a_pt + step
            z_box = a_pt.union(b_pt)
            z_pow = z_box ** n
            # On each side of 1 the square is monotone; endpoints are dyadic.
            distance_squared = ((a_pt - 1) ** 2).union((b_pt - 1) ** 2)
            e_box = z_pow * (distance_squared * (1 + 2 * z_box) + eta * z_box**2)
            num_density = n * step * z_pow * z_box.root(STRIDE) * odd_shape(z_box, q) / z_box
            for idx, h_val in enumerate(heights):
                denom = h_val + e_box
                assert denom > 0
                accum[idx] += h_val * num_density / denom ** 2
            a_pt = b_pt

    # Add rigorous tail bounds on (0, 1/8) and (8, infinity):
    # On (0, 1/8): |B_q(z)| < 2, H + z^n C_eta(z) >= H, z^(1/D) < 1 -> <= (2/H) * 8^(-n).
    # On (8, inf): C_eta(z) >= z^3, |B_q(z)| <= z^3, z^(1/D) <= z -> <= H * 8^(-n).
    for idx, h_val in enumerate(heights):
        tail = (2 / h_val + h_val) * arb(8) ** (-n)
        accum[idx] += arb(0, tail.upper())

    assert [certified_sign(val) for val in accum] == [1, -1, 1, -1, 1]
    target_boxes = (
        ("7/10", "6/5"),
        ("-400000000", "-300000000"),
        ("12/5", "3"),
        ("-29/100", "-23/100"),
        ("1/25", "7/100"),
    )
    print(f"[PASS] n={n}: rigorous Arb enclosures of n H D J(H) on (0, infinity):")
    for w, enc, (lo_s, hi_s) in zip(DEPTH_NODES, accum, target_boxes):
        assert ball(lo_s) < enc < ball(hi_s)
        print(
            f"       w = {str(w):>3}: ({lo_s}, {hi_s}) contains "
            f"[{enc.lower().str(6)}, {enc.upper().str(6)}]"
        )


def verify_supports_and_cycle_ledgers():
    """Verify the support, parity and count arithmetic in the existence proof."""
    expected_ten = [33153, 66305, 66306, 66563, 66822, 67079, 67080, 67081, 67083, 67088]
    for groups_l in (1, 2, 3, 4, 8, 29):
        for j in range(4):
            d = d0 = STRIDE * ORDER
            support = {d0 // 2, d0 - 1, d0}
            support.update(d0 - 1 + 2 * i for i in range(1, j + 1))
            alpha_phys = F(1)
            beta_phys = F((-1) ** j)
            for block in range(groups_l):
                gap = 2 * j if block == 0 else 0
                scale_s = F(2 + block)
                n_step = ORDER + 3 * block
                eta_step = F(2, 3) ** n_step
                q_step = F(3, 4) ** (n_step // 2)
                c_even_2d = -(F(3) - eta_step) * alpha_phys / (scale_s ** (2 * STRIDE))
                c_even_3d = F(2) * alpha_phys / (scale_s ** (3 * STRIDE))
                c_odd_1d = -(F(9, 8) + q_step) * beta_phys / (scale_s ** STRIDE)
                c_odd_3d = beta_phys / (F(8) * (scale_s ** (3 * STRIDE - gap)))
                assert c_even_2d != 0 and c_even_3d > 0
                assert c_odd_1d != 0 and (1 if c_odd_3d > 0 else -1) == (-1) ** j
                added = {
                    d + 2 * STRIDE,
                    d + 3 * STRIDE,
                    d - 1 + gap + STRIDE,
                    d - 1 + 3 * STRIDE,
                }
                assert len(added) == 4 and added.isdisjoint(support)
                support.update(added)
                alpha_phys = c_even_3d
                beta_phys = c_odd_3d
                d += 3 * STRIDE
            assert max(support) == d and (d - 1) in support
            fold_terms = {d + 1, d + 3, d + 8}
            assert fold_terms.isdisjoint(support)
            support.update(fold_terms)

            k = 4 * groups_l + 6 + j
            assert len(support) == k
            assert max(support) == d0 + 3 * STRIDE * groups_l + 8
            even_powers = [p for p in support if p % 2 == 0]
            odd_powers = [p for p in support if p % 2 == 1]
            assert len(even_powers) == 2 * groups_l + 2
            assert len(odd_powers) == 2 * groups_l + 4 + j
            assert min(support) == 33153 and min(support) % 2 == 1

            ord_pos = 6 * groups_l + j
            rem_pos = 6 * groups_l + j
            rem_neg = 4 * groups_l + j
            total_cycles = (2 + ord_pos) + (rem_pos + 2) + (rem_neg + 2) + 4
            assert total_cycles == 16 * groups_l + 10 + 3 * j
            assert total_cycles == 4 * k - 14 - j

            if groups_l == 1 and j == 0:
                assert sorted(support) == expected_ten
                assert max(support) == 67088
                assert total_cycles == 26
                assert len(odd_powers) - 1 == 5

    print("[PASS] Collected supports S_{L,j}, parity counts, ten-power support, and cycle ledgers verified.")


def main():
    if not __debug__:
        raise SystemExit("Run without -O so all certificate assertions are active.")
    ctx.prec = 768
    verify_limit_signs_and_exponents()
    verify_prefix_residue_classes()
    verify_compact_seed_33153()
    eta_257 = F(2, 3) ** ORDER
    q_257 = F(3, 4) ** (ORDER // 2)
    verify_finite_valley_positive_lists(eta_257, q_257)
    verify_finite_valley_depth_integrals(eta_257, q_257)
    verify_supports_and_cycle_ledgers()
    print("ALL POSITIVE-VALLEY CERTIFICATES PASSED.")


if __name__ == "__main__":
    main()
