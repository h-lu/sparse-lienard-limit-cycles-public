"""Compare the numbers printed in the article and supplement with the finite input ledger."""
import json
import re
from fractions import Fraction
from pathlib import Path

from sympy import Rational, expand, solve, symbols

SOURCE = Path(__file__).resolve().parents[1]


def decimal_or_power(token):
    """Read a rational or a signed integer times a power of ten in Table 2."""
    token = token.replace(r"\,", "").strip()
    scaled = re.fullmatch(r"(-?\d+)\\cdot10\^\{?(\d+)\}?", token)
    if scaled:
        return Fraction(int(scaled[1]) * 10**int(scaled[2]))
    return Fraction(token)


def check_rational_margins():
    z, ell, xi = symbols("z ell xi")
    stride = 258
    c = (1 - z)**2 * (1 + 2*z)
    b8 = 8 - 9*z + xi*z**3
    # First differentiate with xi = z^(-ell/D), then bound ell and xi independently.
    derivative = b8.diff(z) - b8.diff(xi) * ell * xi / (stride*z)
    direct = (1 + ell)*c*b8 + stride*z*(c*derivative - c.diff(z)*b8)
    affine = ((8 - 9*z)*((1 + ell)*c - stride*z*c.diff(z)) - 9*stride*z*c
              + xi*z**3*(c + 3*stride*(1 - z**2)))
    assert expand(direct - affine) == 0
    assert Rational(5, 4)**43 > 2**13

    nodes = [Rational(x) for x in ("1/8192", "1/8", "1/2", "3/2", "16", "8192")]
    margins = [Rational(x) for x in ("9/10", "10", "20", "300", "10000000", str(10**22))]
    signs = [1, -1, 1, 1, -1, 1]
    corners, floors = [], []
    for at, margin, sign in zip(nodes, margins, signs, strict=True):
        bounds = [1, Rational(5, 4)] if at < 1 else [Rational(4, 5), 1]
        values = {(i, j): sign*direct.subs({z: at, ell: i, xi: j})/8
                  for i in (0, 6) for j in bounds}
        corner = min(values, key=values.get)
        corners.append(f"({corner[0]},{corner[1]})")
        floors.append(int(100*values[corner]//margin))

    latex = (SOURCE / "tex/valleys.tex").read_text()
    expected_rows = [
        ("z", ["2^{-13}", "1/8", "1/2", "3/2", "16", "2^{13}"]),
        ("s", [f"{item:+d}" for item in signs]),
        (r"(l,\xi)_{\min}", corners),
        ("m_z", ["9/10", "10", "20", "300", "10^7", "10^{22}"]),
        (r"\lfloor 100sW_*/m_z\rfloor", list(map(str, floors))),
    ]
    printed_rows = {line.strip().removesuffix(r"\\") for line in latex.splitlines()}
    for label, values in expected_rows:
        assert "&".join([label, *values]) in printed_rows
    assert all(value > 100 for value in floors)
    example = direct.subs({z: Rational(1, 8), ell: 6, xi: Rational(5, 4)})/8
    assert f"gives $W={example}$" in latex
    print("PASS Symbolic W identity and six exact rational rectangle minima; printed margin floors", *floors)


def check_rising_witness():
    """Differentiate the even primitive and integrate in z to recover the printed witness."""
    z = symbols("z")
    stride, h = 1024, Rational(1, 2**30)
    # With z=r^D write E=r^10 e(z); E'/r=10 r^8 C(z).
    e = 1 + 10*(h-3)*z**2/2058 + Rational(20, 3082)*z**3
    c = expand(e + stride*z*e.diff(z)/10)
    g = 1 + (h-Rational(4, 5))*z + Rational(256, 3375)*z**3
    assert expand(c - (1-z)**2*(1+2*z) - h*z**2) == 0
    # The unnormalized r-integrand is 10 r^16 times this polynomial.
    integrand = expand(c*g + stride*z*(c*g.diff(z)-c.diff(z)*g))
    primitive = sum(integrand.coeff(z, k)*z**k/(17+stride*k) for k in range(7))
    assert expand(17*primitive + stride*z*primitive.diff(z) - integrand) == 0

    example = (SOURCE / "tex/valleys.tex").read_text().split("Sharpness of the six-zero bound.", 1)[1]
    for term in (r"D=1024", r"h=2^{-30}", r"\frac{10(h-3)}{2058}r^{2058}",
                 r"\frac{20}{3082}r^{3082}", r"G(z)=1+(h-4/5)z+\frac{256}{3375}z^3"):
        assert term in example, term
    nodes = [Rational(x) for x in ("1/16", "1/6", "1", "15/8", "4", "32", "4096")]
    floors = [(10**6*primitive.subs(z, x)) // (1+x**6) for x in nodes]
    assert all((-1)**i * value > 0 for i, value in enumerate(floors))
    assert "z&" + "&".join(map(str, nodes)) + r"\\" in example
    assert r"\lfloor 10^6 S(z)/(1+z^6)\rfloor&" + "&".join(map(str, floors)) in example
    print("PASS Article 3.2 rising witness: exact normalized floors", *floors)


def check_lu_term_comparison():
    """Recover Lu's term-count comparisons, including the preclosure ordinary lists."""
    s = symbols("s", integer=True, positive=True)
    m = 8 + 4*s
    b_top = 3 + 2*s
    odd_top = 5 + 2*s
    # Integrating B^2 and the two z^8 controls against du=z dz raises the degrees by two.
    assert 2*b_top + 2 == m and 8 + 2*(2*s - 1) + 2 == m
    # The even exponents are 2 and 8,10,...,m. The odd ones are 3,5,...,5+2s.
    even_count = 2 + (m - 8)/2
    odd_count = 1 + (odd_top - 3)/2
    assert (m - 1 - odd_top).is_positive
    # The shortcut adjoins m-1, realizes with x, then uses the three exterior powers.
    terms = expand(even_count + odd_count + 1 + 1 + 3)
    cycles = 7*s + (3*s + s) + 1
    assert terms == 3*s + 9 and cycles == 11*s + 1
    assert cycles.coeff(s) / terms.coeff(s) == Rational(11, 3)
    for remainder in (0, 1, 2):
        k = terms + remainder
        assert expand((k - 9)/3 - s) == Rational(remainder, 3)
        assert expand(m + 3 - k).is_nonnegative

    # At one radius, rescale z by q1 in the moment integral of Lemma 5.3.
    # The odd height-splitting polynomial is -z^5 + c*z^7/q1^2.
    v = symbols("v")
    weight = v**6 * (1 - v**2)**2
    moment5 = (weight * (v**5).diff(v)).integrate((v, 0, 1))
    moment7 = (weight * (v**7).diff(v)).integrate((v, 0, 1))
    top = moment5 / moment7
    assert top == Rational(85, 77)
    m_one = 12
    b_powers = {3, 5}
    b_squared = {p + q for p in b_powers for q in b_powers}
    evens = {p + 2 for p in {0, 8, 10} | b_squared}
    support8 = {1} | evens | b_powers | {5, 7}
    support9 = support8 | {9}
    assert sorted(support8) == [1, 2, 3, 5, 7, 8, 10, 12] and len(support9) == 9

    # In Lu's depth convention the small launch is positive. The first shoulder
    # multiplies the eight cubic-window signs by s1=-1. In (3.8) the coefficient
    # divided by m*c is b_j*(4*j/(m+j-2)-2), including the moving endpoint.
    sample_signs = [1] + [(-1)**(i + 1) for i in range(8)]
    ordinary_counts = {}
    for support, j, odd_coefficient in [(support8, 7, top), (support9, 9, -1)]:
        coefficient = odd_coefficient * (Rational(4*j, m_one + j - 2) - 2)
        assert coefficient != 0
        tail_sign = 1 if coefficient > 0 else -1
        sample_signs.append(tail_sign)
        # Detune the higher matched reference only after fixing all far tests.
        chain = sample_signs + [-tail_sign]
        assert all(left == -right for left, right in zip(chain, chain[1:]))
        ordinary_counts[len(support)] = len(chain) - 1
    assert ordinary_counts == {8: 10, 9: 11}

    # Closing the nine-term profile as in Proposition 5.1 at e=0 gives degree 14.
    j = m_one - 3
    allowed = support9 | {m_one + 2, j + 2, j + 4}
    assert allowed == set(range(1, 15)) - {4, 6} and len(allowed) == 12
    ordinary = ordinary_counts[9]
    incoming = (3 + 2) + 1 + (1 + 2) + 1
    # The terminal sign is negative, whereas (2-g)b is positive for g=3, b<0.
    assert (ordinary, incoming, ordinary + incoming + 1 + 1) == (11, 10, 23)

    # Row e=1 retains the same even profile and adds the odd power m-1.
    j = m_one - 1
    allowed_e1 = support9 | {j, m_one + 2, j + 2, j + 4}
    assert allowed_e1 == set(range(1, 16)) - {4, 6} and len(allowed_e1) == 13
    ordinary_e1 = ordinary_counts[9] + 1
    # The extra positive-height change comes before the join to the negative list.
    # Here the original table takes zero translation and appends no extra depth.
    incoming_e1 = (3 + 2 + 1) + 1 + (1 + 2)
    assert (ordinary_e1, incoming_e1, ordinary_e1 + incoming_e1 + 2) == (12, 10, 24)

    article = (SOURCE / "main.tex").read_text()
    supplement = (SOURCE / "low-term-supplement.tex").read_text()
    for formula in (r"s=\lfloor(k-9)/3\rfloor", r"H(k)/k\ge11/3",
                    r"H(k)\ge11\lfloor(k-9)/3\rfloor+1", r"$4s+11$"):
        assert formula in article, formula
    assert expand(m + 3) == 4*s + 11
    for formula in (r"7s+4s+1=11s+1", r"(2s+2)+(s+2)+1+1+3=3s+9", r"m+3=4s+11",
                    r"\{1,\ldots,14\}\setminus\{4,6\}", r"H(12)\ge23",
                    r"\{1,\ldots,15\}\setminus\{4,6\}", r"H(13)\ge24"):
        assert formula in supplement, formula
    print("PASS At most 3s+9 powers and at least 11s+1 cycles in the Lu comparison. "
          "All three term-count residues give coefficient 11/3.")
    print("PASS One radius in degree 12: the odd moment gives 85/77. The tail signs are - then +; "
          "ordinary brackets give 10 on eight powers and 11 on nine before exterior closure.")
    print("PASS The original degree-14 example allows powers 1 through 14 except 4 and 6. "
          "11 ordinary cycles, 10 incoming changes and 2 exterior additions give 23.")
    print("PASS The original degree-15 row allows 13 powers, all except 4 and 6. "
          "12 ordinary cycles, 10 incoming changes and 2 exterior additions give 24.")
    return ordinary_counts


def check_matched_turn_variation():
    """Eliminate the turning-level shift in the two first-order action expansions."""
    shift, rising_slope = symbols("shift rising_slope")
    falling_slope = symbols("falling_slope", nonzero=True)
    rising_variation, falling_variation = symbols("rising_variation falling_variation")
    # shift is (minus level - plus level)/perturbation. Both outgoing actions
    # match to first order, with fixed-level variations +V and -V.
    paired_shift, = solve(falling_variation - (falling_slope*shift - falling_variation), shift)
    assert paired_shift == 2*falling_variation/falling_slope
    inward_minus = rising_slope*paired_shift - rising_variation
    lu_log_gap_ratio = inward_minus - rising_variation
    derivative_log_ratio = rising_variation - inward_minus
    corrected = rising_variation - rising_slope*falling_variation/falling_slope
    assert expand(lu_log_gap_ratio + 2*corrected) == 0
    assert expand(derivative_log_ratio - 2*corrected) == 0
    print("PASS Equal upper hits: the shift gives the same descending correction in the "
          "lower-gap and derivative log ratios, with respective prefactors -2 and +2.")


def check_short_profiles():
    """Compare the six printed boxes with the outward-rounded transcript endpoints."""
    section = (SOURCE / "tex/valley-threshold.tex").read_text()
    table = section.split(r"p & J(\eta^2)/10^9 & J(\eta)/10^5 & J(1)", 1)[1]
    table = table.split(r"\end{array}", 1)[0]
    for stride, filename, field in [(2, "stride-two-fixed.txt", " endpoints: "),
                                    (4, "mechanism-replay.txt", ": J lower ")]:
        row = next(line for line in table.splitlines() if line.strip().startswith(f"{stride} &"))
        printed = re.findall(r"\(([-.0-9]+),([-.0-9]+)\)", row)
        log = (SOURCE / "verification" / filename).read_text()
        records = [line.split(field, 1)[1] for line in log.splitlines() if field in line]
        for (low, high), scale, record in zip(printed, (10**9, 10**5, 1), records, strict=True):
            lower, upper = re.findall(r"\[([-+.e0-9]+) \+/- ([-+.e0-9]+)\]", record)
            assert Fraction(low)*scale < Fraction(lower[0]) - Fraction(lower[1])
            assert Fraction(upper[0]) + Fraction(upper[1]) < Fraction(high)*scale
    print("PASS All six printed short-profile boxes contain the corresponding Arb endpoints.")


def check_depth_scales():
    # Squared exponential base of q H / (H+m)^(3/2) for H=w^n, m~(2/3)^n.
    k0, q0 = Fraction(2, 3), Fraction(3, 4)
    points = [Fraction(1, 2), k0, Fraction(7, 8), Fraction(7, 10)]
    bases = [q0*w**2/max(w, k0)**3 for w in points]
    assert bases == [Fraction(81, 128), Fraction(9, 8), Fraction(6, 7), Fraction(15, 14)]
    assert points[-1]/k0 == Fraction(21, 20)
    print("PASS Squared valley bases at w=1/2, 2/3, 7/8 and 7/10:", *bases)


def check_papers(lu_ordinary_counts):
    ledger = json.loads((SOURCE / "inputs/supports.json").read_text())
    main = (SOURCE / "main.tex").read_text()
    supplement = (SOURCE / "low-term-supplement.tex").read_text()
    ten = ledger["ten_power_input"]

    seed = json.loads((SOURCE / "inputs/degree-ten-seed.json").read_text())
    seed_proof = (SOURCE / "tex/balanced-seed.tex").read_text()
    for key in ("radius", "peak"):
        entry = seed[key]
        assert f'I_{entry["digits"]}({entry["left_integer"]})' in seed_proof.replace("I_{15}", "I_15")
    assert "$I_9(N_" in seed_proof and seed["turn_digits"] == 9
    for h, plus, minus, lower, upper in seed["tests"]:
        assert f"{h}&{plus}&{minus}&({lower},{upper})" in seed_proof
    print("PASS Printed degree-ten algebraic radius, peak, six outer radius cells and three action enclosures.")

    proof = (SOURCE / "tex/conclusion.tex").read_text()
    support_block = proof.split(r"\label{eq:ten-support}", 1)[1].split(r"\end{equation}", 1)[0]
    high = [int(s) for s in re.findall(r"\d+", support_block)]
    assert high == ten["exponents"] and len(set(high)) == 10 and max(high) == 67088
    assert ten["n"] == 257 and ten["precision_bits"] == 768

    table = main.split(r"\label{tab:ten-comparison}", 1)[1].split(r"\bottomrule", 1)[0]
    data_rows = table.split(r"\midrule", 1)[1].strip().splitlines()
    entries = []
    degree_column = []
    for row in data_rows:
        _, degree, cycles, unrestricted = row.removesuffix(r"\\").split("&")
        entries.append((degree.strip(), int(cycles)))
        degree_column.append(int(unrestricted))
    nine = next(row for row in ledger["low"] if row["k"] == 9)
    assert nine["project_support"][-1] == "M" and nine["project_hyperbolic"] == 19
    eight_fixed = nine["project_support"][:-1]
    assert sorted(set(eight_fixed)) == [5, 9, 10, 1033, 2058, 3081, 3082, 3083]
    nine_section = (SOURCE / "tex/two-matches.tex").read_text().split(r"\label{subsec:nine-powers}", 1)[1]
    assert "every sufficiently large even" in nine_section.lower() and "$M>3083$" in nine_section
    comparisons = main.split(r"\label{subsec:other-comparisons}", 1)[1].split(
        r"\input{tex/half-returns.tex}", 1)[0]
    assert "adjoins $x^{M+1}$" in nine_section and "at least twenty" in nine_section
    assert nine["project_hyperbolic"] + 1 == 20
    odd_lifts = ledger["public_odd_lift"]
    lifted = odd_lifts["seed"] + odd_lifts["new_powers"]
    assert len(set(lifted)) == 10 and max(lifted) == 23
    palomar = json.loads((SOURCE / "inputs/palomar-samples.json").read_text())
    sextic = list(range(1, len(palomar["coefficients"]) + 1))
    printed = re.search(r"support\s+\$\\\{(1,2,3,4,5,6,[\d,]+)\\\}\$\s+of\s+degree", supplement)
    sextic_lift = list(map(int, printed[1].split(",")))
    appended = sextic_lift[len(sextic):]
    assert sextic_lift[:len(sextic)] == sextic and appended == [7, 9, 11, 13]
    assert len(set(sextic_lift)) == 10 and "$7,9,11,13$" in comparisons
    sextic_count = len(palomar["samples"]) - 1 + len(appended)
    # Published full-support bounds: LMP (deg 3,5), LL (deg 4), DMH (deg >= 6).
    published_full = {3: 1, 4: 1, 5: 2, 6: 4, 7: 5, 8: 6, 9: 7, 10: 8}
    # Preprint exact-degree open sets in Lu, Supplement Proposition S1.2, for k >= 6.
    lu_low = {6: 6, 7: 7, 8: 9, 9: 10, 10: 14}
    assert entries == [("10", published_full[10]), ("10", lu_low[10]), (str(max(sextic_lift)), sextic_count),
                       ("19", 10-1), ("22", 16), (str(max(lifted)), len(lifted)), ("26", 16),
                       ("67088", ten["cycles"])]
    low_degree_bounds = {10: lu_low[10], 13: 22}
    for (degree, _), printed_count in zip(entries, degree_column, strict=True):
        n = int(degree)
        expected = low_degree_bounds[n] if n < 14 else (7*n - 48 + 1)//2
        assert printed_count == expected, (n, printed_count, expected)
    degree_counts = [(7*n - 48 + 1)//2 for n in (22, 26)]
    assert degree_counts == [53, 67]
    assert "gives $%d$ and $%d$" % tuple(degree_counts) in comparisons
    assert "nineteen-cycle nest" in nine_section and "at least twenty" in nine_section
    corollary = (SOURCE / "tex/valley-threshold.tex").read_text().split(
        r"\label{cor:short-lists}", 1)[1].split(r"\end{corollary}", 1)[0]
    short_supports = re.findall(r"\\\{(\d+(?:,\d+)+)\\\}", corollary)
    assert len(short_supports) == 2
    for stride, printed in zip((2, 4), short_supports, strict=True):
        top = 10 + 2 * stride
        powers = sorted([5, 9, 10, 9+stride, 10+stride, 9+2*stride, top, top+1, top+3, top+8])
        assert len(set(powers)) == 10 and list(map(int, printed.split(","))) == powers
        assert (str(max(powers)), 16) in entries
    print("PASS Table 1 and the three valley supports: 26 in degree 67088, 16 each in degrees 22 and 26.")
    print("PASS Article comparison rows: averaging (19,9), published full support (10,8), odd extensions (13,10) and (23,10), "
          "preprint full support (10,14). Pairs are (degree, count).")
    print("PASS Supplement 6.4: the nineteen-cycle nest and one odd extension give twenty on ten terms.")
    print("PASS Exact-degree comparison at 22 and 26: Lu's rounded bounds are 53 and 67.")
    print("PASS Article Table 1, unrestricted exact-degree column:", *degree_column)

    initial, blocks, prefix = symbols("n0 L j", integer=True)
    degree = 258 * (initial + 3 * blocks) + 8
    terms = 4 * blocks + 6 + prefix
    offset = 258 * initial + 8 - Rational(387, 2) * (6 + prefix)
    assert expand(degree - Rational(387, 2) * terms - offset) == 0
    assert degree.subs({initial: ten["n"], blocks: 1}) == max(high)
    print("PASS Degree identity from (4.1) in all residue classes. At n0=257, L=1 the degree is 67088.")

    for snippet in (r"z_p &\in (0.99229,\,0.99230)", r"1-z_c &\in (2.39,\,2.40)\cdot 10^{-44}",
                    r"m &\in (5.55,\,5.56)\cdot 10^{-46}", r"M &\in (2.42,\,2.43)\cdot 10^{-5}",
                    r"z_a\in(0.67022,\,0.67023)", r"7/8<z_p",
                    "(1.76,1.77)", "(-5.11,-5.10)", "(2.07,2.08)", "(-3.71,-3.70)",
                    "(5.71,5.72)", "(-1.22,-1.21)", "(9.84,9.85)"):
        assert snippet in proof, snippet
    boxes = proof.split(r"\label{tab:main-depth-cert}", 1)[1].split(r"\bottomrule", 1)[0]
    rows = re.findall(r"\$([^$]+)\$\s*&\s*\$\(([^,]+),([^)]+)\)\$\s*&\s*\$([+-])\$", boxes)
    for printed, node, bounds in zip(rows, ten["negative_w"], ten["negative_scaled_boxes"], strict=True):
        w, low, high, sign = printed
        assert Fraction(w) == Fraction(node)
        actual = [decimal_or_power(low.strip()), decimal_or_power(high.strip())]
        assert actual == list(map(Fraction, bounds))
        assert actual[0] < actual[1]
        assert (sign == "+" and actual[0] > 0) or (sign == "-" and actual[1] < 0)
    print("PASS Printed critical-point/positive-variation boxes and all five depth enclosures agree with the archived inputs.")

    ordered = (SOURCE / "tex/ordered-returns.tex").read_text().split(r"\midrule")
    added = [int(s) for s in re.findall(r"&\s*(\d+)\\\\", ordered[1])]
    total = int(re.search(r"Total\s*&\s*(\d+)", ordered[2])[1])
    marks = ordered[1].count("$0$ (simple)")
    assert added == [2, 6, 1, 1, 1, 6, 2, 4, 2, 1]
    assert sum(added) == total == ten["cycles"] and marks == 2
    print("PASS Ordered-return table: 24 strict brackets plus two marked orbits.")

    small = supplement.split(r"\label{thm:lower}", 1)[1].split(r"\end{array}", 1)[0]
    rows = re.findall(r"^\s*(\d+)\s*&\s*(\d+)\s*&\s*\\\{(.+?)\\\}", small, flags=re.MULTILINE)
    for (k, cycles, printed), item in zip(rows, ledger["low"], strict=True):
        terms = printed.replace("2^{38}", str(2**38)).split(",")
        assert int(k) == item["k"] and int(cycles) == item["project_hyperbolic"]
        assert terms == [str(power) for power in item["project_support"]]
    comparison = supplement.split(r"\label{tab:low-comparison}", 1)[1].split(r"\bottomrule", 1)[0]
    number_rows = [line.split("&") for line in comparison.splitlines()
                   if re.match(r"^\s*\d+\s*&", line)]
    for cells, item in zip(number_rows, ledger["low"], strict=True):
        k, odd_lift, filled, samples, lu, preclosure, here = [cell.strip() for cell in cells[:7]]
        k = int(k)
        assert (k, int(odd_lift), int(here)) == (
            item["k"], item["seed_or_odd_extension"], item["project_hyperbolic"])
        assert filled == str(published_full[k])
        assert samples == (str((k - 1)//2 + 4) if k >= 6 else "--")
        assert lu == str(lu_low.get(k, "--"))
        assert preclosure == str(lu_ordinary_counts.get(k, "--"))
        top = item["project_support"][-1]
        special_degrees = {"d": r"even $d\ge10$", 2**38: "$2^{38}$", "M": "large even $M$"}
        degree_text = special_degrees.get(top, str(top))
        assert cells[7].strip().removesuffix(r"\\") == degree_text
    support8 = re.search(r"S_8=\\\{([0-9,]+)\\\}", supplement)
    assert list(map(int, support8[1].split(","))) == [1, 2, 3, 5, 7, 8, 10, 12]
    assert r"S_9=S_8\cup\{9\}" in supplement
    assert "$6,7,7,8,8$ brackets" in supplement and (10 - 1)//2 + 4 == 10 - 2
    print("PASS Full supports at degrees 3-9: published 1, 1, 2, 4, 5, 6, 7; "
          "Palomar sample chains at 6-9 give 6, 7, 7, 8. Both give 8 at degree 10.")
    print("PASS Preprint exact-degree open-set comparison on full supports of degrees 6 through 10: "
          "6, 7, 9, 10, 14.")
    print("PASS Supplement Theorem 1.1 and Table 1: all seven supports, cycle counts and degree restrictions agree with the ledger.")


if __name__ == "__main__":
    if not __debug__:
        raise SystemExit("Run with assertions enabled.")
    check_rational_margins()
    check_rising_witness()
    check_depth_scales()
    lu_ordinary_counts = check_lu_term_comparison()
    check_matched_turn_variation()
    check_papers(lu_ordinary_counts)
    check_short_profiles()
