"""Check the eight published depth intervals on finer meshes, independently of the two drivers."""
from fractions import Fraction
from flint import arb, ctx


def point(x):
    return arb(str(x))


def square_at_ends(x, y):
    # The panels split at the zero, so the endpoint squares enclose the range.
    return (x*x).union(y*y)


def within(name, value, lower, upper):
    assert point(lower) < value < point(upper), name
    print(f"{name}: lower {value.lower().str(10)}, upper {value.upper().str(10)}"
          f" in ({lower}, {upper})")


def high_degree():
    ctx.prec = 768
    n = 257
    eta, q = point(Fraction(2, 3)**n), point(Fraction(3, 4)**128)
    # A coarse cell meeting the minimum illustrates cancellation in the expanded cubic.
    cell = point(1-Fraction(1,4096)).union(point(1+Fraction(1,4096)))
    naive = eta + cell**n * (1-(3-eta)*cell**2+2*cell**3)
    assert naive.contains(0)
    print("Expanded denominator on [1-2^-12,1+2^-12] at H=eta contains zero.")

    # Order all dyadic breakpoints from left to right; double the published counts.
    breaks = {Fraction(x) for x in ("1/8","1/4","1/2","1","3/2","2","4","8")}
    for k in range(1,101):
        radius = Fraction(1,2**k)
        breaks.update((1-radius,1+radius))
    ordered = sorted(breaks)
    ws = [Fraction(x) for x in ("1/2","2/3","7/8","3/2","4")]
    depths = [point(w**n) for w in ws]
    totals = [arb(0) for _ in depths]
    cells = 0
    for left, right in zip(ordered, ordered[1:]):
        if left == 1 or right == 1:
            count = 2
        elif right <= Fraction(1,2) or left >= Fraction(3,2):
            count = 2048
        else:
            count = 1024
        cells += count
        width = point((right-left)/count)
        start = point(left)
        for _ in range(count):
            end = start + width
            z = start.union(end)
            sq = square_at_ends(start-1,end-1)
            energy = z**n * (sq*(1+2*z)+eta*z**2)
            odd = (1-z)*(8-z-z**2)/8 - q*z
            density = n * width * z**(n-1) * z.root(258) * odd
            for i,h in enumerate(depths):
                assert h+energy > 0
                totals[i] += density*h/(h+energy)**2
            start = end
    assert cells == 212996
    print(f"Degree 67088, {cells} cells at {ctx.prec} bits; normalized n H D J(H).")
    boxes = [("7/10","6/5"),("-400000000","-300000000"),("12/5","3"),
             ("-29/100","-23/100"),("1/25","7/100")]
    for w,h,total,(lo,hi) in zip(ws,depths,totals,boxes,strict=True):
        tail = (2/h+h) * point(Fraction(1,8**n))
        within(f"w={w}",total+(-tail).union(tail),lo,hi)


def low_degree():
    ctx.prec = 512
    eta, q = arb(2)**-20, arb(2)**-10
    depths = [eta**2,eta,arb(1)]
    cells = 2**17
    width = arb(2)/cells
    totals = [arb(0) for _ in depths]
    for i in range(cells):
        left, right = i*width,(i+1)*width
        r = left.union(right)
        sq = square_at_ends(left**4-1,right**4-1)
        energy = r**10*sq+eta*r**14
        numerator = width*(r**10*sq-q*r**14)
        for j,h in enumerate(depths):
            assert h+energy > 0
            totals[j] += numerator/(h+energy)**2
    tail = point(Fraction(48,17*2**17))
    print(f"Degree 26, {cells} cells at {ctx.prec} bits; unnormalized J(H).")
    for h,total,(lo,hi) in zip(("eta^2","eta","1"),totals,
                                (("6900000000","7100000000"),("-140000","-100000"),("0.044","0.045")),
                                strict=True):
        within(f"H={h}",total+(-tail).union(tail),lo,hi)


if __name__ == "__main__":
    if not __debug__:
        raise SystemExit("Enable assertions.")
    high_degree()
    low_degree()
    print("PASS eight refined whole-cell sums including the two analytic tail estimates.")
