"""Arithmetic and ordinary high-precision checks independent of the certificate code."""
from mpmath import mp


def sign_changes(word):
    return sum(a != b for a, b in zip(word, word[1:]))


def support_and_counts():
    # The plus signs on the two fold test boundaries surround the incoming list.
    folded = '+' + '+-+-+-+' + '-+' + '+-+-+' + '-' + '+'
    ordinary = sign_changes('-+-') + sign_changes('-+-+-+-')
    crest_brackets = 2
    marks = 2
    assert ordinary + sign_changes(folded) + crest_brackets == 24
    assert ordinary + sign_changes(folded) + crest_brackets + marks == 26
    print(f"Degree 67088: {ordinary} ordinary, {sign_changes(folded)} incoming, "
          f"{crest_brackets} fold brackets, {marks} marks.")
    d = 258 * 257
    powers = [d//2, d-1, d, d+257, d+516, d+773, d+774, d+775, d+777, d+782]
    print("L=1, j=0, n=257:", *powers)
    assert powers == [33153, 66305, 66306, 66563, 66822, 67079, 67080, 67081, 67083, 67088]
    short_incoming = '+' + '+-+' + '-+' + '+-+' + '-' + '+'
    assert sign_changes(short_incoming) == 8
    for p in (4, 2):
        top = 10 + 2*p
        short = sorted([5, 9, 10, 9+p, 10+p, 9+2*p, top, top+1, top+3, top+8])
        assert len(set(short)) == 10
        print(f"Stride {p}: {short}, 4 ordinary + 8 incoming + 4 fold = 16.")


def fold_determinant(d):
    # r=exp(u) avoids subtracting the close endpoints 1 and g at large d.
    u_max = mp.log1p(mp.mpf(8)/d) / 8
    def integrand(u):
        f = -d * mp.exp((d-2)*u) * mp.expm1(8*u)
        denominator = 2*mp.cosh(3*u) + 2*mp.cosh(u)
        return f**2 * (6*mp.sinh(3*u) + 2*mp.sinh(u)) / (d * denominator**2)
    return mp.quad(integrand, [0, u_max])


def depth_integral(h, eta, q):
    def integrand(r):
        z = r**4
        valley = (1-z)**2
        return r**10 * (valley-q*z) / (h+r**10*(valley+eta*z))**2
    width = mp.sqrt(h+eta)/4
    near_origin = h**(mp.mpf(1)/10)
    nodes = sorted(set([0, near_origin/2, near_origin, mp.mpf('0.75'),
                        1-width, 1, 1+width, 2, mp.inf]))
    return mp.quad(integrand, nodes)


def numerical_checks():
    mp.dps = 65
    for d in (10, 14, 18, 67080):
        determinant = fold_determinant(d)
        assert determinant > 0
        print(f"d={d}, det={mp.nstr(determinant, 16)}")
    s = mp.mpf(9)/10
    alpha = mp.mpf(3)/2 - s
    t = s/alpha
    constant = 8 / (100*mp.sin(mp.pi/10)) * (1+t)**mp.mpf('1.5')/t**s
    print(f"d=10, p=4: alpha={alpha}, t={t}, K={mp.nstr(constant,16)}")
    eta = mp.mpf(2)**-20
    for h, reference in zip([eta**2, eta, 1],
                            [(mp.mpf('6.9e9'),mp.mpf('7.1e9')),
                             (-140000,-100000),(mp.mpf('.044'),mp.mpf('.045'))]):
        value = depth_integral(h, eta, mp.sqrt(eta))
        assert reference[0] < value < reference[1]
        print("J =", mp.nstr(value, 16))


if __name__ == "__main__":
    if not __debug__:
        raise SystemExit("Enable assertions.")
    support_and_counts()
    numerical_checks()
    print("PASS arithmetic and high-precision comparisons (ordinary quadrature).")
