"""Numerical and symbolic controls for the leading angular Ward kernel."""
import mpmath as mp
import sympy as sp

from pathlib import Path


def main():
    t, beta = sp.symbols("t beta", real=True)
    angular = sp.diff((1 - t**2) * sp.diff(sp.log(1 - beta*t), t), t)
    assert sp.simplify(angular + 1 - (1 - beta**2)/(1 - beta*t)**2) == 0

    mp.mp.dps = 75
    lines = ["Exact symbolic sphere-divergence identity: passed."]
    for beta in map(mp.mpf, ("0.25", "0.6", "0.9")):
        r = beta/(1 + mp.sqrt(1 - beta**2))

        def w(t):
            return (1 - beta**2)/(2*(1 - beta*t)**2)

        def h(t):
            return beta*(1 - t**2)/(2*(1 - beta*t))

        assert abs(mp.quad(w, [-1, r, 1]) - 1) < mp.mpf("1e-70")
        for ell in range(13):
            original = mp.quad(lambda t: w(t)*mp.legendre(ell, t), [-1, r, 1])
            positive = (ell + 1)*mp.quad(lambda t: w(t)*h(t)**ell, [-1, r, 1])
            assert positive > 0 and abs(original/positive - 1) < mp.mpf("1e-50")
        # Divide out the exponential before quadrature at large degree.
        scaled = lambda ell: mp.sqrt(ell + 1)*mp.quad(
            lambda t: w(t)*(h(t)/r)**ell, [-1, r, 1]
        )
        v = [scaled(deg) for deg in (10, 40, 160)]
        limit = mp.sqrt(mp.pi*(1 - r*r)) / 2
        assert abs(v[-1]/limit - 1) < mp.mpf("0.01")
        lines.append(f"beta={beta}: two formulas agree for degrees 0 to 12; "
                     f"normalized degrees 10, 40, 160: {', '.join(mp.nstr(x, 9) for x in v)}; "
                     f"asymptotic constant: {mp.nstr(limit, 9)}")
    report = Path(__file__).resolve().parents[1] / "reports" / "massive-angular.txt"
    report.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
