#!/usr/bin/env python3
"""Exact identities and labelled numerical checks for soft-residue recovery.

These computations support the accompanying analytical arguments. They are
not a formal proof, and a finite frequency grid is not a proof on a whole line.
Requires mpmath and SymPy; neither is downloaded by this script.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
from fractions import Fraction
import json
from math import comb, factorial
from pathlib import Path
import platform
import sys

import mpmath as mp
import sympy as sp


if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")
mp.mp.dps = 60


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def chebyshev_coefficients(max_degree: int) -> list[list[Fraction]]:
    """Coefficients in ascending powers, using T[n+1] = 2*x*T[n]-T[n-1]."""
    polynomials = [[Fraction(1)], [Fraction(0), Fraction(1)]]
    for degree in range(1, max_degree):
        next_poly = [Fraction(0)] * (degree + 2)
        for power, coefficient in enumerate(polynomials[degree]):
            next_poly[power + 1] += 2 * coefficient
        for power, coefficient in enumerate(polynomials[degree - 1]):
            next_poly[power] -= coefficient
        polynomials.append(next_poly)
    return polynomials


def beta_integer_integral(a: int, n: int) -> Fraction:
    """Integrate x**(a-1)*(1-x)**(n-a-1) by polynomial expansion."""
    exponent = n - a - 1
    return sum(
        (Fraction((-1) ** k * comb(exponent, k), a + k)
         for k in range(exponent + 1)),
        Fraction(0),
    )


def exact_checks() -> dict:
    polynomials = chebyshev_coefficients(24)
    for degree in range(1, 25):
        require(polynomials[degree][-1] == 2 ** (degree - 1),
                f"wrong exact leading coefficient for T_{degree}")
        require(len(polynomials[degree]) == degree + 1,
                "wrong polynomial degree")

    beta_rows = []
    for a in range(1, 7):
        for n in range(a + 1, a + 13):
            integrated = beta_integer_integral(a, n)
            factorial_value = Fraction(factorial(a - 1) * factorial(n - a - 1),
                                       factorial(n - 1))
            require(integrated == factorial_value,
                    f"beta polynomial integral failed: a={a}, N={n}")
            beta_rows.append({"a": a, "N": n, "beta": str(integrated)})

    radius = Fraction(3, 2)
    line_norm_rows = []
    for n in range(2, 33):
        norm = radius * beta_integer_integral(1, n)
        require(norm == radius / (n - 1), "wrong a=1 vertical-line norm")
        require(norm != radius / n, "wrong-norm negative control did not fail")
        line_norm_rows.append({"N": n, "R": str(radius),
                               "exact_sup_norm_at_a_1": str(norm)})

    # For integer N: B(z,N-z) = pi/sin(pi*z) * product(k-z)/(N-1)!.
    # The pole at z=0 comes from pi/sin(pi*z), with residue one.
    z = sp.symbols("z")
    residues = []
    for n in (2, 3, 5, 9):
        beta_product = (sp.pi / sp.sin(sp.pi * z)
                        * sp.prod(k - z for k in range(1, n))
                        / sp.factorial(n - 1))
        residue = sp.limit(z * sp.Rational(3, 2) ** z * beta_product, z, 0)
        require(residue == 1, f"unit soft residue failed for N={n}")
        residues.append({"N": n, "residue": str(residue)})

    # Exact algebraic Chebyshev nodes, including radicals for degrees 4--6.
    left, width = sp.Rational(3, 7), sp.Rational(5, 3)
    interpolation = []
    x = sp.symbols("x")
    for degree in range(1, 7):
        nodes = [sp.simplify(left + width * (1 + sp.cos(sp.pi * j / degree)) / 2)
                 for j in range(degree + 1)]
        weights = [sp.simplify(1 / sp.prod(nodes[j] - nodes[k]
                                          for k in range(degree + 1) if k != j))
                   for j in range(degree + 1)]
        for power in range(degree + 1):
            moment = sp.simplify(sum(weights[j] * nodes[j] ** power
                                      for j in range(degree + 1)))
            require(moment == (1 if power == degree else 0),
                    f"exact interpolation moment failed: m={degree}, k={power}")
        expected = 2 ** (2 * degree - 1) / width ** degree
        signed_norm = sp.simplify(sum((-1) ** j * weights[j]
                                       for j in range(degree + 1)))
        require(sp.simplify(signed_norm - expected) == 0,
                f"exact noise norm failed for degree {degree}")
        for j, weight in enumerate(weights):
            require(sp.simplify((-1) ** j * weight).is_positive is True,
                    "weight sign was not proved by exact algebra")
        cheb = sp.chebyshevt(degree, x)
        for j in range(degree + 1):
            value = sp.simplify(cheb.subs(x, sp.cos(sp.pi * j / degree)))
            require(value == (-1) ** j, "exact alternation failed")
        interpolation.append({"degree": degree, "left": str(left),
                              "width": str(width), "sharp_noise_factor": str(expected),
                              "moment_identities": degree + 1})

    # A deliberately inferior design is a control on the claimed optimality.
    degree = 4
    equal_nodes = [left + width * j / degree for j in range(degree + 1)]
    equal_weights = [1 / sp.prod(equal_nodes[j] - equal_nodes[k]
                                 for k in range(degree + 1) if k != j)
                     for j in range(degree + 1)]
    equal_norm = sum(abs(w) for w in equal_weights)
    optimal_norm = 2 ** (2 * degree - 1) / width ** degree
    require(sp.simplify(equal_norm / optimal_norm) == sp.Rational(4, 3),
            "equispaced negative control failed")

    return {
        "chebyshev_exact_leading_coefficients_degrees": [1, 24],
        "beta_integer_polynomial_integrals": beta_rows,
        "beta_unit_residues": residues,
        "whole_vertical_line_integer_a_1_constants": line_norm_rows,
        "whole_line_argument": (
            "For positive A_N, the Mellin integral triangle inequality bounds every "
            "imaginary frequency by its value at zero; zero frequency attains the "
            "bound. Exact beta integrals check that value, not a finite frequency grid."
        ),
        "algebraic_interpolation": interpolation,
        "negative_controls": {
            "wrong_R_over_N_rejected_for_N": [2, 32],
            "equispaced_degree_4_noise_ratio": "4/3",
        },
    }


def numerical_checks() -> dict:
    radius = mp.mpf(3) / 2
    quadrature_rows = []
    for n, real_part, frequency in (
        (3, "1", "0"), (3, "1", "0.75"),
        (5, "1.5", "3"), (7, "1", "5"),
        (12, "2.25", "8"),
    ):
        z = mp.mpc(real_part, frequency)
        # x = omega/(R+omega) maps the original energy integral to (0,1).
        integrand = lambda x: x ** (z - 1) * (1 - x) ** (n - z - 1)
        integral = radius ** z * mp.quad(integrand, [0, mp.mpf(1) / 8,
                                                    mp.mpf(1) / 2,
                                                    mp.mpf(7) / 8, 1])
        gamma_formula = radius ** z * mp.gamma(z) * mp.gamma(n - z) / mp.gamma(n)
        error = abs(integral - gamma_formula)
        require(error < mp.mpf("1e-40"),
                f"60-digit quadrature check failed for N={n}, z={z}: {error}")
        quadrature_rows.append({"N": n, "real_part_above_pole": real_part,
                                "imaginary_frequency": frequency,
                                "absolute_error": mp.nstr(error, 12)})

    frequency_rows = []
    for n in (2, 3, 7, 16, 32):
        exact_norm = radius / (n - 1)
        largest_ratio = mp.mpf(0)
        for frequency in (0, mp.mpf(1) / 2, 2, 5, 10):
            z = mp.mpc(1, frequency)
            value = radius ** z * mp.gamma(z) * mp.gamma(n - z) / mp.gamma(n)
            ratio = abs(value) / exact_norm
            require(ratio <= 1 + mp.mpf("1e-55"), "sample exceeds whole-line norm")
            if frequency:
                require(ratio < 1, "nonzero-frequency strictness control failed")
            largest_ratio = max(largest_ratio, ratio)
        frequency_rows.append({"N": n, "sampled_frequencies": ["0", "0.5", "2", "5", "10"],
                               "largest_sample_to_exact_norm_ratio": mp.nstr(largest_ratio, 15)})

    left, width = mp.mpf(3) / 7, mp.mpf(5) / 3
    cheb_rows = []
    for degree in range(1, 13):
        nodes = [left + width * (1 + mp.cos(mp.pi * j / degree)) / 2
                 for j in range(degree + 1)]
        weights = [1 / mp.fprod(nodes[j] - nodes[k]
                                for k in range(degree + 1) if k != j)
                   for j in range(degree + 1)]
        expected = mp.mpf(2) ** (2 * degree - 1) / width ** degree
        norm_relative_error = abs(mp.fsum(abs(w) for w in weights) / expected - 1)
        moment_error = max(abs(mp.fsum(weights[j] * nodes[j] ** power
                                       for j in range(degree + 1))
                               - (1 if power == degree else 0))
                           for power in range(degree + 1))
        require(norm_relative_error < mp.mpf("1e-45"), "numerical Chebyshev norm failed")
        require(moment_error < mp.mpf("1e-43"), "numerical interpolation moment failed")
        cheb_rows.append({"degree": degree,
                          "noise_factor_relative_error": mp.nstr(norm_relative_error, 12),
                          "maximum_moment_absolute_error": mp.nstr(moment_error, 12)})
    return {"working_decimal_digits": mp.mp.dps,
            "quadrature_absolute_tolerance": "1e-40",
            "quadrature_spot_checks": quadrature_rows,
            "frequency_spot_checks_not_a_whole_line_proof": frequency_rows,
            "chebyshev_interpolation_numerical_degrees_1_to_12": cheb_rows}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path,
                        default=Path(__file__).resolve().parents[1] / "reports" / "recovery.json")
    args = parser.parse_args()
    report = {
        "schema": "celestial-recovery-checks-v1",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "python_version": platform.python_version(),
        "sympy_version": sp.__version__, "mpmath_version": mp.__version__,
        "formal_proof": False,
        "scope": "beta-family identities and fixed-degree Chebyshev recovery only",
        "exact": exact_checks(),
        "numerical": numerical_checks(),
        "passed": True,
        "limits": [
            "Finite computations support, but do not replace, the analytical proofs.",
            "No physical S-matrix, crossing, unitarity or celestial duality is established.",
            "The polynomial-log model is a fixed-degree assumption, not a resummation theorem.",
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(f"PASS: exact beta identities, algebraic Chebyshev m<=6, numerical m<=12; {args.output}")


if __name__ == "__main__":
    main()
