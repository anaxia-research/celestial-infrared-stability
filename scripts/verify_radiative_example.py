#!/usr/bin/env python3
"""Exact and numerical replays for the classical Bondi-radiative example.

This verifies a specified radiative datum and conventional soft functional.
It does not construct a complete Einstein solution or quantum S-matrix.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import platform

import mpmath as mp
import sympy as sp


mp.mp.dps = 40


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def show(value: mp.mpf | mp.mpc) -> str:
    return mp.nstr(value, 24)


def exact_checks() -> dict:
    u = sp.symbols("u", real=True)
    tau, U, B, G = sp.symbols("tau U B G", positive=True)
    pulse = tau / (sp.pi * (u**2 + tau**2))
    shear_fraction = sp.Rational(1, 2) + sp.atan(u / tau) / sp.pi
    require(sp.simplify(sp.diff(shear_fraction, u) - pulse) == 0,
            "news must be the derivative of the shear")
    require(sp.limit(shear_fraction, u, sp.oo) == 1, "wrong final shear")
    require(sp.limit(shear_fraction, u, -sp.oo) == 0, "wrong initial shear")
    full_mass = sp.integrate(pulse, (u, -sp.oo, sp.oo))
    full_flux = sp.integrate(pulse**2, (u, -sp.oo, sp.oo))
    require(full_mass == 1, "Lorentzian memory normalization failed")
    require(sp.simplify(full_flux - 1 / (2 * sp.pi * tau)) == 0,
            "Lorentzian squared-news integral failed")
    window_fraction = sp.simplify(shear_fraction.subs(u, U) - shear_fraction.subs(u, -U))
    require(sp.simplify(window_fraction - 2 * sp.atan(U / tau) / sp.pi) == 0,
            "finite-window memory fraction failed")
    flux_primitive = (u / (u**2 + tau**2) + sp.atan(u / tau) / tau) / (2 * sp.pi**2)
    require(sp.simplify(sp.diff(flux_primitive, u) - pulse**2) == 0,
            "finite-window squared-news primitive failed")
    window_flux = sp.simplify(flux_primitive.subs(u, U) - flux_primitive.subs(u, -U))
    expected_flux = (U / (U**2 + tau**2) + sp.atan(U / tau) / tau) / sp.pi**2
    require(sp.simplify(window_flux - expected_flux) == 0,
            "finite-window news norm failed")

    # Axisymmetric normalized scalar harmonics supply an independent direct
    # coordinate check of the Bochner-derived tensor norm. For u=cos(theta),
    # E_theta_theta=(1-u^2)Y''(u)/2, E_phi_phi=-sin(theta)^2 E_theta_theta.
    # The normalization cancels pi and square roots in the following integral.
    angular_rows = []
    for degree in range(21):
        legendre = sp.legendre(degree, u)
        direct = sp.simplify(sp.Rational(2 * degree + 1, 4)
                             * sp.integrate((1 - u**2)**2 * sp.diff(legendre, u, 2)**2,
                                            (u, -1, 1)))
        lam = degree * (degree + 1)
        expected = sp.Rational(lam * (lam - 2), 2)
        require(direct == expected, f"STF harmonic norm failed at l={degree}")
        require(expected <= sp.Rational(1, 2) * (1 + lam)**2,
                "H2 test-function control failed")
        angular_rows.append({"l": degree, "lambda": lam,
                             "direct_tensor_norm_squared": str(direct),
                             "Bochner_formula": str(expected)})
    require(angular_rows[2]["direct_tensor_norm_squared"] == "12", "quadrupole norm failed")
    ell = sp.symbols("ell", positive=True)
    h_ell = (ell - 1) * ell * (ell + 1) * (ell + 2) / 2
    require(sp.limit(h_ell / ell**4, ell, sp.oo) == sp.Rational(1, 2),
            "large-angular-degree norm growth failed")

    # l=2, unit tensor T=E[Y]/sqrt(12), and memory b=B*T.
    correct_charge = -B * sp.sqrt(12) / (16 * sp.pi * G)
    wrong_charge = -B * sp.sqrt(12) / (4 * sp.pi * G)
    require(sp.simplify(wrong_charge / correct_charge) == 4,
            "misnormalized soft-charge negative control failed")
    correct_energy = B**2 * full_flux / (32 * sp.pi * G)
    require(sp.simplify(correct_energy - B**2 / (64 * sp.pi**2 * G * tau)) == 0,
            "gravitational flux prefactor failed")
    require(sp.simplify(correct_energy - B**2 / (32 * sp.pi**2 * G * tau)) != 0,
            "factor-two energy negative control did not reject")
    eps = sp.symbols("epsilon", positive=True)
    integrated_aligned_noise = sp.sqrt(2 * U) * eps
    paired_noise = sp.sqrt(12) * integrated_aligned_noise / (16 * sp.pi * G)
    operator_bound = sp.sqrt(12) * sp.sqrt(2 * U) * eps / (16 * sp.pi * G)
    require(sp.simplify(paired_noise - operator_bound) == 0,
            "aligned detector noise should attain the charge transfer bound")
    return {"news_integral": str(full_mass), "news_squared_integral": str(full_flux),
            "window_memory_fraction": str(window_fraction), "window_squared_news_integral": str(window_flux),
            "sphere_harmonics": angular_rows, "quadrupole_charge": str(correct_charge),
            "angular_high_degree_control": {
                "charge_operator_norm": "sqrt((l-1)*l*(l+1)*(l+2)/2)/(16*pi*G)",
                "asymptotic": "l^2/(16*pi*G*sqrt(2))",
                "uniform_H2_test_bound": "norm(ell_F)<=norm(F)_H2/(16*pi*G*sqrt(2))",
                "uniform_bound_over_L2_normalized_tests": False},
            "gravitational_energy": str(sp.simplify(correct_energy)),
            "negative_controls": {"wrong_4pi_charge_has_ratio": "4",
                                  "missing_factor_two_in_energy_rejected": True},
            "aligned_detector_noise_attains_dual_bound": True}


def numerical_checks() -> dict:
    fourier_rows = []
    for tau, omega in [(mp.mpf(".7"), mp.mpf(".3")),
                       (mp.mpf(".7"), mp.mpf("1.2")),
                       (mp.mpf("2.5"), mp.mpf(".3")),
                       (mp.mpf("2.5"), mp.mpf("1.2"))]:
        integrand = lambda u: tau * mp.cos(omega * u) / (mp.pi * (u**2 + tau**2))
        measured = 2 * mp.quadosc(integrand, [0, mp.inf], omega=omega)
        expected = mp.exp(-tau * omega)
        require(abs(measured - expected) < mp.mpf("1e-30"), "Lorentzian Fourier check failed")
        fourier_rows.append({"tau": show(tau), "omega": show(omega),
                             "quadrature": show(measured), "formula": show(expected),
                             "absolute_error": show(abs(measured - expected))})

    B, U, G = mp.mpf(3) / 2, mp.mpf(2), mp.mpf(1)
    time_rows = []
    for tau in [mp.mpf(".7"), mp.mpf(2), mp.mpf(10)]:
        pulse = lambda u: tau / (mp.pi * (u**2 + tau**2))
        full = B**2 * mp.quad(lambda u: pulse(u)**2, [-mp.inf, 0, mp.inf])
        expected_full = B**2 / (2 * mp.pi * tau)
        finite = B**2 * mp.quad(lambda u: pulse(u)**2, [-U, 0, U])
        expected_finite = B**2 / mp.pi**2 * (U / (U**2 + tau**2) + mp.atan(U / tau) / tau)
        fraction = mp.quad(pulse, [-U, 0, U])
        expected_fraction = 2 * mp.atan(U / tau) / mp.pi
        require(abs(full - expected_full) < mp.mpf("1e-32"), "full time-flux integral failed")
        require(abs(finite - expected_finite) < mp.mpf("1e-32"), "window time-flux integral failed")
        require(abs(fraction - expected_fraction) < mp.mpf("1e-32"), "window memory integral failed")
        time_rows.append({"tau": show(tau), "full_squared_news_norm": show(full),
                          "window_squared_news_norm": show(finite), "window_memory_fraction": show(fraction)})

    mellin_rows = []
    for a, nu, tau in [(mp.mpf(1), mp.mpf(0), mp.mpf(".7")),
                       (mp.mpf(1), mp.mpf(".8"), mp.mpf(2)),
                       (mp.mpf("1.5"), mp.mpf("2.25"), mp.mpf(".7"))]:
        z = mp.mpc(a, nu)
        measured = B * mp.quad(lambda omega: omega ** (z - 1) * mp.exp(-tau * omega),
                               [0, 1, 4, mp.inf])
        expected = B * mp.gamma(z) * tau ** (-z)
        require(abs(measured - expected) < mp.mpf("1e-30"), "news Mellin quadrature failed")
        norm = B * mp.gamma(a) * tau ** (-a)
        require(abs(expected) <= norm + mp.mpf("1e-35"), "Mellin-line norm exceeded")
        mellin_rows.append({"a": show(a), "nu": show(nu), "tau": show(tau),
                            "quadrature_error": show(abs(measured - expected)),
                            "sample_norm": show(abs(expected)), "whole_line_supremum": show(norm)})

    hiding_rows = []
    fixed_charge = -B * mp.sqrt(12) / (16 * mp.pi * G)
    last_norm, last_energy, last_fraction = mp.inf, mp.inf, mp.inf
    for tau in [mp.mpf(1), mp.mpf(10), mp.mpf(100), mp.mpf(1000)]:
        line_norm = B / tau  # a=1.
        energy = B**2 / (64 * mp.pi**2 * G * tau)
        fraction = 2 * mp.atan(U / tau) / mp.pi
        M = tau * B
        require(line_norm < last_norm and energy < last_energy and fraction < last_fraction,
                "fixed-memory hiding sequence should decrease observed norms")
        require(M >= B, "frequency derivative bound should grow with pulse duration")
        hiding_rows.append({"tau": show(tau), "memory_norm": show(B),
                            "quadrupole_soft_charge": show(fixed_charge), "a_one_Mellin_line_norm": show(line_norm),
                            "gravitational_energy": show(energy), "finite_time_memory_fraction": show(fraction),
                            "required_derivative_bound_M": show(M)})
        last_norm, last_energy, last_fraction = line_norm, energy, fraction

    # The IR-cutoff one-particle contribution has E1(2*tau*epsilon)
    # behavior; its logarithmic divergence is analytical, not a finite-grid claim.
    tau = mp.mpf(2)
    fock_rows = []
    previous = mp.mpf(0)
    for epsilon in [mp.mpf(".1"), mp.mpf(".001"), mp.mpf(".000001")]:
        direct = mp.quad(lambda w: mp.exp(-2 * tau * w) / w, [epsilon, 1, mp.inf])
        expected = mp.e1(2 * tau * epsilon)
        require(abs(direct - expected) < mp.mpf("1e-30"), "cutoff Fock integral failed")
        require(direct > previous, "IR-cutoff norm should grow as cutoff decreases")
        fock_rows.append({"epsilon": show(epsilon), "cutoff_integral_divided_by_B_squared": show(direct)})
        previous = direct

    return {"working_decimal_digits": mp.mp.dps, "fourier_convention": "integral exp(i*omega*u) N(u) du",
            "fourier_quadratures": fourier_rows, "time_quadratures": time_rows,
            "Mellin_quadratures": mellin_rows, "fixed_memory_hiding_sequence": hiding_rows,
            "infrared_Fock_negative_control": {"formula": "integral_epsilon^infinity exp(-2*tau*w)/w dw = E1(2*tau*epsilon)",
                                                "analytical_limit": "+infinity as epsilon goes to zero", "cutoff_checks": fock_rows}}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path,
                        default=Path(__file__).resolve().parents[1] / "reports" / "radiative-example.json")
    args = parser.parse_args()
    report = {"schema": "celestial-radiative-example-v1",
              "generated_utc": datetime.now(timezone.utc).isoformat(),
              "python_version": platform.python_version(), "sympy_version": sp.__version__,
              "mpmath_version": mp.__version__, "formal_proof": False,
              "scope": "A smooth electric-parity Bondi-radiative datum, normalized soft functional, and idealized direct-news detector.",
              "exact": exact_checks(), "numerical": numerical_checks(), "passed": True,
              "limits": ["The analytic proofs establish general statements; finite computations replay identities and numerical examples.",
                         "No complete Einstein scattering solution, antipodal matching, or hard-sector Ward identity is constructed.",
                         "Nonzero memory requires care beyond the ordinary zero-memory Fock representation.",
                         "News Mellin data are not unstripped graviton-amplitude principal-series data.",
                         "The idealized detector noise norm is assumed, not derived from an instrument or quantum dressing prescription."]}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(f"PASS: Bondi normalization, Lorentzian Fourier/Mellin/time integrals, STF harmonic norms, detector and Fock controls; {args.output}")


if __name__ == "__main__":
    main()
