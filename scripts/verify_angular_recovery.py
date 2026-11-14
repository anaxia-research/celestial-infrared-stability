#!/usr/bin/env python3
"""Replay angular Gaussian-recovery identities with labelled numerical scope.

The Hilbert space example uses the two orthonormal real modes
phi_0=1/sqrt(4*pi), phi_1=sqrt(3/(4*pi))*cos(theta) on the sphere.
This verifies a mathematical smeared coefficient, not a physical Ward identity.
The manuscript supplies the infinite-dimensional analytical argument.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import platform

import mpmath as mp
import sympy as sp


mp.mp.dps = 45


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def show(value: mp.mpf | mp.mpc) -> str:
    return mp.nstr(value, 24)


def norm(vector: list[mp.mpc]) -> mp.mpf:
    return mp.sqrt(mp.fsum(abs(value) ** 2 for value in vector))


def charge(test: list[mp.mpc], vector: list[mp.mpc]) -> mp.mpc:
    return mp.fsum(mp.conj(g) * value for g, value in zip(test, vector))


def kernel(nu: mp.mpf, a: mp.mpf, sigma: mp.mpf, h: mp.mpf) -> mp.mpc:
    z = mp.mpc(a, nu)
    return h ** (-z) * mp.exp(sigma**2 * z**2 / 2)


def kernel_norm(a: mp.mpf, sigma: mp.mpf, h: mp.mpf) -> mp.mpf:
    return h ** (-a) * mp.exp(a**2 * sigma**2 / 2) / (mp.sqrt(2 * mp.pi) * sigma)


def scalar_frequency_average(lam: mp.mpf, a: mp.mpf, sigma: mp.mpf,
                             h: mp.mpf, band: mp.mpf = mp.mpf(40)) -> mp.mpf:
    def integrand(nu: mp.mpf) -> mp.mpf:
        z = mp.mpc(a, nu)
        return mp.re(lam ** (-z) * mp.gamma(z) * kernel(nu, a, sigma, h))

    cuts = [mp.mpf(0)] + [mp.mpf(v) for v in (1, 3, 6, 12, 20) if v < band] + [band]
    return mp.quad(integrand, cuts) / mp.pi


def scalar_energy_average(lam: mp.mpf, sigma: mp.mpf, h: mp.mpf) -> mp.mpf:
    def integrand(z: mp.mpf) -> mp.mpf:
        return mp.exp(-lam * h * mp.exp(sigma * z) - z**2 / 2) / mp.sqrt(2 * mp.pi)

    return mp.quad(integrand, [-18, -6, 0, 6, 18])


def exact_angular_checks() -> dict:
    u = sp.symbols("u", real=True)
    modes = [1 / sp.sqrt(4 * sp.pi), sp.sqrt(3 / (4 * sp.pi)) * u]
    gram = sp.Matrix([[sp.simplify(2 * sp.pi * sp.integrate(f * g, (u, -1, 1)))
                       for g in modes] for f in modes])
    require(gram == sp.eye(2), "sphere mode normalization failed")
    r = sp.Matrix([sp.Rational(3, 2), -2 * sp.I / 3])
    g = sp.Matrix([sp.Rational(1, 2), 3 * sp.I / 4])
    target = sp.simplify((sp.conjugate(g).T * r)[0])
    require(target == sp.Rational(1, 4), "smeared coefficient target failed")
    field = sum(r[j] * modes[j] for j in range(2))
    test = sum(g[j] * modes[j] for j in range(2))
    angular_integral = sp.simplify(2 * sp.pi * sp.integrate(sp.conjugate(test) * field,
                                                          (u, -1, 1)))
    require(angular_integral == target, "angular integral disagrees with mode pairing")
    derivative_bound = sp.sqrt(sum(abs(r[j])**2 * [1, 3][j]**2 for j in range(2)))
    require(derivative_bound == sp.Rational(5, 2), "uniform derivative bound failed")
    epsilon = sp.symbols("epsilon", positive=True)
    punctured_norm_squared = sp.integrate(2 * sp.pi / (1 - u),
                                          (u, -1, 1 - epsilon))
    require(sp.simplify(punctured_norm_squared - 2 * sp.pi * sp.log(2 / epsilon)) == 0,
            "punctured collinear norm identity failed")
    require(sp.limit(punctured_norm_squared, epsilon, 0, dir="+") == sp.oo,
            "collinear negative control must diverge in L2")
    return {"orthonormal_gram_matrix": [[str(v) for v in row] for row in gram.tolist()],
            "smeared_target": str(target), "independent_angular_integral": str(angular_integral),
            "uniform_energy_derivative_bound_M": str(derivative_bound),
            "test_function_L2_norm": str(sp.sqrt(13) / 4),
            "collinear_negative_control": {
                "profile": "(1-cos(theta))^(-1/2)",
                "punctured_domain": "cos(theta)<=1-epsilon",
                "squared_L2_norm": "2*pi*log(2/epsilon)",
                "limit_epsilon_to_zero": "+infinity",
                "interpretation": "An angular pole cannot be inserted into an L2 theorem without checking its hypotheses."}}


def numerical_checks() -> dict:
    a, alpha = mp.mpf(1), mp.mpf(1)
    sigma, eta, M = 1 / mp.sqrt(2), mp.mpf("0.0001"), mp.mpf(5) / 2
    h = mp.pi ** (-mp.mpf(1) / 4) * mp.sqrt(eta / M)
    r = [mp.mpc(mp.mpf(3) / 2), mp.mpc(0, -mp.mpf(2) / 3)]
    g = [mp.mpc(mp.mpf(1) / 2), mp.mpc(0, mp.mpf(3) / 4)]
    lambdas = [mp.mpf(1), mp.mpf(3)]
    v = [component / norm(g) for component in g]
    require(abs(norm(v) - 1) < mp.mpf("1e-40"), "adversarial noise direction not unit")
    target = charge(g, r)
    require(abs(target - mp.mpf(1) / 4) < mp.mpf("1e-40"), "wrong numeric angular target")

    # The vector example is f(omega)=sum r_j exp(-lambda_j*omega) phi_j.
    # Each exp remainder is <=lambda_j*omega, yielding the stated M=5/2.
    exact = []
    identity_rows = []
    for lam, coefficient in zip(lambdas, r):
        freq = scalar_frequency_average(lam, a, sigma, h)
        energy = scalar_energy_average(lam, sigma, h)
        require(abs(freq - energy) < mp.mpf("1e-34"), "frequency/lognormal identity failed")
        exact.append(coefficient * freq)
        identity_rows.append({"lambda": show(lam), "frequency_average": show(freq),
                              "energy_average": show(energy),
                              "difference": show(abs(freq - energy))})

    general_kernel_rows = []
    for aa, ss, hh in [(".75", ".4", ".15"), ("1.25", "1.1", ".03")]:
        aa, ss, hh = mp.mpf(aa), mp.mpf(ss), mp.mpf(hh)
        freq = scalar_frequency_average(mp.mpf(2), aa, ss, hh)
        energy = scalar_energy_average(mp.mpf(2), ss, hh)
        require(abs(freq - energy) < mp.mpf("1e-33"), "general-line Gaussian duality failed")
        l1_numeric = mp.quad(lambda nu: abs(kernel(nu, aa, ss, hh)), [0, 2, 5, 10, 20, 40]) / mp.pi
        l1_formula = kernel_norm(aa, ss, hh)
        require(abs(l1_numeric / l1_formula - 1) < mp.mpf("1e-33"),
                "general-line Gaussian kernel norm failed")
        general_kernel_rows.append({"a": show(aa), "sigma": show(ss), "h": show(hh),
                                    "dual_identity_error": show(abs(freq - energy)),
                                    "L1_norm_relative_error": show(abs(l1_numeric / l1_formula - 1))})

    l1 = kernel_norm(a, sigma, h)
    # Noise n(nu)=eta*conj(q(nu))/abs(q(nu))*v has norm eta for every nu.
    # Therefore its exact estimator displacement is eta*||q||_1/(2*pi)*v.
    for nu in [mp.mpf(vv) for vv in (0, 1, 2, 5, 9)]:
        q = kernel(nu, a, sigma, h)
        aligned = [eta * mp.conj(q) / abs(q) * component for component in v]
        require(abs(norm(aligned) - eta) < mp.mpf("1e-40"), "adversarial noise exceeds eta")
        require(norm([q * aligned[j] - eta * abs(q) * v[j] for j in range(2)])
                < mp.mpf("1e-36"), "phase-aligned noise identity failed")
    noise_vector = [eta * l1 * component for component in v]
    noisy = [exact[j] + noise_vector[j] for j in range(2)]
    error_vector = [noisy[j] - r[j] for j in range(2)]
    bias_bound = M * h * mp.exp(sigma**2 / 2)
    upper = bias_bound + eta * l1
    norm_error = norm(error_vector)
    charge_error = abs(charge(g, noisy) - target)
    require(norm_error <= upper, "Hilbert recovery error exceeds analytical bound")
    require(charge_error <= norm(g) * norm_error + mp.mpf("1e-40"),
            "smeared coefficient violates Cauchy-Schwarz")
    aligned_charge = abs(charge(g, noise_vector))
    require(abs(aligned_charge - norm(g) * eta * l1) < mp.mpf("1e-40"),
            "angularly aligned noise does not attain charge operator norm")

    # Independently minimize the already h-optimized noise/bias constant in sigma.
    def optimized_constant(s: mp.mpf) -> mp.mpf:
        b = mp.exp(alpha**2 * s**2 / 2)
        c = mp.exp(a**2 * s**2 / 2) / (mp.sqrt(2 * mp.pi) * s)
        h0 = (a * c / (alpha * b)) ** (1 / (a + alpha))
        return b * h0**alpha + c * h0 ** (-a)

    numeric_sigma = mp.findroot(lambda s: mp.diff(optimized_constant, s), (mp.mpf(".5"), 1))
    expected_constant = 2 * (mp.e / mp.pi) ** (mp.mpf(1) / 4)
    require(abs(numeric_sigma - sigma) < mp.mpf("1e-35"), "bandwidth optimizer failed")
    require(abs(optimized_constant(numeric_sigma) - expected_constant) < mp.mpf("1e-35"),
            "optimized Hölder constant failed")
    require(abs(upper / mp.sqrt(M * eta) - expected_constant) < mp.mpf("1e-35"),
            "optimized dimensional bound mismatch")
    sigma_rows = []
    for s in [mp.mpf(".25"), mp.mpf(".5"), sigma, mp.mpf(1), mp.mpf("1.5")]:
        constant = optimized_constant(s)
        require(constant >= expected_constant - mp.mpf("1e-35"),
                "nonoptimal bandwidth beats exact optimum")
        sigma_rows.append({"sigma": show(s), "optimized_in_h_constant": show(constant)})

    # For this diagonal angular example the line norm is attained at nu=0.
    B = mp.gamma(a) * norm([r[j] * lambdas[j] ** (-a) for j in range(2)])
    bands = []
    for V in [mp.mpf(3), mp.mpf(6), mp.mpf(10)]:
        inside = mp.erf(sigma * V / mp.sqrt(2))
        outside = mp.erfc(sigma * V / mp.sqrt(2))
        numerical_l1_inside = mp.quad(lambda nu: abs(kernel(nu, a, sigma, h)), [0, V]) / mp.pi
        require(abs(numerical_l1_inside - l1 * inside) < mp.mpf("1e-34"),
                "finite-band L1 kernel identity failed")
        band_exact = [r[j] * scalar_frequency_average(lambdas[j], a, sigma, h, V)
                      for j in range(2)]
        tail_error = norm([exact[j] - band_exact[j] for j in range(2)])
        tail_bound = B * l1 * outside
        require(tail_error <= tail_bound + mp.mpf("1e-35"), "finite-band tail bound failed")
        band_noisy = [band_exact[j] + eta * l1 * inside * v[j] for j in range(2)]
        band_error = norm([band_noisy[j] - r[j] for j in range(2)])
        band_bound = bias_bound + l1 * (eta * inside + B * outside)
        require(band_error <= band_bound, "finite-band recovery bound failed")
        bands.append({"V": show(V), "kernel_mass_fraction": show(inside),
                      "kernel_tail_fraction": show(outside),
                      "actual_exact_data_truncation_error": show(tail_error),
                      "truncation_upper_bound": show(tail_bound),
                      "noisy_vector_error": show(band_error),
                      "recovery_upper_bound": show(band_bound)})

    return {"working_decimal_digits": mp.mp.dps,
            "quadrature_control": "Frequency cut at 40; Gaussian expectation cut at |Z|=18. These are numerical replays, not formal integral certificates.",
            "a": show(a), "alpha": show(alpha), "M": show(M), "eta": show(eta),
            "h": show(h), "sigma": show(sigma), "dual_identity_checks": identity_rows,
            "general_line_and_bandwidth_checks": general_kernel_rows,
            "kernel_L1_over_2pi": show(l1), "exact_target_Qg": show(target),
            "noisy_estimate_Qg": show(charge(g, noisy)), "vector_error": show(norm_error),
            "vector_upper_bound": show(upper), "charge_error": show(charge_error),
            "charge_upper_bound": show(norm(g) * upper),
            "aligned_noise_charge_shift": show(aligned_charge),
            "aligned_noise_attains_operator_bound": True,
            "independent_sigma_optimizer": show(numeric_sigma),
            "optimized_constant": show(expected_constant), "sigma_comparison": sigma_rows,
            "whole_line_norm_B": show(B), "finite_band_checks": bands}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path,
                        default=Path(__file__).resolve().parents[1] / "reports" / "angular-recovery.json")
    args = parser.parse_args()
    report = {"schema": "celestial-angular-recovery-v1",
              "generated_utc": datetime.now(timezone.utc).isoformat(),
              "python_version": platform.python_version(), "mpmath_version": mp.__version__,
              "sympy_version": sp.__version__, "formal_proof": False,
              "scope": "Two angular modes, Gaussian direct-Mellin recovery, bounded adversarial noise, and finite-band upper bounds.",
              "exact": exact_angular_checks(), "numerical": numerical_checks(), "passed": True,
              "limits": ["Finite computations support the written analytical theorems.",
                         "The L2 angular assumptions exclude untreated collinear poles.",
                         "No gauge Ward identity, physical soft charge law, or scattering realization is established."]}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(f"PASS: angular Gram identities, Gaussian recovery, adversarial charge noise, bandwidth optimum, finite-band controls; {args.output}")


if __name__ == "__main__":
    main()
