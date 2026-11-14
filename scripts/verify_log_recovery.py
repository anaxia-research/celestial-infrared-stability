#!/usr/bin/env python3
"""Replay exact identities and numerical controls for direct soft-log recovery.

The written theorem supplies the infinite-dimensional analytical proof. These
finite checks distinguish symbolic algebra from numerical quadrature. They do
not prove an S-matrix realization, a loop-resummation theorem, or priority.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from math import comb, factorial
from pathlib import Path
import platform
import sys

import mpmath as mp
import sympy as sp


if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")
mp.mp.dps = 50


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def gaussian_polynomial(power: int, center: sp.Expr,
                        sigma: sp.Expr) -> sp.Expr:
    """Exact E[(center - sigma Z)^power] for standard normal Z."""
    return sp.expand(sum(
        sp.binomial(power, 2 * j) * center ** (power - 2 * j)
        * sigma ** (2 * j) * sp.factorial2(2 * j - 1)
        for j in range(power // 2 + 1)
    ))


def exact_checks() -> dict:
    center, spacing, sigma = sp.symbols("L d sigma", real=True)
    q = sp.symbols("q")
    t = sp.symbols("t", nonnegative=True)
    z, alpha = sp.symbols("z alpha", positive=True)
    rows = []
    wrong_direction_rejected = []
    missing_factorial_rejected = []

    for degree in range(1, 7):
        coefficients = sp.symbols(f"a0:{degree + 1}")
        smoothed = sum(coefficients[k] * gaussian_polynomial(k, center, sigma)
                       for k in range(degree + 1))
        difference = sp.expand(sum(
            (-1) ** (degree - j) * comb(degree, j)
            * smoothed.subs(center, center + j * spacing)
            for j in range(degree + 1)
        ))
        expected = sp.factorial(degree) * coefficients[-1] * spacing ** degree
        require(sp.expand(difference - expected) == 0,
                f"Gaussian nuisance cancellation failed for degree {degree}")

        multiplier = sum((-1) ** (degree - j) * comb(degree, j) * q ** j
                         for j in range(degree + 1))
        require(sp.expand(multiplier - (q - 1) ** degree) == 0,
                "finite-difference frequency multiplier failed")

        reversed_difference = sp.expand(sum(
            (-1) ** (degree - j) * comb(degree, j)
            * smoothed.subs(center, center - j * spacing)
            for j in range(degree + 1)
        ))
        if degree % 2:
            require(sp.expand(reversed_difference - expected) != 0,
                    "wrong-direction negative control was not rejected")
            wrong_direction_rejected.append(degree)
        if degree >= 2:
            require(sp.expand(difference - coefficients[-1] * spacing ** degree) != 0,
                    "missing-factorial negative control was not rejected")
            missing_factorial_rejected.append(degree)

        laplace_moment = sp.integrate(t ** degree * sp.exp(-z * t),
                                     (t, 0, sp.oo))
        require(sp.simplify(laplace_moment
                            - sp.factorial(degree) / z ** (degree + 1)) == 0,
                "cutoff-profile Mellin identity failed")

        envelope = t ** degree * sp.exp(-alpha * t)
        derivative = sp.diff(envelope, t)
        expected_derivative = (t ** (degree - 1) * sp.exp(-alpha * t)
                               * (degree - alpha * t))
        require(sp.simplify(derivative - expected_derivative) == 0,
                "global remainder envelope derivative failed")
        maximum = sp.simplify(envelope.subs(t, degree / alpha))
        require(sp.simplify(maximum
                            - (degree / (alpha * sp.E)) ** degree) == 0,
                "global remainder maximum identity failed")
        require(sp.limit(envelope, t, sp.oo) == 0,
                "global remainder envelope endpoint failed")

        rows.append({
            "degree": degree,
            "gaussian_variance": "sigma^2, symbolic",
            "cancelled_nuisance_degrees": list(range(degree)),
            "difference_divided_by_d_power": f"{factorial(degree)} * a_{degree}",
            "frequency_multiplier": "(exp(d*z)-1)^m",
            "cutoff_Mellin_formula": "c * h0^z * m! / z^(m+1)",
            "global_remainder_supremum": "(m/(alpha*e))^m",
        })

    return {
        "degrees": rows,
        "full_principal_part": exact_jet_checks(),
        "negative_controls": {
            "wrong_scale_shift_direction_rejected_degrees": wrong_direction_rejected,
            "missing_factorial_rejected_degrees": missing_factorial_rejected,
        },
    }


def exact_jet_checks() -> dict:
    """Check inverse-heat weights and amplification bounds using exact algebra."""
    center, spacing, sigma, t = sp.symbols("L d sigma t", real=True)
    rows = []
    wrong_heat_sign = []
    rational_bound_checks = 0
    qvalues = (sp.Rational(2, 3), sp.Rational(3, 2))
    samples = ((sp.Rational(0), sp.Rational(5, 4), sp.Rational(3, 5)),
               (sp.Rational(3, 7), sp.Rational(5, 4), sp.Rational(3, 5)),
               (sp.Rational(4), sp.Rational(5, 4), sp.Rational(3, 5)))

    for degree in range(1, 7):
        coefficients = sp.symbols(f"u0:{degree + 1}")
        polynomial = sum(coefficients[k] * t ** k for k in range(degree + 1))
        smoothed = sp.expand(sum(
            (sigma ** 2 / 2) ** r / sp.factorial(r) * sp.diff(polynomial, t, 2 * r)
            for r in range(degree // 2 + 1)
        ))
        lagrange = [sp.Poly(sp.prod(
            (t - center - i * spacing) / (spacing * (j - i))
            for i in range(degree + 1) if i != j), t)
            for j in range(degree + 1)]

        for k in range(degree + 1):
            n = degree - k
            weights = [sp.expand(sum(
                (-sigma ** 2 / 2) ** r / sp.factorial(r)
                * sp.factorial(k + 2 * r) * ell.coeff_monomial(t ** (k + 2 * r))
                for r in range(n // 2 + 1))) for ell in lagrange]
            reconstructed = sp.expand(sum(
                weights[j] * smoothed.subs(t, center + j * spacing)
                for j in range(degree + 1)))
            require(sp.expand(reconstructed - sp.factorial(k) * coefficients[k]) == 0,
                    f"full principal-part reconstruction failed: m={degree}, k={k}")

            for j, weight in enumerate(weights):
                require(sp.degree(weight, center) <= n,
                        "principal-part weight has an excessive L degree")
                leading = sp.Poly(weight, center).coeff_monomial(center ** n)
                expected = ((-1) ** (n + degree - j) * sp.binomial(degree, j)
                            / (sp.factorial(n) * spacing ** degree))
                require(sp.simplify(leading - expected) == 0,
                        "principal-part leading L coefficient failed")

            wrong_weights = [sp.expand(sum(
                (sigma ** 2 / 2) ** r / sp.factorial(r)
                * sp.factorial(k + 2 * r) * ell.coeff_monomial(t ** (k + 2 * r))
                for r in range(n // 2 + 1))) for ell in lagrange]
            wrong_reconstruction = sp.expand(sum(
                wrong_weights[j] * smoothed.subs(t, center + j * spacing)
                for j in range(degree + 1)))
            if n >= 2:
                require(sp.expand(wrong_reconstruction
                                  - sp.factorial(k) * coefficients[k]) != 0,
                        "wrong inverse-heat sign was not rejected")
                wrong_heat_sign.append({"m": degree, "k": k})

            # Each comparison here is rational arithmetic. The geometric
            # factors q^j represent e^(+/- lambda*j*d), with q>0. This checks
            # the weighted H-bound without numerical exponent evaluation.
            for sample_center, sample_spacing, sample_sigma in samples:
                substitution = {center: sample_center, spacing: sample_spacing,
                                sigma: sample_sigma}
                evaluated_weights = [sp.factor(w.subs(substitution)) for w in weights]
                x = abs(sample_center) + degree * sample_spacing
                hpoly = sum((sample_sigma ** 2 / 2) ** r * x ** (n - 2 * r)
                            / (sp.factorial(r) * sp.factorial(n - 2 * r))
                            for r in range(n // 2 + 1))
                for j, weight in enumerate(evaluated_weights):
                    point_bound = (sp.factorial(degree) * hpoly
                                   / (sample_spacing ** degree
                                      * sp.factorial(j) * sp.factorial(degree - j)))
                    require(sp.simplify(point_bound - abs(weight)).is_nonnegative is True,
                            "exact rational pointwise H amplification bound failed")
                    rational_bound_checks += 1
                for qvalue in qvalues:
                    weighted_norm = sum(abs(evaluated_weights[j]) * qvalue ** j
                                        for j in range(degree + 1))
                    bound = hpoly * (1 + qvalue) ** degree / sample_spacing ** degree
                    require(sp.simplify(bound - weighted_norm).is_nonnegative is True,
                            "exact rational weighted H amplification bound failed")
                    rational_bound_checks += 1

            # The cutoff hard pair's Laurent coefficients exhibit the exact
            # additional factor L^(m-k)/(m-k)! relative to the top coefficient.
            hard_coefficient = sp.factorial(degree) * (-center) ** n / sp.factorial(n)
            polynomial_coefficient = (sp.factorial(k) * sp.binomial(degree, k)
                                      * (-center) ** n)
            require(sp.expand(hard_coefficient - polynomial_coefficient) == 0,
                    "cutoff hard-pair marginal Laurent coefficient failed")
            rows.append({"degree": degree, "coefficient_index_k": k,
                         "pole_power": -(k + 1), "L_degree": n,
                         "target": "k! * u_k",
                         "exact_inverse_heat_identity": True,
                         "hard_pair_log_factor": f"L^{n}/{factorial(n)}"})

    return {"all_coefficients_degrees_1_to_6": rows,
            "exact_rational_H_amplification_comparisons": rational_bound_checks,
            "H_bound_samples": [{"L": str(l), "d": str(d), "sigma": str(s)}
                                for l, d, s in samples],
            "geometric_weights_q": [str(qvalue) for qvalue in qvalues],
            "wrong_inverse_heat_sign_rejected": wrong_heat_sign}


def frequency_kernel(nu, *, a, h, sigma, spacing, degree):
    """Combined frequency kernel; its integral targets m! times a_m."""
    z = mp.mpc(a, nu)
    return (h ** (-z) * mp.exp(sigma ** 2 * z ** 2 / 2)
            * mp.expm1(spacing * z) ** degree
            / (2 * mp.pi * spacing ** degree))


def cutoff_transform(nu, *, a, cutoff, amplitude, degree):
    z = mp.mpc(a, nu)
    return amplitude * cutoff ** z * mp.factorial(degree) / z ** (degree + 1)


def positive_gaussian_moment(v, sigma, degree):
    """Independent closed recurrence for E[(v-sigma Z)_+^m]."""
    cdf = mp.erfc(-v / (mp.sqrt(2) * sigma)) / 2
    if degree == 0:
        return cdf
    density = mp.exp(-(v / sigma) ** 2 / 2) / mp.sqrt(2 * mp.pi)
    previous, current = cdf, v * cdf + sigma * density
    for k in range(2, degree + 1):
        previous, current = current, v * current + (k - 1) * sigma ** 2 * previous
    return current


def format_number(value) -> str:
    return mp.nstr(value, 18)


def jet_weights(degree, k, center, spacing, sigma):
    """Numerical version of the explicit finite-dimensional inverse-heat weights."""
    weights = []
    for j in range(degree + 1):
        polynomial = [mp.mpf(1)]
        denominator = mp.mpf(1)
        for i in range(degree + 1):
            if i == j:
                continue
            root = center + i * spacing
            multiplied = [mp.mpf(0)] * (len(polynomial) + 1)
            for power, coefficient in enumerate(polynomial):
                multiplied[power] -= root * coefficient
                multiplied[power + 1] += coefficient
            polynomial = multiplied
            denominator *= spacing * (j - i)
        weights.append(mp.fsum(
            (-sigma ** 2 / 2) ** r / mp.factorial(r) * mp.factorial(k + 2 * r)
            * polynomial[k + 2 * r] / denominator
            for r in range((degree - k) // 2 + 1)))
    return weights


def h_amplification_polynomial(n, x, sigma):
    return mp.fsum((sigma ** 2 / 2) ** r * x ** (n - 2 * r)
                   / (mp.factorial(r) * mp.factorial(n - 2 * r))
                   for r in range(n // 2 + 1))


def numerical_checks() -> dict:
    a, alpha, modulus = mp.mpf("1.1"), mp.mpf("0.8"), mp.mpf("2")
    sigma, spacing = mp.mpf("0.75"), mp.mpf("0.4")
    cutoff, h = mp.mpf("0.8"), mp.mpf("0.2")
    noise, band = mp.mpf("1e-7"), mp.mpf("3")
    gaussian_bias = mp.exp(alpha ** 2 * sigma ** 2 / 2)
    gaussian_noise = mp.exp(a ** 2 * sigma ** 2 / 2) / (mp.sqrt(2 * mp.pi) * sigma)
    rows = []
    minimax_rows = []
    jet_rows = []
    marginal_minimax_rows = []

    for degree in range(1, 4):
        envelope_max = (degree / (alpha * mp.e)) ** degree
        amplitude = modulus * cutoff ** alpha / envelope_max
        target = mp.factorial(degree) * amplitude
        qnorm = amplitude * cutoff ** a * mp.factorial(degree) / a ** (degree + 1)
        kernel = lambda nu: frequency_kernel(
            nu, a=a, h=h, sigma=sigma, spacing=spacing, degree=degree)
        transform = lambda nu: cutoff_transform(
            nu, a=a, cutoff=cutoff, amplitude=amplitude, degree=degree)

        # Hermitian symmetry makes the exact-data estimator real. The two
        # calculations below use frequency-space quadrature and log-energy
        # Gaussian truncated moments, respectively.
        full_estimate = 2 * mp.quad(lambda nu: mp.re(kernel(nu) * transform(nu)),
                                   [0, 1, 2, 4, 8, 16, mp.inf])
        gaussian_samples = [amplitude * positive_gaussian_moment(
            mp.log(cutoff / h) + j * spacing, sigma, degree)
            for j in range(degree + 1)]
        direct_estimate = mp.fsum((-1) ** (degree - j) * comb(degree, j)
                                  * gaussian_samples[j] for j in range(degree + 1)
                                  ) / spacing ** degree
        inversion_error = abs(full_estimate - direct_estimate)
        require(inversion_error < mp.mpf("1e-35"),
                f"direct Mellin/Gaussian identity failed for m={degree}")

        # A second energy integral checks the independent truncated-normal
        # recurrence used above for several actual shifted sample scales.
        moment_errors = []
        for j in range(degree + 1):
            v = mp.log(cutoff / h) + j * spacing
            energy_integral = mp.quad(
                lambda u: (v - sigma * u) ** degree
                * mp.exp(-u * u / 2) / mp.sqrt(2 * mp.pi),
                [-mp.inf, 0, v / sigma],
            )
            error = abs(energy_integral - positive_gaussian_moment(v, sigma, degree))
            require(error < mp.mpf("1e-40"), "truncated Gaussian moment failed")
            moment_errors.append(error)

        # The Mellin transform is checked against its defining log integral.
        mellin_errors = []
        for nu in (mp.mpf(0), mp.mpf("0.7"), mp.mpf("2")):
            z = mp.mpc(a, nu)
            log_integral = amplitude * cutoff ** z * mp.quad(
                lambda u: u ** degree * mp.exp(-z * u), [0, 1, 4, 12, mp.inf])
            error = abs(log_integral - transform(nu))
            require(error < mp.mpf("1e-35"), "cutoff Mellin integral check failed")
            mellin_errors.append(error)

        upper_bias = (gaussian_bias * modulus * h ** alpha
                      * (1 + mp.exp(-alpha * spacing)) ** degree / spacing ** degree)
        upper_noise_norm = (gaussian_noise * h ** (-a)
                            * (1 + mp.exp(a * spacing)) ** degree / spacing ** degree)
        actual_bias = abs(full_estimate - target)
        require(actual_bias <= upper_bias, "proved log-recovery bias bound failed")

        actual_noise_norm = 2 * mp.quad(lambda nu: abs(kernel(nu)),
                                       [0, 1, 2, 4, 8, 16, mp.inf])
        require(actual_noise_norm <= upper_noise_norm,
                "binomial frequency-noise bound failed")

        # This measurable adversarial noise has magnitude exactly eta and
        # aligns the combined kernel at every frequency.
        adversarial_integral = 2 * mp.quad(
            lambda nu: mp.re(kernel(nu) * noise * mp.conj(kernel(nu)) / abs(kernel(nu))),
            [0, 1, 2, 4, 8, 16, mp.inf],
        )
        require(abs(adversarial_integral - noise * actual_noise_norm) < mp.mpf("1e-35"),
                "adversarial-noise phase control failed")

        band_estimate = 2 * mp.quad(
            lambda nu: mp.re(kernel(nu) * transform(nu)), [0, 1, band])
        band_noise_norm = 2 * mp.quad(lambda nu: abs(kernel(nu)), [0, 1, band])
        omitted_data_error = abs(full_estimate - band_estimate)
        tail_bound = qnorm * upper_noise_norm * mp.erfc(sigma * band / mp.sqrt(2))
        require(omitted_data_error <= tail_bound,
                "finite-band exact-data tail bound failed")
        actual_band_error = abs(band_estimate + noise * band_noise_norm - target)
        full_band_bound = upper_bias + noise * upper_noise_norm + tail_bound
        require(actual_band_error <= full_band_bound,
                "finite-band recovery error bound failed")

        # For x=cutoff*exp(u), the ratio |r(x)|/x^alpha is
        # amplitude*cutoff^(-alpha)*u^m*exp(-alpha*u). Test the
        # proved maximum at u=m/alpha and values on either side.
        maximum_ratio = (amplitude * cutoff ** (-alpha)
                         * (degree / alpha) ** degree * mp.exp(-degree))
        require(abs(maximum_ratio - modulus) < mp.mpf("1e-40"),
                "sharp global remainder envelope normalization failed")
        for multiplier in (mp.mpf("0.1"), mp.mpf("0.5"), mp.mpf(1), 2, 5):
            u = multiplier * degree / alpha
            ratio = amplitude * cutoff ** (-alpha) * u ** degree * mp.exp(-alpha * u)
            require(ratio <= modulus * (1 + mp.mpf("1e-45")),
                    "global remainder envelope sample exceeds exact maximum")

        rows.append({
            "degree": degree,
            "target_highest_pole_coefficient": format_number(target),
            "Mellin_Gaussian_identity_absolute_error": format_number(inversion_error),
            "Gaussian_energy_moment_max_error": format_number(max(moment_errors)),
            "Mellin_defining_integral_max_error": format_number(max(mellin_errors)),
            "actual_bias": format_number(actual_bias),
            "proved_bias_upper_bound": format_number(upper_bias),
            "actual_combined_noise_kernel_L1": format_number(actual_noise_norm),
            "binomial_noise_kernel_upper_bound": format_number(upper_noise_norm),
            "adversarial_noise_error": format_number(adversarial_integral),
            "finite_band_omitted_data_error": format_number(omitted_data_error),
            "finite_band_tail_upper_bound": format_number(tail_bound),
            "finite_band_adversarial_total_error": format_number(actual_band_error),
            "finite_band_total_upper_bound": format_number(full_band_bound),
            "sharp_global_remainder_ratio": format_number(maximum_ratio),
        })

        center = -mp.log(h)
        profile_center = -mp.log(cutoff)
        for k in range(degree + 1):
            weights = jet_weights(degree, k, center, spacing, sigma)
            jet_estimate = mp.fsum(w * sample for w, sample in zip(weights, gaussian_samples))
            jet_target = (amplitude * mp.factorial(k) * comb(degree, k)
                          * (-profile_center) ** (degree - k))
            bias_weight_sum = mp.fsum(abs(weights[j]) * mp.exp(-alpha * j * spacing)
                                     for j in range(degree + 1))
            noise_weight_sum = mp.fsum(abs(weights[j]) * mp.exp(a * j * spacing)
                                      for j in range(degree + 1))
            hpoly = h_amplification_polynomial(degree - k,
                                               abs(center) + degree * spacing, sigma)
            bias_weight_bound = (hpoly / spacing ** degree
                                 * (1 + mp.exp(-alpha * spacing)) ** degree)
            noise_weight_bound = (hpoly / spacing ** degree
                                  * (1 + mp.exp(a * spacing)) ** degree)
            require(bias_weight_sum <= bias_weight_bound * (1 + mp.mpf("1e-45")),
                    "marginal numerical H bias amplification bound failed")
            require(noise_weight_sum <= noise_weight_bound * (1 + mp.mpf("1e-45")),
                    "marginal numerical H noise amplification bound failed")
            jet_bias = abs(jet_estimate - jet_target)
            jet_bias_bound = gaussian_bias * modulus * h ** alpha * bias_weight_bound
            require(jet_bias <= jet_bias_bound,
                    "cutoff profile full principal-part bias bound failed")
            jet_rows.append({
                "degree": degree, "coefficient_index_k": k,
                "target_pole_power": -(k + 1),
                "target": format_number(jet_target),
                "estimate_from_Gaussian_samples": format_number(jet_estimate),
                "actual_bias": format_number(jet_bias),
                "proved_H_bias_bound": format_number(jet_bias_bound),
                "exact_weight_noise_amplification_upper_bound": format_number(
                    gaussian_noise * h ** (-a) * noise_weight_sum),
                "H_noise_amplification_upper_bound": format_number(
                    gaussian_noise * h ** (-a) * noise_weight_bound),
            })

        theta = alpha / (a + alpha)
        previous_eta, previous_risk = None, None
        for eta in (mp.mpf("1e-3"), mp.mpf("1e-5"), mp.mpf("1e-7")):
            hard_cutoff = (2 * eta * a ** (degree + 1) * envelope_max
                           / (mp.factorial(degree) * modulus)) ** (1 / (a + alpha))
            hard_amplitude = modulus * hard_cutoff ** alpha / envelope_max
            hard_qnorm = (hard_amplitude * hard_cutoff ** a * mp.factorial(degree)
                          / a ** (degree + 1))
            risk_lower = mp.factorial(degree) * hard_amplitude / 2
            formula = ((mp.factorial(degree) * modulus / envelope_max)
                       ** (a / (a + alpha))
                       * (2 * eta * a ** (degree + 1)) ** theta / 2)
            require(abs(hard_qnorm / (2 * eta) - 1) < mp.mpf("1e-40"),
                    "hard pair line-norm normalization failed")
            require(abs(risk_lower / formula - 1) < mp.mpf("1e-40"),
                    "minimax hard pair constant failed")
            slope = None
            if previous_eta is not None:
                slope = mp.log(risk_lower / previous_risk) / mp.log(eta / previous_eta)
                require(abs(slope - theta) < mp.mpf("1e-40"),
                        "hard-pair minimax exponent check failed")
            minimax_rows.append({
                "degree": degree, "eta": format_number(eta),
                "line_norm": format_number(hard_qnorm),
                "residue_pair_half_separation": format_number(risk_lower),
                "expected_exponent": format_number(theta),
                "successive_log_slope": None if slope is None else format_number(slope),
            })
            if degree in (2, 3):
                hard_center = -mp.log(hard_cutoff)
                base_rate = modulus ** (a / (a + alpha)) * eta ** theta
                for k in range(degree + 1):
                    n = degree - k
                    marginal_risk = risk_lower * abs(hard_center) ** n / mp.factorial(n)
                    leading_constant = ((mp.factorial(degree) / envelope_max)
                                        ** (a / (a + alpha))
                                        * (2 * a ** (degree + 1)) ** theta
                                        / (2 * mp.factorial(n)))
                    expected_marginal = leading_constant * base_rate * abs(hard_center) ** n
                    require(abs(marginal_risk / expected_marginal - 1) < mp.mpf("1e-40"),
                            "full principal-part marginal hard-pair lower bound failed")
                    marginal_minimax_rows.append({
                        "degree": degree, "coefficient_index_k": k,
                        "eta": format_number(eta),
                        "hard_pair_L": format_number(hard_center),
                        "marginal_half_separation": format_number(marginal_risk),
                        "base_power_rate": format_number(base_rate),
                        "required_logarithmic_power": n,
                        "risk_divided_by_power_rate_and_L_power": format_number(
                            marginal_risk / (base_rate * abs(hard_center) ** n)),
                        "risk_divided_by_power_rate_only": format_number(marginal_risk / base_rate),
                    })
            previous_eta, previous_risk = eta, risk_lower

    return {
        "working_decimal_digits": mp.mp.dps,
        "parameters": {"a": str(a), "alpha": str(alpha), "M": str(modulus),
                       "sigma": str(sigma), "spacing": str(spacing),
                       "profile_cutoff": str(cutoff), "estimator_h": str(h),
                       "noise_eta": str(noise), "frequency_band_Omega": str(band)},
        "degrees_1_to_3": rows,
        "full_principal_part_cutoff_controls": jet_rows,
        "two_profile_minimax_scaling_controls": minimax_rows,
        "marginal_minimax_log_loss_controls": marginal_minimax_rows,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path,
                        default=Path(__file__).resolve().parents[1] / "reports" / "log-recovery.json")
    args = parser.parse_args()
    report = {
        "schema": "celestial-direct-log-recovery-v2",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "python_version": platform.python_version(),
        "sympy_version": sp.__version__, "mpmath_version": mp.__version__,
        "formal_proof": False,
        "scope": ("Exact Gaussian polynomial/finite-difference, full inverse-heat "
                  "principal-part weights, cutoff Mellin identities, and exact rational "
                  "H amplification comparisons for m=1..6; numerical direct-data, "
                  "noise, band, full principal-part, and minimax/log-loss controls "
                  "for m=1..3."),
        "exact": exact_checks(),
        "numerical": numerical_checks(),
        "passed": True,
        "limits": [
            "The manuscript, not finite computations, proves the general analytical statements.",
            "The fixed-degree global remainder prior is assumed, not inferred from data.",
            "The constants are upper/lower bounds; no equality or degree-uniform sharpness is asserted.",
            "These are scalar mathematical kernels, not a proved physical loop amplitude or resummation.",
            "Numerical quadrature values are controls and are not interval-arithmetic certificates.",
            "No Lean formalization of these analytical statements is supplied by this script.",
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print("PASS: symbolic Gaussian and full principal-part weights m=1..6; exact rational "
          "H bounds; direct Mellin, noise, band, full-jet, and marginal log-loss "
          f"controls m=1..3; {args.output}")


if __name__ == "__main__":
    main()
