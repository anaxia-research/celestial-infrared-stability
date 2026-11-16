#!/usr/bin/env python3
"""Numerical and symbolic controls for known-branch Mellin recovery.

The analytical theorem is in the paper (mdpi-preprints/CelestialIRPaper.pdf). Finite executions
support its implementation; numerical quadrature is not a formal or interval
proof, and the prescribed template is not a general physical RG solution.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import platform
import sys

import mpmath as mp
import sympy as sp


if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")
mp.mp.dps = 45


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def number(value):
    return mp.nstr(value, 18)


def template(t, lam, gamma):
    return (1 + lam * t) ** (-gamma)


def transform(z, lam, gamma, lower=mp.mpf(0)):
    """Laplace integral of K(t) on [lower,infinity), principal branch."""
    return (lam ** (-gamma) * mp.exp(z / lam) * z ** (gamma - 1)
            * mp.gammainc(1 - gamma, z * (lower + 1 / lam), mp.inf))


def response(center, sigma, lam, gamma):
    normalizer = mp.sqrt(2 * mp.pi) * sigma
    integrand = lambda t: (mp.exp(-(t - center) ** 2 / (2 * sigma ** 2))
                          * template(t, lam, gamma) / normalizer)
    if center > 0:
        return mp.quad(integrand, [0, center, center + 8 * sigma, mp.inf])
    return mp.quad(integrand, [0, 8 * sigma, mp.inf])


def kernel(nu, a, center, sigma):
    z = mp.mpc(a, nu)
    return mp.exp(center * z + sigma ** 2 * z ** 2 / 2) / (2 * mp.pi)


def exact_checks():
    t, alpha, lam, gamma = sp.symbols("t alpha lambda gamma", positive=True)
    log_envelope = -alpha * t + gamma * sp.log(1 + lam * t)
    second = sp.diff(log_envelope, t, 2)
    expected = -gamma * lam ** 2 / (1 + lam * t) ** 2
    require(sp.simplify(second - expected) == 0,
            "hard-profile logarithmic envelope concavity identity failed")
    require(sp.ask(sp.Q.negative(expected)) is True,
            "hard-profile logarithmic envelope sign was not established")
    return {"log_F_second_derivative": str(expected),
            "negative_for_positive_parameters": True,
            "scope": "Symbolic concavity identity; endpoint argument is proved in the manuscript."}


def numerical_checks():
    a, alpha, modulus = mp.mpf("1.1"), mp.mpf("0.8"), mp.mpf("2")
    lam, gamma, sigma = mp.mpf("0.7"), mp.mpf("0.4"), mp.mpf("0.75")
    amplitude, remainder_amplitude = mp.mpf("1.25"), mp.mpf("0.7")
    center, noise, band = mp.mpf("2"), mp.mpf("1e-6"), mp.mpf("3")
    p0 = mp.erf(1 / mp.sqrt(2))
    abias = mp.exp(alpha ** 2 * sigma ** 2 / 2)
    bnoise = mp.exp(a ** 2 * sigma ** 2 / 2) / (mp.sqrt(2 * mp.pi) * sigma)
    cbranch = lam ** (-gamma) * mp.gamma(1 - gamma)

    transform_rows = []
    for real_part, frequency in (("0.3", "0"), ("1.1", "0.7"), ("1.1", "2")):
        z = mp.mpc(real_part, frequency)
        defining_integral = mp.quad(
            lambda t: mp.exp(-z * t) * template(t, lam, gamma),
            [0, 1, 4, 12, mp.inf])
        error = abs(defining_integral - transform(z, lam, gamma))
        require(error < mp.mpf("1e-32"), "incomplete-gamma Mellin formula failed")
        transform_rows.append({"real_part": real_part, "frequency": frequency,
                               "absolute_error": number(error)})

    branch_rows = []
    for small_z in (mp.mpf("1e-3"), mp.mpf("1e-6"), mp.mpf("1e-9")):
        normalized = small_z ** (1 - gamma) * transform(small_z, lam, gamma)
        relative_error = abs(normalized / cbranch - 1)
        require(relative_error < 2 * small_z ** (1 - gamma),
                "leading branch coefficient control failed")
        branch_rows.append({"z": number(small_z),
                            "normalized_branch_value": number(normalized),
                            "expected_branch_coefficient": number(cbranch),
                            "relative_error": number(relative_error)})

    response_rows = []
    for location in (sigma, mp.mpf(2), mp.mpf(5), mp.mpf(10)):
        phi = response(location, sigma, lam, gamma)
        lower = p0 * template(location + sigma, lam, gamma)
        require(phi >= lower, "known-template Gaussian response lower bound failed")
        response_rows.append({"L": number(location), "Phi": number(phi),
                              "proved_lower_bound": number(lower),
                              "ratio_to_lower_bound": number(phi / lower)})

    phi = response(center, sigma, lam, gamma)
    # Choose a remainder with an independent closed Mellin transform and
    # the exact global bound |R(x)| <= remainder_amplitude*x^alpha <= M*x^alpha.
    # R(x)=r0*x^alpha*exp(-x), including x>1.
    def full_transform(nu):
        z = mp.mpc(a, nu)
        return (amplitude * transform(z, lam, gamma)
                + remainder_amplitude * mp.gamma(z + alpha))

    qnorm = (amplitude * transform(a, lam, gamma)
             + remainder_amplitude * mp.gamma(a + alpha))
    cutoff_band = mp.mpf(14) / sigma
    integration_points = [0, 1, 2, 4, 8, 12, cutoff_band]
    frequency_average = 2 * mp.quad(
        lambda nu: mp.re(kernel(nu, a, center, sigma) * full_transform(nu)),
        integration_points)
    # Finite Gaussian panels avoid evaluating exp(-exp(sigma*z-L))
    # at the extremely large nodes of infinite-interval quadrature.
    z_limit = mp.mpf(18)
    gaussian_remainder = mp.quad(
        lambda z: (remainder_amplitude * mp.exp(alpha * (-center + sigma * z))
                   * mp.exp(-mp.exp(-center + sigma * z))
                   * mp.exp(-z * z / 2) / mp.sqrt(2 * mp.pi)),
        [-z_limit, -8, 0, 8, z_limit])
    tilt = alpha * sigma
    energy_omission_bound = (
        remainder_amplitude * mp.exp(-alpha * center + tilt**2 / 2)
        * (mp.erfc((z_limit - tilt) / mp.sqrt(2))
           + mp.erfc((z_limit + tilt) / mp.sqrt(2))) / 2
    )
    energy_average = amplitude * phi + gaussian_remainder
    frequency_omission_bound = (qnorm * bnoise * mp.exp(a * center)
                                * mp.erfc(sigma * cutoff_band / mp.sqrt(2)))
    inversion_error = abs(frequency_average - energy_average)
    require(inversion_error < mp.mpf("1e-28") + frequency_omission_bound
            + energy_omission_bound,
            "known-branch Mellin/Gaussian average identity failed")

    estimate = frequency_average / phi
    actual_bias = abs(estimate - amplitude)
    bias_bound = abias * modulus * mp.exp(-alpha * center) / phi
    noise_norm = bnoise * mp.exp(a * center) / phi
    require(actual_bias <= bias_bound + frequency_omission_bound / phi,
            "normalized known-branch bias bound failed")
    # Phase-aligned complex noise saturates the complete Gaussian kernel norm.
    adversarial_integral = 2 * mp.quad(
        lambda nu: abs(kernel(nu, a, center, sigma)) * noise / phi,
        [0, 1, 2, 4, 8, 16, mp.inf])
    require(abs(adversarial_integral / (noise * noise_norm) - 1) < mp.mpf("1e-32"),
            "normalized adversarial-noise amplification failed")

    band_estimate = 2 * mp.quad(
        lambda nu: mp.re(kernel(nu, a, center, sigma) * full_transform(nu)),
        [0, 1, band]) / phi
    band_noise = 2 * mp.quad(
        lambda nu: abs(kernel(nu, a, center, sigma)) * noise / phi,
        [0, 1, band])
    band_tail_bound = qnorm * noise_norm * mp.erfc(sigma * band / mp.sqrt(2))
    actual_band_tail = abs(frequency_average / phi - band_estimate)
    require(actual_band_tail <= band_tail_bound + frequency_omission_bound / phi,
            "normalized finite-band tail bound failed")
    band_total_bound = bias_bound + noise * noise_norm + band_tail_bound
    actual_band_error = abs(band_estimate + band_noise - amplitude)
    require(actual_band_error <= band_total_bound,
            "normalized finite-band total recovery bound failed")

    hard_rows = []
    theta = alpha / (a + alpha)
    for eta in (mp.mpf("1e-3"), mp.mpf("1e-5"), mp.mpf("1e-7")):
        location = mp.log(modulus / (2 * a * eta)) / (a + alpha)
        f_end = mp.exp(-alpha * location) / template(location, lam, gamma)
        require(location >= 0 and f_end <= 1, "hard pair small-noise conditions fail")
        b_hard = modulus * f_end
        q_hard = b_hard * transform(a, lam, gamma, lower=location)
        norm_upper = modulus * mp.exp(-(a + alpha) * location) / a
        require(q_hard <= norm_upper * (1 + mp.mpf("1e-35")),
                "hard pair exact norm exceeds its upper bound")
        require(abs(norm_upper / (2 * eta) - 1) < mp.mpf("1e-35"),
                "hard-pair midpoint normalization failed")
        require(q_hard / 2 <= eta, "hard-pair common observation leaves noise balls")
        maximum_ratio = 0
        for fraction in (0, mp.mpf("0.1"), mp.mpf("0.5"), mp.mpf("0.9"), 1):
            t = fraction * location
            ratio = b_hard * template(t, lam, gamma) * mp.exp(alpha * t)
            maximum_ratio = max(maximum_ratio, ratio)
            require(ratio <= modulus * (1 + mp.mpf("1e-35")),
                    "hard pair global remainder bound sample failed")
        tail_error = abs(transform(mp.mpc(a, "0.7"), lam, gamma, lower=location)
                         - mp.quad(lambda u: mp.exp(-mp.mpc(a, "0.7") * (location + u))
                                   * template(location + u, lam, gamma),
                                   [0, 1, 4, 12, mp.inf]))
        require(tail_error < mp.mpf("1e-32"), "hard cutoff Mellin formula failed")
        power_rate = modulus ** (a / (a + alpha)) * eta ** theta
        predicted_lower = ((2 * a) ** theta / 2 * power_rate
                           / template(location, lam, gamma))
        require(abs((b_hard / 2) / predicted_lower - 1) < mp.mpf("1e-35"),
                "hard-pair power times log-gamma lower rate failed")
        hard_rows.append({
            "eta": number(eta), "L": number(location),
            "F_L": number(f_end), "amplitude": number(b_hard),
            "exact_line_norm": number(q_hard),
            "line_norm_upper_bound_2eta": number(norm_upper),
            "common_midpoint_noise_norm": number(q_hard / 2),
            "global_remainder_maximum_sample": number(maximum_ratio),
            "hard_cutoff_Mellin_absolute_error": number(tail_error),
            "risk_lower": number(b_hard / 2), "base_power_rate": number(power_rate),
            "risk_over_power_rate": number(b_hard / (2 * power_rate)),
            "risk_over_power_and_known_log_factor": number(
                b_hard * template(location, lam, gamma) / (2 * power_rate)),
        })

    # Without an amplitude lower bound, different exponents cannot be
    # estimated uniformly. This checks a bounded-noise hard-pair example.
    exponent_noise = mp.mpf("1e-7")
    small_amplitude = a * exponent_noise
    exponents = (mp.mpf("0.25"), mp.mpf("0.75"))
    exponent_differences = [small_amplitude * abs(
        transform(mp.mpc(a, nu), lam, exponents[0])
        - transform(mp.mpc(a, nu), lam, exponents[1]))
        for nu in (0, mp.mpf("0.7"), 2)]
    require(max(exponent_differences) <= 2 * exponent_noise,
            "unknown-exponent vanishing-amplitude control failed")

    return {
        "working_decimal_digits": mp.mp.dps,
        "parameters": {"a": str(a), "alpha": str(alpha), "M": str(modulus),
                       "lambda": str(lam), "gamma": str(gamma), "sigma": str(sigma)},
        "incomplete_gamma_Mellin_controls": transform_rows,
        "leading_fractional_branch_controls": branch_rows,
        "known_template_response_lower_controls": response_rows,
        "direct_recovery_control": {
            "L": number(center), "Phi": number(phi),
            "actual_amplitude": number(amplitude),
            "frequency_integral_band_cutoff": number(cutoff_band),
            "omitted_frequency_bound": number(frequency_omission_bound),
            "omitted_energy_gaussian_tail_bound": number(energy_omission_bound),
            "Mellin_Gaussian_average_absolute_error": number(inversion_error),
            "normalized_estimate": number(estimate),
            "actual_bias": number(actual_bias), "proved_bias_bound": number(bias_bound),
            "noise_eta": number(noise),
            "phase_aligned_noise_error": number(adversarial_integral),
            "exact_full_noise_amplification": number(noise_norm),
            "finite_band_Omega": number(band),
            "finite_band_omitted_data_error": number(actual_band_tail),
            "finite_band_tail_bound": number(band_tail_bound),
            "finite_band_adversarial_total_error": number(actual_band_error),
            "finite_band_total_upper_bound": number(band_total_bound),
        },
        "hard_cutoff_and_log_gamma_scaling_controls": hard_rows,
        "unknown_exponent_control": {
            "gamma_pair": [str(g) for g in exponents],
            "positive_amplitude": number(small_amplitude),
            "noise_eta": number(exponent_noise),
            "sampled_data_differences": [number(v) for v in exponent_differences],
            "proved_whole_line_difference_upper_bound": number(2 * exponent_noise),
            "exact_zero_amplitude_nonidentifiability": True,
        },
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path,
                        default=Path(__file__).resolve().parents[1] / "reports" / "branch-recovery.json")
    args = parser.parse_args()
    report = {
        "schema": "celestial-known-branch-recovery-v1",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "python_version": platform.python_version(),
        "sympy_version": sp.__version__, "mpmath_version": mp.__version__,
        "formal_proof": False,
        "scope": "Known fractional-log template, globally controlled remainder, and bounded-noise Mellin line.",
        "exact": exact_checks(), "numerical": numerical_checks(), "passed": True,
        "limits": [
            "The analytical theorem is proved in the manuscript; finite controls do not replace it.",
            "Numerical quadrature is not an interval-arithmetic or Lean certificate.",
            "The template parameters lambda and gamma are known; no universal unknown-exponent recovery is asserted.",
            "The global remainder prior is assumed and is not derived from a quantum field theory.",
            "No general resummed celestial amplitude, Ward algebra, or physical RG closure is established.",
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print("PASS: known-branch Mellin formula, normalized Gaussian recovery, noise and band bounds, "
          f"hard-pair log-gamma scaling, and exponent negative controls; {args.output}")


if __name__ == "__main__":
    main()
