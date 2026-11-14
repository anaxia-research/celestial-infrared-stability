"""Replay two elementary resummation scope controls using mpmath.

These scalar Laplace--Mellin kernels illustrate established analytic structures.
They are not computed scattering amplitudes or proofs about all physical loops.
Analytic identities and asymptotic arguments are given in paper/angular-logs.tex;
this finite numerical replay checks formula implementations and negative controls.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path

import mpmath as mp


ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "reports" / "resummation.json"
mp.mp.dps = 70
TOLERANCE = mp.mpf("1e-45")


def num(value):
    return mp.nstr(value, 28)


def complex_record(value):
    return {"real": num(mp.re(value)), "imag": num(mp.im(value))}


def relative_error(actual, expected):
    return abs(actual - expected) / max(mp.mpf(1), abs(expected))


def checked_integral(integrand, z, expected):
    # Splitting protects moderate-frequency complex examples from one giant panel.
    actual = mp.quad(lambda t: mp.exp(-z * t) * integrand(t), [0, 1, 4, 16, mp.inf])
    error = relative_error(actual, expected)
    assert error < TOLERANCE, f"Laplace formula mismatch: {error}"
    return {
        "z": complex_record(z),
        "quadrature": complex_record(actual),
        "closed_form": complex_record(expected),
        "scaled_absolute_error": num(error),
        "passed": True,
    }


def exponential_polynomial(t, lam, degree):
    return sum((-lam * t) ** j / mp.factorial(j) for j in range(degree + 1))


def binomial_polynomial(t, lam, gamma, degree):
    return sum(
        (-1) ** j * mp.rf(gamma, j) * (lam * t) ** j / mp.factorial(j)
        for j in range(degree + 1)
    )


def branch_transform(z, lam, gamma):
    return (
        lam ** (-gamma)
        * z ** (gamma - 1)
        * mp.exp(z / lam)
        * mp.gammainc(1 - gamma, z / lam, mp.inf)
    )


def negative_control(kernel, polynomial, alpha):
    times = [10, 20, 40, 80]
    ratios = [abs(kernel(mp.mpf(t)) - polynomial(mp.mpf(t))) * mp.exp(alpha * t) for t in times]
    assert all(right > left for left, right in zip(ratios, ratios[1:]))
    assert ratios[-1] > mp.mpf("1e12")
    return {
        "alpha": num(alpha),
        "t": times,
        "abs_remainder_times_exp_alpha_t": [num(v) for v in ratios],
        "finite_growth_diagnostic_passed": True,
        "analytic_reason": "A nonzero finite polynomial minus a decaying kernel cannot be O(exp(-alpha*t)) as t tends to infinity; finite samples illustrate this obstruction.",
    }


def main():
    lam = mp.mpf("0.7")
    z_values = [mp.mpc("0.7", "0.4"), mp.mpc("1.5", "-0.3"), mp.mpc("0.2", "0.1")]
    shifted = [checked_integral(lambda t: mp.exp(-lam * t), z, 1 / (z + lam)) for z in z_values]

    shifted_truncations = []
    branch_truncations = []
    gamma = mp.mpf("0.4")
    z = mp.mpc("0.8", "0.3")
    tiny_z = mp.mpf("1e-18")
    for degree in range(6):
        exp_closed = lambda zz: sum((-lam) ** j / zz ** (j + 1) for j in range(degree + 1))
        branch_closed = lambda zz: sum(
            (-1) ** j * mp.rf(gamma, j) * lam ** j / zz ** (j + 1)
            for j in range(degree + 1)
        )
        exp_check = checked_integral(lambda t: exponential_polynomial(t, lam, degree), z, exp_closed(z))
        branch_check = checked_integral(lambda t: binomial_polynomial(t, lam, gamma, degree), z, branch_closed(z))
        # Nonzero coefficients of z^(-degree-1) establish the exact pole orders.
        exp_top = (-lam) ** degree
        branch_top = (-1) ** degree * mp.rf(gamma, degree) * lam ** degree
        assert exp_top != 0 and branch_top != 0
        exp_check.update({"degree": degree, "pole_order_at_zero": degree + 1, "top_laurent_coefficient": num(exp_top), "scaled_near_zero_value": num(tiny_z ** (degree + 1) * exp_closed(tiny_z))})
        branch_check.update({"degree": degree, "pole_order_at_zero": degree + 1, "top_laurent_coefficient": num(branch_top), "scaled_near_zero_value": num(tiny_z ** (degree + 1) * branch_closed(tiny_z))})
        assert abs(tiny_z ** (degree + 1) * exp_closed(tiny_z) - exp_top) < mp.mpf("1e-15")
        assert abs(tiny_z ** (degree + 1) * branch_closed(tiny_z) - branch_top) < mp.mpf("1e-15")
        shifted_truncations.append(exp_check)
        branch_truncations.append(branch_check)

    branch_integrals = []
    limits = []
    for gamma in [mp.mpf("0.25"), mp.mpf("0.4"), mp.mpf("0.75")]:
        for z in z_values:
            check = checked_integral(lambda t: (1 + lam * t) ** (-gamma), z, branch_transform(z, lam, gamma))
            check["gamma"] = num(gamma)
            branch_integrals.append(check)
        target = lam ** (-gamma) * mp.gamma(1 - gamma)
        z_sequence = [mp.mpf(10) ** (-power) for power in [2, 4, 8, 16, 24]]
        values = [zz ** (1 - gamma) * branch_transform(zz, lam, gamma) for zz in z_sequence]
        errors = [abs(value / target - 1) for value in values]
        assert all(right < left for left, right in zip(errors, errors[1:]))
        assert errors[-1] < mp.mpf("5e-6")
        limits.append({"gamma": num(gamma), "z": [num(zz) for zz in z_sequence], "scaled_transform": [num(v) for v in values], "limit_target": num(target), "relative_errors": [num(v) for v in errors], "convergence_diagnostic_passed": True})

    alpha = mp.mpf("0.5")
    gamma = mp.mpf("0.4")
    controls = []
    for degree in [0, 1, 3]:
        controls.append({"kernel": "exp(-lambda*t)", "truncation_degree": degree, **negative_control(lambda t: mp.exp(-lam * t), lambda t: exponential_polynomial(t, lam, degree), alpha)})
        controls.append({"kernel": "(1+lambda*t)^(-gamma)", "gamma": num(gamma), "truncation_degree": degree, **negative_control(lambda t: (1 + lam * t) ** (-gamma), lambda t: binomial_polynomial(t, lam, gamma, degree), alpha)})

    report = {
        "schema": "celestial-resummation-controls-v1",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "formal_proof": False,
        "passed": True,
        "status": "passed",
        "scope": "Finite numerical verification of two scalar established-model illustrations; no physical loop amplitude, all-orders soft algebra, or Ward identity is claimed.",
        "precision_decimal_digits": mp.mp.dps,
        "integral_scaled_absolute_tolerance": num(TOLERANCE),
        "parameters": {"lambda": num(lam)},
        "shifted_pole": {"exact_formula": "1/(z+lambda)", "pole_location": num(-lam), "integral_checks": shifted, "finite_taylor_truncations": shifted_truncations},
        "branch_model": {"exact_formula": "lambda^(-gamma)*z^(gamma-1)*exp(z/lambda)*Gamma(1-gamma,z/lambda)", "domain_of_integral": "Re(z)>0, lambda>0, 0<gamma<1", "continuation": "principal branches; fractional branch point at z=0", "integral_checks": branch_integrals, "positive_axis_limit_diagnostics": limits, "finite_binomial_truncations": branch_truncations},
        "invalid_global_exponential_remainder_negative_controls": controls,
        "proof_boundary": "Numerical limits and sampled growth are diagnostics; substitution, finite polynomial integration, domination, and nondecaying polynomial asymptotics justify the corresponding general statements in the manuscript.",
    }
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(f"PASS: {len(shifted) + len(shifted_truncations) + len(branch_integrals) + len(branch_truncations)} integral checks, {len(limits)} limit diagnostics, {len(controls)} invalid-remainder controls")
    print(REPORT)


if __name__ == "__main__":
    main()
