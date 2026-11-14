"""Check the classical two-time-channel Fourier--Mellin calculation.

The theorem uses smooth compact support away from zero. The independent
quadrature control f(omega)=exp(-omega) has an explicit Fourier transform
and satisfies the displayed bridge on 0<Re(s)<1, although it is not compact.
Finite checks support the formulas; the general proofs are in the paper.
"""

from __future__ import annotations

import json
from pathlib import Path

import mpmath as mp


ROOT = Path(__file__).resolve().parents[1]
mp.mp.dps = 70


def number(value):
    return mp.nstr(value, 30)


def record(value):
    return {"real": number(mp.re(value)), "imag": number(mp.im(value))}


def kernels(nu):
    s = mp.mpf("0.5") + 1j * nu
    return s, mp.gamma(s) * mp.exp(-1j * mp.pi * s / 2), mp.gamma(s) * mp.exp(1j * mp.pi * s / 2)


def main():
    identities = []
    sharpness = []
    negative_controls = []
    for nu in map(mp.mpf, ["-5", "-1", "0", "1", "5"]):
        s, kp, km = kernels(nu)
        joint = abs(kp) ** 2 + abs(km) ** 2
        error = abs(joint - 2 * mp.pi)
        assert error < mp.mpf("1e-60")
        identities.append({"nu": number(nu), "joint_squared_norm": number(joint), "error_from_2pi": number(error)})

        # The perturbation lies along the response vector and saturates Cauchy--Schwarz.
        delta = mp.mpf("0.03")
        ep, em = delta * kp / mp.sqrt(2 * mp.pi), delta * km / mp.sqrt(2 * mp.pi)
        target = mp.gamma(1 - s)
        estimate = (mp.conj(kp) * (kp * target + ep) + mp.conj(km) * (km * target + em)) / (2 * mp.pi)
        exact_error = delta / mp.sqrt(2 * mp.pi)
        assert abs(abs(estimate - target) - exact_error) < mp.mpf("1e-60")
        sharpness.append({"nu": number(nu), "joint_noise_norm": number(mp.sqrt(abs(ep) ** 2 + abs(em) ** 2)), "observed_target_error": number(abs(estimate - target)), "sharp_bound": number(exact_error)})

        for sign, kernel, expected in [
            ("positive time", kp, (1 + mp.exp(-2 * mp.pi * nu)) / (2 * mp.pi)),
            ("negative time", km, (1 + mp.exp(2 * mp.pi * nu)) / (2 * mp.pi)),
        ]:
            assert abs(abs(kernel) ** -2 / expected - 1) < mp.mpf("1e-60")
            negative_controls.append({"channel": sign, "nu": number(nu), "inverse_amplification": number(1 / abs(kernel)), "two_channel_amplification": number(1 / mp.sqrt(2 * mp.pi))})

    quadratures = []
    for nu in map(mp.mpf, ["-0.2", "0", "0.2", "2"]):
        s, kp, km = kernels(nu)
        for sign, kernel in [(1, kp), (-1, km)]:
            # C_+(u)=1/(1+iu), C_-(u)=1/(1-iu), F(1-s)=Gamma(1-s).
            actual = mp.quad(lambda u: u ** (s - 1) / (1 + sign * 1j * u), [0, 1, 16, mp.inf])
            expected = kernel * mp.gamma(1 - s)
            error = abs(actual - expected)
            assert error < mp.mpf("1e-28"), error
            quadratures.append({"nu": number(nu), "time_sign": sign, "quadrature": record(actual), "formula": record(expected), "absolute_error": number(error)})

    report = {
        "status": "passed",
        "scope": "Classical kinematic Fourier--Mellin identity on a specified data line, not a new duality or a physical charge theorem.",
        "precision_decimal_digits": mp.mp.dps,
        "theorem_class": "complex-valued C_c^infinity((0,infinity)), support away from zero",
        "quadrature_control": "f(omega)=exp(-omega); explicit C_plus/minus=1/(1 +/- i*u); valid bridge on 0<Re(s)<1, outside compact theorem class",
        "joint_norm_identity": identities,
        "sharp_linear_noise_bound": sharpness,
        "one_channel_negative_controls": negative_controls,
        "independent_quadratures": quadratures,
        "limitations": "Full-line two-time-channel availability is assumed. Re(q)=1/2 is not automatically the physical Re(Delta)=1 line. No quantum state space, Ward identity, OPE, bulk dual, or entanglement formula is established.",
    }
    destination = ROOT / "reports" / "framework-bridge.json"
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(f"PASS: {len(identities)} joint-norm identities, {len(sharpness)} sharpness checks, {len(negative_controls)} single-channel controls, {len(quadratures)} independent quadratures")
    print(destination)


if __name__ == "__main__":
    main()
