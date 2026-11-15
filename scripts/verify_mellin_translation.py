"""Finite-band translation controls using an explicit Mellin wave packet."""

from __future__ import annotations

import json
from pathlib import Path

import mpmath as mp


ROOT = Path(__file__).resolve().parents[1]
mp.mp.dps = 60


def density(nu, decay, time):
    angle = mp.atan(time / decay)
    radius = mp.sqrt(decay**2 + time**2)
    return decay * mp.exp(2 * nu * angle) / (radius * mp.cosh(mp.pi * nu))


def report_record(decay, time, band_half_width, projection_distance):
    retained = mp.sqrt(
        mp.quad(
            lambda nu: density(nu, decay, time),
            [-band_half_width, 0, band_half_width],
        )
    )
    # The displayed density is for the exact full-line formula. Applying the
    # contraction P_B U_tau transfers its retained norm to the normalized
    # exact finite-band projection with error at most projection_distance.
    return {
        "decay_parameter": mp.nstr(decay, 18),
        "time": mp.nstr(time, 18),
        "unprojected_input_retained_norm": mp.nstr(retained, 18),
        "projected_input_lower_bound": mp.nstr(
            max(mp.mpf(0), retained - projection_distance), 18
        ),
        "projected_input_upper_bound": mp.nstr(
            min(mp.mpf(1), retained + projection_distance), 18
        ),
    }


def main():
    quadrature_controls = []
    for nu, time in [(mp.mpf("-0.5"), 1), (mp.mpf("0.5"), 1), (mp.mpf("0.5"), -1)]:
        # h(x)=sqrt(2)*exp(x/2-exp(x)) has unit L² norm. Numerically integrate
        # its translated log-energy Fourier transform in the energy variable.
        amplitude = mp.quad(
            lambda w: w ** (-mp.mpf("0.5") - 1j * nu) * mp.exp(-(1 - 1j * time) * w),
            [0, 1, mp.inf],
        ) / mp.sqrt(mp.pi)
        error = abs(abs(amplitude)**2 - density(nu, mp.mpf(1), time))
        assert error < mp.mpf("1e-27"), error
        quadrature_controls.append(
            {"nu": str(nu), "time": time, "density_absolute_error": mp.nstr(error, 12)}
        )

    band_half_width = mp.mpf(8)
    initial_norm_squared = 2 * mp.atan(mp.sinh(mp.pi * band_half_width)) / mp.pi
    projection_distance = mp.sqrt(2 * (1 - mp.sqrt(initial_norm_squared)))
    low_energy = report_record(10**6, 1, band_half_width, projection_distance)
    high_energy = report_record(mp.mpf("1e-10"), 1, band_half_width, projection_distance)
    assert mp.mpf(low_energy["projected_input_lower_bound"]) > mp.mpf("0.99999")
    assert mp.mpf(high_energy["projected_input_upper_bound"]) < mp.mpf("0.00006")

    report = {
        "status": "passed",
        "scope": "Independent density quadrature and two retained-norm controls for the exact finite-band projection of a massless scalar wave packet.",
        "band_half_width": str(band_half_width),
        "normalized_projection_distance": mp.nstr(projection_distance, 18),
        "low_energy_case": low_energy,
        "high_energy_case": high_energy,
        "quadrature_controls": quadrature_controls,
    }
    destination = ROOT / "reports" / "mellin-translation.json"
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print("PASS: finite-band translation norm controls")


if __name__ == "__main__":
    main()
