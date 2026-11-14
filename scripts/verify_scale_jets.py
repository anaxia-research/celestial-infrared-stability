#!/usr/bin/env python3
"""Exact finite replay of local-unit transformations of Mellin pole jets.

All arithmetic is rational and uses only the Python standard library.  A jet
is ordered (a_-1, ..., a_-r).  Multiplication of its Laurent principal part by
exp(-s*z) gives T(s) = exp(-s*N), where N shifts the coefficient index upward.
The scalar factor exp(-(w+Delta0)*s) is tracked by its formal exponent.

For the transform M(Delta)=integral omega**Delta a(omega) d(omega)/omega,
omega'=exp(-s)*omega gives M'=exp(-(w+Delta)*s)*M.  A dimensionful reference
mu transforms as mu'=exp(-s)*mu.  Thus mu**(-Delta)*M has only weight w:
the compensator cancels the entire nilpotent jet transformation.

These finite identities check bookkeeping; they assert no scattering dynamics,
new gauge interaction, or novelty.  The accompanying mathematical derivation
supplies the general identities, beyond the finite rational cases checked here.
"""

from __future__ import annotations

import argparse
from fractions import Fraction as F
import json
from math import factorial
from pathlib import Path


def rational(value: F) -> str:
    """Serialize exactly; never turn a rational into a floating-point number."""
    return str(value)


def identity(order: int) -> list[list[F]]:
    return [[F(i == j) for j in range(order)] for i in range(order)]


def jet_matrix(order: int, shift: F) -> list[list[F]]:
    """The upper triangular matrix exp(-shift*N)."""
    return [
        [F(0) if j < i else (-shift) ** (j - i) / factorial(j - i)
         for j in range(order)]
        for i in range(order)
    ]


def multiply(left: list[list[F]], right: list[list[F]]) -> list[list[F]]:
    order = len(left)
    return [[sum((left[i][k] * right[k][j] for k in range(order)), F(0))
             for j in range(order)] for i in range(order)]


def apply(matrix: list[list[F]], vector: list[F]) -> list[F]:
    return [sum((a * b for a, b in zip(row, vector)), F(0)) for row in matrix]


def laurent_rescale(vector: list[F], shift: F) -> list[F]:
    """Independent Laurent-polynomial multiplication, keeping negative powers.

    The truncated exponential coefficients are built recursively rather than
    using the closed matrix formula.  Terms of degree >=r cannot affect the
    principal part of a pole of order r.
    """
    order = len(vector)
    exponential = [F(1)]
    for degree in range(1, order):
        exponential.append(exponential[-1] * (-shift) / degree)
    laurent = {-j - 1: coefficient for j, coefficient in enumerate(vector)}
    product: dict[int, F] = {}
    for power, coefficient in laurent.items():
        for degree, exp_coefficient in enumerate(exponential):
            resulting_power = power + degree
            if resulting_power < 0:
                product[resulting_power] = (
                    product.get(resulting_power, F(0))
                    + coefficient * exp_coefficient
                )
    return [product.get(-j - 1, F(0)) for j in range(order)]


def verify(max_order: int) -> dict:
    if not 1 <= max_order <= 50:
        raise ValueError("max_order must be between 1 and 50")
    shifts = (F(-3, 2), F(-1, 5), F(0), F(2, 3), F(5, 2))
    log_references = (F(-7, 4), F(0), F(1, 3), F(9, 5))
    group_cases = inverse_cases = laurent_cases = compensator_cases = 0
    scalar_cases = 0

    for order in range(1, max_order + 1):
        ident = identity(order)
        vectors = [row[:] for row in ident]
        vectors += [[F((-1) ** j * (j + 1), j + 2) for j in range(order)],
                    [F(0)] * order]
        for shift in shifts:
            transform = jet_matrix(order, shift)
            if multiply(transform, jet_matrix(order, -shift)) != ident:
                raise AssertionError(("inverse", order, shift))
            inverse_cases += 1
            for second_shift in shifts:
                if multiply(jet_matrix(order, second_shift), transform) != jet_matrix(
                    order, shift + second_shift
                ):
                    raise AssertionError(("group", order, shift, second_shift))
                group_cases += 1
            for vector in vectors:
                if apply(transform, vector) != laurent_rescale(vector, shift):
                    raise AssertionError(("Laurent product", order, shift, vector))
                laurent_cases += 1
                for ell in log_references:
                    # ell=log(mu); after changing units ell'=ell-shift.
                    normalized_after = apply(
                        jet_matrix(order, ell - shift), apply(transform, vector)
                    )
                    normalized_before = apply(jet_matrix(order, ell), vector)
                    if normalized_after != normalized_before:
                        raise AssertionError(("compensator", order, shift, ell, vector))
                    compensator_cases += 1

    # Remaining scalar factor: e^{-(w+Delta0)s}e^{-Delta0(ell-s)}
    # equals e^{-w*s}e^{-Delta0*ell}.  Comparing rational exponents is exact.
    for weight in (F(-2), F(0), F(3, 4)):
        for pole in (F(-3, 2), F(0), F(5, 3)):
            for shift in shifts:
                for ell in log_references:
                    after = -(weight + pole) * shift - pole * (ell - shift)
                    before = -weight * shift - pole * ell
                    if after != before:
                        raise AssertionError(("scalar compensator", weight, pole, shift, ell))
                    scalar_cases += 1

    # A deliberately false rule must fail: at a neutral double pole, the
    # residue alone is not invariant if the double-pole coefficient is nonzero.
    original = [F(3), F(2)]
    shift = F(5, 7)
    transformed = apply(jet_matrix(2, shift), original)
    if transformed[0] == original[0]:
        raise AssertionError("Negative control did not reject residue invariance")
    if transformed != [F(11, 7), F(2)]:
        raise AssertionError("Unexpected double-pole transformation")
    ell = F(1, 3)
    normalized_before = apply(jet_matrix(2, ell), original)
    normalized_after = apply(jet_matrix(2, ell - shift), transformed)
    if normalized_before != normalized_after:
        raise AssertionError("Explicit compensated double-pole example failed")

    return {
        "status": "passed",
        "arithmetic": "exact rational, Python standard library",
        "max_pole_order": max_order,
        "convention": {
            "raw_mellin": "M(Delta)=integral omega^Delta a(omega) d(omega)/omega",
            "unit_change": "omega'=exp(-sigma)omega, A'=A+d sigma",
            "raw_weight": "w+Delta",
            "normalized_mellin": "Mhat=mu^(-Delta)M; mu'=exp(-sigma)mu",
            "normalized_weight": "w",
        },
        "checks": {
            "group_law_matrix_cases": group_cases,
            "inverse_law_matrix_cases": inverse_cases,
            "independent_laurent_product_vector_cases": laurent_cases,
            "normalized_reference_nilpotent_cancellation_cases": compensator_cases,
            "normalized_reference_scalar_exponent_cases": scalar_cases,
        },
        "negative_control": {
            "incorrect_claim": "The residue of a neutral double pole is invariant",
            "outcome": "rejected",
            "sigma": rational(shift),
            "original_coefficients_a_minus_1_a_minus_2": list(map(rational, original)),
            "transformed_coefficients_a_minus_1_a_minus_2": list(map(rational, transformed)),
            "log_mu_before": rational(ell),
            "normalized_jet_before": list(map(rational, normalized_before)),
            "normalized_jet_after": list(map(rational, normalized_after)),
        },
        "scope": "Finite exact replay of stated scale-jet identities; no physical or novelty claim",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-order", type=int, default=12)
    parser.add_argument(
        "--output", type=Path,
        default=Path(__file__).resolve().parents[1] / "reports" / "scale-jets.json",
    )
    arguments = parser.parse_args()
    report = verify(arguments.max_order)
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(f"Scale jets passed exactly through order {arguments.max_order}; "
          "false residue-invariance rule rejected.")
    print(f"Report: {arguments.output}")


if __name__ == "__main__":
    main()
