#!/usr/bin/env python3
"""Exact local contact-term algebra and a published TT counterexample replay."""
from datetime import datetime, timezone
from math import comb
import json
from pathlib import Path

import sympy as sp


def main():
    # When a smooth polynomial multiplies delta^(n), its term with j
    # derivatives on the polynomial contributes (-1)^j binom(n,j).
    # Combine this with Leibniz and D(1/z)=pi delta. The coefficient of
    # pi P^(k)(bar z_i) delta^(l-k) must be one for every k.
    rows = {}
    for degree in range(1, 9):
        coefficients = [
            sum(comb(degree + 1, i) * (-1)**(k-i) * comb(degree-i, k-i)
                for i in range(k + 1))
            for k in range(degree + 1)
        ]
        assert coefficients == [1] * (degree + 1)
        rows[str(degree)] = coefficients

    # Karan et al. arXiv:2610.07134v1, eqs. 4.26 and 4.27.
    # Divide by the common scalar K and project to TT relative to e_z.
    transverse = sp.diag(1, 1, 0)
    difference = (sp.Rational(125, 16)-sp.Rational(225, 16))*sp.diag(1, 0, 0)
    projected = transverse*difference*transverse
    tt = projected-sp.trace(projected)*transverse/2
    assert tt == sp.diag(sp.Rational(-25, 8), sp.Rational(25, 8), 0)
    assert sp.trace(tt) == 0 and tt[:, 2].is_zero_matrix
    report = {
        'schema': 'celestial-ward-controls-v1',
        'generated_utc': datetime.now(timezone.utc).isoformat(),
        'passed': True,
        'formal_proof': False,
        'contact_coefficients': rows,
        'published_newtonian_TT_difference_divided_by_K': [
            [str(value) for value in row] for row in tt.tolist()
        ],
        'scope': 'Conditional polynomial contact terms and a known TT projection. No quantum loop or Ward identity is derived.',
    }
    output = Path(__file__).resolve().parents[1]/'reports'/'ward-controls.json'
    output.write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')
    print('PASS: degree 1–8 contact coefficients and published TT projection')


if __name__ == '__main__':
    main()
