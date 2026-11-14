#!/usr/bin/env python3
"""Exact ordered-limit identities and labelled double-Mellin spot checks."""
from datetime import datetime, timezone
import json
from pathlib import Path

import mpmath as mp
import sympy as sp

mp.mp.dps = 50
ROOT = Path(__file__).resolve().parents[1]


def main():
    x, y = sp.symbols('x y', positive=True)
    z, w = sp.symbols('z w')
    amplitude = sp.exp(-x-y)/(y*(x+y))
    first_x = sp.limit(x*amplitude, x, 0, dir='+')
    first_y = sp.limit(y*amplitude, y, 0, dir='+')
    xy = sp.limit(y*first_x, y, 0, dir='+')
    yx = sp.limit(x*first_y, x, 0, dir='+')
    assert xy == 0 and yx == 1
    transform = sp.gamma(z+1)*sp.gamma(w)/(z+w)
    rx = sp.residue(transform, z, 0)
    ry = sp.residue(transform, w, 0)
    assert sp.residue(rx, w, 0) == 0
    assert sp.residue(ry, z, 0) == 1
    assert xy != yx  # Negative control: commuting-target assertion is false.
    for q in (sp.Rational(-1), sp.Rational(0), sp.Rational(1, 3),
              sp.Rational(1, 2), sp.Rational(1), sp.Rational(2)):
        assert max(abs(q), abs(q-1)) >= sp.Rational(1, 2)
    rows = []
    for zz, ww in ((mp.mpf('0.7'), mp.mpf('1.2')),
                   (mp.mpc('1.1', '.3'), mp.mpc('.8', '-.2'))):
        # Independent integrations after the proved u=x+y,v=x/u change of variables.
        radial = mp.quad(lambda u: u**(zz+ww-1)*mp.exp(-u), [0, 1, 5, mp.inf])
        angular = mp.quad(lambda v: v**zz*(1-v)**(ww-1), [0, mp.mpf('.5'), 1])
        value = radial*angular
        exact = mp.gamma(zz+1)*mp.gamma(ww)/(zz+ww)
        error = abs(value-exact)
        assert error < mp.mpf('1e-35')
        rows.append({'z': str(zz), 'w': str(ww),
                     'absolute_quadrature_error': mp.nstr(error, 20)})
    report = {'schema': 'celestial-ordered-limits-v1',
              'generated_utc': datetime.now(timezone.utc).isoformat(),
              'formal_proof': False, 'passed': True,
              'exact_sequential_soft_targets': [str(xy), str(yx)],
              'exact_iterated_residues': ['0', '1'],
              'double_Mellin_quadrature_spot_checks': rows,
              'scope': 'Scalar mixed-divisor diagnostic; not a physical OPE or associativity theorem.'}
    (ROOT/'reports'/'ordered-limits.json').write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')
    print('PASS: ordered soft limits, iterated residues, and double-Mellin controls')


if __name__ == '__main__':
    main()
