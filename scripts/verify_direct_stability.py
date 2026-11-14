#!/usr/bin/env python3
"""Independent numerical controls for the analytical stability theorems.

No finite calculation proves distributional convergence or a minimax theorem.
The complete arguments, observation norms, and class assumptions are in the paper.
"""
from __future__ import annotations

from datetime import datetime, timezone
import json
from pathlib import Path

import mpmath as mp

ROOT = Path(__file__).resolve().parents[1]
mp.mp.dps = 45


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def text(x):
    return mp.nstr(x, 18)


def profile(t, n):
    return mp.exp(-n * mp.log1p(mp.exp(t)))


def boundary_norm_squared(n, q):
    pivot = -mp.log(n)
    points = [-mp.inf, pivot - 8, pivot, pivot + 8, 0, 8, mp.inf]
    points = sorted(set(points))
    return 2 * mp.pi * mp.quad(
        lambda t: (1 + t*t)**(-q) * profile(t, n)**2, points)


def main():
    boundary = []
    for q in (mp.mpf('0.75'), mp.mpf(1), mp.mpf('1.5')):
        previous = mp.inf
        for n in (2, 16, 256, 65536):
            norm2 = boundary_norm_squared(n, q)
            require(0 < norm2 < previous, 'boundary norm must decrease with N')
            previous = norm2
            asymptotic = 2 * mp.pi * mp.log(n)**(1-2*q) / (2*q-1)
            boundary.append({'N': n, 'q': text(q), 'norm_squared': text(norm2),
                             'ratio_to_asymptotic_diagnostic': text(norm2/asymptotic)})
    # A separate compact-interval substitution checks normalization and quadrature.
    n = 16
    compact = 2 * mp.pi * mp.quad(
        lambda u: profile(mp.tan(u), n)**2,
        [-mp.pi/2, -mp.atan(mp.log(n)), 0, mp.pi/2])
    require(abs(compact - boundary_norm_squared(n, mp.mpf(1))) < mp.mpf('1e-30'),
            'independent H^-1 quadratures disagree')
    # Fixed Gaussian detector pairings also decrease; this is a spot check only.
    detectors = []
    previous = mp.inf
    for n in (2, 8, 32, 128):
        pairing = mp.sqrt(2*mp.pi) * mp.quad(
            lambda t: mp.exp(-t*t/2) * profile(t, n), [-12, -8, -4, 0, 4, 8, 12])
        require(0 < pairing < previous, 'fixed detector pairing did not decrease')
        previous = pairing
        detectors.append({'N': n, 'gaussian_pairing': text(pairing)})

    lower = []
    for a in (mp.mpf('0.5'), mp.mpf(1), mp.mpf(2)):
        for alpha in (mp.mpf('0.5'), mp.mpf(1), mp.mpf(2)):
            M, eta = mp.mpf('1.7'), mp.mpf('0.0001')
            theta = alpha/(a+alpha)
            c = (2*alpha*eta*M**(a/alpha)/mp.gamma(a/alpha))**theta
            exact_norm = mp.gamma(a/alpha)/alpha*c**(1+a/alpha)*M**(-a/alpha)
            require(abs(exact_norm/(2*eta)-1) < mp.mpf('1e-40'),
                    'two-point observation balls do not meet at noise eta')
            for x in (mp.mpf('0.001'), mp.mpf('0.1'), mp.mpf(1), mp.mpf(10)):
                f = c*mp.exp(-(M/c)*x**alpha)
                require(abs(f-c) <= M*x**alpha, 'lower profile violates modulus')
            z = mp.mpc(a, mp.mpf('0.7'))
            analytic = c/alpha*(c/M)**(z/alpha)*mp.gamma(z/alpha)
            # Independent energy integral after u=(M/c)x^alpha.
            integral = c/alpha*(c/M)**(z/alpha)*mp.quad(
                lambda u: u**(z/alpha-1)*mp.exp(-u), [0, mp.mpf('0.1'), 1, 8, mp.inf])
            require(abs(integral-analytic) < mp.mpf('1e-12'),
                    'lower profile Mellin integral disagrees with gamma formula')
            lower.append({'a': text(a), 'alpha': text(alpha), 'residue_separation': text(c),
                          'line_norm': text(exact_norm), 'minimax_lower_bound': text(c/2)})

    gaussian = []
    a, c, lam, eta = mp.mpf(1), mp.mpf('1.3'), mp.mpf('0.7'), mp.mpf('0.0001')
    B = mp.exp(a*a/2)/mp.sqrt(2*mp.pi)
    for h in (mp.mpf('0.05'), mp.mpf('0.2')):
        def kernel(nu):
            z = mp.mpc(a, nu)
            return h**(-z)*mp.exp(z*z/2)/(2*mp.pi)

        def data(nu):
            z = mp.mpc(a, nu)
            return c*lam**(-z)*mp.gamma(z)

        points = list(range(-14, 15, 2))
        mellin_average = mp.quad(lambda nu: data(nu)*kernel(nu), points)
        energy_average = mp.quad(lambda t: mp.exp(-t*t/2)*c*mp.exp(-lam*h*mp.exp(t)),
                                 [-14, -8, -4, 0, 4, 8, 14])/mp.sqrt(2*mp.pi)
        require(abs(mellin_average-energy_average) < mp.mpf('1e-26'),
                'Gaussian Mellin inversion identity failed')
        noise = lambda nu: eta*mp.exp(-1j*nu*(a-mp.log(h)))
        noise_effect = mp.quad(lambda nu: noise(nu)*kernel(nu), points)
        noise_gain = B*h**(-a)
        require(abs(noise_effect/(eta*noise_gain)-1) < mp.mpf('1e-35'),
                'adversarial noise did not attain the Gaussian operator norm')
        require(abs(mellin_average-c) <= mp.exp(mp.mpf('0.5'))*c*lam*h,
                'bias exceeded the global-modulus bound')
        # Wrong sign of the power in the noise amplification is rejected.
        require(abs(noise_gain-B*h**a) > mp.mpf('0.1'), 'wrong h^a control accepted')
        Omega = mp.mpf(3)
        truncated = mp.quad(lambda nu: data(nu)*kernel(nu), [-Omega, 0, Omega])
        tail_bound = (c/lam)*B*h**(-a)*mp.erfc(Omega/mp.sqrt(2))
        require(abs(truncated-mellin_average) <= tail_bound,
                'finite-band error exceeded its certified analytic upper bound')
        gaussian.append({'h': text(h), 'identity_absolute_error': text(abs(mellin_average-energy_average)),
                         'noise_gain': text(noise_gain), 'bias': text(abs(mellin_average-c)),
                         'bandwidth': text(Omega), 'band_error': text(abs(truncated-mellin_average)),
                         'band_error_upper_bound': text(tail_bound)})

    balance = []
    for a, alpha in ((mp.mpf(1), mp.mpf(1)), (mp.mpf(2), mp.mpf(1)),
                     (mp.mpf('0.5'), mp.mpf('0.5'))):
        A, B = mp.exp(alpha*alpha/2), mp.exp(a*a/2)/mp.sqrt(2*mp.pi)
        M, eta = mp.mpf('1.7'), mp.mpf('0.0001')
        h = (a*B*eta/(alpha*A*M))**(1/(a+alpha))
        bound = lambda x: A*M*x**alpha+B*eta*x**(-a)
        require(bound(h) < bound(h/2) and bound(h) < bound(2*h), 'balance control failed')
        theta = alpha/(a+alpha)
        C = A**(a/(a+alpha))*B**theta*((a/alpha)**theta+(alpha/a)**(a/(a+alpha)))
        require(abs(bound(h)/(C*M**(a/(a+alpha))*eta**theta)-1) < mp.mpf('1e-40'),
                'stated upper constant is not the optimized bound')
        balance.append({'a': text(a), 'alpha': text(alpha), 'h_star': text(h), 'upper_bound': text(bound(h))})

    report = {'schema': 'direct-mellin-stability-controls-v1',
              'generated_utc': datetime.now(timezone.utc).isoformat(),
              'working_decimal_digits': mp.mp.dps, 'formal_proof': False,
              'boundary_norm_quadrature': boundary, 'fixed_detector_spot_checks': detectors,
              'two_point_lower_bounds': lower, 'gaussian_and_band_controls': gaussian,
              'bias_noise_balance': balance, 'passed': True,
              'limits': ['Numerical controls do not prove the infinite analytic statements.',
                         'Continuum frequency data; discrete sampling error is not bounded.',
                         'No physical scattering realization or full celestial duality is established.']}
    output = ROOT/'reports'/'direct-stability.json'
    output.write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')
    print('PASS: boundary norms, Gaussian inversion, adversarial noise, lower bounds, finite band and balance')


if __name__ == '__main__':
    main()
