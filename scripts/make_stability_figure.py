#!/usr/bin/env python3
"""Plot analytical recovery bounds and separately labelled boundary quadratures."""
import json
from pathlib import Path
import math

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
report = json.loads((ROOT/'reports'/'direct-stability.json').read_text())
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10,
                     'axes.spines.top': False, 'axes.spines.right': False})
fig, axes = plt.subplots(1, 2, figsize=(11.6, 4.1), constrained_layout=True)
fig.patch.set_facecolor('#f5fbf8')
for ax in axes:
    ax.set_facecolor('#ffffff')
    ax.grid(alpha=.15)

rows = [r for r in report['boundary_norm_quadrature'] if float(r['q']) == 1]
axes[0].semilogx([r['N'] for r in rows], [math.sqrt(float(r['norm_squared'])) for r in rows],
                'o-', color='#087f6d', label='Boundary H^-1 norm (quadrature)')
axes[0].axhline(1, color='#c2923a', linestyle='--', label='Soft residue (exact)')
axes[0].set(xlabel='Profile index N', ylabel='Norm / residue magnitude',
            title='A fixed residue with vanishing boundary data')
axes[0].legend(fontsize=8)

etas = [10**(-8+i/20) for i in range(141)]
lo = [math.sqrt(e/2) for e in etas]
hi = [2*math.sqrt(math.e)/(2*math.pi)**.25*math.sqrt(e) for e in etas]
axes[1].loglog(etas, hi, color='#087f6d', label='Proved estimator upper bound')
axes[1].loglog(etas, lo, '--', color='#c2923a', label='Proved minimax lower bound')
axes[1].fill_between(etas, lo, hi, color='#087f6d', alpha=.08)
axes[1].set(xlabel='Uniform Mellin-line noise eta', ylabel='Absolute residue error',
            title='Regularity restores the optimal square-root rate')
axes[1].text(.05, .90, 'a = alpha = M = 1; continuum-line data',
             transform=axes[1].transAxes, fontsize=8)
axes[1].legend(loc='lower right', fontsize=8)
destination = ROOT/'assets'/'stability-tradeoff.png'
fig.savefig(destination, dpi=190, facecolor=fig.get_facecolor())
print(destination)
