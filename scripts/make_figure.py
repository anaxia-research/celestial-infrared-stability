"""Plot the exact a=1 Mellin-line norm and its fixed soft residue."""

from pathlib import Path

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


def main():
    root = Path(__file__).resolve().parents[1]
    orders = np.arange(2, 201)
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 11})
    figure, axis = plt.subplots(figsize=(8.2, 4.2), layout='constrained')
    figure.patch.set_facecolor('#ffffff')
    axis.set_facecolor('#f5fafb')
    axis.plot(orders, 1 / (orders - 1), color='#007b83', linewidth=2.5,
              label='Exact Mellin-line supremum')
    axis.axhline(1, color='#c46b22', linewidth=2, linestyle='--', label='Soft residue')
    axis.set_yscale('log')
    axis.set_xlabel('Profile order N')
    axis.set_ylabel('Magnitude, with c = R = a = 1')
    axis.grid(alpha=.2, which='both')
    axis.legend(frameon=False, loc='upper right')
    axis.spines[['top', 'right']].set_visible(False)
    figure.savefig(root / 'assets' / 'residue-instability.png', dpi=190)
    plt.close(figure)
    print('Saved the exact instability figure.')


if __name__ == '__main__':
    main()
