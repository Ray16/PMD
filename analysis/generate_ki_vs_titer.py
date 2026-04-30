#!/usr/bin/env python3
"""Generate Ki vs Titer scatter plot for the 5 variants with measured Ki."""

import numpy as np
from scipy import stats

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT_DIR = '/nfs/lambda_stor_01/homes/rzhu/PMD/analysis'

DATA = {
    'V230E': {'ki': 10,  'titer': 278},
    'WT':    {'ki': 18,  'titer': 475},
    'GKQ':   {'ki': 11,  'titer': 8},
    'HKQ':   {'ki': 80,  'titer': 1079},
    'R74G':  {'ki': 110, 'titer': 975},
}

COLORS = {
    'WT':    '#2ca02c',
    'R74G':  '#1f77b4',
    'GKQ':   '#d62728',
    'HKQ':   '#9467bd',
    'V230E': '#8c564b',
}

MARKERS = {
    'WT': 'o', 'R74G': '^', 'GKQ': 'D', 'HKQ': 'v', 'V230E': 'P',
}

fig, ax = plt.subplots(figsize=(7, 5))

kis, titers, labels = [], [], []
for name, d in DATA.items():
    ki, titer = d['ki'], d['titer']
    kis.append(ki)
    titers.append(titer)
    labels.append(name)
    ax.scatter(ki, titer, s=80, c=COLORS[name], marker=MARKERS[name],
               edgecolors='k', linewidth=1.2, zorder=5, label=name)

offsets = {
    'V230E': (10, -18),
    'WT':    (10, -18),
    'GKQ':   (10, 8),
    'HKQ':   (-12, 10),
    'R74G':  (-15, -20),
}

for name, d in DATA.items():
    ki, titer = d['ki'], d['titer']
    ox, oy = offsets[name]
    ax.annotate(f'{name}\n(Ki={ki}, titer={titer})',
                xy=(ki, titer), xytext=(ox, oy),
                textcoords='offset points',
                fontsize=9, fontweight='bold',
                ha='left' if ox > 0 else 'right',
                va='bottom' if oy > 0 else 'top')

slope, intercept, r, p, se = stats.linregress(kis, titers)
x_fit = np.linspace(0, 120, 100)
y_fit = slope * x_fit + intercept
ax.plot(x_fit, y_fit, '--', color='gray', alpha=0.5, linewidth=1.5, zorder=2)
ax.text(0.05, 0.95, f'r = {r:.2f}, p = {p:.3f}',
        transform=ax.transAxes, fontsize=10, va='top',
        bbox=dict(boxstyle='round,pad=0.3', facecolor='wheat', alpha=0.5))

ax.set_xlabel('Ki (mM)', fontsize=12)
ax.set_ylabel('Isopentenol Titer (mg/L)', fontsize=12)
ax.set_title('Substrate Inhibition Constant vs In Vivo Titer', fontsize=13, fontweight='bold')
ax.tick_params(labelsize=10)
ax.set_xlim(-5, 125)
ax.set_ylim(-80, 1200)

ax.legend(fontsize=9, loc='lower right', framealpha=0.8)

plt.tight_layout()
out_path = f'{OUT_DIR}/ki_vs_titer.png'
fig.savefig(out_path, dpi=200, bbox_inches='tight')
print(f'Saved: {out_path}')
plt.close()
