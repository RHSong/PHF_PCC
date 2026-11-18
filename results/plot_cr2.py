"""Cr2 / cc-pVDZ-DK (sfX2C, frozen Ne core): SUCCSD / SGCCSD with coherent (T) or +(T), against DMRG of
Larsson et al. (arXiv:2206.10738, Table S2): extrapolated E(inf) and unextrapolated E(Dmax).
"""
from pathlib import Path
import csv

import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import CubicSpline

ROOT = Path(__file__).resolve().parent
INPUT = ROOT / 'Cr2.csv'   # the Cr2 sheet of scan_data.xlsx
OUTDIR = ROOT   # figures are written next to this script
HA2EV = 27.211386
E_ATOMS = 2 * -1049.93257
E_INF = {
    'SUCCSD': -2099.8527828, 'SUCCSD_cohT': -2099.8644543, 'SUCCSD_plusT': -2099.8643718,   # 50 A, redo
    'SGCCSD': -2099.8526987, 'SGCCSD_cohT': -2099.8644510, 'SGCCSD_plusT': -2099.8645151,   # 50 A, redo
    'DMRG': E_ATOMS, 'DMRG_raw': E_ATOMS,
}
import sys
ONLY = sys.argv[1:]   # optional: names of the figures to make

RED, BLUE, BLACK = '#e41a1c', '#377eb8', 'black'
ORANGE, GREEN = '#ff7f00', '#4daf4a'   # + (T): warm/cool partner of the family colour
FIGS = {
    'Cr2_DZ_SG_T.pdf': [
        ('SUCCSD', 'SUCCSD', dict(color=RED, linestyle='--', marker='o')),
        ('SUCCSD+coh (T)', 'SUCCSD_cohT', dict(color=RED, linestyle='-', marker='s')),
        ('SUCCSD+(T)', 'SUCCSD_plusT', dict(color=ORANGE, linestyle='-', marker='D')),
        ('SGCCSD', 'SGCCSD', dict(color=BLUE, linestyle='--', marker='o')),
        ('SGCCSD+coh (T)', 'SGCCSD_cohT', dict(color=BLUE, linestyle='-', marker='s')),
        ('SGCCSD+(T)', 'SGCCSD_plusT', dict(color=GREEN, linestyle='-', marker='D')),
    ],
    'Cr2_DZ_SG_plusT.pdf': [
        ('SUCCSD', 'SUCCSD', dict(color=RED, linestyle='--', marker='o')),
        ('SUCCSD+(T)', 'SUCCSD_plusT', dict(color=ORANGE, linestyle='-', marker='s')),
        ('SGCCSD', 'SGCCSD', dict(color=BLUE, linestyle='--', marker='o')),
        ('SGCCSD+(T)', 'SGCCSD_plusT', dict(color=GREEN, linestyle='-', marker='s')),
    ],
}
DMRG = [
    ('DMRG ($D_\\mathrm{max}$)', 'DMRG_raw', dict(color=BLACK, linestyle='--', marker='v')),
    ('DMRG (extrap.)', 'DMRG', dict(color=BLACK, linestyle='-', marker='v')),
]

import yeplot  # noqa: F401
plt.style.use(['yeplot.yeplot', 'yeplot.onecolumn', 'yeplot.times'])


# column headers of the Cr2 sheet -> names used below (the binding-energy columns are not used)
COLUMNS = {'R (Angstrom)': 'R', 'DMRG E(Dmax) (Eh)': 'DMRG_raw', 'DMRG E(extrap) (Eh)': 'DMRG',
           'SUCCSD (Eh)': 'SUCCSD', 'SUCCSD + coherent (T) (Eh)': 'SUCCSD_cohT', 'SUCCSD + (T) (Eh)': 'SUCCSD_plusT',
           'SGCCSD (Eh)': 'SGCCSD', 'SGCCSD + coherent (T) (Eh)': 'SGCCSD_cohT', 'SGCCSD + (T) (Eh)': 'SGCCSD_plusT'}
with INPUT.open(newline='') as handle:
    rows = [{COLUMNS.get(key, key): value for key, value in row.items()} for row in csv.DictReader(handle)]
rows = [row for row in rows if float(row['R']) < 10]   # the 50 A row is the dissociation limit E(inf) below
r = [float(row['R']) for row in rows]

# common y range of the binding-energy panels of all figures (splines included)
_rf = np.linspace(min(r), max(r), 300)
_cols = {c for curves in FIGS.values() for _, c, _ in curves} | {c for _, c, _ in DMRG}
_vals = [CubicSpline(r, (np.array([float(row[c]) for row in rows]) - E_INF[c]) * HA2EV)(_rf) for c in _cols
         if E_INF[c] is not None]
_lo, _hi = min(v.min() for v in _vals), max(v.max() for v in _vals)
YLIM_B = (_lo - 0.05 * (_hi - _lo), _hi + 0.05 * (_hi - _lo))

for name, curves in FIGS.items():
    if ONLY and name not in ONLY:
        continue
    if any(E_INF[c] is None for _, c, _ in curves):
        print(f'skip {name}: an E(inf) is still missing'); continue
    fig, (ax_e, ax_b) = plt.subplots(2, 1, sharex=True, figsize=(3.3, 4.8))
    rf = np.linspace(min(r), max(r), 300)
    for label, column, style in curves + DMRG:
        e = np.array([float(row[column]) for row in rows])
        line = {k: v for k, v in style.items() if k not in ('marker', 'error_only')}
        panels = ((ax_b, (e - E_INF[column]) * HA2EV),) if style.get('error_only') else ((ax_e, e + 2099.0), (ax_b, (e - E_INF[column]) * HA2EV))
        for ax, y in panels:   # error_only: + (T) in the binding-energy panel only, legend entry in the top panel
            # cubic spline through the computed points (line), markers at the points themselves
            ax.plot(rf, CubicSpline(r, y)(rf), linewidth=1.1, **line)
            ax.plot(r, y, linestyle='none', marker=style['marker'], markersize=2.8, color=style['color'])
        ax_e.plot([], [], markersize=2.8, linewidth=1.1, label=label, **{k: v for k, v in style.items() if k != 'error_only'})   # legend entry (top panel): line + marker
    ax_e.set_ylabel(r'$E + 2099\;(E_\mathrm{h})$')
    ax_b.set_ylabel(r'Binding energy (eV)')
    ax_b.set_ylim(*YLIM_B)
    ax_b.set_xlabel(r'$R\;(\mathrm{\AA})$')
    ax_e.set_ylim(-0.93, -0.86)   # cut the steep repulsive wall; compact legend in the free lower-right corner
    for ax, loc in ((ax_e, 'lower right'), (ax_b, 'lower right')):
        h, l = ax_e.get_legend_handles_labels()
        if ax is not ax_e: continue
        ax.legend(h, l, ncol=2, loc=loc, fontsize=6.0, columnspacing=0.5, handlelength=1.0, handletextpad=0.3, borderaxespad=0.15, borderpad=0.25, labelspacing=0.35)
    for ax, tag in ((ax_e, '(a)'), (ax_b, '(b)')):
        ax.grid(True, linestyle='--', alpha=0.5)
        # (a) bottom left; (b) top RIGHT (the steep 1.6 A points fill the top-left corner of the binding panel)
        ax.text(0.03 if ax is ax_e else 0.97, 0.05 if ax is ax_e else 0.95, tag, transform=ax.transAxes,
                ha='left' if ax is ax_e else 'right', va='bottom' if ax is ax_e else 'top')
    fig.savefig(OUTDIR / name, format='pdf')
    plt.close(fig)
    print(f'wrote {OUTDIR / name}')
