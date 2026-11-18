"""N2 / cc-pVDZ (all electrons): SUCCSD / SGCCSD with coherent (T) or +(T), against FCI.
Usage: python plot_n2.py   (figures are written next to this script)
"""
from pathlib import Path
import csv

import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import CubicSpline

ROOT = Path(__file__).resolve().parent
INPUT = ROOT / 'N2.csv'   # the N2 sheet of scan_data.xlsx
OUTDIR = ROOT   # figures are written next to this script
HA2EV = 27.211386
E_INF = {
    'SUCCSD': -108.9586862, 'SUCCSD_cohT': -108.9599861, 'SUCCSD_plusT': -108.9599570,   # redo, 20 A
    'SGCCSD': -108.9583377, 'SGCCSD_cohT': -108.9602950, 'SGCCSD_plusT': -108.9597573,
    'FCI': 2 * -54.4801150544,
}
RED, BLUE, BLACK = '#e41a1c', '#377eb8', 'black'
ORANGE, GREEN = '#ff7f00', '#4daf4a'   # + (T): warm/cool partner of the family colour
FIGS = {
    'N2_DZ_T.pdf': [
        ('SUCCSD', 'SUCCSD', dict(color=RED, linestyle='none', marker='o')),
        ('SUCCSD+coh (T)', 'SUCCSD_cohT', dict(color=RED, linestyle='-', marker='s')),
        ('SUCCSD+(T)', 'SUCCSD_plusT', dict(color=ORANGE, linestyle='-', marker='D', error_only=True)),
        ('SGCCSD', 'SGCCSD', dict(color=BLUE, linestyle='none', marker='o')),
        ('SGCCSD+coh (T)', 'SGCCSD_cohT', dict(color=BLUE, linestyle='-', marker='s')),
        ('SGCCSD+(T)', 'SGCCSD_plusT', dict(color=GREEN, linestyle='-', marker='D', error_only=True)),
    ],
    'N2_DZ_plusT.pdf': [
        ('SUCCSD', 'SUCCSD', dict(color=RED, linestyle='none', marker='o')),
        ('SUCCSD+(T)', 'SUCCSD_plusT', dict(color=ORANGE, linestyle='-', marker='s')),
        ('SGCCSD', 'SGCCSD', dict(color=BLUE, linestyle='none', marker='o')),
        ('SGCCSD+(T)', 'SGCCSD_plusT', dict(color=GREEN, linestyle='-', marker='s')),
    ],
}

import yeplot  # noqa: F401
plt.style.use(['yeplot.yeplot', 'yeplot.onecolumn', 'yeplot.times'])


# column headers of the N2 sheet -> names used below
COLUMNS = {'R (bohr)': 'R_bohr', 'FCI (Eh)': 'FCI', 'SUHF (Eh)': 'SUHF', 'SGHF (Eh)': 'SGHF',
           'SUCCSD (Eh)': 'SUCCSD', 'SUCCSD + coherent (T) (Eh)': 'SUCCSD_cohT', 'SUCCSD + (T) (Eh)': 'SUCCSD_plusT',
           'SGCCSD (Eh)': 'SGCCSD', 'SGCCSD + coherent (T) (Eh)': 'SGCCSD_cohT', 'SGCCSD + (T) (Eh)': 'SGCCSD_plusT'}
with INPUT.open(newline='') as handle:
    rows = [{COLUMNS[key]: value for key, value in row.items()} for row in csv.DictReader(handle)]
r = np.array([float(row['R_bohr']) for row in rows])
rf = np.linspace(r.min(), r.max(), 400)
col = lambda c: np.array([float(row[c]) for row in rows])

fci = col('FCI')

# SUHF / SGHF energies on the same R grid
phf_r = r
PHF_CURVES = [   # label, energies, colour
    ('SUHF', col('SUHF'), RED),
    ('SGHF', col('SGHF'), BLUE),
]
err_e = lambda c: (col(c) - fci) * 1e3
err_b = lambda c: ((col(c) - E_INF[c]) - (fci - E_INF['FCI'])) * 1e3

def common_ylim(f):   # same y range in both figures (splines included)
    cs = {c for curves in FIGS.values() for _, c, _ in curves}
    v = np.concatenate([CubicSpline(r, f(c))(rf) for c in cs]); lo, hi = v.min(), v.max()
    return lo - 0.05 * (hi - lo), hi + 0.05 * (hi - lo)
YLIM_E, YLIM_B = common_ylim(err_e), common_ylim(err_b)

for name, curves in FIGS.items():
    fig, (ax_e, ax_b) = plt.subplots(2, 1, sharex=True, figsize=(3.3, 4.8))
    ax_e.plot(rf, CubicSpline(r, fci)(rf), color='black', linewidth=1.2, label='CDFCI')
    phf_fine = np.linspace(phf_r.min(), phf_r.max(), 400)
    for phf_label, phf_energy, phf_colour in PHF_CURVES:   # spline line + small markers at the PHF points
        ax_e.plot(phf_fine, CubicSpline(phf_r, phf_energy)(phf_fine), color=phf_colour, linewidth=1.0)
        ax_e.plot(phf_r, phf_energy, linestyle='none', marker='o', markersize=2.2, color=phf_colour)
        ax_e.plot([], [], color=phf_colour, linewidth=1.0, marker='o', markersize=2.2, label=phf_label)
    for label, column, style in curves:
        line = {k: v for k, v in style.items() if k != 'marker'}
        dots = style['linestyle'] == 'none'   # CCSD: open circles only, no line
        panels = ((ax_b, err_e(column)),) if style.get('error_only') else ((ax_e, col(column)), (ax_b, err_e(column)))
        for ax, y in panels:   # error_only: + (T) shown in the error panel only, legend entry in the top panel
            if dots:
                ax.plot(r, y, linestyle='none', marker='o', markersize=3.0, markerfacecolor='none',
                        markeredgewidth=0.8, color=style['color'])
            else:
                ax.plot(r, y, linestyle='none', marker=style['marker'], markersize=3.0, color=style['color'])
        if dots:
            ax_e.plot([], [], linestyle='none', marker='o', markersize=3.5, markerfacecolor='none', markeredgewidth=0.8,
                      color=style['color'], label=label)
        else:
            ax_e.plot([], [], linestyle='none', marker=style['marker'], markersize=3.5, color=style['color'], label=label)
    ax_b.axhline(0.0, color='0.4', linewidth=0.6)
    ax_e.set_ylabel(r'$E\;(E_\mathrm{h})$')
    ax_b.set_ylabel(r'$E - E_\mathrm{CDFCI}\;(\mathrm{m}E_\mathrm{h})$')
    ax_e.set_ylim(-109.32, -108.74)
    ax_b.set_ylim(-10.0, YLIM_E[1] + 6.0)   # room for the CC marker legend
    ax_b.set_xlabel(r'$R\;(\mathrm{bohr})$')
    # panel (a) legend: the curves (CDFCI, SUHF, SGHF); panel (b) legend: the CC markers
    handles, labels = ax_e.get_legend_handles_labels()
    is_line = [('CCSD' not in lab) for lab in labels]
    ax_e.legend([h for h, k in zip(handles, is_line) if k], [l for l, k in zip(labels, is_line) if k], ncol=1, loc='lower right',
                fontsize='x-small', handlelength=1.8, handletextpad=0.4, borderaxespad=0.3, labelspacing=0.25)
    ax_b.legend([h for h, k in zip(handles, is_line) if not k], [l for l, k in zip(labels, is_line) if not k], ncol=2, loc='upper right',
                fontsize='x-small', handlelength=1.2, handletextpad=0.3, columnspacing=0.6, borderaxespad=0.3, labelspacing=0.25)
    for ax in (ax_e, ax_b):
        ax.grid(True, linestyle='--', alpha=0.5)
    ax_e.text(0.03, 0.04, '(a)', transform=ax_e.transAxes, ha='left', va='bottom')   # bottom left: PHF curves fill the top left
    ax_b.text(0.03, 0.95, '(b)', transform=ax_b.transAxes, ha='left', va='top')
    name = name[:-4] + '_PHF.pdf'
    fig.savefig(OUTDIR / name, format='pdf')
    plt.close(fig)
    print(f'wrote {OUTDIR / name}')
