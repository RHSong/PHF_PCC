"""C2 / cc-pVDZ (all electrons): C4 x SUHF-based PCCSD with coherent (T) or +(T), against iCI.
Usage: python plot_c2.py   (figures are written next to this script)
"""
import os
from pathlib import Path
import csv
import numpy as np
import matplotlib.pyplot as plt
import yeplot  # noqa: F401  (registers the styles)
from scipy.interpolate import CubicSpline

plt.style.use(['yeplot.yeplot', 'yeplot.onecolumn', 'yeplot.times'])
ROOT = Path(__file__).resolve().parent
RED, BLUE, ORANGE, GREEN = '#e41a1c', '#377eb8', '#ff7f00', '#4daf4a'   # coherent (T): red/blue; + (T): orange/green

# column headers of the C2 sheet -> names used below
COLUMNS = {'R (bohr)': 'R_bohr', '1Sigma_g+ iCI (Eh)': 'iCI_1Sigma_g', '1Delta_g iCI (Eh)': 'iCI_1Delta_g',
           'C4A SUHF1 (pi4) (Eh)': 'SUHF_C4A_pi4', 'C4A SUHF2 (pz+px) (Eh)': 'SUHF_C4A_pzpx', 'C4B SUHF (Eh)': 'SUHF_C4B',
           'C4A (best) PCCSD (Eh)': 'C4A_PCCSD', 'C4A (best) PCCSD + coherent (T) (Eh)': 'C4A_cohT',
           'C4A (best) PCCSD + (T) (Eh)': 'C4A_plusT', 'C4B PCCSD (Eh)': 'C4B_PCCSD',
           'C4B PCCSD + coherent (T) (Eh)': 'C4B_cohT', 'C4B PCCSD + (T) (Eh)': 'C4B_plusT'}
rows = [{COLUMNS.get(key, key): value for key, value in row.items()} for row in csv.DictReader(open(ROOT / 'C2.csv'))]
col = lambda name: np.array([float(row[name]) for row in rows])

# iCI extrapolated energies
fr = col('R_bohr')
fS = col('iCI_1Sigma_g')
fD = col('iCI_1Delta_g')
m = (fr >= 1.95) & (fr <= 4.05)

# CC data (C4(A) = lower-PCCSD branch at each R)
R = col('R_bohr')
S_cc, S_coh, S_plus = col('C4A_PCCSD'), col('C4A_cohT'), col('C4A_plusT')
D_cc, D_coh, D_plus = col('C4B_PCCSD'), col('C4B_cohT'), col('C4B_plusT')

# PHF data: one column per SUHF branch, empty where a branch has no point
phf_rows = rows
phf_r = R
def phf_column(name):
    return np.array([float(row[name]) if row[name] else np.nan for row in phf_rows])
PHF_CURVES = [   # label, energies, colour, linestyle
    (r'$C_4(A)$ SUHF1', phf_column('SUHF_C4A_pi4'), RED, '-'),
    (r'$C_4(A)$ SUHF2', phf_column('SUHF_C4A_pzpx'), RED, ':'),
    (r'$C_4(B)$ SUHF', phf_column('SUHF_C4B'), BLUE, '-'),
]

fS_at = np.array([dict(zip(np.round(fr, 2), fS))[round(x, 2)] for x in R]); fD_at = np.array([dict(zip(np.round(fr, 2), fD))[round(x, 2)] for x in R])
FIGS = {'C2_ccpvdz_T.pdf': ('+ coh (T)', S_coh, D_coh, RED, BLUE), 'C2_ccpvdz_plusT.pdf': ('+ (T)', S_plus, D_plus, ORANGE, GREEN)}
for name, (lab, S_t, D_t, cS, cD) in FIGS.items():
    fig, (ax_e, ax_b) = plt.subplots(2, 1, sharex=True, figsize=(3.3, 4.8))
    tl = lab.replace('(T)', r'$(T)$')
    ax_e.plot(fr[m], np.array(fS)[m], color='black', linestyle='-', label=r'iCI $^1\Sigma_g^+$')
    ax_e.plot(fr[m], np.array(fD)[m], color='black', linestyle='--', label=r'iCI $^1\Delta_g$')
    for phf_label, phf_energy, phf_colour, phf_style in PHF_CURVES:   # lines only, spline through the PHF points
        have = ~np.isnan(phf_energy)
        r_fine = np.linspace(phf_r[have].min(), phf_r[have].max(), 400)
        ax_e.plot(r_fine, CubicSpline(phf_r[have], phf_energy[have])(r_fine), color=phf_colour, linestyle=phf_style,
                  linewidth=1.0)
        ax_e.plot(phf_r[have], phf_energy[have], linestyle='none', marker='o', markersize=2.2, color=phf_colour)
        ax_e.plot([], [], color=phf_colour, linestyle=phf_style, linewidth=1.0, marker='o', markersize=2.2, label=phf_label)
    for ax, (ys, yt, yds, ydt) in ((ax_e, (S_cc, S_t, D_cc, D_t)), (ax_b, ((S_cc - fS_at) * 1e3, (S_t - fS_at) * 1e3, (D_cc - fD_at) * 1e3, (D_t - fD_at) * 1e3))):
        ax.plot(R, ys, linestyle='none', marker='o', mfc='none', color=RED, label=r'$C_4(A)$ SUCCSD')
        ax.plot(R, yt, linestyle='none', marker='s', color=cS, label=r'$C_4(A)$ SUCCSD %s' % tl)
        ax.plot(R, yds, linestyle='none', marker='o', mfc='none', color=BLUE, label=r'$C_4(B)$ SUCCSD')
        ax.plot(R, ydt, linestyle='none', marker='s', color=cD, label=r'$C_4(B)$ SUCCSD %s' % tl)
        ax.grid(linestyle='--', alpha=0.5)
    if name == 'C2_ccpvdz_T.pdf':   # error panel also shows + (T) (plain Fock); legend entries go to the top panel
        for yp, c, lab2 in (((S_plus - fS_at) * 1e3, ORANGE, r'$C_4(A)$ SUCCSD + $(T)$'), ((D_plus - fD_at) * 1e3, GREEN, r'$C_4(B)$ SUCCSD + $(T)$')):
            ax_b.plot(R, yp, linestyle='none', marker='D', markersize=3.0, color=c)
            ax_e.plot([], [], linestyle='none', marker='D', markersize=3.0, color=c, label=lab2)
    ax_b.axhline(0.0, color='0.4', linewidth=0.6)
    ax_e.set_ylabel(r'$E$ ($E_\mathrm{h}$)'); ax_b.set_ylabel(r'$E - E_\mathrm{iCI}$ (m$E_\mathrm{h}$)')
    ax_e.set_ylim(-75.775, -75.30)   # headroom for the legend above the 1Delta_g curve at 2.0 bohr
    # panel (a) legend: the curves (iCI, SUHF); panel (b) legend: the CC markers
    handles, labels = ax_e.get_legend_handles_labels()
    is_line = [('SUCCSD' not in lab) for lab in labels]
    ax_e.legend([h for h, k in zip(handles, is_line) if k], [l for l, k in zip(labels, is_line) if k], loc='upper right', ncol=2,
                fontsize=6.3, columnspacing=0.6, handletextpad=0.3, handlelength=1.8, borderaxespad=0.25, labelspacing=0.2)
    ax_b.legend([h for h, k in zip(handles, is_line) if not k], [l for l, k in zip(labels, is_line) if not k], loc='upper right', ncol=2,
                fontsize=6.3, columnspacing=0.4, handletextpad=0.3, handlelength=1.2, borderaxespad=0.25, labelspacing=0.2)
    ax_b.set_ylim(top=ax_b.get_ylim()[1] + (4.5 if name == 'C2_ccpvdz_T.pdf' else 2.5))   # room for the marker legend
    ax_e.text(0.03, 0.04, '(a)', transform=ax_e.transAxes, va='bottom')
    ax_b.text(0.03, 0.96, '(b)', transform=ax_b.transAxes, va='top')   # top left: the 2.0-bohr points sit at the bottom
    ax_b.set_xlabel(r'$R$ (bohr)'); ax_b.set_xlim(1.95, 4.05)
    fig.savefig(ROOT / (name[:-4] + '_PHF.pdf'))   # figures are written next to this script
    print('saved', name)
