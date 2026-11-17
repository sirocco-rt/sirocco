#!/usr/bin/env python
'''
Plot the spectrum in a wind cell: the fine histogram of J_nu, the band model that
Sirocco fitted to it, the band boundaries and the photoionization edges.


Command line usage:
    band_plot.py [options] run [run2 ...]

    run is a run directory or a FITS file written by windsave2fits (root_cellspec.fits).
    For a run directory, root_cellspec.fits must exist: make it with
        windsave2fits -thresh root     (if the run used band model cell_spectrum or -pi_rate_table)
        windsave2fits root             (otherwise)

Options:
    -h              print this documentation
    -cell N         the plasma cell to plot (default 0); several as -cell 0,10,40
    -labels a,b     labels for the runs (default: directory or file names)
    -emin E         lowest photon energy shown, eV (default 10)
    -emax E         highest photon energy shown, eV (default 300)
    -ion "O III"    add a lower panel with the photoionization-rate integrand of this ion,
                    4 pi sigma J_nu / (h nu), from the histogram and from the band model,
                    and print the two rates
    -fmin f         mark the edges of ions with at least this fraction of their element
                    (default 1e-3)
    -nuJ            plot nu J_nu instead of J_nu
    -o file         output file (default band_plot.png)

Description:
    Each column is one run (for example the same model with two band models); each row a cell.
    - grey steps: the fine cell histogram of J_nu (what band model cell_spectrum uses)
    - coloured lines: the band model, a power law or an exponential in each band, as used
      for the photoionization rates (pi_rates.c)
    - dotted vertical lines: the band boundaries
    - labelled ticks at the top: photoionization edges of the ions present in the cell
    With -ion, the lower panel shows where the rate comes from, and the ratio of the band-model
    rate to the histogram rate (as written by sirocco -pi_rate_table).

Primary routines:
    read_cellspec   reads a windsave2fits file
    band_model_jnu  evaluates the band model of one cell on a frequency grid
    steer           processes the command line

Notes:
    Written 2026-11-15 for the banding comparison on the common test wind (Sirocco-2611)
    and for PR #1250.  The band model is evaluated over each whole band; Sirocco limits it to
    the frequencies where photons were seen (fmin_mod, fmax_mod), which windsave2fits does not
    save, so the curves can extend slightly beyond what Sirocco uses.
'''

import sys
import os
from glob import glob
import numpy as np
from astropy.io import fits
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

H = 6.62607015e-27
K = 1.380649e-16
HEV = 4.135667696e-15
ELEMENT = {1: 'H', 2: 'He', 6: 'C', 7: 'N', 8: 'O', 10: 'Ne', 12: 'Mg', 14: 'Si', 16: 'S', 18: 'Ar', 20: 'Ca', 26: 'Fe'}
ROMAN = ['I', 'II', 'III', 'IV', 'V', 'VI', 'VII', 'VIII', 'IX', 'X', 'XI', 'XII', 'XIII', 'XIV', 'XV', 'XVI',
         'XVII', 'XVIII', 'XIX', 'XX', 'XXI', 'XXII', 'XXIII', 'XXIV', 'XXV', 'XXVI', 'XXVII']
SPEC_MOD_PL = 1
SPEC_MOD_EXP = 2
BAND_COLOURS = ['#1f4e79', '#c0504d']


def find_fits(name):
    if os.path.isfile(name):
        return name
    files = glob(os.path.join(name, '*_cellspec.fits'))
    if len(files) != 1:
        print('Error: %s should contain one *_cellspec.fits file (run windsave2fits), not %d' % (name, len(files)))
        sys.exit(1)
    return files[0]


def read_cellspec(filename):
    h = fits.open(filename)
    d = {
        'jnu': np.array(h['Jnu'].data, dtype=float),
        'edges': np.array(h['nu_edges'].data, dtype=float),
        'nu_model': np.array(h['nu_model'].data, dtype=float),
        'model': h['spec_model'].data,
        'z': np.array(h['ions'].data['z']),
        'istate': np.array(h['ions'].data['istate']),
        'sigma': np.array(h['sigma_outer'].data, dtype=float),
        'density': np.array(h['density'].data, dtype=float),
    }
    d['centre'] = np.sqrt(d['edges'][:-1] * d['edges'][1:])
    h.close()
    return d


def band_model_jnu(d, cell, nu):
    '''the band model of one plasma cell on the frequencies nu (nan where there is no model)'''
    m = d['model'][d['model']['nplasma'] == cell]
    out = np.full(len(nu), np.nan)
    band_of = np.searchsorted(d['nu_model'], nu) - 1
    for row in m:
        b = row['band']
        sel = band_of == b
        if row['spec_mod_type'] == SPEC_MOD_PL:
            out[sel] = 10.0 ** (row['pl_log_w'] + row['pl_alpha'] * np.log10(nu[sel]))
        elif row['spec_mod_type'] == SPEC_MOD_EXP:
            out[sel] = row['exp_w'] * np.exp(-H * nu[sel] / (K * row['exp_temp']))
    return out


def ion_edges(d, cell, fmin):
    '''(threshold in eV, name) of the ions present in the cell'''
    out = []
    dens = d['density'][cell]
    for k in range(len(d['z'])):
        z = d['z'][k]
        total = np.sum(dens[d['z'] == z])
        if total <= 0 or dens[k] / total < fmin:
            continue
        nz = np.where(d['sigma'][k] > 0)[0]
        if len(nz) == 0:
            continue
        out.append((d['edges'][nz[0]] * HEV, '%s %s' % (ELEMENT.get(z, 'Z%d' % z), ROMAN[d['istate'][k] - 1])))
    return sorted(out)


def find_ion(d, ion):
    words = ion.split()
    for k in range(len(d['z'])):
        if ELEMENT.get(d['z'][k]) == words[0] and ROMAN[d['istate'][k] - 1] == words[1]:
            return k
    return None


def draw_cell(ax, d, cell, emin, emax, fmin, nuj, colour):
    nu = d['centre']
    e = nu * HEV
    jnu = d['jnu'][cell]
    fine = np.logspace(np.log10(emin / HEV), np.log10(emax / HEV), 3000)
    model = band_model_jnu(d, cell, fine)
    yh = jnu
    ym = model
    if nuj:
        yh = jnu * nu
        ym = model * fine
    ax.step(d['edges'][:-1] * HEV, yh, where='post', color='0.55', lw=1.0, label='cell histogram')
    ax.plot(fine * HEV, ym, color=colour, lw=2.0, label='band model')
    for b in d['nu_model']:
        if emin < b * HEV < emax:
            ax.axvline(b * HEV, color=colour, ls=':', lw=0.8)
    sel = (e > emin) & (e < emax) & (yh > 0)
    if sel.any():
        ax.set_ylim(np.min(yh[sel]) / 3, np.max(yh[sel]) * 3)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(emin, emax)
    for energy, name in ion_edges(d, cell, fmin):
        if emin < energy < emax:
            ax.axvline(energy, color='k', lw=0.6, alpha=0.5)
            ax.text(energy, 1.0, ' ' + name, transform=ax.get_xaxis_transform(), rotation=90, va='top', ha='left', fontsize=7)
    ax.grid(alpha=0.25, which='both')


def draw_integrand(ax, d, cell, k, emin, emax, colour, ion):
    nu = d['centre']
    width = np.diff(d['edges'])
    sigma = d['sigma'][k]
    jnu = d['jnu'][cell]
    integrand_hist = 4 * np.pi * sigma * jnu / (H * nu)
    rate_hist = np.sum(integrand_hist * width)
    model = band_model_jnu(d, cell, nu)
    model = np.where(np.isfinite(model), model, 0.0)
    integrand_model = 4 * np.pi * sigma * model / (H * nu)
    rate_model = np.sum(integrand_model * width)
    ax.step(d['edges'][:-1] * HEV, integrand_hist * nu, where='post', color='0.55', lw=1.0, label='histogram')
    ax.plot(nu * HEV, integrand_model * nu, color=colour, lw=1.6, label='band model')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(emin, emax)
    top = max(np.max(integrand_hist * nu), np.max(integrand_model * nu))
    if top > 0:
        ax.set_ylim(top / 1e4, top * 3)
    ratio = np.nan
    if rate_hist > 0:
        ratio = rate_model / rate_hist
    ax.set_ylabel('%s rate per ln nu' % ion, fontsize=9)
    ax.text(0.98, 0.92, 'band / hist = %.2f (%+.2f dex)' % (ratio, np.log10(ratio)), transform=ax.transAxes,
            ha='right', va='top', fontsize=9)
    ax.grid(alpha=0.25, which='both')
    return rate_model, rate_hist


def steer(argv):
    cells = [0]
    labels = []
    emin = 10.0
    emax = 300.0
    ion = ''
    fmin = 1e-3
    nuj = False
    outfile = 'band_plot.png'
    runs = []

    i = 1
    while i < len(argv):
        if argv[i] == '-h':
            print(__doc__)
            return
        elif argv[i] == '-cell':
            i += 1
            cells = [int(c) for c in argv[i].split(',')]
        elif argv[i] == '-labels':
            i += 1
            labels = argv[i].split(',')
        elif argv[i] == '-emin':
            i += 1
            emin = float(argv[i])
        elif argv[i] == '-emax':
            i += 1
            emax = float(argv[i])
        elif argv[i] == '-ion':
            i += 1
            ion = argv[i]
        elif argv[i] == '-fmin':
            i += 1
            fmin = float(argv[i])
        elif argv[i] == '-nuJ':
            nuj = True
        elif argv[i] == '-o':
            i += 1
            outfile = argv[i]
        elif argv[i].startswith('-'):
            print('Error: unknown option %s' % argv[i])
            print(__doc__)
            return
        else:
            runs.append(argv[i])
        i += 1

    if len(runs) == 0:
        print(__doc__)
        return
    if len(labels) != len(runs):
        labels = [os.path.basename(os.path.normpath(r)) for r in runs]

    data = [read_cellspec(find_fits(r)) for r in runs]
    nrow_per_cell = 1
    if ion:
        nrow_per_cell = 2
    ncol = len(runs)
    nrow = nrow_per_cell * len(cells)
    ratios = [3] * nrow
    if ion:
        ratios = [3, 1.6] * len(cells)
    fig, axes = plt.subplots(nrow, ncol, figsize=(6.5 * ncol, sum(ratios) * 1.25 + 0.8), squeeze=False,
                             gridspec_kw={'height_ratios': ratios}, sharex=True)
    for c, cell in enumerate(cells):
        for n, d in enumerate(data):
            colour = BAND_COLOURS[n % len(BAND_COLOURS)]
            ax = axes[c * nrow_per_cell, n]
            draw_cell(ax, d, cell, emin, emax, fmin, nuj, colour)
            ylabel = 'J$_\\nu$'
            if nuj:
                ylabel = '$\\nu$ J$_\\nu$'
            ax.set_ylabel('cell %d: %s' % (cell, ylabel))
            if c == 0:
                ax.set_title(labels[n], fontsize=11)
                ax.legend(fontsize=8, loc='lower left')
            if ion:
                k = find_ion(d, ion)
                if k is None:
                    print('Error: %s is not in %s' % (ion, runs[n]))
                    sys.exit(1)
                rm, rh = draw_integrand(axes[c * nrow_per_cell + 1, n], d, cell, k, emin, emax, colour, ion)
                print('%s cell %d %s: band-model rate %.4e, histogram rate %.4e, ratio %.3f' % (labels[n], cell, ion, rm, rh, rm / rh))
    for ax in axes[-1]:
        ax.set_xlabel('photon energy (eV)')
    fig.tight_layout()
    fig.savefig(outfile, dpi=100)
    print('Wrote %s' % outfile)


if __name__ == '__main__':
    steer(sys.argv)
