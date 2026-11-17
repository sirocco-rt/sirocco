#!/usr/bin/env python
'''
Monitor how the ionization changes from cycle to cycle, alongside Sirocco's own
convergence tests (prototype for run_check.py).


Command line usage:
    ion_monitor.py [options] run [run2 ...]

    run is a run directory or a root name; the per-cycle tables in diag_<root>/
    (sirocco -d, then windsave2table) are used: root.NN.master.txt, .converge.txt,
    .heat.txt, .spec.txt and .<element>.frac.txt.

Options:
    -h              print this documentation
    -tol x          tolerance for the ion changes (default 0.1)
    -nmin n         a cell has information on an element if, had nothing been absorbed, at least
                    n photons would have crossed it in each band holding the threshold of an ion
                    with at least 10% of the element (default 10)
    -cycle N        the cycle for the maps (default the last)
    -labels a;b     labels for the runs, separated by semicolons
    -o root         root for the output files (default ion_monitor)

Description:
    For each run, two figures.

    <root>_cycles.png  one panel per run, fraction of wind cells against ionization cycle:
        solid: Sirocco's tests, as in run_check.py (t_r and t_e change less than
               2 epsilon = 10%, heating and cooling within epsilon = 5%; Converged = all
               three, from the converge column)
        dashed: the ion changes delta_X < tol for X = H, He, N and the abundance-weighted
               C, N, O ("heavy"), where for one element
                   delta = sqrt( sum over ions 2 (f_n - f_(n-1))^2 / (f_n + f_(n-1)) )
               with f the ion fractions in this cycle and the previous one
        dotted: the fraction of cells that pass delta < tol without information (a "vacuous"
               pass), for He and for C/N/O: the cell would have received fewer than nmin photons,
               even if nothing had been absorbed, in the band holding the threshold of an ion with
               at least 10% of its element (thresholds inside the photon range only). The expected
               number is the photons through the cell times the band's share in the least-absorbed
               well-sampled cells. Few photons because the gas absorbed them is information (J is
               really low); few because few were generated in that band is not.

    <root>_maps.png  one row per run, maps for one cycle:
        tests failed (Sirocco's converge column); delta He and delta heavy; and for each a
        three-way map: fails / passes with information / passes without information

    The ion changes are monitored, not added to Sirocco's convergence: their tolerance
    and whether they belong in the criterion are open questions.

Notes:
    Written 2026-11-15 as a prototype for run_check.py (Sirocco-2611). The thresholds are
    those of the atomic data used in the common-wind tests (hhecno_simple), for H, He, C,
    N and O only.
'''

import sys
import os
from glob import glob
import copy
import numpy as np
from astropy.io import ascii
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from mpl_toolkits.axes_grid1 import make_axes_locatable

HEV = 4.135667696e-15
EPSILON = 0.05
ELEMENTS = ['H', 'He', 'C', 'N', 'O']
ABUN = {'C': 10 ** (8.56 - 12), 'N': 10 ** (8.05 - 12), 'O': 10 ** (8.93 - 12)}
# ground-state thresholds (eV) of ion istate = 1, 2, ... (hhecno_simple; the highest from standard values)
THRESH = {'H': [13.60], 'He': [24.31, 54.42], 'C': [11.45, 26.59, 48.06, 64.43, 392.1, 490.0],
          'N': [16.26, 29.14, 47.35, 75.15, 97.79, 552.1, 667.0],
          'O': [14.91, 35.92, 53.70, 75.08, 110.30, 137.97, 739.3, 871.4]}
COLOURS = {'t_r': '#1f77b4', 't_e': '#ff7f0e', 'hc': '#2ca02c', 'conv': '#d62728',
           'H': '#7f7f7f', 'He': '#9467bd', 'N': '#8c564b', 'heavy': '#17becf'}


def find_root(name):
    if os.path.isdir(name):
        dirs = glob(os.path.join(name, 'diag_*'))
        roots = [os.path.basename(d)[5:] for d in dirs if os.path.isdir(d)]
        if len(roots) != 1:
            print('Error: %s should contain one diag_<root> directory' % name)
            sys.exit(1)
        return name, roots[0]
    directory = os.path.dirname(name) or '.'
    root = os.path.basename(name)
    if root.endswith('.pf'):
        root = root[:-3]
    return directory, root


def ncycles(directory, root):
    return len(glob(os.path.join(directory, 'diag_%s' % root, '%s.*.master.txt' % root)))


def table(directory, root, cycle, kind):
    return ascii.read(os.path.join(directory, 'diag_%s' % root, '%s.%02d.%s.txt' % (root, cycle, kind)))


def read_cycle(directory, root, cycle):
    '''one cycle: arrays over wind cells, in the order of the master table'''
    m = table(directory, root, cycle, 'master')
    w = np.array(m['inwind']) == 0
    c = table(directory, root, cycle, 'converge')
    h = table(directory, root, cycle, 'heat')
    out = {'i': np.array(m['i'])[w], 'j': np.array(m['j'])[w], 'converge': np.array(m['converge'])[w],
           't_e': np.array(c['t_e'])[w], 't_e_old': np.array(c['t_e_old'])[w],
           't_r': np.array(c['t_r'])[w], 't_r_old': np.array(c['t_r_old'])[w],
           'heat': np.array(h['heat_tot'])[w] + np.array(h['heat_shoc'])[w], 'cool': np.array(h['cool_tot'])[w],
           'ntot': np.array(m['ntot'])[w].astype(float)}
    for el in ELEMENTS:
        t = table(directory, root, cycle, '%s.frac' % el)
        cols = [k for k in t.colnames if len(k) == 3 and k[0] == 'i' and k[1:].isdigit()]
        out[el] = np.array([np.array(t[k])[w] for k in cols])
    # photons in the band containing each threshold, per cell
    s = table(directory, root, cycle, 'spec')
    key = {(i, j): n for n, (i, j) in enumerate(zip(out['i'], out['j']))}
    bands = sorted(set(s['nband']))
    ncell = len(out['i'])
    nx = np.zeros((len(bands), ncell))
    f1 = np.zeros(len(bands))
    f2 = np.zeros(len(bands))
    for row in s:
        n = key.get((int(row['i']), int(row['j'])))
        if n is None:
            continue
        b = int(row['nband'])
        nx[b, n] = row['nxtot']
        if row['fmax'] > row['fmin']:
            f1[b] = row['fmin'] if f1[b] == 0 else min(f1[b], row['fmin'])
            f2[b] = max(f2[b], row['fmax'])
    out['band_nx'] = nx
    out['band_f1'] = f1
    out['band_f2'] = f2
    # each band's share of the photons where the least was absorbed: the largest share among well-sampled cells
    total = np.maximum(nx.sum(axis=0), 1.0)
    good = total >= 1000
    if not good.any():
        good = total >= np.percentile(total, 90)
    out['band_share'] = np.max(nx[:, good] / total[good], axis=1)
    return out


def change(a, b):
    tot = np.abs(a) + np.abs(b)
    out = np.zeros_like(tot, dtype=float)
    ok = tot > 0
    out[ok] = 2 * np.abs(a[ok] - b[ok]) / tot[ok]
    return out


def delta(p, q):
    tot = p + q
    term = np.zeros_like(tot)
    ok = tot > 0
    term[ok] = 2 * (p[ok] - q[ok]) ** 2 / tot[ok]
    return np.sqrt(term.sum(axis=0))


def photons_at_threshold(cur, el, fion=0.1):
    """
    per cell, the fewest photons the cell would have received, if nothing had been absorbed, in the band holding the
    threshold of any ion of the element with at least fion of the element (thresholds inside the photon range only).
    Expected photons = (photons through the cell) x (the band's share in the least-absorbed, well-sampled cells).
    Few photons because the gas absorbed them is information (J really is low); few because few were generated
    in that band is not. Cells with no such ion get a large number.
    """
    f = cur[el]
    nion = f.shape[0]
    top = np.max(cur['band_f2'])
    lower = np.where(cur['band_f1'] > 0, cur['band_f1'], np.inf)
    out = np.full(f.shape[1], 1e30)
    for k in range(min(nion - 1, len(THRESH[el]))):     # the bare ion has no threshold
        th = THRESH[el][k] / HEV
        if th >= top:
            continue
        below = np.where(lower <= th)[0]
        if len(below) == 0:
            continue
        b = below[np.argmax(lower[below])]
        expected = cur['ntot'] * cur['band_share'][b]
        use = f[k] >= fion
        out[use] = np.minimum(out[use], expected[use])
    return out


def measures(cur, prev, tol, nmin):
    m = {'t_r': change(cur['t_r'], cur['t_r_old']) < 2 * EPSILON,
         't_e': change(cur['t_e'], cur['t_e_old']) < 2 * EPSILON,
         'hc': change(cur['heat'], cur['cool']) < EPSILON,
         'conv': cur['converge'] == 0}
    d = {}
    for el in ('H', 'He', 'N'):
        d[el] = delta(cur[el], prev[el])
    d['heavy'] = np.sqrt(sum(ABUN[e] * delta(cur[e], prev[e]) ** 2 for e in ABUN) / sum(ABUN.values()))
    for k in d:
        m['d_' + k] = d[k] < tol
    info_he = photons_at_threshold(cur, 'He') >= nmin
    info_heavy = np.all([photons_at_threshold(cur, e) >= nmin for e in ('C', 'N', 'O')], axis=0)
    m['vac_He'] = m['d_He'] & ~info_he
    m['vac_heavy'] = m['d_heavy'] & ~info_heavy
    m['delta'] = d
    m['info_He'] = info_he
    m['info_heavy'] = info_heavy
    return m


def by_cycle(directory, root, tol, nmin):
    n = ncycles(directory, root)
    rows = []
    prev = read_cycle(directory, root, 1)
    for k in range(2, n + 1):
        cur = read_cycle(directory, root, k)
        m = measures(cur, prev, tol, nmin)
        rows.append((k, {q: np.mean(v) for q, v in m.items() if q not in ('delta',) and v.dtype == bool}))
        prev = cur
    return rows


def plot_cycles(runs, labels, tol, nmin, outroot):
    fig, axes = plt.subplots(1, len(runs), figsize=(5.2 * len(runs), 4.8), sharey=True, squeeze=False)
    for a, (run, label) in enumerate(zip(runs, labels)):
        ax = axes[0, a]
        rows = by_cycle(*run, tol, nmin)
        x = [r[0] for r in rows]
        for q, name in (('t_r', 't_r'), ('t_e', 't_e'), ('hc', 'heating:cooling'), ('conv', 'Converged')):
            ax.plot(x, [r[1][q] for r in rows], '-', color=COLOURS[q], lw=2.2 if q == 'conv' else 1.5, label=name)
        for q in ('H', 'He', 'N', 'heavy'):
            name = {'H': 'H', 'He': 'He', 'N': 'N', 'heavy': 'C, N, O'}[q]
            ax.plot(x, [r[1]['d_' + q] for r in rows], '--', color=COLOURS[q], lw=1.6,
                    label='$\\delta$%s < %.2f' % (name, tol))
        for q in ('He', 'heavy'):
            name = {'He': 'He', 'heavy': 'C, N, O'}[q]
            ax.plot(x, [r[1]['vac_' + q] for r in rows], ':', color=COLOURS[q], lw=2.0,
                    label='%s passes without information' % name)
        ax.set_title(label, fontsize=10)
        ax.set_xlabel('ionization cycle')
        ax.set_ylim(-0.02, 1.02)
        ax.grid(alpha=0.3)
    axes[0, 0].set_ylabel('fraction of cells in the wind')
    axes[0, -1].legend(fontsize=7, loc='center right')
    fig.suptitle("Sirocco's convergence tests (solid) and the ion changes (dashed, monitored, not part of 'Converged'); "
                 "dotted: passes with fewer than %d photons expected at the threshold" % nmin, fontsize=10)
    fig.tight_layout()
    fig.savefig('%s_cycles.png' % outroot, dpi=100)
    print('Wrote %s_cycles.png' % outroot)


def grid(cur, values):
    nx = int(np.max(cur['i'])) + 1
    nz = int(np.max(cur['j'])) + 1
    g = np.full((nx, nz), np.nan)
    g[cur['i'], cur['j']] = values
    return np.ma.masked_invalid(g)


def draw(ax, cur, values, title, cmap, norm=None):
    g = grid(cur, values)
    cm = copy.copy(cmap)
    cm.set_bad('white')
    im = ax.pcolormesh(np.arange(g.shape[0] + 1) - 0.5, np.arange(g.shape[1] + 1) - 0.5, g.T, cmap=cm, norm=norm,
                       shading='flat')
    ax.set_title(title, fontsize=9)
    ax.set_xlabel('i')
    ax.set_ylabel('j')
    div = make_axes_locatable(ax)
    cax = div.append_axes('right', size='5%', pad=0.05)
    return im, cax


def plot_maps(runs, labels, cycle, tol, nmin, outroot):
    fig, axes = plt.subplots(len(runs), 5, figsize=(22, 4.2 * len(runs)), squeeze=False)
    three = ListedColormap(['#d62728', '#2ca02c', '#ffbf00'])
    for r, (run, label) in enumerate(zip(runs, labels)):
        directory, root = run
        k = cycle if cycle else ncycles(directory, root)
        cur = read_cycle(directory, root, k)
        prev = read_cycle(directory, root, k - 1)
        m = measures(cur, prev, tol, nmin)
        im, cax = draw(axes[r, 0], cur, cur['converge'], '%s, cycle %d\nSirocco: tests failed (%.0f%% converged)' %
                       (label, k, 100 * np.mean(cur['converge'] == 0)), plt.get_cmap('viridis', 4),
                       plt.Normalize(-0.5, 3.5))
        plt.colorbar(im, cax=cax, ticks=[0, 1, 2, 3])
        for c, (q, name) in enumerate((('He', 'He'), ('heavy', 'C, N, O'))):
            dv = m['delta'][q]
            im, cax = draw(axes[r, 1 + 2 * c], cur, np.maximum(dv, 1e-4), '$\\delta$%s since last cycle' % name,
                           plt.cm.coolwarm, matplotlib.colors.LogNorm(1e-3, 1))
            plt.colorbar(im, cax=cax)
            state = np.where(~m['d_' + q], 0, np.where(m['info_' + q], 1, 2))
            frac = [np.mean(state == s) for s in range(3)]
            im, cax = draw(axes[r, 2 + 2 * c], cur, state,
                           '%s: fails %.0f%%, passes %.0f%%,\npasses without information %.0f%%' %
                           (name, 100 * frac[0], 100 * frac[1], 100 * frac[2]), three, plt.Normalize(-0.5, 2.5))
            cb = plt.colorbar(im, cax=cax, ticks=[0, 1, 2])
            cb.ax.set_yticklabels(['fails', 'passes', 'no info'], fontsize=8)
    fig.suptitle('Maps for one cycle: tolerance %.2f on the ion changes; "no information" = fewer than %d photons '
                 'expected (unabsorbed) at the threshold of an ion with at least 10%% of its element' % (tol, nmin), fontsize=11)
    fig.tight_layout()
    fig.savefig('%s_maps.png' % outroot, dpi=90)
    print('Wrote %s_maps.png' % outroot)


def steer(argv):
    tol = 0.1
    nmin = 10
    cycle = 0
    labels = []
    outroot = 'ion_monitor'
    names = []
    i = 1
    while i < len(argv):
        if argv[i] == '-h':
            print(__doc__)
            return
        elif argv[i] == '-tol':
            i += 1
            tol = float(argv[i])
        elif argv[i] == '-nmin':
            i += 1
            nmin = int(argv[i])
        elif argv[i] == '-cycle':
            i += 1
            cycle = int(argv[i])
        elif argv[i] == '-labels':
            i += 1
            labels = argv[i].split(';')
        elif argv[i] == '-o':
            i += 1
            outroot = argv[i]
        elif argv[i].startswith('-'):
            print('Error: unknown option %s' % argv[i])
            return
        else:
            names.append(argv[i])
        i += 1
    if not names:
        print(__doc__)
        return
    runs = [find_root(n) for n in names]
    if len(labels) != len(runs):
        labels = [os.path.basename(os.path.normpath(n)) for n in names]
    plot_cycles(runs, labels, tol, nmin, outroot)
    plot_maps(runs, labels, cycle, tol, nmin, outroot)


if __name__ == '__main__':
    steer(sys.argv)
