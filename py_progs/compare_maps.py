#!/usr/bin/env python
'''
Maps of wind quantities, side by side, for one run or for two runs and their difference.


Command line usage:
    compare_maps.py [options] root
    compare_maps.py [options] rootA rootB

    root may be a root name (model, or model.pf), or a run directory containing a single
    .pf file. The tables are those written by windsave2table; with -cycle N the per-cycle
    tables in diag_<root>/ written with sirocco -d (root.NN.master.txt etc.) are used.

Options:
    -h              print this documentation
    -vars a,b,c     the quantities to show (default t_e,t_r,converge). Any column of the
                    master, heat or converge tables, plus the derived quantities below.
    -tests          one run: show Sirocco's convergence tests cell by cell instead, i.e.
                    the change in t_r and t_e since the previous cycle, the heating-cooling
                    imbalance, and the number of tests failed (converge)
    -cycle N        use ionization cycle N (the per-cycle tables in diag_<root>/)
    -cycleB M       for two runs, the cycle of the second run (default the same as -cycle);
                    with one root and -cycleB, the two cycles of one run are compared
    -grid ij|xz|log axes: cell indices (default, as in run_check.py), x and z, or log x, log z
    -frac           show the difference of positive quantities as 2 (B - A) / (A + B) instead of
                    log10(B / A); best for small differences, e.g. two successive cycles
    -o file         output file (default compare_maps.png, or <root>_tests.png with -tests)

Description:
    One run: a row of maps, one per quantity.
    Two runs (or two cycles): one row per quantity, with columns A, B and the difference.
    A and B share one colour scale. The difference is
        log10(B / A)          for positive quantities (t_e, t_r, ne, ...), in dex;
        2 (B - A) / (A + B)   instead with -frac: the fractional change Sirocco uses in its
                              convergence test, times 2 (it saturates at 2 for large ratios);
        B - A                 for converge and other quantities that can be zero or negative.

    Derived quantities (computed from the converge table of a single cycle):
        dt_r   2 |t_r - t_r_old| / (t_r + t_r_old)     change in t_r since the previous cycle
        dt_e   2 |t_e - t_e_old| / (t_e + t_e_old)     change in t_e since the previous cycle
        hc     2 |heat_tot - cool_tot| / (heat_tot + cool_tot)   heating-cooling imbalance
    Sirocco calls a cell converged when dt_r and dt_e are below 2 epsilon (epsilon = 0.05)
    and hc is below epsilon (see convergence.c); converge is the number of these tests failed
    (0 = converged, 3 = all failed).  With -tests the colour scale for the change maps has the
    tolerance at its middle, so cells that pass are blue and cells that fail are red.

    Cells outside the wind are left white.

Primary routines:
    get_cells   reads the tables for one run (and cycle) into one dictionary keyed by column
    one_run     maps of one run
    two_runs    maps of two runs and their difference
    steer       processes the command line

Notes:
    Written 2026-11-15 to compare band models and to show the convergence tests on the
    common test wind (Sirocco-2611).  It complements plot_wind.py (one quantity) and
    run_check.py (convergence by cycle and the converge map of one run).
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
from matplotlib.colors import LogNorm, TwoSlopeNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable

EPSILON = 0.05

LABELS = {
    't_e': 'electron temperature (K)',
    't_r': 'radiation temperature (K)',
    'ne': 'electron density (cm$^{-3}$)',
    'converge': 'tests failed (0 = converged)',
    'ntot': 'photons in cell',
    'dt_r': 'change in t_r since last cycle',
    'dt_e': 'change in t_e since last cycle',
    'hc': 'heating-cooling imbalance',
}

TOLERANCE = {'dt_r': 2 * EPSILON, 'dt_e': 2 * EPSILON, 'hc': EPSILON}

SIGNED = ('converge', 'inwind', 'i', 'j', 'v_x', 'v_y', 'v_z', 'dmo_dt_x', 'dmo_dt_y', 'dmo_dt_z')


def find_root(name):
    '''
    Return (directory, root) for a root name, a .pf file, or a run directory
    containing one .pf file (*.out.pf files are ignored)
    '''
    if os.path.isdir(name):
        pfs = [f for f in glob(os.path.join(name, '*.pf')) if not f.endswith('.out.pf')]
        roots = [os.path.basename(f)[:-3] for f in pfs]
        # several .pf files: keep those that have tables
        if len(roots) > 1:
            roots = [r for r in roots if os.path.isfile(os.path.join(name, '%s.master.txt' % r))
                     or os.path.isdir(os.path.join(name, 'diag_%s' % r))]
        if len(roots) != 1:
            print('Error: cannot tell which model in %s to use (%s); give the root name' % (name, ' '.join(roots)))
            sys.exit(1)
        return name, roots[0]
    directory = os.path.dirname(name)
    if directory == '':
        directory = '.'
    root = os.path.basename(name)
    if root.endswith('.pf'):
        root = root[:-3]
    return directory, root


def table_name(directory, root, kind, cycle):
    if cycle is None:
        return os.path.join(directory, '%s.%s.txt' % (root, kind))
    return os.path.join(directory, 'diag_%s' % root, '%s.%02d.%s.txt' % (root, cycle, kind))


def get_cells(name, cycle=None):
    '''
    Read the master, heat and converge tables of one run (and cycle) into one dictionary
    of 2d arrays [i, j], and add the derived quantities dt_r, dt_e and hc.
    '''
    directory, root = find_root(name)
    cells = {}
    for kind in ('master', 'heat', 'converge'):
        filename = table_name(directory, root, kind, cycle)
        if not os.path.isfile(filename):
            if kind == 'master':
                print('Error: %s not found' % filename)
                sys.exit(1)
            continue
        t = ascii.read(filename)
        for col in t.colnames:
            if col not in cells:
                cells[col] = np.array(t[col])

    if 'j' not in cells:
        print('Error: %s is not a 2d model; compare_maps.py makes maps of 2d models' % name)
        sys.exit(1)

    nx = int(np.max(cells['i'])) + 1
    nz = int(np.max(cells['j'])) + 1
    for col in list(cells.keys()):
        if len(cells[col]) == nx * nz:
            cells[col] = cells[col].reshape(nx, nz)
        else:
            del cells[col]

    if 't_r_old' in cells:
        cells['dt_r'] = change(cells['t_r'], cells['t_r_old'])
        cells['dt_e'] = change(cells['t_e'], cells['t_e_old'])
    if 'heat_tot' in cells and 'cool_tot' in cells:
        cells['hc'] = change(cells['heat_tot'], cells['cool_tot'])

    label = os.path.basename(os.path.normpath(name))
    if cycle is not None:
        label = '%s, cycle %d' % (label, cycle)
    cells['label'] = label
    return cells


def change(a, b):
    '''2 |a - b| / (a + b), zero where both are zero'''
    a = np.abs(np.asarray(a, dtype=float))
    b = np.abs(np.asarray(b, dtype=float))
    total = a + b
    out = np.zeros_like(total)
    ok = total > 0
    out[ok] = 2.0 * np.abs(a[ok] - b[ok]) / total[ok]
    return out


def axes_for(cells, grid):
    if grid == 'xz':
        return cells['x'], cells['z'], 'x (cm)', 'z (cm)'
    if grid == 'log':
        x = np.log10(np.maximum(cells['x'], np.min(cells['x'][cells['x'] > 0])))
        z = np.log10(np.maximum(cells['z'], np.min(cells['z'][cells['z'] > 0])))
        return x, z, 'log x (cm)', 'log z (cm)'
    return cells['i'], cells['j'], 'i', 'j'


def masked(cells, var):
    if var not in cells:
        print('Error: %s is not in the tables of %s' % (var, cells['label']))
        sys.exit(1)
    return np.ma.array(np.asarray(cells[var], dtype=float), mask=(cells['inwind'] != 0))


def colour_scale(values, var):
    '''a Normalize for one quantity: log if it is positive and spans more than a factor 30'''
    good = values.compressed()
    if len(good) == 0:
        return None
    vmin = np.min(good)
    vmax = np.max(good)
    if var in TOLERANCE:
        tol = TOLERANCE[var]
        top = max(vmax, 2 * tol)
        return TwoSlopeNorm(vcenter=tol, vmin=0.0, vmax=top)
    if var in SIGNED or vmin <= 0:
        return plt.Normalize(vmin=vmin, vmax=vmax)
    if vmax / vmin > 30:
        return LogNorm(vmin=vmin, vmax=vmax)
    return plt.Normalize(vmin=vmin, vmax=vmax)


def colour_map(var):
    if var in TOLERANCE:
        cmap = copy.copy(plt.cm.coolwarm)
    elif var == 'converge':
        cmap = copy.copy(plt.get_cmap('viridis', 4))
    else:
        cmap = copy.copy(plt.cm.viridis)
    cmap.set_bad(color='white')
    return cmap


def draw(ax, x, y, values, norm, cmap, title, xlabel, ylabel):
    im = ax.pcolormesh(x, y, values, norm=norm, cmap=cmap, shading='auto')
    ax.set_title(title, fontsize=10)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    plt.colorbar(im, cax=cax)
    return im


def one_run(cells, variables, grid, outfile):
    '''a row of maps of one run, one per quantity'''
    n = len(variables)
    fig, axes = plt.subplots(1, n, figsize=(4.6 * n, 4.4), squeeze=False)
    x, y, xlabel, ylabel = axes_for(cells, grid)
    for k, var in enumerate(variables):
        values = masked(cells, var)
        title = LABELS.get(var, var)
        if var in TOLERANCE:
            frac = np.mean(values.compressed() < TOLERANCE[var])
            title = '%s\n(%.0f%% of cells below %.2f)' % (title, 100 * frac, TOLERANCE[var])
        if var == 'converge':
            frac = np.mean(values.compressed() == 0)
            title = '%s\n(%.0f%% converged)' % (title, 100 * frac)
        draw(axes[0, k], x, y, values, colour_scale(values, var), colour_map(var), title, xlabel, ylabel)
    fig.suptitle(cells['label'], fontsize=12)
    fig.tight_layout()
    fig.savefig(outfile, dpi=90)
    print('Wrote %s' % outfile)


def two_runs(a, b, variables, grid, outfile, frac=False):
    '''one row per quantity: run A, run B, and the difference'''
    n = len(variables)
    fig, axes = plt.subplots(n, 3, figsize=(14, 4.2 * n), squeeze=False)
    x, y, xlabel, ylabel = axes_for(a, grid)
    for k, var in enumerate(variables):
        va = masked(a, var)
        vb = masked(b, var)
        if va.shape != vb.shape:
            print('Error: the two runs have different grids')
            sys.exit(1)
        both = np.ma.concatenate([va.flatten(), vb.flatten()])
        norm = colour_scale(both, var)
        cmap = colour_map(var)
        label = LABELS.get(var, var)
        draw(axes[k, 0], x, y, va, norm, cmap, 'A: %s\n%s' % (a['label'], label), xlabel, ylabel)
        draw(axes[k, 1], x, y, vb, norm, cmap, 'B: %s\n%s' % (b['label'], label), xlabel, ylabel)

        if var in SIGNED or np.ma.min(both) <= 0:
            diff = vb - va
            dlabel = 'B - A'
        elif not frac:
            diff = np.ma.log10(vb / va)
            dlabel = 'log10(B / A)'
        else:
            total = va + vb
            diff = np.ma.where(total > 0, 2.0 * (vb - va) / np.ma.where(total > 0, total, 1.0), 0.0)
            dlabel = '2 (B - A) / (A + B)'
        good = diff.compressed()
        top = np.max(np.abs(good)) if len(good) else 1.0
        if top == 0:
            top = 1.0
        dmap = copy.copy(plt.cm.RdBu_r)
        dmap.set_bad(color='white')
        title = '%s, %s\n(median |diff| %.3g, max %.3g)' % (label, dlabel, np.median(np.abs(good)), top)
        draw(axes[k, 2], x, y, diff, plt.Normalize(vmin=-top, vmax=top), dmap, title, xlabel, ylabel)
    fig.tight_layout()
    fig.savefig(outfile, dpi=90)
    print('Wrote %s' % outfile)


def steer(argv):
    variables = ['t_e', 't_r', 'converge']
    tests = False
    cycle = None
    cycle_b = None
    grid = 'ij'
    outfile = ''
    frac = False
    roots = []

    i = 1
    while i < len(argv):
        if argv[i] == '-h':
            print(__doc__)
            return
        elif argv[i] == '-vars':
            i += 1
            variables = argv[i].split(',')
        elif argv[i] == '-frac':
            frac = True
        elif argv[i] == '-tests':
            tests = True
        elif argv[i] == '-cycle':
            i += 1
            cycle = int(argv[i])
        elif argv[i] == '-cycleB':
            i += 1
            cycle_b = int(argv[i])
        elif argv[i] == '-grid':
            i += 1
            grid = argv[i]
        elif argv[i] == '-o':
            i += 1
            outfile = argv[i]
        elif argv[i].startswith('-'):
            print('Error: unknown option %s' % argv[i])
            print(__doc__)
            return
        else:
            roots.append(argv[i])
        i += 1

    if len(roots) == 0 or len(roots) > 2:
        print(__doc__)
        return

    if cycle_b is None:
        cycle_b = cycle

    if len(roots) == 1 and (cycle_b == cycle or tests):
        cells = get_cells(roots[0], cycle)
        if tests:
            variables = ['dt_r', 'dt_e', 'hc', 'converge']
            if outfile == '':
                outfile = '%s_tests.png' % find_root(roots[0])[1]
        if outfile == '':
            outfile = 'compare_maps.png'
        one_run(cells, variables, grid, outfile)
        return

    if len(roots) == 1:
        a = get_cells(roots[0], cycle)
        b = get_cells(roots[0], cycle_b)
    else:
        a = get_cells(roots[0], cycle)
        b = get_cells(roots[1], cycle_b)
    if outfile == '':
        outfile = 'compare_maps.png'
    two_runs(a, b, variables, grid, outfile, frac)


if __name__ == '__main__':
    steer(sys.argv)
