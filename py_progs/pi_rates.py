#!/usr/bin/env python
'''
Compare the photoionization rates from the band model with those from the fine cell
histogram and from the direct estimators, as written by sirocco -pi_rate_table.


Command line usage:
    pi_rates.py [options] run [run2 ...]

    run is a run directory, or a root name (model or model.pf) in the current directory.
    The rates are read from root.pi_rates.*.txt (one file per MPI rank).

Options:
    -h              print this documentation
    -labels a,b     labels for the runs in the figures (default: the directory names)
    -cycles n1:n2   the ionization cycles to use, counting from 1 as Sirocco does in its
                    output (default: the second half of the cycles in the file)
    -fmin f         only ions that are at least this fraction of their element (default 1e-3)
    -elements H,He  only these elements (default all in the file)
    -ion "O III"    also plot this ion's rate ratios against radius
    -o root         root for the output files (default pi_rates)

Description:
    sirocco -pi_rate_table writes, for every cell, ion and ionization cycle, three
    photoionization rates computed from the same photons:
        rate_bands   from the band model (Wind.ionization.band_model)
        rate_hist    from the fine histogram of the cell spectrum (-cell_spec_dim bins, plus
                     the photoionization thresholds), with J_nu taken as constant in each bin
        rate_direct  directly from the photon path lengths and the cross section
    With enough photons in each cell, rate_hist and rate_direct agree, and rate_bands/rate_hist
    measures how well the band model describes the spectrum (its representation error).

    pi_rates.py writes:
        <root>.png      one row per ion: the median of log10(bands/hist) over cells and cycles
                        (dot, with the 16-84% range as a bar), and the median of
                        log10(direct/hist) (cross); one colour per run
        <root>.txt      the same numbers as a table
        <root>_<ion>.png    with -ion, log10(bands/hist) and log10(direct/hist) against radius
                        (median over the cycles), one colour per run

    The ion fractions used for -fmin come from the per-cycle tables in diag_<root>/
    (sirocco -d, then windsave2table, which writes root.NN.<element>.frac.txt), or, if
    those are missing, from the final root.<element>.frac.txt. The table root.master.txt
    (or the per-cycle master) gives the radius of each plasma cell.

Primary routines:
    read_rates      reads all the rate files of one run
    ion_fractions   attaches the ion fraction and radius of each row
    summary         statistics of the ratios for each ion
    steer           processes the command line

Notes:
    Written 2026-11-15 for the banding comparison on the common test wind (Sirocco-2611)
    and for PR #1250. Cycles in the rate file (geo.wcycle) count from 0; the diag tables
    count from 1, and so do the -cycles option and the output.
'''

import sys
import os
from glob import glob
import numpy as np
from astropy.io import ascii
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ELEMENT = {1: 'H', 2: 'He', 3: 'Li', 4: 'Be', 5: 'B', 6: 'C', 7: 'N', 8: 'O', 9: 'F', 10: 'Ne', 11: 'Na', 12: 'Mg',
           13: 'Al', 14: 'Si', 15: 'P', 16: 'S', 17: 'Cl', 18: 'Ar', 19: 'K', 20: 'Ca', 26: 'Fe', 28: 'Ni'}
ROMAN = ['I', 'II', 'III', 'IV', 'V', 'VI', 'VII', 'VIII', 'IX', 'X', 'XI', 'XII', 'XIII', 'XIV', 'XV', 'XVI',
         'XVII', 'XVIII', 'XIX', 'XX', 'XXI', 'XXII', 'XXIII', 'XXIV', 'XXV', 'XXVI', 'XXVII']
COLOURS = ['#1f4e79', '#c0504d', '#4f8f3a', '#e69f00', '#7f5f9f', '#56b4e9', '#444444']


def ion_name(z, istate):
    return '%s %s' % (ELEMENT.get(z, 'Z%d' % z), ROMAN[istate - 1])


def find_root(name):
    '''(directory, root) for a run directory, a root name, or a .pf file'''
    if os.path.isdir(name):
        files = glob(os.path.join(name, '*.pi_rates.*.txt'))
        roots = sorted(set(os.path.basename(f).split('.pi_rates.')[0] for f in files))
        if len(roots) != 1:
            print('Error: %s should contain the rate files of exactly one model, not %d' % (name, len(roots)))
            sys.exit(1)
        return name, roots[0]
    directory = os.path.dirname(name)
    if directory == '':
        directory = '.'
    root = os.path.basename(name)
    if root.endswith('.pf'):
        root = root[:-3]
    return directory, root


def read_rates(directory, root):
    '''all rows of root.pi_rates.*.txt, as a dictionary of arrays'''
    rows = []
    files = sorted(glob(os.path.join(directory, '%s.pi_rates.*.txt' % root)))
    if len(files) == 0:
        print('Error: no %s.pi_rates.*.txt in %s (run sirocco with -pi_rate_table)' % (root, directory))
        sys.exit(1)
    for filename in files:
        for line in open(filename):
            words = line.split()
            if len(words) == 0 or not words[0].isdigit():
                continue
            rows.append([float(w) for w in words])
    r = np.array(rows)
    return {'cycle': r[:, 0].astype(int) + 1, 'nplasma': r[:, 1].astype(int), 'nion': r[:, 2].astype(int),
            'z': r[:, 3].astype(int), 'istate': r[:, 4].astype(int), 'density': r[:, 5],
            'bands': r[:, 6], 'hist': r[:, 7], 'direct': r[:, 8]}


def read_table(directory, root, kind, cycle):
    '''a windsave2table table for one cycle (diag_<root>/), or the final one; None if missing'''
    names = [os.path.join(directory, 'diag_%s' % root, '%s.%02d.%s.txt' % (root, cycle, kind)),
             os.path.join(directory, '%s.%s.txt' % (root, kind))]
    for filename in names:
        if os.path.isfile(filename):
            return ascii.read(filename)
    return None


def ion_fractions(rates, directory, root, cycles):
    '''
    Attach to each row the fraction of its element in that ion, and the radius of its cell.
    Rows of cycles whose tables are missing get the final tables.
    '''
    nrow = len(rates['cycle'])
    frac = np.full(nrow, np.nan)
    radius = np.full(nrow, np.nan)
    for cycle in cycles:
        sel = np.where(rates['cycle'] == cycle)[0]
        if len(sel) == 0:
            continue
        master = read_table(directory, root, 'master', cycle)
        if master is None:
            print('Error: no master table for %s; run windsave2table (and sirocco -d for each cycle)' % root)
            sys.exit(1)
        plasma_col = 'nplasma'
        if plasma_col not in master.colnames:
            plasma_col = 'iplasma'
        inwind = np.array(master['inwind']) >= 0
        row_of_plasma = {}
        for k in np.where(inwind)[0]:
            row_of_plasma[int(master[plasma_col][k])] = k
        if 'r' in master.colnames:
            rcol = np.array(master['r'])
        else:
            rcol = np.sqrt(np.array(master['x']) ** 2 + np.array(master['z']) ** 2)
        for z in sorted(set(rates['z'][sel])):
            name = ELEMENT.get(z, '')
            table = read_table(directory, root, '%s.frac' % name, cycle)
            if table is None:
                continue
            for k in sel[rates['z'][sel] == z]:
                m = row_of_plasma.get(rates['nplasma'][k])
                if m is None:
                    continue
                col = 'i%02d' % rates['istate'][k]
                if col in table.colnames:
                    frac[k] = table[col][m]
                radius[k] = rcol[m]
    rates['frac'] = frac
    rates['r'] = radius
    return rates


def summary(rates, cycles, fmin, elements):
    '''for each ion: median and 16-84% of log10(bands/hist), median of log10(direct/hist), count'''
    use = np.isin(rates['cycle'], cycles) & (rates['frac'] >= fmin) & (rates['hist'] > 0) & (rates['bands'] > 0)
    if elements:
        use &= np.isin([ELEMENT.get(z, '') for z in rates['z']], elements)
    out = {}
    for z, istate in sorted(set(zip(rates['z'][use], rates['istate'][use]))):
        s = use & (rates['z'] == z) & (rates['istate'] == istate)
        if s.sum() < 3:
            continue
        lb = np.log10(rates['bands'][s] / rates['hist'][s])
        good = rates['direct'][s] > 0
        if good.any():
            ld = np.median(np.log10(rates['direct'][s][good] / rates['hist'][s][good]))
        else:
            ld = np.nan
        out[(z, istate)] = (np.median(lb), np.percentile(lb, 16), np.percentile(lb, 84), ld, int(s.sum()))
    return out


def plot_summary(results, labels, outroot):
    ions = sorted(set(k for res in results for k in res))
    nrun = len(results)
    fig, ax = plt.subplots(figsize=(7.5, 0.45 * len(ions) + 1.8))
    ax.axvspan(-0.05, 0.05, color='0.9', zorder=0)
    ax.axvline(0, color='k', lw=0.8)
    for n, res in enumerate(results):
        offset = 0.0
        if nrun > 1:
            offset = 0.6 * (n / (nrun - 1) - 0.5)
        for k, key in enumerate(ions):
            if key not in res:
                continue
            med, lo, hi, ld, count = res[key]
            y = k + offset
            first = (k == min(i for i, kk in enumerate(ions) if kk in res))
            ax.plot([lo, hi], [y, y], color=COLOURS[n % len(COLOURS)], lw=2)
            ax.plot(med, y, 'o', color=COLOURS[n % len(COLOURS)], ms=6, label=(labels[n] if first else None))
            if np.isfinite(ld):
                ax.plot(ld, y, 'x', color=COLOURS[n % len(COLOURS)], ms=6,
                        label=('direct / hist' if (first and n == 0) else None))
    ax.set_yticks(range(len(ions)))
    ax.set_yticklabels([ion_name(z, i) for z, i in ions])
    ax.set_ylim(-0.7, len(ions) - 0.3)
    ax.set_xlabel('log$_{10}$ (rate / rate from the cell histogram)')
    ax.set_title('Photoionization rates: band model (dots, 16-84% bars) and direct (crosses)\n'
                 'against the fine histogram; grey = within 0.05 dex', fontsize=10)
    ax.grid(alpha=0.3, axis='x')
    ax.legend(fontsize=8, loc='best')
    fig.tight_layout()
    fig.savefig('%s.png' % outroot, dpi=100)
    print('Wrote %s.png' % outroot)


def plot_ion(all_rates, labels, cycles, ion, outroot):
    words = ion.split()
    zmatch = [z for z, e in ELEMENT.items() if e == words[0]]
    if len(words) != 2 or len(zmatch) != 1 or words[1] not in ROMAN:
        print('Error: give the ion as e.g. "O III", not %s' % ion)
        return
    z = zmatch[0]
    istate = ROMAN.index(words[1]) + 1
    fig, ax = plt.subplots(figsize=(7.5, 4.5))
    ax.axhspan(-0.05, 0.05, color='0.9', zorder=0)
    ax.axhline(0, color='k', lw=0.8)
    for n, rates in enumerate(all_rates):
        s = np.isin(rates['cycle'], cycles) & (rates['z'] == z) & (rates['istate'] == istate) & (rates['hist'] > 0)
        cells = sorted(set(rates['nplasma'][s]))
        r = []
        lb = []
        ld = []
        for c in cells:
            k = s & (rates['nplasma'] == c)
            r.append(np.median(rates['r'][k]))
            lb.append(np.median(np.log10(np.maximum(rates['bands'][k], 1e-300) / rates['hist'][k])))
            ld.append(np.median(np.log10(np.maximum(rates['direct'][k], 1e-300) / rates['hist'][k])))
        order = np.argsort(r)
        r = np.array(r)[order]
        colour = COLOURS[n % len(COLOURS)]
        ax.plot(r, np.array(lb)[order], '-o', color=colour, ms=3, label='%s: bands / hist' % labels[n])
        ax.plot(r, np.array(ld)[order], ':', color=colour, label='%s: direct / hist' % labels[n])
    ax.set_xscale('log')
    ax.set_xlabel('radius (cm)')
    ax.set_ylabel('log$_{10}$ (rate / rate from histogram)')
    ax.set_title('%s photoionization rate against radius (median over cycles %d-%d)' % (ion, cycles[0], cycles[-1]),
                 fontsize=10)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    name = '%s_%s.png' % (outroot, ion.replace(' ', ''))
    fig.savefig(name, dpi=100)
    print('Wrote %s' % name)


def steer(argv):
    labels = []
    cycle_range = ''
    fmin = 1e-3
    elements = []
    ion = ''
    outroot = 'pi_rates'
    runs = []

    i = 1
    while i < len(argv):
        if argv[i] == '-h':
            print(__doc__)
            return
        elif argv[i] == '-labels':
            i += 1
            labels = argv[i].split(',')
        elif argv[i] == '-cycles':
            i += 1
            cycle_range = argv[i]
        elif argv[i] == '-fmin':
            i += 1
            fmin = float(argv[i])
        elif argv[i] == '-elements':
            i += 1
            elements = argv[i].split(',')
        elif argv[i] == '-ion':
            i += 1
            ion = argv[i]
        elif argv[i] == '-o':
            i += 1
            outroot = argv[i]
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

    all_rates = []
    results = []
    lines = ['# pi_rates.py: log10(rate_bands/rate_hist) median, 16%, 84%; log10(rate_direct/rate_hist) median; n',
             '# fmin %g' % fmin]
    cycles = None
    for run, label in zip(runs, labels):
        directory, root = find_root(run)
        rates = read_rates(directory, root)
        ncycle = int(np.max(rates['cycle']))
        if cycle_range:
            a, b = cycle_range.split(':')
            run_cycles = list(range(int(a), int(b) + 1))
        else:
            run_cycles = list(range(ncycle // 2 + 1, ncycle + 1))
        if cycles is None:
            cycles = run_cycles
        rates = ion_fractions(rates, directory, root, run_cycles)
        res = summary(rates, run_cycles, fmin, elements)
        all_rates.append(rates)
        results.append(res)
        lines.append('## %s (%s, cycles %d-%d)' % (label, run, run_cycles[0], run_cycles[-1]))
        for (z, istate), (med, lo, hi, ld, count) in sorted(res.items()):
            lines.append('%-8s %+.3f %+.3f %+.3f   direct %+.3f   n %d' % (ion_name(z, istate), med, lo, hi, ld, count))

    open('%s.txt' % outroot, 'w').write('\n'.join(lines) + '\n')
    print('Wrote %s.txt' % outroot)
    plot_summary(results, labels, outroot)
    if ion:
        plot_ion(all_rates, labels, cycles, ion, outroot)


if __name__ == '__main__':
    steer(sys.argv)
