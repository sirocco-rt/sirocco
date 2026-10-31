#!/usr/bin/env python
'''
Read the Q and U spectra that Sirocco writes when polarization is on, and
summarize the polarization integrated over a band.


Command line usage (if any):

    usage::

        c_spectra.py [-h] [-wmin WMIN] [-wmax WMAX] [-compare FILE]
                     [-plot FILE] [-nbin NBIN] spec [spec ...]

    spec        Sirocco spectrum files (.spec, .spec_wind, .log_spec, ...)
    -wmin       shortest wavelength (A) to include (default: all)
    -wmax       longest wavelength (A) to include (default: all)
    -compare    a step4_regress.py results table; matching models and angles
                are printed alongside
    -plot       write a plot of I, q and u against wavelength to FILE
                (only for the first spectrum file)
    -nbin       combine this many wavelength bins in the plot (default 50)

Description:

    When Sirocco is run with @Spectrum.polarization yes (which needs -d),
    the spectrum files have two extra columns for each observer, e.g.
    A80P0.50_Q and A80P0.50_U, after all of the usual columns.  Q and U are
    in the same units as I and are referred to North on the sky, the
    projected +z axis (IAU convention), so Q > 0 is an electric vector
    parallel to the projected symmetry axis.

    The band-integrated values are sums of L_nu dnu over the bins, so the
    flambda or fnu units of the files are first converted back to L_nu dnu.
    q = Q/I and u = U/I.  For an axisymmetric model U should be zero to
    within the noise, which gives a rough measure of the noise in q.

Primary routines:

    read_spec, band_pol, rebin_pol

Notes:

    Requires numpy, and matplotlib for -plot.

History:

    261013 ksl Coding begun, with Claude Code

'''

import argparse
import os
import sys

import numpy as np


def read_spec(fname):
    '''
    Read a Sirocco spectrum file.

    Returns (units, freq, wave, cols), where units is 'flambda', 'fnu' or
    'lnu', and cols is a dict from column name to array.
    '''
    units = None
    names = None
    rows = []
    with open(fname) as f:
        for line in f:
            if line.startswith('#'):
                if 'Units:' in line:
                    if 'flambda' in line:
                        units = 'flambda'
                    elif 'Fnu' in line:
                        units = 'fnu'
                    else:
                        units = 'lnu'
                continue
            words = line.split()
            if not words:
                continue
            if words[0] == 'Freq.':
                names = words
                continue
            if names is not None:
                rows.append([float(w) for w in words])
    if names is None:
        raise ValueError(f'{fname}: no column header found')
    data = np.array(rows)
    cols = {name: data[:, i] for i, name in enumerate(names)}
    return units, cols['Freq.'], cols['Lambda'], cols


def lnu_dnu(units, freq, logbins):
    '''The factor that turns a column in these units into L_nu dnu per bin (up to a constant).'''
    if units == 'flambda':
        x = 1.0 / freq ** 2
    else:
        x = np.ones_like(freq)
    if logbins:
        x = x * freq
    return x


def observers(cols):
    '''The names of the observer spectra that have Q and U columns.'''
    return [name[:-2] for name in cols if name.endswith('_Q') and name[:-2] in cols]


def band_pol(fname, wmin=None, wmax=None):
    '''
    Band-integrated polarization for each observer in a spectrum file.

    Returns a list of dicts with name, angle, I, q, u (in %), P (%) and PA (deg).
    '''
    units, freq, wave, cols = read_spec(fname)
    logbins = 'log_spec' in os.path.basename(fname)
    sel = np.ones(len(freq), dtype=bool)
    if wmin is not None:
        sel &= wave >= wmin
    if wmax is not None:
        sel &= wave <= wmax
    fac = lnu_dnu(units, freq, logbins)[sel]

    out = []
    for name in observers(cols):
        I = np.sum(cols[name][sel] * fac)
        Q = np.sum(cols[name + '_Q'][sel] * fac)
        U = np.sum(cols[name + '_U'][sel] * fac)
        q = 100 * Q / I if I > 0 else np.nan
        u = 100 * U / I if I > 0 else np.nan
        angle = float(name[1:].split('P')[0])
        out.append(dict(name=name, angle=angle, I=I, q=q, u=u,
                        P=np.hypot(q, u), PA=np.degrees(0.5 * np.arctan2(u, q)) % 180))
    return out


def rebin_pol(fname, name, nbin=50):
    '''Wavelength, I, q and u (%) for one observer, combining nbin bins at a time.'''
    units, freq, wave, cols = read_spec(fname)
    fac = lnu_dnu(units, freq, 'log_spec' in os.path.basename(fname))
    n = (len(freq) // nbin) * nbin

    def comb(x):
        return x[:n].reshape(-1, nbin).sum(axis=1)

    I = comb(cols[name] * fac)
    Q = comb(cols[name + '_Q'] * fac)
    U = comb(cols[name + '_U'] * fac)
    w = wave[:n].reshape(-1, nbin).mean(axis=1)
    with np.errstate(invalid='ignore', divide='ignore'):
        return w, I, 100 * Q / I, 100 * U / I


def read_step4(fname):
    '''Results of step4_regress.py: {(model, angle): (q, sigma_q, u, sigma_u)} in %.'''
    res = {}
    model = None
    with open(fname) as f:
        for line in f:
            words = line.split()
            if not words or line.startswith('#'):
                continue
            if '(grid:' in line:
                model = words[0]
                continue
            if model and len(words) >= 5 and '±' in words[3]:
                try:
                    angle = float(words[0])
                except ValueError:
                    continue
                q, sq = (float(x) for x in words[3].split('±'))
                u, su = (float(x) for x in words[4].split('±'))
                res[(model, angle)] = (q, sq, u, su)
    return res


def plot(fname, outfile, nbin):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    units, freq, wave, cols = read_spec(fname)
    names = observers(cols)
    fig, ax = plt.subplots(3, 1, figsize=(8, 9), sharex=True)
    for name in names:
        w, I, q, u = rebin_pol(fname, name, nbin)
        ax[0].plot(w, I / np.nanmax(I), label=name)
        ax[1].plot(w, q, label=name)
        ax[2].plot(w, u, label=name)
    ax[0].set_ylabel('I (normalized)')
    ax[1].set_ylabel('q (%)')
    ax[2].set_ylabel('u (%)')
    ax[2].set_xlabel('wavelength (A)')
    for a in ax[1:]:
        a.axhline(0, color='k', lw=0.5)
    ax[0].legend(fontsize=8)
    ax[0].set_title(os.path.basename(fname))
    fig.tight_layout()
    fig.savefig(outfile, dpi=120)


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[1].strip())
    parser.add_argument('spec', nargs='+')
    parser.add_argument('-wmin', type=float, default=None)
    parser.add_argument('-wmax', type=float, default=None)
    parser.add_argument('-compare', default=None)
    parser.add_argument('-plot', default=None)
    parser.add_argument('-nbin', type=int, default=50)
    args = parser.parse_args()

    step4 = read_step4(args.compare) if args.compare else {}

    for fname in args.spec:
        model = os.path.basename(fname).split('.')[0]
        print(f'\n{fname}')
        head = f'  {"observer":<10} {"q (%)":>8} {"u (%)":>8} {"P (%)":>7} {"PA":>6}'
        if step4:
            head += f'   {"Step 4 q (%)":>16} {"u (%)":>14}'
        print(head)
        for r in band_pol(fname, args.wmin, args.wmax):
            line = f'  {r["name"]:<10} {r["q"]:+8.3f} {r["u"]:+8.3f} {r["P"]:7.3f} {r["PA"]:6.1f}'
            key = (model, r['angle'])
            if key in step4:
                q, sq, u, su = step4[key]
                line += f'   {q:+8.3f}±{sq:<6.3f} {u:+7.3f}±{su:.3f}'
            print(line)

    if args.plot:
        plot(args.spec[0], args.plot, args.nbin)


if __name__ == '__main__':
    sys.exit(main())
