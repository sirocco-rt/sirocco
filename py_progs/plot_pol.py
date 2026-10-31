#!/usr/bin/env python
'''
Plot the spectrum and its linear polarization for each observer of a
Sirocco run made with polarization switched on.


Command line usage (if any):

    usage::

        plot_pol.py [-h] [-wmin WMIN] [-wmax WMAX] [-emin EMIN] [-emax EMAX]
                    [-nbin NBIN] [-max_obs N] [-out FILE] [-linear]
                    [-separate] [-xray | -optical] [-table]
                    rootname

        rootname    the root name of the run, or the spectrum file itself
                    (rootname.spec is read by default, rootname.log_spec
                    for X-rays)
        -wmin       shortest wavelength (A) to plot (default: all)
        -wmax       longest wavelength (A) to plot (default: all)
        -emin       lowest energy (keV) to plot; sets -wmax
        -emax       highest energy (keV) to plot; sets -wmin
        -nbin       combine this many wavelength bins for the polarization
                    (default: about 100 bins across the plot)
        -max_obs    plot at most this many inclinations, spread evenly over
                    those in the file (default 4; 0 for all)
        -out        name of the plot (default rootname_pol.png)
        -linear     plot the flux on a linear axis (default logarithmic)
        -separate   one figure for each inclination, e.g. rootname_pol_A80.png,
                    instead of all of them in one
        -xray       plot against energy, as for X-rays (see below)
        -optical    plot against wavelength (by default the choice is made
                    from the wavelengths in the spectrum)
        -h          print this documentation
        -table      make no plot; print q, u, P and PA integrated over
                    -wmin to -wmax for every observer, including those
                    selected by the number of scatters

Description:

    When Sirocco is run with @Spectrum.polarization yes (which needs -d),
    the spectrum files have two extra columns for each observer, e.g.
    A80P0.50_Q and A80P0.50_U.  Q and U are in the same units as the
    spectrum and are referred to North on the sky, the projected +z axis
    (IAU convention), so Q > 0 means the electric vector is parallel to
    the projected symmetry axis of the model.

    For each observer the plot has three panels: the spectrum F_lambda
    (logarithmic unless -linear),
    the degree of polarization P = sqrt(Q^2 + U^2)/I in per cent, and the
    position angle PA = 0.5 atan2(U, Q) of the electric vector, measured
    from the projected symmetry axis (PA = 0 or 180: E along the axis;
    PA = 90: E parallel to the disk).  The polarization is much noisier
    than the spectrum, so it is calculated in coarse bins, each combining
    nbin bins of the spectrum.  P and PA integrated over the plotted band
    are given in the top panel.

    Sirocco's models are axisymmetric, so U should be zero to within the
    noise.  The rms of u = U/I about zero is a measure of the noise in each
    coarse bin.  It is calculated over 11 neighbouring bins, since the noise
    varies with wavelength, and is drawn on the P panel, since noise biases
    P upwards by about sigma when P is small.  The value in the top panel
    is the rms over the whole band.  Position angles are grey where
    P < 3 sigma and are not meaningful there.

    If the whole plotted range is shortward of 124 A (above 0.1 keV), or
    with -xray, the plot is made as is usual for X-rays: against energy
    in keV on a logarithmic axis, with the flux as E F_E (= nu F_nu).
    Then, unless a spectrum file is named, rootname.log_spec is read
    rather than rootname.spec, since Sirocco's .spec files have bins of
    equal width in frequency, which are very coarse at the low-energy
    end of a broad X-ray band, while the .log_spec files have bins of
    equal width in log frequency.

    If there are several spectra at the same inclination, selected by the
    number of scatters (Spectrum.select_specific_no_of_scatters_in_spectra,
    e.g. A80P0.50_sc:0, _sc:1, _sc:>1), they are drawn together: the
    spectrum of all the photons in black and the selections in colour in
    the spectrum panel, and the P of each selection in an extra panel below
    the P of all the photons (the P of the scattered light can be much
    larger).  The position angle, the noise and the band values are for all
    the photons; the band values for the selections are printed.

Primary routines:

    steer
        processes the command line
    doit
        makes the plot (plot_pol.doit(root) is the routine to call from
        a script, e.g. run_check.py)
    has_pol
        says whether a spectrum file contains polarization
    band_pol
        band-integrated q, u, P and PA for each observer, summing L_nu dnu
        over the band (also used by the polarization tests in
        py_progs/Tests/polarization)

Notes:

    Requires numpy and matplotlib.

History:

    261016 ksl Coding begun, with Claude Code, from py_progs/polarization/c_spectra.py

    261018 ksl Added band_pol and -table from c_spectra.py, which this replaces;
    spectra at the same inclination drawn together; -separate

    261019 ksl Plots against energy for X-rays (-xray, -optical, -emin, -emax)

'''

import os

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator, NullFormatter


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
        raise ValueError('%s: no column header found' % fname)
    data = np.array(rows)
    cols = {name: data[:, i] for i, name in enumerate(names)}
    return units, cols['Freq.'], cols['Lambda'], cols


def observers(cols):
    '''
    The names of the observer spectra that have Q and U columns.  If some
    are spectra of all the photons and others are selected by the number
    of scatters (_sc:...), only the former are returned.
    '''
    names = [name[:-2] for name in cols if name.endswith('_Q') and name[:-2] in cols]
    every = [name for name in names if '_sc:' not in name]
    return every if every else names


def has_pol(fname):
    '''True if the spectrum file has Q and U columns.'''
    try:
        with open(fname) as f:
            for line in f:
                if line.split()[:1] == ['Freq.']:
                    return '_Q' in line
    except (IOError, IndexError):
        pass
    return False


def lnu_dnu(units, freq, logbins):
    '''The factor that turns a column in these units into L_nu dnu per bin (up to a constant).'''
    if units == 'flambda':
        x = 1.0 / freq ** 2
    else:
        x = np.ones_like(freq)
    if logbins:
        x = x * freq
    return x


def groups(cols):
    '''
    The observer spectra with Q and U columns, grouped by observer: a list of
    (main, [names]), where names are the spectra of all the photons and of
    those selected by the number of scatters (e.g. A80P0.50, A80P0.50_sc:0,
    A80P0.50_sc:1), and main is the one of all the photons if there is one.
    '''
    out = {}
    for name in cols:
        if name.endswith('_Q') and name[:-2] in cols:
            out.setdefault(name[:-2].split('_')[0], []).append(name[:-2])
    return [(key if key in names else names[0], names) for key, names in out.items()]


def label(name):
    '''A short label for an observer column, e.g. A80P0.50 -> 80 deg'''
    angle = name[1:].split('P')[0]
    extra = name.split('P0.50')[-1] if 'P0.50' in name else ''
    return '%s$^\\circ$%s' % (angle, (' ' + extra.strip('_')) if extra else '')


def selection(name):
    '''Which photons a column contains: all, or e.g. sc:1 for those scattered once'''
    return name.split('_', 1)[1] if '_' in name else 'all'


HC_KEV = 12.398420              # hc in keV A
XRAY_WMAX = 124.0               # A; a spectrum entirely shortward of this (0.1 keV) is plotted as X-rays


def is_xray(wave, wmin=0, wmax=0):
    '''True if the wavelengths to be plotted are all shortward of XRAY_WMAX'''
    longest = wmax if wmax > 0 else np.max(wave)
    return longest <= XRAY_WMAX


def doit(root, wmin=0, wmax=0, nbin=0, max_obs=4, outfile='', logflux=True, separate=False, only='', mode=''):
    '''
    Plot the spectrum and q, u for the observers in root.spec.

    Returns the name of the plot, or 'None' if there is no polarization.
    With separate=True, each inclination gets a figure of its own, named
    after outfile with the observer added (e.g. root_pol_A80.png), and the
    list of their names is returned.  only restricts the plot to one
    observer (e.g. A80P0.50).  mode is 'xray' (against energy, from
    root.log_spec), 'optical' (against wavelength, from root.spec), or
    '' to choose from the wavelengths in the spectrum.
    '''
    named = os.path.isfile(root) and not os.path.isdir(root) and '.' in os.path.basename(root)
    fname = root if named else root + '.spec'
    if not has_pol(fname):
        print('plot_pol: %s has no Q and U columns' % fname)
        return 'None'
    if outfile == '':
        outfile = os.path.basename(fname).split('.')[0] + '_pol.png'

    units, freq, wave, cols = read_spec(fname)
    if mode == '':
        mode = 'xray' if is_xray(wave, wmin, wmax) else 'optical'
    xray = mode == 'xray'
    if xray and not named and has_pol(root + '.log_spec'):
        fname = root + '.log_spec'
        units, freq, wave, cols = read_spec(fname)
    fac = lnu_dnu(units, freq, 'log_spec' in os.path.basename(fname))
    # the quantity plotted along x, and the flux in the top panel (nu F_nu for X-rays)
    xval = HC_KEV / wave if xray else wave
    if xray:
        fscale = wave if units == 'flambda' else freq
    else:
        fscale = np.ones_like(wave)

    sel = np.ones(len(wave), dtype=bool)
    if wmin > 0:
        sel &= wave >= wmin
    if wmax > 0:
        sel &= wave <= wmax
    idx = np.where(sel)[0]
    if len(idx) == 0:
        print('plot_pol: no wavelengths between %g and %g' % (wmin, wmax))
        return 'None'
    if nbin <= 0:
        nbin = max(1, len(idx) // 100)
    nkeep = (len(idx) // nbin) * nbin
    idx = idx[:nkeep]

    def comb(x):
        return x[idx].reshape(-1, nbin).sum(axis=1)

    obs = groups(cols)
    if max_obs > 0 and len(obs) > max_obs:
        obs = [obs[i] for i in np.round(np.linspace(0, len(obs) - 1, max_obs)).astype(int)]
    if separate:
        base = os.path.splitext(outfile)[0]
        return [doit(root, wmin, wmax, nbin, 0, '%s_%s.png' % (base, name.split('_')[0].replace('P0.50', '')),
                     logflux, only=name.split('_')[0], mode=mode) for name, variants in obs]
    if only:
        obs = [one for one in obs if one[0].split('_')[0] == only]

    nobs = len(obs)
    nsel = max(len(variants) for name, variants in obs)
    fig = plt.figure(figsize=(9, (4.6 if nsel == 1 else 5.8) * nobs))
    outer = fig.add_gridspec(nobs, 1, hspace=0.22)
    nsmooth = max(1, nbin // 10)
    wb = xval[idx].reshape(-1, nbin).mean(axis=1)

    for k, (name, variants) in enumerate(obs):
        I = comb(cols[name] * fac)
        Q = comb(cols[name + '_Q'] * fac)
        U = comb(cols[name + '_U'] * fac)
        ok = I > 0
        with np.errstate(invalid='ignore', divide='ignore'):
            q = np.where(ok, 100 * Q / I, np.nan)
            u = np.where(ok, 100 * U / I, np.nan)
        p = np.hypot(q, u)
        qband = 100 * Q.sum() / I.sum()
        uband = 100 * U.sum() / I.sum()
        pband = np.hypot(qband, uband)
        paband = np.degrees(0.5 * np.arctan2(uband, qband)) % 180
        # the models are axisymmetric, so u is noise; its rms is the noise in q, u (and P) per bin,
        # overall (sigma) and locally, over 11 neighbouring bins (sig_loc), since the noise varies
        sigma = np.sqrt(np.nanmean(u ** 2))
        u2 = np.nan_to_num(u ** 2)
        nwin = np.convolve(np.isfinite(u).astype(float), np.ones(11), mode='same')
        sig_loc = np.sqrt(np.convolve(u2, np.ones(11), mode='same') / np.maximum(nwin, 1))
        # position angle, wrapped into the 180 deg range centred on the band value
        pa = np.degrees(0.5 * np.arctan2(u, q))
        pa = (pa - paband + 90) % 180 + paband - 90
        good = p > 3 * sig_loc

        # the spectrum of all the photons (or the first selection) in black, with the selections by number
        # of scatters, if any, in colour on the same axes, and their P in a panel of their own, since P of
        # the scattered light can be much larger than P of the total
        others = [v for v in variants if v != name]
        if others:
            gs = outer[k].subgridspec(4, 1, height_ratios=[1.6, 1, 1, 0.8], hspace=0)
        else:
            gs = outer[k].subgridspec(3, 1, height_ratios=[1.6, 1, 0.8], hspace=0)
        ax1 = fig.add_subplot(gs[0])
        ax2 = fig.add_subplot(gs[1], sharex=ax1)
        ax2s = fig.add_subplot(gs[2], sharex=ax1) if others else None
        ax3 = fig.add_subplot(gs[-1], sharex=ax1)

        flux = np.convolve(cols[name] * fscale, np.ones(nsmooth) / nsmooth, mode='same')
        ax1.plot(xval[idx], flux[idx], color='0.65', lw=0.5)
        ax1.plot(wb, comb(cols[name] * fscale) / nbin, color='k', lw=1.0, label=selection(name))
        for j, v in enumerate(others):
            ax1.plot(wb, comb(cols[v] * fscale) / nbin, color='C%d' % (j + 1), lw=0.9, label=selection(v))
        if others:
            ax1.legend(loc='upper left', fontsize=7, ncol=len(variants))
        if logflux:
            fb = np.concatenate([comb(cols[v] * fscale) / nbin for v in variants])
            top = np.nanmax(fb[fb > 0]) if np.any(fb > 0) else 1.0
            ax1.set_yscale('log')
            ax1.set_ylim(max(np.nanmin(fb[fb > 0]) / 2, top / 1e4), top * 2)
        else:
            ax1.set_ylim(bottom=0)
            ax1.ticklabel_format(axis='y', style='sci', scilimits=(-2, 3), useMathText=True)
            ax1.yaxis.set_major_locator(MaxNLocator(nbins=5, prune='lower'))
        if xray:
            ax1.set_ylabel('E F$_E$ ($\\nu$ F$_\\nu$)')
            ax1.set_xscale('log')
        else:
            ax1.set_ylabel('F$_\\lambda$' if units == 'flambda' else 'flux')
        ax1.text(0.98, 0.92, label(name.split('_')[0]), transform=ax1.transAxes, ha='right', va='top')
        ax1.text(0.02, 0.06, 'band: P = %.2f%%, PA = %.0f$^\\circ$;  $\\sigma$(bin) = %.2f%%' % (pband, paband, sigma),
                 transform=ax1.transAxes, ha='left', va='bottom', fontsize=9,
                 bbox=dict(facecolor='white', edgecolor='none', alpha=0.8))

        ax2.step(wb, p, where='mid', color='C0' if not others else 'k')
        ax2.plot(wb, sig_loc, color='C3', lw=0.8, ls=':', label='$\\sigma$ (noise)')
        ax2.set_ylim(bottom=0)
        ax2.set_ylabel('P (%)' if not others else 'P, all (%)')
        ax2.yaxis.set_major_locator(MaxNLocator(nbins=4, prune='upper'))
        if k == 0:
            ax2.legend(loc='upper right', fontsize=8)

        if others:
            for j, v in enumerate(others):
                Iv = comb(cols[v] * fac)
                with np.errstate(invalid='ignore', divide='ignore'):
                    pv = np.where(Iv > 0, 100 * np.hypot(comb(cols[v + '_Q'] * fac), comb(cols[v + '_U'] * fac)) / Iv,
                                  np.nan)
                ax2s.step(wb, pv, where='mid', color='C%d' % (j + 1), lw=0.9)
            ax2s.set_ylim(bottom=0)
            ax2s.set_ylabel('P, selected (%)')
            ax2s.yaxis.set_major_locator(MaxNLocator(nbins=4, prune='upper'))

        ax3.plot(wb[good], pa[good], '.', color='C0', ms=4)
        ax3.plot(wb[~good], pa[~good], '.', color='0.75', ms=3)
        for x, txt in ((0, 'axis'), (90, 'disk'), (180, 'axis')):
            if paband - 90 <= x <= paband + 90:
                ax3.axhline(x, color='0.6', lw=0.5, ls='--')
                ax3.text(1.005, x, txt, transform=ax3.get_yaxis_transform(), fontsize=7, va='center', color='0.4')
        ax3.set_ylim(paband - 90, paband + 90)
        # ticks well inside the panel, so they do not run into the P labels above
        pc = 45 * np.round(paband / 45) if abs(paband - 45 * np.round(paband / 45)) < 10 else np.round(paband)
        ax3.set_yticks([pc - 45, pc, pc + 45])
        ax3.set_xlim(xval[idx].min(), xval[idx].max())      # the actual range; shared by this observer's panels
        if xray:
            # plain numbers (0.5, 1, 2, ... keV) rather than powers of ten on the log axis
            lo, hi = ax3.get_xlim()
            ticks = [t for t in (0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 100, 200, 500) if lo <= t <= hi]
            ax3.set_xticks(ticks)
            ax3.set_xticklabels(['%g' % t for t in ticks])
            ax3.xaxis.set_minor_formatter(NullFormatter())
        ax3.set_ylabel('PA ($^\\circ$)')

        for a in (ax1, ax2, ax2s):
            if a is not None:
                plt.setp(a.get_xticklabels(), visible=False)
        if k < nobs - 1:
            plt.setp(ax3.get_xticklabels(), visible=False)
        print('%-16s P = %6.3f%%  PA = %5.1f  (q = %+7.3f%%  u = %+7.3f%%)  sigma(bin) = %.3f%%'
              % (name, pband, paband, qband, uband, sigma))
        for v in others:
            Iv, Qv, Uv = (comb(cols[v + t] * fac).sum() for t in ('', '_Q', '_U'))
            if Iv > 0:
                print('  %-14s P = %6.3f%%  (q = %+7.3f%%  u = %+7.3f%%)'
                      % (selection(v), 100 * np.hypot(Qv, Uv) / Iv, 100 * Qv / Iv, 100 * Uv / Iv))

    ax3.set_xlabel('Energy (keV)' if xray else 'Wavelength ($\\AA$)')
    fig.axes[0].set_title('%s: spectrum, polarization and position angle (%d bins per point; grey: P < 3$\\sigma$)'
                          % (os.path.basename(fname), nbin), fontsize=10)
    fig.savefig(outfile, dpi=100, bbox_inches='tight')
    plt.close(fig)
    return outfile


def band_pol(fname, wmin=0, wmax=0):
    '''
    Band-integrated polarization for every observer in a spectrum file,
    including those selected by the number of scatters.

    Returns a list of dicts with name, angle, I, q, u (in %), P (%) and PA (deg).
    '''
    units, freq, wave, cols = read_spec(fname)
    sel = np.ones(len(freq), dtype=bool)
    if wmin > 0:
        sel &= wave >= wmin
    if wmax > 0:
        sel &= wave <= wmax
    fac = lnu_dnu(units, freq, 'log_spec' in os.path.basename(fname))[sel]

    out = []
    for name in [one[:-2] for one in cols if one.endswith('_Q') and one[:-2] in cols]:
        I = np.sum(cols[name][sel] * fac)
        Q = np.sum(cols[name + '_Q'][sel] * fac)
        U = np.sum(cols[name + '_U'][sel] * fac)
        q = 100 * Q / I if I > 0 else np.nan
        u = 100 * U / I if I > 0 else np.nan
        angle = float(name[1:].split('P')[0])
        out.append(dict(name=name, angle=angle, I=I, q=q, u=u,
                        P=np.hypot(q, u), PA=np.degrees(0.5 * np.arctan2(u, q)) % 180))
    return out


def table(root, wmin=0, wmax=0):
    '''Print the band-integrated q, u, P and PA of every observer in root.spec'''
    fname = root if os.path.isfile(root) and not os.path.isdir(root) and '.' in os.path.basename(root) \
        else root + '.spec'
    if not has_pol(fname):
        print('plot_pol: %s has no Q and U columns' % fname)
        return
    print('%s' % fname)
    print('  %-16s %8s %8s %7s %6s' % ('observer', 'q (%)', 'u (%)', 'P (%)', 'PA'))
    for r in band_pol(fname, wmin, wmax):
        print('  %-16s %+8.3f %+8.3f %7.3f %6.1f' % (r['name'], r['q'], r['u'], r['P'], r['PA']))


def steer(argv):
    '''
    Process the command line
    '''

    root = ''
    wmin = 0
    wmax = 0
    emin = 0
    emax = 0
    nbin = 0
    max_obs = 4
    outfile = ''
    logflux = True
    separate = False
    mode = ''
    xtable = False

    i = 1
    while i < len(argv):
        if argv[i] == '-h':
            print(__doc__)
            return
        elif argv[i] == '-wmin':
            i += 1
            wmin = float(argv[i])
        elif argv[i] == '-wmax':
            i += 1
            wmax = float(argv[i])
        elif argv[i] == '-emin':
            i += 1
            emin = float(argv[i])
        elif argv[i] == '-emax':
            i += 1
            emax = float(argv[i])
        elif argv[i] == '-nbin':
            i += 1
            nbin = int(argv[i])
        elif argv[i] == '-max_obs':
            i += 1
            max_obs = int(argv[i])
        elif argv[i] == '-out':
            i += 1
            outfile = argv[i]
        elif argv[i] == '-linear':
            logflux = False
        elif argv[i] == '-separate':
            separate = True
        elif argv[i] == '-xray':
            mode = 'xray'
        elif argv[i] == '-optical':
            mode = 'optical'
        elif argv[i] == '-table':
            xtable = True
        elif argv[i][0] == '-':
            print('Error: Unknown switch ---  %s' % argv[i])
            return
        else:
            root = argv[i]
        i += 1

    if root == '':
        print('Error: no rootname given')
        return

    if emin > 0:
        wmax = HC_KEV / emin
    if emax > 0:
        wmin = HC_KEV / emax

    if xtable:
        table(root, wmin, wmax)
    else:
        print(doit(root, wmin, wmax, nbin, max_obs, outfile, logflux, separate, mode=mode))


# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1:
        steer(sys.argv)
    else:
        print(__doc__)
