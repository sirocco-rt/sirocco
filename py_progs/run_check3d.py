#!/usr/bin/env python

"""
run_check3d.py -- Comprehensive run diagnostics using interactive Plotly plots.

Handles 1D, 2D, and 3D wind models.  All plots are embedded directly in the
HTML summary page.  The 3D wind plot (large) is kept in a separate file and
shown via an iframe.

Requires Plotly (pip install plotly).

Usage
-----
    run_check3d.py root [root ...]
    run_check3d.py -all          # process every .wind_save in current directory
"""

import os
import sys
import subprocess
import numpy as np
from glob import glob
from astropy.io import ascii
from astropy.table import Table

import xhtml
import plot_wind_3d

try:
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
    PLOTLY_AVAILABLE = True
except ImportError:
    PLOTLY_AVAILABLE = False
    print('Error: Plotly not found — run_check3d.py requires Plotly.')
    print('       Install with:  pip install plotly')
    sys.exit(1)


# ---------------------------------------------------------------------------
# Utility functions
# ---------------------------------------------------------------------------

def windsave2table(root):
    """Run windsave2table; fall back to versioned binary on failure."""

    def _xwindsave2table(root):
        sfiles = glob('%s*spec' % root)
        if sfiles:
            with open(sfiles[0]) as f:
                line = f.readline()
            words = line.split()
            if words[1] in ('Python', 'Sirocco'):
                xver = words[3]
                command = 'windsave2table-%s %s' % (xver, root)
            else:
                return True
        else:
            return True
        print('Trying fallback: %s' % command)
        proc = subprocess.Popen(command, shell=True,
                                stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = proc.communicate()
        if proc.returncode or len(stderr):
            print('Error: fallback also failed')
            return True
        return False

    command = 'windsave2table %s' % root
    proc = subprocess.Popen(command, shell=True,
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = proc.communicate()
    if proc.returncode:
        print('Error: windsave2table returned', proc.returncode)
        return _xwindsave2table(root)
    elif len(stderr):
        print('Error: windsave2table stderr:')
        print(stderr.decode())
        return True
    return False


def check_completion(root):
    """Parse .sig file and return a list of summary strings."""
    try:
        with open(root + '.sig') as f:
            lines = f.readlines()
    except IOError:
        print('Error: %s.sig not found' % root)
        return []

    ion_string = spec_string = ''
    restart = False
    tot = 0
    for i, line in enumerate(lines):
        if 'Finished' in line:
            if 'ionization' in line:
                ion_string = line
            elif 'spectrum' in line:
                spec_string = line
        elif 'RESTART' in line:
            restart = True
            ion_string = spec_string = ''
            word = lines[i - 1].split()
            tot += eval(word[5])

    complete_string = lines[-1]
    msgs = []
    if restart:
        msgs.append('%s was restarted; earlier run(s) took %.1f s.' % (root, tot))
        msgs.append('Times below are since the last restart.')

    ion_time = 0
    if 'COMPLETE' in complete_string:
        word = complete_string.split()
        msgs.append('%s ran to completion in %s s.' % (root, word[5]))
        try:
            word = ion_string.split()
            msgs.append('%s ionization cycles in %s s.' % (word[8], word[5]))
            ion_time = eval(word[5])
        except Exception:
            msgs.append('No ionization cycles.')
        try:
            word = spec_string.split()
            msgs.append('%s spectrum cycles in %s s.' % (word[8], eval(word[5]) - ion_time))
        except Exception:
            msgs.append('No spectrum cycles.')
    else:
        msgs.append('WARNING: %s HAS NOT COMPLETED.' % root)
        word = complete_string.split()
        try:
            msgs.append('Stopped after ~%s s in %s of %s %s cycles.' %
                        (word[5], word[8], word[10], word[11]))
        except IndexError:
            msgs.append('Could not read cycle information.')
    return msgs


def read_diag(root):
    """Extract per-cycle convergence statistics from the diag file."""
    filename = 'diag_%s/%s_00.diag' % (root, root)
    if not os.path.exists(filename):
        filename = 'diag_%s/%s_0.diag' % (root, root)

    command = "grep 'Check_convergence' %s" % filename
    proc = subprocess.Popen(command, shell=True,
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = proc.communicate()
    if len(stderr):
        return [], [], [], [], []

    converged = []
    converging = []
    t_r = []
    t_e = []
    hc = []
    ncells = 1
    for line in stdout.decode().split('\n'):
        if 'converged' in line:
            x = line.replace('(', ' ').replace(')', ' ').split()
            converged.append(eval(x[2]))
            converging.append(eval(x[6]))
            ncells = eval(x[9])
        elif 'hc(real' in line:
            x = line.split()
            t_r.append(eval(x[2]))
            t_e.append(eval(x[4]))
            hc.append(eval(x[8]))

    if t_r:
        t_r = np.array(t_r) / ncells
        t_e = np.array(t_e) / ncells
        hc  = np.array(hc)  / ncells
        ncycle = np.arange(len(converged))
        Table([ncycle, converged, converging, t_r, t_e, hc],
              names=['Ncycle', 'Converged', 'Converging',
                     'T_r_converged', 'T_e_converged', 'hc_converged']
              ).write('%s.convergence.txt' % root,
                      format='ascii.fixed_width_two_line', overwrite=True)
    else:
        t_r = t_e = hc = []
    return converged, converging, t_r, t_e, hc


def py_error(root):
    """Run py_error.py and return its output lines."""
    proc = subprocess.Popen('py_error.py %s' % root, shell=True,
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = proc.communicate()
    if len(stderr):
        return []
    return stdout.decode().split('\n')


def how_many_dimensions(filename):
    """Return 1, 2, or 3 based on whether the master table has j and/or k columns."""
    x = ascii.read(filename)
    if 'k' in x.colnames:
        return 3
    if 'j' in x.colnames:
        return 2
    return 1


# ---------------------------------------------------------------------------
# Interactive Plotly plot builders — each returns an HTML div string or None
# ---------------------------------------------------------------------------

def _smooth(x, n=21):
    """Boxcar smooth."""
    if n <= 1 or len(x) <= n:
        return x
    kernel = np.ones(n) / n
    return np.convolve(x, kernel, mode='same')


def make_convergence_div(root, converged, converging, t_r, t_e, hc,
                         include_plotlyjs=True):
    """Interactive convergence-by-cycle plot."""
    ncycle = list(range(len(converged)))
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=ncycle, y=list(t_r), mode='lines', name='t_r',
                             line=dict(dash='dash')))
    fig.add_trace(go.Scatter(x=ncycle, y=list(t_e), mode='lines', name='t_e',
                             line=dict(dash='dash')))
    fig.add_trace(go.Scatter(x=ncycle, y=list(hc), mode='lines',
                             name='heating:cooling', line=dict(dash='dash')))
    fig.add_trace(go.Scatter(x=ncycle, y=converged, mode='lines',
                             name='Converged', line=dict(width=2)))
    fig.add_trace(go.Scatter(x=ncycle, y=converging, mode='lines',
                             name='Converging', line=dict(width=2)))
    fig.update_layout(
        title='%s — convergence by cycle' % root,
        xaxis_title='Ionization cycle',
        yaxis_title='Fraction of cells in wind',
        height=420, width=750,
        legend=dict(x=0.01, y=0.99),
    )
    return fig.to_html(full_html=False, include_plotlyjs=include_plotlyjs)


def make_spec_tot_div(root, include_plotlyjs=False):
    """Interactive nuLnu vs frequency from <root>.log_spec_tot."""
    filename = root + '.log_spec_tot'
    try:
        data = ascii.read(filename)
    except IOError:
        print('Warning: %s not found; skipping spec_tot plot' % filename)
        return None

    freq = np.array(data['Freq.'], dtype=float)
    fig = go.Figure()

    for col, label, dash in [
        ('Created',  'Created (ext)',   'solid'),
        ('Emitted',  'Observed (tot)',  'solid'),
        ('Wind',     'Observed (wind)', 'dash'),
        ('HitSurf',  'Hit surface',     'dot'),
    ]:
        try:
            y = freq * _smooth(np.array(data[col], dtype=float))
            fig.add_trace(go.Scatter(x=freq, y=y, mode='lines',
                                     name=label, line=dict(dash=dash)))
        except Exception:
            pass
    try:
        y = freq * _smooth(np.array(data['WCreated'], dtype=float))
        fig.add_trace(go.Scatter(x=freq, y=y, mode='lines',
                                 name='Created (wind)', line=dict(dash='dashdot')))
    except Exception:
        pass

    for freq_val, label in [(3.29e15, 'Ly limit'), (2.47e15, 'Lyα'),
                             (1.31e16, 'HeII'), (6.91e14, 'Hβ')]:
        fig.add_vline(x=freq_val, line_width=1, line_dash='dot', line_color='grey',
                      annotation_text=label, annotation_position='top right')

    fig.update_xaxes(type='log', title='Frequency (Hz)')
    fig.update_yaxes(type='log', title='νL<sub>ν</sub>')
    fig.update_layout(title='%s — spectral energy distribution' % root,
                      height=500, width=950, legend=dict(x=0.01, y=0.99))
    return fig.to_html(full_html=False, include_plotlyjs=include_plotlyjs)


def make_spec_div(root, include_plotlyjs=False):
    """Interactive flux vs wavelength from <root>.spec; angles togglable."""
    filename = root + '.spec'
    try:
        data = ascii.read(filename)
    except IOError:
        print('Warning: %s not found; skipping angle spectra plot' % filename)
        return None

    wave = np.array(data['Lambda'], dtype=float)
    angle_cols = [c for c in data.colnames if c.startswith('A')]
    if not angle_cols:
        return None

    fig = go.Figure()
    for col in angle_cols:
        label = col.replace('P0.50', '').lstrip('A') + '°'
        flux = _smooth(np.array(data[col], dtype=float))
        fig.add_trace(go.Scatter(x=wave, y=flux, mode='lines', name=label))

    fig.update_xaxes(title='Wavelength (Å)')
    fig.update_yaxes(title='Flux')
    fig.update_layout(
        title='%s — angle spectra (click legend to toggle)' % root,
        height=500, width=950, legend=dict(x=0.01, y=0.99),
    )
    return fig.to_html(full_html=False, include_plotlyjs=include_plotlyjs)


def make_wind_2d_div(root, master_file,
                     vars=('t_e', 't_r', 'ne', 'rho'),
                     include_plotlyjs=False):
    """
    Interactive heatmap subplots for a 2D (CYLIND or RTHETA) wind model.
    Uses physical coordinates (log rho vs z) with the same masking as plot_wind.
    """
    try:
        data = ascii.read(master_file)
    except IOError:
        print('Warning: cannot read %s' % master_file)
        return None

    is_sph = 'rcen' in data.colnames
    xcol = 'rcen' if is_sph else 'xcen'
    ycol = 'thetacen' if is_sph else 'zcen'
    # Fall back to x/z for plain 2D tables
    if xcol not in data.colnames:
        xcol, ycol = 'x', 'z'

    ii = np.array(data['i'])
    jj = np.array(data['j'])
    ndim = int(np.max(ii)) + 1
    mdim = int(np.max(jj)) + 1
    inwind = np.array(data['inwind']).reshape(ndim, mdim)

    x_raw = np.array(data[xcol]).reshape(ndim, mdim)
    y_raw = np.array(data[ycol]).reshape(ndim, mdim)

    # Log-scale the radial axis
    xpos = x_raw[x_raw > 1]
    xlogmin = np.log10(np.min(xpos) / 10) if len(xpos) else 0.0
    x_plot = np.where(x_raw > 1, np.log10(x_raw), xlogmin)

    nvars = len(vars)
    ncols = min(2, nvars)
    nrows = (nvars + 1) // 2
    fig = make_subplots(rows=nrows, cols=ncols,
                        subplot_titles=list(vars),
                        horizontal_spacing=0.12,
                        vertical_spacing=0.15)

    for idx, var in enumerate(vars):
        row = idx // ncols + 1
        col = idx % ncols + 1
        if var not in data.colnames:
            continue
        z = np.array(data[var], dtype=float).reshape(ndim, mdim)
        mask = inwind < 0
        z[mask] = np.nan
        with np.errstate(divide='ignore', invalid='ignore'):
            zlog = np.where(z > 0, np.log10(z), np.nan)
        fig.add_trace(
            go.Heatmap(z=zlog.T, x=x_plot[:, 0], y=y_raw[0, :],
                       colorscale='Viridis', showscale=True,
                       name='log(%s)' % var,
                       hovertemplate='log(ρ)=%{x:.2g}<br>z=%{y:.2g}<br>'
                                     + 'log(%s)' % var + '=%{z:.3f}<extra></extra>'),
            row=row, col=col,
        )
        fig.update_xaxes(title_text='log(ρ)' if not is_sph else 'log(r)',
                         row=row, col=col)
        fig.update_yaxes(title_text='z' if not is_sph else 'θ (°)',
                         row=row, col=col)

    fig.update_layout(
        title='%s — wind structure' % root,
        height=420 * nrows, width=950,
    )
    return fig.to_html(full_html=False, include_plotlyjs=include_plotlyjs)


def make_wind_1d_div(root, master_file,
                     vars=('t_e', 't_r', 'ne', 'rho'),
                     include_plotlyjs=False):
    """
    Interactive line plots for a 1D (spherical/shell) wind model.
    """
    try:
        data = ascii.read(master_file)
    except IOError:
        print('Warning: cannot read %s' % master_file)
        return None

    rcol = 'rcen' if 'rcen' in data.colnames else 'xcen'
    if rcol not in data.colnames:
        rcol = 'x'

    inwind = np.array(data['inwind'])
    r = np.array(data[rcol], dtype=float)
    mask = inwind >= 0

    fig = go.Figure()
    for var in vars:
        if var not in data.colnames:
            continue
        y = np.array(data[var], dtype=float)
        y[~mask] = np.nan
        with np.errstate(divide='ignore', invalid='ignore'):
            ylog = np.where(y > 0, np.log10(y), np.nan)
        fig.add_trace(go.Scatter(x=r[mask], y=ylog[mask], mode='lines',
                                 name='log(%s)' % var))

    fig.update_xaxes(type='log', title='r (cm)')
    fig.update_yaxes(title='log(value)')
    fig.update_layout(title='%s — wind structure' % root,
                      height=450, width=750, legend=dict(x=0.01, y=0.99))
    return fig.to_html(full_html=False, include_plotlyjs=include_plotlyjs)


def make_wind_3d_file(root, master_file, outdir):
    """Write the 3D wind interactive HTML; returns the file path or None."""
    outfile = os.path.join(outdir, '%s_wind_interactive.html' % root)
    return plot_wind_3d.doit(master_file, var=['t_e', 't_r', 'ne', 'rho'],
                              scale='log', outfile=outfile)


# ---------------------------------------------------------------------------
# HTML summary
# ---------------------------------------------------------------------------

def make_html(root, complete_message, errors, divs, wind_3d_path=None):
    """
    Write the HTML summary.

    divs : ordered list of (heading, html_div_string) pairs; None entries skipped.
    wind_3d_path : relative path to 3D wind iframe file, or None.
    """
    string = xhtml.begin('%s: Run summary' % root)
    string += xhtml.paragraph('Summary for %s' % root)
    string += xhtml.hline()

    string += xhtml.h2('Run status')
    for msg in complete_message:
        string += xhtml.paragraph(msg)

    for heading, div in divs:
        if div is None:
            continue
        string += xhtml.h2(heading)
        string += div + '\n'

    if wind_3d_path:
        string += xhtml.h2('3D wind structure (interactive)')
        string += xhtml.paragraph(
            'Use the sliders to explore phi and z/theta slices.')
        string += ('<iframe src="%s" width="100%%" height="720px" '
                   'frameborder="0" style="border:1px solid #ccc;"></iframe>\n'
                   % wind_3d_path)

    string += xhtml.h2('Error summary')
    string += xhtml.preformat(errors)
    string += xhtml.end()

    htmlfile = root + '.html'
    with open(htmlfile, 'w') as f:
        f.write(string)
    print('Written: %s' % htmlfile)
    return htmlfile


# ---------------------------------------------------------------------------
# Main driver
# ---------------------------------------------------------------------------

def doit(root='test'):
    """Run diagnostics for a single model root."""

    print('\nEvaluating %s\n' % root)

    if windsave2table(root):
        print('Exiting: windsave2table failed for %s' % root)
        return

    complete_message = check_completion(root)
    for msg in complete_message:
        print(msg)

    master_file = '%s.master.txt' % root
    if not os.path.exists(master_file):
        master_file = master_file.replace('master', '0.master')
        if os.path.exists(master_file):
            print('Warning: multiple domains detected; reporting domain 0')
        else:
            print('Error: cannot find master file for %s' % root)
            return

    xdim = how_many_dimensions(master_file)
    print('Grid dimensionality: %dD' % xdim)

    # ----------------------------------------------------------------
    # Convergence (first div — embeds Plotly.js for whole page)
    # ----------------------------------------------------------------
    converged, converging, t_r, t_e, hc = read_diag(root)
    if len(converged) > 1:
        conv_div = make_convergence_div(root, converged, converging, t_r, t_e, hc,
                                        include_plotlyjs=True)
    else:
        print('Not enough cycles to plot convergence history')
        conv_div = None

    # ----------------------------------------------------------------
    # Spectral plots (reuse already-loaded Plotly.js)
    # ----------------------------------------------------------------
    spec_tot_div = make_spec_tot_div(root, include_plotlyjs=(conv_div is None))
    spec_div = make_spec_div(root,
                             include_plotlyjs=(conv_div is None and
                                              spec_tot_div is None))

    # ----------------------------------------------------------------
    # Wind structure plots
    # ----------------------------------------------------------------
    wind_3d_path = None
    wind_div = None

    if xdim == 3:
        interactive_dir = './diag_%s/interactive' % root
        os.makedirs(interactive_dir, exist_ok=True)
        wind_path = make_wind_3d_file(root, master_file, interactive_dir)
        if wind_path:
            wind_3d_path = os.path.relpath(wind_path)
    elif xdim == 2:
        wind_div = make_wind_2d_div(root, master_file, include_plotlyjs=False)
    else:
        wind_div = make_wind_1d_div(root, master_file, include_plotlyjs=False)

    # ----------------------------------------------------------------
    # HTML summary
    # ----------------------------------------------------------------
    errors = py_error(root)

    divs = [
        ('Convergence by cycle',              conv_div),
        ('Spectral energy distribution',      spec_tot_div),
        ('Angle spectra',                     spec_div),
        ('Wind structure',                    wind_div),
    ]

    make_html(root, complete_message, errors, divs, wind_3d_path=wind_3d_path)


# ---------------------------------------------------------------------------
# Command-line interface
# ---------------------------------------------------------------------------

def steer(argv):
    xall = False
    files = []
    i = 1
    while i < len(argv):
        if argv[i] == '-h':
            print(__doc__)
            return
        elif argv[i] == '-all':
            xall = True
            break
        elif argv[i].endswith('.out.pf'):
            pass
        elif argv[i].endswith('.pf'):
            files.append(argv[i].replace('.pf', ''))
        elif argv[i].startswith('-'):
            print('Unknown option: %s' % argv[i])
            return
        else:
            files.append(argv[i])
        i += 1

    if xall:
        files = [f.replace('.wind_save', '') for f in glob('*.wind_save')]

    for root in files:
        if os.path.isfile(root + '.wind_save'):
            doit(root)
        else:
            print('Warning: %s.wind_save not found; skipping' % root)


if __name__ == '__main__':
    steer(sys.argv)
