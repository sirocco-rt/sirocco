#!/usr/bin/env python

"""
run_check3d.py -- Comprehensive run diagnostics using interactive Plotly plots.

Intended as a future replacement for run_check.py.  Handles 1D, 2D, and 3D
wind models.  All plots are embedded directly as interactive Plotly figures
in the HTML summary; no static matplotlib images are produced.

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
# Text strings — identical to run_check.py
# ---------------------------------------------------------------------------

convergence_message = '''
The plot shows the fraction of cells which satisfy various convergence criteria, comparing the radiation
temperature, t_r, the electron temperature, t_e between cycles, or the balance between heating and cooling
in a given cycle.  If the fractional variation in, for example, the electron temperature has changed by
less than 5%, then the electron temperature criterion has been passed.  If a cell passes all three tests,
it is said to have "converged".  The plot also indicates the number of cells for which the electron
temperature is evolving to a higher or lower value, rather than oscillating up and down in electron
temperature.  Such cells are said to be "converging".
'''

convergence2 = '''
The plot below indicates which cells have converged.  The values given for each cell indicate the number
of convergence tests which have failed, so a value of 0 indicates that a cell is "converged".
'''

error_description = '''
Python accumulates a fairly large number of errors and warnings.  Most are benign, especially if they
only occur a few times, but one should beware if an error message occurs many times, or if messages
that look unusual start to appear.  A summary of the errors for this run of the program is shown below.
More information about where the errors occurred can be found in the diag files directory.  This summary
presents the total number of times an error message of a particular type was generated; many times the
same error message will occur in each of the threads.
'''

tot_plot_description = '''
This plot is based on the .spec_tot file.  It shows the photons that were created and emitted in the
last ionization cycle.  The photons created by a central object or disk (the external sources) are
plotted separately from the wind (the internal source).  For a simple atom model, the wind photons are
generated based on physical conditions in the wind at the beginning of the cycle.  For a macro-atom
model, the wind photons are generated as photons pass through the wind and interact with macro-atoms.
The observed flux from photons that escape to infinity and a subset of this, the emitted flux of photons
arising in the wind are also shown.  Finally, the flux of the photons that hit the disk or star is also
shown (as they would be seen at a distance of 100 pc).
'''

spec_plot_description = '''
This plot is from the .spec file, and shows the expected spectra (at 100 pc) as a function of the
various inclination angles requested.  Note that the y-axis scale varies for different inclination angles.
'''

ne_plot_description = '''
This plot shows the electron density (ne) throughout the wind.  The electron density traces the density
structure of the outflow and is directly related to the optical depth of the wind.
'''


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
        msgs.append('%s was restarted. The earlier run(s) took %.1f s.' % (root, tot))
        msgs.append('Here we report times and progress since the (last) restart.')

    ion_time = 0
    if 'COMPLETE' in complete_string:
        word = complete_string.split()
        msgs.append('%s ran to completion in %s s.' % (root, word[5]))
        try:
            word = ion_string.split()
            msgs.append('%s ionization cycles were completed in %s s.' % (word[8], word[5]))
            ion_time = eval(word[5])
        except Exception:
            msgs.append('There were no ionization cycles for this run.')
        try:
            word = spec_string.split()
            msgs.append('%s spectrum cycles were completed in %s s.' %
                        (word[8], eval(word[5]) - ion_time))
        except Exception:
            msgs.append('There were no spectrum cycles for this run.')
    else:
        msgs.append('WARNING: RUN %s HAS NOT COMPLETED SUCCESSFULLY.' % root)
        word = complete_string.split()
        try:
            msgs.append('If not running, it stopped after about %s s in %s of %s %s cycles.' %
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
# Interactive Plotly builders — each returns an HTML div string or None
# ---------------------------------------------------------------------------

def _smooth(x, n=21):
    if n <= 1 or len(x) <= n:
        return x
    kernel = np.ones(n) / n
    return np.convolve(x, kernel, mode='same')


def _first_true(include_plotlyjs):
    """Return the include_plotlyjs value for the first div that isn't None."""
    return include_plotlyjs


def make_convergence_div(root, converged, converging, t_r, t_e, hc,
                         include_plotlyjs=True):
    """Convergence fractions vs ionization cycle."""
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
        xaxis_title='Ionization Cycle',
        yaxis_title='Fraction of Cells in Wind',
        height=420, width=750,
        legend=dict(x=0.01, y=0.99),
    )
    return fig.to_html(full_html=False, include_plotlyjs=include_plotlyjs)


def make_wind_3d_div(master_file, var, scale='log', inwind='',
                     include_plotlyjs=False):
    """
    Inline interactive 3D wind plot for a single variable using plot_wind_3d.
    Returns an HTML div string or None.
    """
    fig = plot_wind_3d._build_figure(master_file, var, scale=scale, inwind=inwind)
    if fig is None:
        return None
    return fig.to_html(full_html=False, include_plotlyjs=include_plotlyjs)


def _wind_heatmap_2d(data, xcol, ycol, var, scale, inwind_arr, ndim, mdim,
                     xlabel, ylabel):
    """Return (x_plot, y_plot, zlog, title) arrays for a single 2D wind variable."""
    x_raw = np.array(data[xcol]).reshape(ndim, mdim)
    y_raw = np.array(data[ycol]).reshape(ndim, mdim)
    inwind = inwind_arr.reshape(ndim, mdim)

    xpos = x_raw[x_raw > 1]
    xlogmin = np.log10(np.min(xpos) / 10) if len(xpos) else 0.0
    x_plot = np.where(x_raw > 1, np.log10(x_raw), xlogmin)

    z = np.array(data[var], dtype=float).reshape(ndim, mdim)
    if inwind != 'all':
        mask = inwind_arr.reshape(ndim, mdim) < 0
        z[mask] = np.nan

    if scale == 'log':
        with np.errstate(divide='ignore', invalid='ignore'):
            z = np.where(z > 0, np.log10(z), np.nan)
        title = 'log(%s)' % var
    else:
        title = var

    return x_plot[:, 0], y_raw[0, :], z.T, title, xlabel, ylabel


def make_wind_2d_div(root, master_file, vars=('t_e',), scale='log', inwind='',
                     include_plotlyjs=False):
    """Plotly heatmap subplots for a 2D wind model; one panel per variable."""
    try:
        data = ascii.read(master_file)
    except IOError:
        print('Warning: cannot read %s' % master_file)
        return None

    is_sph = 'rcen' in data.colnames
    xcol = 'rcen' if is_sph else ('xcen' if 'xcen' in data.colnames else 'x')
    ycol = 'thetacen' if is_sph else ('zcen' if 'zcen' in data.colnames else 'z')
    xlabel = 'log(r)' if is_sph else 'log(ρ)'
    ylabel = 'θ (°)' if is_sph else 'z'

    ii = np.array(data['i'])
    jj = np.array(data['j'])
    ndim = int(np.max(ii)) + 1
    mdim = int(np.max(jj)) + 1
    inwind_arr = np.array(data['inwind'])

    vars = [v for v in vars if v in data.colnames]
    if not vars:
        return None

    ncols = min(2, len(vars))
    nrows = (len(vars) + 1) // 2
    fig = make_subplots(rows=nrows, cols=ncols,
                        subplot_titles=list(vars),
                        horizontal_spacing=0.12,
                        vertical_spacing=0.15)

    for idx, var in enumerate(vars):
        row = idx // ncols + 1
        col = idx % ncols + 1
        var_scale = 'lin' if var == 'converge' else scale
        x, y, z, title, xl, yl = _wind_heatmap_2d(
            data, xcol, ycol, var, var_scale,
            inwind_arr if var != 'converge' else np.zeros_like(inwind_arr),
            ndim, mdim, xlabel, ylabel)
        fig.add_trace(
            go.Heatmap(z=z, x=x, y=y, colorscale='Viridis', showscale=True,
                       name=title,
                       hovertemplate=xl + '=%{x:.2g}<br>'
                                     + yl + '=%{y:.2g}<br>'
                                     + title + '=%{z:.3f}<extra></extra>'),
            row=row, col=col,
        )
        fig.update_xaxes(title_text=xl, row=row, col=col)
        fig.update_yaxes(title_text=yl, row=row, col=col)

    fig.update_layout(title='%s — %s' % (root, ', '.join(vars)),
                      height=420 * nrows, width=950)
    return fig.to_html(full_html=False, include_plotlyjs=include_plotlyjs)


def make_wind_1d_div(root, master_file, vars=('t_e',), scale='log', inwind='',
                     include_plotlyjs=False):
    """Plotly line plots for a 1D wind model."""
    try:
        data = ascii.read(master_file)
    except IOError:
        print('Warning: cannot read %s' % master_file)
        return None

    rcol = 'rcen' if 'rcen' in data.colnames else ('xcen' if 'xcen' in data.colnames else 'x')
    inwind_arr = np.array(data['inwind'])
    r = np.array(data[rcol], dtype=float)
    mask = inwind_arr >= 0

    fig = go.Figure()
    for var in vars:
        if var not in data.colnames:
            continue
        y = np.array(data[var], dtype=float)
        y[~mask] = np.nan
        if scale == 'log' and var != 'converge':
            with np.errstate(divide='ignore', invalid='ignore'):
                y = np.where(y > 0, np.log10(y), np.nan)
            label = 'log(%s)' % var
        else:
            label = var
        fig.add_trace(go.Scatter(x=r[mask], y=y[mask], mode='lines', name=label))

    fig.update_xaxes(type='log', title='r (cm)')
    fig.update_yaxes(title='log(value)')
    fig.update_layout(title='%s — %s' % (root, ', '.join(vars)),
                      height=450, width=750, legend=dict(x=0.01, y=0.99))
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
    fig.update_layout(title='%s — total spectra' % root,
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


# ---------------------------------------------------------------------------
# HTML builder
# ---------------------------------------------------------------------------

def make_html(root, complete_message, errors,
              conv_by_cycle_div, converge_map_div,
              te_div, tr_div, ne_div,
              spec_tot_div, spec_div):
    """Write the HTML summary, matching run_check.py structure and text."""

    string = xhtml.begin('%s: How well did the run go?' % root)
    string += xhtml.paragraph('Provide an overview of whether the run of %s has succeeded.' % root)
    string += xhtml.add_list(complete_message)
    string += xhtml.hline()

    # --- Did the run converge? ---
    string += xhtml.h2('Did the run converge?')
    if conv_by_cycle_div:
        string += conv_by_cycle_div + '\n'
    else:
        string += xhtml.paragraph('There is no convergence plot. OK if no ionization cycles.')
    string += xhtml.paragraph(convergence_message)
    string += xhtml.paragraph(convergence2)
    if converge_map_div:
        string += converge_map_div + '\n'
    string += xhtml.hline()

    # --- Temperatures ---
    string += xhtml.h2('What do the temperatures look like?')
    if te_div:
        string += te_div + '\n'
    if tr_div:
        string += tr_div + '\n'
    string += xhtml.hline()

    # --- Electron density ---
    string += xhtml.h2('What does the electron density look like?')
    if ne_div:
        string += ne_div + '\n'
    string += xhtml.paragraph(ne_plot_description)
    string += xhtml.hline()

    # --- Total spectra ---
    string += xhtml.h2('What do the total spectra look like (somewhat smoothed)?')
    if spec_tot_div:
        string += spec_tot_div + '\n'
        string += xhtml.paragraph(tot_plot_description)
    else:
        string += xhtml.paragraph('There is no total spectrum plot. OK if no ionization cycles.')
    string += xhtml.hline()

    # --- Angle spectra ---
    string += xhtml.h2('What do the final spectra look like (somewhat smoothed)?')
    if spec_div:
        string += spec_div + '\n'
        string += xhtml.paragraph(spec_plot_description)
    else:
        string += xhtml.paragraph(
            'There is no plot of a detailed spectrum, probably because '
            'detailed spectra were not created.')
    string += xhtml.hline()

    # --- Errors ---
    string += xhtml.h2('Errors and Warnings')
    string += xhtml.paragraph(error_description)
    string += xhtml.preformat(errors)
    string += xhtml.hline()

    # --- Parameter file ---
    string += xhtml.h2('The parameter file used:')
    try:
        with open(root + '.pf') as f:
            lines = f.readlines()
        string += xhtml.add_list(lines)
    except IOError:
        string += xhtml.paragraph('Parameter file %s.pf not found.' % root)
    string += xhtml.hline()

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

    # --- convergence by cycle (embeds Plotly.js once for the whole page) ---
    converged, converging, t_r, t_e, hc = read_diag(root)
    if len(converged) > 1:
        conv_by_cycle_div = make_convergence_div(
            root, converged, converging, t_r, t_e, hc, include_plotlyjs=True)
    else:
        print('Not enough cycles to plot convergence history.')
        conv_by_cycle_div = None

    # All subsequent divs reuse the already-loaded Plotly.js
    kw = dict(include_plotlyjs=False)

    # --- wind structure plots ---
    print('Generating wind plots ...')
    if xdim == 3:
        converge_map_div = make_wind_3d_div(master_file, 'converge',
                                            scale='lin', inwind='all', **kw)
        te_div = make_wind_3d_div(master_file, 't_e', **kw)
        tr_div = make_wind_3d_div(master_file, 't_r', **kw)
        ne_div = make_wind_3d_div(master_file, 'ne',  **kw)
    elif xdim == 2:
        converge_map_div = make_wind_2d_div(root, master_file, vars=('converge',),
                                            scale='lin', inwind='all', **kw)
        te_div = make_wind_2d_div(root, master_file, vars=('t_e',), **kw)
        tr_div = make_wind_2d_div(root, master_file, vars=('t_r',), **kw)
        ne_div = make_wind_2d_div(root, master_file, vars=('ne',),  **kw)
    else:
        converge_map_div = make_wind_1d_div(root, master_file, vars=('converge',),
                                            scale='lin', inwind='all', **kw)
        te_div = make_wind_1d_div(root, master_file, vars=('t_e',), **kw)
        tr_div = make_wind_1d_div(root, master_file, vars=('t_r',), **kw)
        ne_div = make_wind_1d_div(root, master_file, vars=('ne',),  **kw)

    # --- spectral plots ---
    print('Generating spectral plots ...')
    spec_tot_div = make_spec_tot_div(root, **kw)
    spec_div     = make_spec_div(root,     **kw)

    # --- write HTML ---
    errors = py_error(root)
    make_html(root, complete_message, errors,
              conv_by_cycle_div, converge_map_div,
              te_div, tr_div, ne_div,
              spec_tot_div, spec_div)


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
