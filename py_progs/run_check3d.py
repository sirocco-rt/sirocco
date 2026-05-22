#!/usr/bin/env python

"""
run_check3d.py -- Comprehensive run diagnostics for all model dimensionalities.

Intended as a future replacement for run_check.py.  Handles 1D, 2D, and 3D
wind models and adds interactive Plotly plots alongside the traditional static
matplotlib summaries.

Interactive HTML files are stored in diag_<root>/interactive/ and linked
relatively from the main HTML summary page.

If Plotly is not available the interactive section is skipped gracefully and
the static plots are still produced.

Usage
-----
    run_check3d.py root [root ...]
    run_check3d.py -all          # process every .wind_save in current directory
"""

import os
import sys
import subprocess
import numpy
import numpy as np
from glob import glob
from astropy.io import ascii
from astropy.table import Table
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as pylab

import xhtml
import plot_wind
import plot_wind_1d
import plot_spec
import plot_tot

# ---------------------------------------------------------------------------
# Optional Plotly imports — graceful fallback
# ---------------------------------------------------------------------------
try:
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
    import plot_wind_3d
    PLOTLY_AVAILABLE = True
except ImportError:
    PLOTLY_AVAILABLE = False

if not PLOTLY_AVAILABLE:
    print('Note: Plotly not found — interactive plots will not be generated.')
    print('      Install with:  pip install plotly')


# ---------------------------------------------------------------------------
# Utility functions  (copied from run_check.py so this file is self-contained)
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
# Static matplotlib plots  (same as run_check.py)
# ---------------------------------------------------------------------------

def plot_converged(root, converged, converging, t_r, t_e, hc):
    """Save a convergence-by-cycle plot to diag_<root>/convergence.png."""
    pylab.close(1)
    pylab.figure(1, (6, 6))
    pylab.plot(t_r, '--', label='t_r')
    pylab.plot(t_e, '--', label='t_e')
    pylab.plot(hc,  '--', label='heating:cooling')
    pylab.plot(converged,  lw=2, label='Converged')
    pylab.plot(converging, lw=2, label='Converging')
    pylab.legend(loc='best')
    pylab.xlabel('Ionization Cycle')
    pylab.ylabel('Fraction of Cells in Wind')
    pylab.savefig('diag_%s/convergence.png' % root)


# ---------------------------------------------------------------------------
# Interactive Plotly plots
# ---------------------------------------------------------------------------

def _smooth(x, n=21):
    """Simple boxcar smooth (same algorithm as plot_tot.xsmooth)."""
    if n <= 1 or len(x) <= n:
        return x
    kernel = np.ones(n) / n
    return np.convolve(x, kernel, mode='same')


def make_interactive_spec_tot(root, outdir):
    """
    Create a Plotly nuLnu vs frequency plot from <root>.log_spec_tot.
    Returns the output path on success, or None.
    """
    filename = root + '.log_spec_tot'
    try:
        data = ascii.read(filename)
    except IOError:
        print('Warning: %s not found; skipping interactive spec_tot' % filename)
        return None

    freq = np.array(data['Freq.'], dtype=float)

    fig = go.Figure()

    for col, label, dash in [
        ('Created',  'Created (ext)',    'solid'),
        ('Emitted',  'Observed (tot)',   'solid'),
        ('Wind',     'Observed (wind)',  'dash'),
        ('HitSurf',  'Hit surface',      'dot'),
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

    # Mark key spectral lines
    for freq_val, label in [(3.29e15, 'Ly limit'), (2.47e15, 'Lyα'),
                             (1.31e16, 'HeII'), (6.91e14, 'Hβ')]:
        fig.add_vline(x=freq_val, line_width=1, line_dash='dot', line_color='grey',
                      annotation_text=label, annotation_position='top right')

    fig.update_xaxes(type='log', title='Frequency (Hz)')
    fig.update_yaxes(type='log', title='νL<sub>ν</sub>')
    fig.update_layout(title='%s — spec_tot' % root, height=500, width=950,
                      legend=dict(x=0.01, y=0.99))

    outfile = os.path.join(outdir, '%s_spec_tot_interactive.html' % root)
    fig.write_html(outfile, include_plotlyjs=True)
    return outfile


def make_interactive_spec(root, outdir):
    """
    Create a Plotly flux vs wavelength plot from <root>.spec.
    All inclination angles are shown as separate traces, togglable via the legend.
    Returns the output path on success, or None.
    """
    filename = root + '.spec'
    try:
        data = ascii.read(filename)
    except IOError:
        print('Warning: %s not found; skipping interactive spec' % filename)
        return None

    wave = np.array(data['Lambda'], dtype=float)

    # Find angle columns (start with 'A')
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
        height=500, width=950,
        legend=dict(x=0.01, y=0.99),
    )

    outfile = os.path.join(outdir, '%s_spec_interactive.html' % root)
    fig.write_html(outfile, include_plotlyjs=True)
    return outfile


def make_interactive_wind(root, master_file, outdir):
    """
    Create a Plotly 3D wind plot via plot_wind_3d for a 3D model.
    Returns the output path on success, or None.
    """
    outfile = os.path.join(outdir, '%s_wind_interactive.html' % root)
    result = plot_wind_3d.doit(master_file, var=['t_e', 't_r', 'ne', 'rho'],
                                scale='log', outfile=outfile)
    return result


# ---------------------------------------------------------------------------
# HTML summary
# ---------------------------------------------------------------------------

_convergence_description = '''
The plot shows the fraction of cells satisfying convergence criteria for the
radiation temperature (t_r), electron temperature (t_e), and
heating-to-cooling ratio (hc).  A cell "converges" when all three criteria
are met; "converging" cells are evolving monotonically rather than
oscillating.
'''

_ne_description = '''
Electron density throughout the wind (k=0 phi slice for 3D models).
Traces the density structure; directly related to wind optical depth.
'''

_spec_tot_description = '''
Photon budget from the last ionization cycle: photons created externally
(star/disk), observed total, observed from wind, and those hitting the disk
or star surface (nuLnu vs frequency).
'''

_spec_description = '''
Angle-resolved spectra from the .spec file (flux vs wavelength).
'''

_interactive_description = '''
Interactive plots allow zooming, hovering for exact values, and (for 3D
wind models) sliders to explore different phi and z/theta slices.
'''


def make_html(root, converge_plot, te_plot, tr_plot, ne_plot,
              spec_tot_plot, spec_plot, nspectra,
              complete_message, errors,
              interactive_spec_tot=None,
              interactive_spec=None,
              interactive_wind=None):
    """Write the HTML summary page."""

    string = xhtml.begin('%s: Run summary' % root)
    string += xhtml.paragraph('Summary for %s' % root)
    string += xhtml.hline()

    # Completion status
    string += xhtml.h2('Run status')
    for msg in complete_message:
        string += xhtml.paragraph(msg)

    # Convergence plot
    if os.path.isfile('diag_%s/convergence.png' % root):
        string += xhtml.h2('Convergence by cycle')
        string += xhtml.paragraph(_convergence_description)
        string += xhtml.image('file:./diag_%s/convergence.png' % root)

    # Wind structure plots
    string += xhtml.h2('Wind structure (static)')
    for plot, label in [(converge_plot, 'Convergence map'),
                        (te_plot,       'Electron temperature'),
                        (tr_plot,       'Radiation temperature'),
                        (ne_plot,       'Electron density')]:
        if plot and plot != 'none' and os.path.isfile(plot):
            string += xhtml.paragraph(label)
            string += xhtml.image('file:%s' % plot)
    string += xhtml.paragraph(_ne_description)

    # Spectral plots (static)
    string += xhtml.h2('Spectra (static)')
    if spec_tot_plot and os.path.isfile(spec_tot_plot):
        string += xhtml.paragraph(_spec_tot_description)
        string += xhtml.image('file:%s' % spec_tot_plot, width=800)
    if spec_plot and spec_plot != 'None' and spec_plot != 'none':
        if os.path.isfile(str(spec_plot)):
            string += xhtml.paragraph(_spec_description)
            string += xhtml.image('file:%s' % spec_plot, width=800)

    # Interactive plots section
    have_interactive = any([interactive_spec_tot, interactive_spec, interactive_wind])
    if have_interactive:
        string += xhtml.h2('Interactive plots')
        string += xhtml.paragraph(_interactive_description)
        links = []
        if interactive_spec_tot:
            links.append(xhtml.link('Spectral energy distribution (interactive)',
                                    href=interactive_spec_tot))
        if interactive_spec:
            links.append(xhtml.link('Angle spectra (interactive)',
                                    href=interactive_spec))
        if interactive_wind:
            links.append(xhtml.link('3D wind structure (interactive)',
                                    href=interactive_wind))
        string += xhtml.add_list(links)
    elif not PLOTLY_AVAILABLE:
        string += xhtml.h2('Interactive plots')
        string += xhtml.paragraph(
            'Plotly is not installed; interactive plots were not generated. '
            'Install with: pip install plotly')

    # Errors
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

    # Locate master file
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

    # Ensure the interactive output directory exists
    interactive_dir = './diag_%s/interactive' % root
    if PLOTLY_AVAILABLE:
        os.makedirs(interactive_dir, exist_ok=True)

    # ----------------------------------------------------------------
    # Static wind plots (matplotlib)
    # ----------------------------------------------------------------
    plot_dir = './diag_%s' % root
    if xdim == 3:
        print('3D grid: generating k=0 phi-slice static wind plots')
        converge_plot = plot_wind.doit(master_file, 'converge', plot_dir=plot_dir, k=0)
        te_plot       = plot_wind.doit(master_file, 't_e',      plot_dir=plot_dir, k=0)
        tr_plot       = plot_wind.doit(master_file, 't_r',      plot_dir=plot_dir, k=0)
        ne_plot       = plot_wind.doit(master_file, 'ne',       plot_dir=plot_dir, k=0)
    elif xdim == 2:
        converge_plot = plot_wind.doit(master_file, 'converge', plot_dir=plot_dir)
        te_plot       = plot_wind.doit(master_file, 't_e',      plot_dir=plot_dir)
        tr_plot       = plot_wind.doit(master_file, 't_r',      plot_dir=plot_dir)
        ne_plot       = plot_wind.doit(master_file, 'ne',       plot_dir=plot_dir)
    else:
        converge_plot = plot_wind_1d.doit(master_file, 'converge', plot_dir=plot_dir)
        te_plot       = plot_wind_1d.doit(master_file, 't_e',      plot_dir=plot_dir)
        tr_plot       = plot_wind_1d.doit(master_file, 't_r',      plot_dir=plot_dir)
        ne_plot       = plot_wind_1d.doit(master_file, 'ne',       plot_dir=plot_dir)

    # ----------------------------------------------------------------
    # Convergence cycle plot
    # ----------------------------------------------------------------
    converged, converging, t_r, t_e, hc = read_diag(root)
    if len(converged) > 1:
        plot_converged(root, converged, converging, t_r, t_e, hc)
    else:
        print('Not enough cycles to plot convergence history')

    # ----------------------------------------------------------------
    # Static spectral plots
    # ----------------------------------------------------------------
    plot_tot.doit(root)
    spec_tot_plot = root + '.spec_tot.png'
    try:
        spec_plot, nspectra = plot_spec.do_mosaic(root, wmin=0, wmax=0)
    except Exception:
        print('Warning: could not generate static spectrum plot')
        spec_plot = 'none'
        nspectra = 0

    # ----------------------------------------------------------------
    # Interactive Plotly plots
    # ----------------------------------------------------------------
    interactive_spec_tot = None
    interactive_spec     = None
    interactive_wind     = None

    if PLOTLY_AVAILABLE:
        print('Generating interactive Plotly plots ...')
        interactive_spec_tot = make_interactive_spec_tot(root, interactive_dir)
        interactive_spec     = make_interactive_spec(root, interactive_dir)
        if xdim == 3:
            interactive_wind = make_interactive_wind(root, master_file, interactive_dir)

        # Convert absolute paths to relative paths for HTML links
        def _rel(path):
            if path and os.path.isabs(path):
                return os.path.relpath(path)
            return path
        interactive_spec_tot = _rel(interactive_spec_tot)
        interactive_spec     = _rel(interactive_spec)
        interactive_wind     = _rel(interactive_wind)
    else:
        print('Skipping interactive plots (Plotly not available)')

    # ----------------------------------------------------------------
    # Error summary
    # ----------------------------------------------------------------
    errors = py_error(root)

    # ----------------------------------------------------------------
    # Write HTML summary
    # ----------------------------------------------------------------
    make_html(root, converge_plot, te_plot, tr_plot, ne_plot,
              spec_tot_plot, spec_plot, nspectra,
              complete_message, errors,
              interactive_spec_tot=interactive_spec_tot,
              interactive_spec=interactive_spec,
              interactive_wind=interactive_wind)


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
