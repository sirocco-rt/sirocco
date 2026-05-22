#!/usr/bin/env python

"""
plot_wind_3d.py -- Interactive 3D wind visualisation using Plotly.

For CYLIND3D / SPH3D master tables (those with i, j, k grid-index columns)
this script produces a self-contained HTML file containing two panels:

  Left  : rho-z heatmap for each phi slice, browsable via a slider.
  Right : rho-phi heatmap at a chosen z slice (default: smallest |z|).

Both panels share the same colour scale and use physical coordinates.
The rho axis is logarithmic.

Command-line usage
------------------
    plot_wind_3d.py  filename  [variable]  [options]

    variable   column to plot (default: t_e)

    -log       log10 colour scale (default)
    -lin       linear colour scale
    -all       show all cells, not just in-wind ones
    -j N       z-slice index for right panel (default: smallest |z|)
    -o file    output HTML path (default: <root>_<var>_3d.html)

Python usage
------------
    import plot_wind_3d
    plot_wind_3d.doit('agn_cyl3d.master.txt', 't_e')
    plot_wind_3d.doit('agn_cyl3d.master.txt', 'ne', scale='lin', j_slice=25)
"""

import os
import sys
import numpy as np
from astropy.io import ascii
import plotly.graph_objects as go
from plotly.subplots import make_subplots


# ---------------------------------------------------------------------------
# Data preparation
# ---------------------------------------------------------------------------

def get_data(filename, var, scale='log', inwind=''):
    """
    Read a 3D master.txt table and return reshaped arrays.

    Returns a dict with the arrays needed by doit(), or None on error.
    Storage order is i*(mdim*pdim) + j*pdim + k, so a plain reshape works.
    """
    try:
        data = ascii.read(filename)
    except IOError:
        print('Error: cannot read %s' % filename)
        return None

    if 'k' not in data.colnames:
        print('Error: %s has no k column — not a 3D table' % filename)
        return None

    if var not in data.colnames:
        print('Error: variable "%s" not in %s' % (var, filename))
        print('Available columns: ' + ', '.join(data.colnames))
        return None

    ndim = int(np.max(data['i'])) + 1
    mdim = int(np.max(data['j'])) + 1
    pdim = int(np.max(data['k'])) + 1

    inwind_arr = np.array(data['inwind'])
    xvar = np.array(data[var], dtype=float)

    if inwind != 'all':
        xvar[inwind_arr < 0] = np.nan

    # Reshape: storage order is i slowest, k fastest
    xvar_3d = xvar.reshape(ndim, mdim, pdim)

    # Coordinate centres — extract from k=0 slice and i=0,j=0 row
    k0 = np.array(data['k']) == 0
    rho_cen = np.array(data['xcen'][k0]).reshape(ndim, mdim)[:, 0]
    z_cen   = np.array(data['zcen'][k0]).reshape(ndim, mdim)[0, :]

    ij0 = (np.array(data['i']) == 0) & (np.array(data['j']) == 0)
    phi_cen = np.array(data['phicen'][ij0])

    title = var
    if scale == 'log':
        with np.errstate(divide='ignore', invalid='ignore'):
            xvar_3d = np.where(xvar_3d > 0, np.log10(xvar_3d), np.nan)
        title = 'log₁₀(' + var + ')'

    return dict(xvar_3d=xvar_3d, rho_cen=rho_cen, z_cen=z_cen,
                phi_cen=phi_cen, ndim=ndim, mdim=mdim, pdim=pdim,
                title=title)


# ---------------------------------------------------------------------------
# Main plotting function
# ---------------------------------------------------------------------------

def doit(filename, var='t_e', scale='log', outfile=None, inwind='', j_slice=None):
    """
    Create an interactive HTML visualisation for a 3D wind variable.

    Parameters
    ----------
    filename : str   Path to a master.txt file with i, j, k columns.
    var      : str   Variable to plot (default 't_e').
    scale    : str   'log' (default) or 'lin'.
    outfile  : str   Output HTML path.  Default: <root>_<var>_3d.html.
    inwind   : str   '' = in-wind cells only (default); 'all' = show all.
    j_slice  : int   j index for the rho-phi panel (default: |z| minimum).
    """
    d = get_data(filename, var, scale=scale, inwind=inwind)
    if d is None:
        return None

    xvar_3d = d['xvar_3d']
    rho_cen = d['rho_cen']
    z_cen   = d['z_cen']
    phi_cen = d['phi_cen']
    ndim, mdim, pdim = d['ndim'], d['mdim'], d['pdim']
    title   = d['title']

    # Shared colour limits
    finite = xvar_3d[np.isfinite(xvar_3d)]
    if len(finite) == 0:
        print('Warning: no finite in-wind values for %s; try -all' % var)
        zmin, zmax = 0.0, 1.0
    else:
        zmin = float(np.nanmin(finite))
        zmax = float(np.nanmax(finite))

    # j index for the rho-phi panel
    if j_slice is None:
        j_slice = int(np.argmin(np.abs(z_cen)))

    # ------------------------------------------------------------------
    # Build figure with two side-by-side subplots
    # ------------------------------------------------------------------
    phi_deg = np.degrees(phi_cen)
    z_label = 'z = %.2e cm  (j=%d)' % (z_cen[j_slice], j_slice)

    fig = make_subplots(
        rows=1, cols=2,
        subplot_titles=[
            'rho-z plane  (use slider to change ϕ)',
            'rho-ϕ plane  (' + z_label + ')',
        ],
        horizontal_spacing=0.18,
    )

    # --- Left panel: one Heatmap trace per phi slice --------------------
    for k in range(pdim):
        # Heatmap expects z[row, col] with row = y-axis, col = x-axis
        z_slice = xvar_3d[:, :, k].T    # (mdim, ndim): rows=z, cols=rho
        fig.add_trace(
            go.Heatmap(
                z=z_slice,
                x=rho_cen,
                y=z_cen,
                colorscale='Viridis',
                zmin=zmin, zmax=zmax,
                showscale=True,
                colorbar=dict(x=0.40, thickness=15, title=dict(text=title, side='right')),
                visible=(k == 0),
                name='ϕ=%.0f°' % phi_deg[k],
                hovertemplate='rho=%{x:.2e}<br>z=%{y:.2e}<br>' + title + '=%{z:.3f}<extra></extra>',
            ),
            row=1, col=1,
        )

    # --- Right panel: rho-phi slice at chosen j (always visible) -------
    eq_slice = xvar_3d[:, j_slice, :].T    # (pdim, ndim): rows=phi, cols=rho
    fig.add_trace(
        go.Heatmap(
            z=eq_slice,
            x=rho_cen,
            y=phi_deg,
            colorscale='Viridis',
            zmin=zmin, zmax=zmax,
            showscale=True,
            colorbar=dict(x=1.01, thickness=15, title=dict(text=title, side='right')),
            name='rho-phi',
            hovertemplate='rho=%{x:.2e}<br>ϕ=%{y:.0f}°<br>' + title + '=%{z:.3f}<extra></extra>',
        ),
        row=1, col=2,
    )

    # --- Phi slider (controls left panel only) --------------------------
    n_left = pdim   # indices 0..pdim-1 are the left traces; pdim is the right trace
    steps = []
    for k in range(pdim):
        visible = [kk == k for kk in range(n_left)] + [True]
        steps.append(dict(
            method='update',
            args=[{'visible': visible}],
            label='%.0f°' % phi_deg[k],
        ))

    sliders = [dict(
        active=0,
        steps=steps,
        currentvalue=dict(prefix='ϕ = ', font=dict(size=14)),
        pad=dict(t=55),
        x=0.0, len=0.42,
    )]

    # --- Axes and layout ------------------------------------------------
    fig.update_xaxes(title_text='rho (cm)', type='log', row=1, col=1)
    fig.update_yaxes(title_text='z (cm)',               row=1, col=1)
    fig.update_xaxes(title_text='rho (cm)', type='log', row=1, col=2)
    fig.update_yaxes(title_text='ϕ (degrees)',     row=1, col=2)

    fig.update_layout(
        title_text='%s  —  %s' % (os.path.basename(filename), title),
        sliders=sliders,
        height=560,
        width=1150,
        margin=dict(t=100, b=80),
    )

    # --- Save -----------------------------------------------------------
    if outfile is None:
        root = os.path.basename(filename)
        for ext in ('.master.txt', '.txt'):
            root = root.replace(ext, '')
        outfile = '%s_%s_3d.html' % (root, var)

    fig.write_html(outfile)
    print('Written: %s' % outfile)
    return outfile


# ---------------------------------------------------------------------------
# Command-line interface
# ---------------------------------------------------------------------------

def steer(argv):
    filename = ''
    var = 't_e'
    scale = 'log'
    outfile = None
    inwind = ''
    j_slice = None

    i = 1
    while i < len(argv):
        if argv[i] == '-h':
            print(__doc__)
            return
        elif argv[i] == '-log':
            scale = 'log'
        elif argv[i] == '-lin':
            scale = 'lin'
        elif argv[i] == '-all':
            inwind = 'all'
        elif argv[i] == '-j':
            i += 1
            j_slice = int(argv[i])
        elif argv[i] == '-o':
            i += 1
            outfile = argv[i]
        elif filename == '':
            filename = argv[i]
        elif var == 't_e':
            var = argv[i]
        else:
            print('Unrecognised argument: %s' % argv[i])
            return
        i += 1

    if filename == '':
        print('Usage: plot_wind_3d.py filename [variable] [-log|-lin] [-all] [-j N] [-o file.html]')
        return

    doit(filename, var=var, scale=scale, outfile=outfile, inwind=inwind, j_slice=j_slice)


if __name__ == '__main__':
    steer(sys.argv)
