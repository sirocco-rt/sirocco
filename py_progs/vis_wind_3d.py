#!/usr/bin/env python

"""
vis_wind_3d.py -- Rotatable 3D scatter plots from a CYLIND3D or SPH3D master table.

The table must have xcen, ycen, zcen columns (written by windsave2table for
CYLIND3D and SPH3D models).  Only cells with inwind >= 0 are shown.  Each
variable is coloured by log10(value).

Command-line usage
------------------
    vis_wind_3d.py  master.txt  var1  [var2  var3 ...]  [options]

    master.txt   path to a windsave2table master.txt file
    var1 ...     one or more column names, e.g. rho ne t_e

Options
-------
    -q           quantile (rank-based) colour normalisation: equal numbers of
                 cells fall in each colour band.  The colorbar tick labels
                 show the actual log10 values at selected quantiles.
    -p lo hi     clip the colour scale to the lo-th and hi-th percentiles of
                 the log10 values (e.g. -p 2 98).  Ignored when -q is used.
    -a alpha     marker opacity, 0–1 (default 1.0).
    -s size      marker size in pixels (default 3).
    -c scale     Plotly colorscale name, e.g. Plasma, Inferno, Hot (default Viridis).
                 Run `python -c "import plotly.express as px; print(px.colors.named_colorscales())"``
                 for the full list of ~70 named scales.
    -o file      output HTML path (default: <root>_<vars>.html beside input).

Output
------
    A self-contained HTML file.  Open in any browser; rotating one panel
    rotates all panels to the same viewpoint.

Python usage
------------
    from vis_wind_3d import vis_wind_3d
    vis_wind_3d('run.master.txt', ['rho', 'ne', 't_e'],
                quantile=True, alpha=0.7)
"""

import sys
import os
import argparse
import numpy as np
from astropy.io import ascii
from scipy.stats import rankdata
from plotly.subplots import make_subplots
import plotly.graph_objects as go


# ---------------------------------------------------------------------------

def _strip_root(path):
    """Return the file root, stripping .master.txt or other extensions."""
    base = os.path.basename(path)
    if base.endswith('.master.txt'):
        return base[:-11]
    root, _ = os.path.splitext(base)
    if root.endswith('.master'):
        root = root[:-7]
    return root


def _quantile_color(log_vals):
    """
    Map log_vals to rank-normalised colour values in [0, 1].

    Returns
    -------
    color_vals : ndarray
        Normalised ranks in [0, 1]; NaN where log_vals is NaN.
    tick_ranks : list of float
        Rank positions for colorbar ticks.
    tick_text : list of str
        Actual log10 values at those ranks, for colorbar labels.
    """
    valid = np.isfinite(log_vals)
    color_vals = np.full_like(log_vals, np.nan)
    n = int(np.sum(valid))
    if n == 0:
        return color_vals, [], []

    # Map to (0, 1) using midpoint ranks so no point sits exactly at 0 or 1
    color_vals[valid] = (rankdata(log_vals[valid]) - 0.5) / n

    # Tick positions at 9 evenly spaced quantiles
    pcts = [0, 12.5, 25, 37.5, 50, 62.5, 75, 87.5, 100]
    tick_ranks = [p / 100.0 for p in pcts]
    tick_logvals = np.nanpercentile(log_vals, pcts)
    tick_text = [f'{v:.2g}' for v in tick_logvals]

    return color_vals, tick_ranks, tick_text


# ---------------------------------------------------------------------------

def vis_wind_3d(tablefile, variables, outfile=None,
                quantile=False, plo=None, phi=None,
                alpha=1.0, marker_size=3, colorscale='Viridis'):
    """
    Create rotatable 3D scatter HTML panels from a master.txt table.

    Parameters
    ----------
    tablefile : str
        Path to the windsave2table master.txt file.
    variables : list of str
        Column names to plot (one panel per variable).
    outfile : str or None
        Output HTML path.  Defaults to <root>_<vars>.html beside tablefile.
    quantile : bool
        If True, use rank-based colour normalisation.
    plo, phi : float or None
        Percentile clip for colour scale (e.g. 2, 98).  Ignored if quantile=True.
    alpha : float
        Marker opacity in [0, 1].
    marker_size : int
        Marker size in pixels.
    colorscale : str
        Plotly named colorscale (default 'Viridis').
    """

    t = ascii.read(tablefile)
    t = t[t['inwind'] >= 0]

    for col in ('xcen', 'ycen', 'zcen'):
        if col not in t.colnames:
            raise ValueError(
                f"Column '{col}' not found in {tablefile}. "
                "Is this a CYLIND3D or SPH3D master.txt?"
            )

    X = np.array(t['xcen'], dtype=float)
    Y = np.array(t['ycen'], dtype=float)
    Z = np.array(t['zcen'], dtype=float)

    valid_vars = []
    for var in variables:
        if var in t.colnames:
            valid_vars.append(var)
        else:
            print(f"Warning: '{var}' not found in table, skipping")
    if not valid_vars:
        raise ValueError("No valid variables to plot.")

    ncols = len(valid_vars)
    gap = 0.05
    subplot_w = (1.0 - (ncols - 1) * gap) / ncols

    fig = make_subplots(
        rows=1,
        cols=ncols,
        specs=[[{'type': 'scatter3d'}] * ncols],
        subplot_titles=valid_vars,
        horizontal_spacing=gap,
    )

    for idx, var in enumerate(valid_vars):
        col = idx + 1

        vals = np.array(t[var], dtype=float)
        with np.errstate(divide='ignore', invalid='ignore'):
            log_vals = np.where(vals > 0, np.log10(vals), np.nan)

        cb_x = idx * (subplot_w + gap) + subplot_w + 0.01

        if quantile:
            color_vals, tick_ranks, tick_text = _quantile_color(log_vals)
            colorbar = dict(
                title=f'log {var}',
                x=cb_x,
                len=0.75,
                thickness=15,
                tickvals=tick_ranks,
                ticktext=tick_text,
            )
            cmin, cmax = 0.0, 1.0
        else:
            color_vals = log_vals
            lo = plo if plo is not None else 0
            hi = phi if phi is not None else 100
            cmin = np.nanpercentile(log_vals, lo)
            cmax = np.nanpercentile(log_vals, hi)
            colorbar = dict(
                title=f'log {var}',
                x=cb_x,
                len=0.75,
                thickness=15,
            )

        fig.add_trace(
            go.Scatter3d(
                x=X, y=Y, z=Z,
                mode='markers',
                marker=dict(
                    size=marker_size,
                    color=color_vals,
                    colorscale=colorscale,
                    cmin=cmin,
                    cmax=cmax,
                    opacity=alpha,
                    colorbar=colorbar,
                    showscale=True,
                ),
                name=var,
            ),
            row=1,
            col=col,
        )

        scene_key = 'scene' if col == 1 else f'scene{col}'
        fig.update_layout(**{
            scene_key: dict(
                xaxis_title='x (cm)',
                yaxis_title='y (cm)',
                zaxis_title='z (cm)',
                aspectmode='data',
            )
        })

    fig.update_layout(
        title_text=os.path.basename(tablefile),
        height=700,
        width=max(700, 550 * ncols),
    )

    if outfile is None:
        dirname = os.path.dirname(tablefile) or '.'
        root = _strip_root(tablefile)
        varstr = '_'.join(valid_vars)
        outfile = os.path.join(dirname, f"{root}_{varstr}.html")

    plot_div = fig.to_html(full_html=False, include_plotlyjs='cdn', div_id='vw3d')

    sync_js = (
        "<script>\n"
        "(function () {\n"
        "    var plot = document.getElementById('vw3d');\n"
        "    var _busy = false;\n"
        "    plot.on('plotly_relayout', function (ev) {\n"
        "        if (_busy) return;\n"
        "        var camera = null;\n"
        "        for (var k of Object.keys(ev)) {\n"
        r"            if (/^scene\d*\.camera$/.test(k)) { camera = ev[k]; break; }" "\n"
        "        }\n"
        "        if (!camera) return;\n"
        "        var update = {};\n"
        "        for (var k of Object.keys(plot.layout)) {\n"
        r"            if (/^scene\d*$/.test(k)) update[k + '.camera'] = camera;" "\n"
        "        }\n"
        "        _busy = true;\n"
        "        Plotly.relayout(plot, update).then(function () { _busy = false; });\n"
        "    });\n"
        "})();\n"
        "</script>\n"
    )

    title = os.path.basename(tablefile)
    html = (
        "<!DOCTYPE html>\n<html>\n"
        f"<head><meta charset='utf-8'><title>{title}</title></head>\n"
        f"<body>\n{plot_div}\n{sync_js}\n</body>\n</html>\n"
    )

    with open(outfile, 'w') as f:
        f.write(html)

    print(f"Written {outfile}")
    return outfile


# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='Rotatable 3D wind scatter plots from a master.txt table.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument('tablefile', help='windsave2table master.txt file')
    parser.add_argument('variables', nargs='+', help='column names to plot')
    parser.add_argument('-q', '--quantile', action='store_true',
                        help='quantile (rank-based) colour normalisation')
    parser.add_argument('-p', nargs=2, type=float, metavar=('LO', 'HI'),
                        help='percentile clip for colour scale (e.g. -p 2 98)')
    parser.add_argument('-a', '--alpha', type=float, default=1.0,
                        help='marker opacity 0-1 (default 1.0)')
    parser.add_argument('-s', '--size', type=int, default=3,
                        help='marker size in pixels (default 3)')
    parser.add_argument('-c', '--colorscale', default='Viridis',
                        help='Plotly colorscale name (default Viridis)')
    parser.add_argument('-o', '--outfile', default=None,
                        help='output HTML path')

    args = parser.parse_args()

    plo = args.p[0] if args.p else None
    phi = args.p[1] if args.p else None

    vis_wind_3d(
        args.tablefile,
        args.variables,
        outfile=args.outfile,
        quantile=args.quantile,
        plo=plo,
        phi=phi,
        alpha=args.alpha,
        marker_size=args.size,
        colorscale=args.colorscale,
    )


if __name__ == '__main__':
    main()
