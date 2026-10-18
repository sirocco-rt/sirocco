#!/usr/bin/env python
'''
A cartoon of the geometries of three of the polarization tests, drawn in the
x-z plane (all are axisymmetric and mirrored about z = 0), with the density
of each model from the formula its generator uses:

    brown_mclean   n_e = n0 (r_in/r)^2 sin^2(theta), r_in < r < 10 r_in,
                   around a star of 0.1 r_in (make_bm.py)
    hillier        chi = chi0 (R_min/r)^4 (1 + 10 cos^2(theta)),
                   R_min < r < 30 R_min, around a star of 0.2 R_min
                   (make_hillier.py)
    sr_shell       the Hillier shell, optically thin, moving radially at
                   beta c, around a star of 0.01 R_min (make_srshell.py)

The observers' directions and the E-vector of the scattered light are
marked: a bar along the projected axis for q > 0, a circle (E perpendicular
to the figure, parallel to the equator) for q < 0.  The q values are from the
runs on which the tests' tolerances are based.

usage: shell_cartoons.py            writes shell_cartoons.png
'''
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle


def density_map(ax, extent, rmin, rmax, ne, ndec=3, cmap='Blues'):
    '''Shade log10 of ne(r, theta) between rmin and rmax, over ndec decades below its maximum'''
    x = np.linspace(-extent, extent, 600)
    X, Z = np.meshgrid(x, x)
    R = np.hypot(X, Z)
    mu = np.divide(Z, R, out=np.zeros_like(R), where=R > 0)
    with np.errstate(divide='ignore', invalid='ignore'):
        # shade out to rmax or the edge of the plotted region, whichever is smaller, so the map ends in a circle
        d = np.where((R >= rmin) & (R <= min(rmax, extent)), ne(R, mu), np.nan)
        logd = np.log10(d)
    top = np.nanmax(logd)
    ax.imshow(logd, extent=(-extent, extent, -extent, extent), origin='lower', cmap=cmap,
              vmin=top - ndec, vmax=top, interpolation='bilinear')


def observers(ax, incs, q, radius, label_radius):
    '''Lines of sight in the upper right quadrant, with the E-vector of the scattered light'''
    for inc in incs:
        t = np.radians(inc)
        d = np.array([np.sin(t), np.cos(t)])
        p = radius * d
        ax.plot([0, 0.92 * radius * d[0]], [0, 0.92 * radius * d[1]], ls=':', color='0.55', lw=0.7)
        ax.annotate('', xy=0.97 * p, xytext=0.90 * p, arrowprops=dict(arrowstyle='->', color='k', lw=1))
        north = np.array([-np.cos(t), np.sin(t)])          # the projected +z axis on the sky, in this plane
        s = 0.04 * radius
        if q[inc] > 0:
            ax.plot(*np.array([p + 0.05 * radius * d - s * north, p + 0.05 * radius * d + s * north]).T,
                    color='#cb181d', lw=2.2)
        else:
            ax.plot(*(p + 0.05 * radius * d), marker='o', ms=7, mfc='none', mec='#cb181d', mew=2)
        ax.text(*(label_radius * d), '%g$^\\circ$\nq = %+.2g%%' % (inc, q[inc]), fontsize=8, ha='center',
                va='center')


fig, axes = plt.subplots(1, 3, figsize=(17, 6.4))

# Brown & McLean: an equatorially concentrated envelope, r_in = 1, r_out = 10 (units of r_in = 1e11 cm)
ax = axes[0]
density_map(ax, 10.5, 1.0, 10.0, lambda r, mu: r ** -2 * (1 - mu ** 2))
ax.add_patch(Circle((0, 0), 0.1, color='#fd8d3c', zorder=5))
observers(ax, (10, 30, 50, 70, 90), {10: 0.022, 30: 0.22, 50: 0.54, 70: 0.83, 90: 0.95}, 11.0, 14.0)
ax.set_title('brown_mclean: n$_e \\propto$ r$^{-2}$ sin$^2\\theta$, r$_{in}$ < r < 10 r$_{in}$', fontsize=10)
ax.text(0.02, 0.02, 'densest at the equator, empty on the axis\nstar 0.1 r$_{in}$; pure H, electron scattering only\n'
        'q > 0: E along the axis (q for $\\tau_{eq}$ = 0.1)', transform=ax.transAxes, fontsize=8.5, va='bottom')
ax.set_xlabel('x (r$_{in}$ = 10$^{11}$ cm)')
ax.set_ylabel('z (r$_{in}$)')

# Hillier: a detached prolate shell (units of R_min), drawn to 6 R_min
ax = axes[1]
shell = lambda r, mu: r ** -4 * (1 + 10 * mu ** 2)
density_map(ax, 6.3, 1.0, 30.0, shell)
ax.add_patch(Circle((0, 0), 0.2, color='#fd8d3c', zorder=5))
observers(ax, (11.25, 33.75, 56.25, 78.75), {11.25: -0.21, 33.75: -1.5, 56.25: -2.6, 78.75: -3.1}, 6.6, 8.5)
ax.set_title('hillier: $\\chi \\propto$ r$^{-4}$ (1 + 10 cos$^2\\theta$), R$_{min}$ < r < 30 R$_{min}$', fontsize=10)
ax.text(0.02, 0.02, '11$\\times$ denser on the axis than at the equator\nstar 0.2 R$_{min}$; drawn to 6 R$_{min}$\n'
        'q < 0: E parallel to the equator (q for $\\chi_0$ = 0.1;\nq turns positive at large i for $\\chi_0 \\geq$ 1)',
        transform=ax.transAxes, fontsize=8.5, va='bottom')
ax.set_xlabel('x (R$_{min}$)')
ax.set_ylabel('z (R$_{min}$)')

# Relativistic shell: the Hillier shape, optically thin, moving radially at beta c
ax = axes[2]
density_map(ax, 6.3, 1.0, 30.0, shell, cmap='Purples')
ax.add_patch(Circle((0, 0), 0.12, color='#fd8d3c', zorder=5))
for th in np.radians(np.arange(0, 360, 30)):
    d = np.array([np.sin(th), np.cos(th)])
    ax.annotate('', xy=3.6 * d, xytext=1.6 * d, arrowprops=dict(arrowstyle='->', color='#54278f', lw=1.1))
observers(ax, (22.5, 45, 67.5, 90), {22.5: -2.6, 45: -9.9, 67.5: -19.5, 90: -24.2}, 6.6, 8.5)
ax.set_title('sr_shell: the Hillier shell, thin ($\\chi_0$ = 0.003), v = $\\beta$c radial', fontsize=10)
ax.text(0.02, 0.02, '$\\beta$ = 3$\\times$10$^{-4}$, 0.1, 0.3, and 0.3 with -nonrel; also a spherical shell\n'
        'star 0.01 R$_{min}$ (shown enlarged); drawn to 6 R$_{min}$\n'
        'q < 0 (q of the light scattered once, $\\beta$ = 0.3)', transform=ax.transAxes, fontsize=8.5, va='bottom')
ax.set_xlabel('x (R$_{min}$)')
ax.set_ylabel('z (R$_{min}$)')

for ax, l in zip(axes, (16, 10, 10)):
    ax.set_xlim(-l * 0.70, l)
    ax.set_ylim(-l * 0.70, l)
    ax.set_aspect('equal')
    ax.axhline(0, color='0.85', lw=0.5, zorder=0)
    ax.axvline(0, color='0.85', lw=0.5, zorder=0)
fig.suptitle('Geometries of the brown_mclean, hillier and sr_shell polarization tests (shading: log density over 3 '
             'decades; orange: the star)\nred bar: E-vector in the plane of the figure, along the projected axis '
             '(q > 0); red circle: E-vector perpendicular to the figure, parallel to the equator (q < 0)',
             fontsize=9.5)
fig.tight_layout(rect=(0, 0, 1, 0.92))
fig.savefig('shell_cartoons.png', dpi=130, bbox_inches='tight')
print('shell_cartoons.png')
