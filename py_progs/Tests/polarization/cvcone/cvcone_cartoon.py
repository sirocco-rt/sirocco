#!/usr/bin/env python
'''
A cartoon of the two CV-like wind models of the cvcone polarization test:
a polar cone (pbl_es, 0-20 deg) and an equatorial wedge (eq_es, 60-88 deg),
drawn to scale in the x-z plane (the models are axisymmetric and mirrored
about the disk plane).

The q values are those of the light scattered once in the 261017 runs.

usage: cvcone_cartoon.py            writes cvcone_cartoon.png
'''
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Wedge, Circle

R_IN, R_OUT, R_DISK = 0.28, 10.0, 2.4          # 1e10 cm (make_cvcone.py defaults)
MODELS = [('pbl_es: polar cone, 0-20$^\\circ$', 0, 20, 'q < 0: E parallel to the disk',
           {10: -1.9, 45: -30.9, 80: -83.3}),
          ('eq_es: equatorial wedge, 60-88$^\\circ$', 60, 88, 'q > 0: E along the axis',
           {10: 3.1, 45: 14.5, 80: 22.2})]
INCS = (10, 45, 80)


def doit():
    '''
    Draw the cartoon
    '''

    fig, axes = plt.subplots(1, 2, figsize=(13, 6.6))
    for ax, (title, tmin, tmax, sense, q1) in zip(axes, MODELS):
        # the wind, in both hemispheres; theta is measured from the +z axis
        if tmin == 0:
            # a cone about the axis: one wedge in each hemisphere
            spans = [(90 - tmax, 90 + tmax), (270 - tmax, 270 + tmax)]
        else:
            spans = []
            for hemi in (1, -1):
                for side in (1, -1):
                    a1, a2 = (90 - tmax, 90 - tmin) if hemi == 1 else (-90 + tmin, -90 + tmax)
                    if side == -1:
                        a1, a2 = 180 - a2, 180 - a1
                    spans.append((a1, a2))
        for a1, a2 in spans:
            ax.add_patch(Wedge((0, 0), R_OUT, a1, a2, width=R_OUT - R_IN, facecolor='#9ecae1',
                               edgecolor='#3182bd', lw=0.8, alpha=0.8))
        # radial outflow arrows in the wind
        tc = np.radians(0.5 * (tmin + tmax))
        for hemi in (1, -1):
            for side in (1, -1):
                x0, z0 = side * np.sin(tc), hemi * np.cos(tc)
                ax.annotate('', xy=(7.5 * x0, 7.5 * z0), xytext=(4.0 * x0, 4.0 * z0),
                            arrowprops=dict(arrowstyle='->', color='#08519c', lw=1.2))
        # the disk: flat, absorbing, not radiating
        ax.plot([-R_DISK, R_DISK], [0, 0], color='0.25', lw=4, solid_capstyle='butt')
        ax.text(R_DISK + 0.2, -0.35, 'disk (absorbs,\nno emission)', fontsize=8, color='0.25', va='top')
        # the central source
        ax.add_patch(Circle((0, 0), 0.35, color='#fd8d3c', zorder=5))
        ax.text(0.3, 0.55, 'central source\n(10$^{37}$ erg s$^{-1}$, 10$^5$ K)', fontsize=8, color='#d94801')
        # observers, in the upper right quadrant, with the E-vector of the scattered light
        for inc in INCS:
            t = np.radians(inc)
            d = np.array([np.sin(t), np.cos(t)])
            p = 11.6 * d
            ax.annotate('', xy=10.9 * d, xytext=10.1 * d, arrowprops=dict(arrowstyle='->', color='k', lw=1))
            ax.plot([0, 10.1 * d[0]], [0, 10.1 * d[1]], ls=':', color='0.6', lw=0.7)
            ax.text(*(14.0 * d), '%d$^\\circ$\nq$_1$ = %+.0f%%' % (inc, q1[inc]), fontsize=8.5, ha='center', va='center')
            # E-vector of the scattered light on the sky: along the projected axis (q > 0) or perpendicular (q < 0)
            north = np.array([-np.cos(t), np.sin(t)])          # projected +z on the sky, in this plane
            e = north if q1[inc] > 0 else None
            if e is not None:
                ax.plot(*np.array([p - 0.45 * e, p + 0.45 * e]).T, color='#cb181d', lw=2.2)
            else:
                ax.plot(p[0], p[1], marker='o', ms=7, mfc='none', mec='#cb181d', mew=2)
        ax.set_title(title, fontsize=11)
        ax.text(0.02, 0.02, 'constant density, radial $\\tau_{es}$ = 1\nsimple atoms, no wind emission\n' + sense,
                transform=ax.transAxes, fontsize=9, va='bottom')
        ax.set_xlim(-11, 16.5)
        ax.set_ylim(-11, 15)
        ax.set_aspect('equal')
        ax.set_xlabel('x (10$^{10}$ cm)')
        ax.set_ylabel('z (10$^{10}$ cm)')
        ax.axhline(0, color='0.85', lw=0.5, zorder=0)
        ax.axvline(0, color='0.85', lw=0.5, zorder=0)
    fig.suptitle('The two CV-like wind models of the polarization tests (to scale; the wind extends from 2.8$\\times$10$^9$ '
                 'to 10$^{11}$ cm)\nred bar: E-vector in the plane of the figure (q > 0); red circle: E-vector perpendicular '
                 'to it, parallel to the disk (q < 0).  q$_1$: Sirocco, light scattered once', fontsize=9.5)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig('cvcone_cartoon.png', dpi=130, bbox_inches='tight')
    print('cvcone_cartoon.png')


def steer(argv):
    '''
    Process the command line; the routine takes no arguments except -h
    '''

    i = 1
    while i < len(argv):
        if argv[i] == '-h':
            print(__doc__)
            return
        else:
            print('Error: Unknown switch ---  %s' % argv[i])
            return
        i += 1

    doit()


# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    steer(sys.argv)
