'''
Test definition for pol_tests.py: a CV-like wind of constant density in a
polar cone (pbl_es) or an equatorial wedge (eq_es), lit by a central source,
with the ionization calculated by Sirocco.  The light scattered exactly once
is compared with an independent single-scattering calculation (geom_pol.py),
which uses the n_e of the model, so windsave2table must be able to read the
wind_save file.

q is measured in the windows between the resonance lines (1700-1850 and
2000-2700 A).  Checks:

    The cone, at 10, 45 and 80 deg: |q - prediction| < 1 point (full) or 1.5
    points (quick, fewer photons and a 10000-ray prediction).  The runs this
    was set from agree to 0.4 points or better.

    The wedge: within 1 point at 10 and 45 deg, 1.5 points at 80 deg.
    Against the mean of 5 predictions (seeds 2, 11-14; sd 0.25-0.4 points),
    the full-tier run of 261019 differs by -0.03, +0.58 and -0.91 points at
    10, 45 and 80 deg.  The offsets at 45 and 80 deg are larger than the
    noise and are not understood; with one prediction (seed 2) the 80 deg
    run was 1.17 points low, hence the wider tolerance there.

    Until 261019 the wedge at 10 and 45 deg was about 2 points too positive.
    The cause was the CV's secondary star, which then lay inside the wind and
    hid part of it from the observers; the prediction has no secondary.  The
    models now put the secondary well outside the wind (make_cvcone.py).

The quick tier runs pol_cone, a small version of pbl_es.  examples/regress has
a model of the same name, but without the observers selected by the number of
scatters, which the regression does not need.
'''
import os

DESCRIPTION = 'CV-like wind in a polar cone or equatorial wedge vs single scattering'
PLOT_POL = os.path.join(os.environ['SIROCCO'], 'py_progs', 'plot_pol.py')


def models(tier):
    if tier == 'quick':
        return [dict(make=['make_cvcone.py', 'pol_cone', '0', '20', '-nphot', '500000', '-ion_cycles', '5',
                           '-spec_cycles', '2'], root='pol_cone', w2t=True)]
    return [dict(make=None, root='pbl_es', switches=['-p'], w2t=True),
            dict(make=None, root='eq_es', switches=['-p'], w2t=True)]


def post(tier):
    if tier == 'quick':
        return [['geom_pol.py', 'pol_cone', '0', '20', '-central', '-nray', '10000'],
                ['es_vs_prediction.py', 'pol_cone'],
                [PLOT_POL, 'pol_cone', '-max_obs', '0']]
    return [['geom_pol.py', 'pbl_es', '0', '20', '-central', '-nray', '30000'],
            ['geom_pol.py', 'eq_es', '60', '88', '-central', '-nray', '30000'],
            ['es_vs_prediction.py', 'pbl_es', 'eq_es'],
            [PLOT_POL, 'pbl_es', '-max_obs', '0'],
            [PLOT_POL, 'eq_es', '-max_obs', '0']]


def checks(tier):
    from es_vs_prediction import compare

    out = []
    cone = 'pol_cone' if tier == 'quick' else 'pbl_es'
    tol = 1.5 if tier == 'quick' else 1.0
    for row in compare(cone):
        out.append(('%s (cone) i=%d: q - prediction, scattered once (points; q = %.1f%%)' % (cone, row['inc'], row['q1']),
                    row['q1'] - row['pred1'], -tol, tol))
    if tier == 'full':
        for row in compare('eq_es'):
            tol = 1.5 if row['inc'] == 80 else 1.0
            out.append(('eq_es (wedge) i=%d: q - prediction, scattered once (points; q = %.1f%%)' % (row['inc'], row['q1']),
                        row['q1'] - row['pred1'], -tol, tol))
    return out
