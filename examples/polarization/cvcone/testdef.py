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

    The wedge at 80 deg: within 1 point (runs: 0.2).

    The wedge at 10 and 45 deg is reported but not checked: Sirocco is about 2
    points more positive than the prediction there, which is not understood.

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
            label = 'eq_es (wedge) i=%d: q - prediction, scattered once (points; q = %.1f%%)' % (row['inc'], row['q1'])
            if row['inc'] == 80:
                out.append((label, row['q1'] - row['pred1'], -1.0, 1.0))
            else:
                out.append((label + ' [known discrepancy, not checked]', row['q1'] - row['pred1'], None, None))
    return out
