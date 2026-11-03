'''
Test definition for pol_tests.py: an optically thin prolate shell moving
radially, compared with an exact special-relativistic single-scattering
calculation (sr_single.py).

The prediction leaves out the attenuation of the scattered light on its way
out (tau of 0.2-2%), so Sirocco should be slightly low.  Checks, on the
bolometric light scattered exactly once, relative to the direct light:

    I_sc ratio (Sirocco / prediction) at every angle, and Q_sc ratio at 45,
    67.5 and 90 deg, between 0.96 and 1.01.  (At 22.5 deg Q is too small to
    compare.)  The runs this was set from gave 0.971-0.999.  The -nonrel run
    is compared with the static prediction.

    The spherical shell: max |q_sc| < 0.3% (runs: 0.11%).

The quick tier runs only the shell at 0.3c, with a quarter of the photons;
there the noise in each ratio is about 0.7%, so the range is 0.95-1.03 (a
quick run gave 0.989-1.001).
'''
import numpy as np

DESCRIPTION = 'Thin shell moving at up to 0.3c vs exact relativistic single scattering'
FULL = [('sb0', '3e-4', []), ('sb0.1', '0.1', []), ('sb0.3', '0.3', []), ('sb0.3n', '0.3', ['-nonrel']),
        ('ss0.3', '0.3', [])]


def models(tier):
    if tier == 'quick':
        return [dict(make=['make_srshell.py', 'sb0.3', '0.003', '0.3', '-nphot', '1000000'], root='sb0.3')]
    out = []
    for root, beta, sw in FULL:
        make = ['make_srshell.py', root, '0.003', beta] + (['-sphere'] if root.startswith('ss') else [])
        out.append(dict(make=make, root=root, switches=sw))
    return out


def post(tier):
    return [['compare_sr.py', 'sb0.3']] if tier == 'quick' else [['compare_sr.py'] + [root for root, beta, sw in FULL]]


def checks(tier):
    from compare_sr import compare

    lo, hi = (0.95, 1.03) if tier == 'quick' else (0.96, 1.01)
    out = []
    for root in (['sb0.3'] if tier == 'quick' else [r[0] for r in FULL]):
        beta, nonrel, rows = compare(root)
        label = '%s (beta %g%s)' % (root, beta, ', -nonrel' if nonrel else '')
        if root.startswith('ss'):
            out.append(('%s sphere: max |q_sc| (%%)' % label, max(abs(r['q']) for r in rows), 0.0, 0.3))
            continue
        ri = [r['ratio_i'] for r in rows]
        rq = [r['ratio_q'] for r in rows if r['inc'] >= 45]
        out.append(('%s: I_sc / prediction, lowest' % label, min(ri), lo, hi))
        out.append(('%s: I_sc / prediction, highest' % label, max(ri), lo, hi))
        out.append(('%s: Q_sc / prediction (45-90 deg), lowest' % label, np.nanmin(rq), lo, hi))
        out.append(('%s: Q_sc / prediction (45-90 deg), highest' % label, np.nanmax(rq), lo, hi))
    return out
