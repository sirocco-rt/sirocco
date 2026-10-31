'''
Test definition for pol_tests.py: the Hillier (1994) prolate shell, the test
used by Bulla et al. (2015, Sect. 3.1).

Hillier's values (hillier_ref.py) were read by eye from his figures, to
about +-0.05%, so the tolerances are absolute differences in q (%):

    q at each inclination (hx*): within 0.2% of Hillier.  The runs this was
    set from differ by up to 0.08% for chi0 <= 0.16 and 0.56, 0.11% for
    chi0 = 1 and 0.16% for chi0 = 0.32 (probably the reading of his Fig 3).

    q against albedo at 90 deg (hal90): within 0.1% (runs: up to 0.03%).

    q against albedo at 56.25 deg (hal56): within 0.3% (runs: up to 0.22%).
    This is where Hillier's own Monte Carlo and transfer results differ
    most, by about 0.1%.

    max |u| < 0.1%: U should be 0 for an axisymmetric model.

The quick tier runs only chi0 = 0.1 (tau_ave = 0.29), with fewer photons.
'''
import io
import contextlib

DESCRIPTION = 'Hillier (1994) prolate shell: polarization vs optical depth and albedo'
CHI0_FULL = [0.01, 0.0316, 0.056, 0.1, 0.16, 0.32, 0.56, 1.0]


def models(tier):
    if tier == 'quick':
        return [dict(make=['make_hillier.py', 'hx0.1', '0.1', '-nphot', '1000000', '-ncycles', '2'], root='hx0.1')]
    out = [dict(make=['make_hillier.py', 'hal90', '1.0', '-nscat', '-angle', '90', '-nphot', '2000000', '-ncycles', '2'],
                root='hal90'),
           dict(make=['make_hillier.py', 'hal56', '1.0', '-nscat', '-angle', '56.25', '-nphot', '2000000', '-ncycles', '2'],
                root='hal56')]
    for c in CHI0_FULL:
        out.append(dict(make=['make_hillier.py', 'hx%s' % c, str(c), '-nphot', '2000000', '-ncycles', '2'], root='hx%s' % c))
    return out


def post(tier):
    return [] if tier == 'quick' else [['compare_hillier.py']]


def checks(tier):
    import hillier_pol as hp
    import hillier_ref as ref

    out = []
    for c in ([0.1] if tier == 'quick' else CHI0_FULL):
        b = hp.band('hx%s.spec' % c)
        q = [100 * Q / I for o, I, Q, U in b]
        u = [100 * U / I for o, I, Q, U in b]
        d = [abs(x - y) for x, y in zip(q, ref.P_INC[c]) if y is not None]
        out.append(('hx%s (tau_ave %.3f): max |q - Hillier| (%%)' % (c, 2.888 * c), max(d), 0.0, 0.2))
        out.append(('hx%s: max |u| (%%)' % c, max(abs(x) for x in u), 0.0, 0.1))
    if tier == 'full':
        for root, ang, tol in (('hal90', 90.0, 0.1), ('hal56', 56.25, 0.3)):
            with contextlib.redirect_stdout(io.StringIO()):
                res = hp.albedo(root)
            q = {round(a, 2): x for a, x, y in res}
            d = [abs(q[round(a, 2)] - p) for a, p in zip(ref.ALB, ref.P_ALB[ang])]
            out.append(('%s: max |q - Hillier| over albedo 0.1-1 (%%)' % root, max(d), 0.0, tol))
    return out
