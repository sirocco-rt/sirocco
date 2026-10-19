'''
Test definition for pol_tests.py: thin, flattened electron-scattering
envelopes around a point source, compared with single scattering
(single_scatter.py), which is Brown & McLean (1977) with attenuation.

q is measured over 4500-6500 A, away from the edges of the band, where the
hot (1e6 K) electrons Doppler-shift scattered light in or out of the band.

Checks, at 50, 70 and 90 deg (q is small at lower inclinations and noisy):

    q / q(single scattering) between 0.96 and 1 + tau_eq + 0.04.
    Multiple scattering only adds polarized light, by about 0.8 tau_eq, so
    the ratio should be at least 1; 0.96 allows for noise.  The runs this
    was set from gave 1.09 (tau_eq 0.1), 1.02-1.03 (0.03) and 0.97-0.985
    (0.01): the thinnest envelope is a few per cent low, close to the limit.

    max |u| / q(90) < 0.02: U should be 0 for an axisymmetric model.
'''
DESCRIPTION = 'Thin envelopes around a point source vs Brown & McLean (1977) with attenuation'
CASES = {'quick': [('bm01', 0.1, 4000000)],
         'full': [('bm01', 0.1, 8000000), ('bm003', 0.03, 8000000), ('bm001', 0.01, 16000000)]}


def models(tier):
    return [dict(make=['make_bm.py', root, str(tau), '-nphot', str(nphot)], root=root)
            for root, tau, nphot in CASES[tier]]


def post(tier):
    return [['compare_bm.py'] + ['%s:%g' % (root, tau) for root, tau, nphot in CASES[tier]]]


def checks(tier):
    from compare_bm import compare

    out = []
    for root, tau, nphot in CASES[tier]:
        rows = {r['inc']: r for r in compare(root, tau)}
        for inc in (50, 70, 90):
            out.append(('%s (tau_eq %g) i=%d: q / single scattering' % (root, tau, inc), rows[inc]['ratio'],
                        0.96, 1 + tau + 0.04))
        out.append(('%s: max |u| / q(90)' % root, max(r['u'] for r in rows.values()) / rows[90]['q'], 0.0, 0.02))
    return out
