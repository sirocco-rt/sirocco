#!/usr/bin/env python
'''
Diagnostics for a make_cvcone.py run made with -split:

  - for each observer and band, the fraction of the flux that has been
    scattered, the polarization of the scattered light, and the total P
  - the continuum optical depth from the centre out along an inclination
    that passes through the wind (default 10 deg; needs sirocco_optd;
    written to optd/root.spec_tau)

Reads root.spec; prints the table; writes optd/root.spec_tau.

usage: analyze_split.py root [inclination]
'''
import os
import subprocess
import sys
import numpy as np

sys.path.insert(0, os.path.join(os.environ['SIROCCO'], 'py_progs'))
import plot_pol

root = sys.argv[1]
incl = sys.argv[2] if len(sys.argv) > 2 else '10'

units, freq, wave, cols = plot_pol.read_spec(root + '.spec')
fac = plot_pol.lnu_dnu(units, freq, False)
# bands: thirds of the wavelength range of the spectrum
edges = np.linspace(wave.min(), wave.max(), 4)
BANDS = tuple((round(edges[i]), round(edges[i + 1])) for i in range(3))
print('%s: frac_sc  P_sc(%%)  P_tot(%%) in bands %s' % (root, BANDS))
angles = sorted({int(k[1:].split('P')[0]) for k in cols if k.endswith('_sc:>1')})  # the observers in the file
for a in angles:
    n0, n1 = 'A%dP0.50_sc:0' % a, 'A%dP0.50_sc:>1' % a
    s = []
    for lo, hi in BANDS:
        m = (wave > lo) & (wave < hi)
        I0, I1 = (cols[n0] * fac)[m].sum(), (cols[n1] * fac)[m].sum()
        Q = (cols[n0 + '_Q'] * fac)[m].sum() + (cols[n1 + '_Q'] * fac)[m].sum()
        Q1 = (cols[n1 + '_Q'] * fac)[m].sum()
        s.append('%.3f %6.1f %6.2f' % (I1 / (I0 + I1), 100 * Q1 / I1, 100 * Q / (I0 + I1)))
    print('%2d  ' % a + '  |  '.join(s))

os.makedirs('optd', exist_ok=True)
for ext in ('.wind_save', '.pf'):
    dst = 'optd/' + root + ext
    if not os.path.exists(dst):
        os.symlink('../' + root + ext, dst)
if not os.path.exists('optd/data'):
    os.symlink('../data', 'optd/data')
out = subprocess.run([os.path.join(os.environ['SIROCCO'], 'bin', 'sirocco_optd'), '-i', incl, root], cwd='optd',
                     capture_output=True, text=True).stdout
for line in out.splitlines():
    if 'column density' in line or 'tau_H' in line:
        print(line.strip())
d = np.loadtxt('optd/' + root + '.spec_tau', comments='#', skiprows=6)
w, t = d[:, 1], d[:, 2]
print('continuum tau from the centre along %s deg:' % incl)
tedges = np.linspace(w.min(), w.max(), 7)                     # six windows across the range
print('  ' + '  '.join('%d-%d: %.2f' % (lo, hi, np.median(t[(w > lo) & (w < hi)]))
                       for lo, hi in zip(tedges[:-1], tedges[1:])))
