#!/usr/bin/env python

'''
Create a SPH3D import file from a windsave2table master file.

Read the master file produced by windsave2table for a SPH3D model
and write a file that can be imported back into Sirocco.

Command line usage::

    import_sph3d.py rootname [outputfile]

    e.g.
        import_sph3d.py cv_sph3d
        reads  cv_sph3d.master.txt
        writes cv_sph3d.import.txt

Output columns (space-separated, readable by import_sph3d.c)::

    i  j  k  inwind  r  theta  phi  v_x  v_y  v_z  rho  t_e  t_r

Column descriptions:

- ``r``       inner radial boundary of the cell (cm)
- ``theta``   inner polar-angle boundary (degrees, 0–180)
- ``phi``     inner azimuthal boundary (degrees, 0–360)
- ``v_x, v_y, v_z``  Cartesian velocity components (cm/s)
- ``rho``     observer-frame mass density (CMF rho * Lorentz gamma, g/cm³)
- ``t_e, t_r``  electron and radiation temperatures (K)

The SPH3D master.txt corner coordinates are labelled r, theta, phi
(distinct from the cell-centre columns rcen, thetacen, phicen).
'''

import sys
import numpy as np
from astropy.io import ascii
from astropy.table import Table


def doit(root='model', outputfile=''):
    '''
    Read an SPH3D master.txt and write a Sirocco 3D spherical polar import file.

    Parameters
    ----------
    root : str
        Rootname.  Reads <root>.master.txt.
    outputfile : str, optional
        Output path.  Defaults to <root>.import.txt.

    Returns
    -------
    astropy.table.Table
        The table that was written, or None on error.
    '''

    filename = root + '.master.txt'
    if outputfile == '':
        outputfile = root + '.import.txt'

    try:
        data = ascii.read(filename)
    except Exception as e:
        print('Error reading %s: %s' % (filename, e))
        return None

    # Verify this looks like a SPH3D master file
    required = {'r', 'theta', 'phi', 'i', 'j', 'k', 'inwind',
                'v_x', 'v_y', 'v_z', 'rho', 't_e', 't_r'}
    missing = required - set(data.colnames)
    if missing:
        print('Error: master file is missing columns: %s' % ', '.join(sorted(missing)))
        print('Expected an SPH3D master.txt; got columns: %s' % ', '.join(data.colnames))
        return None

    C = 2.997925e10

    v2 = data['v_x']**2 + data['v_y']**2 + data['v_z']**2
    gamma = 1.0 / np.sqrt(1.0 - v2 / C**2)

    out = Table()
    out['i'] = data['i']
    out['j'] = data['j']
    out['k'] = data['k']
    out['inwind'] = data['inwind']
    out['r'] = data['r']
    out['theta'] = data['theta']
    out['phi'] = data['phi']
    out['v_x'] = data['v_x']
    out['v_y'] = data['v_y']
    out['v_z'] = data['v_z']
    out['rho'] = data['rho'] * gamma
    out['t_e'] = data['t_e']
    out['t_r'] = data['t_r']

    for col in ('r', 'theta', 'phi', 'v_x', 'v_y', 'v_z', 'rho', 't_e', 't_r'):
        out[col].format = '.4e'

    ascii.write(out, outputfile, format='basic', overwrite=True)

    ndim = int(np.max(data['i'])) + 1
    mdim = int(np.max(data['j'])) + 1
    pdim = int(np.max(data['k'])) + 1
    print('Grid dimensions: ndim=%d  mdim=%d  pdim=%d  (total %d cells)' % (ndim, mdim, pdim, len(out)))
    print('Read  %d rows from %s' % (len(out), filename))
    print('Wrote %d rows to   %s' % (len(out), outputfile))

    return out


if __name__ == '__main__':
    if len(sys.argv) > 2:
        doit(sys.argv[1], sys.argv[2])
    elif len(sys.argv) > 1:
        doit(sys.argv[1])
    else:
        print(__doc__)
