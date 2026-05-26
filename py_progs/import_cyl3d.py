#!/usr/bin/env python

'''
Create a CYLIND3D import file from a windsave2table master file.

Read the master file produced by windsave2table for a CYLIND3D model
and write a file that can be imported back into Sirocco.

Command line usage::

    import_cyl3d.py rootname [outputfile]

    e.g.
        import_cyl3d.py agn_cyl3d
        reads  agn_cyl3d.master.txt
        writes agn_cyl3d.import.txt

Output columns (space-separated, readable by import_cylindrical3d.c)::

    i  j  k  inwind  x  z  phi  v_x  v_y  v_z  rho  t_e  t_r

where x and z are in cm, phi is the azimuthal lower boundary in degrees (0–360),
rho is the observer-frame mass density (CMF rho * Lorentz gamma), and
t_e, t_r are the electron and radiation temperatures.

The CYLIND3D master.txt corner coordinates are labelled x, z, phi
(distinct from the cell-centre columns xcen, zcen, phicen).
'''

import sys
import numpy as np
from astropy.io import ascii
from astropy.table import Table


def doit(root='model', outputfile=''):
    '''
    Read a CYLIND3D master.txt and write a Sirocco 3D cylindrical import file.

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

    # Verify this looks like a CYLIND3D master file
    required = {'x', 'z', 'phi', 'i', 'j', 'k', 'inwind', 'v_x', 'v_y', 'v_z', 'rho', 't_e', 't_r'}
    missing = required - set(data.colnames)
    if missing:
        print('Error: master file is missing columns: %s' % ', '.join(sorted(missing)))
        print('Expected a CYLIND3D master.txt; got columns: %s' % ', '.join(data.colnames))
        return None

    C = 2.997925e10

    v2 = data['v_x']**2 + data['v_y']**2 + data['v_z']**2
    gamma = 1.0 / np.sqrt(1.0 - v2 / C**2)

    out = Table()
    out['i'] = data['i']
    out['j'] = data['j']
    out['k'] = data['k']
    out['inwind'] = data['inwind']
    out['x'] = data['x']
    out['z'] = data['z']
    out['phi'] = data['phi']
    out['v_x'] = data['v_x']
    out['v_y'] = data['v_y']
    out['v_z'] = data['v_z']
    out['rho'] = data['rho'] * gamma
    out['t_e'] = data['t_e']
    out['t_r'] = data['t_r']

    for col in ('x', 'z', 'phi', 'v_x', 'v_y', 'v_z', 'rho', 't_e', 't_r'):
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
