#!/usr/bin/env python

'''
Create a polar-coordinate import file from a windsave2table master file.

Read the master file produced by windsave2table for an rtheta (polar)
coordinate model and produce a file that can be imported into Sirocco.

Command line usage:
    import_rtheta.py rootname [outputfile]

    e.g.
        import_rtheta.py agn_rtheta_small
    reads  agn_rtheta_small.master.txt
    writes agn_rtheta_small.import.txt

Output columns (space-separated, readable by import_rtheta.c):
    i  j  inwind  r  theta  v_x  v_y  v_z  rho  t_e  t_r

Notes:
    windsave2table stores rho in the CMF (co-moving) frame.  This script
    converts it to the observer frame by multiplying by the Lorentz factor
    gamma = 1/sqrt(1 - v^2/c^2).

    For xmem3d two-hemisphere models theta runs from 0 to 180 deg and the
    output file covers both hemispheres.  import_rtheta.c detects the
    two-hemisphere case automatically when the last theta value exceeds 90 deg.

    For single-hemisphere (memory-branch) models theta runs from 0 to 90 deg
    and the file is identical in format.
'''

import sys
import numpy as np
from astropy.io import ascii
from astropy.table import Table


def doit(root='rtheta', outputfile=''):
    '''
    Read a master.txt file and write a Sirocco import file.

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

    C = 2.997925e10

    v2 = data['v_x']**2 + data['v_y']**2 + data['v_z']**2
    gamma = 1.0 / np.sqrt(1.0 - v2 / C**2)

    out = Table()
    out['i'] = data['i']
    out['j'] = data['j']
    out['inwind'] = data['inwind']
    out['r'] = data['r']
    out['theta'] = data['theta']
    out['v_x'] = data['v_x']
    out['v_y'] = data['v_y']
    out['v_z'] = data['v_z']
    out['rho'] = data['rho'] * gamma
    out['t_e'] = data['t_e']
    out['t_r'] = data['t_r']

    for col in ('r', 'theta', 'v_x', 'v_y', 'v_z', 'rho', 't_e', 't_r'):
        out[col].format = '.4e'

    ascii.write(out, outputfile, format='basic', overwrite=True)

    theta_max = float(np.max(data['theta']))
    hemi = 'two-hemisphere (theta 0-180)' if theta_max > 90 else 'single-hemisphere (theta 0-90)'
    print('Read  %d rows from %s  [%s]' % (len(out), filename, hemi))
    print('Wrote %d rows to   %s' % (len(out), outputfile))

    return out


if __name__ == '__main__':
    if len(sys.argv) > 2:
        doit(sys.argv[1], sys.argv[2])
    elif len(sys.argv) > 1:
        doit(sys.argv[1])
    else:
        print(__doc__)
