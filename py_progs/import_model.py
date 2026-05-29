#!/usr/bin/env python

'''
Create a Sirocco import file from a windsave2table master file.

Reads the coord_system comment written by windsave2table and produces
the appropriate import file automatically.  Supports all five coordinate
systems: spherical, cylindrical, polar, cyl3d, sph3d.

Command line usage::

    import_model.py rootname [outputfile] [--keep-partial]

    e.g.
        import_model.py cv_polar
        reads  cv_polar.master.txt  (detects: coord_system: polar)
        writes cv_polar.import.txt

Old master.txt files without a coord_system header are handled via
column-based detection with a printed warning.

By default, cells with ``inwind = 1`` (W_PART_INWIND, partially in the
wind) are converted to ``inwind = -1`` before writing.  Sirocco's import
code converts both W_NOT_INWIND (−1) and W_PART_INWIND (1) to W_IGNORE
(−2) on import anyway, so this produces a cleaner file.  Pass
``--keep-partial`` on the command line (or ``keep_partial=True`` when
calling :func:`doit`) to pass W_PART_INWIND cells through unchanged.

After writing the import file the script automatically runs
check_import_model.py and prints a validation summary.
'''

import sys
import numpy as np
from astropy.io import ascii
from astropy.table import Table
from check_import_model import check as _check_output


def _read_coord_system(filename):
    '''Return the coord_system string from the first line of a master.txt, or None.'''
    try:
        with open(filename) as f:
            line = f.readline().strip()
        if line.startswith('# coord_system:'):
            return line.split(':', 1)[1].strip()
    except IOError:
        pass
    return None


def _gamma(data):
    C = 2.997925e10
    v2 = data['v_x']**2 + data['v_y']**2 + data['v_z']**2
    return 1.0 / np.sqrt(1.0 - v2 / C**2)


def _doit_spherical(data, outputfile):
    C = 2.997925e10
    v = np.sqrt(data['v_x']**2 + data['v_y']**2 + data['v_z']**2)
    gamma = 1.0 / np.sqrt(1.0 - (v / C)**2)

    out = Table()
    out['i'] = data['i']
    out['inwind'] = data['inwind']
    out['r'] = data['r']
    out['v'] = v
    out['rho'] = data['rho'] * gamma
    out['t_e'] = data['t_e']
    out['t_r'] = data['t_r']

    for col in ('r', 'v', 'rho', 't_e', 't_r'):
        out[col].format = '.4e'

    ascii.write(out, outputfile, format='fixed_width_two_line', overwrite=True)
    return out


def _doit_cylindrical(data, outputfile):
    gamma = _gamma(data)

    out = Table()
    out['i'] = data['i']
    out['j'] = data['j']
    out['inwind'] = data['inwind']
    out['x'] = data['x']
    out['z'] = data['z']
    out['v_x'] = data['v_x']
    out['v_y'] = data['v_y']
    out['v_z'] = data['v_z']
    out['rho'] = data['rho'] * gamma
    out['t_e'] = data['t_e']
    out['t_r'] = data['t_r']

    for col in ('x', 'z', 'v_x', 'v_y', 'v_z', 'rho', 't_e', 't_r'):
        out[col].format = '.4e'

    ascii.write(out, outputfile, format='basic', overwrite=True)
    return out


def _doit_polar(data, outputfile):
    gamma = _gamma(data)

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

    ndim = int(np.max(out['i'])) + 1
    mdim = int(np.max(out['j'])) + 1
    theta_max = float(np.max(data['theta']))
    two_hemi = theta_max > 90.0

    guard_mask = (
        (out['i'] >= ndim - 2) |
        (out['j'] <= 1) |
        (two_hemi & (out['j'] >= mdim - 2))
    )
    n_forced = int(np.sum(guard_mask & (out['inwind'] != -1)))
    if n_forced:
        out['inwind'][guard_mask] = -1
        print('Forced %d boundary cells to inwind=-1 (guard cells)' % n_forced)

    for col in ('r', 'theta', 'v_x', 'v_y', 'v_z', 'rho', 't_e', 't_r'):
        out[col].format = '.4e'

    ascii.write(out, outputfile, format='basic', overwrite=True)
    return out


def _doit_cyl3d(data, outputfile):
    gamma = _gamma(data)

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
    return out


def _doit_sph3d(data, outputfile):
    gamma = _gamma(data)

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

    ndim = int(np.max(out['i'])) + 1
    mdim = int(np.max(out['j'])) + 1
    pdim = int(np.max(out['k'])) + 1
    theta_max = float(np.max(data['theta']))
    two_hemi = theta_max > 90.0

    guard_mask = (
        (out['i'] >= ndim - 2) |
        (out['j'] <= 1) |
        (two_hemi & (out['j'] >= mdim - 2))
    )
    n_forced = int(np.sum(guard_mask & (out['inwind'] != -1)))
    if n_forced:
        out['inwind'][guard_mask] = -1
        print('Forced %d boundary cells to inwind=-1 (guard cells)' % n_forced)

    for col in ('r', 'theta', 'phi', 'v_x', 'v_y', 'v_z', 'rho', 't_e', 't_r'):
        out[col].format = '.4e'

    ascii.write(out, outputfile, format='basic', overwrite=True)
    return out


_DISPATCH = {
    'spherical':   _doit_spherical,
    'cylindrical': _doit_cylindrical,
    'polar':       _doit_polar,
    'cyl3d':       _doit_cyl3d,
    'sph3d':       _doit_sph3d,
}

# Fallback detection for old master.txt files without a coord_system header.
# Each entry: (name, required_cols, forbidden_cols).  Checked in order;
# more-specific entries (3D) must come before less-specific ones (2D/1D).
_FALLBACK = [
    ('sph3d',       {'r', 'theta', 'phi', 'k'}, set()),
    ('cyl3d',       {'x', 'z', 'phi', 'k'},     set()),
    ('polar',       {'r', 'theta', 'j'},         {'phi', 'k'}),
    ('cylindrical', {'x', 'z', 'j'},             {'phi', 'k'}),
    ('spherical',   {'r'},                       {'j', 'k'}),
]


def _detect_coord(colnames):
    cols = set(colnames)
    for name, required, forbidden in _FALLBACK:
        if required <= cols and not (forbidden & cols):
            return name
    return None


def doit(root='model', outputfile='', keep_partial=False):
    '''
    Read a master.txt and write a Sirocco import file.

    Parameters
    ----------
    root : str
        Rootname.  Reads <root>.master.txt.
    outputfile : str, optional
        Output path.  Defaults to <root>.import.txt.
    keep_partial : bool, optional
        If False (default), W_PART_INWIND cells (inwind=1) are converted to
        inwind=-1 before writing.  Sirocco converts them to W_IGNORE on import
        anyway, so the default produces a cleaner file.  Set True to pass them
        through unchanged.

    Returns
    -------
    astropy.table.Table or None on error.
    '''

    filename = root + '.master.txt'
    if outputfile == '':
        outputfile = root + '.import.txt'

    coord = _read_coord_system(filename)

    try:
        data = ascii.read(filename)
    except Exception as e:
        print('Error reading %s: %s' % (filename, e))
        return None

    if coord is None:
        coord = _detect_coord(data.colnames)
        if coord is None:
            print('Error: cannot determine coordinate system from %s' % filename)
            print('Columns: %s' % ', '.join(data.colnames))
            return None
        print('Warning: no coord_system header; detected %s from columns' % coord)

    fn = _DISPATCH.get(coord)
    if fn is None:
        print('Error: unknown coord_system "%s"' % coord)
        return None

    if not keep_partial:
        n_partial = int(np.sum(data['inwind'] == 1))
        if n_partial:
            data['inwind'][data['inwind'] == 1] = -1
            print('Converted %d W_PART_INWIND (inwind=1) cells to inwind=-1' % n_partial)

    out = fn(data, outputfile)
    print('Read  %d rows from %s' % (len(out), filename))
    print()
    _check_output(outputfile, apply_fix=False)
    return out


if __name__ == '__main__':
    args = [a for a in sys.argv[1:] if not a.startswith('-')]
    flags = [a for a in sys.argv[1:] if a.startswith('-')]
    keep_partial = '--keep-partial' in flags

    if len(args) > 1:
        doit(args[0], args[1], keep_partial=keep_partial)
    elif len(args) == 1:
        doit(args[0], keep_partial=keep_partial)
    else:
        print(__doc__)
