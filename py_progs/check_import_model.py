#!/usr/bin/env python

'''
Validate (and fix) a Sirocco import file.

Reports basic grid statistics, checks for guard-cell problems, and by
default writes a corrected copy with any fixable issues resolved.

Guard-cell violations (wrong inwind values at grid boundaries) are fixed
automatically.  Physical issues (rho <= 0, T <= 0, v >= c) are reported
but cannot be fixed automatically.

Command line usage::

    check_import_model.py filename [--no-fix]

    filename   path to an import file, or a rootname (tries <root>.import.txt)
    --no-fix   report only; do not write a corrected file

Output file (when fixes are applied):
    <name>.fixed.txt   e.g. polar.import.txt -> polar.import.fixed.txt
'''

import sys
import os
import numpy as np
from astropy.io import ascii
from astropy.table import Table


_DETECT = [
    ('sph3d',       {'r', 'theta', 'phi', 'k'}, set()),
    ('cyl3d',       {'x', 'z', 'phi', 'k'},     set()),
    ('polar',       {'r', 'theta', 'j'},         {'phi', 'k'}),
    ('cylindrical', {'x', 'z', 'j'},             {'phi', 'k'}),
    ('spherical',   {'r'},                       {'j', 'k'}),
]


def _detect_coord(colnames):
    cols = set(colnames)
    for name, required, forbidden in _DETECT:
        if required <= cols and not (forbidden & cols):
            return name
    return None


def _make_fixed_filename(filename):
    if filename.endswith('.txt'):
        return filename[:-4] + '.fixed.txt'
    return filename + '.fixed.txt'


def _check_guard(data, mask, label):
    '''Return (n_bad, label) for cells in mask that are not inwind=-1.'''
    bad = data[mask & (data['inwind'] != -1)]
    n = len(bad)
    if n == 0:
        print('  %-42s OK' % label)
    else:
        i_min, i_max = int(np.min(bad['i'])), int(np.max(bad['i']))
        if 'j' in bad.colnames:
            j_min, j_max = int(np.min(bad['j'])), int(np.max(bad['j']))
            print('  %-42s FAIL  %d cells  (i=%d..%d, j=%d..%d)' %
                  (label, n, i_min, i_max, j_min, j_max))
        else:
            print('  %-42s FAIL  %d cells  (i=%d..%d)' %
                  (label, n, i_min, i_max))
    return n


def _print_basic_info(data, coord):
    ndim = int(np.max(data['i'])) + 1
    mdim = int(np.max(data['j'])) + 1 if 'j' in data.colnames else None
    pdim = int(np.max(data['k'])) + 1 if 'k' in data.colnames else None

    print('Coord system     : %s  (detected from columns)' % coord)
    dims = 'ndim=%d' % ndim
    if mdim:
        dims += '  mdim=%d' % mdim
    if pdim:
        dims += '  pdim=%d' % pdim
    print('Grid             : %s' % dims)

    total     = len(data)
    n_inwind  = int(np.sum(data['inwind'] == 0))
    n_not     = int(np.sum(data['inwind'] == -1))
    n_part    = int(np.sum(data['inwind'] == -2))
    n_other   = total - n_inwind - n_not - n_part
    print('Total cells      : %d' % total)
    print('In wind          : %d' % n_inwind)
    print('Not in wind      : %d' % n_not)
    if n_part:
        print('Part in wind     : %d' % n_part)
    if n_other:
        print('Other inwind     : %d  (unexpected values)' % n_other)

    return ndim, mdim, pdim


def _two_hemi(data, coord):
    if coord in ('polar', 'sph3d'):
        return float(np.max(data['theta'])) > 90.0
    if coord in ('cylindrical', 'cyl3d'):
        return float(np.min(data['z'])) < 0.0
    return False


def _guard_masks(data, coord, ndim, mdim, two_hemi_flag):
    '''Return list of (mask, label) pairs for guard-cell boundaries.'''
    guards = []
    if coord == 'spherical':
        guards.append((data['i'] == ndim - 1,
                       'Outer radial (i=%d)' % (ndim - 1)))
    else:
        guards.append((data['i'] >= ndim - 2,
                       'Outer radial (i=%d,%d)' % (ndim - 2, ndim - 1)))
        if coord in ('polar', 'sph3d') and mdim is not None:
            guards.append((data['j'] <= 1, 'North pole (j=0,1)'))
            if two_hemi_flag:
                guards.append((data['j'] >= mdim - 2,
                               'South pole (j=%d,%d)' % (mdim - 2, mdim - 1)))
    return guards


def _run_guard_checks(data, coord, ndim, mdim, two_hemi_flag):
    guards = _guard_masks(data, coord, ndim, mdim, two_hemi_flag)
    total_bad = 0
    for mask, label in guards:
        total_bad += _check_guard(data, mask, label)
    return total_bad


def _run_physical_checks(data):
    C = 2.997925e10
    v2 = data['v_x']**2 + data['v_y']**2 + data['v_z']**2
    v_over_c = float(np.sqrt(np.max(v2))) / C
    flag = '  WARNING: v/c >= 1!' if v_over_c >= 1.0 else 'OK'
    print('  Max v/c                                    %.4f  %s' % (v_over_c, flag))

    n_inwind = int(np.sum(data['inwind'] == 0))
    phys_problems = 0
    if n_inwind > 0:
        in_mask = data['inwind'] == 0
        rho = data['rho'][in_mask]
        n_zero_rho = int(np.sum(rho <= 0.0))
        flag = '  WARNING: %d cells have rho <= 0' % n_zero_rho if n_zero_rho else 'OK'
        print('  Min rho (in-wind)                          %.3e  %s' %
              (float(np.min(rho)), flag))
        phys_problems += n_zero_rho

        te = data['t_e'][in_mask]
        tr = data['t_r'][in_mask]
        n_cold = int(np.sum(te <= 0.0)) + int(np.sum(tr <= 0.0))
        flag = '  WARNING: %d cells have T <= 0' % n_cold if n_cold else 'OK'
        print('  Min t_e / t_r (in-wind)                    %.3e / %.3e  %s' %
              (float(np.min(te)), float(np.min(tr)), flag))
        phys_problems += n_cold

    return phys_problems


def check(filename, apply_fix=True):
    '''
    Validate (and optionally fix) a Sirocco import file.

    Parameters
    ----------
    filename : str
        Path to the import file.
    apply_fix : bool
        If True (default) write a corrected copy when guard-cell problems
        are found.

    Returns
    -------
    int
        0 if the file is clean (or was successfully fixed), 1 otherwise.
    '''

    try:
        data = ascii.read(filename)
    except Exception as e:
        print('Error reading %s: %s' % (filename, e))
        return 1

    coord = _detect_coord(data.colnames)
    if coord is None:
        print('Error: cannot determine coordinate system from columns:')
        print('  %s' % ', '.join(data.colnames))
        return 1

    print('File             : %s' % filename)
    ndim, mdim, pdim = _print_basic_info(data, coord)

    th = _two_hemi(data, coord)
    if coord in ('polar', 'sph3d'):
        print('Hemisphere       : %s' %
              ('two (theta 0-180)' if th else 'single (theta 0-90)'))
    elif coord in ('cylindrical', 'cyl3d'):
        print('Hemisphere       : %s' %
              ('two (z < 0 to z > 0)' if th else 'single (z >= 0)'))

    print('\nGuard cell checks:')
    n_guard_bad = _run_guard_checks(data, coord, ndim, mdim, th)

    print('\nPhysical checks:')
    n_phys_bad = _run_physical_checks(data)

    # ---- Fix and write ----------------------------------------------------
    if n_guard_bad > 0 and apply_fix:
        out = Table(data)           # copy

        guards = _guard_masks(out, coord, ndim, mdim, th)
        n_fixed = 0
        for mask, _ in guards:
            bad = mask & (out['inwind'] != -1)
            n_fixed += int(np.sum(bad))
            out['inwind'][bad] = -1

        fixed_filename = _make_fixed_filename(filename)

        # Preserve numeric formats from original columns
        float_cols = [c for c in out.colnames
                      if out[c].dtype.kind == 'f' and c != 'inwind']
        for col in float_cols:
            out[col].format = '.4e'

        ascii.write(out, fixed_filename, format='basic', overwrite=True)

        print('\n--- Fixed file: %s ---' % fixed_filename)
        print('Guard cells corrected: %d cells set to inwind=-1' % n_fixed)
        print('\nGuard cell checks (fixed file):')
        _run_guard_checks(out, coord, ndim, mdim, th)

        n_inwind_fixed = int(np.sum(out['inwind'] == 0))
        n_not_fixed    = int(np.sum(out['inwind'] == -1))
        print('\nIn wind (fixed)  : %d' % n_inwind_fixed)
        print('Not in wind      : %d' % n_not_fixed)

    # ---- Final summary ----------------------------------------------------
    print()
    if n_guard_bad == 0 and n_phys_bad == 0:
        print('Summary: CLEAN')
        return 0
    elif n_guard_bad > 0 and apply_fix:
        if n_phys_bad > 0:
            print('Summary: guard cells fixed; physical issues remain (see above)')
        else:
            print('Summary: guard cells fixed and written to %s' %
                  _make_fixed_filename(filename))
        return 0
    else:
        parts = []
        if n_guard_bad > 0:
            parts.append('%d guard-cell problem(s)' % n_guard_bad)
        if n_phys_bad > 0:
            parts.append('%d physical issue(s)' % n_phys_bad)
        print('Summary: %s found' % ' and '.join(parts))
        return 1


def doit(arg, apply_fix=True):
    if os.path.exists(arg):
        return check(arg, apply_fix)
    candidate = arg + '.import.txt'
    if os.path.exists(candidate):
        return check(candidate, apply_fix)
    print('Error: cannot find %s or %s' % (arg, candidate))
    return 1


if __name__ == '__main__':
    args = [a for a in sys.argv[1:] if not a.startswith('-')]
    flags = [a for a in sys.argv[1:] if a.startswith('-')]

    if not args:
        print(__doc__)
        sys.exit(0)

    apply_fix = '--no-fix' not in flags
    rc = doit(args[0], apply_fix)
    sys.exit(0 if rc == 0 else 1)
