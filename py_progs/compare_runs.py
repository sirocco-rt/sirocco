#!/usr/bin/env python
'''
Compare the ASCII outputs of two Sirocco runs, by default bitwise.

Synopsis:

    compare_runs.py [-tol x] [-files ext1,ext2,...] run1/root run2/root

Description:

    run1/root and run2/root are the root names of two runs, including the
    directory, e.g. shm/agn_macro and noshm/agn_macro.  The two runs may have
    different root names.

    By default the files compared are the spectra (.spec, .spec_wind,
    .log_spec, .log_spec_wind, .log_spec_tot, .log_spec_tot_wind) and any
    windsave2table outputs (root.*.txt) that exist for either run.  A file
    that exists for one run but not the other counts as a difference.

    Comment lines (starting with #) are ignored, since they hold the date,
    the version and the git hash.  Everything else must be identical, unless
    -tol is given, in which case numbers may differ by up to that relative
    amount.

    For each file that differs, the columns that differ are listed, with the
    number of rows affected and the largest relative difference.

Arguments:

    -tol x      Allow a relative difference of x between numbers (default 0,
                i.e. bitwise identical)
    -files ...  Comma-separated list of extensions to compare instead of the
                defaults, e.g. -files spec,log_spec_tot,master.txt
    -h          Print this help

Returns:

    Exit status 0 if all the files agree, 1 if any differ, 2 on a usage error.

Notes:

    Typical uses are checking that a change does not alter results, and that
    runs with different MPI layouts (e.g. with and without -no-shared-memory)
    give the same answer.  Runs with different numbers of MPI ranks use
    different random number streams, so they are not expected to agree.

    Only the python standard library is used.
'''

import sys
import os
import glob

SPEC_EXTENSIONS = ['spec', 'spec_wind', 'log_spec', 'log_spec_wind', 'log_spec_tot', 'log_spec_tot_wind']


def read_lines(filename):
    '''Return the non-comment, non-blank lines of a file, each split into words'''
    lines = []
    with open(filename) as f:
        for line in f:
            if line.startswith('#') or line.strip() == '':
                continue
            lines.append(line.split())
    return lines


def rel_diff(a, b):
    '''Relative difference between two words, or None if they are not both numbers'''
    try:
        x = float(a)
        y = float(b)
    except ValueError:
        return None
    if x == y:
        return 0.0
    return abs(x - y) / max(abs(x), abs(y))


def compare_file(file1, file2, tol):
    '''
    Compare two files, returning (agree, a list of lines describing the differences)
    '''
    lines1 = read_lines(file1)
    lines2 = read_lines(file2)

    if lines1 == lines2:
        return True, []

    if len(lines1) != len(lines2):
        return False, ['    different numbers of lines: %d and %d' % (len(lines1), len(lines2))]

    # Take column names from the first line if it is a header (any word not a number)
    header = None
    if lines1 and any(rel_diff(w, w) is None for w in lines1[0]):
        header = lines1[0]

    ndiff = {}  # column -> rows that differ
    maxrel = {}  # column -> largest relative difference (None if not numeric)
    nrow_shape = 0

    for row1, row2 in zip(lines1, lines2):
        if len(row1) != len(row2):
            nrow_shape += 1
            continue
        for j, (a, b) in enumerate(zip(row1, row2)):
            if a == b:
                continue
            r = rel_diff(a, b)
            if r is not None and r <= tol:
                continue
            ndiff[j] = ndiff.get(j, 0) + 1
            if r is None:
                maxrel[j] = None
            elif maxrel.get(j, 0.0) is not None:
                maxrel[j] = max(maxrel.get(j, 0.0), r)

    if not ndiff and nrow_shape == 0:
        return True, []

    report = []
    if nrow_shape:
        report.append('    %d rows have different numbers of columns' % nrow_shape)
    for j in sorted(ndiff):
        name = header[j] if header and j < len(header) else 'column %d' % (j + 1)
        rel = 'non-numeric' if maxrel[j] is None else 'max rel diff %.2e' % maxrel[j]
        report.append('    %-16s %6d of %d rows differ, %s' % (name, ndiff[j], len(lines1), rel))
    return False, report


def files_to_compare(root1, root2, extensions):
    '''Return a sorted list of the extensions to compare'''
    if extensions:
        return extensions
    exts = list(SPEC_EXTENSIONS)
    for root in (root1, root2):
        for f in glob.glob(root + '.*.txt'):
            exts.append(f[len(root) + 1:])
    return sorted(set(exts), key=lambda e: (e not in SPEC_EXTENSIONS, e))


def steer(argv):
    '''Parse the command line'''
    tol = 0.0
    extensions = None
    roots = []
    i = 1
    while i < len(argv):
        if argv[i] == '-h':
            print(__doc__)
            sys.exit(0)
        elif argv[i] == '-tol':
            i += 1
            tol = float(argv[i])
        elif argv[i] == '-files':
            i += 1
            extensions = [e.strip() for e in argv[i].split(',') if e.strip()]
        elif argv[i].startswith('-'):
            print('Unknown option %s' % argv[i])
            print(__doc__)
            sys.exit(2)
        else:
            roots.append(argv[i])
        i += 1

    if len(roots) != 2:
        print(__doc__)
        sys.exit(2)
    return roots[0], roots[1], tol, extensions


def main(argv):
    root1, root2, tol, extensions = steer(argv)

    what = 'bitwise' if tol == 0 else 'to a relative tolerance of %g' % tol
    print('Comparing %s and %s (%s, ignoring comment lines)' % (root1, root2, what))

    nfail = 0
    ncompared = 0
    for ext in files_to_compare(root1, root2, extensions):
        f1 = '%s.%s' % (root1, ext)
        f2 = '%s.%s' % (root2, ext)
        e1 = os.path.isfile(f1)
        e2 = os.path.isfile(f2)
        if not e1 and not e2:
            if extensions:
                print('  %-22s MISSING from both runs' % ext)
                nfail += 1
            continue
        if not (e1 and e2):
            print('  %-22s MISSING from %s' % (ext, root2 if e1 else root1))
            nfail += 1
            continue
        ncompared += 1
        agree, report = compare_file(f1, f2, tol)
        if agree:
            print('  %-22s identical' % ext)
        else:
            print('  %-22s DIFFERS' % ext)
            for line in report:
                print(line)
            nfail += 1

    if nfail:
        print('FAILED: %d file(s) differ or are missing' % nfail)
        return 1
    if ncompared == 0:
        print('FAILED: no files found for either run')
        return 1
    print('PASSED: %d file(s) agree' % ncompared)
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv))
