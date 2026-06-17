#!/usr/bin/env python
"""
ConvertMacro2Simple.py -- Convert a Sirocco macro-atom atomic data set to simple-atom format.

The script reads a macro-atom master data file (e.g. data/h20_hetop_standard80.dat),
identifies the component files that contain macro-atom specific keywords (IonM,
LevMacro, LinMacro), converts them to their simple-atom equivalents, and writes
a new master file that assembles the converted and unchanged files.

Conversions applied:
    IonM   -> IonV  (elem_ions files)   nlte field set to 0
    LevMacro -> LevTop (levels files)   -1 inserted for islp and qqnum fields
    LinMacro -> Line   (lines files)    filtered by lower level index and oscillator strength

Files that require no conversion (topbase phot files, collision data, etc.) are
referenced in-place in the new master file -- no copying.

Usage:
    macro2simple.py input_master.dat [options]

    input_master.dat  Path to macro-atom master file, relative to $SIROCCO
                      (e.g. data/h20_hetop_standard80.dat) or absolute.

Options:
    --max-lower-level N   Keep LinMacro lines with ll <= N  [default: 1]
    --min-osc-strength F  Keep lines with f >= F            [default: 0.0]
    --outdir DIR          Directory for converted files     [default: $SIROCCO/data/atomic_simple]

Output:
    $SIROCCO/data/<stem>_simple.dat          new master file
    $SIROCCO/data/atomic_simple/<f>_simple.dat  converted component files
"""

import os
import sys
import argparse
import re

# ---------------------------------------------------------------------------
# Path helpers
# ---------------------------------------------------------------------------

SIROCCO = os.environ.get('SIROCCO', '')


def resolve_data_path(data_path):
    """Resolve a data/-relative path to an absolute path.

    In a Sirocco run directory 'data/' is a symlink to $SIROCCO/xdata/, so
    paths in master files look like 'data/atomic_macro/h20_lines.dat' and
    resolve to $SIROCCO/xdata/atomic_macro/h20_lines.dat.
    """
    if os.path.isabs(data_path):
        return data_path
    if SIROCCO:
        # data/ -> $SIROCCO/xdata/
        xdata_path = data_path.replace('data/', 'xdata/', 1) if data_path.startswith('data/') else data_path
        candidate = os.path.join(SIROCCO, xdata_path)
        if os.path.isfile(candidate):
            return candidate
        # Also try without xdata substitution (e.g. if called with an xdata/ path)
        candidate2 = os.path.join(SIROCCO, data_path)
        if os.path.isfile(candidate2):
            return candidate2
    # Fall back to cwd-relative (handles run-directory usage where data/ symlink exists)
    local = os.path.abspath(data_path)
    if os.path.isfile(local):
        return local
    return local


def to_output_rel(abs_path):
    """Return a CWD-relative path for a converted output file.

    Since Sirocco resolves all paths in master files relative to the run
    directory, we express output paths relative to CWD.  For a master in
    data/ this yields data/atomic_simple/<name>.
    """
    return os.path.relpath(abs_path, os.getcwd())


# ---------------------------------------------------------------------------
# File-type detection
# ---------------------------------------------------------------------------

def file_keywords(filepath):
    """Return the set of macro-atom keywords present in a file."""
    found = set()
    try:
        with open(filepath) as f:
            for line in f:
                word = line.split()[0] if line.split() else ''
                if word in ('IonM', 'LevMacro', 'LinMacro', 'PhotMacS'):
                    found.add(word)
                if len(found) == 4:
                    break
    except (IOError, OSError):
        pass
    return found


# ---------------------------------------------------------------------------
# Conversion routines -- each returns (lines_written, lines_skipped)
# ---------------------------------------------------------------------------

def convert_elem_ions(inpath, outpath):
    """Convert IonM -> IonV and set nlte to 0."""
    n_converted = 0
    with open(inpath) as fin, open(outpath, 'w') as fout:
        for raw in fin:
            stripped = raw.strip()
            if not stripped.startswith('IonM'):
                fout.write(raw)
                continue

            # IonM  name  z  istate  g  ip  nmax  nlte  label
            # Replace IonM with IonV and set nlte (second-to-last numeric field) to 0.
            # We parse conservatively to preserve spacing of the label.
            parts = stripped.split()
            # parts[0]=IonM parts[1]=name parts[2]=z parts[3]=istate
            # parts[4]=g parts[5]=ip parts[6]=nmax parts[7]=nlte parts[8+]=label
            if len(parts) < 8:
                fout.write(raw)
                continue
            parts[0] = 'IonV'
            # Keep nlte (parts[7]) at its original value so that level nden
            # slots are allocated and PhotTopS cross-sections can be matched.
            # Reconstruct with consistent spacing; preserve label verbatim
            label = ' '.join(parts[8:])
            fout.write('%-8s %-4s %3s %3s %3s %12s %5s %3s     %s\n' % (
                parts[0], parts[1], parts[2], parts[3],
                parts[4], parts[5], parts[6], parts[7], label))
            n_converted += 1

    return n_converted


def convert_levels(inpath, outpath):
    """Convert LevMacro -> LevTop, inserting -1 for islp and qqnum fields."""
    n_converted = 0
    with open(inpath) as fin, open(outpath, 'w') as fout:
        for raw in fin:
            stripped = raw.strip()
            if not stripped.startswith('LevMacro'):
                fout.write(raw)
                continue

            # LevMacro  z  ion  ilv  ion_pot  ex_energy  g  rad_rate  ()  label
            # LevTop    z  ion  islp  ilv  ion_pot  ex_energy  g  qqnum  rad_rate  ()  label
            parts = stripped.split()
            # parts[0]=LevMacro [1]=z [2]=ion [3]=ilv [4]=e [5]=exx [6]=g [7]=rl [8]=() [9+]=label
            if len(parts) < 9:
                fout.write(raw)
                continue
            label = ' '.join(parts[9:])
            fout.write('LevTop %3s %3s %4d %3s %12s %12s %5s %6s %11s %s %s\n' % (
                parts[1], parts[2],
                -1,          # islp  -- not used by simple atoms
                parts[3],    # ilv
                parts[4],    # ion_pot
                parts[5],    # ex_energy
                parts[6],    # g
                -1,          # qqnum -- not used by simple atoms
                parts[7],    # rad_rate
                parts[8],    # ()
                label))
            n_converted += 1

    return n_converted


def convert_lines(inpath, outpath, max_ll, min_f):
    """Convert LinMacro -> Line, keeping only ll <= max_ll and f >= min_f."""
    n_total = 0
    n_kept = 0
    with open(inpath) as fin, open(outpath, 'w') as fout:
        for raw in fin:
            stripped = raw.strip()
            if not stripped.startswith('LinMacro'):
                fout.write(raw)
                continue

            n_total += 1
            # LinMacro  z  ion  lambda  f  gl  gu  el  eu  ll  lu
            parts = stripped.split()
            if len(parts) < 11:
                fout.write(raw)
                n_kept += 1
                continue

            ll = int(parts[9])
            f_val = float(parts[4])

            if ll > max_ll or f_val < min_f:
                continue

            # Replace keyword only; preserve all field values
            fout.write('Line' + raw[len('LinMacro'):])
            n_kept += 1

    return n_total, n_kept


def convert_phot(inpath, outpath):
    """Convert PhotMacS/PhotMac -> PhotTopS/PhotTop.

    PhotMacS format: PhotMacS z istate levl levu exx np
    PhotTopS format: PhotTopS z istate islp  ilv  exx np

    We use levl as ilv and -1 as islp, which matches the -1 islp values
    written by convert_levels for the LevTop entries.  The upper level
    index (levu) is dropped since PhotTopS only references the lower level.
    """
    n_converted = 0
    with open(inpath) as fin, open(outpath, 'w') as fout:
        for raw in fin:
            stripped = raw.strip()
            if stripped.startswith('PhotMacS'):
                # PhotMacS z istate levl levu exx np
                parts = stripped.split()
                if len(parts) >= 7:
                    z, istate, levl, levu, exx, np_ = parts[1:7]
                    fout.write(f'PhotTopS {z} {istate}  -1  {levl}  {exx}  {np_}\n')
                    n_converted += 1
                else:
                    fout.write(raw)
            elif stripped.startswith('PhotMac'):
                # PhotMac freq xsection -> PhotTop freq xsection
                fout.write('PhotTop' + raw[len('PhotMac'):])
            else:
                fout.write(raw)
    return n_converted


# ---------------------------------------------------------------------------
# Output filename helper
# ---------------------------------------------------------------------------

def simple_outpath(inpath, simple_dir):
    """Derive the output path for a converted file in simple_dir."""
    base = os.path.basename(inpath)
    root, ext = os.path.splitext(base)
    return os.path.join(simple_dir, root + '_simple' + ext)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='Convert a macro-atom Sirocco atomic data set to simple-atom format')
    parser.add_argument('master',
                        help='Master data file (e.g. data/h20_hetop_standard80.dat or absolute path)')
    parser.add_argument('--max-lower-level', type=int, default=1, metavar='N',
                        help='Keep LinMacro lines with ll <= N (default: 1, ground-state only)')
    parser.add_argument('--min-osc-strength', type=float, default=0.0, metavar='F',
                        help='Keep lines with oscillator strength >= F (default: 0.0)')
    parser.add_argument('--outdir', default=None,
                        help='Directory for converted component files '
                             '(default: $SIROCCO/data/atomic_simple)')
    args = parser.parse_args()

    # Resolve the master file
    master_abs = resolve_data_path(args.master)
    if not os.path.isfile(master_abs):
        print(f'Error: cannot find master file: {args.master}', file=sys.stderr)
        if not SIROCCO:
            print('  Hint: $SIROCCO is not set', file=sys.stderr)
        sys.exit(1)

    # Output goes alongside the input master file, without resolving symlinks.
    # os.path.abspath preserves symlinks; os.path.realpath would follow them.
    # This means that if the user passes 'data/h20_hetop_standard80.dat' where
    # data/ is a local test directory, output lands in that local directory.
    # Converted component files go in atomic_simple/ beside the master.
    # Unchanged file references keep their original data/... paths.

    master_input_dir = os.path.dirname(os.path.abspath(args.master))

    # Output directory for converted component files
    if args.outdir:
        simple_dir = os.path.abspath(args.outdir)
    else:
        simple_dir = os.path.join(master_input_dir, 'atomic_simple')

    os.makedirs(simple_dir, exist_ok=True)

    # Output master file: <stem>_simple.dat in the same directory as input
    base = os.path.basename(master_abs)
    root, ext = os.path.splitext(base)
    out_master_path = os.path.join(master_input_dir, root + '_simple' + ext)

    print(f'Input master : {master_abs}')
    print(f'Output master: {out_master_path}')
    print(f'Component dir: {simple_dir}')
    print(f'Line filter  : ll <= {args.max_lower_level}, f >= {args.min_osc_strength}')
    print()

    # Read and process the master file line by line
    new_master_lines = [
        f'# Simple-atom version of {os.path.basename(master_abs)}\n',
        f'# Generated by macro2simple.py\n',
        f'# LinMacro lines filtered: ll <= {args.max_lower_level}',
    ]
    if args.min_osc_strength > 0:
        new_master_lines[-1] += f', f >= {args.min_osc_strength}'
    new_master_lines[-1] += '\n'
    new_master_lines.append(
        '# IonM->IonV (nlte=0), LevMacro->LevTop, LinMacro->Line.\n')
    new_master_lines.append(
        '# Phot and all atomic/ files referenced in-place.\n')
    new_master_lines.append('#\n')

    with open(master_abs) as f:
        raw_lines = f.readlines()

    for raw in raw_lines:
        stripped = raw.strip()

        # Blank lines and comments pass through
        if not stripped or stripped.startswith('#'):
            new_master_lines.append(raw)
            continue

        # It's a component file reference
        data_path = stripped
        abs_path = resolve_data_path(data_path)

        if not os.path.isfile(abs_path):
            print(f'  Warning: cannot find {data_path} -- keeping reference unchanged')
            new_master_lines.append(raw)
            continue

        keywords = file_keywords(abs_path)

        if not keywords:
            # Nothing to convert -- reference in-place
            new_master_lines.append(raw)
            continue

        # File needs conversion
        out_path = simple_outpath(abs_path, simple_dir)
        out_rel = to_output_rel(out_path)

        if keywords == {'IonM'}:
            n = convert_elem_ions(abs_path, out_path)
            print(f'  IonM->IonV      : {data_path}')
            print(f'                 -> {out_rel}  ({n} ions converted)')

        elif keywords == {'LevMacro'}:
            n = convert_levels(abs_path, out_path)
            print(f'  LevMacro->LevTop: {data_path}')
            print(f'                 -> {out_rel}  ({n} levels converted)')

        elif keywords == {'LinMacro'}:
            n_total, n_kept = convert_lines(abs_path, out_path,
                                            args.max_lower_level,
                                            args.min_osc_strength)
            print(f'  LinMacro->Line  : {data_path}')
            print(f'                 -> {out_rel}  ({n_kept}/{n_total} lines kept)')

        elif keywords == {'PhotMacS'}:
            n = convert_phot(abs_path, out_path)
            print(f'  PhotMacS->PhotTopS: {data_path}')
            print(f'                 -> {out_rel}  ({n} xsection blocks converted)')

        else:
            # Mixed file -- apply conversions sequentially via a temp file
            import shutil
            tmp = out_path + '.tmp'
            tmp2 = out_path + '.tmp2'
            src = abs_path
            if 'IonM' in keywords:
                convert_elem_ions(src, tmp);  src = tmp
            if 'LevMacro' in keywords:
                convert_levels(src, tmp2);    src = tmp2
            if 'LinMacro' in keywords:
                n_total, n_kept = convert_lines(src, out_path,
                                                args.max_lower_level,
                                                args.min_osc_strength)
            elif 'PhotMacS' in keywords:
                convert_phot(src, out_path)
            else:
                shutil.copy(src, out_path)
            for t in (tmp, tmp2):
                if os.path.exists(t):
                    os.remove(t)
            print(f'  Mixed           : {data_path}  -> {out_rel}')

        new_master_lines.append(out_rel + '\n')

    # Write the new master file
    with open(out_master_path, 'w') as f:
        f.writelines(new_master_lines)

    print(f'\nWrote {out_master_path}')


if __name__ == '__main__':
    main()
