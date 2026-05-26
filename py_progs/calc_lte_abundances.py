#!/usr/bin/env python
"""
                    Space Telescope Science Institute

Synopsis:

Calculate LTE ion abundances using Sirocco atomic data and the Saha equation.


Command line usage::

    calc_lte_abundances.py masterfile temperature nh
    calc_lte_abundances.py data/standard80.dat 10000 1e10

    calc_lte_abundances.py -many [options]
    calc_lte_abundances.py -many -masterfile data/standard80.dat -nh 1e10
    calc_lte_abundances.py -many -log_tmin 4 -log_tmax 5 -n_per_decade 20
    calc_lte_abundances.py -many -ion_tables

Options for ``-many`` mode:

- ``-masterfile FILE``  Atomic data masterfile (default: data/h20_hetop_standard80.dat)
- ``-nh VALUE``         Hydrogen number density in cm^-3 (default: 1e9)
- ``-log_tmin VALUE``   Log10 of minimum temperature (default: 3)
- ``-log_tmax VALUE``   Log10 of maximum temperature (default: 6)
- ``-n_per_decade N``   Number of temperature points per decade (default: 30)
- ``-outfile FILE``     Output filename (default: auto-generated from masterfile)
- ``-ion_tables``       Also produce per-element ion tables (like windsave2table format)

Description:

This script calculates ionization equilibrium using the Saha equation,
following the approach in Sirocco's saha.c. It reads atomic data from
a masterfile and computes ion fractions for a given temperature and
hydrogen number density.

Before running, ensure that Setup_Sirocco_Dir has been run to create
the necessary symlinks so that file paths in the masterfile are valid.

The -ion_tables option produces additional output files for each element,
with ion fractions as columns (i01, i02, etc.) similar to the format
produced by windsave2table.

Primary routines:

- ``steer``            Processes command line options and calls do_one or do_many
- ``do_one``           Loads data, runs calculation for a single temperature
- ``do_many``          Runs calculations over a range of temperatures, writes table
- ``get_ion``          Extracts per-element tables from do_many output
- ``calc_ionization``  Performs the ionization calculation
- ``read_atomicdata``  Reads atomic data from a masterfile

Notes:

    Results are returned as a list of dictionaries with keys:
        element, z, ion, istate, density, fraction, ne

    This can be converted to an astropy Table:
        from astropy.table import Table
        t = Table(results)

    To select specific ions:
        t[t['ion'] == 'H II']
        t[(t['z'] == 6) & (t['istate'] == 3)]  # C III

History:

2024    Generated for Sirocco project

"""

import os
import sys
import math
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

import numpy as np
from astropy.table import Table, join
from astropy.io import ascii

# Physical constants (matching Sirocco's constants.h)
BOLTZMANN = 1.38062e-16      # erg/K
EV2ERGS = 1.602192e-12       # eV to ergs conversion
SAHA_CONST = 4.82907e15      # 2*(2*pi*MELEC*k)^1.5 / h^3

# Numerical parameters (matching Sirocco's saha.c)
MIN_TEMP = 100.0             # Minimum temperature (K)
MAXITERATIONS = 200          # Max iterations for ne convergence
FRACTIONAL_ERROR = 0.03      # Convergence criterion for ne
THETAMAX = 1e4               # Ionization parameter limit
DENSITY_MIN = 1e-20          # Minimum density floor

# Roman numerals for ion notation
ROMAN = ['I', 'II', 'III', 'IV', 'V', 'VI', 'VII', 'VIII', 'IX', 'X',
         'XI', 'XII', 'XIII', 'XIV', 'XV', 'XVI', 'XVII', 'XVIII',
         'XIX', 'XX', 'XXI', 'XXII', 'XXIII', 'XXIV', 'XXV', 'XXVI', 'XXVII']


@dataclass
class Element:
    """Represents an element with its properties."""
    name: str
    z: int                   # Atomic number
    abun: float              # Abundance (linear, relative to H = 1)
    atomic_weight: float
    firstion: int = -1       # Index to first ion in ion list
    nions: int = 0           # Number of ions for this element


@dataclass
class Ion:
    """Represents an ion with its properties."""
    name: str
    z: int                   # Atomic number
    istate: int              # Ionization state (1=neutral, 2=singly ionized, etc.)
    g: float                 # Ground state statistical weight
    ip: float                # Ionization potential in ergs
    ip_ev: float             # Ionization potential in eV (for reference)
    is_macro: bool = False   # True if this is a macro atom (IonM)
    firstlevel: int = -1     # Index to first level
    nlevels: int = 0         # Number of levels
    partition: float = 1.0   # Partition function


@dataclass
class Level:
    """Represents an energy level."""
    z: int                   # Atomic number
    istate: int              # Ionization state
    g: float                 # Statistical weight
    ex: float                # Excitation energy in ergs
    ex_ev: float             # Excitation energy in eV (for reference)
    is_macro: bool = False   # True if this is a macro atom level (LevMacro)


class AtomicData:
    """Container for all atomic data."""

    def __init__(self):
        self.elements: Dict[int, Element] = {}  # keyed by atomic number z
        self.ions: List[Ion] = []
        self.levels: List[Level] = []
        self.ion_index: Dict[Tuple[int, int], int] = {}  # (z, istate) -> ion index

    def parse_masterfile(self, masterfile: str) -> List[str]:
        """
        Parse masterfile to get list of data files.

        Args:
            masterfile: Path to masterfile

        Returns:
            List of data file paths
        """
        data_files = []
        missing_files = []

        # Get directory of masterfile to resolve relative paths
        masterfile_dir = os.path.dirname(os.path.abspath(masterfile))

        with open(masterfile, 'r') as f:
            for line in f:
                line = line.strip()
                # Skip comments and empty lines
                if not line or line.startswith('#'):
                    continue

                # Try the path as-is first, then relative to masterfile directory
                if os.path.exists(line):
                    data_files.append(line)
                else:
                    alt_path = os.path.join(masterfile_dir, line)
                    if os.path.exists(alt_path):
                        data_files.append(alt_path)
                    else:
                        missing_files.append(line)

        if missing_files:
            print(f"Warning: Could not find the following files:")
            for f in missing_files:
                print(f"  {f}")

        return data_files

    def _read_data_file(self, filepath: str):
        """Read a single atomic data file."""
        with open(filepath, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue

                parts = line.split()
                if len(parts) < 2:
                    continue

                keyword = parts[0]

                if keyword == 'Element':
                    self._parse_element(parts)
                elif keyword == 'IonV':
                    self._parse_ion(parts, is_macro=False)
                elif keyword == 'IonM':
                    self._parse_ion(parts, is_macro=True)
                elif keyword == 'LevTop':
                    self._parse_level_topbase(parts)
                elif keyword == 'LevMacro':
                    self._parse_level_macro(parts)
                elif keyword == 'Level':
                    self._parse_level(parts)

    def _parse_element(self, parts: List[str]):
        """Parse Element line."""
        if len(parts) < 5:
            return

        z = int(parts[1])
        name = parts[2]
        log_abun = float(parts[3])
        atomic_weight = float(parts[4])
        abun = 10.0 ** (log_abun - 12.0)

        self.elements[z] = Element(
            name=name, z=z, abun=abun, atomic_weight=atomic_weight
        )

    def _parse_ion(self, parts: List[str], is_macro: bool = False):
        """Parse IonV or IonM line.

        Args:
            parts: Split line from data file
            is_macro: True if this is a macro atom (IonM), False for simple atom (IonV)
        """
        if len(parts) < 6:
            return

        name = parts[1]
        z = int(parts[2])
        istate = int(parts[3])
        g = float(parts[4])
        ip_ev = float(parts[5])
        ip = ip_ev * EV2ERGS

        if ip_ev > 1e15:
            ip = 0.0

        ion = Ion(name=name, z=z, istate=istate, g=g, ip=ip, ip_ev=ip_ev, is_macro=is_macro)
        self.ions.append(ion)

    def _parse_level_topbase(self, parts: List[str]):
        """Parse LevTop line (Topbase format)."""
        if len(parts) < 8:
            return

        z = int(parts[1])
        istate = int(parts[2])
        ex_ev = float(parts[6])
        g = float(parts[7])
        ex = ex_ev * EV2ERGS

        level = Level(z=z, istate=istate, g=g, ex=ex, ex_ev=ex_ev)
        self.levels.append(level)

    def _parse_level_macro(self, parts: List[str]):
        """Parse LevMacro line (macro atom format).

        Format: LevMacro z istate lvl energy ex_energy g rate description
        Example: LevMacro 1 1 1 -13.59843 0.000000 2 1.00e+21 () n=1
        """
        if len(parts) < 7:
            return

        z = int(parts[1])
        istate = int(parts[2])
        ex_ev = float(parts[5])  # excitation energy in eV
        g = float(parts[6])      # statistical weight
        ex = ex_ev * EV2ERGS

        level = Level(z=z, istate=istate, g=g, ex=ex, ex_ev=ex_ev, is_macro=True)
        self.levels.append(level)

    def _parse_level(self, parts: List[str]):
        """Parse Level line."""
        if len(parts) < 5:
            return

        z = int(parts[1])
        istate = int(parts[2])
        g = float(parts[4])
        ex_ev = float(parts[5]) if len(parts) > 5 else 0.0
        ex = ex_ev * EV2ERGS

        level = Level(z=z, istate=istate, g=g, ex=ex, ex_ev=ex_ev)
        self.levels.append(level)

    def _build_indices(self):
        """Build indices for fast lookup."""
        # Sort ions by (z, istate) to ensure they are contiguous by element
        # and in order of ionization state
        self.ions.sort(key=lambda ion: (ion.z, ion.istate))

        for i, ion in enumerate(self.ions):
            self.ion_index[(ion.z, ion.istate)] = i

            if ion.z in self.elements:
                elem = self.elements[ion.z]
                if elem.firstion < 0:
                    elem.firstion = i
                elem.nions += 1

        for level in self.levels:
            key = (level.z, level.istate)
            if key in self.ion_index:
                ion_idx = self.ion_index[key]
                ion = self.ions[ion_idx]
                if ion.firstlevel < 0:
                    ion.firstlevel = self.levels.index(level)
                ion.nlevels += 1


def read_atomicdata(masterfile: str) -> AtomicData:
    """
    Read atomic data from a masterfile.

    Args:
        masterfile: Path to masterfile (e.g., 'data/standard80.dat')

    Returns:
        AtomicData object with loaded data
    """
    data = AtomicData()
    data_files = data.parse_masterfile(masterfile)

    for filepath in data_files:
        data._read_data_file(filepath)

    data._build_indices()
    return data


def _calculate_partition_functions(data: AtomicData, temperature: float):
    """Calculate partition functions for all ions.

    For macro atoms (IonM), only LevMacro levels are used.
    For simple atoms (IonV), only LevTop/Level levels are used.
    """
    kt = BOLTZMANN * temperature

    for i, ion in enumerate(data.ions):
        z_partition = ion.g
        ground_ex = 0.0

        # Filter levels: macro ions use macro levels, simple ions use simple levels
        for level in data.levels:
            if level.z == ion.z and level.istate == ion.istate:
                # Only use levels that match the ion's macro/simple status
                if level.is_macro != ion.is_macro:
                    continue
                if level.ex < ground_ex or ground_ex == 0.0:
                    ground_ex = level.ex

        for level in data.levels:
            if level.z == ion.z and level.istate == ion.istate:
                # Only use levels that match the ion's macro/simple status
                if level.is_macro != ion.is_macro:
                    continue
                if level.ex > ground_ex:
                    delta_e = level.ex - ground_ex
                    if delta_e < 10 * kt:
                        z_partition += level.g * math.exp(-delta_e / kt)

        ion.partition = max(z_partition, ion.g)


def _saha_equation(data: AtomicData, ne: float, temperature: float, nh: float) -> List[float]:
    """Calculate ion densities using the Saha equation."""
    xsaha = SAHA_CONST * pow(temperature, 1.5)
    kt = BOLTZMANN * temperature

    densities = [0.0] * len(data.ions)

    for z, elem in data.elements.items():
        if elem.firstion < 0:
            continue

        first = elem.firstion
        last = first + elem.nions

        densities[first] = 1.0
        total = 1.0
        big = pow(10.0, 250.0 / max(elem.nions, 1))

        for nion in range(first + 1, last):
            prev_ion = data.ions[nion - 1]
            curr_ion = data.ions[nion]

            if prev_ion.partition > 0 and ne > 0:
                b = (xsaha * curr_ion.partition *
                     math.exp(-prev_ion.ip / kt) /
                     (ne * prev_ion.partition))
                b = min(b, big)
            else:
                b = big

            densities[nion] = densities[nion - 1] * b
            total += densities[nion]

        if total > 0:
            norm = nh * elem.abun / total
            for nion in range(first, last):
                densities[nion] *= norm
                densities[nion] = max(densities[nion], DENSITY_MIN)

    return densities


def _get_electron_density(data: AtomicData, densities: List[float]) -> float:
    """Calculate electron density from ion densities."""
    ne = 0.0
    for i, ion in enumerate(data.ions):
        charge = ion.istate - 1
        ne += densities[i] * charge
    return max(ne, DENSITY_MIN)


def calc_ionization(data: AtomicData, temperature: float, nh: float) -> List[Dict]:
    """
    Calculate ionization equilibrium for given conditions.

    This is the main calculation routine that can be called from scripts.

    Args:
        data: AtomicData object from read_atomicdata()
        temperature: Temperature in Kelvin
        nh: Hydrogen number density in cm^-3

    Returns:
        List of dictionaries with keys:
            element: Element name (e.g., 'H', 'He', 'C')
            z: Atomic number
            ion: Ion name (e.g., 'H I', 'H II', 'C IV')
            istate: Ionization state (1=neutral, 2=singly ionized, etc.)
            density: Number density in cm^-3
            fraction: Ion fraction for this element
            ne: Electron density in cm^-3

        Can be converted to astropy Table:
            from astropy.table import Table
            t = Table(results)
    """
    t = max(temperature, MIN_TEMP)

    # Calculate partition functions
    _calculate_partition_functions(data, t)

    # Initial electron density estimate from H ionization
    xsaha = SAHA_CONST * pow(t, 1.5)

    h1_idx = data.ion_index.get((1, 1))
    if h1_idx is not None:
        h1_ip = data.ions[h1_idx].ip
        theta = xsaha * math.exp(-h1_ip / (BOLTZMANN * t)) / nh

        if theta < THETAMAX:
            x = (-theta + math.sqrt(theta * theta + 4 * theta)) / 2.0
            ne = x * nh
        else:
            ne = nh
    else:
        ne = 0.5 * nh

    ne = max(ne, 1e-6)

    # Iterate to convergence
    for iteration in range(MAXITERATIONS):
        densities = _saha_equation(data, ne, t, nh)
        ne_new = _get_electron_density(data, densities)
        ne_new = max(ne_new, DENSITY_MIN)

        rel_error = abs(ne - ne_new) / ne_new if ne_new > 0 else 1.0

        if rel_error < FRACTIONAL_ERROR or ne_new < 1e-6:
            break

        ne = (ne + ne_new) / 2.0
    else:
        print(f"Warning: Failed to converge after {MAXITERATIONS} iterations")

    # Build results list
    results = []

    for z in sorted(data.elements.keys()):
        elem = data.elements[z]
        if elem.firstion < 0:
            continue

        # Calculate total density for this element
        total = 0.0
        for nion in range(elem.firstion, elem.firstion + elem.nions):
            total += densities[nion]

        for nion in range(elem.firstion, elem.firstion + elem.nions):
            ion = data.ions[nion]
            density = densities[nion]
            fraction = density / total if total > 0 else 0.0

            # Ion notation
            if ion.istate <= len(ROMAN):
                ion_name = f"{elem.name} {ROMAN[ion.istate - 1]}"
            else:
                ion_name = f"{elem.name}^{ion.istate - 1}+"

            results.append({
                'element': elem.name,
                'z': z,
                'ion': ion_name,
                'istate': ion.istate,
                'density': density,
                'fraction': fraction,
                'ne': ne_new
            })

    return results


def print_results(results: List[Dict], temperature: float, nh: float):
    """
    Print formatted results to stdout.

    Args:
        results: List of dictionaries from calc_ionization()
        temperature: Temperature in K
        nh: Hydrogen number density in cm^-3
    """
    if not results:
        print("No results to display")
        return

    ne = results[0]['ne']

    print("\n" + "=" * 70)
    print("IONIZATION EQUILIBRIUM RESULTS")
    print("=" * 70)
    print(f"Temperature:        {temperature:.2e} K")
    print(f"H number density:   {nh:.2e} cm^-3")
    print(f"Electron density:   {ne:.2e} cm^-3")
    print(f"ne/nh ratio:        {ne/nh:.3f}")
    print("=" * 70)

    current_element = None
    for row in results:
        if row['element'] != current_element:
            current_element = row['element']
            # Find abundance from first row of this element
            print(f"\n{current_element} (Z={row['z']})")
            print("-" * 50)
            print(f"{'Ion':<10} {'State':<8} {'Density (cm^-3)':<18} {'Fraction':<12}")
            print("-" * 50)

        # Always print hydrogen, otherwise filter by fraction
        if row['z'] == 1 or row['fraction'] > 1e-10:
            print(f"{row['ion']:<10} {row['istate']:<8} {row['density']:<18.3e} {row['fraction']:<12.4e}")


def do_one(masterfile, temperature, nh):
    """
    Calculate ionization equilibrium for given conditions.

    This is the main routine to call from scripts.

    Args:
        masterfile: Path to masterfile (e.g., 'data/standard80.dat')
        temperature: Temperature in Kelvin
        nh: Hydrogen number density in cm^-3

    Returns:
        List of dictionaries with keys:
            element, z, ion, istate, density, fraction, ne

        Can be converted to astropy Table:
            from astropy.table import Table
            t = Table(results)
    """
    data = read_atomicdata(masterfile)
    results = calc_ionization(data, temperature, nh)
    return results


def do_many(masterfile='data/h20_hetop_standard80.dat', nh=1e9, log_tmin=3, log_tmax=6, n_per_decade=30, outfile=''):
    """
    Calculate LTE ionization over a range of temperatures and write to a table.

    Args:
        masterfile: Path to atomic data masterfile
        nh: Hydrogen number density in cm^-3
        log_tmin: Log10 of minimum temperature in K
        log_tmax: Log10 of maximum temperature in K
        n_per_decade: Number of temperature points per decade
        outfile: Output filename (auto-generated if empty)

    Returns:
        Astropy Table with columns: element, z, ion, istate, density,
        fraction, ne, t_e, nh

    The output file uses ascii.fixed_width_two_line format.
    """
    delta=1./n_per_decade
    log_t=log_tmin
    results=[]
    while log_t<=log_tmax:
        # print(log_t)
        result=do_one(masterfile,10**log_t,nh)
        for r in result:
            r['t_e']     = 10.0**log_t
            r['nh']    = nh
        results.extend(result)
        log_t+=delta
    all_results=Table(rows=results)
    if outfile=='':
        word=masterfile.split('/')[-1]
        outfile=word.split('.')[0]
        outfile='%s_nh_%.2e.txt' % (outfile,nh)
    all_results.write(outfile,format='ascii.fixed_width_two_line',overwrite=True)
    return all_results


def get_ion(element, filename='lte_from_data.txt', xsave=False):
    """
    Extract per-element ion table from do_many output, similar to windsave2table format.

    Args:
        element: Element name (e.g., 'H', 'He', 'C')
        filename: Input file from do_many
        xsave: If True, write the result to a file

    Returns:
        Astropy Table with columns: t_e, ne, i01, i02, ... (ion fractions)
        Returns empty list if element not found.

    The output format resembles windsave2table ion files, with temperature
    and electron density followed by ion fraction columns.
    """
    xtab = ascii.read(filename)
    xelem = xtab[xtab['element'] == element]
    if len(xelem) == 0:
        print('Could not find %s' % element)
        print(np.unique(xtab['element']))
        return []

    istate = np.unique(xelem['istate'])
    start = True
    for one_state in istate:
        xx = xelem[xelem['istate'] == one_state]
        ztab = xx['t_e', 'ne', 'fraction']
        ztab.rename_column('fraction', 'i%02d' % one_state)
        if start:
            result = ztab.copy()
            start = False
        else:
            result = join(result, ztab['t_e', 'i%02d' % one_state], join_type='left')

    if xsave:
        # Derive output filename from input filename
        outfile = filename.replace('.txt', '.%s.frac.txt' % element)
        if outfile == filename:
            outfile = '%s.%s.frac.txt' % (filename, element)
        result.write(outfile, format='ascii.fixed_width_two_line', overwrite=True)
        print(f"  Wrote {element} table to {outfile}")

    return result

import matplotlib.pyplot as plt
from astropy.io import ascii

def plot_one(
    sirr='sirocco_bb/Sum_sirroco_He_bb.txt',
    cloud='cloudy_bb/Sum_Cloudy_helium.txt',
    outfile=''
):
    """
    This is just a plotting routine, intended to make it possible
    to compate two different sets of data of the same element.

    The figure has two panels:
      - Top: log-log comparison of fractional abundances
      - Bottom: difference (sirocco - cloudy) on a linear scale,
        with cloudy interpolated onto the sirocco temperature grid
    """
    cloudy  = ascii.read(cloud)
    sirocco = ascii.read(sirr)

    # Find ion abundance columns (pattern: i01, i02, etc.)
    import re
    ion_pattern = re.compile(r'^i(\d+)$')
    cols = [c for c in sirocco.colnames if ion_pattern.match(c)]
    # Sort by ion number to ensure correct order
    cols.sort(key=lambda c: int(ion_pattern.match(c).group(1)))

    fig, (ax, ax_diff) = plt.subplots(2, 1, figsize=(8, 10),
                                       gridspec_kw={'height_ratios': [2, 1]})

    # Get temperature arrays
    t_sirocco = np.array(sirocco['t_e'])
    t_cloudy = np.array(cloudy['t_e'])

    # Use one color per ion stage, reused for both files
    colors = {}
    for col in cols:
        ion = int(col.replace('i', ''))
        roman_label = ROMAN[ion - 1] if ion <= len(ROMAN) else f'{ion}+'

        # Sirocco: points
        ax.loglog(
            sirocco['t_e'], sirocco[col],
            marker='o', linestyle='none',
            label=roman_label
        )
        colors[col] = ax.lines[-1].get_color()

        # Cloudy: lines, same color
        ax.loglog(
            cloudy['t_e'], cloudy[col],
            linestyle=':', alpha=0.6,
            color=colors[col]
        )

        # Interpolate cloudy onto sirocco temperature grid for difference
        # Use log of temperature for interpolation since grid is logarithmic
        cloudy_interp = np.interp(
            np.log10(t_sirocco),
            np.log10(t_cloudy),
            np.array(cloudy[col])
        )
        diff = np.array(sirocco[col]) - cloudy_interp

        # Plot difference in bottom panel
        ax_diff.semilogx(
            t_sirocco, diff,
            marker='o', linestyle='-', markersize=4,
            color=colors[col],
            label=roman_label
        )

    ax.set_xlim(1e3, 1e6)
    ax.set_ylim(1e-4, 2)

    ax.set_xlabel(r'T$_e$', fontsize=16)
    ax.set_ylabel('Fractional Abundance', fontsize=16)
    ax.set_title('Sirocco vs Cloudy', fontsize=14)
    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.tick_params(axis='both', which='minor', labelsize=14)

    first=sirr.split('/')[-1]
    second=cloud.split('/')[-1]
    first_element=first.split('.')[-3]
    second_element=second.split('.')[-3]
    first=first.replace('.txt','')
    second=second.replace('.txt','')

    if first_element==second_element:
        ax.set_title(first_element, fontsize=14)

    # Determine legend layout based on number of ions
    n_ions = len(cols)
    if n_ions <= 6:
        legend_ncol = 2
        legend_fontsize = 9
    elif n_ions <= 12:
        legend_ncol = 4
        legend_fontsize = 8
    else:
        # Many ions (e.g., Fe with 27): use more columns
        legend_ncol = 9
        legend_fontsize = 7

    # Secondary legend explaining styles (add first so ion legend is on top)
    from matplotlib.lines import Line2D
    style_legend = [
        Line2D([0], [0], marker='o', linestyle='none', label=first),
        Line2D([0], [0], linestyle=':', label=second)
    ]
    style_leg = ax.legend(
        handles=style_legend,
        loc='lower left',
        frameon=False
    )
    ax.add_artist(style_leg)

    # Main legend: ion stages
    if n_ions > 12:
        # Place legend below the bottom panel
        ax.legend(ncol=legend_ncol, fontsize=legend_fontsize, loc='upper right')
    else:
        ax.legend(ncol=legend_ncol, fontsize=legend_fontsize, loc='best')

    # Configure difference panel
    ax_diff.set_xlim(1e3, 1e6)
    ax_diff.axhline(y=0, color='gray', linestyle='--', linewidth=0.8)
    ax_diff.set_xlabel(r'T$_e$', fontsize=16)
    ax_diff.set_ylabel('Difference', fontsize=14)
    ax_diff.tick_params(axis='both', which='major', labelsize=14)
    ax_diff.tick_params(axis='both', which='minor', labelsize=12)

    # For the difference panel, place legend below when many ions
    if n_ions > 12:
        ax_diff.legend(ncol=legend_ncol, fontsize=legend_fontsize,
                       loc='upper center', bbox_to_anchor=(0.5, -0.15))
    else:
        ax_diff.legend(ncol=legend_ncol, fontsize=legend_fontsize, loc='best')

    if n_ions > 12:
        # Make room for legend below the plot
        fig.tight_layout(rect=[0, 0.08, 1, 1])
    else:
        fig.tight_layout()
    if outfile!='':
        plt.savefig(outfile)

def steer(argv):
    """
    Parse command line arguments and call do_one or do_many.

    Args:
        argv: Command line arguments (sys.argv)
    """
    # Check for -many mode first
    if len(argv) > 1 and argv[1] == '-many':
        steer_many(argv[2:])
        return

    if len(argv) < 4:
        print(__doc__)
        return

    # Parse arguments for single calculation
    positional = []
    i = 1
    while i < len(argv):
        if argv[i] == '-h' or argv[i] == '--help':
            print(__doc__)
            return
        elif argv[i][0] == '-':
            print(f'Error: Unknown switch --- {argv[i]}')
            return
        else:
            positional.append(argv[i])
        i += 1

    # Expect: masterfile temperature nh
    if len(positional) >= 3:
        masterfile = positional[0]
        temperature = float(positional[1])
        nh = float(positional[2])
    else:
        print("Error: Need masterfile, temperature, and nh")
        print(__doc__)
        return

    results = do_one(masterfile, temperature, nh)
    print_results(results, temperature, nh)


def steer_many(argv):
    """
    Parse command line arguments for do_many.

    Args:
        argv: Command line arguments after '-many'
    """
    # Default values matching do_many signature
    masterfile = 'data/h20_hetop_standard80.dat'
    nh = 1e9
    log_tmin = 3
    log_tmax = 6
    n_per_decade = 30
    outfile = ''
    ion_tables = False

    i = 0
    while i < len(argv):
        if argv[i] == '-h' or argv[i] == '--help':
            print("Usage: calc_lte_abundances.py -many [options]")
            print("\nOptions:")
            print("  -masterfile FILE   Atomic data masterfile (default: data/h20_hetop_standard80.dat)")
            print("  -nh VALUE          Hydrogen number density in cm^-3 (default: 1e9)")
            print("  -log_tmin VALUE    Log10 of minimum temperature (default: 3)")
            print("  -log_tmax VALUE    Log10 of maximum temperature (default: 6)")
            print("  -n_per_decade N    Number of temperature points per decade (default: 30)")
            print("  -outfile FILE      Output filename (default: auto-generated)")
            print("  -ion_tables        Also produce per-element ion tables (like windsave2table)")
            return
        elif argv[i] == '-masterfile' and i + 1 < len(argv):
            masterfile = argv[i + 1]
            i += 2
        elif argv[i] == '-nh' and i + 1 < len(argv):
            nh = float(argv[i + 1])
            i += 2
        elif argv[i] == '-log_tmin' and i + 1 < len(argv):
            log_tmin = float(argv[i + 1])
            i += 2
        elif argv[i] == '-log_tmax' and i + 1 < len(argv):
            log_tmax = float(argv[i + 1])
            i += 2
        elif argv[i] == '-n_per_decade' and i + 1 < len(argv):
            n_per_decade = int(argv[i + 1])
            i += 2
        elif argv[i] == '-outfile' and i + 1 < len(argv):
            outfile = argv[i + 1]
            i += 2
        elif argv[i] == '-ion_tables':
            ion_tables = True
            i += 1
        else:
            print(f'Error: Unknown option or missing value --- {argv[i]}')
            return

    print(f"Running do_many with:")
    print(f"  masterfile:   {masterfile}")
    print(f"  nh:           {nh:.2e}")
    print(f"  log_tmin:     {log_tmin}")
    print(f"  log_tmax:     {log_tmax}")
    print(f"  n_per_decade: {n_per_decade}")
    if outfile:
        print(f"  outfile:      {outfile}")
    if ion_tables:
        print(f"  ion_tables:   True")

    results = do_many(masterfile, nh, log_tmin, log_tmax, n_per_decade, outfile)

    # Determine actual output filename for ion tables
    if outfile == '':
        word = masterfile.split('/')[-1]
        actual_outfile = word.split('.')[0]
        actual_outfile = '%s_nh_%.2e.txt' % (actual_outfile, nh)
    else:
        actual_outfile = outfile

    print(f"\nWrote {len(results)} rows to {actual_outfile}")

    if ion_tables:
        print("\nGenerating per-element ion tables:")
        elements = np.unique(results['element'])
        for elem in elements:
            get_ion(elem, filename=actual_outfile, xsave=True)


# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    if len(sys.argv) > 1:
        steer(sys.argv)
    else:
        print(__doc__)
