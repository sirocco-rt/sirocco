#!/usr/bin/env python
"""
                    Space Telescope Science Institute

Synopsis:

Calculate LTE ion abundances using Sirocco atomic data and the Saha equation.


Command line usage:

    usage: calc_lte_abundances.py [-h] [-v] masterfile temperature nh

    Examples:
        calc_lte_abundances.py data/standard80.dat 10000 1e10
        calc_lte_abundances.py data/standard80.dat 50000 1e8 -v

Description:

    This script calculates ionization equilibrium using the Saha equation,
    following the approach in Sirocco's saha.c. It reads atomic data from
    a masterfile and computes ion fractions for a given temperature and
    hydrogen number density.

    Before running, ensure that Setup_Sirocco_Dir has been run to create
    the necessary symlinks so that file paths in the masterfile are valid.

Primary routines:

    steer           Processes command line options and calls do_one
    do_one          Loads data, runs calculation, prints and returns results
    calc_ionization Performs the ionization calculation, returns results as
                    a list of dictionaries suitable for conversion to a Table
    read_atomicdata Reads atomic data from a masterfile

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

        with open(masterfile, 'r') as f:
            for line in f:
                line = line.strip()
                # Skip comments and empty lines
                if not line or line.startswith('#'):
                    continue

                if os.path.exists(line):
                    data_files.append(line)
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
                    self._parse_ion(parts)
                elif keyword == 'LevTop':
                    self._parse_level_topbase(parts)
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

    def _parse_ion(self, parts: List[str]):
        """Parse IonV line."""
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

        ion = Ion(name=name, z=z, istate=istate, g=g, ip=ip, ip_ev=ip_ev)
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
    """Calculate partition functions for all ions."""
    kt = BOLTZMANN * temperature

    for i, ion in enumerate(data.ions):
        z_partition = ion.g
        ground_ex = 0.0

        for level in data.levels:
            if level.z == ion.z and level.istate == ion.istate:
                if level.ex < ground_ex or ground_ex == 0.0:
                    ground_ex = level.ex

        for level in data.levels:
            if level.z == ion.z and level.istate == ion.istate:
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


def calc_ionization(data: AtomicData, temperature: float, nh: float,
                    verbose: bool = False) -> List[Dict]:
    """
    Calculate ionization equilibrium for given conditions.

    This is the main calculation routine that can be called from scripts.

    Args:
        data: AtomicData object from read_atomicdata()
        temperature: Temperature in Kelvin
        nh: Hydrogen number density in cm^-3
        verbose: Print convergence information

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

    if verbose:
        print(f"Initial ne estimate: {ne:.3e} cm^-3")

    # Iterate to convergence
    for iteration in range(MAXITERATIONS):
        densities = _saha_equation(data, ne, t, nh)
        ne_new = _get_electron_density(data, densities)
        ne_new = max(ne_new, DENSITY_MIN)

        rel_error = abs(ne - ne_new) / ne_new if ne_new > 0 else 1.0

        if verbose and iteration % 10 == 0:
            print(f"  Iteration {iteration}: ne = {ne_new:.3e}, rel_error = {rel_error:.3e}")

        if rel_error < FRACTIONAL_ERROR or ne_new < 1e-6:
            if verbose:
                print(f"Converged after {iteration + 1} iterations")
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

        if row['fraction'] > 1e-10:
            print(f"{row['ion']:<10} {row['istate']:<8} {row['density']:<18.3e} {row['fraction']:<12.4e}")


def do_one(masterfile, temperature, nh, verbose=False):
    """
    Calculate ionization equilibrium for given conditions.

    This is the main routine to call from scripts.

    Args:
        masterfile: Path to masterfile (e.g., 'data/standard80.dat')
        temperature: Temperature in Kelvin
        nh: Hydrogen number density in cm^-3
        verbose: Print convergence and loading information

    Returns:
        List of dictionaries with keys:
            element, z, ion, istate, density, fraction, ne

        Can be converted to astropy Table:
            from astropy.table import Table
            t = Table(results)
    """
    # Load atomic data
    if verbose:
        print(f"Loading atomic data from {masterfile}")

    data = read_atomicdata(masterfile)

    if verbose:
        print(f"Loaded {len(data.elements)} elements, {len(data.ions)} ions, {len(data.levels)} levels")

    # Run calculation
    results = calc_ionization(data, temperature, nh, verbose=verbose)

    # Print results
    if verbose:
        print_results(results, temperature, nh)

    return results


def steer(argv):
    """
    Parse command line arguments and call do_one.

    Args:
        argv: Command line arguments (sys.argv)
    """
    verbose = False

    if len(argv) < 4:
        print(__doc__)
        return

    # Parse arguments
    positional = []
    i = 1
    while i < len(argv):
        if argv[i] == '-h' or argv[i] == '--help':
            print(__doc__)
            return
        elif argv[i] == '-v' or argv[i] == '--verbose':
            verbose = True
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

    results=do_one(masterfile, temperature, nh, verbose=verbose)


# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    if len(sys.argv) > 1:
        steer(sys.argv)
    else:
        print(__doc__)
