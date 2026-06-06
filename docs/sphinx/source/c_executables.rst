C Executable Programs
#####################

In addition to the main ``sirocco`` program, the distribution includes several
standalone C executables that read a wind save file (``root.wind_save``) and
perform additional analysis or export data for use with other codes.  All of
these programs are compiled and installed alongside ``sirocco`` by the standard
``make install`` build.

.. contents:: Programs
   :local:
   :depth: 1

----

swind
=====

``swind`` is an **interactive** program for inspecting the contents of a wind
save file.  It is the primary tool for examining the ionization state,
temperatures, densities, velocities, and other per-cell quantities produced by
a Sirocco run.

Usage::

    swind [-h] [-s] [-d] [-p parameter_file] [root]

Options:

``-h``
    Print a brief help message and exit.

``-s``
    Write a standard set of per-cell quantities to individual ASCII files and
    exit (non-interactive).

``-d``
    When per-cycle wind save files exist, write ASCII files for each
    ionization cycle as well as for the final state (equivalent to ``-s`` when
    no per-cycle saves are present).

``-p parameter_file``
    Read the interactive choices from *parameter_file* instead of the command
    line.  This allows the same sequence of queries to be replayed on a
    different wind save file.  The commands executed interactively are saved to
    ``swind.pf`` at the end of a session (if the user quits with ``q``).

``root``
    Root name of the wind save file.  If omitted, the program prompts for it.

Running ``swind`` without ``-s`` or ``-d`` starts an interactive session in
which the user selects quantities to display from a menu.  Output can
optionally be written to ASCII files.  The files produced use either the
original Sirocco grid (prefix ``x.``) or a re-gridded linear array (prefix
``z.``), which is convenient for creating contour plots.

----

windsave2table
==============

``windsave2table`` writes a **standard set of ASCII tables** from a wind save
file.  Unlike ``swind``, it is entirely non-interactive: the output files are
fixed and the program exits immediately after writing them.

Usage::

    windsave2table [-d] [-s] [-a] [-edge] [-x windcell] [-xall]
                   [--version] [-h] root

Options:

``-d``
    Report ion densities rather than ion fractions in the ion tables.

``-s``
    Report the number of scatters per unit volume rather than ion fractions.

``-a``
    Write additional tables with extended ion information.

``-edge``
    Include edge (boundary) cells in the output tables.

``-x windcell``
    In addition to the standard tables, write the detailed cell spectrum for
    the wind cell numbered *windcell*.

``-xall``
    Write detailed cell spectra for every in-wind cell to a single large file.

``--version``
    Print version information and exit.

``-h``
    Print this help and exit.

Output files all begin with the root name and a domain number (e.g.
``root.0.master.txt``).  Multiple domains produce separate files.  Every file
starts with columns ``x``, ``z``, ``i``, ``j``, ``inwind`` so that plotting
routines can identify each cell unambiguously.  The *master* file records
:math:`n_e`, velocity, :math:`\rho`, :math:`T_e`, :math:`T_r`, and ion
fractions for a set of key ions.

----

sirocco_optd
============

``sirocco_optd`` (formerly ``py_optical_depth``) computes **optical depth
spectra and photospheric surfaces** through the Sirocco wind model.

Usage::

    sirocco_optd [-h] [-d ndom] [-p tau_stop] [-cion nion]
                 [-freq_min min] [-freq_max max] [-i i1 i2 ...]
                 [--nonrel] [--smax frac] [--no-es] [--version]
                 root

Options:

``-h``
    Print the full help message and exit.

``-d ndom``
    Set the domain from which photons are launched.

``-p tau_stop``
    Instead of tracing photons to escape, integrate outwards from the
    inner boundary to find the surface of constant electron-scattering
    optical depth *tau_stop* (i.e. the photosphere).

``-cion nion``
    Extract the column density for the ion indexed *nion*.

``-freq_min min``, ``-freq_max max``
    Frequency boundaries for the optical-depth spectrum.

``-i i1 i2 ...``
    Calculate optical depths along the listed sight-line inclinations
    (degrees from the pole).  Overrides the inclinations defined in the
    model.

``--nonrel``
    Use linear frequency transforms.  Use this when Sirocco was run in
    non-relativistic mode.

``--smax frac``
    Set the maximum fraction of a cell width that a photon may travel in
    a single step.

``--no-es``
    Exclude electron scattering from the opacity.

``--version``
    Print version information and exit.

By default the program integrates the continuum optical depth along every
observer line of sight defined in the model.  If no observers were defined,
a set of default sight lines is used.

----

windsave2fits
=============

``windsave2fits`` exports selected wind-save quantities to **FITS files**.  It
requires the optional ``cfitsio`` library and is the most convenient way to
export large two-dimensional data arrays (e.g. cell spectra) for analysis in
Python or other tools that read FITS.

.. note::
   ``windsave2fits`` is built only when ``cfitsio`` is detected during
   ``./configure``.  If the executable is absent, ``cfitsio`` is likely not
   installed.

Usage::

    windsave2fits [options] root

The primary use case is exporting the per-cell spectral model data (the
spectral bands used to estimate ionization rates) in a form that can be
compared with the actual cell spectra recorded during transport.

----

rad_hydro_files
===============

``rad_hydro_files`` is a **post-processing tool for radiation-hydrodynamics
(rad-hydro) coupling**.  It reads a Sirocco wind save file and writes a
collection of ASCII files that communicate the radiative heating, cooling,
driving forces, and ionization state of the wind to a hydrodynamics code.

The outputs are used by the `PLUTO–Sirocco
<https://github.com/sirocco-rt/pluto-sirocco>`_ coupled rad-hydro framework,
in which PLUTO advances the hydrodynamics and Sirocco provides the
radiative-transfer physics at each exchange step.

Usage::

    rad_hydro_files root

where ``root`` is the root name of a Sirocco wind save file.  The parameter
file is **not** read; all configuration is taken from the wind save.  Because
several of the derived quantities (line cooling, recombination emissivities)
are not stored directly in the wind save, ``rad_hydro_files`` recomputes them
by calling the standard Sirocco cooling and luminosity routines before writing
the output.

Output files
~~~~~~~~~~~~

All files are written to the current directory:

``py_heatcool.dat``
    Per-cell heating and cooling rates: position, :math:`T_e`, ionisation
    parameter :math:`\xi`, :math:`n_e`, X-ray photoionisation heating,
    Compton heating, line heating, free–free heating, Compton cooling, line
    and recombination cooling, free–free cooling, :math:`\rho`, and hydrogen
    number density.

``py_driving.dat``
    Per-cell radiation driving forces: geometry, density, electron density,
    directional flux in three bands (optical, UV, X-ray), and the electron
    scattering and bound-free radiation force vectors.

``py_ion_data.dat``
    Per-cell ion density for every ion in the atomic data set.

``py_spec_data.dat``
    Per-cell spectral model data (band boundaries and model parameters for
    the spectral estimator used in the ionization calculation).

``py_pcon_data.dat``
    Per-cell Sobolev parameter data: :math:`T_e`, :math:`\rho`, hydrogen
    number density, :math:`n_e`, and dimensionless optical-depth parameters
    for the optical, UV, and X-ray bands.

``py_fluxes.dat``
    Per-cell directional flux in each of the three photometric bands.

``directional_flux_theta.dat``, ``directional_flux_phi.dat``, ``directional_flux_r.dat``
    Per-cell UV flux decomposed into angular bins in the :math:`\theta`,
    :math:`\phi`, and :math:`r` directions respectively, intended for
    computing anisotropic radiation forces.

``py_debug_data.dat``
    Diagnostic file: position, thermal velocity, velocity gradient, and
    mean intensity per cell.

MPI
~~~

``rad_hydro_files`` is MPI-enabled and may be run in parallel to accelerate
the cooling-rate recomputation::

    mpirun -n 8 rad_hydro_files root

----

modify_wind
===========

``modify_wind`` is a **developer/diagnostic tool** for modifying the contents
of a wind save file without rerunning Sirocco.  Its primary current use is
overwriting ion densities in selected cells with prescribed values, making it
useful for constructing controlled test cases.

.. note::
   ``modify_wind`` is a prototype.  The modifications to be applied must be
   hard-coded in the source file ``source/modify_wind.c`` before recompiling.
   It is not intended for routine use.

Usage::

    modify_wind root

The program reads the wind save ``root.wind_save``, applies the coded
modifications, and writes the result to a new wind save file.

----

inspect_wind
============

``inspect_wind`` is a **diagnostic program** that prints selected internal
variables from a wind save file to ASCII files.  It was originally written to
allow detailed inspection of macro-atom variables in parallel-mode runs (where
direct memory inspection is not straightforward) and is intended to be
customised in source for each diagnostic task.

.. note::
   Like ``modify_wind``, ``inspect_wind`` requires the user to edit
   ``source/inspect_wind.c`` to select which variables to output, then
   recompile.

Usage::

    inspect_wind root
