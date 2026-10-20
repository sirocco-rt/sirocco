FITS Diagnostic Tools
#####################

Three programs are provided for exporting wind-save data to FITS binary tables,
which can be read by Python/astropy, TOPCAT, or any FITS-aware tool.
They are all invoked with a model root name and write one or more FITS extensions
to a single output file.

.. contents:: Tools
   :local:
   :depth: 1

windsave2fits
=============

**Usage**::

    windsave2fits root

Reads ``root.wind_save`` and writes ``root_cellspec.fits``.

This tool exports the per-cell spectral model that SIROCCO uses during the
ionization cycle — the band-limited radiation field that drives the ionization
balance.  Each plasma cell stores a coarse spectrum (``geo.nbins_in_cell_spec``
wavelength bins) representing the local radiation field.  The output FITS file
contains one binary-table extension per domain with one row per plasma cell and
one column per wavelength bin.

This is primarily useful for diagnosing how well the cell-spectrum approximation
captures the true radiation field, for instance by comparing the cell spectrum
against the full extracted spectrum at the same position.

**Output file:** ``root_cellspec.fits``

**Extensions:**

``cellspec``
   Binary table; one row per plasma cell, columns are the band fluxes
   (erg s\ :sup:`-1` cm\ :sup:`-2` sr\ :sup:`-1`).

windsave2table
==============

**Usage**::

    windsave2table [options] root

Reads ``root.wind_save`` and writes a set of ASCII tables, one per quantity
family.

This is the standard tool for inspecting the ionization and thermal state of
the wind.  Each output file contains a grid of values, one line per in-wind
cell, in the native coordinate system of the domain (cylindrical, rtheta, or
spherical).  Columns common to most tables are:

``x``, ``z``
   Cylindrical R and z coordinates of the inner vertex of the cell (cm).
``xcen``, ``zcen``
   Cylindrical R and z coordinates of the cell centre (cm).
``i``, ``j``
   Radial and axial grid indices.
``nwind``
   Absolute index of the cell in the wind array.
``nplasma``
   Index of the associated plasma cell.
``inwind``
   In-wind flag (0 = fully in-wind, >0 = partially or not in-wind).
``converge``
   Convergence flag (0 = converged).

Standard output files include (``root.master.txt``, ``root.heat.txt``,
``root.ion.*.txt``, etc.).  Run with ``-h`` for a full list of options.

windstruct2fits
===============

**Usage**::

    windstruct2fits root

Reads ``root.wind_save`` and writes ``root_windstruct.fits``.

This tool covers the full transport grid — all 2 × NDIM2 wind cells including
lower-hemisphere cells — not just the upper-hemisphere plasma cells exposed by
``windsave2table``.  It is designed for diagnosing hemisphere asymmetries,
inspecting velocity gradient tensors, and cross-checking cell geometry.

Because the lower-hemisphere cells share the same plasma state as their
upper-hemisphere counterparts (same ``nplasma`` index) but undergo independent
photon transport, this program is the only way to compare transport-level
diagnostics between the two hemispheres directly.

**Output file:** ``root_windstruct.fits``

**Extensions:**

``windstruct`` — wind transport cells (2 × NDIM2 rows)
   One row per wind cell, including both upper (index ``n < NDIM2``) and lower
   (index ``n ≥ NDIM2``) hemisphere cells.

   .. list-table::
      :header-rows: 1
      :widths: 18 8 55

      * - Column
        - Type
        - Description
      * - ``nwind``
        - J
        - Cell index in wmain
      * - ``ndom``
        - J
        - Domain number
      * - ``nplasma``
        - J
        - Plasma cell index (upper and lower hemisphere cells share the same value)
      * - ``inwind``
        - J
        - In-wind flag
      * - ``x``
        - E
        - Cylindrical R of inner vertex (cm)
      * - ``z``
        - E
        - Cylindrical z of inner vertex (cm)
      * - ``xcen``
        - E
        - Cylindrical R of cell centre (cm)
      * - ``zcen``
        - E
        - Cylindrical z of cell centre (cm)
      * - ``i``
        - J
        - Radial grid index
      * - ``j``
        - J
        - Axial grid index
      * - ``x_3d``
        - 3E
        - Full 3-vector inner vertex position (cm)
      * - ``xcen_3d``
        - 3E
        - Full 3-vector cell centre position (cm)
      * - ``v``
        - 3E
        - Velocity at inner vertex (cm/s)
      * - ``v_grad``
        - 9E
        - Velocity gradient tensor, row-major (1/s)
      * - ``cell_z_min``
        - E
        - Signed lower z boundary (cm)
      * - ``cell_z_max``
        - E
        - Signed upper z boundary (cm)
      * - ``dvds_ave``
        - E
        - Direction-averaged \|dv/ds\| (1/s)
      * - ``dvds_max``
        - E
        - Maximum \|dv/ds\| at vertex (1/s)
      * - ``nscat``
        - J
        - Scatter count (hemisphere-specific: ``nscat_upper`` / ``nscat_lower``)
      * - ``nrad``
        - J
        - Emission count (hemisphere-specific: ``nrad_upper`` / ``nrad_lower``)
      * - ``ntot``
        - K
        - Total photon passages through this transport cell
      * - ``ntot_star``
        - K
        - Photon passages originating from the central star
      * - ``ntot_bl``
        - K
        - Photon passages originating from the boundary layer
      * - ``ntot_disk``
        - K
        - Photon passages originating from the disk
      * - ``ntot_wind``
        - K
        - Photon passages originating from wind re-emission
      * - ``ntot_agn``
        - K
        - Photon passages originating from the AGN continuum
      * - ``nioniz``
        - K
        - H-ionizing photon passages (hν > 13.6 eV)

``plasma`` — plasma cells (NPLASMA rows)
   One row per plasma cell (upper-hemisphere in-wind cells only).  This mirrors
   the data available through ``windsave2table`` but in a single FITS binary
   table, making it convenient for Python analysis.

   .. list-table::
      :header-rows: 1
      :widths: 18 8 55

      * - Column
        - Type
        - Description
      * - ``nwind``
        - J
        - Upper-hemisphere wind cell index
      * - ``nplasma``
        - J
        - Self index
      * - ``x``
        - E
        - Cylindrical R of inner vertex (cm)
      * - ``z``
        - E
        - Cylindrical z of inner vertex (cm)
      * - ``xcen``
        - E
        - Cylindrical R of cell centre (cm)
      * - ``zcen``
        - E
        - Cylindrical z of cell centre (cm)
      * - ``i``
        - J
        - Radial grid index
      * - ``jj``
        - J
        - Axial grid index (named ``jj`` to avoid conflict with mean intensity ``j``)
      * - ``ne``
        - E
        - Electron number density (cm\ :sup:`-3`)
      * - ``rho``
        - E
        - Mass density (g cm\ :sup:`-3`)
      * - ``vol``
        - E
        - Cell volume (cm\ :sup:`3`)
      * - ``t_r``
        - E
        - Radiation temperature (K)
      * - ``t_e``
        - E
        - Electron temperature (K)
      * - ``t_r_old``
        - E
        - Radiation temperature from previous cycle (K)
      * - ``t_e_old``
        - E
        - Electron temperature from previous cycle (K)
      * - ``w``
        - E
        - Radiation dilution factor
      * - ``j``
        - E
        - Mean intensity (erg s\ :sup:`-1` cm\ :sup:`-2` sr\ :sup:`-1`)
      * - ``j_direct``
        - E
        - Mean intensity, direct photons
      * - ``j_scatt``
        - E
        - Mean intensity, scattered photons
      * - ``ave_freq``
        - E
        - Intensity-weighted mean frequency (Hz)
      * - ``ip``
        - E
        - Ionization parameter, Cloudy definition
      * - ``heat_tot``
        - E
        - Total heating rate (erg s\ :sup:`-1`)
      * - ``heat_lines`` / ``heat_ff`` / ``heat_photo`` / ``heat_auger`` / ``heat_comp``
        - E
        - Heating by process (erg s\ :sup:`-1`)
      * - ``cool_tot``
        - E
        - Total cooling rate (erg s\ :sup:`-1`)
      * - ``ntot``
        - J
        - Total photon passages (plasma estimator)
      * - ``nioniz``
        - J
        - H-ionizing photon passages (plasma estimator)
      * - ``xj``
        - 20E
        - Band-limited mean intensities
      * - ``xave_freq``
        - 20E
        - Band-limited mean frequencies (Hz)
      * - ``rad_force_es`` / ``rad_force_bf`` / ``rad_force_ff``
        - 4E
        - Radiative force components (dyne)
      * - ``F_vis`` / ``F_UV`` / ``F_Xray``
        - 4E
        - Band fluxes (erg s\ :sup:`-1` cm\ :sup:`-2`)
      * - ``lum_tot`` / ``lum_lines`` / ``lum_ff``
        - E
        - Radiative luminosities (erg s\ :sup:`-1`)
      * - ``cool_rr`` / ``cool_comp`` / ``cool_di`` / ``cool_dr``
        - E
        - Cooling rates by process (erg s\ :sup:`-1`)
      * - ``dt_e``
        - E
        - Change in electron temperature last iteration (K)
      * - ``gain``
        - E
        - Iteration gain factor
      * - ``xi``
        - E
        - Ionization parameter, TTS definition
      * - ``converge_t_r`` / ``converge_t_e`` / ``converge_hc``
        - E
        - Convergence fractional changes
      * - ``converge_whole``
        - J
        - Sum of convergence check flags
      * - ``converging``
        - J
        - Converging flag
      * - ``nscat_es`` / ``nscat_res``
        - J
        - Electron and resonant scatters
      * - ``nscat_upper`` / ``nscat_lower``
        - J
        - Scatters in upper / lower hemisphere transport cells
      * - ``nrad`` / ``nrad_upper`` / ``nrad_lower``
        - J
        - Wind photon emissions, total and per hemisphere

**Hemisphere indexing note:**  Upper-hemisphere cells have indices
``0 .. NDIM2−1``; lower-hemisphere cells have indices ``NDIM2 .. 2×NDIM2−1``.
Cells with the same ``i``, ``j`` in each hemisphere share a ``nplasma`` index
(same ionization state) but accumulate independent photon-passage counters
(``ntot``, ``ntot_star``, etc.) via ``wind_counts_main``.

**Example Python snippet:**

.. code-block:: python

    from astropy.io import fits
    import numpy as np

    with fits.open("mymodel_windstruct.fits") as hdul:
        wt = hdul["windstruct"].data
        pt = hdul["plasma"].data

    # Compare ntot between upper and lower hemisphere for cell (i=3, j=5)
    ndim2 = len(pt)           # NPLASMA == NDIM2 for a single-domain model
    mask_up = (wt["i"] == 3) & (wt["j"] == 5) & (wt["nwind"] < ndim2)
    mask_lo = (wt["i"] == 3) & (wt["j"] == 5) & (wt["nwind"] >= ndim2)
    print("upper ntot:", wt["ntot"][mask_up])
    print("lower ntot:", wt["ntot"][mask_lo])
