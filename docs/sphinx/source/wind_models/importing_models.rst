.. imported:

Importing Models
################

SIROCCO can read 1-D, 2.5-D, or full 3-D grids of density and velocity
instead of setting up the model from an analytic prescription.  Caution should
be exercised with this mode, as it is still in a development phase, and the
mode requires the user to ensure that things like mass and angular momentum
conservation are enforced.

This mode is activated via the wind type ``imported``, which triggers an extra
question for the file to import, e.g.

.. code::

   Wind.type(SV,star,hydro,corona,kwd,homologous,shell,imported)             imported
   Wind.coord_system(spherical,cylindrical,polar,cyl3d,sph3d)               cylindrical
   Wind.model2import                    cv.import.txt

An example in cylindrical geometry, :code:`cv_import.pf`, is given with a
supplementary grid file in :code:`examples/beta/`.

General rules for all import formats
=====================================

* All physical units are CGS (velocities in cm/s, lengths in cm, densities in
  g/cm\ :sup:`3`, temperatures in K) unless stated otherwise.
* Cell position coordinates are supplied at the **inner corner** (lower-left)
  of each cell, not at cell centres.
* Guard cells **must** be included at the outer edges of the grid and at the
  poles (for polar-type grids).  Guard cells should carry a plausible velocity
  but have density and temperature set to zero.
* The ``inwind`` flag marks whether a cell is in the wind:

  .. code:: c

     W_IGNORE      = -2   // ignore this grid cell (transparent to photons)
     W_NOT_INWIND  = -1   // this cell is not in the wind (treated as W_IGNORE on import)
     W_ALL_INWIND  =  0   // this cell is in the wind

  Cells with ``inwind = 1`` (partially in wind) are treated as ``W_IGNORE``.

  After the grid is built, SIROCCO checks that the required boundary slices
  contain no active (``inwind = 0``) cells.  All violations are reported to
  the error log, and then the program exits.  Every boundary that is missing
  guard cells is listed before the exit so the full extent of the problem is
  visible in one run.

* Temperature columns are optional.  If one temperature value per cell is
  provided, SIROCCO treats it as the electron temperature and sets
  :math:`T_r = 1.1\,T_e`.  If two values are provided, the first is taken as
  :math:`T_e` and the second as :math:`T_r`.  If no temperature is provided,
  both are initialised from the ``Wind.t.init`` parameter.

Spherical Grids
---------------

A 1-D spherically symmetric model.

**Required columns** (one row per cell, i increasing outward):

* ``i``              — cell index
* ``inwind``         — in-wind flag
* ``r``              — inner radial boundary (cm)
* ``v_r``            — radial velocity (cm/s)
* ``rho``            — mass density (g/cm\ :sup:`3`)
* ``t_e`` (optional) — electron temperature (K)
* ``t_r`` (optional) — radiation temperature (K)

.. admonition:: Guard cells

   Three guard cells are expected: one at the inner edge and two at the outer
   edge.  The radial coordinate must be strictly increasing.

Cylindrical Grids
-----------------

A 2.5-D axisymmetric model on a (rho, z) grid.

**Required columns** (one row per cell):

* ``i``              — radial (rho) index
* ``j``              — vertical (z) index
* ``inwind``         — in-wind flag
* ``x``              — cylindrical rho coordinate (cm)
* ``z``              — signed z coordinate (cm); negative for lower hemisphere
* ``v_x``            — x velocity (cm/s)
* ``v_y``            — y velocity (cm/s)
* ``v_z``            — z velocity (cm/s)
* ``rho``            — mass density (g/cm\ :sup:`3`)
* ``t_e`` (optional) — electron temperature (K)
* ``t_r`` (optional) — radiation temperature (K)

.. admonition:: Hemisphere coverage

   If all ``z`` values are ≥ 0 (**single-hemisphere** file), lower-hemisphere
   photons are automatically folded to the mirror upper-hemisphere cell.
   If any ``z`` value is < 0 (**two-hemisphere** file), each cell is used for
   the hemisphere it belongs to, allowing asymmetric winds.

.. admonition:: Guard cells

   Two guard cells at the outer radial boundary (large x) and two at the outer
   z boundary of each hemisphere.  For single-hemisphere files guard cells at
   z = 0 are not required.

Polar Grids
-----------

A 2.5-D model on a (r, theta) grid.

**Required columns** (one row per cell):

* ``i``              — radial index
* ``j``              — polar-angle (theta) index
* ``inwind``         — in-wind flag
* ``r``              — inner radial boundary (cm)
* ``theta``          — inner polar angle (degrees, 0–180)
* ``v_x``            — x velocity (cm/s)
* ``v_y``            — y velocity (cm/s)
* ``v_z``            — z velocity (cm/s)
* ``rho``            — mass density (g/cm\ :sup:`3`)
* ``t_e`` (optional) — electron temperature (K)
* ``t_r`` (optional) — radiation temperature (K)

.. admonition:: Hemisphere coverage

   If the maximum theta value in the file is ≤ 90°, the grid covers only the
   upper hemisphere and SIROCCO folds lower-hemisphere photons symmetrically.
   If any theta exceeds 90°, the grid covers both hemispheres explicitly.

.. admonition:: Guard cells

   One guard cell at the north pole (theta near 0°).  For two-hemisphere files,
   two guard cells at the south pole (theta near 180°) are also required.

3-D Cylindrical Grids (CYLIND3D)
---------------------------------

A full 3-D model on a (rho, z, phi) grid, with phi running from 0 to 360°.

**Required columns** (one row per cell, k varies fastest, then j, then i):

* ``i``              — radial (rho) index
* ``j``              — vertical (z) index
* ``k``              — azimuthal (phi) index
* ``inwind``         — in-wind flag
* ``x``              — cylindrical rho coordinate at lower-rho corner (cm)
* ``z``              — z coordinate at lower-z corner (cm)
* ``phi``            — azimuthal lower boundary (degrees, 0–360)
* ``v_x``            — x velocity (cm/s)
* ``v_y``            — y velocity (cm/s)
* ``v_z``            — z velocity (cm/s)
* ``rho``            — mass density (g/cm\ :sup:`3`)
* ``t_e`` (optional) — electron temperature (K)
* ``t_r`` (optional) — radiation temperature (K)

The grid must be **complete**: all ``ndim × mdim × pdim`` cells must be
present.

.. admonition:: Phi ordering

   Phi values must be **strictly increasing** in k order and must span the
   full azimuth (0–360°).  The phi lookup grid is built from the ``phi`` column
   of the cells at (i=0, j=0, k=0 … pdim-1), so all cells with the same k
   must carry the same ``phi`` boundary value.  The phi spacing may be
   non-uniform (e.g. finer sampling near a feature of interest), but the
   boundaries must be monotonically increasing.  Sirocco converts phi to
   radians internally; the import file always uses degrees.

.. admonition:: Guard cells

   Same rules as the 2-D cylindrical case applied independently in rho and z.
   No guard cells are required in the phi direction.

.. admonition:: Single phi-slice (axisymmetric)

   Setting ``pdim = 1`` in the import file (all cells have ``k = 0``) gives a
   single phi slice covering the full 360° azimuth, equivalent to the 2-D
   cylindrical case.

3-D Spherical Polar Grids (SPH3D)
-----------------------------------

A full 3-D model on a (r, theta, phi) grid, with phi running from 0 to 360°.

**Required columns** (one row per cell, k varies fastest, then j, then i):

* ``i``              — radial index
* ``j``              — polar-angle (theta) index
* ``k``              — azimuthal (phi) index
* ``inwind``         — in-wind flag
* ``r``              — inner radial boundary (cm)
* ``theta``          — inner polar-angle boundary (degrees, 0–180)
* ``phi``            — inner azimuthal boundary (degrees, 0–360)
* ``v_x``            — x velocity (cm/s)
* ``v_y``            — y velocity (cm/s)
* ``v_z``            — z velocity (cm/s)
* ``rho``            — mass density (g/cm\ :sup:`3`)
* ``t_e`` (optional) — electron temperature (K)
* ``t_r`` (optional) — radiation temperature (K)

Velocities must be given in **Cartesian** coordinates (not spherical
components), consistent with all other SIROCCO coordinate systems.

The grid must be **complete**: all ``ndim × mdim × pdim`` cells must be
present.

.. admonition:: Phi ordering

   Phi values must be **strictly increasing** in k order and must span the
   full azimuth (0–360°).  The phi lookup grid is built from the ``phi`` column
   of the cells at (i=0, j=0, k=0 … pdim-1), so all cells with the same k
   must carry the same ``phi`` boundary value.  The phi spacing may be
   non-uniform, but the boundaries must be monotonically increasing.
   Sirocco converts phi to radians internally; the import file always uses degrees.

.. admonition:: Hemisphere coverage

   As for the 2-D polar case: if the maximum theta in the file is ≤ 90°,
   only upper-hemisphere cells are present and lower-hemisphere photons are
   folded symmetrically.  If theta extends beyond 90°, both hemispheres are
   covered explicitly.

.. admonition:: Guard cells

   Same rules as the 2-D polar case applied in (r, theta).  No guard cells are
   required in the phi direction.

.. admonition:: Single phi-slice (axisymmetric)

   Setting ``pdim = 1`` (all ``k = 0``) gives a single phi slice equivalent to
   the 2-D polar case.

Example format (first two cells of a small SPH3D grid):

.. code::

   i  j  k  inwind  r       theta   phi     v_x     v_y     v_z     rho     t_e    t_r
   0  0  0    -1    1.4e9   0.0     0.0     0.0     0.0     0.0     0.0     0.0    0.0
   0  0  1    -1    1.4e9   0.0     90.0    0.0     0.0     0.0     0.0     0.0    0.0

Setting Wind Temperatures
-------------------------

Reading in a temperature is optional when importing a model.  If one
temperature value per cell is provided, SIROCCO assumes it is the electron
temperature and sets:

.. math::

   T_{r} = 1.1\,T_{e}

If two values are provided, the first is :math:`T_e` and the second is
:math:`T_r`.  If no temperature is provided, both are initialised from:

`Wind.t.init 40000`

and the electron temperature uses the Lucy approximation:

.. math::

   T_{e} = 0.9\,T_{r}

Maximum and Minimum Wind Radius
--------------------------------

The maximum and minimum spherical extent of the wind is calculated
automatically by SIROCCO from cells with ``inwind ≥ 0``.  Guard cells are
excluded from this calculation.

Python tools for import files
=============================

Two scripts in ``py_progs/`` support the import workflow:
``import_model.py`` for creating import files from SIROCCO output, and
``check_import_model.py`` for validating any import file.

Generating an import file from a SIROCCO run
---------------------------------------------

The easiest way to learn the import format for a particular coordinate system
is to run SIROCCO with an analytic wind model, then convert the output with
``import_model.py``.

.. code:: bash

   windsave2table test          # produces test.master.txt
   import_model.py test         # produces test.import.txt

``windsave2table`` writes a ``# coord_system:`` comment as the first line of
every table it produces.  ``import_model.py`` reads that comment and selects
the correct format automatically, supporting all five coordinate systems:
spherical, cylindrical, polar, cyl3d, and sph3d.

``import_model.py`` also:

* Applies a Lorentz factor correction to the density column:

  .. math::

     \rho_{\mathrm{import}} = \rho_{\mathrm{CMF}} \times \gamma

  For typical sub-relativistic winds (v/c << 1) this correction is negligible.

* Enforces guard cells at all required boundaries (outer radial rows,
  north- and south-pole rows for polar-type grids), correcting any cells that
  were not already marked ``inwind = -1``.

* Automatically runs ``check_import_model.py`` on the output file and prints a
  validation summary.

Validating an import file
--------------------------

``check_import_model.py`` can be used to validate any import file, whether
produced by ``import_model.py`` or created externally.  By default it also
fixes any guard-cell violations it finds and writes the corrected file as
``<name>.fixed.txt``:

.. code:: bash

   check_import_model.py test.import.txt        # fix + write test.import.fixed.txt
   check_import_model.py test.import.txt --no-fix   # report only

The script reports:

* Coordinate system (detected from column names)
* Grid dimensions and cell counts (total, in-wind, not-in-wind)
* Guard-cell check results for each boundary
* Physical sanity checks (velocities, density, temperature)
* A final ``CLEAN`` or problem summary

Physical issues (rho ≤ 0, T ≤ 0, v ≥ c) are reported but cannot be fixed
automatically.

.. warning::

   If you use astropy or other tools to modify an import file, check that all
   variables remain reasonable.  In particular, astropy tables can silently set
   ``rho`` to 0.0 if the output format is not specified.  The provided scripts
   set explicit formats for all floating-point columns to avoid this.

.. note::

   Spectra produced by the original analytic-model run and by the imported
   re-run will not be identical.  Two effects contribute:

   1. **Partial cells**: the original run classifies cells as partially in
      the wind and handles them specially; imported models accept only fully
      in-wind or fully out-of-wind cells.
   2. **Density placement**: SIROCCO assumes velocities at cell corners and
      densities at cell centres.  If an import file places densities at the
      same position as the corner coordinates there will be a small systematic
      offset.
