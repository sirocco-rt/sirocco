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
* Ghost cells **must** be included at the outer edges of the grid and at the
  poles (for polar-type grids).  Ghost cells should carry a plausible velocity
  but have density and temperature set to zero.
* The ``inwind`` flag marks whether a cell is in the wind:

  .. code:: c

     W_IGNORE      = -2   // ignore this grid cell
     W_NOT_INWIND  = -1   // this cell is not in the wind
     W_ALL_INWIND  =  0   // this cell is in the wind

  Cells with ``inwind = 1`` (partially in wind) are treated as ``W_IGNORE``.

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

A full 3-D model on a (rho, z, phi) grid, with phi running from 0 to 2 pi.

**Required columns** (one row per cell, k varies fastest, then j, then i):

* ``i``              — radial (rho) index
* ``j``              — vertical (z) index
* ``k``              — azimuthal (phi) index
* ``inwind``         — in-wind flag
* ``x``              — cylindrical rho coordinate at lower-rho corner (cm)
* ``z``              — z coordinate at lower-z corner (cm)
* ``phi``            — azimuthal lower boundary (radians, 0–2 pi)
* ``v_x``            — x velocity (cm/s)
* ``v_y``            — y velocity (cm/s)
* ``v_z``            — z velocity (cm/s)
* ``rho``            — mass density (g/cm\ :sup:`3`)
* ``t_e`` (optional) — electron temperature (K)
* ``t_r`` (optional) — radiation temperature (K)

The grid must be **complete**: all ``ndim × mdim × pdim`` cells must be
present.  The phi grid must be **uniformly spaced** from 0 to 2 pi; the last
phi boundary is set automatically to 2 pi and need not be supplied.

.. admonition:: Guard cells

   Same rules as the 2-D cylindrical case applied independently in rho and z.
   No guard cells are required in the phi direction.

.. admonition:: Single phi-slice (axisymmetric)

   Setting ``pdim = 1`` in the import file (all cells have ``k = 0``) gives a
   single phi slice covering the full 0–2 pi azimuth, equivalent to the 2-D
   cylindrical case.

3-D Spherical Polar Grids (SPH3D)
-----------------------------------

A full 3-D model on a (r, theta, phi) grid, with phi running from 0 to 2 pi.

**Required columns** (one row per cell, k varies fastest, then j, then i):

* ``i``              — radial index
* ``j``              — polar-angle (theta) index
* ``k``              — azimuthal (phi) index
* ``inwind``         — in-wind flag
* ``r``              — inner radial boundary (cm)
* ``theta``          — inner polar-angle boundary (degrees, 0–180)
* ``phi``            — inner azimuthal boundary (radians, 0–2 pi)
* ``v_x``            — x velocity (cm/s)
* ``v_y``            — y velocity (cm/s)
* ``v_z``            — z velocity (cm/s)
* ``rho``            — mass density (g/cm\ :sup:`3`)
* ``t_e`` (optional) — electron temperature (K)
* ``t_r`` (optional) — radiation temperature (K)

Velocities must be given in **Cartesian** coordinates (not spherical
components), consistent with all other SIROCCO coordinate systems.

The grid must be **complete**: all ``ndim × mdim × pdim`` cells must be
present.  The phi grid must be **uniformly spaced** from 0 to 2 pi.

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
   0  0  0    -1    1.4e9   0.0     0.0000  0.0     0.0     0.0     0.0     0.0    0.0
   0  0  1    -1    1.4e9   0.0     1.5708  0.0     0.0     0.0     0.0     0.0    0.0

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

Generating import files from a SIROCCO run
==========================================

The easiest way to learn the import format for a particular coordinate system
is to first run SIROCCO with an analytic wind model in the same coordinate
system, then convert the output to an import file using the provided Python
scripts.

For example, to create an import file for a 3-D spherical polar model with
root name ``test``:

.. code:: bash

   windsave2table test          # produces test.master.txt
   import_sph3d.py test         # produces test.import.txt

The import file can then be used in a new run with ``Wind.type = imported``.

The available conversion scripts in ``py_progs/`` are:

* ``import_1d.py``      — 1-D spherical models
* ``import_cyl.py``     — 2.5-D cylindrical models
* ``import_rtheta.py``  — 2.5-D polar (r-theta) models
* ``import_cyl3d.py``   — 3-D cylindrical models (CYLIND3D)
* ``import_sph3d.py``   — 3-D spherical polar models (SPH3D)

Each script reads ``<root>.master.txt`` and writes ``<root>.import.txt``.  The
scripts also apply a Lorentz factor correction to the density column:

.. math::

   \rho_{\mathrm{import}} = \rho_{\mathrm{CMF}} \times \gamma

For typical sub-relativistic winds (v/c << 1) this correction is negligible.

.. warning::

   If you use astropy or other tools to modify an import file, check that all
   variables remain reasonable.  In particular, we have observed that astropy
   tables can silently set ``rho`` to 0.0 if the format is not specified.
   The provided scripts set explicit output formats for all floating-point
   columns to avoid this.

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
