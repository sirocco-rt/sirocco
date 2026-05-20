Coordinate grids
----------------

SIROCCO supports five coordinate gridding schemes:

* 1-D spherical
* 2-D cylindrical
* 2-D polar (r-theta)
* 3-D cylindrical (rho, z, phi) — ``cyl3d``
* 3-D spherical polar (r, theta, phi) — ``sph3d``

These options are controlled by the
:doc:`/input/parameters/wind/Wind/Wind.coord_system` keyword.

Although SIROCCO incorporates several analytic wind models (e.g. SV, KWD,
Knigge), the velocities and other properties are placed on a discrete grid
during program setup.  It is up to the user to choose an appropriate coordinate
system and the number of grid points for any particular run.

For the 2-D and 3-D coordinate systems, cells are spaced logarithmically in
the radial direction so that cells grow larger with distance from the central
source.  For polar (r-theta) coordinates the angular spacing is linear; for
cylindrical coordinates the z-spacing is also logarithmic.  For imported models
the exact coordinate gridding is set by the user.

Two-hemisphere grids
====================

For the 2-D cylindrical, 2-D polar, 3-D cylindrical, and 3-D spherical polar
coordinate systems, SIROCCO uses a **full two-hemisphere grid** in which the
upper hemisphere (z > 0 or theta < 90 deg) and the lower hemisphere (z < 0 or
theta > 90 deg) each have their own independent set of cells.  This allows the
wind to be asymmetric above and below the disk plane, and ensures that
estimators accumulated in one hemisphere do not contaminate the other.

The user specifies the number of grid cells in the z (or theta) direction
**per hemisphere** via
:doc:`/input/parameters/wind/Wind/Wind.dim.in.z_or_theta.direction`.
The code doubles this internally so that the full grid covers both hemispheres.
For example, specifying 30 theta-cells gives an internal grid of 60 theta-cells:
indices 0–29 cover theta 0–90 deg (upper hemisphere) and indices 30–59 cover
theta 90–180 deg (lower hemisphere).

For the 3-D coordinate systems the same hemisphere doubling is applied in the
theta (or z) direction; the additional phi dimension is specified separately via
:doc:`/input/parameters/wind/Wind/Wind.dim.in.phi.direction` and is **not**
doubled.

Guard cells
-----------

A few cells at the outer boundary of each hemisphere are reserved as guard
(boundary buffer) cells and are unconditionally set to ``W_NOT_INWIND`` without
consulting the wind-cone geometry.  Their purpose is to give the grid a
well-defined outer edge and to prevent the interpolation stencil from reaching
beyond the last real cell.

For **cylindrical** grids, ``cylind_is_cell_in_wind`` excludes any cell for
which the within-hemisphere z-index satisfies ``j_hemi >= mdim_half - 2``,
i.e. the two outermost cells of both the upper and lower hemispheres are guard
cells.  In the radial direction, the two outermost r-cells (``i >= ndim - 2``)
are also guard cells.

For **polar** grids, ``rtheta_is_cell_in_wind`` excludes:

* The two outermost radial cells (``i >= ndim - 2``), as in the cylindrical case.
* In the upper hemisphere (j = 0 .. mdim/2-1): the cell at j = 0 is a
  north-pole guard (theta near 0).  There is **no equatorial guard** — the cell
  immediately above the disk plane (j = mdim/2 - 1, theta just below 90 deg) is
  classified by the wind-cone geometry, not pre-excluded.
* In the lower hemisphere (j = mdim/2 .. mdim-1): the two cells closest to the
  south pole (``j_hemi >= mdim_half - 2``, theta near 180 deg) are guard cells.

For **3-D cylindrical** (CYLIND3D) and **3-D spherical polar** (SPH3D) grids,
the same guard-cell rules apply independently in (rho or r) and (z or theta).
The phi dimension covers 0 to 2 pi without guard cells; azimuthal extent is
determined purely by the ``inwind`` flags in the import file.

.. note::
   The ``cyl_var`` coordinate system (cylindrical with radially-varying z
   boundaries) is not supported in this version of SIROCCO.

3-D Cylindrical coordinates (CYLIND3D)
=======================================

The 3-D cylindrical coordinate system extends the 2-D cylindrical grid by
adding an azimuthal dimension.  Cells are indexed by (i, j, k) where:

* ``i`` is the radial (rho) index, increasing outward
* ``j`` is the vertical (z) index, with the lower hemisphere occupying
  ``j = 0 .. mdim/2-1`` and the upper hemisphere ``j = mdim/2 .. mdim-1``
* ``k`` is the azimuthal (phi) index, running from 0 to ``pdim-1`` and
  covering 0 to 2 pi uniformly

The number of phi cells is set by
:doc:`/input/parameters/wind/Wind/Wind.dim.in.phi.direction`.
For a purely axisymmetric model, ``pdim = 1`` gives a single phi slice
covering the full 0–2 pi azimuth; this is equivalent to the 2-D cylindrical
geometry with a small computational overhead.

3-D Spherical polar coordinates (SPH3D)
=========================================

The 3-D spherical polar coordinate system extends the 2-D polar (r-theta) grid
by adding an azimuthal (phi) dimension.  Cells are indexed by (i, j, k) where:

* ``i`` is the radial index, increasing outward
* ``j`` is the polar-angle (theta) index, with the upper hemisphere (theta 0–90
  deg) occupying ``j = 0 .. mdim/2-1`` and the lower hemisphere (theta 90–180
  deg) occupying ``j = mdim/2 .. mdim-1``
* ``k`` is the azimuthal (phi) index, running from 0 to ``pdim-1`` and
  covering 0 to 2 pi uniformly

The number of phi cells is set by
:doc:`/input/parameters/wind/Wind/Wind.dim.in.phi.direction`.
As for CYLIND3D, ``pdim = 1`` gives a single phi slice covering the full
azimuth and is equivalent to the 2-D polar geometry.

Cell volume for SPH3D is

.. math::

   V = \frac{\Delta\phi}{3}(r_{\max}^3 - r_{\min}^3)(\cos\theta_{\min} - \cos\theta_{\max})

Velocity interpolation
======================

**2-D grids (spherical, cylindrical, polar)**

Within each 2-D grid cell SIROCCO performs a bilinear interpolation of the
wind velocity using the four cell corners.  For cylindrical cells the corners
are at :math:`(\rho_{\min}, z_{\min})`, :math:`(\rho_{\max}, z_{\min})`,
:math:`(\rho_{\min}, z_{\max})`, and :math:`(\rho_{\max}, z_{\max})`.  For
polar cells the corners are at :math:`(r_{\min}, \theta_{\min})`,
:math:`(r_{\max}, \theta_{\min})`, :math:`(r_{\min}, \theta_{\max})`, and
:math:`(r_{\max}, \theta_{\max})`.  For 1-D spherical cells a linear
interpolation in r between the inner and outer boundaries is used.

The four corner velocities are stored in the per-cell fields ``v``,
``v_rmax``, ``v_thetamax``, and ``vmax`` of the ``WindPtr`` array, computed
during grid initialisation.

**3-D grids (CYLIND3D, SPH3D)**

For 3-D cells SIROCCO performs a **trilinear** interpolation using all eight
corners of the cell.  The phi-min face uses the same four corners as the 2-D
case (``v``, ``v_rmax``, ``v_thetamax``, ``vmax``); the phi-max face uses
four additional fields ``v_phimax``, ``v_rmax_phimax``,
``v_thetamax_phimax``, and ``vmax_phimax``.

For axisymmetric analytic models all eight corners have the same (r, theta)
velocity, so the trilinear interpolation reduces to bilinear in (r, theta) or
(rho, z).  For imported non-axisymmetric models the phi-face corner velocities
carry genuine phi variation.

Grid size and memory
====================

The fidelity of a model improves with the number of grid cells, but larger
grids consume more memory and more compute time.  The most important quantity
is the number of cells that are **in the wind**, since Sirocco maintains full
plasma and radiation-field information only for those cells.  A narrow wind on
a 100×100 grid may require far less memory than a wide-angle wind on the same
grid.

For 3-D grids the total cell count is ``ndim × mdim × pdim`` (before the
hemisphere doubling of mdim), and memory requirements grow accordingly.

Partial cells
=============

Parameterized wind models often have regions of space that are in the wind and
regions that are not.  When a coordinate grid is overlaid on such a model some
cells cross the wind boundary — these are *partial cells*.

SIROCCO treats partial cells as follows: by default they are excluded from the
calculation and their densities are set to zero.  Because densities are
interpolated from cell corners, this can affect the first fully in-wind cell
adjacent to the boundary.

Two alternatives are available:

* Partial cells can be included.  This is most reasonable for the KWD model,
  where the velocity law extends naturally outside the wind.  Partial cells
  tend to converge more slowly than fully in-wind cells because of their
  smaller effective volume.

* As an advanced option, partial cells can be excluded but assigned the density
  of the nearest fully in-wind cell, so that photon transport near the boundary
  does not interpolate toward zero.

In most cases the treatment of partial cells has only a minor effect on the
predicted spectrum, but it is worth checking in situations where a significant
fraction of photons traverse very few wind cells (e.g. a narrow wind with a
small opening angle).
