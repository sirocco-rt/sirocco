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
=============
**3-D grids (CYLIND3D, SPH3D)**

For 3-D cells SIROCCO performs a **trilinear** interpolation in (rho, z, phi)
for CYLIND3D, or (r, theta, phi) for SPH3D.  The bilinear weights in the first two
coordinates are combined with linear weights between the two neighbouring phi
slices, so eight values are used.  The interpolation is periodic in phi: the last
slice wraps round to the first.  The per-cell corner fields used by the 2-D grids
are not used for 3-D grids.

In a 3-D grid the velocity stored for each cell, ``v``, is the velocity at the
cell's inner vertex at the cell's own azimuth, expressed in **local**
(v_rho, v_phi, v_z) components at that azimuth.  The velocity gradient
``v_grad`` is stored in the same local frame.  Because the components are local,
values from neighbouring phi slices can be interpolated directly; the result is
then rotated to the Cartesian frame at the photon's position.  Imported 3-D
models use the same convention (see :doc:`wind_models/importing_models`).

For axisymmetric models every phi slice holds the same local velocities, so the
interpolation reduces to the bilinear 2-D case.

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
