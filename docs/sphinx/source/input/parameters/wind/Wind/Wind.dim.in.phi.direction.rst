Wind.dim.in.phi.direction
=========================
The number of grid cells in the azimuthal (phi) direction.  This parameter
is only requested when the coordinate system is ``cyl3d`` (3-D cylindrical)
or ``sph3d`` (3-D spherical polar).

For a purely axisymmetric model, set ``pdim = 1`` to use a single phi slice
covering the full 0–2π azimuth; this reproduces the corresponding 2-D
geometry with negligible overhead.  For non-axisymmetric imported models,
set ``pdim`` to the number of phi divisions in the import file.

Unlike the z/theta direction, the phi dimension is **not** doubled for
hemispheres: the grid always covers 0 to 2π in a single pass.

For analytic wind models the phi cells are always uniformly spaced.  For
imported models (``Wind.type = imported``) the phi spacing may be
non-uniform, but the phi boundaries must be **strictly increasing** in k
order and must span 0 to 2π.  The phi lookup grid is built from the ``phi``
column of the cells at (i=0, j=0, k=0 … pdim-1), so every cell with the
same k must carry the same phi boundary value.

Type
  Integer

Values
  Greater than or equal to 1.

File
  `setup_domains.c <https://github.com/sirocco-rt/sirocco/blob/master/source/setup_domains.c>`_


Parent(s)
  * :ref:`Wind.number_of_components`: Greater than or equal to 0. Once per wind.

  * :ref:`Wind.coord_system`: ``cyl3d`` or ``sph3d``
