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

Type
  Integer

Values
  Greater than or equal to 1.

File
  `setup_domains.c <https://github.com/sirocco-rt/sirocco/blob/master/source/setup_domains.c>`_


Parent(s)
  * :ref:`Wind.number_of_components`: Greater than or equal to 0. Once per wind.

  * :ref:`Wind.coord_system`: ``cyl3d`` or ``sph3d``
