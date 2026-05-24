Wind.dim.in.x_or_r.direction
============================
The number of grid cells in the radial direction (rho for cylindrical
grids, r for spherical or polar grids).  Applies to all coordinate
systems: spherical, cylindrical, polar, cyl3d, and sph3d.

Because some cells at the outer boundary are reserved as guard (buffer)
cells, the number of active wind cells is slightly smaller than the value
given here.

Note that in some situations there may be more than one wind
component, known technically as a domain.  In that case the user
will be queried for this value multiple times, one for each domain.

Type
  Integer

Values
  Greater than or equal to 4, to allow for boundaries.

File
  `setup_domains.c <https://github.com/sirocco-rt/sirocco/blob/master/source/setup_domains.c>`_


Parent(s)
  * :ref:`Wind.number_of_components`: Greater than or equal to 0. Once per wind.

  * :ref:`Wind.type`: Not imported


