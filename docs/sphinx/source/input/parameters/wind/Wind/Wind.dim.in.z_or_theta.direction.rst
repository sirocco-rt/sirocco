Wind.dim.in.z_or_theta.direction
================================
The number of grid cells in the vertical or polar-angle direction
(z for cylindrical grids, theta for polar grids).  Applies to all
non-spherical coordinate systems: cylindrical, polar, cyl3d, and sph3d.
Not used for the 1-D spherical coordinate system.

For **cylindrical**, **polar**, **cyl3d**, and **sph3d** coordinate
systems this value is the number of grid cells **per hemisphere**.
The code doubles it internally so that both hemispheres (above and below
the disk plane, or upper and lower theta) are covered by independent
cells.  For example, a value of 30 produces an internal grid of 60 cells
in the z or theta direction.  A few cells at each end of the range serve
as boundary buffer (guard) cells, so the number of active wind cells is
slightly smaller than the total grid size.

Note that in some situations there may be more than one wind
component, known technically as a domain.  In that case the user
will be queried for this value multiple times, one for each domain.

Type
  Integer

Values
  Greater than 0

File
  `setup_domains.c <https://github.com/sirocco-rt/sirocco/blob/master/source/setup_domains.c>`_


Parent(s)
  * :ref:`Wind.number_of_components`: Greater than 0. Once per wind.

  * :ref:`Wind.type`: Not imported


