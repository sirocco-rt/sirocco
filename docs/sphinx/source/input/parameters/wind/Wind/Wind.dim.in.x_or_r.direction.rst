Wind.dim.in.x_or_r.direction
============================
Winds are calulated on spherical, cylindrical, or polar grids.
This input variable gives the size of the grid in the x or r
direction.  Because some grid cells are used as a buffer, the
actual wind cells are contained in a slightly smaller grid than
the number given.

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


