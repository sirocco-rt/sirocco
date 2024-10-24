Wind.dim.in.z_or_theta.direction
================================
Winds are calulated on spherical, cylindrical, or polar grids.
This input variable gives the size of the grid in the z or theta
direction.  Because some grid cells are used as a buffer, the
actual wind cells are contained in a slightly smaller grid than
the number given.

Note that in some situations there may be more than one wind
component, known technically as a domain.  In that case the user
will be queried for this value mulitple times, one for each domain

Type
  Integer

Values
  Greater than 0

File
  `setup_domains.c <https://github.com/sirocco-rt/sirocco/blob/master/source/setup_domains.c>`_


Parent(s)
  * :ref:`Wind.number_of_components`: Greater than 0. Once per wind.

  * :ref:`Wind.type`: Not imported


