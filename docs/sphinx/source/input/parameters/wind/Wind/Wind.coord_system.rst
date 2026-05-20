Wind.coord_system
=================
The coordinate system used for a describing a component of the wind.

Type
  Enumerator

Values
  spherical
    1-D spherical (r only)

  cylindrical
    2-D cylindrical (rho, z)

  polar
    2-D spherical polar (r, theta)

  cyl3d
    3-D cylindrical (rho, z, phi); requires
    :doc:`Wind.dim.in.phi.direction <Wind.dim.in.phi.direction>`

  sph3d
    3-D spherical polar (r, theta, phi); requires
    :doc:`Wind.dim.in.phi.direction <Wind.dim.in.phi.direction>`


.. note::
   The ``cyl_var`` option (cylindrical with radially-varying z boundaries)
   is not supported in this version of SIROCCO.


File
  `setup_domains.c <https://github.com/sirocco-rt/sirocco/blob/master/source/setup_domains.c>`_


Parent(s)
  * :ref:`Wind.number_of_components`: Greater than 0. Once per wind.


