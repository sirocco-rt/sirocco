Boundary_layer.power_law_cutoff
===============================
This is a low frequency cutoff for an AGN-style power law spectrum
of a form :math:`L_\nu=K\nu^\alpha`, as applied to the boundary layer of a star.
It prevents the power-law being applied to low frequencies and giving an odd SED.

Type
  Double

Unit
  Hz

Values
  Greater than 0

File
  `setup_star_bh.c <https://github.com/sirocco-rt/sirocco/blob/master/source/setup_star_bh.c>`_


Parent(s)
  * :ref:`Boundary_layer.rad_type_to_make_wind`: ``power_law``
