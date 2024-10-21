Boundary_layer.rad_type_to_make_wind
====================================
Determines the luminosity and SED of the boundary layer.
The code can cause a source to radiate differently in the ionisation and spectral cycles.
This variable allows a boundary layer source to radiate differently from :ref:`Boundary_layer.rad_type_in_final_spectrum`
during the cycles used to calculate the wind ionisation state and temperature.

Type
  Enumerator

Values
  bb
    Black-body radiation. The boundary layer radiates as a black-body source with surface luminosity set by its
    effective temperature (:ref:`Boundary_layer.temp`) and resulting in a total luminosity
    proportional to its surface area.

  models
    Radiate according to a model. SIROCCO can support tabulated models that output with a binned luminosity distribution
    depending on system properties like temperature and gravity. See :ref:`Input_spectra.model_file`. The total
    luminosity will be set by :ref:`Boundary_layer.luminosity`.

  power
    Radiate following a power-law model as :math:`L_\nu=K\nu^\alpha`. The total luminosity will be set by :ref:`Boundary_layer.luminosity`.


File
  `setup_star_bh.c <https://github.com/agnwinds/python/blob/master/source/setup_star_bh.c>`_


Parent(s)
  * :ref:`Boundary_layer.radiation`: ``True``


Child(ren)
  * :ref:`Boundary_layer.power_law_index`

  * :ref:`Input_spectra.model_file`

  * :ref:`Boundary_layer.luminosity`

  * :ref:`Boundary_layer.power_law_cutoff`

  * :ref:`Boundary_layer.temp`

