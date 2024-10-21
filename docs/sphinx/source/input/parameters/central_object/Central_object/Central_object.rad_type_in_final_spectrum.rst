Central_object.rad_type_in_final_spectrum
=========================================

Determines the spectral energy distribution of the central object in the spectral cycles. The luminosity is set by the options for the
ionisation cycles, however.

Type
  Enumerator

Values
  bb
    Available for :ref:`System_type` of ``star`` or ``cv``.
    Black-body radiation. The boundary layer radiates as a black-body source with surface luminosity set by its
    effective temperature (:ref:`Central_object.temp`) and resulting in a total luminosity
    proportional to its surface area.

  models
    Available for :ref:`System_type` of ``star`` or ``cv``.
    Radiate according to a model. SIROCCO can support tabulated models that output with a binned luminosity distribution
    depending on system properties like temperature and gravity. See :ref:`Input_spectra.model_file`. The total
    luminosity will be set by :ref:`Central_object.luminosity`.

  uniform
    Available for :ref:`System_type` of ``star`` or ``cv``.
    Photons are generated with a random, uniformly-distributed wavelength between
    :ref:`Spectrum.wavemin` and :ref:`Spectrum.wavemax`. Can in some cases substitute for a Kurcuz spectrum.
    This mode is only available when generating final spectra.

  brems
    Available for :ref:`System_type` of ``agn`` or ``bh``.
    Central object radiates with spectral energy distribution of a brehmsstralung spectrum as :math:`L_\nu=\nu^{\alpha}e^{-h\nu/kT}`.
    This was originally developed to allow comparison to spectra generated
    according to Blondin heating and cooling rates.

  cloudy
    Available for :ref:`System_type` of ``agn`` or ``bh``.
    Central object radiates with a 'broken' power law, intended largely for testing purposes against Cloudy.
    The spectral energy distribution form is :math:`L_\nu=K\nu^\alpha``. However, beyond the provided high and low energy
    breakpoints, the luminosity falls off sharply.

  power
    Available for :ref:`System_type` of ``agn`` or ``bh``.
    Radiate following a power-law model as :math:`L_\nu=K\nu^\alpha`.
    The total luminosity will be set by :ref:`Boundary_layer.luminosity`.


File
  `python.c <https://github.com/agnwinds/python/blob/master/source/python.c>`_


Parent(s)
  * :ref:`Central_object.radiation`: ``True``


Child(ren)
  * :ref:`Input_spectra.model_file`

