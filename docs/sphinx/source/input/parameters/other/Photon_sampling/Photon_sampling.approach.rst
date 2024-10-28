Photon_sampling.approach
========================
Choice of whether and how to use stratified sampling in creating photons during the
ionization stage of the calculation.

Type
  Enumerator

Values
  T_star
    Sets a single band based on the temperature given

  cv
    Traditional cv setup

  yso
    YSO setup

  AGN
    Test for balance matching the bands we have been using for AGN runs

  min_max_freq
    Mode 1 sets a single wide band defined by f1 and f2

  user_bands
    User-defined bands

  cloudy_test
    Set up to compare with cloudy power law table command note
    that this also sets up the weight and photon index for the PL, to ensure a continuous distribution

  wide
    Test for balance to have a really wide frequency range

  logarithmic
    Generalized method to set up logarithmic bands


File
  `bands.c <https://github.com/sirocco-rt/sirocco/blob/master/source/bands.c>`_


Child(ren)
  * :ref:`Photon_sampling.nbands`

  * :ref:`Photon_sampling.high_energy_limit`

  * :ref:`Photon_sampling.low_energy_limit`

