Spectrum_cycles
===============

In SIROCCO, the detailed spectra are calculate with flights of photons that have
the same number of photon bundles as in ionization cycles. The detailed spectra
are calculated between a minimum and maximum wavelenth, and at specific inclination
angles. Various other conditions can be placed on each spectrum that is accumulated.

The spectra are actually output at the end of each spectrum cycle with a normalization
that is correct for a distance of 100 pc, so one can look at the detailed spectra as
it builds up, if one wishes.

Type
  Integer

File
  `setup.c <https://github.com/agnwinds/python/blob/master/source/setup.c>`_


Child(ren)
  * :ref:`Spectrum.orbit_phase`

  * :ref:`Spectrum.no_observers`

  * :ref:`Spectrum.wavemin`

  * :ref:`Spectrum.select_photons_by_position`

  * :ref:`Spectrum.type`

  * :ref:`Spectrum.live_or_die`

  * :ref:`Spectrum.select_specific_no_of_scatters_in_spectra`

  * :ref:`Spectrum.wavemax`

