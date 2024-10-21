Wind.ionization
===============
The approach used by SIROCCO to calculate the ionization
of the wind during ionization cycles.  A number of these
modes are historical or included for diagnostic purposes.

Type
  Enumerator

Values
  on.the.spot
    Use a simple on the spot approximation to calculated the ionization.

  LTE_te
    Calculate ionization based on the Saha equation using
    the electron temperature.  (This is intended as a diagnostic
    mode.)

  LTE_tr
    Calculate ionization based on the Saha equation using
    the radiation temperature. (This is intended as a diagnstic mode)

  ML93
    Use the modified on the spot approimation  described by 
    `Mazzli & Lucy 1993 <https://ui.adsabs.harvard.edu/abs/1993A%26A...279..447M/abstract>`_  

  fixed
    Read the ion aboundances in from a file.  All cells will have
    the same abundances. (This is intended
    as a diagnostic mode, mainly to investigate the details of raditive transrfer.
    It should be used with caution.  In particular, if the elements for which
    abundances are provided differ from the elements to be used as described in the
    elements/ions portion of the atomic data, then one should not expect the calculated
    electron density to be that that comes simply from the fixd concetnrations file.)

  matrix_bb
    Estimate photoionization rates by approximating the spectrum in
    each cell based on the radiation temperature and an effective
    weight.  Invert the rate matrix equations to calculate the ionization

  matrix_pow
    Estimate photionization rates by approximating the spectrum in a cell by a piecewise
    approximation, usually a power law.  Invert the rate matrix equation to
    calculate the ionization. (This is the preferred ionization mode for most
    calculations)

  matrix_est
    Estimate photoionization rates by calculating rates directly from the photons that pass
    through a cell.  There is no attempt to model the spectrum. Invert the rate matrix equation to
    calculate the ionization.


File
  `setup.c <https://github.com/agnwinds/python/blob/master/source/setup.c>`_


Child(ren)
  * :ref:`Wind.fixed_concentrations_file`

