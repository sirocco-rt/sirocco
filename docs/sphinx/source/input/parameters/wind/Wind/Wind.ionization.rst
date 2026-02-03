Wind.ionization
===============
The approach used by SIROCCO to calculate the ionization
of the wind during ionization cycles.  A number of these
modes are historical or included for diagnostic purposes.

Type
  Enumerator

Values
  on.the.spot
    Use a simple on-the-spot approximation to calculate the ionization.

  LTE_te
    Calculate ionization assuming Local Thermodynamic Equilibrium (LTE)
    using the Saha equation with the electron temperature t_e. The electron
    temperature is set equal to the initial radiation temperature (twind),
    or from the imported model if available, and remains fixed throughout
    the simulation. The partition functions are evaluated at t_e. This mode
    forces Saha ionization for all species, including macro atoms.
    (This is intended as a diagnostic mode.)

  LTE_tr
    Calculate ionization assuming Local Thermodynamic Equilibrium (LTE)
    using the Saha equation with the radiation temperature t_r. The
    radiation temperature is updated each ionization cycle based on the
    mean frequency of the radiation field. The electron temperature is
    set to t_e = 0.9 * t_r each cycle. The partition functions are
    evaluated at t_r. This mode forces Saha ionization for all species,
    including macro atoms. (This is intended as a diagnostic mode.)

  ML93
    Use the modified on-the-spot approximation described by
    `Mazzali & Lucy 1993 <https://ui.adsabs.harvard.edu/abs/1993A%26A...279..447M/abstract>`_.
    At each ionization cycle, the electron temperature is set to
    t_e = 0.9 * t_r. Ionization fractions are first computed using
    the Saha equation at t_r, then corrected using the Lucy-Mazzali
    estimators which account for the dilute, non-Planckian radiation field.
    The partition functions are evaluated using the dilution factor W.
    This mode does not attempt to balance heating and cooling to determine
    t_e self-consistently.  

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

  matrix_multi
    Similar to matrix_pow, except in this case a more aggessive approach to reaching
    ion state/temperature balance is attempted.  With matrix_pow and the various other
    matrix methods, a single attempt is made to balance heating and cooling and then
    a new ionization state is calculate based on the derived electron temperature.  With
    this method, an interation is performed, with this process being carried out multiple
    times. Currently number of iterations is fixed and set by the parameter NEBULARMODE_MATRIX_MULTISHOT
    which can be found in sirrocco.h.  This mode was developed to accelerate convergence in
    dense plasmas, but is still regarded as somewhat experimental as of 2025 Sept.  


File
  `setup.c <https://github.com/sirocco-rt/sirocco/blob/master/source/setup.c>`_


Child(ren)
  * :ref:`Wind.fixed_concentrations_file`

