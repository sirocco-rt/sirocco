Wind_heating.extra_processes
============================

The user can selected setup options associated with non radiative processes that can affect the thermal balance of the wind. These are adiabatic processes and an extra heating term explicitly implemented for FU Ori stars.  The default is set to adiabatic. 

Type
  Enumerator

Values
  none
    This mode disables any extra heating processes.

  adiabatic
    This mode enables adiabatic cooling processes which contributes to thermal balance within the wind.

  nonthermal
   This is a very special option put in place for modelling FU Ori stars, and should be used with extreme caution.

   This mode enables shock heating into the thermal balance of the wind. The shock heating is defined initally as a luminosity to be added to wind but is immediately converted to a luminosity per unit volume (:ref:`Wind_heating.extra_luminosity`). Since nearly all systems that we are dealing with have a star, we initialize the amount of extra heating as a fraction of the stellar luminosity. See `cooling.c <https://github.com/sirocco-rt/sirocco/blob/master/source/>`_ :code:`shock_heating`.

  both
    This is a very special option put in place for modelling FU Ori stars, and should be used with extreme caution.

    The adiabatic and nonthermal cases from above combined.


File
  `setup.c <https://github.com/sirocco-rt/sirocco/blob/master/source/setup.c>`_


Child(ren)
  * :ref:`Wind_heating.extra_luminosity`
  * :ref:`Wind_heating.kpacket_frac` - Requires macro-atom modes (:ref:`Line_transfer` ``macro_atoms``, ``macro_atoms_thermal_trapping``).


