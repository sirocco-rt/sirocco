Wind_heating.extra_luminosity
=============================

This options allows the user to determine the shock factor. The shock heating is defined initally as a luminosity to be added to wind but is immediately converted to a luminosity per unit volume. Since nearly all systems that we are dealing with have a star, we initialize the amount of extra heating as a fraction of the stellar luminosity. See `cooling.c <https://github.com/sirocco-rt/sirocco/blob/master/source/>`_ :code:`shock_heating`.

This is a very special option put in place for modelling FU Ori stars, and should be used with extreme caution.
The default value is calculated as a function of the star's radius and temperature. 

Type
  Double

Units
  ergs/s

Values
  Greater than 0

File
  `setup.c <https://github.com/sirocco-rt/sirocco/blob/master/source/setup.c>`_


Parent(s)
  * :ref:`Wind_heating.extra_processes`: ``nonthermal``, ``both``


