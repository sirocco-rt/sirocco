Wind_heating.kpacket_frac
=========================
The user can set the fraction of additional k-packets (energy into the thermal pool) added during the macro-atom calculations. This subsequently reduces the number of photons emitted at each cycle. The default fraction is set to 0.1.

Type
  Double

Unit
  None

Values
  Greater than 0

File
  `setup.c <https://github.com/sirocco-rt/sirocco/blob/master/source/setup.c>`_


Parent(s)
  * :ref:`Wind_heating.extra_processes`: ``nonthermal``, ``both``

  * :ref:`Line_transfer`: ``macro_atoms``, ``macro_atoms_thermal_trapping``


