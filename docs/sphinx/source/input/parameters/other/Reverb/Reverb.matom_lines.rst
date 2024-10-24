Reverb.matom_lines
==================
Number of macro-atom lines to track paths for individually. This many
reverb.matom_line entries are required, and the line associated with each has
the path of photons deexciting into it recorded in its own array. Note: This
doesn't give rise to any noticable differences to the pure wind mode in most
simulations.

Type
  Integer

Values
  Greater than or equal to 0

File
  `setup_reverb.c <https://github.com/sirocco-rt/sirocco/blob/master/source/setup_reverb.c>`_


Parent(s)
  * :ref:`Reverb.type`: matom

  * :ref:`Line_transfer`: ``macro_atoms``, ``macro_atoms_thermal_trapping``


Child(ren)
  * :ref:`Reverb.matom_line`

