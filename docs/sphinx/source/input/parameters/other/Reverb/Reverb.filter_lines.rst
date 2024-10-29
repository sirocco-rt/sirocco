Reverb.filter_lines
===================
Whether or not to filter any lines out of the output file. This is used to keep output
file sizes down, and avoid them overwhelming the user.

Type
  Int

Values
  0
    *No filtering*
    
    Include *all* photons that contribute to the spectra in the output
    file. Not recommended as it leads to gargantuan file sizes.

  -1
    *Filter continuum*
    
    Include all photons whose last interaction was scatter
    or emission in a line. Recommended setting for exploratory runs where you'd
    like to identify which lines are the easiest to process.

  N
    *Filter lines*
    
    Include N :ref:`Reverb.filter_line` entries, each specifying one
    line to keep in the output file. If :ref:`reverb.matom_lines` is >0, all macro-atom
    lines of interest are automatically included in the filter list.


File
  `setup_reverb.c <https://github.com/sirocco-rt/sirocco/blob/master/source/setup_reverb.c>`_


Parent(s)
  * :ref:`Reverb.type`: ``wind``, ``matom``


Child(ren)
  * :ref:`Reverb.filter_line`

