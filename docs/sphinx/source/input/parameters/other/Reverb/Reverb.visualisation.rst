Reverb.visualisation
====================
Which type of visualisation to output, if any. Reverb modes that keep arrays
of photon paths per cell can output them either as averages in a 3d model, or
as a selection of flat text files with full bin-by-bin breakdowns. Useful for
diagnostics.

Type
  Enumerator

Values
  none
    No visualisation.

  vtk
    Mesh visualisation. Outputs mean incident path per cell, photon count per cell, and mean
    observed delay to '.vtk' format, readable using a range of programs including
    (my preferred option) VisIt, available at https://visit.llnl.gov/.

  dump
    Outputs distributions of paths for continuum heating and each line to a range of 'dump cells'
    specified by X & Z position.

  both
    Outputs both vtk and dump.


File
  `setup_reverb.c <https://github.com/sirocco-rt/sirocco/blob/master/source/setup_reverb.c>`_


Parent(s)
  * :ref:`Reverb.type`: ``wind``, ``matom``


Child(ren)
  * :ref:`Reverb.dump_cells`

  * :ref:`Reverb.angle_bins`

