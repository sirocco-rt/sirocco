Diag.track_resonant_scatters
============================

Saves output to ``<file_name>.ext.txt``. This is a diagnostic that allows one to print out information about a
photon that has gone through a resonant scatter during transport through the wind. 
track_resonant_scatters uses the ``track_scatters`` function.

The ``track_scatters`` columns in the output are assigned in order as follows:

.. code:: 

   "Scatter",
    p->np,
    p->x[0], 
    p->x[1], 
    p->x[2], 
    p->grid, 
    p->freq, 
    p->w, 
    nplasma, 
    comment

Type
  Boolean (yes/no)

File
  `diag.c <https://github.com/sirocco-rt/sirocco/blob/master/source/diag.c>`_


Parent(s)
  * :ref:`Diag.extra`: ``True``


