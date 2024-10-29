Diag.save_extract_photons
=========================
Saves output to ``<file_name>.ext.txt``. This is a diagnostic that allows one to print out information about a photon at any time. 
Extract is called when a photon begins its flight and every time that photon scatters.
Photons are saved before and after each Compton step determining the direction for a photon undergoing Compton scattering. 
Save_extract_photons uses the ``save_photons`` function.

The ``save_photons`` columns in the output are assigned in order as follows:

.. code:: 

  geo.pcycle, 
  p->np, 
  p->freq_orig,
  p->freq, 
  p->w_orig,
  p->w,
  p->x[0],
  p->x[1],
  p->x[2],
  p->lmn[0],
  p->lmn[1],
  p->lmn[2],
  p->ds,
  p->tau,
  p->grid,
  p->istat,
  p->origin,
  p->nscat,
  p->nres,
  p->frame,
  comment;

Type
  Boolean (yes/no)

File
  `diag.c <https://github.com/sirocco-rt/sirocco/blob/master/source/diag.c>`_


Parent(s)
  * :ref:`Diag.extra`: ``True``


