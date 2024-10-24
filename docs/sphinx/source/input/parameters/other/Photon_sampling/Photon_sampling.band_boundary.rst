Photon_sampling.band_boundary
=============================
When the user specifies what bands are used for stratfied sampling, this parameter specifies the boundaries
between energy bands in which a minimum fraction of photons will be generated.  The number of times this
parameter is request depends upon the number of energies bands being used.

Type
  Double

Unit
  eV

Values
  Greater than 0, monotonically increasing

File
  `bands.c <https://github.com/sirocco-rt/sirocco/blob/master/source/bands.c>`_


Parent(s)
  * :ref:`Photon_sampling.nbands`: Greater than 0, once per band


