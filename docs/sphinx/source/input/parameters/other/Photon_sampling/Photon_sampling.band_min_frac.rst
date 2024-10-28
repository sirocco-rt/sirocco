Photon_sampling.band_min_frac
=============================
When specifying manually the bands used for generating photons during the ionization phase, this
parameter specifies the The minimum fraction of photons to be generated in this energy band.
The number of times this parameter will be reqested depends upon the number of bands.  The summ
of the fractions need not sum to 1, in which case the remaining photons will be distributed according
to the luminosity in the energy bands

Type
  Double

Values
  Greater than 0 and should sum to less than 1 over all bands

File
  `bands.c <https://github.com/sirocco-rt/sirocco/blob/master/source/bands.c>`_


Parent(s)
  * :ref:`Photon_sampling.nbands`: Greater than 0, once per band


