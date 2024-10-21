Photon_sampling.nbands
======================
SIROCCO uses stratified samplign to generate photons during the ionization phase.  This
parameter allows the user to define the number of bands for stratified sampling, if s/he
wants to customize the bands used for the generation of photons

Type
  Integer

Values
  Greater than 0

File
  `bands.c <https://github.com/agnwinds/python/blob/master/source/bands.c>`_


Parent(s)
  * :ref:`Photon_sampling.approach`: user_bands


Child(ren)
  * :ref:`Photon_sampling.band_min_frac`

  * :ref:`Photon_sampling.band_boundary`

