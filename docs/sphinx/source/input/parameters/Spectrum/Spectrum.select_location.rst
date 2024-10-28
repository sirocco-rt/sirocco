Spectrum.select_location
========================
One of several related parameters that permit one to apply
additional conditions on the location of photons extracted in
the detailed spectrum. The location refers here to the either
where the photons was created or where it last scattered

Type
  Enumerator

Values
  all
    Select photons regardless of where they are generated

  below_disk
    Select only photons generated from below (-z) the disk

  above_disk
    Select only photons orginating above the disk

  spherical_region
    Select photons by defining a spherical region


File
  `setup.c <https://github.com/sirocco-rt/sirocco/blob/master/source/setup.c>`_


Parent(s)
  * :ref:`Spectrum.select_photons_by_position`: ``True``


Child(ren)
  * :ref:`Spectrum.select_r`

  * :ref:`Spectrum.select_rho`

  * :ref:`Spectrum.select_azimuth`

  * :ref:`Spectrum.select_z`

