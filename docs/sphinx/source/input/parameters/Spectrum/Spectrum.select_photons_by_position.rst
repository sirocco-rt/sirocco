Spectrum.select_photons_by_position
===================================
Advanced command associated with adding conditions for
the detailed spectra that are extracted.  This command simply
asks whether one would like to construct spectra from photons
that originate or last scattered from a certain regions of space.

If yes, then one will be asked to specify the regions for all
extraction angles.

This option is useful for diagnostic purposes, such as differentiating
between photons that read the observer from the near or far side of
the disk.

*Note: This option is only available in extract mode. If one attempts to select
photons by position in live or die mode. The SIROCCO will warn the user and exit.*

Type
  Boolean (yes/no)

File
  `setup.c <https://github.com/agnwinds/python/blob/master/source/setup.c>`_


Parent(s)
  * :ref:`Spectrum_cycles`: Greater than or equal to 0


Child(ren)
  * :ref:`Spectrum.select_location`

