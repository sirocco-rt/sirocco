Input_spectra.model_file
========================
In addition to being able to generate several types of spectra, such
as blackbodies and power laws, SIROCCO can read in a series of spectra
which are tabulated and are describable in terms of (usually) temperature
and gravity). This parameter defines the name of the file which gives the
location of the individual spectra and the temperate and gravity associated
with each spectrum. (One may wish to use the same files for several radiation sources, viz the disk and the star).
SIROCCO actually only reads in the data the first time.

Type
  String

File
  `setup.c <https://github.com/agnwinds/python/blob/master/source/setup.c>`_


Parent(s)
  * :ref:`Central_object.rad_type_to_make_wind`: models

  * :ref:`Central_object.rad_type_in_final_spectrum`: models

  * :ref:`Disk.rad_type_to_make_wind`: models

  * :ref:`Disk.rad_type_in_final_spectrum`: models

  * :ref:`Boundary_layer.rad_type_to_make_wind`: models

  * :ref:`Boundary_layer.rad_type_in_final_spectrum`: models


