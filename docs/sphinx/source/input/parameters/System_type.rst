System_type
===========
The parameter is provides the program with a broad
overview of the type of system that will be simulated. Depending upon the user's choice of system type, the subsequent variables asked for generating the .pf file changes. The options are:

Type
  Enumerator

Values
  star
    A system in which the central object is a star.

  cv
    A system with central star and a companion star. The companion can occult the central object and disk, depending on the phase specified by the user. This will affect the outputted spectrum for an observer's position.

  bh
    A system with a black hole binary. The binary object, similarly to cv, can occult the central object and disk.

  agn
    A system in which the central object is a supermassive black hole.

  previous
    In this case, the user is starting from a previous SIROCCO run. The user can either continue the
    run or change some parameters associated with radiation sources.


File
  `sirocco.c <https://github.com/sirocco-rt/sirocco/blob/master/source/sirocco.c>`_


Child(ren)
  * :ref:`Atomic data`
  * :ref:`Binary.period`  
  * :ref:`Binary.mass_sec` 
  * :ref:`Boundary_layer.radiation`
  * :ref:`Central_object.geometry_for_source`  
  * :ref:`Central_object.temp`  
  * :ref:`Central_object.blackbody_temp`  
  * :ref:`Central_object.luminosity`  
  * :ref:`Spectrum.orbit_phase`
  * :ref:`Wind.number_of_components`
  * :ref:`Wind.old_windfile`