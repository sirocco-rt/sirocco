Wind Model Parameters
######################

SIROCCO allows for various types of wind models, which are defined by the following parameters. This page focuses on the actual parameters in the file, but further description of the wind models and instructions on how to import models can be found under :doc:`Wind Models <../wind_models>`.

.. code::

  ### Parameters describing the various winds or coronae in the system
  Wind.number_of_components                         1
  Wind.type(SV,star,hydro,corona,kwd,homologous,shell,imported)                   SV
  Wind.coord_system(spherical,cylindrical,polar,cyl_var)          cylindrical
  Wind.dim.in.x_or_r.direction                     30
  Wind.dim.in.z_or_theta.direction                   30

:ref:`Wind.number_of_components` is usually 1, but can be greater if one wishes to construct a wind from a combination of several wind models,
for example a fast flow near the poles of a system, and a slow for near the disk.
If the number of components exceeds 1, then the remaining questions relating to the wind will be posed multiple times.

:ref:`Wind.type`: The wind models incorporated into SIROCCO currently are:

:ref:`SV`
  The Shlosman and Vitello parameterization of a bi-conical flow.

:ref:`Stellar_wind`
  A stellar-wind model. A fairly standard parameterization of a spherical outflow for a hot star.

:ref:`hydro`
  A special purpose mode used by the python collaboration for importing models from Zeus and Pluto. (Depreciated)

:ref:`Corona`
  A simple model for a corona above the disk.

:ref:`KWD`
   The Knigge Woods and Drew parameterization of a bi-conical flow.

:ref:`Homologous`
  A homologous expansion law useful for simulating SNe.

:ref:`Shell`
  A model of a thin shell useful for diagnostic studies.

:ref:`Imported <Importing models>`
  A general purpose mode for importing a wind from an ascii file (see also :doc:`SIROCCO Script documentation <../wind_models/importing_models>`).


:ref:`Wind.coord_system` is the coordinate system in which the wind is defined.

:ref:`Wind.dim.in.x_or_r.direction` is the number of grid cells in the x or r direction. 

:ref:`Wind.dim.in.z_or_theta.direction` is the number of grid cells in the z or theta direction.
