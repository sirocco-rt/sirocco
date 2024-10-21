System Description
##################

The first set of parameters which SIROCCO needs are information about the overall system.

.. code::
   
   System_type(star,cv,bh,agn,previous)                   bh

   ### Parameters for the Central Object
   Central_object.mass(msol)                        10
   Central_object.radius(cm)               8.85667e+06
   Binary.mass_sec(msol)                           0.4
   Binary.period(hr)                               3.2

   ### Parameters for the Disk (if there is one)
   Disk.type(none,flat,vertically.extended,rmin>central.obj.rad)                 flat
   Disk.radiation(yes,no)                          yes
   Disk.rad_type_to_make_wind(bb,models,mod_bb)                   bb
   Disk.temperature.profile(standard,readin)             standard
   Disk.mdot(msol/yr)                            1e-08
   Disk.radmax(cm)                           2.657e+08

   ### Parameters for Boundary Layer or the compact object in an X-ray Binary or AGN
   Central_object.radiation(yes,no)                  yes
   Central_object.rad_type_to_make_wind(bb,models,power,cloudy,brems,mono)                power
   Central_object.luminosity(ergs/s)          4.72063e+37
   Central_object.power_law_index                 -1.5
   Central_object.geometry_for_source(sphere,lamp_post,bubble)               sphere

:ref:`System_type` is starting point, a basic classification of the type of object one is trying to model.
This is used to guide further questions about the object and to set defaults.

Most of the other parameters are fairly self-explanatory, and are documented fully in the various :ref:`Parameters` entries.
