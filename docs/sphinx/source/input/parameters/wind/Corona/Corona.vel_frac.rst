Corona.vel_frac
===============
For the coronal model, the azimuthal velocity is
set by the velocity of the underlying disk. The user
can also give the corona a radial velocity component. 
``vel_frac`` asks for the fraction of the disk's velocity to be used for the radial velocity.
(As coded, if this number is positive, the velocity is in the r direction
towards the central object).

Type
  Double

Unit
  :math:`v_{disk}` (in units of the disk's velocity)

Values
  Any, 0 implies no radial velocity.

File
  `corona.c <https://github.com/sirocco-rt/sirocco/blob/master/source/corona.c>`_


Parent(s)
  * :ref:`Wind.type`: ``corona``


