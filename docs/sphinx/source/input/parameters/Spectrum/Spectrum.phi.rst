Spectrum.phi
============
The azimuthal angle (in degrees) at which the spectrum observer is placed,
for non-binary systems.  This controls the phi component of the observer
line-of-sight direction and is relevant for non-axisymmetric 3-D models
(CYLIND3D or SPH3D) where the spectrum can vary with phi.

For binary systems the equivalent parameter is
:ref:`Spectrum.orbit_phase`, which sets the observer position via the
orbital phase (0–1).

A value of 0 places the observer in the x-z plane (phi = 0 rad), looking
in the −x direction.

Type
  Double

Values
  0 to 360 (degrees). Default 0.

File
  `setup.c <https://github.com/sirocco-rt/sirocco/blob/master/source/setup.c>`_


Parent(s)
  * :ref:`Spectrum.angle`: Once per observer angle, for non-binary systems only.

  * :ref:`System_type`: Not binary (i.e. not ``cv`` or ``bh``)
