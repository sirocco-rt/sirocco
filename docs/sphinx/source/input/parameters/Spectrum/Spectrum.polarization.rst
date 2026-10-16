Spectrum.polarization
=====================
Advanced command that switches on the calculation of linear polarization in the
spectral cycles.  When it is ``yes``, photon packets carry the Stokes parameters
q = Q/I and u = U/I, electron scattering is treated with the polarized
(Thomson or Klein-Nishina) scattering matrix, and Q and U spectra are produced
for each observer.  Two columns, e.g. ``A80P0.50_Q`` and ``A80P0.50_U``, are
added for each observer to the ``.spec``, ``.spec_wind``, ``.log_spec`` and
``.log_spec_wind`` files.  Q and U are referred to North on the sky, the
projected symmetry axis of the system, so Q > 0 means the electric vector is
parallel to that axis.

The ionization cycles are not affected.  At present only electron scattering
polarizes the light; line scattering and all other emission processes are
treated as unpolarized.  See :doc:`/physics/polarization`.

Type
  Boolean (yes/no)

Values
  ``no`` (the default) or ``yes``

File
  `setup.c <https://github.com/sirocco-rt/sirocco/blob/master/source/setup.c>`_


Parent(s)
  * :ref:`Spectrum_cycles`: Greater than or equal to 0


