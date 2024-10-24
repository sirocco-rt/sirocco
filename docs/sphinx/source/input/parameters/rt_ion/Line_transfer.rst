Line_transfer
=============
The way in which line transfer and scattering is dealt with
in the code. Governs whether we adopt any approximations
for radiative transfer, whether to use the indivisible packet
and macro-atom machinery, and whether to use isotropic or
anisotropic scattering.

Thermal trapping mode is recommended for non-macro atom runs,
while thermal trapping macro-atom mode is recommended for macro-atom runs.

Type
  Enumerator

Values
  pure_abs
    *Pure absorption*
    
    The pure absortion approximation.

  pure_scat
    *Pure scattering*
    
    The pure scattering approximation.

  sing_scat
    *Single scattering*
    
    The single scattering approximation.

  escape_prob
    *Escape probability*
    
    Resonance scattering and electron scattering is dealt with isotropically.
    free-free, compton and bound-free opacities attenuate the weight of the photon
    wind emission produces additional photons, which have their directions chosen isotropically.
    The amount of radiation produced is attenuated by the escape probability.

  thermal_trapping
    *Escape probability + anisotropic scattering*
    
    We use the 'thermal trapping method' to choose an
    anistropic direction when an r-packet deactivation
    or scatter occurs.

  macro_atoms
    *Macro-atoms*
    
    use macro-atom line transfer.
    Packets are indivisible and thus all opacities are dealt with by activate a macro-atom, scattering,
    or creating a k-packet.
    the new direction following electron scattering or deactivation of
    a macro atom is chosen isotropically.

  macro_atoms_thermal_trapping
    *Macro-atoms + anisotropic scattering*
    
    as macro_atoms, but we use the 'thermal trapping method' to choose an anistropic direction
    when an r-packet deactivation or scatter occurs.


File
  `setup_line_transfer.c <https://github.com/sirocco-rt/sirocco/blob/master/source/setup_line_transfer.c>`_


Child(ren)
  * :ref:`Reverb.matom_lines`

  * :ref:`Wind_heating.kpacket_frac`

