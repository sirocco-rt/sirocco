Reverb.type
===========
Whether to perform reverberation mapping. Reverberation mapping tracks the
path of photons emitted in the simulation as they travel through the geometry,
assuming that any delays from recombination etc. are negligible and all delays
are due to light travel time. For each final spectrum, all contributing
photons are output to a '.delay_dump' file that can then be processed using
our 'tfpy' SIROCCO (no relation) library.

Type
  Enumerator

Values
  none
    **Off**

  photon
    Each photon is assigned an initial path based on its distance from the
    central source (assuming emission in the disk and wind is correlated with
    emission from the CO).

  wind
    CO photons are assigned paths as in Photon mode, disk photons are assigned
    paths as set by the reverb.disk_type parameter. Photons generated in the
    wind are assigned a path based on the *distribution* of paths of photons
    that have contributed to continuum absorption in that cell.

  matom
    This works as wind mode, but for a number of specified macro-atom lines
    paths are tracked for those photons who cause a deexcitation into a given
    line. When a photon is emitted in one of those lines, the path is drawn from
    that specific distribution. This distribution is build up not just from the
    last cycle of the simulation, but from all cycles after the wind achieves
    >90% convergence. This is necessary as some lines are poorly-sampled.
    
    This mode gives pretty much identical results to wind mode, but at least we
    made it to check rather than just assuming it would be fine.
    
    This requires that :ref:`Line_transfer` is either ``macro_atoms`` or 
    ``macro_atoms_thermal_trapping``


File
  `setup_reverb.c <https://github.com/agnwinds/python/blob/master/source/setup_reverb.c>`_


Child(ren)
  * :ref:`Reverb.matom_lines`

  * :ref:`Reverb.filter_lines`

  * :ref:`Reverb.path_bins`

  * :ref:`Reverb.visualisation`

  * :ref:`Reverb.disk_type`

