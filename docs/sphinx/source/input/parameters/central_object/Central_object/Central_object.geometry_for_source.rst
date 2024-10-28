Central_object.geometry_for_source
==================================


Setting the geometry model of the AGN/BH's radiation. This is only specified if :ref:`Central_object.radiation` is set to ``yes``. This is applicable even for black-body sources, where the *luminosity* depends on :ref:`Central_object.radius`.

Type
  Enumerator

Values
  lamp_post
    The central source radiates from two point sources
    located on the system axis above and below the disk plane.
    Emission is completely isotropic.

  sphere
    The central source radiates from a spherical surface with radius :ref:`Central_object.radius`.
    Emission is cosine-weighted in the direction of the sphere normal at the point of emission.
    Photons that would be spawned in an extended disk (if :ref:`Disk.type` is `vertically.extended`)
    are re-generated.

  bubble
    The central source radiates from random positions within a sphere of a given radius. 
    This model generates isotropic photons.

File
  `setup_star_bh.c <https://github.com/sirocco-rt/sirocco/blob/master/source/setup_star_bh.c>`_


Parent(s)
  * :ref:`System_type`: ``agn``, ``bh``

  * :ref:`Central_object.radiation`: ``True``


Child(ren)
  * :ref:`Central_object.lamp_post_height`

