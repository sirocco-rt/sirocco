Spectrum.live_or_die
====================
Normally in creating detailed spectrum SIROCCO "extracts" photons in a certain
direction reweighting them to account for the fact that they have been extracted
in a certain direction.  It is possible to just count the photons that are emitted
in a single angle range. The two methods should yield the same or very similar results
but the extraction method is much more efficient and live or die is basically a
diagnostic mode.

*Note: Small differences between the two methods can occur if a wind model has an 
edge a particular inclination angle, because in the live-or-die method, one must
accept photons in a small angle range, whereas the "extract" method involves photon packets
that emerge from the system at a single angle.*

Type
  Enumerator

Values
  live.or.die
    Count only those photons that escape within a small angle range towards the observer

  extract
    Extract a component of all photons that scatter towards the observer


File
  `setup.c <https://github.com/agnwinds/python/blob/master/source/setup.c>`_


Parent(s)
  * :ref:`Spectrum_cycles`: Greater than or equal to 0


