SV.v_zero_mode
==============

Option for how to set the wind velocity at the base of the wind. To be fixed or determined from the speed of sound in the wind.
The values for both these modes are set by :ref:`SV.v_zero`.

Type
  Enumerator

Values
  fixed
    A fixed speed for the initial velocity of the wind.

  sound_speed
    The initial velocity of the wind corresponds to the speed of sound in the wind.


File
  `sv.c <https://github.com/sirocco-rt/sirocco/blob/master/source/sv.c>`_


Parent(s)
  * :ref:`Wind.type`: ``SV``


Child(ren)
  * :ref:`SV.v_zero`

