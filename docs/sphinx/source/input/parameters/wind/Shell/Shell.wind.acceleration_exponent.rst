Shell.wind.acceleration_exponent
================================

This value is the :math:`\beta` exponent from the Caster & Lamers velocity equation of a stellar wind,

.. math::
  v(r) = v_0 + (v_\infty-v_0) * (1 - R_s/r)^{\beta}.

Type
  Double

Values
  Greater than or equal to 0

File
  `shell_wind.c <https://github.com/sirocco-rt/sirocco/blob/master/source/shell_wind.c>`_


Parent(s)
  * :ref:`Wind.type`: ``shell``


