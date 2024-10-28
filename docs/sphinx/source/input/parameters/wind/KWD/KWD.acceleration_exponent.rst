KWD.acceleration_exponent
=========================

The acceleration_exponent sets the length scale over which the accleration to :math:`v_\infty` is accomplished.
This value is the :math:`\beta` exponent from the Caster & Lamers equation of a stellar wind,

.. math::
  v(r) = v_0 + (v_\infty-v_0) * (1 - R_s/r)^{\beta}.

Type
  Double

Values
  Greater than 0

File
  `knigge.c <https://github.com/sirocco-rt/sirocco/blob/master/source/knigge.c>`_


Parent(s)
  * :ref:`Wind.type`: ``kwd``


