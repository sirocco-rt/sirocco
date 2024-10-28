KWD.mdot_r_exponent
===================
The exponent (:math:`\alpha`) for the mass loss rate as defined in the KWD model,

.. Math::

  \dot{m}(r) \propto F(r)^\alpha \propto T(r)^{4\alpha}

F is the local luminous flux and T is the local temperature at a radius R. A
value of 0 sets a uniform mass loss rate.

Type
  Double

Values
  Greater than or equal to 0

File
  `knigge.c <https://github.com/sirocco-rt/sirocco/blob/master/source/knigge.c>`_


Parent(s)
  * :ref:`Wind.type`: ``kwd``


