SV.acceleration_exponent
========================
The power-law acceleration exponent (i.e. :math:`\alpha`) of a line driven wind in a Shlosman & Vitello (SV) CV disk wind model.
Sets the length scale over which the acceleration to :math:`v_{\infty}` is accomplished.
When equal to 1, the results resemble those of a linear velocity law.
Typically for an SV type wind, this power law exponent is 1.5.

.. math::
    v_l = v_o + (v_{\infty}(r_o)-v_o) \left[\frac {(l/R_v)^{\alpha}}{(l/R_v)^{\alpha}+1}\right]

Equation (2) Shlosman & Vitello 1993, ApJ 409, 372.

Type
  Double

Values
  Greater than 0

File
  `sv.c <https://github.com/sirocco-rt/sirocco/blob/master/source/sv.c>`_


Parent(s)
  * :ref:`Wind.type`: ``SV``


