SV.acceleration_length
======================
The wind acceleration height scale (:math:`R_v`) for a disk wind described by the
Shlosman Vitelo model.

.. math::
    v_l = v_o + (v_{\infty}(r_o)-v_o) \left[\frac {(l/R_v)^{\alpha}}{(l/R_v)^{\alpha}+1}\right]

Equation (2) Shlosman & Vitello 1993, ApJ 409, 372.

Type
  Double

Unit
  cm

Values
  Greater than 0

File
  `sv.c <https://github.com/sirocco-rt/sirocco/blob/master/source/sv.c>`_


Parent(s)
  * :ref:`Wind.type`: ``SV``


