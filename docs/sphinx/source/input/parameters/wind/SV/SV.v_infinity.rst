SV.v_infinity
=============
Asymptotic (i.e. final) velocity of a line driven wind in a Shlosman & Vitello CV disk wind model.
Assumed to scale with the local velocity at the base of the streamline.

.. math::
    v_l = v_o + (v_{\infty}(r_o)-v_o) \left[\frac {(l/R_v)^{\alpha}}{(l/R_v)^{\alpha}+1}\right]

Equation (2) Shlosman & Vitello 1993, ApJ 409, 372.

Type
  Double

Unit
  :math:`v_{esc}` (in units of escape velocity)

Values
  Greater than 0

File
  `sv.c <https://github.com/sirocco-rt/sirocco/blob/master/source/sv.c>`_


Parent(s)
  * :ref:`Wind.type`: ``SV``


