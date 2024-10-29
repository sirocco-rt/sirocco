KWD.v_infinity
==============

A parameter controlling the radial streamline velocities of a steller wind at large distances. Described by the KWD model and given in terms of Castor & Lamers equation,

.. math:: 
  v(r) = v_0 + (v_\infty - v_0) * (1 - R_s/r)^\beta.


Type
  Double

Unit
  :math:`v_{esc}` (in units of escape velocity)

Values
  Greater than 0

File
  `knigge.c <https://github.com/sirocco-rt/sirocco/blob/master/source/knigge.c>`_


Parent(s)
  * :ref:`Wind.type`: ``kwd``


