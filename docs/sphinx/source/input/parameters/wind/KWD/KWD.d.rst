KWD.d
=====
The ratio :math:`d/d_{min}` is used to describe the degree of geometric collimation of
the disk wind in the KWD model. However, only d (the distance to the focal point in
central object radii) is set by the user as this provides a more natural parameter.

For high values of d, the wind is highly collimated, parallel to the rotational axis. 
For intermediate values of d (:math:`d\approx d_{min}`), the wind is more spherically symmetric.
For low values of d, the wind is 'sheet-like' equatorially. 

Type
  Double

Unit
  :ref:`Central_object.radius`

Values
  Greater than 0

File
  `knigge.c <https://github.com/sirocco-rt/sirocco/blob/master/source/knigge.c>`_


Parent(s)
  * :ref:`Wind.type`: ``kwd``


