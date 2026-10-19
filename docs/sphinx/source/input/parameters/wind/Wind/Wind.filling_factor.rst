Wind.filling_factor
===================
The volume filling factor of the outflow. The implementation
of clumping (microclumping) is described in
Matthews et al. (2016), 2016MNRAS.458..293M. Asked once per domain.

The density given by the wind model (or in an imported model) is the
volume-averaged density.  SIROCCO divides it by the filling factor to obtain
the density inside the clumps, which is the density it stores and uses, and
multiplies the cell volume by the filling factor.

Type
  Double

Values
  :math:`0\lt f\le1`, where 1 is a fully smooth wind.

File
  `setup_domains.c <https://github.com/sirocco-rt/sirocco/blob/master/source/setup_domains.c>`_


Parent(s)
  * :ref:`Wind.number_of_components`: Greater than 0. Once per domain.


