SV.mdot_r_exponent
==================
The exponent :math:`\lambda` for the mass loss rate per unit area as defined in the Shlosman Vitelo model. 

.. math::        
    \frac{\delta\dot{m}}{\delta A} \propto \dot{m}_{wind} r_o^{\lambda} cos(\theta(r_o))

Equation (4) Shlosman & Vitelo,ApJ,1993,409,372.

Type
  Double

Values
  Greater than or equal to 0. 0 sets a uniform mass loss rate.

File
  `sv.c <https://github.com/sirocco-rt/sirocco/blob/master/source/sv.c>`_


Parent(s)
  * :ref:`Wind.type`: ``SV``


