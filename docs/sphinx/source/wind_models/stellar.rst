The stellar wind model
############################################################


The stellar wind models implements the common 
`Caster & Larmers <https://ui.adsabs.harvard.edu/abs/1979ApJS...39..481C/abstract>`_
velocity law , where

.. math::
    v(r)=V_o + (V_{\infty}-V_o) (1-R_o/r)^{\beta}

Evidently, if :math:`\beta`
is 1, then the velocty will 
expand uniformly between :math:`V_o` and :math:`V_{\infty}`
