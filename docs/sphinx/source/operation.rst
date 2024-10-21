Code Operation
#####################

The basic code operation of SIROCCO is split into different cycles; First, the ionization state is calculated (:doc:`Ionization Cycles <operation/ionization_cycles>`). As these photons pass through the simulation grid, their heating and ionizing effect on the plasma is recorded through the use of Monte Carlo estimators. This process continues until the code converges on a solution in which the heating and cooling processes are balanced and the temperature stops changing significantly (see :doc:`Convergence & Errors <output/evaluation>`). Once the ionization and temperature structure of the outflow has been calculated, the spectrum is synthesized by tracking photons through the plasma until sufficient signal-to-noise is achieved in the output spectrum for lines to be easily identified (:doc:`Spectral Cycles <operation/spectral_cycles>`).

.. toctree::
   :glob:

   operation/*
