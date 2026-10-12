early_stopping
==============
Parameters controlling convergence-based early stopping of ionization cycles (accessed using the
``-early_stopping`` flag, see :ref:`Running SIROCCO`). When enabled, sirocco monitors the
convergence fraction across cycles and can terminate ionization early when both a stability
criterion and a minimum convergence floor are satisfied. The parameters are prefixed with
``@estop.`` in the parameter file.

.. toctree::
   :glob:

   early_stopping/*
