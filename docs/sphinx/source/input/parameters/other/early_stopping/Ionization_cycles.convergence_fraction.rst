@estop.Ionization_cycles.convergence_fraction
===============================================
The minimum percentage of cells that must be converged before early stopping is allowed.
This acts as a floor to prevent the code from stopping early on a model that is stably
unconverged (e.g. only 50% of cells converged). For example, a value of 80 means at least
80% of cells must be converged for the early stopping check to trigger.

Type
  Double

Values
  0 to 100

Default
  80

File
  `setup.c <https://github.com/sirocco-rt/sirocco/blob/dev/source/setup.c>`_

Parent(s)
  * ``-early_stopping`` switch
