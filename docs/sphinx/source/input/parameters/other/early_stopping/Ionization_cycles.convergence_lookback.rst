@estop.Ionization_cycles.convergence_lookback
===============================================
The number of consecutive cycles over which the convergence fraction must be stable
(average absolute change less than the tolerance) before early stopping triggers.
Larger values are more conservative.

Type
  Integer

Values
  2 to 100

Default
  5

File
  `setup.c <https://github.com/sirocco-rt/sirocco/blob/dev/source/setup.c>`_

Parent(s)
  * ``-early_stopping`` switch
