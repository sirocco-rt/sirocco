@estop.Ionization_cycles.min_cycles
=====================================
The minimum number of ionization cycles to complete before early stopping is allowed.
This prevents premature termination before the wind has had a chance to settle.

Type
  Integer

Values
  Greater than or equal to 0

Default
  0

File
  `setup.c <https://github.com/sirocco-rt/sirocco/blob/dev/source/setup.c>`_

Parent(s)
  * ``-early_stopping`` switch
