@estop.Ionization_cycles.convergence_tolerance
================================================
Percentage tolerance for early stopping of ionization cycles. If the average cycle-to-cycle
change in the convergence fraction over the lookback window falls below this percentage *and*
the convergence fraction is at or above the floor set by
:ref:`@estop.Ionization_cycles.convergence_fraction`, ionization stops early and proceeds to
spectrum extraction. For example, a value of 2 means stop when the average change is less than 2%.

Type
  Double

Values
  Greater than 0

Default
  2

File
  `setup.c <https://github.com/sirocco-rt/sirocco/blob/dev/source/setup.c>`_

Parent(s)
  * ``-early_stopping`` switch
