Photons_per_cycle
=================
This parameter allows the user to set the number of photons per cycle. The photon value given is used for both ionization and spectral cycles. The total number of photons passed through a simulation is the photons_per_cycle multiplied by the number of cycles. Higher numbers of photons_per_cycle will result in higher memory usage. However, counteracting a reduction in photons_per_cycle by adding more cycles will result in more overhead and CPU time from non-parallelized code. The user should find a balance between the photon number per cycle and the number of cycles. The default value given is 100,000.

Type
  Double

Values
  Greater than 0

File
  `setup.c <https://github.com/sirocco-rt/sirocco/blob/master/source/setup.c>`_


