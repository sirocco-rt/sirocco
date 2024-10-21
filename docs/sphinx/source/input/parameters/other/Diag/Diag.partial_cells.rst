Diag.partial_cells
================================

Additional options for how to deal with cells that are partially filled by wind.
Somewhat degenerate with the ``-include_partial_cells`` flag under :ref:`Running SIROCCO`.


Type
  Enumerator

Values
  include
    Include wind cells that are only partially filled by the wind   

  zero_densities
    Ignore wind cells that are only partially filled by the wind by zeroing their density

  extend_full_cells
    Experimental model that extends full cells