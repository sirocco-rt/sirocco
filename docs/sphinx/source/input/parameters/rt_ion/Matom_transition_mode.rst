Matom_transition_mode
================================

This option allows the user to determine which mode to use for the macro-atom transition probabilities. The macro atom is not a single atom; it represents an ensemble of atoms of a particular type. The mc_jumps method is a stochastic treatment and the matrix method is a state machine. The default is set to mc_jumping. More information on macro atoms can be found in the :ref:`macro atoms` section.

Matrix mode stores one matrix per wind cell, with :math:`(N+1)^2` entries for :math:`N` macro-atom levels in total, so its memory grows quickly with the size of the atomic data. The matrices are shared by the processes on a node, but even so, for 1000 cells they take about 60 MB with the 85-level H/He data (``h20_hetop_standard80``) and about 1.1 GB with the 377-level CNO data (``master_cno``). For large atoms or large grids, use ``mc_jumps``. The command-line option ``-no-matrix-storage`` avoids the memory, but then the matrix is recomputed every time a macro atom is activated, which is usually much slower than ``mc_jumps``.

Type
  Enumerator

Values
  mc_jumps
    A straightforward implementation of the macro-atom scheme. The models follow the interactions step by step, drawing a random number at each stage of the process until an interaction that excites a macro atom (or creates a k-packet directly) generates a new r-packet.

  matrix
    Matrix mode treats the macro atoms as a state machine with a finite number of states. State machines allow the computation of transition outcomes through matrix inversions rather than step-by-step tracking. Transitions to absorbing states
    correspond to transitions that generate k-packets. Transitions between levels of the macro atom via collisions with the thermal pool correspond to internal transitions of the absorbing state machine. The advantage of an absorbing state machine is that the new r-packet is not generated probabilistically but through several matrix inversions. All of the macro atoms and the thermal pool are treated simultaneously.

File
  `setup_line_transfer.c <https://github.com/sirocco-rt/sirocco/blob/master/source/setup_line_transfer.c>`_

Parent(s)
  * :ref:`Line_transfer`

